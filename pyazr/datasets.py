"""Dataset and extrapolation provenance for an AZURE2 model, from the ``.azr``.

The ``<segmentsData>`` block records, per data segment, where the data came from
(the data file), which reaction channel it measures (entrance / exit particle
pairs), its energy and angle range, the observable type, and the normalization
systematic error.  The ``<segmentsTest>`` block declares the *extrapolations*:
grids of energies / angles AZURE2 evaluates the model on with no data attached.
AZURE2's API exposes the *values* (energies, cross sections) but not this
provenance, so -- exactly as the fitting examples already read the per-segment
systematic error straight from the file -- this module parses both blocks in
Python.  No C++/API change is needed.

Layout of a ``<segmentsData>`` line (see ``include/SegLine.h``)::

    isActive entranceKey exitKey minE maxE minA maxA isDiff
        [phaseJ phaseL if isDiff==2]
        dataNorm varyNorm dataNormError
        [energyShift energyShiftError varyEnergyShift]      (optional)
        dataFile  [advanced flags ...]

Layout of a ``<segmentsTest>`` line (see ``include/ExtrapLine.h``)::

    isActive entranceKey exitKey minE maxE eStep minA maxA aStep isDiff
        [phaseJ phaseL      if isDiff==2]
        [maxAngDistOrder    if isDiff==3]
        [isAdvanced ...]                                    (optional)

Note the two blocks differ in more than the step columns: an extrapolation line
carries ``eStep``/``aStep``, and the ``isDiff`` codes are *not* the same set --
``3`` is total capture in a data segment but an angular distribution in an
extrapolation, and cm-differential is ``4`` vs ``5``.  Both maps below are taken
from ``ESegment``'s two constructors (``src/ESegment.cpp``), which is what
AZURE2 itself acts on.
"""

import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


def _isfloat(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


# isDiff codes for <segmentsData> -- ESegment::ESegment(SegLine).
_OBSERVABLE = {
    0: "angle-integrated",
    1: "differential",
    2: "phase-shift",
    3: "total-capture",
    4: "differential-cm",
    5: "angle-integrated-E1",
    6: "angle-integrated-E2",
    7: "analyzing-power",
}

# isDiff codes for <segmentsTest> -- ESegment::ESegment(ExtrapLine).
_EXTRAP_OBSERVABLE = {
    0: "angle-integrated",
    1: "differential",
    2: "phase-shift",
    3: "angular-distribution",
    4: "total-capture",
    5: "differential-cm",
    7: "analyzing-power",
}


@dataclass
class Segment:
    """One data segment: its provenance and what it measures."""

    key: int                       # 1-based segment number (file order)
    active: bool
    entrance_key: int              # entrance particle-pair key
    exit_key: int                  # exit pair key (-1 = total / summed)
    energy_min: float
    energy_max: float
    angle_min: float
    angle_max: float
    observable: str                # 'differential' | 'angle-integrated' | ...
    norm: float                    # applied normalization
    vary_norm: bool                # is the normalization a fit parameter?
    norm_error: float              # normalization systematic, fractional
    data_file: str                 # path the data came from  <-- provenance
    energy_shift: float = 0.0      # applied beam-energy shift (MeV)
    energy_shift_error: float = 0.0    # its systematic (MeV)
    vary_shift: bool = False       # is the energy shift a fit parameter?
    operation: Optional[str] = None    # 'sum' | 'ratio' for a composite segment
    components: tuple = ()             # (entrance, exit, angle, scaling) each

    @property
    def composite(self) -> bool:
        """True if this segment is a sum or ratio of several pathways."""
        return bool(self.components)

    @property
    def name(self) -> str:
        """A short label: the data file's base name without extension."""
        return os.path.splitext(os.path.basename(self.data_file))[0]

    def reaction(self, pairs=None) -> str:
        """Human-readable reaction, e.g. ``'pair1 -> pair2'`` (or with pair
        spins if a :class:`~pyazr.PairSet` is supplied)."""
        ex = "total" if self.exit_key == -1 else f"pair{self.exit_key}"
        return f"pair{self.entrance_key} -> {ex}"

    def describe(self) -> str:
        """The reaction, spelled out for a composite segment.

        The segment's own entrance/exit pair is the first term: for a ratio it
        is the numerator and the components are the denominator, for a sum it is
        simply the first addend.
        """
        def term(e, x, a=-999.0, s=1.0):
            out = f"pair{e}->" + ("total" if x == -1 else f"pair{x}")
            if a > -900:
                out += f"@{a:g}deg"
            # -999 is the "absent" sentinel, the same one the angle uses.
            if s is not None and s > -900.0 and s != 1.0:
                out += f" x{s:g}"
            return out

        head = term(self.entrance_key, self.exit_key)
        if not self.composite:
            return head
        rest = [term(*c) for c in self.components]
        if self.operation == "ratio":
            return head + " / " + " / ".join(rest)
        return " + ".join([head] + rest)

    def __repr__(self):
        extra = f", {self.operation} of {len(self.components)}" if self.composite else ""
        return (f"Segment(#{self.key} {self.name!r}, "
                f"{self.entrance_key}->{self.exit_key}, {self.observable}, "
                f"E={self.energy_min:g}-{self.energy_max:g}, "
                f"norm_err={self.norm_error:g}{extra})")


class SegmentSet(list):
    """All data segments of a model, with provenance and convenience views."""

    @classmethod
    def from_file(cls, path):
        with open(path) as f:
            lines = f.readlines()
        try:
            start = lines.index("<segmentsData>\n") + 1
            end = lines.index("</segmentsData>\n")
        except ValueError:
            raise ValueError(f"{path}: no <segmentsData> block.")
        segs = cls()
        key = 0
        for raw in lines[start:end]:
            if not raw.strip():
                continue
            key += 1
            segs.append(cls._parse(raw.split(), key))
        return segs

    @staticmethod
    def _parse(t, key):
        active = int(float(t[0])) == 1
        entrance_key = int(float(t[1]))
        exit_key = int(float(t[2]))
        eMin, eMax, aMin, aMax = (float(t[3]), float(t[4]),
                                  float(t[5]), float(t[6]))
        isDiff = int(float(t[7]))
        i = 8
        if isDiff == 2:                 # phase-shift carries J, L
            i += 2
        norm = float(t[i]); vary_norm = int(float(t[i + 1])) == 1
        norm_error = float(t[i + 2]); i += 3
        # optional energy-shift triple, then the (non-numeric) data file
        shift = shift_error = 0.0
        vary_shift = False
        if i + 2 < len(t) and all(_isfloat(tok) for tok in t[i:i + 3]):
            shift, shift_error = float(t[i]), float(t[i + 1])
            vary_shift = int(float(t[i + 2])) == 1
            i += 3
        while i < len(t) and _isfloat(t[i]):
            i += 1
        data_file = t[i] if i < len(t) else ""

        # A composite ("advanced") segment continues after the data file with
        #   isAdvanced operationType nComponents [entrance exit angle [scaling]]*
        # and a leading -1 in the count position marks the newer layout that
        # carries a per-component scaling factor.  See SegLine.h.
        operation, components = None, []
        tail = t[i + 1:]
        if tail and _isfloat(tail[0]) and int(float(tail[0])) == 1:
            operation = "ratio" if (len(tail) > 1 and
                                    int(float(tail[1])) == 1) else "sum"
            j, scaled = 2, False
            n = int(float(tail[j])) if len(tail) > j else 0
            j += 1
            if n == -1:                      # marker: per-component scaling
                scaled = True
                n = int(float(tail[j])) if len(tail) > j else 0
                j += 1
            for _ in range(n):
                if j + 2 >= len(tail) + 1 and not scaled:
                    break
                try:
                    e, x, a = (int(float(tail[j])), int(float(tail[j + 1])),
                               float(tail[j + 2]))
                except (IndexError, ValueError):
                    break
                j += 3
                s = 1.0
                if scaled and j < len(tail) and _isfloat(tail[j]):
                    s = float(tail[j]); j += 1
                components.append((e, x, a, s))

        return Segment(
            key=key, active=active, entrance_key=entrance_key,
            exit_key=exit_key, energy_min=eMin, energy_max=eMax,
            angle_min=aMin, angle_max=aMax,
            observable=_OBSERVABLE.get(isDiff, f"code{isDiff}"),
            norm=norm, vary_norm=vary_norm,
            norm_error=norm_error, data_file=data_file,
            energy_shift=shift, energy_shift_error=shift_error,
            vary_shift=vary_shift,
            operation=operation, components=tuple(components))

    # -- views ----------------------------------------------------------------

    @property
    def active(self):
        """Only the active data segments."""
        return SegmentSet(s for s in self if s.active)

    def by_reaction(self, entrance_key=None, exit_key=None):
        """The data segments for one entrance and exit pair."""
        return SegmentSet(
            s for s in self
            if (entrance_key is None or s.entrance_key == entrance_key)
            and (exit_key is None or s.exit_key == exit_key))

    def by_key(self, key):
        """The segment with this segment key, or None.

        A key is the segment's position in the input file counting the
        *inactive* ones too, so it is not an index into this set -- which is
        what makes ``datasets[key - 1]`` wrong for any project that has one.
        A ``Parameter``'s ``segment_key`` is a key in this sense.
        """
        for s in self:
            if s.key == key:
                return s
        return None

    def sys_errors(self, active_only=True, vary_only=False, fractional=True):
        """Per-segment normalization systematic errors, in file order.

        ``fractional`` divides the file's percent value by 100 (matching the
        convention the fitting examples use).  ``vary_only`` restricts to
        segments whose normalization is actually a fit parameter.
        """
        src = self.active if active_only else self
        if vary_only:
            src = SegmentSet(s for s in src if s.vary_norm)
        f = 1e-2 if fractional else 1.0
        return [s.norm_error * f for s in src]

    def shift_errors(self, active_only=True, vary_only=False):
        """Per-segment beam-energy-shift errors (MeV), in file order.

        The energy-shift counterpart of :meth:`sys_errors`; ``vary_only``
        restricts to segments whose shift is actually a fit parameter.  Unlike
        the normalization systematic, the file value is already absolute (MeV),
        so there is nothing to convert.
        """
        src = self.active if active_only else self
        if vary_only:
            src = SegmentSet(s for s in src if s.vary_shift)
        return [s.energy_shift_error for s in src]

    def files(self):
        """The data files the segments read."""
        return [s.data_file for s in self]

    def table(self, pairs=None):
        """A printable table of the data segments."""
        rows = [("#", "data file", "reaction", "observable", "E range",
                 "norm_err%", "vary")]
        for s in self:
            rows.append((
                str(s.key), s.name, s.reaction(pairs), s.observable,
                f"{s.energy_min:g}-{s.energy_max:g}",
                f"{s.norm_error:g}", "*" if s.vary_norm else ""))
        w = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        return "\n".join("  ".join(cell.ljust(w[c]) for c, cell in enumerate(r))
                         for r in rows)

    def __repr__(self):
        return f"SegmentSet({len(self)} segments, {len(self.active)} active)"


@dataclass
class TestSegment:
    """One extrapolation segment: a grid to evaluate the model on, no data.

    The energy / angle grids are declared as ``min``, ``max``, ``step``; AZURE2
    generates the points itself, so the calculated arrays come back from the API
    (``azure2.calculate_energies`` and friends in extrapolation mode) rather than
    from any file.
    """

    key: int                       # 1-based segment number (file order)
    active: bool
    entrance_key: int              # entrance particle-pair key
    exit_key: int                  # exit pair key (-1 = total / summed)
    energy_min: float
    energy_max: float
    energy_step: float
    angle_min: float
    angle_max: float
    angle_step: float
    observable: str                # 'differential' | 'total-capture' | ...
    phase_J: Optional[float] = None        # isDiff == 2 only
    phase_L: Optional[int] = None          # isDiff == 2 only
    max_ang_dist_order: Optional[int] = None   # isDiff == 3 only

    @property
    def name(self) -> str:
        """A short label, e.g. ``'1->2'`` or ``'1->total'``."""
        ex = "total" if self.exit_key == -1 else str(self.exit_key)
        return f"{self.entrance_key}->{ex}"

    @property
    def is_angle_integrated(self) -> bool:
        """Is the observable angle-integrated rather than differential?"""
        return self.observable in ("angle-integrated", "total-capture")

    @property
    def n_energies(self) -> int:
        """Number of energy points AZURE2 will generate for this grid."""
        if self.energy_step <= 0:
            return 1
        return int((self.energy_max - self.energy_min) / self.energy_step) + 1

    def reaction(self, pairs=None) -> str:
        """Human-readable reaction, e.g. ``'pair1 -> pair2'``."""
        ex = "total" if self.exit_key == -1 else f"pair{self.exit_key}"
        return f"pair{self.entrance_key} -> {ex}"

    def __repr__(self):
        return (f"TestSegment(#{self.key} {self.entrance_key}->{self.exit_key}, "
                f"{self.observable}, "
                f"E={self.energy_min:g}-{self.energy_max:g}"
                f"/{self.energy_step:g})")


class TestSegmentSet(list):
    """All extrapolation segments of a model, parsed from ``<segmentsTest>``."""

    @classmethod
    def from_file(cls, path):
        with open(path) as f:
            lines = f.readlines()
        try:
            start = lines.index("<segmentsTest>\n") + 1
            end = lines.index("</segmentsTest>\n")
        except ValueError:
            # A model need not declare any extrapolations.
            return cls()
        segs = cls()
        key = 0
        for raw in lines[start:end]:
            if not raw.strip():
                continue
            key += 1
            segs.append(cls._parse(raw.split(), key))
        return segs

    @staticmethod
    def _parse(t, key):
        active = int(float(t[0])) == 1
        entrance_key = int(float(t[1]))
        exit_key = int(float(t[2]))
        eMin, eMax, eStep = float(t[3]), float(t[4]), float(t[5])
        aMin, aMax, aStep = float(t[6]), float(t[7]), float(t[8])
        isDiff = int(float(t[9]))
        phase_J = phase_L = max_order = None
        i = 10
        if isDiff == 2:                 # phase shift carries J, L
            phase_J = float(t[i]); phase_L = int(float(t[i + 1])); i += 2
        elif isDiff == 3:               # angular distribution carries its order
            max_order = int(float(t[i])); i += 1
        return TestSegment(
            key=key, active=active, entrance_key=entrance_key,
            exit_key=exit_key,
            energy_min=eMin, energy_max=eMax, energy_step=eStep,
            angle_min=aMin, angle_max=aMax, angle_step=aStep,
            observable=_EXTRAP_OBSERVABLE.get(isDiff, f"code{isDiff}"),
            phase_J=phase_J, phase_L=phase_L, max_ang_dist_order=max_order)

    # -- views ----------------------------------------------------------------

    @property
    def active(self):
        """Only the active test segments."""
        return TestSegmentSet(s for s in self if s.active)

    def by_reaction(self, entrance_key=None, exit_key=None):
        """The test segments for one entrance and exit pair."""
        return TestSegmentSet(
            s for s in self
            if (entrance_key is None or s.entrance_key == entrance_key)
            and (exit_key is None or s.exit_key == exit_key))

    def by_key(self, key):
        """The test segment with this segment key, or None."""
        for s in self:
            if s.key == key:
                return s
        return None

    def table(self, pairs=None):
        """A printable table of the test segments."""
        rows = [("#", "reaction", "observable", "E range", "E step",
                 "angle", "active")]
        for s in self:
            angle = ("-" if s.is_angle_integrated
                     else f"{s.angle_min:g}-{s.angle_max:g}")
            rows.append((
                str(s.key), s.reaction(pairs), s.observable,
                f"{s.energy_min:g}-{s.energy_max:g}", f"{s.energy_step:g}",
                angle, "*" if s.active else ""))
        w = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        return "\n".join("  ".join(cell.ljust(w[c]) for c, cell in enumerate(r))
                         for r in rows)

    def __repr__(self):
        return (f"TestSegmentSet({len(self)} segments, "
                f"{len(self.active)} active)")


def fitted_norms(path="output/chiSquared.out", nsegments=None):
    """Per-segment normalizations from a fit, indexed by 0-based segment.

    ``chiSquared.out`` is one comma-separated line per segment --
    ``key, chi2, N, norm,`` -- followed by a totals line that is *not* CSV, so
    it has to be parsed by segment key rather than by row position.  (Reading
    it with ``skiprows=2`` drops segment 1 and shifts every normalization onto
    the wrong segment.)

    AZURE2's residual is ``(fit - data*n)/(err*n)``: the normalization scales
    the *data*, and its own ``AZUREOut_*.out`` files report the scaled values.
    Multiply measured points by ``fitted_norms()[i]`` to put them on the same
    footing as the model curve.
    """
    norms = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 4:
                continue                       # header, blank, or totals line
            try:
                key = int(parts[0])
                norms[key] = float(parts[3])
            except ValueError:
                continue
    if not norms:
        raise ValueError(f"{path}: no per-segment lines found.")
    n = nsegments if nsegments is not None else max(norms)
    missing = [k for k in range(1, n + 1) if k not in norms]
    if missing:
        raise ValueError(f"{path} has no normalization for segment(s) "
                         f"{missing[:5]}{'...' if len(missing) > 5 else ''}.")
    return np.array([norms[k] for k in range(1, n + 1)])
