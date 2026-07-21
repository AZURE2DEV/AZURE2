"""Dataset provenance for an AZURE2 model, read from the ``.azr`` file.

The ``<segmentsData>`` block records, per data segment, where the data came from
(the data file), which reaction channel it measures (entrance / exit particle
pairs), its energy and angle range, the observable type, and the normalization
systematic error.  AZURE2's API exposes the *values* (energies, cross sections)
but not this provenance, so -- exactly as the fitting examples already read the
per-segment systematic error straight from the file -- this module parses the
block in Python.  No C++/API change is needed.

Layout of a ``<segmentsData>`` line (see ``include/SegLine.h``)::

    isActive entranceKey exitKey minE maxE minA maxA isDiff
        [phaseJ phaseL if isDiff%10==2]
        dataNorm varyNorm dataNormError
        [energyShift energyShiftError varyEnergyShift]      (optional)
        dataFile  [advanced flags ...]
"""

import os
from dataclasses import dataclass
from typing import List, Optional


def _isfloat(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


_OBSERVABLE = {0: "angle-integrated", 1: "differential", 2: "phase-shift"}


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
    is_thm: bool                   # Trojan-Horse (HOES) modified segment?
    norm: float                    # applied normalization
    vary_norm: bool                # is the normalization a fit parameter?
    norm_error: float              # normalization systematic, fractional
    data_file: str                 # path the data came from  <-- provenance

    @property
    def name(self) -> str:
        """A short label: the data file's base name without extension."""
        return os.path.splitext(os.path.basename(self.data_file))[0]

    def reaction(self, pairs=None) -> str:
        """Human-readable reaction, e.g. ``'pair1 -> pair2'`` (or with pair
        spins if a :class:`~pyazr.PairSet` is supplied)."""
        ex = "total" if self.exit_key == -1 else f"pair{self.exit_key}"
        return f"pair{self.entrance_key} -> {ex}"

    def __repr__(self):
        return (f"Segment(#{self.key} {self.name!r}, "
                f"{self.entrance_key}->{self.exit_key}, {self.observable}, "
                f"E={self.energy_min:g}-{self.energy_max:g}, "
                f"norm_err={self.norm_error:g})")


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
        if isDiff % 10 == 2:            # phase-shift carries J, L
            i += 2
        norm = float(t[i]); vary_norm = int(float(t[i + 1])) == 1
        norm_error = float(t[i + 2]); i += 3
        # optional energy-shift triple, then the (non-numeric) data file
        while i < len(t) and _isfloat(t[i]):
            i += 1
        data_file = t[i] if i < len(t) else ""
        return Segment(
            key=key, active=active, entrance_key=entrance_key,
            exit_key=exit_key, energy_min=eMin, energy_max=eMax,
            angle_min=aMin, angle_max=aMax,
            observable=_OBSERVABLE.get(isDiff % 10, f"code{isDiff % 10}"),
            is_thm=isDiff >= 10, norm=norm, vary_norm=vary_norm,
            norm_error=norm_error, data_file=data_file)

    # -- views ----------------------------------------------------------------

    @property
    def active(self):
        return SegmentSet(s for s in self if s.active)

    def by_reaction(self, entrance_key=None, exit_key=None):
        return SegmentSet(
            s for s in self
            if (entrance_key is None or s.entrance_key == entrance_key)
            and (exit_key is None or s.exit_key == exit_key))

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

    def files(self):
        return [s.data_file for s in self]

    def table(self, pairs=None):
        rows = [("#", "data file", "reaction", "observable", "E range",
                 "norm_err%", "vary")]
        for s in self:
            rows.append((
                str(s.key), s.name, s.reaction(pairs),
                s.observable + (" THM" if s.is_thm else ""),
                f"{s.energy_min:g}-{s.energy_max:g}",
                f"{s.norm_error:g}", "*" if s.vary_norm else ""))
        w = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        return "\n".join("  ".join(cell.ljust(w[c]) for c, cell in enumerate(r))
                         for r in rows)

    def __repr__(self):
        return f"SegmentSet({len(self)} segments, {len(self.active)} active)"
