"""Editable, file-backed model of an AZURE2 ``.azr`` level scheme.

``AZURE2`` builds its R-matrix model from the ``<levels>`` block of a ``.azr``
file: one whitespace-separated line per *channel*, grouped into levels.  This
module parses that block into structured, mutable objects so a caller can add or
remove levels and channels from Python and write the result to a new file --
**without touching the original**.  Everything outside ``<levels>`` (config,
potential, data segments, target integration, ...) is preserved verbatim.

Because AZURE2 reads its model from a file, an edited scheme is *applied* by
writing a file and launching a fresh instance from it::

    from pyazr import AzrModel, azure2

    model = AzrModel.from_file("13N.azr")
    print(model)                                  # inspect the scheme
    model.remove_level(jpi="1/2+", energy=20)     # drop a background pole
    model.add_level(J=1.5, parity=+1, energy=4.1,
                    channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0),
                              dict(pair=2, L=1, S=0.5, gamma=0.1)])
    path = model.write("13N_edited.azr")          # original file untouched
    azr = azure2(path)                            # run with the edited scheme

The 31 fields of a channel line match ``NucLine`` in the AZURE2 source
(``include/NucLine.h``); note the file stores ``2*S`` and ``2*L`` as integers,
which the ``S`` / ``L`` properties convert.
"""

import os
import re
import tempfile
from typing import List, Optional

# Field order of a <levels> channel line, matching NucLine's read order.
_FIELDS = [
    "levelJ", "levelPi", "levelE", "levelFix", "aa", "ir", "s2", "l2",
    "levelID", "isActive", "channelFix", "gamma", "j1", "pi1", "j2", "pi2",
    "e2", "m1", "m2", "z1", "z2", "entranceSepE", "sepE", "j3", "pi3", "e3",
    "pType", "chRad", "g1", "g2", "ecMultMask",
]
_IDX = {name: i for i, name in enumerate(_FIELDS)}
_NFIELDS = len(_FIELDS)


def _token_spans(line):
    """(start, end) of every whitespace-separated token in ``line``."""
    return [m.span() for m in re.finditer(r"\S+", line)]


def _isnum(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _fmt(x):
    """Round-trippable, compact token for a number."""
    if isinstance(x, int):
        return str(x)
    if float(x).is_integer() and abs(x) < 1e15:
        return str(int(x))
    return repr(float(x))


class AzrChannel:
    """One channel line (31 fields).  Level fields are shared by a level's lines.

    Keeps the line exactly as it was read, so an untouched channel round-trips
    byte for byte; typed accessors parse on demand and mutators rewrite only the
    field they change, in place, keeping the column the file had it in.
    """

    def __init__(self, tokens, raw=None):
        if len(tokens) != _NFIELDS:
            raise ValueError(
                f"a <levels> line needs {_NFIELDS} fields, got {len(tokens)}: "
                f"{tokens}")
        self.tokens = list(tokens)
        #: The line as read, or None for a channel built from tokens alone.
        self._raw = raw
        #: Field names whose token no longer matches ``_raw``.
        self._dirty = set()

    # -- typed field access ---------------------------------------------------

    def _get(self, name, cast):
        return cast(self.tokens[_IDX[name]])

    def _set(self, name, value):
        formatted = _fmt(value)
        if self.tokens[_IDX[name]] == formatted:
            return                      # unchanged: do not disturb the raw line
        self.tokens[_IDX[name]] = formatted
        self._dirty.add(name)

    # channel identity

    @property
    def pair(self):
        """1-based particle-pair number this channel decays to."""
        return self._get("ir", lambda v: int(float(v)))

    @property
    def entrance_key(self):
        """Entrance-pair key stored on the line."""
        return self._get("aa", lambda v: int(float(v)))

    @property
    def L(self):
        """Orbital angular momentum of the channel."""
        return self._get("l2", lambda v: int(float(v))) // 2

    @L.setter
    def L(self, v):        self._set("l2", int(2 * v))

    @property
    def S(self):
        """Channel spin."""
        return self._get("s2", lambda v: int(float(v))) / 2.0

    @S.setter
    def S(self, v):        self._set("s2", int(round(2 * v)))

    @property
    def gamma(self):
        """The <levels> width field: a partial width in eV for an open channel, an ANC in fm^-1/2 for a closed one -- not a reduced-width amplitude."""
        return self._get("gamma", float)

    @gamma.setter
    def gamma(self, v):    self._set("gamma", float(v))

    @property
    def channel_fixed(self):
        """Is this channel's width held fixed in the fit?"""
        return self._get("channelFix", lambda v: int(float(v))) != 0

    @channel_fixed.setter
    def channel_fixed(self, v): self._set("channelFix", 1 if v else 0)

    @property
    def ptype(self):
        """Particle type code: 0 for a particle channel, nonzero for a photon."""
        return self._get("pType", lambda v: int(float(v)))

    @property
    def is_photon(self):
        """Is this a photon channel?"""
        return self.ptype != 0

    @property
    def channel_radius(self):
        """Channel radius in fm."""
        return self._get("chRad", float)

    @channel_radius.setter
    def channel_radius(self, v): self._set("chRad", float(v))

    @property
    def active(self):
        """Is the level active?"""
        return self._get("isActive", lambda v: int(float(v))) != 0

    # pair physics -- the same quantities GET_PAIRS_INFO reports at runtime, but
    # readable straight from the file, so a caller can identify a model's
    # channels without launching AZURE2.

    @property
    def Z1(self):
        """Charge number of the light particle."""
        return self._get("z1", lambda v: int(float(v)))

    @property
    def Z2(self):
        """Charge number of the heavy particle."""
        return self._get("z2", lambda v: int(float(v)))

    @property
    def M1(self):
        """Mass of the light particle, in u."""
        return self._get("m1", float)

    @property
    def M2(self):
        """Mass of the heavy particle, in u."""
        return self._get("m2", float)

    @property
    def J1(self):
        """Intrinsic spin of the light particle."""
        return self._get("j1", float)

    @property
    def parity1(self):
        """Parity of the light particle."""
        return self._get("pi1", lambda v: int(float(v)))

    @property
    def J2(self):
        """Intrinsic spin of the heavy particle."""
        return self._get("j2", float)

    @property
    def parity2(self):
        """Parity of the heavy particle."""
        return self._get("pi2", lambda v: int(float(v)))

    @property
    def excitation(self):
        """Excitation energy of the pair's residual nucleus (MeV).

        Zero for a ground-state pair; this is what distinguishes the capture
        channels of a multi-transition model (gamma_0, gamma_1, ...).
        """
        return self._get("e2", float)

    @property
    def sep_energy(self):
        """Separation energy of the pair, in MeV."""
        return self._get("sepE", float)

    # level fields (shared across a level's channel lines)

    @property
    def levelJ(self):
        """Total angular momentum of the level."""
        return self._get("levelJ", float)

    @property
    def levelPi(self):
        """Parity of the level."""
        return self._get("levelPi", lambda v: int(float(v)))

    @property
    def levelE(self):
        """Level energy in MeV -- an excitation energy of the compound nucleus."""
        return self._get("levelE", float)

    @property
    def levelID(self):
        """The level's number, shared by all of its channel lines."""
        return self._get("levelID", lambda v: int(float(v)))

    @property
    def level_fixed(self):
        """Is the level energy held fixed in the fit?"""
        return self._get("levelFix", lambda v: int(float(v))) != 0

    def _set_level(self, J, parity, energy, level_fixed, level_id):
        self._set("levelJ", float(J))
        self._set("levelPi", int(parity))
        self._set("levelE", float(energy))
        self._set("levelFix", 1 if level_fixed else 0)
        self._set("levelID", int(level_id))

    def clone(self):
        """A copy of this channel line, sharing nothing with the original.

        The raw line comes with it, so a cloned channel whose couplings are then
        changed keeps the column layout of the one it was cloned from.
        """
        copy = AzrChannel(self.tokens, self._raw)
        copy._dirty = set(self._dirty)
        return copy

    def to_line(self):
        """The channel as one line of the ``<levels>`` block.

        An untouched line is returned exactly as it was read.  An edited one is
        the same line with the changed fields substituted in place, right
        aligned in the width the file used, so a single edit does not reflow the
        whole block.  A channel with no raw line behind it -- built from tokens,
        not parsed -- falls back to single spaces.
        """
        if self._raw is None:
            return " ".join(self.tokens)
        if not self._dirty:
            return self._raw
        spans = _token_spans(self._raw)
        if len(spans) != _NFIELDS:      # not the line we parsed; do not guess
            return " ".join(self.tokens)
        out = self._raw
        for name in sorted(self._dirty, key=lambda n: -_IDX[n]):
            start, end = spans[_IDX[name]]
            new = self.tokens[_IDX[name]]
            width = end - start
            out = out[:start] + (new.rjust(width) if len(new) <= width else new) + out[end:]
        return out

    @classmethod
    def from_line(cls, line):
        """Parse one .azr channel line, keeping it verbatim for round-tripping."""
        return cls(line.split(), raw=line)


class AzrLevel:
    """A level: shared (J, parity, energy, fixed) + its channels."""

    def __init__(self, channels: List[AzrChannel]):
        if not channels:
            raise ValueError("a level needs at least one channel.")
        self.channels = list(channels)

    @property
    def J(self):
        """Total angular momentum of the level."""
        return self.channels[0].levelJ

    @property
    def parity(self):
        """Parity of the level."""
        return self.channels[0].levelPi

    @property
    def energy(self):
        """Level energy in MeV."""
        return self.channels[0].levelE

    @property
    def level_id(self):
        """The level's number."""
        return self.channels[0].levelID

    @property
    def fixed(self):
        """Is the level energy held fixed?"""
        return self.channels[0].level_fixed

    @property
    def jpi(self):
        """J^pi as text, e.g. "3/2-"."""
        j = int(self.J) if float(self.J).is_integer() else f"{int(round(2*self.J))}/2"
        return f"{j}{'+' if self.parity > 0 else '-'}"

    def set_energy(self, energy):
        """Set the level energy on every channel line of the level."""
        for c in self.channels:
            c._set("levelE", float(energy))

    def set_fixed(self, fixed):
        """Fix or free the level energy on every channel line."""
        for c in self.channels:
            c._set("levelFix", 1 if fixed else 0)

    def _renumber(self, level_id):
        for c in self.channels:
            c._set("levelID", int(level_id))

    def __repr__(self):
        return (f"AzrLevel(J^pi={self.jpi}, E={self.energy:g} MeV, "
                f"{len(self.channels)} channels"
                f"{', fixed' if self.fixed else ''})")


class AzrModel:
    """A ``.azr`` file parsed into an editable level scheme.

    Only the ``<levels>`` block is interpreted; the surrounding text is kept
    verbatim and re-emitted unchanged by :meth:`write`.
    """

    def __init__(self, prefix, levels, suffix, source=None):
        self._prefix = prefix       # text up to and including "<levels>\n"
        self.levels: List[AzrLevel] = levels
        self._suffix = suffix       # text from "</levels>" onward
        self.source = source

    # -- parsing / serialization ---------------------------------------------

    @classmethod
    def from_file(cls, path):
        """Parse a .azr file; only <levels> is interpreted, the rest is kept verbatim."""
        with open(path) as f:
            text = f.read()
        start = text.find("<levels>")
        end = text.find("</levels>")
        if start == -1 or end == -1:
            raise ValueError(f"{path}: no <levels> ... </levels> block.")
        body_start = text.index("\n", start) + 1
        prefix = text[:body_start]
        body = text[body_start:end]
        suffix = text[end:]

        # group channel lines into levels by their levelID
        levels, cur, cur_id = [], [], None
        for line in body.splitlines():
            if not line.strip():
                continue
            ch = AzrChannel.from_line(line)
            if cur and ch.levelID != cur_id:
                levels.append(AzrLevel(cur))
                cur = []
            cur.append(ch)
            cur_id = ch.levelID
        if cur:
            levels.append(AzrLevel(cur))
        return cls(prefix, levels, suffix, source=path)

    def set_output_dir(self, path):
        """Point the model's ``<config>`` output directory somewhere else.

        AZURE2 writes its results *and* its external-capture integral caches
        (``intEC.dat`` / ``intEC.extrap``) here.  Giving an edited model its own
        directory keeps it from overwriting the main run's outputs -- and, since
        the cache is only valid for the grid it was built on, keeps a variant's
        integrals from being read back into a run with different segments.

        The directory is created if it does not exist; AZURE2 does not create it
        and fails quietly otherwise.
        """
        path = str(path)
        if not path.endswith("/"):
            path += "/"
        os.makedirs(path, exist_ok=True)
        lines = self._prefix.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if line.rstrip().endswith("#Full Path to Output Directory"):
                comment = line[line.index("#"):]
                lines[i] = f"{path:<100}{comment}"
                break
        else:
            raise ValueError("no output-directory line in the <config> block.")
        self._prefix = "".join(lines)
        return self

    def to_text(self):
        """The whole file as text, with <levels> re-emitted and everything else unchanged."""
        self._renumber()
        blocks = ["\n".join(c.to_line() for c in lv.channels)
                  for lv in self.levels]
        body = "\n\n".join(blocks)
        return self._prefix + body + "\n\n" + self._suffix

    def write(self, path):
        """Write the (edited) model to ``path`` and return it.  The original
        file is never modified."""
        with open(path, "w") as f:
            f.write(self.to_text())
        return path

    def to_tempfile(self, suffix=".azr", dir=None):
        """Write to a fresh temp file (for launching an edited instance) and
        return its path.  The caller owns the file."""
        fd, path = tempfile.mkstemp(suffix=suffix, dir=dir)
        os.close(fd)
        return self.write(path)

    # -- queries --------------------------------------------------------------

    def find(self, jpi=None, energy=None, tol=1e-3, index=None):
        """Return the levels matching ``jpi`` and/or ``energy`` (or by index)."""
        if index is not None:
            return [self.levels[index]]
        out = []
        for lv in self.levels:
            if jpi is not None and lv.jpi != jpi:
                continue
            if energy is not None and abs(lv.energy - energy) > tol:
                continue
            out.append(lv)
        return out

    def _pair_template(self, pair):
        """An existing channel that uses ``pair`` (for cloning its pair data)."""
        for lv in self.levels:
            for c in lv.channels:
                if c.pair == pair:
                    return c
        raise ValueError(
            f"pair {pair} is not used by any existing channel; adding a brand-"
            f"new particle pair is not supported -- reference an existing pair.")

    def _renumber(self):
        """Give every level a fresh consecutive levelID (1-based)."""
        for i, lv in enumerate(self.levels, start=1):
            lv._renumber(i)

    # -- editing --------------------------------------------------------------

    def remove_level(self, jpi=None, energy=None, index=None, tol=1e-3):
        """Remove the matching level(s).  Returns the removed :class:`AzrLevel`
        objects.  Raises if the selector matches nothing."""
        victims = self.find(jpi=jpi, energy=energy, index=index, tol=tol)
        if not victims:
            raise KeyError(f"no level matches jpi={jpi} energy={energy} "
                           f"index={index}.")
        self.levels = [lv for lv in self.levels if lv not in victims]
        self._renumber()
        return victims

    def add_level(self, J, parity, energy, channels, level_fixed=True,
                  at=None):
        """Add a level with the given channels.

        ``channels`` is a list of dicts, each ``{pair, L, S, gamma[, fixed]}``,
        referencing an existing particle pair by number; the pair's physical
        data (masses, charges, separation energy, radius) is copied from an
        existing channel that uses that pair.

        **All levels of a given J^pi share one channel set** (an R-matrix
        J-group requirement).  If a level of this ``(J, parity)`` already exists,
        this method clones that group's exact channel structure and only applies
        your ``gamma`` / ``fixed`` values, matched by ``(pair, L, S)`` -- any
        channel you don't mention is added inert (``gamma=0``, fixed).  A spec
        that matches none of the group's channels raises, listing the allowed
        ones.  For a brand-new J^pi the channels you give define the group.

        Returns the new :class:`AzrLevel`.
        """
        new_id = max((lv.level_id for lv in self.levels), default=0) + 1
        existing = [lv for lv in self.levels
                    if lv.J == J and lv.parity == parity]

        if existing:
            # Clone the group's channel structure; start every channel inert.
            template = existing[0]
            chans = [c.clone() for c in template.channels]
            for ch in chans:
                ch._set_level(J, parity, energy, level_fixed, new_id)
                ch._set("isActive", 1)
                ch.gamma = 0.0
                ch.channel_fixed = True
            allowed = [(ch.pair, ch.L, ch.S) for ch in chans]
            for spec in channels:
                match = [ch for ch in chans
                         if ch.pair == spec["pair"] and ch.L == spec["L"]
                         and abs(ch.S - spec["S"]) < 1e-6]
                if not match:
                    raise ValueError(
                        f"J^pi {template.jpi} already exists with a fixed "
                        f"channel set {allowed}; the spec "
                        f"(pair={spec['pair']}, L={spec['L']}, S={spec['S']}) "
                        f"matches none of them.")
                match[0].gamma = spec.get("gamma", 0.0)
                match[0].channel_fixed = spec.get("fixed", False)
        else:
            chans = []
            for spec in channels:
                ch = self._pair_template(spec["pair"]).clone()
                ch._set_level(J, parity, energy, level_fixed, new_id)
                ch._set("isActive", 1)
                ch.L = spec["L"]
                ch.S = spec["S"]
                ch.gamma = spec.get("gamma", 0.0)
                ch.channel_fixed = spec.get("fixed", False)
                chans.append(ch)

        level = AzrLevel(chans)
        if at is None:
            self.levels.append(level)
        else:
            self.levels.insert(at, level)
        self._renumber()
        return level

    def remove_channel(self, jpi, pair, L, S, tol=1e-6):
        """Remove a channel ``(pair, L, S)`` from *every* level of a J^pi group.

        All levels of a J-group share one channel set, so a channel is removed
        from the whole group at once (removing it from a single level would make
        AZURE2 reject the file).  Returns the number of channel lines removed.
        Raises if no level of ``jpi`` has that channel, or if it is a level's
        only channel.
        """
        removed = 0
        for lv in self.levels:
            if lv.jpi != jpi:
                continue
            keep = [c for c in lv.channels
                    if not (c.pair == pair and c.L == L
                            and abs(c.S - S) < tol)]
            if len(keep) == len(lv.channels):
                continue
            if not keep:
                raise ValueError(
                    f"channel (pair={pair}, L={L}, S={S}) is the only channel "
                    f"of a {jpi} level; cannot remove it.")
            removed += len(lv.channels) - len(keep)
            lv.channels = keep
        if removed == 0:
            raise KeyError(f"no {jpi} level has channel "
                           f"(pair={pair}, L={L}, S={S}).")
        return removed

    def add_channel(self, level, pair, L, S, gamma=0.0, fixed=False):
        """Add one channel to an existing level (referencing an existing pair)."""
        template = self._pair_template(pair)
        ch = template.clone()
        ch._set_level(level.J, level.parity, level.energy, level.fixed,
                      level.level_id)
        ch._set("isActive", 1)
        ch.L, ch.S, ch.gamma = L, S, gamma
        ch.channel_fixed = fixed
        level.channels.append(ch)
        return ch

    # -- channel radius -------------------------------------------------------

    def channel_radii(self):
        """``{pair: radius}`` in fm, for every particle pair in the file."""
        out = {}
        for lv in self.levels:
            for c in lv.channels:
                if not c.is_photon:
                    out[c.pair] = c.channel_radius
        return dict(sorted(out.items()))

    def set_channel_radius(self, pair, radius):
        """Set the channel radius (fm) of one particle pair, on every line.

        AZURE2 stores the radius per channel line, so it must be written to all
        of them or the model is inconsistent; this does that and returns the
        number of lines changed.

        The radius is the matching surface between the internal R-matrix region
        and the external Coulomb solutions, so changing it changes the
        penetrabilities, shift functions, boundary conditions, Wigner limits and
        the lower limit of every external-capture integral. A reduced width
        therefore *means* something different afterwards, and the level scheme
        is no longer fitted -- **refit before reading anything off the result**.
        Delete ``output/intEC.dat`` / ``output/intEC.extrap`` (or write the new
        model into its own output directory) so the cached external-capture
        integrals, which belong to the old radius, cannot be reused.

        >>> mdl = AzrModel.from_file("7Be.azr")
        >>> mdl.channel_radii()                 # {1: 4.24151, 2: 3.94396, 3: 3.94396}
        >>> mdl.set_channel_radius(1, 5.0)      # 3He+alpha
        >>> path = mdl.write("a5.0.azr")
        """
        radius = float(radius)
        if not radius > 0:
            raise ValueError(f"channel radius must be positive, got {radius}.")
        pair = int(pair)
        changed = 0
        for lv in self.levels:
            for c in lv.channels:
                if c.pair == pair and not c.is_photon:
                    c.channel_radius = radius
                    changed += 1
        if not changed:
            raise KeyError(
                f"pair {pair} has no particle channel in this model "
                f"(radii: {self.channel_radii()}).")
        return changed

    # -- level activation -----------------------------------------------------

    def deactivate_level(self, jpi=None, energy=None, index=None, tol=1e-3):
        """Switch off the matching level(s) without removing them.

        Every channel of a matched level has its reduced width (``gamma``) set
        to zero and is fixed, so the level stays in the file (and in its
        J-group's shared channel set) but couples to nothing -- it contributes
        exactly nothing to the calculation.  This is the file-level counterpart
        of :meth:`azure2.without_level`; use it when you want to persist a
        "level removed" model and refit it.  Returns the affected levels.
        """
        victims = self.find(jpi=jpi, energy=energy, index=index, tol=tol)
        if not victims:
            raise KeyError(f"no level matches jpi={jpi} energy={energy} "
                           f"index={index}.")
        for lv in victims:
            for c in lv.channels:
                c.gamma = 0.0
                c.channel_fixed = True
        return victims

    # -- extrapolations (edits the <segmentsTest> block) ----------------------

    # observable name -> isDiff code for a <segmentsTest> line
    # (ESegment::ESegment(ExtrapLine); mirrors datasets._EXTRAP_OBSERVABLE).
    _EXTRAP_CODE = {
        "angle-integrated": 0, "differential": 1, "phase-shift": 2,
        "angular-distribution": 3, "total-capture": 4, "differential-cm": 5,
        # Vector analyzing power A_y. Differential in the centre-of-mass frame
        # like code 5, but the quantity is a dimensionless ratio, not a cross
        # section.
        "analyzing-power": 7,
    }

    def _splice_segments_test(self, new_lines):
        """Replace the body of the <segmentsTest> block in the verbatim tail.

        ``new_lines`` is a list of pre-formatted extrapolation lines (no
        trailing newline).  Creates the block if the file has none.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsTest>")
            end = lines.index("</segmentsTest>")
        except ValueError:
            # No block yet: append one at the end of the tail.
            self._suffix = (self._suffix.rstrip("\n") + "\n\n<segmentsTest>\n"
                            + "\n".join(new_lines) + "\n</segmentsTest>\n")
            return
        self._suffix = "\n".join(lines[:start + 1] + new_lines
                                 + lines[end:])

    def clear_extrapolations(self):
        """Remove every ``<segmentsTest>`` line (leave the block empty)."""
        self._splice_segments_test([])
        return self

    def add_extrapolation(self, entrance, exit, e_min, e_max, e_step,
                          observable="angle-integrated", angle=None,
                          angle_min=0.0, angle_max=0.0, angle_step=0.0,
                          phase_J=None, phase_L=None, order=None, active=True):
        """Append one extrapolation segment (a grid to evaluate the model on).

        ``entrance`` / ``exit`` are particle-pair keys (``exit=-1`` for a summed
        / total observable such as capture).  The model is evaluated on the
        energy grid ``e_min : e_max : e_step`` (entrance-channel c.m. MeV).

        ``observable`` is one of ``angle-integrated``, ``differential``,
        ``differential-cm``, ``total-capture``, ``angular-distribution``,
        ``phase-shift``.  For a single-angle differential give ``angle=`` (sets
        min=max, step 0); for an angular grid give ``angle_min/max/step``.
        ``phase_J``/``phase_L`` are required for a phase-shift, ``order`` for an
        angular distribution.
        """
        if observable not in self._EXTRAP_CODE:
            raise ValueError(f"unknown observable {observable!r}; expected one "
                             f"of {sorted(self._EXTRAP_CODE)}.")
        isDiff = self._EXTRAP_CODE[observable]
        if angle is not None:
            angle_min = angle_max = angle
            angle_step = 0.0
        toks = [1 if active else 0, int(entrance), int(exit),
                _fmt(e_min), _fmt(e_max), _fmt(e_step),
                _fmt(angle_min), _fmt(angle_max), _fmt(angle_step), isDiff]
        if isDiff == 2:                       # phase shift carries J, L
            if phase_J is None or phase_L is None:
                raise ValueError("a phase-shift extrapolation needs phase_J "
                                 "and phase_L.")
            toks += [_fmt(phase_J), int(phase_L)]
        elif isDiff == 3:                     # angular distribution carries order
            if order is None:
                raise ValueError("an angular-distribution extrapolation needs "
                                 "order.")
            toks += [int(order)]
        toks += [0]                           # trailing isAdvanced flag
        line = "  ".join(_fmt(t) if not isinstance(t, str) else t for t in toks)
        lines = self._suffix.splitlines()
        if "<segmentsTest>" in lines:
            end = lines.index("</segmentsTest>")
            self._suffix = "\n".join(lines[:end] + [line] + lines[end:])
        else:
            self._splice_segments_test([line])
        return self

    def set_extrapolations(self, specs):
        """Replace the whole ``<segmentsTest>`` block with ``specs``.

        ``specs`` is an iterable of dicts, each a keyword bundle for
        :meth:`add_extrapolation` (e.g. ``dict(entrance=1, exit=-1, e_min=0.05,
        e_max=5, e_step=0.05, observable="total-capture")``).
        """
        self.clear_extrapolations()
        for spec in specs:
            self.add_extrapolation(**spec)
        return self

    def keep_extrapolations(self, keys):
        """Drop every ``<segmentsTest>`` line except the given 1-based ``keys``,
        in the order they are listed here.

        The surviving lines are copied verbatim, so the grids are exactly the
        ones the file declares.  AZURE2 re-evaluates *every* active segment on
        each forward pass, so trimming the block to the handful you actually
        want is what makes a finite-difference uncertainty band affordable --
        see :func:`pyazr.bands.extrapolation_bands`.

        Note this invalidates ``output/intEC.extrap``: AZURE2 caches the
        external-capture integrals there and reuses them on a changed grid,
        silently corrupting the result.  Delete it before running the edited
        model.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsTest>") + 1
            end = lines.index("</segmentsTest>")
        except ValueError:
            raise ValueError("no <segmentsTest> block to trim.")
        body = [ln for ln in lines[start:end] if ln.strip()]
        missing = [k for k in keys if not 1 <= k <= len(body)]
        if missing:
            raise KeyError(f"<segmentsTest> has {len(body)} segments; no "
                           f"{missing}.")
        self._splice_segments_test([body[k - 1] for k in keys])
        return self

    # -- segment normalizations (edits the <segmentsData> block) --------------

    def set_segment_norm(self, file_substr, vary=None, sys_error=None):
        """Edit the normalization of every ``<segmentsData>`` line whose text
        matches ``file_substr`` (e.g. a data-file name).

        ``vary=True/False`` frees/fixes the normalization; ``sys_error`` sets its
        systematic (percent, as stored in the file).  Returns the number of
        segment lines changed.  Operates on the verbatim tail text, so the
        levels editing is unaffected.
        """
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                isDiff = int(float(t[7]))
                i = 8 + (2 if isDiff % 10 == 2 else 0)   # dataNorm index
                if vary is not None:
                    t[i + 1] = "1" if vary else "0"
                if sys_error is not None:
                    t[i + 2] = _fmt(sys_error)
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._suffix = "\n".join(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    def engine_level_keys(self):
        """``{(jgroup, level): AzrLevel}`` -- the numbering AZURE2 itself uses.

        The engine builds a J-group the first time it meets a ``(J, parity)``
        while reading ``<levels>``, and numbers levels within a group in file
        order.  Reproducing that here gives the same ``(jgroup, level)`` a
        :class:`~pyazr.parameters.Parameter` reports, which is the only reliable
        way to say *which* level a fitted value belongs to: a level at
        ``Ex = 0`` comes back from the API with ``level_energy = None``, so an
        energy-based key cannot identify it.

        Both indices are 1-based, as the API reports them.
        """
        order, seen, out = [], {}, {}
        for lv in self.levels:
            k = (int(round(2 * lv.J)), int(lv.parity))
            if k not in seen:
                order.append(k)
                seen[k] = 0
            seen[k] += 1
            out[(order.index(k) + 1, seen[k])] = lv
        return out

    def apply_fit(self, parameters, x, transform=None, physical=False,
                  pairs=None, strict=True):
        """Write a fitted parameter vector into the levels block.

        **The ``gamma`` field of a ``<levels>`` line is not a reduced-width
        amplitude.**  It holds the same *physical* value AZURE2 prints in
        ``parameters.out`` and shows in the GUI: a partial width in eV for an
        open particle channel, an ANC in fm^-1/2 for a closed (sub-threshold)
        one, and a partial width in eV for a photon channel.  Writing an rwa
        there produces a file that loads without complaint and is wrong -- for
        the 7Be model the two differ by factors of 10^2 to 10^7.

        So ``x`` (a free vector in ``params_rwa`` order) must be converted
        first.  Either hand over the transform and let this do it:

        >>> mdl.apply_fit(m.parameters, x_best, transform=m.transform_rwa)

        or convert yourself and say so:

        >>> mdl.apply_fit(m.parameters, m.transform_rwa(x_best), physical=True)

        Passing neither raises, rather than silently writing an rwa.

        Levels are matched on the ``(jgroup, level)`` the engine reports, via
        :meth:`engine_level_keys`, and channels on ``(pair, L, S)`` within the
        matched level.

        **Pass ``pairs=m.pairs``.**  A ``Parameter``'s ``pair`` is the engine's
        number, which counts particle pairs in the order ``<levels>`` first
        mentions them -- not the pair key the file writes.  The two differ
        whenever the levels do not introduce the pairs in key order: on the 8Be
        model engine pair 1 is file key 2, and file key 1 is engine pair 6, so
        matching one against the other writes every width to the wrong channel.
        The :class:`~pyazr.parameters.PairSet` carries both, and is the only
        thing that can translate.

        ``strict`` (the default) raises if any free parameter finds no home,
        rather than skipping it: a partial write produces a file that loads
        cleanly and is a mixture of two fits.

        Note this covers ``<levels>`` only -- normalizations and energy shifts
        are not in that block, so a fit that moved them is only half saved.
        :meth:`pyazr.azure2.azure2.save_fit` writes the companion
        ``param.sav`` and verifies the result; prefer it to calling this
        directly.  The number of values written is left in :attr:`applied`.
        """
        if transform is None and not physical:
            raise ValueError(
                "apply_fit needs physical values, not reduced-width amplitudes: "
                "pass transform=m.transform_rwa, or convert with "
                "m.transform_rwa(x) yourself and pass physical=True.")
        if transform is not None:
            if physical:
                raise ValueError("pass transform= or physical=True, not both.")
            x = transform(x)

        lvlmap = self.engine_level_keys()
        # engine pair number -> the pair key the file writes
        pairkey = {p.number: p.key for p in pairs} if pairs is not None else {}
        written, unplaced = 0, []
        for p in parameters:
            if p.fixed or p.free_index is None:
                continue
            if p.kind not in ("energy", "width"):
                continue            # norms and shifts do not live in <levels>
            if p.free_index >= len(x):
                unplaced.append(f"{p.name} (free_index {p.free_index} beyond the vector)")
                continue
            lv = lvlmap.get((p.jgroup, p.level))
            if lv is None:
                unplaced.append(f"{p.name} (no level at jgroup {p.jgroup}, level {p.level})")
                continue
            v = float(x[p.free_index])
            if p.kind == "energy":
                lv.set_energy(v)
                written += 1
                continue
            want_pair = pairkey.get(p.pair, p.pair)
            for c in lv.channels:
                if (c.pair == want_pair and c.L == p.L
                        and abs(c.S - (p.S or 0.0)) < 1e-6):
                    c.gamma = v
                    written += 1
                    break
            else:
                unplaced.append(
                    f"{p.name} (no channel pair={want_pair} L={p.L} S={p.S} "
                    f"in {lv.jpi} at {lv.energy} MeV)")

        if unplaced and strict:
            raise ValueError(
                f"apply_fit could not place {len(unplaced)} of the free "
                f"parameters, so the result would be a mixture of two fits:\n  "
                + "\n  ".join(unplaced[:8])
                + (f"\n  ... and {len(unplaced) - 8} more" if len(unplaced) > 8 else "")
                + "\nThe .azr and the parameter set do not describe the same model. "
                  "Pass strict=False to write the rest anyway."
                + ("\nNote apply_fit was given no pairs=, so it assumed the "
                   "engine's pair numbers are the file's pair keys. Pass "
                   "pairs=m.pairs if they are not." if pairs is None else ""))
        self.applied = written
        self.unplaced = unplaced
        return self

    def set_segment_datafile(self, file_substr, new_path):
        """Repoint matching ``<segmentsData>`` lines to ``new_path`` (the data
        file is the first non-numeric token).  Returns the number changed."""
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                for i, tok in enumerate(t):
                    if not _isnum(tok):
                        t[i] = new_path
                        break
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._suffix = "\n".join(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    def set_segment_active(self, file_substr, active):
        """Activate/deactivate every ``<segmentsData>`` line matching
        ``file_substr`` (sets the leading isActive field).  Deactivated segments
        are ignored by AZURE2.  Returns the number of lines changed."""
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                t[0] = "1" if active else "0"
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._suffix = "\n".join(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    # -- adding / removing whole data segments --------------------------------

    # observable name -> isDiff code for a <segmentsData> line
    # (ESegment::ESegment(SegLine); mirrors datasets._OBSERVABLE).
    _DATA_CODE = {
        "angle-integrated": 0, "differential": 1, "phase-shift": 2,
        "total-capture": 3, "differential-cm": 4, "angle-integrated-E1": 5,
        "angle-integrated-E2": 6, "analyzing-power": 7,
    }

    def add_data_segment(self, data_file, entrance, exit,
                         observable="angle-integrated",
                         energy_min=0.0, energy_max=5.0,
                         angle_min=0.0, angle_max=180.0,
                         norm=1.0, vary_norm=False, norm_error=0.0,
                         energy_shift=0.0, energy_shift_error=0.0,
                         vary_shift=False, phase_J=None, phase_L=None,
                         active=True):
        """Append one data segment (a ``<segmentsData>`` line).

        ``data_file`` is the data file path (relative to the run directory,
        which for a model living in its own folder is usually ``data/...``).
        ``entrance`` / ``exit`` are particle-pair keys (``exit=-1`` for a
        summed/total observable).  ``observable`` is one of
        ``angle-integrated``, ``differential``, ``differential-cm``,
        ``total-capture``, ``phase-shift``, ``angle-integrated-E1``,
        ``angle-integrated-E2``.

        ``norm`` is the normalization applied to the data, ``norm_error`` its
        systematic error (percent, as stored in the file), ``energy_shift``
        the beam-energy shift (MeV) with its ``_error``.  ``vary_norm`` /
        ``vary_shift`` free the corresponding parameter.

        Returns ``self`` so calls chain.
        """
        if observable not in self._DATA_CODE:
            raise ValueError(f"unknown observable {observable!r}; expected one "
                             f"of {sorted(self._DATA_CODE)}.")
        isDiff = self._DATA_CODE[observable]
        toks = [1 if active else 0, int(entrance), int(exit),
                _fmt(energy_min), _fmt(energy_max),
                _fmt(angle_min), _fmt(angle_max), isDiff]
        if isDiff == 2:                       # phase shift carries J, L
            if phase_J is None or phase_L is None:
                raise ValueError("a phase-shift data segment needs phase_J "
                                 "and phase_L.")
            toks += [_fmt(phase_J), int(phase_L)]
        toks += [_fmt(norm), 1 if vary_norm else 0, _fmt(norm_error),
                 _fmt(energy_shift), _fmt(energy_shift_error),
                 1 if vary_shift else 0, str(data_file), 0, 0]
        line = "  ".join(t if isinstance(t, str) else _fmt(t) for t in toks)
        lines = self._suffix.splitlines()
        if "<segmentsData>" not in lines:
            raise ValueError("no <segmentsData> block to add to.")
        end = lines.index("</segmentsData>")
        self._suffix = "\n".join(lines[:end] + [line] + lines[end:])
        return self

    def remove_data_segments(self, file_substr):
        """Remove every ``<segmentsData>`` line whose text matches
        ``file_substr`` (e.g. a data-file name).  Returns the number of
        segments removed.  Raises if nothing matches.

        Removing data changes which energies AZURE2 evaluates, so the
        external-capture integrals must be recalculated before the next run --
        delete ``output/intEC.dat`` / ``output/intEC.extrap`` (or write the
        model into its own output directory) so the stale cache cannot be
        reused.
        """
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, removed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                removed += 1
                continue
            out.append(line)
        self._suffix = "\n".join(out)
        if removed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return removed

    def clear_data_segments(self):
        """Remove every ``<segmentsData>`` line (leave the block empty).

        The external-capture caches ``output/intEC.dat`` / ``output/intEC.extrap``
        belong to the removed grids and must be deleted before the next run.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsData>")
            end = lines.index("</segmentsData>")
        except ValueError:
            raise ValueError("no <segmentsData> block to clear.")
        self._suffix = "\n".join(lines[:start + 1] + lines[end:])
        return self

    # -- rendering ------------------------------------------------------------

    def __str__(self):
        head = f"AzrModel"
        if self.source:
            head += f"  [{os.path.basename(self.source)}]"
        head += f"   {len(self.levels)} levels"
        lines = [head, "=" * len(head)]
        seen = []
        for lv in self.levels:
            if lv.jpi not in seen:
                seen.append(lv.jpi)
                lines.append(f"\nJ^pi = {lv.jpi}")
            fix = "fixed" if lv.fixed else "FREE"
            lines.append(f"  E = {lv.energy:>9.4g} MeV  [{fix}]")
            for c in lv.channels:
                kind = f"{'photon' if c.is_photon else 'particle'} pair{c.pair}"
                cfix = "fixed" if c.channel_fixed else "FREE"
                lines.append(f"      {kind:<16} L={c.L} S={c.S:g}  "
                             f"Gamma={c.gamma:>11.5g}  [{cfix}]")
        return "\n".join(lines)

    def __repr__(self):
        return f"AzrModel({len(self.levels)} levels, source={self.source!r})"
