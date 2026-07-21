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

    Stores the raw string tokens so untouched channels round-trip exactly;
    typed accessors parse on demand and mutators re-format only the field they
    change.
    """

    def __init__(self, tokens):
        if len(tokens) != _NFIELDS:
            raise ValueError(
                f"a <levels> line needs {_NFIELDS} fields, got {len(tokens)}: "
                f"{tokens}")
        self.tokens = list(tokens)

    # -- typed field access ---------------------------------------------------
    def _get(self, name, cast):
        return cast(self.tokens[_IDX[name]])

    def _set(self, name, value):
        self.tokens[_IDX[name]] = _fmt(value)

    # channel identity
    @property
    def pair(self):        return self._get("ir", lambda v: int(float(v)))
    @property
    def entrance_key(self): return self._get("aa", lambda v: int(float(v)))
    @property
    def L(self):           return self._get("l2", lambda v: int(float(v))) // 2
    @L.setter
    def L(self, v):        self._set("l2", int(2 * v))
    @property
    def S(self):           return self._get("s2", lambda v: int(float(v))) / 2.0
    @S.setter
    def S(self, v):        self._set("s2", int(round(2 * v)))
    @property
    def gamma(self):       return self._get("gamma", float)
    @gamma.setter
    def gamma(self, v):    self._set("gamma", float(v))
    @property
    def channel_fixed(self): return self._get("channelFix", lambda v: int(float(v))) != 0
    @channel_fixed.setter
    def channel_fixed(self, v): self._set("channelFix", 1 if v else 0)
    @property
    def ptype(self):       return self._get("pType", lambda v: int(float(v)))
    @property
    def is_photon(self):   return self.ptype != 0
    @property
    def active(self):      return self._get("isActive", lambda v: int(float(v))) != 0

    # pair physics -- the same quantities GET_PAIRS_INFO reports at runtime, but
    # readable straight from the file, so a caller can identify a model's
    # channels without launching AZURE2.
    @property
    def Z1(self):          return self._get("z1", lambda v: int(float(v)))
    @property
    def Z2(self):          return self._get("z2", lambda v: int(float(v)))
    @property
    def M1(self):          return self._get("m1", float)
    @property
    def M2(self):          return self._get("m2", float)
    @property
    def excitation(self):
        """Excitation energy of the pair's residual nucleus (MeV).

        Zero for a ground-state pair; this is what distinguishes the capture
        channels of a multi-transition model (gamma_0, gamma_1, ...).
        """
        return self._get("e2", float)
    @property
    def sep_energy(self):  return self._get("sepE", float)

    # level fields (shared across a level's channel lines)
    @property
    def levelJ(self):      return self._get("levelJ", float)
    @property
    def levelPi(self):     return self._get("levelPi", lambda v: int(float(v)))
    @property
    def levelE(self):      return self._get("levelE", float)
    @property
    def levelID(self):     return self._get("levelID", lambda v: int(float(v)))
    @property
    def level_fixed(self): return self._get("levelFix", lambda v: int(float(v))) != 0

    def _set_level(self, J, parity, energy, level_fixed, level_id):
        self._set("levelJ", float(J))
        self._set("levelPi", int(parity))
        self._set("levelE", float(energy))
        self._set("levelFix", 1 if level_fixed else 0)
        self._set("levelID", int(level_id))

    def clone(self):
        return AzrChannel(self.tokens)

    def to_line(self):
        return " ".join(self.tokens)

    @classmethod
    def from_line(cls, line):
        return cls(line.split())


class AzrLevel:
    """A level: shared (J, parity, energy, fixed) + its channels."""

    def __init__(self, channels: List[AzrChannel]):
        if not channels:
            raise ValueError("a level needs at least one channel.")
        self.channels = list(channels)

    @property
    def J(self):        return self.channels[0].levelJ
    @property
    def parity(self):   return self.channels[0].levelPi
    @property
    def energy(self):   return self.channels[0].levelE
    @property
    def level_id(self): return self.channels[0].levelID
    @property
    def fixed(self):    return self.channels[0].level_fixed

    @property
    def jpi(self):
        j = int(self.J) if float(self.J).is_integer() else f"{int(round(2*self.J))}/2"
        return f"{j}{'+' if self.parity > 0 else '-'}"

    def set_energy(self, energy):
        for c in self.channels:
            c._set("levelE", float(energy))

    def set_fixed(self, fixed):
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

    def to_text(self):
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

    def apply_fit(self, parameters, x):
        """Write a fitted free-parameter vector into the levels block.

        ``parameters`` is a model's ``azure2.parameters`` (metadata) and ``x`` is
        the fitted free-parameter vector (``params_rwa`` order).  Each free
        level energy and reduced width is written into the matching level /
        channel (matched by 2J, parity, input energy, and channel pair/L/S), so
        the file becomes a snapshot of the fit -- the reference to reuse.
        """
        def key(J, parity, E):
            return (int(round(2 * J)), int(parity),
                    round(E, 2) if E is not None else None)

        lvlmap = {}
        for lv in self.levels:
            lvlmap.setdefault(key(lv.J, lv.parity, lv.energy), lv)

        for p in parameters:
            if p.fixed or p.free_index is None or p.J is None:
                continue
            lv = lvlmap.get(key(p.J, p.parity, p.level_energy))
            if lv is None:
                continue
            v = float(x[p.free_index])
            if p.kind == "energy":
                lv.set_energy(v)
            elif p.kind == "width":
                for c in lv.channels:
                    if (c.pair == p.pair and c.L == p.L
                            and abs(c.S - (p.S or 0.0)) < 1e-6):
                        c.gamma = v
                        break
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
