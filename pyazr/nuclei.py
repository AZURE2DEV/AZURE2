"""Compound nuclei of a project: which nucleus a pair forms, and merging several
single-nucleus projects into one file (natural-abundance targets).

A project can hold the level schemes of several compound nuclei.  No column
says which nucleus a level belongs to: the pair columns of its channel lines
do, by the rule the engine applies (``PPair``):

* particle pair (pType 0):  Z = Z1 + Z2,  A = round(M1 + M2)
* gamma pair (pType 10):    Z = Z2,       A = round(M2)
* beta pair (pType 20):     Z = Z2 - Z1,  A = round(M2)   (Z1 = -1 for beta-minus)

Pairs with the same (Z, A) form one nucleus; the engine keeps their J^pi
groups apart, so the level matrix is block-diagonal.  A segment line ending in
``sample 1`` (older files: ``natural 1``) weights its own term and each composite component by the
isotopic abundance of the term's entrance pair's nucleus, read from the
``<abundances>`` block (one line per nucleus: ``Z A abundance vary sigma``).

The identity in the file is the COMPOUND nucleus: n + 20Ne forms 21Ne, n + 22Ne
forms 23Ne, n + 21Ne forms 22Ne.  The abundance belongs to the target isotope,
so a label shows both -- :meth:`Nucleus.label` gives the compound, and
:func:`AzrModel.nuclei` the pairs that form it.
"""
import os
from typing import NamedTuple

__all__ = ["SYMBOLS", "Nucleus", "nucleus_of_channel", "nucleus_of_pair", "merge_projects", "MergeReport"]

#: Element symbols by charge number (0 = the neutron).
SYMBOLS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


class Nucleus(NamedTuple):
    """A compound nucleus, identified by its charge and mass numbers."""
    Z: int
    A: int

    @property
    def symbol(self) -> str:
        return SYMBOLS[self.Z] if 0 <= self.Z < len(SYMBOLS) else f"Z{self.Z}"

    @property
    def label(self) -> str:
        """``"23Ne"`` style label of the compound nucleus."""
        return f"{self.A}{self.symbol}"

    @classmethod
    def parse(cls, spec):
        """``Nucleus``, ``(Z, A)`` or a label such as ``"23Ne"`` / ``"Ne23"``."""
        if isinstance(spec, Nucleus):
            return spec
        if isinstance(spec, (tuple, list)) and len(spec) == 2:
            return cls(int(spec[0]), int(spec[1]))
        s = str(spec).strip()
        digits = "".join(ch for ch in s if ch.isdigit())
        letters = "".join(ch for ch in s if ch.isalpha())
        if not digits or not letters:
            raise ValueError(f"cannot read a nucleus from {spec!r}; give (Z, A) or e.g. '23Ne'")
        sym = letters[0].upper() + letters[1:].lower()
        if sym not in SYMBOLS:
            raise ValueError(f"unknown element symbol {letters!r} in {spec!r}")
        return cls(SYMBOLS.index(sym), int(digits))

    def __str__(self):
        return self.label


def nucleus_of(ptype, Z1, Z2, M1, M2) -> Nucleus:
    """The engine's rule, from a pair's type, charges and masses."""
    if ptype == 20:
        return Nucleus(int(Z2) - int(Z1), int(round(M2)))
    if ptype == 10:
        return Nucleus(int(Z2), int(round(M2)))
    return Nucleus(int(Z1) + int(Z2), int(round(M1 + M2)))


def nucleus_of_channel(channel) -> Nucleus:
    """The compound nucleus of an :class:`pyazr.azrfile.AzrChannel`'s pair."""
    return nucleus_of(channel.ptype, channel.Z1, channel.Z2, channel.M1, channel.M2)


def nucleus_of_pair(pair) -> Nucleus:
    """The compound nucleus of a runtime :class:`pyazr.parameters.Pair` (the
    engine reports it; this recomputes it from the pair's columns)."""
    return nucleus_of(pair.ptype, pair.Z1, pair.Z2, pair.M1, pair.M2)


class MergeReport(NamedTuple):
    """What :func:`merge_projects` did: the pair-key, data-segment and
    test-segment maps per member (member index -> {old: new}) and warnings."""
    pair_map: list
    segment_map: list
    test_map: list
    warnings: list


def merge_projects(models, out, abundances, sum_sigma=0.0, shared=(), natural=True):
    """Merge several single-nucleus projects into one file.

    ``models`` are :class:`pyazr.azrfile.AzrModel` instances (or paths); every
    member's pair keys are offset so they stay unique, its segment lines, their
    composite components, the ``<targetInt>`` segment lists and the
    ``<potential>`` ``pair=`` sections are renumbered accordingly, the
    ``<levels>`` blocks are concatenated and ``<abundances>`` is written from
    ``abundances`` ({nucleus label or (Z, A): (value, vary, sigma)} or a plain
    value).  Data paths are rewritten relative to ``out``'s directory.

    ``shared`` names data sets the members have in common, one dict per set
    mapping member index -> that member's data-segment key: the first member's
    line is kept as the carrier, turned into a SUM composite whose components
    are the other members' main terms (and their components), marked
    ``sample 1`` when ``natural`` is set, and the other members' lines are
    dropped.

    Not merged, reported in the returned :class:`MergeReport`: the members'
    ``<parameterSettings>`` (names refer to GUI rows) and ``<mcmc>`` blocks,
    members' ``<lastRun>`` blocks other than member 0's (kept verbatim: its
    pair keys do not move), and member engine flags (one engine has one set).
    """
    from .azrfile import AzrModel
    import shlex
    from .datasets import _split_keyword_tail

    def _quoted_tokens(line):
        # non-POSIX mode keeps the quotes on a quoted token (a data path with a
        # space, a <targetInt> stopping-power equation), so the line is written
        # back as the GUI and the engine read it
        return shlex.split(line, posix=False)

    out = os.path.abspath(out)
    out_dir = os.path.dirname(out)
    members = [m if isinstance(m, AzrModel) else AzrModel.from_file(m) for m in models]
    sources = [m.source if getattr(m, "source", None) else None for m in members]
    warnings = []

    # 1. pair-key offsets: each member's keys are shifted past the previous members' largest key
    pair_map, offset = [], 0
    for m in members:
        keys = sorted({c.pair for lv in m.levels for c in lv.channels})
        pair_map.append({k: k + offset for k in keys})
        offset += max(keys) if keys else 0

    # 2. levels: every channel line with its pair key offset, level IDs renumbered on write
    merged_levels = []
    for m, pm in zip(members, pair_map):
        for lv in m.levels:
            new = []
            for c in lv.channels:
                cc = c.clone()
                cc._set("ir", pm[c.pair])
                new.append(cc)
            merged_levels.append(type(lv)(new))

    # 3. segment lines, renumbered by position; shared groups folded into one carrier line
    def remap_line(tokens, pm, is_test):
        """Offset the entrance/exit keys and the component keys of one tokenized segment line."""
        t = list(tokens)
        for j in (1, 2):
            if int(float(t[j])) != -1:
                t[j] = str(pm[int(float(t[j]))])
        # the composite block: isAdvanced op [-1] n {e x a [s]}* after the data file (data) / isDiff fields (test)
        tail_start = None
        if is_test:
            isDiff = int(float(t[9]))
            tail_start = 10 + (2 if isDiff == 2 else 1 if isDiff == 3 else 0)
        else:
            isDiff = int(float(t[7]))
            i = 8 + (2 if isDiff == 2 else 0) + 3   # norm, varyNorm, normErr
            if i + 2 < len(t) and all(_isnum(x) for x in t[i:i + 3]):
                i += 3                                  # shift, shiftErr, varyShift
            while i < len(t) and _isnum(t[i]):
                i += 1
            tail_start = i + 1                           # after the data file token
        positional, blocks = _split_keyword_tail(t[tail_start:])
        if positional and _isnum(positional[0]) and int(float(positional[0])) == 1 and len(positional) >= 3:
            j, scaled = 2, False
            n = int(float(positional[j])); j += 1
            if n == -1:
                scaled = True
                n = int(float(positional[j])); j += 1
            for _ in range(n):
                positional[j] = str(pm[int(float(positional[j]))])
                positional[j + 1] = str(pm[int(float(positional[j + 1]))])
                j += 3 + (1 if scaled else 0)
        tail = positional + [tok for kw, vals in blocks.items() for tok in [kw] + list(vals)]
        return t[:tail_start] + tail, blocks

    def _isnum(x):
        try:
            float(x); return True
        except ValueError:
            return False

    def block_lines(text, tag):
        lines = text.splitlines()
        try:
            a = lines.index(f"<{tag}>"); b = lines.index(f"</{tag}>")
        except ValueError:
            return []
        return [ln for ln in lines[a + 1:b] if ln.strip()]

    # data segments
    data_lines, segment_map = [], [dict() for _ in members]
    carriers = {}   # (member, key) -> group index
    for g, grp in enumerate(shared):
        for mi, key in grp.items():
            carriers[(int(mi), int(key))] = g
    group_lines = [None] * len(shared)
    group_extra = [[] for _ in shared]
    for mi, (m, pm) in enumerate(zip(members, pair_map)):
        member_dir = os.path.dirname(os.path.abspath(sources[mi])) if sources[mi] else out_dir
        for key, raw in enumerate(block_lines(m._suffix, "segmentsData"), start=1):
            toks, blocks = remap_line(_quoted_tokens(raw), pm, False)
            # the data file, relative to the merged file
            isDiff = int(float(toks[7]))
            i = 8 + (2 if isDiff == 2 else 0) + 3
            if i + 2 < len(toks) and all(_isnum(x) for x in toks[i:i + 3]):
                i += 3
            while i < len(toks) and _isnum(toks[i]):
                i += 1
            path = toks[i].strip('"')
            if not os.path.isabs(path):
                path = os.path.relpath(os.path.join(member_dir, path), out_dir)
            toks[i] = f'"{path}"' if " " in path else path
            g = carriers.get((mi, key))
            if g is None:
                data_lines.append(" ".join(toks))
                segment_map[mi][key] = len(data_lines)
            elif group_lines[g] is None:
                group_lines[g] = (mi, key, toks)
                data_lines.append(None)                      # placeholder, filled below
                segment_map[mi][key] = len(data_lines)
            else:
                # another member's term: becomes components of the carrier
                positional, _ = _split_keyword_tail(toks[i + 1:])
                comps = [(toks[1], toks[2], "-999")]
                if positional and int(float(positional[0])) == 1:
                    j, scaled = 2, False
                    n = int(float(positional[j])); j += 1
                    if n == -1:
                        scaled = True; n = int(float(positional[j])); j += 1
                    for _ in range(n):
                        comps.append((positional[j], positional[j + 1], positional[j + 2]))
                        j += 3 + (1 if scaled else 0)
                group_extra[g].extend(comps)
                segment_map[mi][key] = segment_map[group_lines[g][0]][group_lines[g][1]]
    for g, carrier in enumerate(group_lines):
        if carrier is None:
            continue
        mi, key, toks = carrier
        isDiff = int(float(toks[7]))
        i = 8 + (2 if isDiff == 2 else 0) + 3
        if i + 2 < len(toks) and all(_isnum(x) for x in toks[i:i + 3]):
            i += 3
        while i < len(toks) and _isnum(toks[i]):
            i += 1
        positional, blocks = _split_keyword_tail(toks[i + 1:])
        comps = []
        if positional and int(float(positional[0])) == 1:
            if int(float(positional[1])) == 1:
                raise ValueError(f"shared data set {g}: the carrier line (member {mi}, key {key}) is a RATIO")
            j, scaled = 2, False
            n = int(float(positional[j])); j += 1
            if n == -1:
                scaled = True; n = int(float(positional[j])); j += 1
            for _ in range(n):
                comps.append((positional[j], positional[j + 1], positional[j + 2]))
                j += 3 + (1 if scaled else 0)
        comps.extend(group_extra[g])
        new_tail = ["1", "0", str(len(comps))] + [tok for c in comps for tok in c] + ["0"]
        if natural:
            blocks.pop("natural", None)
            blocks["sample"] = ["1"]
        new_tail += [tok for kw, vals in blocks.items() for tok in [kw] + list(vals)]
        data_lines[segment_map[mi][key] - 1] = " ".join(toks[:i + 1] + new_tail)

    # test segments
    test_lines, test_map = [], [dict() for _ in members]
    for mi, (m, pm) in enumerate(zip(members, pair_map)):
        for key, raw in enumerate(block_lines(m._suffix, "segmentsTest"), start=1):
            toks, _ = remap_line(raw.split(), pm, True)
            test_lines.append(" ".join(toks))
            test_map[mi][key] = len(test_lines)

    # 4. target effects: segment lists remapped per member; effects whose lines vanished are dropped
    target_lines = []
    for mi, m in enumerate(members):
        for raw in block_lines(m._suffix, "targetInt"):
            toks = _quoted_tokens(raw)
            keys = _expand_ranges(toks[1].strip('"'))
            new_keys = sorted({segment_map[mi][k] for k in keys if k in segment_map[mi]})
            if not new_keys:
                warnings.append(f"member {mi}: target effect on segments {toks[1]} dropped (no surviving line)")
                continue
            toks[1] = '"' + ",".join(str(k) for k in new_keys) + '"'
            target_lines.append(" ".join(toks))

    # 5. potential block: member 0's defaults, every member's pair sections offset
    potential_lines = []
    for mi, m in enumerate(members):
        for raw in block_lines(m._prefix + m._suffix, "potential"):
            s = raw.strip()
            if s.startswith("pair="):
                k = int(s.split("=", 1)[1].split()[0])
                potential_lines.append(f"pair={pair_map[mi][k]}")
            elif mi == 0:
                potential_lines.append(s)
            elif s and s not in potential_lines:
                warnings.append(f"member {mi}: potential setting {s!r} ignored (member 0's defaults apply)")

    # 5b. measured scattering lengths: every member's lines, pair keys offset
    scatlen_lines = []
    for mi, m in enumerate(members):
        for e in m.scattering_lengths:
            scatlen_lines.append(e._replace(pair=pair_map[mi][e.pair]).line())

    # 6. assemble: member 0's config (output dir to out's), merged levels, abundances, segments, targetInt
    base = members[0]
    text = base._prefix
    text = text[:text.index("<levels>")]
    config = text[text.index("<config>"):text.index("</config>")].rstrip("\n").split("\n")
    config = [ln for ln in config if not ln.strip().startswith("abundanceSumSigma=")]
    if sum_sigma > 0:
        config.append(f"abundanceSumSigma={sum_sigma:g}")
    merged = AzrModel("\n".join(config) + "\n</config>\n" + ("<potential>\n" + "\n".join(potential_lines) + "\n</potential>\n" if potential_lines else "") + "<levels>\n",
                      merged_levels, "")
    abund_lines = []
    nuclei = merged.nuclei()
    for nuc in nuclei:
        spec = abundances.get(nuc.label, abundances.get((nuc.Z, nuc.A), abundances.get(nuc))) if isinstance(abundances, dict) else None
        if spec is None:
            value, vary, sigma = 1.0, 0, 0.0
            warnings.append(f"no abundance given for {nuc.label}: 1.0, fixed")
        elif isinstance(spec, (tuple, list)):
            value, vary, sigma = (list(spec) + [0, 0.0])[:3]
        else:
            value, vary, sigma = float(spec), 0, 0.0
        abund_lines.append(f"{nuc.Z} {nuc.A} {value:g} {int(bool(vary))} {sigma:g}")
    # member 0's <lastRun> block is carried verbatim (its pair keys are unchanged); an EMPTY
    # <lastRun> or <mcmc> block is not valid for the GUI's readers, so neither is written
    # when there is nothing to put in it.
    last_run = block_lines(base._suffix, "lastRun")
    suffix = ["</levels>", "<abundances>"] + abund_lines + ["</abundances>", "<segmentsData>"] + \
             [ln for ln in data_lines if ln is not None] + ["</segmentsData>", "<segmentsTest>"] + test_lines + \
             ["</segmentsTest>"] + \
             (["<scatteringLengths>"] + scatlen_lines + ["</scatteringLengths>"] if scatlen_lines else []) + \
             ["<targetInt>"] + target_lines + ["</targetInt>",
             "<parameterSettings>", "</parameterSettings>"] + \
             (["<lastRun>"] + last_run + ["</lastRun>"] if last_run else [])
    merged._suffix = "\n".join(suffix) + "\n"
    merged.set_output_dir(os.path.join(out_dir, "output"))
    for mi, m in enumerate(members):
        if block_lines(m._suffix, "parameterSettings"):
            warnings.append(f"member {mi}: <parameterSettings> dropped (its names refer to GUI rows)")
        if block_lines(m._suffix, "mcmc"):
            warnings.append(f"member {mi}: <mcmc> block dropped (its parameters are named by position)")
        if mi > 0 and block_lines(m._suffix, "lastRun"):
            warnings.append(f"member {mi}: <lastRun> block dropped (member 0's is kept)")
    merged.write(out)
    return MergeReport(pair_map, segment_map, test_map, warnings)


def _expand_ranges(spec):
    """``"3-5,8"`` -> [3, 4, 5, 8]."""
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out
