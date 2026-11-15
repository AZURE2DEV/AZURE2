"""Structured identities for free parameters, stable across model variants.

pyazr hands back a free vector ``params_rwa`` whose order is the ``.azr``
file's, and a ``Parameter`` table that names each entry by the engine's own
numbering: J-group index, level ordinal within the group, engine pair number,
segment key.  Two of those change between variants of one model (a shard
holds a subset of segments, so its segment keys are renumbered; a level added
mid-file shifts the ordinals after it), so rmfit never moves values between
variants by position.  It moves them by :class:`ParamKey`:

- an energy is ``('energy', level)``,
- a width is ``('width', level, pair, L, S)`` with the *file's* pair key,
- a norm or shift is ``('norm'|'shift', seg)`` with the *base model's* segment key; the sqrt(E)
  shift coefficient is ``('shift_sqrt', seg)`` and shift-function parameter a_j is
  ``('shift_fn', seg, term=j)`` (2026-10-22),

where ``level`` is ``"<J^pi>#<ordinal>"``, the ordinal counting the *active*
levels of that J^pi in file order (which is exactly the engine's numbering,
see ``CNuc::Fill``).  Keep that numbering stable by only ever appending new
levels at the end of the file and by "removing" levels with
``AzrModel.deactivate_level`` (widths zeroed and fixed, line kept) rather than
deleting lines mid-campaign.

A project with several compound nuclei (natural targets, AZURE2 >= the multi-CN work) has one J-group per (J^pi,
nucleus), each numbered from 1, so ``"1/2-#2"`` alone names two levels there (21Ne and 23Ne in a natural-neon target, 2026-10-30).
In such a model every level label carries its nucleus, ``"<J^pi>#<ordinal>@<nucleus>"`` (``"1/2-#2@21Ne"``); a
single-nucleus model keeps the plain label, byte for byte.  :func:`level_label` is the one place labels are made;
:func:`level_group` gives the same-J^pi group of a label (``"1/2-@21Ne"``, or ``"1/2-"``).

The canonical sign pattern modded out by the model's exact sign symmetries
lives here too, because it is a function of the keys.
"""
from dataclasses import dataclass
from typing import Tuple, Dict, List, Optional, Sequence

import numpy as np


def jpi_string(J, parity) -> str:
    """``"5/2+"`` style, the same text :class:`pyazr.azrfile.AzrLevel.jpi` uses."""
    J = float(J)
    j = int(J) if J.is_integer() else f"{int(round(2 * J))}/2"
    return f"{j}{'+' if float(parity) > 0 else '-'}"


@dataclass(frozen=True, order=True)
class ParamKey:
    kind: str                       # 'energy' | 'width' | 'norm' | 'shift' | 'shift_sqrt' | 'shift_fn' | 'abundance' | 'density'
                                    # | 'doppler' | 'thickness' | 'screening'
    level: Optional[str] = None     # "5/2+#3"
    pair: Optional[int] = None      # file pair key (not the engine's number)
    L: Optional[int] = None
    S: Optional[float] = None
    seg: Optional[int] = None       # base-model <segmentsData> key
    member: str = ""                # joint fits (rmfit.mixture): the member project the key belongs to
    nucleus: Optional[Tuple[int, int]] = None   # 'abundance': (Z, A) of the compound nucleus (<abundances> block)
    term: Optional[int] = None      # 'shift_fn': j of the shift-function parameter a_j
    effect: Optional[int] = None    # 'density' / 'doppler' / 'screening': the <targetInt> effect holding it (place among the active lines;
                                    # a densitygroup's first line), as in AZURE2's density_<effect>

    @property
    def jpi(self) -> Optional[str]:
        # a member-tagged label reads "22Ne:1/2+#3"; a multi-nucleus one "1/2+#3@23Ne"
        return None if self.level is None else self.level.split("#")[0].split(":")[-1]

    @property
    def ordinal(self) -> Optional[int]:
        return None if self.level is None else int(self.level.split("#")[1].split("@")[0])

    @property
    def nucleus_tag(self) -> str:
        """``"23Ne"`` for a level label of a multi-nucleus model, else ``""``."""
        return "" if self.level is None or "@" not in self.level else self.level.split("@", 1)[1]

    def __str__(self):
        if self.kind == "energy":
            return f"E[{self.level}]"
        if self.kind == "width":
            return f"g[{self.level} p{self.pair} L{self.L} S{self.S:g}]"
        if self.kind == "weight":
            return f"w[{self.level}]"
        if self.kind == "abundance":
            return f"a[{self.nucleus[0]}_{self.nucleus[1]}]"
        if self.kind == "shift_fn":
            return f"shift_fn[{self.member + ':' if self.member else ''}seg{self.seg} a{self.term}]"
        if self.kind == "density":
            return f"density[t{self.effect}]"
        if self.kind == "doppler":
            return f"doppler[t{self.effect}]"
        if self.kind == "screening":
            return f"screening[t{self.effect}]"
        return f"{self.kind}[{self.member + ':' if self.member else ''}seg{self.seg}]"


def level_group(label: str) -> str:
    """The same-J^pi group of a level label: ``"1/2-#2@21Ne"`` -> ``"1/2-@21Ne"``, ``"5/2+#3"`` -> ``"5/2+"``
    (a member prefix ``"22Ne:"`` is kept).  Levels of one group share a J-group in the engine; levels of two nuclei
    never do, whatever their J^pi."""
    label = str(label)
    head, _, rest = label.partition("#")
    _, at, nuc = rest.partition("@")
    return head + (at + nuc if at else "")


def jgroup_nucleus_tags(m) -> Dict[int, str]:
    """``{engine J-group: "@<nucleus>"}`` when the model's pairs form more than one compound nucleus, else ``{}``.

    A J-group is nucleus-pure in the engine, so the pair of any of its widths names its nucleus (pyazr
    ``Pair.nucleus_Z/A``; older builds report none, which reads as one nucleus).  A J-group without any width
    parameter gets no tag; :func:`build_keys` refuses the duplicate key that could then arise."""
    nuc = {int(p.number): (int(getattr(p, "nucleus_Z", 0) or 0), int(getattr(p, "nucleus_A", 0) or 0))
           for p in m.pairs}
    if len(set(nuc.values())) <= 1:
        return {}
    try:
        from pyazr.nuclei import Nucleus
        name = lambda za: Nucleus(*za).label
    except ImportError:                                  # pyazr without nuclei.py
        name = lambda za: f"Z{za[0]}A{za[1]}"
    tags: Dict[int, str] = {}
    for p in m.parameters:
        if p.kind == "width" and p.jgroup is not None and p.pair is not None and int(p.pair) in nuc:
            tags.setdefault(int(p.jgroup), "@" + name(nuc[int(p.pair)]))
    return tags


def level_label(J, parity, level, jgroup=None, tags: Optional[Dict[int, str]] = None) -> str:
    """``"<J^pi>#<ordinal>"``, plus ``"@<nucleus>"`` from ``tags`` (:func:`jgroup_nucleus_tags`) in a
    multi-nucleus model."""
    tag = (tags or {}).get(int(jgroup), "") if jgroup is not None else ""
    return f"{jpi_string(J, parity)}#{int(level)}{tag}"


def build_keys(m, seg_map: Optional[Dict[int, int]] = None) -> List[ParamKey]:
    """Keys of ``m.params_rwa`` in order.

    ``seg_map`` translates this model's segment keys to the base model's (a
    shard's renumbered lines back to the lines they came from); identity if
    omitted.
    """
    pairkey = {p.number: p.key for p in m.pairs}
    tags = jgroup_nucleus_tags(m)
    n = len(m.params_rwa)
    keys: List[Optional[ParamKey]] = [None] * n
    for p in m.parameters:
        if p.fixed or p.free_index is None:
            continue
        if p.kind in ("energy", "width"):
            level = level_label(p.J, p.parity, p.level, p.jgroup, tags)
            if p.kind == "energy":
                k = ParamKey("energy", level=level)
            else:
                k = ParamKey("width", level=level, pair=int(pairkey.get(p.pair, p.pair)),
                             L=int(p.L), S=float(p.S if p.S is not None else 0.0))
        elif p.kind in ("norm", "shift", "shift_sqrt", "shift_fn"):
            sk = int(p.segment_key)
            if seg_map is not None:
                sk = int(seg_map[sk])
            k = ParamKey(p.kind, seg=sk, term=(int(p.fn_term) if p.kind == "shift_fn" else None))
        elif p.kind == "abundance":
            # one per compound nucleus of a project with an <abundances> block (AZURE2 >= feature/multi-compound-nucleus)
            k = ParamKey("abundance", nucleus=(int(p.nucleus_Z), int(p.nucleus_A)))
        elif p.kind == "density":
            # a fitted target density (AZURE2 >= feature/fit-target-density, 2026-10-30): one per densityfit line or
            # densitygroup, keyed by its <targetInt> effect number.  A shard renumbers nothing in <targetInt>, but
            # whether a shard keeps every line is untested: SingleSession only for now.
            if seg_map is not None:
                raise NotImplementedError(f"{p.name}: fitted target densities are supported in SingleSession only")
            k = ParamKey("density", effect=int(p.effect))
        elif p.kind == "doppler":
            # a fitted free-gas Doppler temperature (AZURE2 >= 2026-11-05, 'doppler T vary sigma [Tmax]' on a
            # <targetInt> line): keyed by its effect number, like a density; SingleSession only for now.
            if seg_map is not None:
                raise NotImplementedError(f"{p.name}: fitted Doppler temperatures are supported in SingleSession only")
            k = ParamKey("doppler", effect=int(p.effect))
        elif p.kind == "screening":
            # an electron-screening potential U_e (eV) of a <targetInt> line ('screening Ue vary sigma', AZURE2 dev
            # since 2026-10): keyed by its effect number; SingleSession only, as a density (shard <targetInt> untested).
            # Before 2026-11-05 build_keys stopped at this kind ("unknown parameter kind").
            if seg_map is not None:
                raise NotImplementedError(f"{p.name}: fitted screening potentials are supported in SingleSession only")
            k = ParamKey("screening", effect=int(p.effect))
        elif p.kind == "thickness":
            # a fitted sample thickness of a transmission / self-shielded-yield segment (AZURE2 >= 2026-11-05,
            # 'transmission|yield n vary sigma', 'thicknessgroup'): keyed by the segment.  Not in shards: the
            # transform exp(-n sigma) of a split segment and a thicknessgroup across shards are not designed yet.
            if seg_map is not None:
                raise NotImplementedError(f"{p.name}: fitted sample thicknesses are supported in SingleSession only")
            k = ParamKey("thickness", seg=int(p.segment_key))
        else:
            raise ValueError(f"unknown parameter kind {p.kind!r} for {p.name}")
        if p.free_index >= n:
            raise RuntimeError(f"{p.name}: free_index {p.free_index} beyond the free vector ({n})")
        keys[p.free_index] = k
    missing = [i for i, k in enumerate(keys) if k is None]
    if missing:
        raise RuntimeError(f"no key for free-vector positions {missing[:8]}")
    seen, dup = set(), []
    for k in keys:
        if k.kind not in ("energy", "width"):
            continue                    # a shard may carry one segment's norm twice (pieces of a split segment)
        if k in seen:
            dup.append(str(k))
        seen.add(k)
    if dup:
        # 2026-10-30: two free R-matrix parameters with one key would share a column everywhere keys are matched (take,
        # index_map, bounds); refuse instead of mixing two levels up
        raise RuntimeError(f"duplicate parameter keys {sorted(set(dup))[:8]}: two levels share a label "
                           f"(several compound nuclei with a J-group that has no width?)")
    return keys  # type: ignore[return-value]


def index_map(src: Sequence[ParamKey], dst: Sequence[ParamKey]) -> np.ndarray:
    """For every key of ``src`` its position in ``dst`` (-1 if absent).

    Many-to-one is allowed: two shard pieces of one split segment both carry
    that segment's norm key and both map to the single global column.
    """
    where = {}
    for j, k in enumerate(dst):
        if k in where:
            raise ValueError(f"duplicate key in destination: {k}")
        where[k] = j
    return np.array([where.get(k, -1) for k in src], dtype=int)


def take(x_src, src_keys, dst_keys, x_dst_default):
    """A copy of ``x_dst_default`` with every key shared with ``src`` taken from ``x_src``."""
    out = np.array(x_dst_default, float, copy=True)
    idx = index_map(src_keys, dst_keys)
    for i, j in enumerate(idx):
        if j >= 0:
            out[j] = x_src[i]
    return out


def level_order(keys: Sequence[ParamKey]) -> List[str]:
    """Levels in order of first appearance among the keys (file order)."""
    seen, out = set(), []
    for k in keys:
        if k.level is not None and k.level not in seen:
            seen.add(k.level)
            out.append(k.level)
    return out


def sign_pattern(x, keys: Sequence[ParamKey]) -> str:
    """Raw signs of the width parameters, one character each (``+``/``-``/``0``)."""
    return "".join(("+" if x[i] > 0 else "-" if x[i] < 0 else "0")
                   for i, k in enumerate(keys) if k.kind == "width")


def canonical_signs(x, keys: Sequence[ParamKey],
                    fixed: Sequence[ParamKey] = ()) -> str:
    """The width sign pattern modulo the model's exact sign symmetries.

    Two transformations leave every R-matrix cross section unchanged when all
    channels are particle channels: flipping every amplitude of one level
    (``R = sum gamma_c gamma_c' / (E - E_l)`` is bilinear per level), and
    flipping every amplitude of one particle pair in every level of every
    J-group (a channel phase convention, ``U -> D U D``; the cross-J
    interference terms of differential data pick up the phase twice).  So the
    gauge freedom is one sign per level and one per pair.  Flipping one channel
    inside a single J-group is *not* a symmetry of differential data and is
    deliberately not modded out.

    ``fixed`` lists the *fixed, nonzero* amplitudes of the model (they are not
    in ``x``).  A gauge transformation may not change them, which ties the
    level sign to the pair sign for every such amplitude: on one model
    flipping every free pair-1 amplitude is *not* a symmetry (13 %
    relative), because some levels carry fixed pair-1 amplitudes, while
    flipping every free pair-1 *and* pair-2 amplitude is one to 1e-16.
    Levels and pairs joined by fixed amplitudes therefore share one gauge
    sign (connected components of that graph).  Zero amplitudes carry no sign.

    Canonical form (2026-11-12): the gauge components are nodes and every
    nonzero amplitude an edge.  In key order, an amplitude that joins two
    so-far separate trees is a spanning-forest edge: it is written ``+`` and
    the two trees are merged with the relative gauge that makes it so.  Any
    other amplitude closes a cycle and is written as its gauge-invariant
    product with the forest path between its ends.  The forest depends only on
    which amplitudes are nonzero and on the key order, so every pattern of one
    orbit gives the same string, and different orbits give different strings
    (the string is the pattern itself in one particular gauge).  Before
    2026-11-12 an amplitude joining two trees was written with the trees'
    unreconciled gauges, so about 8 % of random patterns gave a string that a
    gauge-equivalent pattern did not (missed duplicates; never a false match).
    """
    x = np.asarray(x, float)
    widths = [(i, k) for i, k in enumerate(keys) if k.kind == "width"]
    # gauge components: level and pair nodes joined by fixed nonzero amplitudes
    comp: Dict[object, object] = {}

    def croot(u):
        comp.setdefault(u, u)
        while comp[u] != u:
            comp[u] = comp[comp[u]]
            u = comp[u]
        return u

    for k in fixed:
        if k.kind == "width" and k.level is not None and k.pair is not None:
            ru, rv = croot(("L", k.level)), croot(("P", k.member, k.pair))
            if ru != rv:
                comp[ru] = rv
    # spanning forest over the components, with each node's sign relative to its tree root
    parent: Dict[object, object] = {}
    par: Dict[object, int] = {}

    def find(u):
        if u not in parent:
            parent[u], par[u] = u, 1
            return u, 1
        path = []
        while parent[u] != u:
            path.append(u)
            u = parent[u]
        root, acc = u, 1
        for v in reversed(path):            # compress: parity of v to the root
            acc *= par[v]
            parent[v], par[v] = root, acc
        return root, (par[path[0]] if path else 1)

    out = []
    for i, k in widths:
        if x[i] == 0:
            out.append("0")
            continue
        s = 1 if x[i] > 0 else -1
        A, B = croot(("L", k.level)), croot(("P", k.member, k.pair))
        if A == B:
            out.append("+" if s > 0 else "-")     # gauge-invariant
            continue
        (ra, pa), (rb, pb) = find(A), find(B)
        if ra != rb:
            parent[ra], par[ra] = rb, pa * pb * s     # makes pa' * pb * s = +1 for this edge
            out.append("+")
        else:
            out.append("+" if pa * pb * s > 0 else "-")
    return "".join(out)


def flip(x, keys: Sequence[ParamKey], which: Sequence[ParamKey]) -> np.ndarray:
    """A copy of ``x`` with the sign of each listed width parameter reversed."""
    out = np.array(x, float, copy=True)
    pos = {k: i for i, k in enumerate(keys)}
    for k in which:
        out[pos[k]] = -out[pos[k]]
    return out


def sign_distance(x1, x2, keys: Sequence[ParamKey], fixed: Sequence[ParamKey] = (),
                  exact_max: int = 20, starts: int = 64):
    """Minimal number of single-amplitude sign flips that turn the width sign
    pattern of ``x1`` into that of ``x2`` modulo the exact gauge symmetries
    (one sign per level, one per particle pair, tied by fixed nonzero
    amplitudes exactly as in :func:`canonical_signs`).  Amplitudes that are zero
    in either vector are ignored.  Returns ``(distance, differing_keys)``, the
    keys in key order; flipping ``differing_keys`` in ``x1`` gives a pattern
    gauge-equivalent to ``x2`` (also when the distance is not minimal, below).

    The problem is the frustration index of a signed graph: the gauge
    components (connected components of the level/pair graph joined by fixed
    amplitudes) are the nodes, every amplitude is an edge between its level's
    and its pair's component with sign +1 (x1, x2 agree) or -1, and a gauge is
    one sign per node.  Every edge has a pair end, so the components holding a
    pair end ("core") cover all edges and each other component is optimal in
    closed form, ``(n_v - |sum_e d_e h_pair(e)|) / 2``.  With at most
    ``exact_max`` core components the core gauges are enumerated exactly
    (vectorised, one global sign fixed); above that the distance is the best of
    ``starts`` deterministic multistart local searches (single and pair flips
    of core gauges, from the identity, a spanning-tree parity gauge and random
    gauges) -- an upper bound that was exact on every test instance, and
    distance 0 is always detected exactly (a balanced graph's parity gauge has
    no frustrated edge).

    Until 2026-11-12 the pair gauges were enumerated in a plain Python loop over
    all 2^(pairs) choices, each rescanning every amplitude per level: fine for
    the handful of pairs of the earlier campaigns, but a model with 27 pairs
    (134M choices) hung a population round for hours
    (an archived campaign readme).  That version also gave several
    levels tied to the same amplitude-less pair independent signs.
    """
    x1 = np.asarray(x1, float)
    x2 = np.asarray(x2, float)
    parent: Dict[object, object] = {}

    def find(u):
        parent.setdefault(u, u)
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    def union(u, v):
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv

    for k in fixed:
        if k.kind == "width" and k.level is not None and k.pair is not None:
            union(("L", k.level), ("P", k.member, k.pair))
    amps = [(i, k) for i, k in enumerate(keys)
            if k.kind == "width" and x1[i] != 0 and x2[i] != 0]
    if not amps:
        return (0, [])
    # pair nodes are (member, pair), as the unions above tag them (a joint fit
    # has one gauge per member's pair); plain fits have member "" throughout.
    ea = [find(("L", k.level)) for _, k in amps]
    eb = [find(("P", k.member, k.pair)) for _, k in amps]
    dd = np.array([1 if (x1[i] > 0) == (x2[i] > 0) else -1 for i, _ in amps], int)
    core_ids = {r: j for j, r in enumerate(sorted(set(eb), key=str))}
    free_ids = {r: j for j, r in enumerate(sorted({a for a in ea if a not in core_ids}, key=str))}
    R, F = len(core_ids), len(free_ids)
    const = np.array([a == b for a, b in zip(ea, eb)])          # both ends in one component: gauge-invariant
    cc = np.array([(a in core_ids) and a != b for a, b in zip(ea, eb)])
    fr = np.array([a in free_ids for a in ea])
    ia = np.array([core_ids.get(a, -1) for a in ea])
    ib = np.array([core_ids[b] for b in eb])
    fa = np.array([free_ids.get(a, -1) for a in ea])
    c_ia, c_ib, c_d = ia[cc], ib[cc], dd[cc]
    M = np.zeros((F, R))                                      # free component x core component: sum of d
    np.add.at(M, (fa[fr], ib[fr]), dd[fr])
    n_free = np.bincount(fa[fr], minlength=F).astype(float)
    base = float(np.sum(dd[const] < 0))

    def cost(H):                                              # H: (m, R) of +-1 -> (m,) frustrated edges
        H = np.atleast_2d(H)
        c = np.full(H.shape[0], base)
        if c_d.size:
            c += np.sum(H[:, c_ia] * H[:, c_ib] * c_d < 0, axis=1)
        if F:
            c += np.sum(n_free - np.abs(H @ M.T), axis=1) / 2.0
        return c

    if R <= exact_max:
        best_h, best_c = None, np.inf
        n = max(R - 1, 0)                                     # h_0 = +1: the global flip is a symmetry
        chunk = 1 << min(n, 16)
        bits = np.arange(n)
        for lo in range(0, 1 << n, chunk):
            idx = np.arange(lo, min(lo + chunk, 1 << n))
            H = np.ones((idx.size, R))
            if n:
                H[:, 1:] = 1 - 2 * ((idx[:, None] >> bits) & 1)
            c = cost(H)
            j = int(np.argmin(c))
            if c[j] < best_c:
                best_c, best_h = c[j], H[j].copy()
    else:
        # spanning-tree parity gauge over all components (zero cost if balanced)
        nodes = list(core_ids) + list(free_ids)
        adj: Dict[object, list] = {u: [] for u in nodes}
        for a, b, s in zip(ea, eb, dd):
            if a != b:
                adj[a].append((b, s))
                adj[b].append((a, s))
        sg: Dict[object, int] = {}
        for u0 in nodes:
            if u0 in sg:
                continue
            sg[u0], stack = 1, [u0]
            while stack:
                u = stack.pop()
                for v, s in adj[u]:
                    if v not in sg:
                        sg[v] = sg[u] * s
                        stack.append(v)
        rng = np.random.default_rng(0)
        inits = [np.ones(R), np.array([sg[r] for r in core_ids], float)]
        inits += [rng.choice([-1.0, 1.0], R) for _ in range(max(int(starts) - 2, 0))]
        flips = [np.eye(R)[i] for i in range(R)] + [np.eye(R)[i] + np.eye(R)[j] for i in range(R) for j in range(i + 1, R)]
        Fm = 1 - 2 * np.array(flips)                          # rows: multiply h by these to flip one or two gauges
        best_h, best_c = None, np.inf
        for h in inits:
            c0 = cost(h)[0]
            while True:
                cand = cost(Fm * h)
                j = int(np.argmin(cand))
                if cand[j] >= c0 - 1e-9:
                    break
                h, c0 = Fm[j] * h, cand[j]
            if c0 < best_c:
                best_c, best_h = c0, h.copy()
            if best_c == 0:
                break
    h = best_h
    # the gauge of each free component given the core gauges, then the frustrated edges
    g_free = np.sign(M @ h) if F else np.zeros(0)
    g_free[g_free == 0] = 1
    diff = []
    for e, (i, k) in enumerate(amps):
        if const[e]:
            bad = dd[e] < 0
        else:
            ga = g_free[fa[e]] if fr[e] else h[ia[e]]
            bad = dd[e] * ga * h[ib[e]] < 0
        if bad:
            diff.append(k)
    assert len(diff) == int(round(best_c)), (len(diff), best_c)
    return (len(diff), diff)
