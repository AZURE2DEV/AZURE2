"""Physics-informed proposals: what to try next from an incumbent vector.

Every move is a small, explicit change to a parameter vector (a sign flip, a
magnitude kick, an energy kick, a J-group transplant); structure moves (level
add / remove, dataset off) live in :mod:`rmfit.structure` because they change
the model file.  Proposals are ranked by *interference exposure* -- how much a
level's amplitude in a channel can interfere with the other levels of its
J^pi group in that channel -- so the archive's lesson (perturb one to a few
physically relevant channels, never everything) is built in.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .keys import ParamKey, canonical_signs, flip, level_order, sign_distance


@dataclass
class Move:
    kind: str                       # flip | pair_flip | magnitude | energy | transplant
    targets: List[ParamKey]
    params: dict = field(default_factory=dict)
    description: str = ""

    def to_json(self):
        import dataclasses
        return dict(kind=self.kind, targets=[dataclasses.asdict(k) for k in self.targets],
                    params=self.params, description=self.description)

    @property
    def jpi(self) -> Optional[str]:
        return self.targets[0].jpi if self.targets else None

    @property
    def level(self) -> Optional[str]:
        return self.targets[0].level if self.targets else None


def apply(x, keys: Sequence[ParamKey], move: Move) -> np.ndarray:
    pos = {k: i for i, k in enumerate(keys)}
    out = np.array(x, float, copy=True)
    if move.kind in ("flip", "pair_flip", "multi_flip"):
        for k in move.targets:
            out[pos[k]] = -out[pos[k]]
    elif move.kind == "magnitude":
        for k in move.targets:
            out[pos[k]] *= move.params["factor"]
    elif move.kind == "energy":
        out[pos[move.targets[0]]] = move.params["energy"]
    elif move.kind == "transplant":
        src = np.asarray(move.params["x_from"], float)
        for k in move.targets:
            out[pos[k]] = src[pos[k]]
    else:
        raise ValueError(move.kind)
    return out


# -- exposure -------------------------------------------------------------------

def level_energies(ev, x) -> Dict[str, float]:
    """Current energy of every level: the free value if the energy is free,
    the seed otherwise."""
    e = dict(ev.level_energy)
    for j, k in enumerate(ev.keys):
        if k.kind == "energy":
            e[k.level] = float(x[j])
    return e


def exposure(ev, x) -> Dict[ParamKey, float]:
    """Interference exposure of every free width amplitude.

    For amplitude ``gamma_{l c}`` (level ``l``, channel ``c``) the sum over the
    other levels ``l'`` of the same J^pi that have a nonzero amplitude in the
    same channel of ``|gamma_{l c} gamma_{l' c}| / (|E_l - E_l'| + 0.05 MeV)``
    -- the size of the off-diagonal R-matrix term the two levels share,
    relative to their separation.  Amplitudes with zero exposure cannot change
    the fit by a sign flip.
    """
    x = np.asarray(x, float)
    E = level_energies(ev, x)
    by_chan: Dict[tuple, List[tuple]] = {}
    for j, k in enumerate(ev.keys):
        if k.kind == "width" and x[j] != 0:
            by_chan.setdefault((k.jpi, k.pair, k.L, k.S), []).append((k, x[j]))
    out: Dict[ParamKey, float] = {}
    for chan, amps in by_chan.items():
        for k, g in amps:
            s = 0.0
            for k2, g2 in amps:
                if k2.level == k.level:
                    continue
                s += abs(g * g2) / (abs(E.get(k.level, 0.0) - E.get(k2.level, 0.0)) + 0.05)
            out[k] = s
    for j, k in enumerate(ev.keys):
        if k.kind == "width" and k not in out:
            out[k] = 0.0
    return out


# -- proposal generators ----------------------------------------------------------

def single_flips(x, keys: Sequence[ParamKey], jpi: Optional[Sequence[str]] = None,
                 exposure_map: Optional[Dict[ParamKey, float]] = None,
                 skip: Sequence[str] = (), fixed: Sequence[ParamKey] = ()) -> List[Move]:
    """Every single-amplitude sign flip that is a genuinely different sign
    pattern (modulo the exact symmetries), optionally restricted to some J^pi
    groups, ordered by exposure (highest first).  ``skip`` lists ``str(key)``
    of flips already known to revert or worsen."""
    x = np.asarray(x, float)
    jpi = set(jpi) if jpi else None
    skip = set(skip)
    c0 = canonical_signs(x, keys, fixed)
    seen = {c0}
    # a skipped flip's equivalent patterns are skipped too (flipping either
    # amplitude of a two-amplitude level is the same move)
    for j, k in enumerate(keys):
        if k.kind == "width" and x[j] != 0 and str(k) in skip:
            seen.add(canonical_signs(flip(x, keys, [k]), keys, fixed))
    moves = []
    for j, k in enumerate(keys):
        if k.kind != "width" or x[j] == 0:
            continue
        if jpi and k.jpi not in jpi:
            continue
        if str(k) in skip:
            continue
        c = canonical_signs(flip(x, keys, [k]), keys, fixed)
        if c in seen:
            continue
        seen.add(c)
        moves.append(Move("flip", [k], description=f"flip {k}"))
    if exposure_map:
        moves.sort(key=lambda m: -exposure_map.get(m.targets[0], 0.0))
    return moves



def random_flip_set(moves: Sequence[Move], rng, n: int, per_group: int = 3) -> List[Move]:
    """``n`` single flips drawn at random (``rng``: numpy Generator) from ``moves`` -- normally the highest-exposure
    single_flips -- with at most ``per_group`` per J^pi group.  Used by Campaign.run_population (2026-11-08): the
    restarts of a population must differ from one another, which perturb's deterministic top-exposure choice does not."""
    chosen: List[Move] = []
    groups: Dict[Optional[str], int] = {}
    for j in rng.permutation(len(moves)):
        mv = moves[int(j)]
        if groups.get(mv.jpi, 0) >= per_group:
            continue
        chosen.append(mv)
        groups[mv.jpi] = groups.get(mv.jpi, 0) + 1
        if len(chosen) >= n:
            break
    return chosen


def pair_flips(x, keys: Sequence[ParamKey], jpi: str, exposure_map: Dict[ParamKey, float],
               top: int = 20, fixed: Sequence[ParamKey] = ()) -> List[Move]:
    """Two-amplitude flips inside one J^pi group, the ``top`` pairs by the
    product of their exposures, excluding pairs equivalent to a single flip."""
    x = np.asarray(x, float)
    amps = [k for j, k in enumerate(keys) if k.kind == "width" and k.jpi == jpi and x[j] != 0]
    amps.sort(key=lambda k: -exposure_map.get(k, 0.0))
    amps = amps[:max(top, 12)]
    seen = {canonical_signs(x, keys, fixed)}
    for k in amps:
        seen.add(canonical_signs(flip(x, keys, [k]), keys, fixed))
    cands = []
    for i, a in enumerate(amps):
        for b in amps[i + 1:]:
            c = canonical_signs(flip(x, keys, [a, b]), keys, fixed)
            if c in seen:
                continue
            seen.add(c)
            cands.append((exposure_map.get(a, 0.0) * exposure_map.get(b, 0.0), a, b))
    cands.sort(key=lambda t: -t[0])
    return [Move("pair_flip", [a, b], description=f"flip {a} & {b}") for _, a, b in cands[:top]]


def magnitude_kicks(x, keys: Sequence[ParamKey], level: str, factors=(0.3, 3.0)) -> List[Move]:
    ks = [k for j, k in enumerate(keys) if k.kind == "width" and k.level == level and x[j] != 0]
    return [Move("magnitude", ks, params=dict(factor=f), description=f"scale {level} widths x{f:g}")
            for f in factors]


def energy_kick(keys: Sequence[ParamKey], level: str, energy: float) -> Move:
    k = [k for k in keys if k.kind == "energy" and k.level == level]
    if not k:
        raise ValueError(f"the energy of {level} is not free")
    return Move("energy", k, params=dict(energy=float(energy)), description=f"move {level} to {energy:.4f} MeV")


def transplant(keys: Sequence[ParamKey], jpi: str, x_from) -> Move:
    ks = [k for k in keys if k.kind in ("energy", "width") and k.jpi == jpi]
    return Move("transplant", ks, params=dict(x_from=np.asarray(x_from, float).tolist()),
                description=f"transplant J^pi {jpi} block from another basin")


def group_of(keys: Sequence[ParamKey], jpi: str) -> np.ndarray:
    return np.array([k.kind in ("energy", "width") and k.jpi == jpi for k in keys])


def channel_pair_flips(x, keys: Sequence[ParamKey], jpi: str, fixed: Sequence[ParamKey] = (),
                       skip_patterns: Sequence[str] = ()) -> List[Move]:
    """Every flip of two amplitudes of the *same channel* in two levels of one
    J^pi group (the pair whose relative sign sets their interference in that
    channel), excluding patterns equivalent to a single flip or already seen."""
    x = np.asarray(x, float)
    by_chan: Dict[tuple, List[ParamKey]] = {}
    for j, k in enumerate(keys):
        if k.kind == "width" and k.jpi == jpi and x[j] != 0:
            by_chan.setdefault((k.pair, k.L, k.S), []).append(k)
    seen = {canonical_signs(x, keys, fixed)} | set(skip_patterns)
    for j, k in enumerate(keys):
        if k.kind == "width" and k.jpi == jpi and x[j] != 0:
            seen.add(canonical_signs(flip(x, keys, [k]), keys, fixed))
    out = []
    for chan, ks in by_chan.items():
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                c = canonical_signs(flip(x, keys, [a, b]), keys, fixed)
                if c in seen:
                    continue
                seen.add(c)
                out.append(Move("pair_flip", [a, b], description=f"flip {a} & {b}"))
    return out


def channel_enumeration(x, keys: Sequence[ParamKey], jpi: str, fixed: Sequence[ParamKey] = (),
                        max_patterns: int = 256, skip_patterns: Sequence[str] = ()) -> List[Move]:
    """For each channel of a J^pi group with n >= 3 nonzero amplitudes (one per
    level) and 2^(n-1) <= ``max_patterns``: every relative-sign pattern of
    that channel across its levels, as a multi-amplitude flip from the current
    pattern -- the archive's "always test both relative signs of interfering
    same-J^pi levels", done exhaustively per channel.  Patterns equivalent to
    the incumbent, to a single flip, or already seen are dropped."""
    x = np.asarray(x, float)
    by_chan: Dict[tuple, List[ParamKey]] = {}
    for j, k in enumerate(keys):
        if k.kind == "width" and k.jpi == jpi and x[j] != 0:
            by_chan.setdefault((k.pair, k.L, k.S), []).append(k)
    seen = {canonical_signs(x, keys, fixed)} | set(skip_patterns)
    for j, k in enumerate(keys):
        if k.kind == "width" and k.jpi == jpi and x[j] != 0:
            seen.add(canonical_signs(flip(x, keys, [k]), keys, fixed))
    out = []
    for chan, ks in by_chan.items():
        n = len(ks)
        if n < 3 or 2 ** (n - 1) > max_patterns:
            continue
        for mask in range(1, 2 ** (n - 1)):
            targets = [ks[i + 1] for i in range(n - 1) if (mask >> i) & 1]
            if len(targets) < 2:
                continue                       # singles are proposed elsewhere
            c = canonical_signs(flip(x, keys, targets), keys, fixed)
            if c in seen:
                continue
            seen.add(c)
            out.append(Move("multi_flip", targets,
                            description=f"flip {len(targets)} of {jpi} p{chan[0]} L{chan[1]} S{chan[2]:g}: "
                                        + ", ".join(k.level for k in targets)))
    return out


def crossover_moves(x_inc, x_donor, keys: Sequence[ParamKey], fixed: Sequence[ParamKey] = (),
                    label: str = "donor") -> List[Move]:
    """Where the donor's sign pattern differs from the incumbent's (minimal flip
    set modulo the gauge symmetries): the donor's pattern of each differing
    channel, of each J^pi group that differs in more than one channel, and of
    everything at once, as multi-amplitude flips of the incumbent."""
    dist, diff = sign_distance(x_inc, x_donor, keys, fixed)
    by_chan: Dict[tuple, List[ParamKey]] = {}
    by_group: Dict[str, List[ParamKey]] = {}
    for k in diff:
        by_chan.setdefault((k.jpi, k.pair, k.L, k.S), []).append(k)
        by_group.setdefault(k.jpi, []).append(k)
    out = []
    for ch, ks in sorted(by_chan.items(), key=lambda t: -len(t[1])):
        out.append(Move("multi_flip", list(ks),
                        description=f"{label} pattern of {ch[0]} p{ch[1]} L{ch[2]} S{ch[3]:g} ({len(ks)} flips)"))
    if len(by_chan) > 1:
        for g, ks in sorted(by_group.items(), key=lambda t: -len(t[1])):
            biggest = max(len(v) for c, v in by_chan.items() if c[0] == g)
            if len(ks) > biggest:
                out.append(Move("multi_flip", list(ks),
                                description=f"{label} pattern of all {g} channels ({len(ks)} flips)"))
        if len(by_group) > 1:              # otherwise identical to that group's move
            out.append(Move("multi_flip", list(diff), description=f"{label} pattern of all channels ({len(diff)} flips)"))
    return out
