"""Sign rounds: the interference-sign search an evaluator does by hand.

One round on an incumbent vector:

1. **screen** -- for every canonical single flip (and the top pair flips) of
   the requested J^pi groups, flip and take one damped Gauss-Newton step on the
   group's columns (+ norms): ~4 sharded evaluations each; rank by the
   predicted objective;
2. **polish** the top-k: stage 1 polishes the group + norms with the flipped
   amplitude *sign-locked* away from zero (in rwa space gamma = 0 is a smooth
   point, so an unlocked refit walks straight back -- the "reverts" the
   archive kept seeing); stage 2 releases everything;
3. **greedy union** -- accepted flips are added one at a time in order of gain
   to the incumbent, each followed by a global polish, keeping only the ones
   that still help (a naive union of individually good flips made a large
   fit much worse);
4. **verify** the result on a fresh single session and record everything.

The round's *improving-flip count* is the convergence signal: the archive's
rounds went 16 -> 4 -> 3 -> 1.
"""
import time
import dataclasses
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .engine import Evaluator, EvalError
from .keys import ParamKey, canonical_signs
from .local import (Bounds, BoundsPolicy, FitResult, default_bounds, free_mask,
                    gn_step, polish, sanity)
from .moves import (Move, apply, exposure, pair_flips, single_flips, channel_pair_flips, channel_enumeration,
                    crossover_moves)


@dataclass
class SearchConfig:
    jpi: Optional[Sequence[str]] = None     # J^pi groups to search (all if None)
    top_k: int = 10                         # candidates polished per round
    pairs_per_group: int = 6                # exposure-ranked pair flips GN-screened per J^pi group
    channel_pairs: bool = True              # also every same-channel pair flip, frozen-prescreened
    enumerate_channels: bool = False        # ... and every relative-sign pattern of small channels (<= max_patterns)
    max_patterns: int = 256
    frozen_top: int = 15                    # ... of which this many *per J^pi group* get the GN screen
    gn_iters_multi: int = 3                 # GN iterations (with backtracking) when screening multi-amplitude patterns
    multi_top: int = 4                      # best multi-level channel patterns polished in addition to top_k
    flip_floor: float = 0.01                # a flip must move the frozen objective by this x tau to earn a polish
    screen_damping: float = 1e-2
    stage1_nfev: int = 40
    stage2_nfev: int = 120
    union_nfev: int = 150
    tau_rel: float = 1e-3                   # acceptance: gain >= max(tau_rel*obj, tau_sigma*s2*ln N)
    tau_sigma: float = 5.0
    dup_rel: float = 5e-4
    verbose: int = 1
    # 2026-10-02 local release: >0 = the screen's GN step and the stage-1 polish free every level of ANY J^pi within
    # +-local_release MeV of the flipped level, all background-level widths (>= the policy's bg_energy) and the norms,
    # instead of the flipped level's J^pi group; 0 = the old group-only behaviour
    local_release: float = 0.0


@dataclass
class ScreenResult:
    move: Move
    obj_flip: float          # objective right after the flip (no refit)
    obj_gn: float            # after one damped GN step on the group
    x_gn: np.ndarray

    @property
    def delta(self):
        return self.obj_gn


@dataclass
class FlipTrial:
    move: Move
    stage1: FitResult
    stage2: Optional[FitResult]
    outcome: str             # accept | marginal | worse | reverts | drift | unphysical | crashed
    x: np.ndarray
    objective: float
    gain: float


@dataclass
class RoundStat:
    n_proposed: int = 0
    n_screened: int = 0
    n_polished: int = 0
    n_improving_screen: int = 0
    n_accepted: int = 0
    n_reverted: int = 0
    n_union_kept: int = 0
    control_after: float = np.nan
    best_before: float = np.nan
    best_after: float = np.nan
    verified_after: float = np.nan
    wall: float = 0.0
    trials: List[dict] = field(default_factory=list)
    screen: List[dict] = field(default_factory=list)


def tau_flip(objective, npts, cfg: SearchConfig) -> float:
    s2 = objective / max(npts, 1)
    return max(cfg.tau_rel * objective, cfg.tau_sigma * s2 * np.log(max(npts, 2)))


def group_cols(ev: Evaluator, jpi: str, with_norms=True) -> np.ndarray:
    return np.where(free_mask(ev, jpi=[jpi], with_norms=with_norms))[0]


# -- screening ----------------------------------------------------------------------

def move_energies(ev: Evaluator, x, mv: Move) -> List[float]:
    """Current energies of the levels a move touches."""
    from .local import level_energies_now
    E = level_energies_now(ev, x)
    return sorted({E[k.level] for k in mv.targets if k.level is not None and k.level in E})


def move_cols(ev: Evaluator, x, mv: Move, cfg: SearchConfig, policy: Optional[BoundsPolicy] = None) -> np.ndarray:
    """Columns freed for a move's screen step / stage-1 polish: its J^pi group, or the local release."""
    if getattr(cfg, "local_release", 0.0) and cfg.local_release > 0:
        from .local import local_mask
        return np.where(local_mask(ev, x, move_energies(ev, x, mv), cfg.local_release,
                                   bg_energy=(policy.bg_energy if policy is not None else None)))[0]
    return group_cols(ev, mv.jpi)


def screen(ev: Evaluator, x, moves: Sequence[Move], bounds: Bounds, cfg: SearchConfig,
           log=print, policy: Optional[BoundsPolicy] = None) -> List[ScreenResult]:
    """Flip + one damped Gauss-Newton step on the group's columns per move."""
    x = np.asarray(x, float)
    out = []
    t0 = time.time()
    for i, mv in enumerate(moves):
        xf = bounds.clip(apply(x, ev.keys, mv))
        cols = move_cols(ev, xf, mv, cfg, policy)
        n_it = cfg.gn_iters_multi if mv.kind == "multi_flip" else 1
        try:
            of = ev.objective(xf)
            xg, og = xf, of
            damping = cfg.screen_damping
            for _ in range(max(n_it, 1)):
                # a coordinated flip of several amplitudes starts far from the
                # pattern's own minimum: one linear step may overshoot, so the
                # step is kept only if it improves, else retried more damped
                x1, o1, _ = gn_step(ev, xg, cols, damping=damping, bounds=bounds)
                if o1 < og:
                    xg, og = x1, o1
                    continue
                x1, o1, _ = gn_step(ev, xg, cols, damping=damping * 100, bounds=bounds)
                if o1 < og:
                    xg, og = x1, o1
                    damping *= 10
                else:
                    break
        except EvalError as err:
            log(f"    screen {mv.description}: engine failure ({err})")
            continue
        out.append(ScreenResult(mv, of, og, xg))
        if cfg.verbose > 1 or (i + 1) % 10 == 0:
            log(f"    screened {i+1}/{len(moves)} ({time.time()-t0:.0f}s) last: {mv.description}: "
                f"flip {of:.1f} -> gn {og:.1f}")
    out.sort(key=lambda s: s.obj_gn)
    return out


# -- polishing one flip ---------------------------------------------------------------

def polish_flip(ev: Evaluator, x_ref, mv: Move, policy: BoundsPolicy, cfg: SearchConfig,
                x_start=None, log=print) -> FlipTrial:
    """Stage 1 (group + norms, sign-locked) then stage 2 (everything free)."""
    x_ref = np.asarray(x_ref, float)
    x0 = apply(x_ref, ev.keys, mv) if x_start is None else np.asarray(x_start, float)
    pos = {k: i for i, k in enumerate(ev.keys)}
    lock = {k: (1 if x0[pos[k]] > 0 else -1) for k in mv.targets}
    b1 = default_bounds(ev, x_ref, policy, sign_lock=lock)
    x0 = b1.clip(x0)
    free1 = free_mask(ev, jpi=[mv.jpi], with_norms=True)
    if getattr(cfg, "local_release", 0.0) and cfg.local_release > 0:
        free1 = np.zeros(len(ev.keys), bool); free1[move_cols(ev, x0, mv, cfg, policy)] = True
    f1 = polish(ev, x0, free=free1, bounds=b1, policy=policy, max_nfev=cfg.stage1_nfev, log=log)
    if f1.status == "crashed":
        return FlipTrial(mv, f1, None, "crashed", x_ref, np.inf, -np.inf)
    b2 = default_bounds(ev, x_ref, policy)
    f2 = polish(ev, f1.x, bounds=b2, policy=policy, max_nfev=cfg.stage2_nfev, log=log)
    if f2.status == "crashed":
        f2 = f1
    reverted = any((f2.x[pos[k]] > 0) != (lock[k] > 0) for k in mv.targets)
    outcome = "reverts" if reverted else "candidate"
    return FlipTrial(mv, f1, f2, outcome, f2.x, f2.objective, np.nan)


# -- the round --------------------------------------------------------------------------

def sign_round(ev: Evaluator, x_inc, obj_inc: float, cfg: SearchConfig, policy: BoundsPolicy,
               skip: Sequence[str] = (), log=print) -> tuple:
    """One sign round from the incumbent ``x_inc`` (objective ``obj_inc``).

    Returns ``(x_best, obj_best, RoundStat, trials)`` on the *sharded*
    evaluator; the caller verifies ``x_best`` on a single session.
    """
    t_round = time.time()
    x_inc = np.asarray(x_inc, float)
    stat = RoundStat(best_before=obj_inc)
    bounds = default_bounds(ev, x_inc, policy)
    expo = exposure(ev, x_inc)
    fixed = getattr(ev, "fixed_amplitudes", ())
    # 2026-09-14: amplitudes a polish left at ~0 (|rwa| < sign_eps * gamma_W: open channels of
    # a freshly added level that the data do not use) are not flip candidates -- flipping
    # ~0 is a no-op that "improves" by noise and then eats a staged polish (the 7/2+ and 5/2-
    # post-add rounds spent 3.5 h each on such flips).  The generators skip exact zeros, so
    # they get a proposal vector with those amplitudes zeroed; screening still starts
    # from x_inc.
    x_prop = x_inc.copy()
    try:
        gW = ev.wigner_amplitudes()
        tiny = [j for j, k in enumerate(ev.keys) if k.kind == "width" and j in gW and gW[j] > 0
                and 0 < abs(x_inc[j]) < policy.sign_eps * gW[j]]
        for j in tiny:
            x_prop[j] = 0.0
        if tiny:
            log(f"  {len(tiny)} amplitude(s) below sign_eps * gamma_W excluded from the flip proposals: "
                + ", ".join(str(ev.keys[j]) for j in tiny[:8]) + (" ..." if len(tiny) > 8 else ""))
    except Exception as err:
        log(f"  tiny-amplitude exclusion skipped ({err})")
    moves = single_flips(x_prop, ev.keys, jpi=cfg.jpi, exposure_map=expo, skip=skip, fixed=fixed)
    groups = sorted({m.jpi for m in moves})
    for g in groups:
        moves += pair_flips(x_prop, ev.keys, g, expo, top=cfg.pairs_per_group, fixed=fixed)
    # same-channel pair flips are too many for the GN screen (hundreds per
    # group): a frozen evaluation each (1 forward pass), then the GN screen
    # for the best few -- the frozen number is not predictive of a refit, but
    # it is a usable filter for "which pairs are worth a Jacobian"
    extra = []
    if cfg.channel_pairs:
        from .keys import canonical_signs as _cs
        known = {_cs(apply(x_prop, ev.keys, m), ev.keys, fixed) for m in moves}
        for g in groups:
            cp = channel_pair_flips(x_prop, ev.keys, g, fixed=fixed, skip_patterns=known)
            extra += cp
            if cfg.enumerate_channels:
                # a two-level pattern of a channel *is* one of the pair flips above
                known_g = known | {_cs(apply(x_prop, ev.keys, m), ev.keys, fixed) for m in cp}
                extra += channel_enumeration(x_prop, ev.keys, g, fixed=fixed, max_patterns=cfg.max_patterns,
                                             skip_patterns=known_g)
        if extra:
            t0 = time.time()
            frozen = []
            for mv in extra:
                try:
                    frozen.append((ev.objective(bounds.clip(apply(x_inc, ev.keys, mv))), mv))
                except EvalError:
                    continue
            frozen.sort(key=lambda t: t[0])
            # per J^pi group: half of the quota by frozen objective (which
            # favours flips of small amplitudes -- they change little), half
            # by interference exposure of the flipped amplitudes, so the big
            # coordinated patterns are screened too
            import math
            chosen = []
            n_obj = (cfg.frozen_top + 1) // 2
            for g in groups:
                fg = [t for t in frozen if t[1].jpi == g]
                pick = fg[:n_obj]
                rest = fg[n_obj:]
                rest.sort(key=lambda t: -sum(math.log(expo.get(k, 0.0) + 1e-30) for k in t[1].targets))
                chosen += pick + rest[:cfg.frozen_top - n_obj]
            n_enum = sum(1 for mv in extra if mv.kind == "multi_flip")
            log(f"  frozen pre-screen of {len(extra)} same-channel patterns ({len(extra) - n_enum} pair flips, "
                f"{n_enum} enumerated) in {time.time()-t0:.0f}s; {len(chosen)} kept (per group <= {cfg.frozen_top}: "
                f"half by frozen objective, half by exposure); best frozen "
                + ", ".join(f"{mv.description} {o:.0f}" for o, mv in frozen[:3]))
            moves += [mv for _, mv in chosen]
    stat.n_proposed = len(moves) + len(extra)
    log(f"  round: {len(moves)} proposals GN-screened ({len(extra)} channel pairs pre-screened) over "
        f"{len(groups)} J^pi groups; incumbent {obj_inc:.3f}")

    return finish_round(ev, x_inc, obj_inc, moves, bounds, cfg, policy, stat, t_round, log=log)


def finish_round(ev: Evaluator, x_inc, obj_inc: float, moves: Sequence[Move], bounds: Bounds,
                 cfg: SearchConfig, policy: BoundsPolicy, stat: RoundStat, t_round: float, log=print) -> tuple:
    """The second half of any round: GN-screen the proposals, control-polish the
    incumbent, polish the best proposals (top_k plus multi_top multi-level
    patterns) sign-locked then released, classify, greedy union.  Returns
    ``(x_best, obj_best, stat, trials)`` on the sharded evaluator."""
    if getattr(cfg, "local_release", 0.0) and cfg.local_release > 0:
        log(f"  local release: screen and stage-1 polish free every level within +-{1000 * cfg.local_release:.0f} keV "
            f"of the flipped one, all background widths (>= {policy.bg_energy} MeV) and the norms")
    scr = screen(ev, x_inc, moves, bounds, cfg, log=log, policy=policy)
    stat.n_screened = len(scr)
    stat.n_improving_screen = sum(1 for s in scr if s.obj_gn < obj_inc)
    stat.screen = [dict(move=s.move.description, kind=s.move.kind, frozen=float(s.obj_flip), gn=float(s.obj_gn))
                   for s in scr]
    log(f"  screen: {stat.n_improving_screen} of {len(scr)} improve after one GN step; best "
        + ", ".join(f"{s.move.description} {s.obj_gn:.1f}" for s in scr[:5]))
    multi = [s for s in scr if s.move.kind == "multi_flip"]
    if multi:
        log(f"  best multi-level patterns after {cfg.gn_iters_multi} GN iterations: "
            + ", ".join(f"{s.move.description} {s.obj_gn:.1f}" for s in multi[:4]))
    # 2026-09-16: a flip whose FROZEN response is numerical noise cannot be worth a
    # 12-minute staged polish.  The tiny-amplitude exclusion above does not catch these:
    # measured on the 9/2+ 10.909 round, `g[9/2+#6 p4 L4 S0.5]` has theta^2 = 0.98 and
    # flipping it moves the objective by +0.24, because pairs 4/6/7/8/9 appear in NO
    # active data segment -- an unobserved channel's sign only reaches the observables
    # through second-order multiple-scattering terms.  Rather than bookkeeping which
    # pairs are observed, use the number the screen already has: require |frozen - obj|
    # to clear `flip_floor` (default 1% of tau) before a candidate can take
    # a polish slot.  The GN step can only exploit what the flip itself changed.
    tau = tau_flip(obj_inc, ev.npts, cfg)
    noise = [s for s in scr if abs(s.obj_flip - obj_inc) < cfg.flip_floor * tau]
    if noise:
        log(f"  {len(noise)} proposal(s) whose frozen response is below {cfg.flip_floor:g} x tau are not "
            f"polished (sign of a channel the data do not see): "
            + ", ".join(s.move.description for s in noise[:5]) + (" ..." if len(noise) > 5 else ""))
        scr = [s for s in scr if s not in noise]
    # the polish list: the best top_k overall, plus the best multi_top
    # multi-level channel patterns (their screened objective is rarely
    # competitive -- the magnitudes of a coordinated flip need a real refit)
    to_polish = list(scr[:cfg.top_k])
    for s in multi:
        if len([t for t in to_polish if t.move.kind == "multi_flip"]) >= cfg.multi_top:
            break
        if s not in to_polish:
            to_polish.append(s)

    # control: the incumbent itself polished with the same budget a candidate
    # gets, so a flip is credited only with what it gains *beyond* the plain
    # continuation of the incumbent (an incumbent left at max_nfev keeps
    # improving on its own by hundreds of chi2 units per polish)
    ctrl = polish(ev, x_inc, bounds=default_bounds(ev, x_inc, policy), policy=policy,
                  max_nfev=cfg.stage2_nfev, log=log)
    obj_ref = min(obj_inc, ctrl.objective)
    x_ctrl = ctrl.x if ctrl.objective < obj_inc else x_inc
    stat.control_after = float(ctrl.objective)
    log(f"  control polish of the incumbent: {obj_inc:.3f} -> {ctrl.objective:.3f} ({ctrl.status}); "
        f"flips are judged against {obj_ref:.3f}")
    trials: List[FlipTrial] = []
    for s in to_polish:
        log(f"  polish {s.move.description} (screen {s.obj_gn:.1f})")
        tr = polish_flip(ev, x_inc, s.move, policy, cfg, x_start=s.x_gn, log=log)
        stat.n_polished += 1
        if tr.outcome == "candidate":
            tr.gain = obj_ref - tr.objective
            rep = sanity(ev, tr.x, policy, default_bounds(ev, x_inc, policy), transform=False, x_ref=x_inc)
            if not rep.ok:
                tr.outcome = "unphysical"
            elif tr.gain >= tau:
                tr.outcome = "accept"
                stat.n_accepted += 1
            else:
                tr.outcome = "worse" if tr.gain <= 0 else "marginal"
        elif tr.outcome == "reverts":
            stat.n_reverted += 1
            tr.gain = obj_ref - tr.objective
            # the sign walked back, but the released polish that did so may
            # have found a genuinely better vector on the way (once a clear gain
            # from a reverted 4-level pattern): keep it as a "drift" result
            if tr.gain >= tau:
                rep = sanity(ev, tr.x, policy, default_bounds(ev, x_inc, policy), transform=False, x_ref=x_inc)
                if rep.ok:
                    tr.outcome = "drift"
                    stat.n_accepted += 1
        log(f"    -> {tr.outcome}: objective {tr.objective:.3f} (gain {tr.gain:+.1f}, tau {tau:.0f})")
        trials.append(tr)
        stat.trials.append(dict(move=s.move.description, screen=float(s.obj_gn), outcome=tr.outcome,
                                objective=float(tr.objective), gain=float(tr.gain)))

    # greedy sequential union of the accepted flips (from the control if it
    # improved the incumbent, so a round without accepted flips still returns
    # the continuation)
    x_best, obj_best = x_ctrl, obj_ref
    accepted = sorted([t for t in trials if t.outcome in ("accept", "drift")], key=lambda t: -t.gain)
    if len(accepted) == 1:
        x_best, obj_best = accepted[0].x, accepted[0].objective
        stat.n_union_kept = 1
    elif accepted:
        log(f"  greedy union of {len(accepted)} accepted flips")
        x_best, obj_best = accepted[0].x, accepted[0].objective
        stat.n_union_kept = 1
        for t in accepted[1:]:
            tr = polish_flip(ev, x_best, t.move, policy, cfg, log=log)
            if tr.outcome in ("candidate", "reverts") and tr.objective < obj_best - 0.1 * tau:
                x_best, obj_best = tr.x, tr.objective
                stat.n_union_kept += 1
                log(f"    + {t.move.description}: {obj_best:.3f}")
            else:
                log(f"    - {t.move.description} adds nothing ({tr.outcome}, {tr.objective:.3f})")
        if stat.n_union_kept > 1:
            fu = polish(ev, x_best, bounds=default_bounds(ev, x_best, policy), policy=policy,
                        max_nfev=cfg.union_nfev, log=log)
            if fu.objective < obj_best:
                x_best, obj_best = fu.x, fu.objective
    stat.best_after = obj_best
    stat.wall = time.time() - t_round
    log(f"  round done in {stat.wall/60:.1f} min: {obj_inc:.3f} -> {obj_best:.3f} "
        f"({stat.n_accepted} accepted, {stat.n_reverted} reverted, {stat.n_union_kept} kept)")
    return x_best, obj_best, stat, trials


def crossover_round(ev: Evaluator, x_inc, obj_inc: float, x_donor, cfg: SearchConfig, policy: BoundsPolicy,
                    log=print, label="donor") -> tuple:
    """Basin crossover: where the donor vector's interference-sign pattern differs
    from the incumbent's (modulo the gauge symmetries), apply the donor's pattern
    per channel, per J^pi group and all at once, and polish every proposal.  The
    oracle experiment on a large model showed that a channel pattern taken from another
    good basin can gain a great deal from an incumbent that the ordinary
    proposals (whose screens never rank such patterns) leave alone."""
    t_round = time.time()
    x_inc = np.asarray(x_inc, float)
    x_donor = np.asarray(x_donor, float)
    stat = RoundStat(best_before=obj_inc)
    bounds = default_bounds(ev, x_inc, policy)
    fixed = getattr(ev, "fixed_amplitudes", ())
    moves = crossover_moves(x_inc, x_donor, ev.keys, fixed=fixed, label=label)
    stat.n_proposed = len(moves)
    log(f"  crossover with {label}: {len(moves)} donor patterns (per channel, per J^pi group, all) "
        f"GN-screened and all polished; incumbent {obj_inc:.3f}")
    if not moves:
        stat.best_after = obj_inc
        return x_inc, obj_inc, stat, []
    cfg2 = dataclasses.replace(cfg, top_k=len(moves), multi_top=len(moves))
    return finish_round(ev, x_inc, obj_inc, moves, bounds, cfg2, policy, stat, t_round, log=log)
