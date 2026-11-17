#!/usr/bin/env python3
"""BoundsPolicy.photon_width_max_ev / bg_photon_width_max_ev (2026-11-10): photon amplitudes have no Wigner limit and were
unbounded (a from-scratch demo: a background pole's Gamma_gamma polished to 1e12 eV).  Now |x| <= sqrt(Gamma_max / c) with
Gamma_gamma = c x^2 from the transform at x_ref (a probe at 1e-4 for a channel sitting at exactly zero).  Part 1 on a stub,
part 2 on the 13N fixture with the engine (seconds).  Run:  python3 rmfit/tests/test_photon_bound.py"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT); os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_photon_bound.py)")
    sys.exit(77)
import numpy as np
from rmfit.keys import ParamKey
from rmfit.local import BoundsPolicy, default_bounds
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)


class Stub:
    """A physical 1/2+ at 2.37 (p0 and a photon channel) and a 1/2- background at 20 (p0 and a photon channel at zero)."""
    keys = [ParamKey("width", level="1/2+#1", pair=1, L=0, S=0.5), ParamKey("width", level="1/2+#1", pair=2, L=1, S=0.5),
            ParamKey("width", level="1/2-#2", pair=1, L=1, S=0.5), ParamKey("width", level="1/2-#2", pair=2, L=1, S=0.5),
            ParamKey("energy", level="1/2+#1")]
    level_energy = {"1/2+#1": 2.37, "1/2-#2": 20.0}
    n_rmatrix = 5
    azr = None
    photon_pairs = frozenset([2])
    C = {1: 4.0, 3: 25.0}          # Gamma_gamma = C x^2 (eV) for the two photon channels

    def wigner_amplitudes(self):
        return {0: 2.0, 2: 1.0}     # particle channels only

    def wigner_widths(self, x):
        return {0: 1.0e6, 2: 1.0e6}

    def transform(self, x):
        out = np.array(x, float)
        for j, c in self.C.items():
            out[j] = c * x[j] ** 2 * (1 if x[j] >= 0 else -1)
        return out


def main():
    ev = Stub()
    x = np.array([0.1, 0.5, 0.1, 0.0, 2.37])
    pol = BoundsPolicy(theta2_cap=0.5, bg_theta2_cap=3.0, bg_energy=15.0, threshold_margin_mev=0.0)
    b = default_bounds(ev, x, pol)
    check("particle channels keep the theta2 caps", np.isclose(b.ub[0], np.sqrt(0.5) * 2.0) and np.isclose(b.ub[2], np.sqrt(3.0) * 1.0), b.ub)
    check("physical-level photon channel: |x| <= sqrt(1e4 eV / 4)", np.isclose(b.ub[1], np.sqrt(1e4 / 4.0)) and np.isclose(b.lb[1], -b.ub[1]), b.ub[1])
    check("background photon channel AT ZERO (probed): |x| <= sqrt(1e5 eV / 25)", np.isclose(b.ub[3], np.sqrt(1e5 / 25.0)), b.ub[3])
    b2 = default_bounds(ev, x, BoundsPolicy(theta2_cap=0.5, bg_energy=15.0, threshold_margin_mev=0.0, photon_width_max_ev=0.0, bg_photon_width_max_ev=0.0))
    check("0 = unbounded (the old behaviour)", not np.isfinite(b2.ub[1]) and not np.isfinite(b2.ub[3]))
    b3 = default_bounds(ev, x, pol, sign_lock={ev.keys[1]: -1})
    check("sign lock on a photon channel: one side of zero", b3.ub[1] < 0 and np.isclose(b3.lb[1], -np.sqrt(1e4 / 4.0)), (b3.lb[1], b3.ub[1]))
    b4 = default_bounds(ev, x, pol, cap_scale=4.0)
    check("cap_scale 4: a start inside the bound keeps it (the relaxation only admits a violating start)", np.isclose(b4.ub[1], np.sqrt(1e4 / 4.0)), b4.ub[1])
    x5 = x.copy(); x5[1] = 80.0
    b5 = default_bounds(ev, x5, pol, cap_scale=4.0)
    check("cap_scale 4: a start at 80 (bound 50) is admitted up to min(80, 2 x 50)", np.isclose(b5.ub[1], 80.0), b5.ub[1])
    x6 = x.copy(); x6[1] = 150.0
    b6 = default_bounds(ev, x6, pol, cap_scale=4.0)
    check("cap_scale 4: a start at 150 is walked in to 2 x 50", np.isclose(b6.ub[1], 100.0), b6.ub[1])
    # -- part 2: the engine on 13N (12C(p,g): pair 2 is the photon pair) --------------------------------------------
    from rmfit import AZURE2_ROOT
    from rmfit.engine import SingleSession
    from rmfit.model import single_variant
    from rmfit.azr3 import Model as AzrModel
    tmp = tempfile.mkdtemp(prefix="rmfit_photon_"); work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    bg = [lv for lv in mdl.active_levels if float(lv.energy) >= 15.0]
    ch = next(c for c in bg[0].channels if c.is_photon)
    ch.gamma = 1e-6; ch.channel_fixed = False                     # a (nearly) zero free photon channel on a background pole
    base = mdl.write(os.path.join(work, "base.azr"))
    with SingleSession(single_variant(base, "b", tmp), cwd=tmp, name="b") as ev:
        x = ev.x0.copy()
        ph = [j for j, k in enumerate(ev.keys[:ev.n_rmatrix]) if k.kind == "width" and k.pair == 2]
        pa = [j for j, k in enumerate(ev.keys[:ev.n_rmatrix]) if k.kind == "width" and k.pair == 1]
        pol = BoundsPolicy(bg_energy=15.0, e_window=0.5)
        b = default_bounds(ev, x, pol)
        b0 = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0, e_window=0.5, photon_width_max_ev=0.0, bg_photon_width_max_ev=0.0))
        print("  photon keys:", [(str(ev.keys[j]), round(float(x[j]), 6), round(float(b.ub[j]), 6)) for j in ph])
        check("13N: every free photon width has a finite bound", ph and all(np.isfinite(b.ub[j]) and b.ub[j] > 0 for j in ph), str([b.ub[j] for j in ph]))
        check("13N: particle widths unchanged by the photon cap", all(b.ub[j] == b0.ub[j] and b.lb[j] == b0.lb[j] for j in pa))
        check("13N: photon widths were unbounded before", all(not np.isfinite(b0.ub[j]) for j in ph))
        phys = np.asarray(ev.transform(x), float)
        inside = all(abs(x[j]) <= b.ub[j] + 1e-12 for j in ph)
        check("13N: the start lies inside the bounds", inside, str([(x[j], b.ub[j]) for j in ph]))
        # the bound maps to Gamma_max: push one physical-level photon amplitude to its bound and read the width back
        j = next(j for j in ph if ev.level_energy.get(ev.keys[j].level, 0.0) < 15.0)
        y = x.copy(); y[j] = b.ub[j]
        g = abs(float(np.asarray(ev.transform(y), float)[j]))
        check("13N: at the bound Gamma_gamma is within 20 % of photon_width_max_ev (asymptotic factor)", abs(g - 1e4) < 0.2 * 1e4, f"{g:.1f} eV")
    shutil.rmtree(tmp, ignore_errors=True)
    if failures:
        print(f"FAILED {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all photon-bound checks passed")


if __name__ == "__main__":
    main()
