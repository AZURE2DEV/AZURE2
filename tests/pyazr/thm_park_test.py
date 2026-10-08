#!/usr/bin/env python3
"""THM under Park's parametrization: Brune equivalence and the THM derivatives.

`--use-park` (pyazr `use_park=True`) fits Park's amplitudes,
gamma_Park = gamma_Brune sqrt(J), J = 1 - sum_c gamma_Park^2 dS_c/dE.  Park's
level matrix is Brune's with diag(sqrt J) on both sides, and the HOES amplitude
sum gamma_f A v (v = sum gamma M_l) is bilinear in the amplitudes, so a THM
model has to be the same function of the physics in both parametrizations.
The THM rows of the Jacobian are central differences of the HOES model in the
fit's own amplitudes, so under Park they must carry the chain rule through J
(which depends on every amplitude of the level and, through dS/dE, on its
energy).  Checked on tests/18O_p_a_thm (two interfering 1/2+ levels, every
channel entered as a Brune amplitude, gammaIsRWA), unfolded, the first level's
energy and all four amplitudes free, in two variants:

  A. a THM experiment with a linear background (profiled norm and background)
     and a coherent background cbackground=1/2+:2 (two more free parameters);
  B. the Coulomb line shape (lineshape=on), whose N_C takes each level's total
     width from the fit amplitudes (ThmLevelWidth).

For each variant:

  1. the Park and Brune sessions from the same .azr give the same chi2 and
     residuals (rel 1e-9), their free vectors are related by
     gamma_Brune = gamma_Park / sqrt(J) with park_norms (rel 1e-12), and
     transform_rwa (what save_fit writes) gives the same physical values, the
     amplitude inputs as Brune's amplitudes;
  2. residual_jacobian (Park) agrees with central differences of the
     residuals (rel 1e-4 of the largest entry);
  3. the chain rule: the Park Jacobian equals the Brune session's Jacobian
     times d x_Brune / d x_Park (rel 1e-5) -- the derivative with respect to a
     Park amplitude includes the J factor;
  4. chi2_and_grad (Park) equals 2 J^T r and central differences of chi2.

Then the J > 0 wall on a THM project: a level's amplitudes scaled until
J < 0 raise calculate_chi2_rwa above the data chi2 by sum (J/1e-3)^2, which
park_norms and penalties report, and the gradient there matches differences
of the chi2 (1e-2: see park_gradient_test.py).

Needs the compiled engine; exits 77 (skipped) without it.

Run from anywhere:  python3 tests/pyazr/thm_park_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
KIN = "beam=18O target=3He spectator=d Ebeam=115"
VARIANTS = {
    "A (experiment, linear and coherent background)":
        "experiment[A] segments=1,2 background=linear cbackground=1/2+:2=0.3,-0.2",
    "B (lineshape=on)": f"experiment[A] segments=1,2 {KIN} lineshape=on",
}

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    sys.stdout.reconfigure(line_buffering=True)
    import numpy as np
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def stage(tmp, name, block):
    """The project unfolded, level 1 (8.6026 MeV) energy free and 10 keV high,
    every amplitude free, with BLOCK as its <thm> block."""
    work = os.path.join(tmp, name)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    out, skip, levels = [], False, False
    with open(os.path.join(SOURCE, "18O_p_a_thm.azr")) as f:
        for line in f.read().replace("\r", "").split("\n"):
            if line.startswith("<targetInt>"):
                out.append(line)
                skip = True
                continue
            if line.startswith("</targetInt>"):
                skip = False
            if skip:
                continue
            if line.startswith("<levels>"):
                levels = True
            elif line.startswith("</levels>"):
                levels = False
            t = line.split()
            if levels and len(t) > 30:
                if t[2] == "8.602600":
                    t[2], t[3] = "8.612600", "0"
                t[10] = "0"
                line = " ".join(t)
            out.append(line)
    out.append(f"<thm>\n{block}\n</thm>\n")
    path = os.path.join(work, "run.azr")
    with open(path, "w") as f:
        f.write("\n".join(out))
    return path


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    scale = max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300)
    return float(np.max(np.abs(a - b)) / scale)


def level_of(m):
    """Per free parameter: the index of its level in park_norms() order for a
    width, None otherwise."""
    order = {key: i for i, key in enumerate(m.physical_levels())}
    out = [None] * len(m.params_rwa)
    for key, ps in m.parameters.by_physical_level().items():
        for p in ps:
            if p.kind == "width" and not p.fixed and p.free_index is not None:
                out[p.free_index] = order[key]
    return out


def to_brune(m, levels, x):
    """gamma_Brune = gamma_Park / sqrt(J) for the widths of x (Park)."""
    J = m.park_norms(x)
    y = np.array(x, float)
    for i, lv in enumerate(levels):
        if lv is not None:
            y[i] = x[i] / np.sqrt(J[lv])
    return y


def residual_fd(m, x, h):
    J = np.zeros((m.residual_jacobian(x)[0].size, x.size))
    for i in range(x.size):
        xp = x.copy(); xp[i] += h[i]
        xm = x.copy(); xm[i] -= h[i]
        J[:, i] = (m.residual_jacobian(xp)[0] - m.residual_jacobian(xm)[0]) / (2.0 * h[i])
    return J


def chi2(m, x):
    return float(np.sum(m.calculate_chi2_rwa(x)))


with tempfile.TemporaryDirectory() as tmp:
    for k, (label, block) in enumerate(VARIANTS.items()):
        print(f"\n{label}")
        azr_b = stage(tmp, f"v{k}.brune", block)
        azr_p = stage(tmp, f"v{k}.park", block)
        with azure2(azr_b, use_brune=True) as mb, azure2(azr_p, use_park=True) as mp:
            xb = np.asarray(mb.params_rwa, float)
            xp = np.asarray(mp.params_rwa, float)
            levels = level_of(mp)
            nwidth = sum(lv is not None for lv in levels)
            kinds = [p.kind for p in sorted((p for p in mp.parameters if p.free_index is not None),
                                            key=lambda p: p.free_index)]
            check(f"free parameters: {len(xp)} ({nwidth} amplitudes, kinds {sorted(set(kinds))})",
                  len(xp) == len(xb) and nwidth == 4 and (k != 0 or "cbkg" in kinds))

            # 1. the same model from the same file
            cb, cp = chi2(mb, xb), chi2(mp, xp)
            check(f"chi2 Park {cp:.10g} = Brune {cb:.10g} (rel 1e-9)", abs(cp - cb) <= 1e-9 * cb)
            rb, rp = mb.residual_jacobian(xb)[0], mp.residual_jacobian(xp)[0]
            check(f"residuals agree (rel {rel(rp, rb):.1e})", rel(rp, rb) < 1e-9)
            J = mp.park_norms(xp)
            check(f"0 < J < 1 for both levels (J = {np.array2string(J, precision=4)})",
                  bool(np.all((J > 0) & (J < 1))))
            d = rel(to_brune(mp, levels, xp), xb)
            check(f"x_Brune = x_Park / sqrt(J) for the amplitudes (rel {d:.1e})", d < 1e-12)
            # What save_fit writes: the amplitude-input channels as Brune's
            # amplitudes (the file's convention) in both modes.
            tb = np.asarray(mb.transform_rwa(xb), float)
            tp = np.asarray(mp.transform_rwa(xp), float)
            check(f"physical values (save_fit) agree, amplitude inputs as Brune's (rel {rel(tp, tb):.1e})",
                  rel(tp, tb) < 1e-12)

            # 2. the Park Jacobian against differences of its own residuals
            h = 1e-5 * (np.abs(xp) + 1.0)
            r, Jp = mp.residual_jacobian(xp)
            Jfd = residual_fd(mp, xp, h)
            d = rel(Jp, Jfd)
            check(f"residual_jacobian (Park) = central differences (rel {d:.1e})", d < 1e-4)

            # 3. the chain rule through J: dr/dx_Park = dr/dx_Brune . dx_Brune/dx_Park
            _, Jb = mb.residual_jacobian(xb)
            T = np.zeros((xp.size, xp.size))
            for i in range(xp.size):
                a = xp.copy(); a[i] += h[i]
                b = xp.copy(); b[i] -= h[i]
                T[:, i] = (to_brune(mp, levels, a) - to_brune(mp, levels, b)) / (2.0 * h[i])
            off = T - np.diag(np.diag(T))
            check("x_Brune depends on the level's other amplitudes and energy (J chain)",
                  float(np.max(np.abs(off))) > 1e-3)
            d = rel(Jp, Jb @ T)
            check(f"Park Jacobian = Brune Jacobian x dx_Brune/dx_Park (rel {d:.1e})", d < 1e-5)

            # 4. the gradient
            val, g = mp.chi2_and_grad(xp)
            check("chi2_and_grad value = chi2", abs(val - cp) <= 1e-9 * cp)
            d = rel(g, 2.0 * Jp.T @ r)
            check(f"chi2_and_grad = 2 J^T r (rel {d:.1e})", d < 1e-6)
            gfd = np.array([(chi2(mp, xp + h[i] * np.eye(xp.size)[i]) -
                             chi2(mp, xp - h[i] * np.eye(xp.size)[i])) / (2.0 * h[i])
                            for i in range(xp.size)])
            d = rel(g, gfd)
            check(f"chi2_and_grad = central differences of chi2 (rel {d:.1e})", d < 1e-4)

            if k != 0:
                continue
            # The J > 0 wall: scale level 2's amplitudes until its J < 0.
            idx = [i for i, lv in enumerate(levels) if lv == 1]
            for scale in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
                xw = xp.copy()
                xw[idx] *= scale
                Jw = mp.park_norms(xw)
                if np.any(Jw < 0):
                    break
            check(f"level 2 driven to J < 0 (scale {scale:g}, J = {np.array2string(Jw, precision=3)})",
                  bool(np.any(Jw < 0)))
            pen = float(mp.penalties(xw)["park"])
            want = float(np.sum((Jw[Jw < 0] / 1e-3) ** 2))
            check(f"penalties()['park'] = sum (J/1e-3)^2 = {want:.6g}", abs(pen - want) <= 1e-9 * want)
            data = float(np.sum(mp.segment_chi2(xw)))
            check("calculate_chi2_rwa = data chi2 + penalty", abs(chi2(mp, xw) - data - pen) <= 1e-6 * pen,
                  f"{chi2(mp, xw)} vs {data} + {pen}")
            _, gw = mp.chi2_and_grad(xw)
            hw = 1e-6 * (np.abs(xw) + 1.0)
            gfdw = np.array([(chi2(mp, xw + hw[i] * np.eye(xw.size)[i]) -
                              chi2(mp, xw - hw[i] * np.eye(xw.size)[i])) / (2.0 * hw[i])
                             for i in range(xw.size)])
            d = rel(gw, gfdw)
            check(f"gradient with the wall active = central differences (rel {d:.1e})", d < 1e-2)

print()
if failures:
    print(f"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
