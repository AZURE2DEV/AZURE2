"""sign_distance (2026-11-12 rewrite): exact path vs an independent brute force
over every gauge component, the multistart path vs the exact one, flip-set
validity, agreement with the pre-2026-11-12 enumeration, and speed at
large-model size (27 pairs).  Pure python, no AZURE2.  Run:
python3 rmfit/tests/test_sign_distance.py [path/to/old_keys.py]
"""
import itertools, os, sys, time, importlib.util
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from rmfit.keys import ParamKey, canonical_signs, flip, sign_distance

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)

def instance(rng, nlev, npair, p_amp=0.5, nfix=0, p_zero=0.05):
    keys, fixed = [], []
    for l in range(nlev):
        for p in range(1, npair + 1):
            if rng.random() < p_amp:
                keys.append(ParamKey("width", f"1/2+#{l+1}", p, 0, 0.5))
    for _ in range(nfix):
        fixed.append(ParamKey("width", f"1/2+#{rng.integers(nlev)+1}", int(rng.integers(1, npair + 1)), 9, 0.5))
    keys.append(ParamKey("norm", seg=1))
    x1 = rng.choice([-1.0, 1.0], len(keys)) * rng.uniform(0.5, 2, len(keys))
    x2 = x1 * rng.choice([-1.0, 1.0], len(keys), p=[0.3, 0.7])
    z = rng.random(len(keys)) < p_zero
    x1[z] = 0.0
    x2[z & (rng.random(len(keys)) < 0.5)] = 0.0       # zero in one or both vectors: carries no sign
    return keys, fixed, x1, x2

def brute(x1, x2, keys, fixed):
    """Every sign of every gauge component (levels and pairs joined by fixed amplitudes)."""
    parent = {}
    def find(u):
        parent.setdefault(u, u)
        while parent[u] != u:
            u = parent[u]
        return u
    for k in fixed:
        a, b = find(("L", k.level)), find(("P", k.member, k.pair))
        if a != b: parent[a] = b
    amps = [(i, k) for i, k in enumerate(keys) if k.kind == "width" and x1[i] != 0 and x2[i] != 0]
    ed = [(find(("L", k.level)), find(("P", k.member, k.pair)), 1 if (x1[i] > 0) == (x2[i] > 0) else -1) for i, k in amps]
    nodes = sorted({u for a, b, _ in ed for u in (a, b)}, key=str)
    best = None
    for s in itertools.product((1, -1), repeat=len(nodes)):
        g = dict(zip(nodes, s))
        c = sum(1 for a, b, d in ed if d * g[a] * g[b] < 0)
        best = c if best is None else min(best, c)
    return best or 0

def valid(x1, x2, keys, fixed, diff):
    """Flipping diff in x1 must give x2's pattern up to a gauge: sign_distance == 0 (exact)
    and, since the canonical_signs fix of 2026-11-12, equal canonical strings on the common support."""
    y = flip(x1, keys, diff)
    both = (np.asarray(x1) != 0) & (np.asarray(x2) != 0)
    return (sign_distance(y, x2, keys, fixed)[0] == 0
            and canonical_signs(np.where(both, y, 0), keys, fixed) == canonical_signs(np.where(both, x2, 0), keys, fixed))

rng = np.random.default_rng(12)
# 1. exact path == brute force, with ties from fixed amplitudes (incl. several levels on one amplitude-less pair)
bad = []
for t in range(300):
    keys, fixed, x1, x2 = instance(rng, int(rng.integers(1, 6)), int(rng.integers(1, 5)), nfix=int(rng.integers(0, 4)))
    if len({k.level for k in keys if k.level} | {(k.pair) for k in keys if k.pair}) + len(fixed) > 12: continue
    d, diff = sign_distance(x1, x2, keys, fixed)
    b = brute(x1, x2, keys, fixed)
    y = flip(x1, keys, diff)
    if d != b or brute(y, x2, keys, fixed) != 0: bad.append((t, d, b))
check("exact path == brute force over all components (300 random)", not bad, str(bad[:5]))
# 2. multistart path (exact_max=0) == exact path, never below, flip set valid
gaps, inval, n = [], 0, 0
for t in range(200):
    keys, fixed, x1, x2 = instance(rng, int(rng.integers(5, 40)), int(rng.integers(4, 16)), p_amp=0.4, nfix=int(rng.integers(0, 6)))
    de, _ = sign_distance(x1, x2, keys, fixed)
    dh, diff = sign_distance(x1, x2, keys, fixed, exact_max=0)
    n += 1; gaps.append(dh - de); inval += not valid(x1, x2, keys, fixed, diff)
gaps = np.array(gaps)
check(f"multistart never below exact ({n} instances, up to 15 pairs)", gaps.min() >= 0, str(gaps.min()))
check(f"multistart == exact on {np.sum(gaps == 0)}/{n}", np.sum(gaps == 0) >= 0.98 * n, f"gaps {gaps[gaps>0]}")
check("multistart flip sets valid", inval == 0, str(inval))
# 3. distance 0 found exactly by the multistart path (a pure gauge transform of x1)
z = 0
for t in range(100):
    keys, fixed, x1, _ = instance(rng, 40, 27, p_amp=0.5, nfix=int(rng.integers(0, 8)))
    lv = sorted({k.level for k in keys if k.level}); pr = sorted({k.pair for k in keys if k.pair})
    gl = dict(zip(lv, rng.choice([-1, 1], len(lv)))); gp = dict(zip(pr, rng.choice([-1, 1], len(pr))))
    tied = {(k.level, k.pair) for k in fixed}
    x2 = x1.copy()
    for i, k in enumerate(keys):
        if k.kind == "width": x2[i] *= gl[k.level] * gp[k.pair]
    # a gauge transform must leave fixed amplitudes unchanged: only test untied cases or tie-consistent ones
    if any(gl[l] * gp[p] < 0 for l, p in tied): continue
    z += sign_distance(x1, x2, keys, fixed)[0] != 0
check("27 pairs: a gauge transform has distance 0", z == 0, str(z))
# 4. agreement with the pre-2026-11-12 enumeration (no fixed ties, so its multi-level tie bug cannot enter)
if len(sys.argv) > 1:
    spec = importlib.util.spec_from_file_location("keys_old", sys.argv[1]); old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    dis = []
    for t in range(200):
        keys, fixed, x1, x2 = instance(rng, int(rng.integers(1, 30)), int(rng.integers(1, 11)), nfix=0)
        a = sign_distance(x1, x2, keys)[0]; b = old.sign_distance(x1, x2, keys)[0]
        if a != b: dis.append((t, a, b))
    check("== old enumeration on 200 untied instances (<= 10 pairs)", not dis, str(dis[:5]))
# 5. speed at a large model's size: 27 pairs, 60 levels, ~600 amplitudes, ties
keys, fixed, x1, x2 = instance(rng, 60, 27, p_amp=0.37, nfix=12)
t0 = time.time(); d, diff = sign_distance(x1, x2, keys, fixed); dt = time.time() - t0
check(f"27 pairs / {sum(k.kind=='width' for k in keys)} amplitudes: {dt:.2f} s (distance {d})", dt < 30 and valid(x1, x2, keys, fixed, diff))
keys, fixed, x1, x2 = instance(rng, 60, 20, p_amp=0.5, nfix=0)
t0 = time.time(); d, diff = sign_distance(x1, x2, keys, fixed); dt = time.time() - t0
check(f"20 pairs exact path: {dt:.2f} s", dt < 30)
# 6. canonical_signs (2026-11-12 fix): invariant under every gauge transform that keeps the fixed amplitudes,
#    and equal strings <=> distance 0 (different orbits never share a string)
def gauge(rng, x, keys, fixed):
    lv = sorted({k.level for k in keys if k.level}); pr = sorted({(k.member, k.pair) for k in keys if k.pair})
    for _ in range(50):
        gl = dict(zip(lv, rng.choice([-1, 1], len(lv)))); gp = dict(zip(pr, rng.choice([-1, 1], len(pr))))
        if all(gl.get(k.level, 1) * gp.get((k.member, k.pair), 1) > 0 for k in fixed):
            return np.array([x[i] * gl[k.level] * gp[(k.member, k.pair)] if k.kind == "width" else x[i]
                             for i, k in enumerate(keys)])
    return x.copy()
ninv = nsep = nsame = 0
for t in range(1000):
    keys, fixed, x1, x2 = instance(rng, int(rng.integers(1, 30)), int(rng.integers(1, 10)), p_amp=0.5, nfix=int(rng.integers(0, 5)))
    ninv += canonical_signs(x1, keys, fixed) != canonical_signs(gauge(rng, x1, keys, fixed), keys, fixed)
    x2 = np.where(x1 != 0, x2, 0)                         # same support
    same = sign_distance(x1, x2, keys, fixed)[0] == 0
    nsame += same
    nsep += same != (canonical_signs(x1, keys, fixed) == canonical_signs(x2, keys, fixed))
check("canonical_signs invariant under 1000 random gauge transforms (ties, zeros)", ninv == 0, str(ninv))
check(f"canonical_signs equal <=> sign distance 0 (1000 pairs, {nsame} gauge-equivalent)", nsep == 0, str(nsep))
print("test_sign_distance:", "FAILED " + ", ".join(failures) if failures else "ok")
sys.exit(1 if failures else 0)
