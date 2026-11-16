#!/usr/bin/env python3
"""rmfit level labels in a model with several compound nuclei (2026-10-30).

The engine has one J-group per (J^pi, compound nucleus), each numbered from 1, so in a natural-target project with two
odd-A nuclei "1/2+#1" names two levels (21Ne and 23Ne in a natural-neon target).  Before this fix rmfit's labels carried no
nucleus: SingleSession.level_energy kept ONE energy per label (default_bounds then capped widths, fixed background
poles and set width floors from the other nucleus's level) and two free energies shared the key E[1/2+#1].
Now a multi-nucleus model labels every level "<J^pi>#<n>@<nucleus>"; a single-nucleus model is unchanged.

Fixture tests/fixtures/multicn_labels/mcn.azr: p+12C -> 13N and p+14C -> 15N, levels 1/2+, 3/2-, 5/2+ in each
(13N at 2.365 / 3.502 / 3.547 MeV, 15N at 11.0074 / 11.5074 / 11.9074 MeV), energies and widths free.
  (a) label helpers (pure python): level_group, ordinal, nucleus_tag, jpi;
  (b) SingleSession: six distinct energy keys, level_energy has six entries with each level's own energy;
  (c) default_bounds with bg_energy = 10 MeV: the 15N energies are fixed as background poles, the 13N ones are not,
      and the widths take the background / physical caps of their own level;
  (d) Sharded (2 shards, one nucleus's data each): the same keys, objective = SingleSession's;
  (e) the removal-variant helpers resolve a level of the right nucleus (campaign._level_energy,
      structure.deactivate_level_variant(nucleus=...));
  (f) single-nucleus model (AZURE2 tests/13N): keys and level_energy byte-identical to the pre-fix labels.
    AZURE2_ROOT=<tree> python3 rmfit/tests/test_level_labels_multicn.py      (about a minute)"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
FIX = os.path.join(HERE, "fixtures", "multicn_labels")
failures = []


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        failures.append(what)


def main():
    from rmfit import AZURE2_ROOT  # noqa: F401  (puts the AZURE2 tree's pyazr on the path)
    from rmfit.keys import ParamKey, jpi_string
    from rmfit.engine import SingleSession, Sharded
    from rmfit.local import default_bounds, BoundsPolicy

    # (a) label helpers (guarded so that the rest also runs on a pre-fix rmfit, which has none of them)
    try:
        from rmfit.keys import level_group
        k = ParamKey("energy", level="1/2-#2@21Ne")
        check((k.jpi, k.ordinal, k.nucleus_tag) == ("1/2-", 2, "21Ne"),
              "(a) 1/2-#2@21Ne -> jpi 1/2-, ordinal 2, nucleus 21Ne")
        p = ParamKey("energy", level="5/2+#3")
        check((p.jpi, p.ordinal, p.nucleus_tag) == ("5/2+", 3, ""), "(a) a plain label parses as before")
        check(level_group("1/2-#2@21Ne") == "1/2-@21Ne" and level_group("5/2+#3") == "5/2+"
              and level_group("22Ne:1/2+#3") == "22Ne:1/2+", "(a) level_group keeps the nucleus and a member prefix")
        m = ParamKey("energy", level="22Ne:1/2+#3")
        check((m.jpi, m.ordinal) == ("1/2+", 3), "(a) a member-tagged label (rmfit.mixture) parses as before")
    except (ImportError, AttributeError, ValueError) as err:
        check(False, f"(a) label helpers: {err.__class__.__name__}: {err}")

    work = tempfile.mkdtemp(prefix="rmfit_mcn_")
    try:
        fix = os.path.join(work, "mcn")
        shutil.copytree(FIX, fix)
        os.makedirs(os.path.join(fix, "output"), exist_ok=True)
        os.makedirs(os.path.join(fix, "checks"), exist_ok=True)
        azr = os.path.join(fix, "mcn.azr")
        # absolute data paths: the shard files live in another directory
        txt = open(azr).read()
        for f in ("data_13n.dat", "data_15n.dat"):
            txt = txt.replace(" " + f, " " + os.path.join(fix, f))
        open(azr, "w").write(txt)
        flags = dict(use_brune=True)
        S = SingleSession(azr, cwd=fix, flags=flags, name="mcn")
        try:
            ekeys = [str(k) for k in S.keys if k.kind == "energy"]
            want = {f"E[{j}#1@{n}]" for j in ("1/2+", "3/2-", "5/2+") for n in ("13N", "15N")}
            check(set(ekeys) == want and len(ekeys) == 6, f"(b) six distinct energy keys {sorted(ekeys)}")
            check(len({str(k) for k in S.keys}) == len(S.keys), "(b) no duplicate key")
            le = S.level_energy
            exp = {"1/2+#1@13N": 2.365, "3/2-#1@13N": 3.502, "5/2+#1@13N": 3.547,
                   "1/2+#1@15N": 11.0074, "3/2-#1@15N": 11.5074, "5/2+#1@15N": 11.9074}
            check(set(le) == set(exp) and all(abs(le[l] - e) < 1e-6 for l, e in exp.items()),
                  f"(b) level_energy has each level's own energy: {le}")
            wl = {k.level for k in S.keys if k.kind == "width"}
            check(wl == set(exp), "(b) width keys carry the same labels")

            # (c) bounds from each level's own energy
            pol = BoundsPolicy(bg_energy=10.0)
            b = default_bounds(S, S.x0, pol)
            fx = {str(S.keys[j]) for j in range(len(S.keys)) if b.fixed[j]}
            check(fx == {f"E[{j}#1@15N]" for j in ("1/2+", "3/2-", "5/2+")},
                  f"(c) bg_energy 10 MeV fixes the 15N energies only: {sorted(fx)}")
            gW = {c: float(np.sqrt(w)) for c, w in S.wigner_limits.items()}     # gamma_W^2 -> gamma_W
            ok = True
            for j, kk in enumerate(S.keys):
                if kk.kind != "width" or j not in gW:
                    continue
                cap = pol.cap_for(exp[kk.level])
                ok &= abs(b.ub[j] - np.sqrt(cap) * gW[j]) <= 1e-9 * b.ub[j]
            check(ok, f"(c) every width cap = sqrt(cap(E of its own level)) * gamma_W "
                      f"(theta^2 caps {pol.theta2_cap} / background {pol.bg_theta2_cap})")

            # (d) Sharded: one shard per nucleus's data, same keys and objective
            H = Sharded(azr, os.path.join(work, "sh"), nshards=2, flags=flags)
            try:
                check([str(k) for k in H.keys] == [str(k) for k in S.keys], "(d) Sharded keys = SingleSession keys")
                check(H.level_energy == S.level_energy, "(d) Sharded level_energy = SingleSession's")
                rng = np.random.default_rng(3)
                x1 = S.x0 * (1 + 0.01 * rng.standard_normal(S.x0.size))
                oS, oH = S.objective(x1), H.objective(x1)
                check(abs(oS - oH) <= 1e-9 * oS, f"(d) objective Single {oS:.8f} = Sharded {oH:.8f}")
            finally:
                H.close()
        finally:
            S.close()

        # (e) the removal-variant helpers find the level of the right nucleus
        from rmfit.campaign import Campaign
        from rmfit.structure import deactivate_level_variant
        from rmfit.azr3 import Model as AzrModel
        fake = Campaign.__new__(Campaign)
        fake.base_azr = azr
        e13, e15 = fake._level_energy("3/2-#1@13N"), fake._level_energy("3/2-#1@15N")
        check(e13 == 3.502 and e15 == 11.5074, f"(e) campaign._level_energy: 13N {e13}, 15N {e15}")
        out = os.path.join(fix, "rm.azr")
        info = deactivate_level_variant(azr, out, "3/2-", 11.5074, nucleus="15N")
        mdl = AzrModel.from_file(out)
        off = [(lv.jpi, lv.energy) for lv in mdl.levels if all(c.gamma == 0.0 for c in lv.channels)]
        check(off == [("3/2-", 11.5074)], f"(e) deactivate_level_variant(nucleus='15N') switches off {off}")
        try:
            deactivate_level_variant(azr, out, "3/2-", 11.5074, nucleus="13N")
            check(False, "(e) a level of the other nucleus is not matched")
        except KeyError:
            check(True, "(e) a level of the other nucleus is not matched")

        # (f) one nucleus: labels as before the fix
        root = os.environ.get("AZURE2_ROOT", os.path.expanduser("~/AZURE2"))
        one = os.path.join(work, "13N")
        shutil.copytree(os.path.join(root, "tests", "13N"), one)
        os.makedirs(os.path.join(one, "output"), exist_ok=True)
        os.makedirs(os.path.join(one, "checks"), exist_ok=True)
        from pyazr import azure2
        S1 = SingleSession(os.path.join(one, "13N.azr"), cwd=one, flags=flags, name="one")
        E1 = azure2(os.path.join(one, "13N.azr"), cwd=one, **flags)
        try:
            old_level = lambda p: f"{jpi_string(p.J, p.parity)}#{int(p.level)}"     # the pre-fix label
            old = {p.free_index: old_level(p) for p in E1.parameters
                   if p.kind in ("energy", "width") and not p.fixed and p.free_index is not None}
            check(all(S1.keys[j].level == lab for j, lab in old.items()) and old,
                  f"(f) single nucleus: the {len(old)} level keys are the pre-fix labels")
            old_le = {old_level(k): (0.0 if k.energy is None else float(k.energy))
                      for k in E1.parameters.by_physical_level()}
            check(S1.level_energy == old_le, "(f) single nucleus: level_energy identical to the pre-fix table")
            check(not any("@" in l for l in S1.level_energy), "(f) single nucleus: no '@' anywhere")
        finally:
            S1.close()
            E1.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("FAIL" if failures else "PASS", "test_level_labels_multicn", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
