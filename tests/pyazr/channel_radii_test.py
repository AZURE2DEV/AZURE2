#!/usr/bin/env python3
"""Per-channel radii through pyazr (AZURE2 dev 2026-11-10): the <channelRadii> block in AzrModel and the
session's set_channel_radius(..., L=, S=, J=) / channel_radius_table().

Model: tests/channel_space/channel_space.azr -- p+12C (pair 1), p+12C(2+) (pair 2, channel spins 3/2 and 5/2, so a
J-group holds several channels of one pair and l) and a gamma pair (3), every radius 5.0 fm.

  1. AzrModel: set_channel_radius(2, 4.0, L=2, S=3/2) writes "<channelRadii> / 2 2 4 s=3/2 / </channelRadii>" right
     after </levels>; a second call for the same (pair, L, S, J) replaces the line; entries read back; the pair form
     (no L) still rewrites column 28; S without L, a photon pair and a non-positive radius are refused;
     remove_channel_radius drops the line and the empty block; channel_radius_of follows the engine's precedence.
  2. A session on the file with the block: channel_radius_table() gives 4.0 fm to exactly the channels (pair 2, l = 2,
     s = 3/2) and 5.0 to all others; chi2 differs from the plain file (the radius matters).
  3. A session on the plain file with set_channel_radius(2, 4.0, L=2, S=1.5): the same radii and the same chi2 (1e-12),
     and it survives a mode switch and a later pair radius change (pair 2 -> 6.0 moves only the channels without an
     entry).
  4. Precedence: "2 2 4.0" and "2 2 4.4 s=5/2 J=1/2" -> (2, 2, 5/2, 1/2) at 4.4, (2, 2, 3/2, 1/2) at 4.0, l != 2 at 5.0,
     in the engine and in AzrModel.channel_radius_of alike.
  5. transform.levels_from_azr and the level scheme carry the per-channel radius.
  6. An s entry and a J entry matching one channel without an l+s+J entry: RuntimeError naming the engine's refusal.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/channel_radii_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def raises(exc, fn):
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import azure2, AzrModel, ChannelRadiusEntry
        from pyazr.transform import levels_from_azr
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 0

    work = tempfile.mkdtemp(prefix="pyazr_channel_radii_")
    try:
        src = os.path.join(ROOT, "tests", "channel_space")

        def project(name, edit=None):
            d = os.path.join(work, name)
            for sub in ("output", "checks"):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
            shutil.copytree(os.path.join(src, "data"), os.path.join(d, "data"))
            mdl = AzrModel.from_file(os.path.join(src, "channel_space.azr"))
            if edit:
                edit(mdl)
            return mdl.write(os.path.join(d, "m.azr")), d

        def chi(s):
            return float(np.sum(s.calculate_chi2_rwa(np.asarray(s.params_rwa, float))))

        def radii(s):
            return {(r.pair, r.L, r.S, r.J): r.radius for r in s.channel_radius_table()}

        print("1. AzrModel")
        mdl = AzrModel.from_file(os.path.join(src, "channel_space.azr"))
        check("no entries in the plain file", mdl.channel_radius_entries() == [])
        mdl.set_channel_radius(2, 4.0, L=2, S=1.5)
        text = mdl.text() if hasattr(mdl, "text") else open(mdl.write(os.path.join(work, "t1.azr"))).read()
        check("block written after </levels>", "</levels>\n<channelRadii>\n2 2 4 s=3/2\n</channelRadii>" in text,
              text[text.find("</levels>"):text.find("</levels>") + 80])
        mdl.set_channel_radius(2, 4.25, L=2, S="3/2")
        check("same (pair, L, S, J) replaces the line",
              mdl.channel_radius_entries() == [ChannelRadiusEntry(2, 2, 4.25, 1.5, None)], f"{mdl.channel_radius_entries()}")
        mdl.set_channel_radius(2, 4.4, L=2, S=2.5, J=0.5)
        check("two entries", len(mdl.channel_radius_entries()) == 2, f"{mdl.channel_radius_entries()}")
        check("channel_radius_of: l+s+J entry", mdl.channel_radius_of(2, 2, 2.5, 0.5) == 4.4)
        check("channel_radius_of: l+s entry", mdl.channel_radius_of(2, 2, 1.5, 0.5) == 4.25)
        check("channel_radius_of: no entry -> pair radius", mdl.channel_radius_of(2, 1, 1.5, 1.5) == 5.0)
        n = mdl.set_channel_radius(1, 4.8)
        check("pair form still rewrites column 28", n > 0 and mdl.channel_radii()[1] == 4.8, f"{n} {mdl.channel_radii()}")
        check("S without L refused", raises(ValueError, lambda: mdl.set_channel_radius(2, 4.0, S=1.5)))
        check("photon pair refused", raises(KeyError, lambda: mdl.set_channel_radius(3, 4.0, L=1)))
        check("non-positive radius refused", raises(ValueError, lambda: mdl.set_channel_radius(2, 0.0, L=1)))
        mdl.remove_channel_radius(2, 2, S=2.5, J=0.5)
        mdl.remove_channel_radius(2, 2, S=1.5)
        t2 = open(mdl.write(os.path.join(work, "t2.azr"))).read()
        check("removing every entry drops the block", "<channelRadii>" not in t2)
        check("removing a missing entry raises KeyError", raises(KeyError, lambda: mdl.remove_channel_radius(2, 2)))

        print("\n2. session on a file with the block")
        plain, plain_dir = project("plain")
        blk, blk_dir = project("block", lambda m: m.set_channel_radius(2, 4.0, L=2, S=1.5))
        with azure2(plain, cwd=plain_dir) as p:
            chi_plain = chi(p)
            r_plain = radii(p)
        check("plain file: every particle channel at 5.0 fm",
              all(v == 5.0 for k, v in r_plain.items() if k[0] in (1, 2)), f"{r_plain}")
        with azure2(blk, cwd=blk_dir) as b:
            chi_blk = chi(b)
            r_blk = radii(b)
            scheme = b.level_scheme
        moved = {k for k, v in r_blk.items() if v == 4.0}
        want = {k for k in r_blk if k[0] == 2 and k[1] == 2 and abs(k[2] - 1.5) < 1e-9}
        check("4.0 fm on exactly (pair 2, l = 2, s = 3/2)", moved == want and len(want) > 0, f"{sorted(moved)} vs {sorted(want)}")
        check("every other particle channel at 5.0 fm",
              all(v == 5.0 for k, v in r_blk.items() if k not in want and k[0] in (1, 2)), f"{r_blk}")
        check("the radius matters (chi2 differs from the plain file)", rel(chi_blk, chi_plain) > 1e-6,
              f"{chi_blk} vs {chi_plain}")

        print("\n3. the same through the session")
        with azure2(plain, cwd=plain_dir) as s:
            got = s.set_channel_radius(2, 4.0, L=2, S=1.5)
            check("set_channel_radius returns the matched channels", len(got) == len(want) and all(r.radius == 4.0 for r in got),
                  f"{got}")
            check("same radii as the file's block", radii(s) == r_blk, f"{radii(s)}")
            c = chi(s)
            check("same chi2 as the file's block (1e-12)", rel(c, chi_blk) <= 1e-12, f"{c} vs {chi_blk}")
            s.data_mode()
            s.configure()
            check("kept through a mode switch", radii(s) == r_blk)
            s.set_channel_radius(2, 6.0)
            r6 = radii(s)
            check("a later pair radius moves only the channels without an entry",
                  all((v == 4.0) if k in want else (v == 6.0) for k, v in r6.items() if k[0] == 2) and
                  all(v == 5.0 for k, v in r6.items() if k[0] == 1), f"{r6}")

        print("\n4. precedence")

        def two(m):
            m.set_channel_radius(2, 4.0, L=2)
            m.set_channel_radius(2, 4.4, L=2, S=2.5, J=0.5)
        prec, prec_dir = project("prec", two)
        pm = AzrModel.from_file(prec)
        with azure2(prec, cwd=prec_dir) as q:
            rq = radii(q)
        check("engine: (2, 2, 5/2, 1/2) at 4.4", rq.get((2, 2, 2.5, 0.5)) == 4.4, f"{rq}")
        check("engine: (2, 2, 3/2, 1/2) at 4.0", rq.get((2, 2, 1.5, 0.5)) == 4.0, f"{rq}")
        check("engine: pair 2 l != 2 at 5.0", all(v == 5.0 for k, v in rq.items() if k[0] == 2 and k[1] != 2), f"{rq}")
        check("AzrModel.channel_radius_of agrees with the engine on every channel",
              all(abs(pm.channel_radius_of(*k) - v) < 1e-12 for k, v in rq.items() if k[0] in (1, 2)))

        print("\n5. transform and level scheme")
        lv = levels_from_azr(AzrModel.from_file(blk))
        rr = {(c.pair, c.L, c.S, l.J): c.channel_radius for l in lv for c in l.channels if c.pair in (1, 2)}
        check("levels_from_azr: per-channel radius", all((v == 4.0) == (k in want) for k, v in rr.items()), f"{rr}")
        sch = {(c.pair, c.L, c.S, lvl.J): c.radius for lvl in scheme.levels for c in lvl.channels if c.pair in (1, 2)}
        check("level scheme: per-channel radius", all((v == 4.0) == (k in want) for k, v in sch.items()), f"{sch}")

        print("\n6. an s entry and a J entry on one channel")
        amb_dir = os.path.join(work, "plain")
        with azure2(plain, cwd=amb_dir) as a:
            a.set_channel_radius(2, 4.1, L=2, S=1.5)
            try:
                a.set_channel_radius(2, 4.2, L=2, J=0.5)
                check("RuntimeError on an ambiguous channel", False, "no exception")
            except RuntimeError as e:
                check("RuntimeError on an ambiguous channel", "could not rebuild" in str(e), str(e))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
