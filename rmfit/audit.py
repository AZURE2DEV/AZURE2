"""``rmfit audit``: pre-flight checks of a model file and its data before (or without) a campaign.  No engine session,
nothing is changed; seconds on the login node.  Writes AUDIT.md next to campaign.json (or next to the .azr).

Why (DESIGN.md, findings 2.5 and 3.2-3.5): input errors were found by hand audits or by
the evaluator, never by the fit -- an EXFOR FLAG column read as errors (for months), normErr written as a
fraction where AZURE2 reads percent (penalties 100x too tight), hundreds of free channels exactly 0 and several
nominal norms != 1, legacy .sav names silently ignored (thousands of files), absolute output paths writing into the
archive, <targetInt> lines folding the test segments.  Every check here is one of those, generalized.

Checks (code -> severity):
  engine-pin            campaign.json["engine"] vs the engine this process imports            HIGH on mismatch
  output-path           <config> output directory absolute and outside the model's directory    HIGH
  pair-mass             pair masses vs AME2020 (pyazr.masses), integer masses, S vs table        MEDIUM / HIGH (> 1 keV)
  zero-width            free (unfixed) channel widths exactly 0 on active levels                 MEDIUM (count)
  level-merge           two active levels of one J^pi read as one (same ID within 1 keV, or same energy)  HIGH
  level-close           two active levels of one J^pi within 1 keV, different IDs (two levels since dev e239a12)  MEDIUM
  norm-nominal          active segment with nominal norm != 1                                   MEDIUM
  norm-unpenalized      vary_norm with normError 0                                              MEDIUM
  norm-error-units      0 < normError < 1 (percent expected: a fraction gives 100x the weight)   HIGH
  shift-unpenalized     vary_shift with shift error 0                                           MEDIUM
  data-file             missing file, parse problems, non-finite rows                           HIGH
  data-values           value <= 0 (not A_y), error <= 0, error > |value|                        MEDIUM
  data-duplicates       duplicate (E, angle) rows inside one file                               LOW
  data-rows-unused      rows of a file outside every active window                              INFO
  data-below-threshold  rows below the exit channel's lab threshold                             MEDIUM
  data-spikes           isolated point > 1.8x / < 0.55x its neighbours at > 5 sigma             LOW
  data-crossdup         rows of one file repeated in another (same E, value, error)             MEDIUM
  data-error-pattern    error / value takes <= 5 distinct exact percent values (FLAG column?)    HIGH
  data-error-tiny       rows with error < 0.1 % of the value (a typo outweighs thousands of points) MEDIUM / HIGH
  targetint-test        active <targetInt> line naming a <segmentsTest> key (folded extrapolation) HIGH
  legacy-sav            j=.._la=.. names in a .sav/.par next to the model (ignored by AZURE2)   HIGH
  intec-stale           output/intEC.* older than the .azr                                      LOW
  frames                A_y (isDiff 7) or differential data with angles > 180 / < 0             MEDIUM
The three idle-time reviews (input file, data, literature) are subagent PROMPTS, not checks here:
`rmfit <campaign> audit --kind input|data|literature|all` (audit_subagent.py, templates/audit_*.md).
"""
import glob
import os
import re
import time
from collections import defaultdict
from typing import List, Optional

import numpy as np

from .spec import Finding, write_findings_md, pair_lines

LEGACY_RE = re.compile(r"^j=(\d+)_la=(\d+)_(energy|ch=\d+_rwa)$")


def _block(text, tag):
    m = re.search(r"<%s>(.*?)</%s>" % (tag, tag), text, re.S)
    return m.group(1) if m else ""


def read_points(path):
    rows, probs = [], []
    with open(path) as fh:
        for k, line in enumerate(fh, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            t = s.split()
            if len(t) < 4:
                probs.append(f"line {k}: {len(t)} columns")
                continue
            try:
                rows.append([float(x) for x in t[:4]])
            except ValueError:
                probs.append(f"line {k}: non-numeric")
    return (np.array(rows) if rows else np.zeros((0, 4))), probs


def in_window(d, seg):
    m = (d[:, 0] >= seg.energy_min - 1e-9) & (d[:, 0] <= seg.energy_max + 1e-9)
    if seg.observable in ("differential", "analyzing-power") or seg.angle_max > seg.angle_min:
        if seg.observable not in ("angle-integrated", "total-capture", "total"):
            m &= (d[:, 1] >= seg.angle_min - 1e-9) & (d[:, 1] <= seg.angle_max + 1e-9)
    return m


def audit_model(azr_path: str, campaign_dir: Optional[str] = None, campaign_cfg: Optional[dict] = None,
                run_target_energies: bool = True) -> List[Finding]:
    azr_path = os.path.abspath(azr_path)
    base = os.path.dirname(azr_path)
    out: List[Finding] = []
    text = open(azr_path).read()

    # -- engine pin ------------------------------------------------------------------------------
    if campaign_cfg is not None:
        from . import enginepin
        rec = campaign_cfg.get("engine")
        if rec and rec.get("so_md5"):
            cur = enginepin.fingerprint()
            if cur["so_md5"] != rec["so_md5"]:
                out.append(Finding("HIGH", "engine-pin", f"campaign pinned {enginepin.short(rec)}; this process imports {enginepin.short(cur)}"))
            else:
                out.append(Finding("INFO", "engine-pin", f"pinned {enginepin.short(rec)}"))
        else:
            out.append(Finding("LOW", "engine-pin", "no engine pin in campaign.json (campaign made before 2026-11-07)"))

    # -- output path -------------------------------------------------------------------------------
    cfg = _block(text, "config")
    outdir = None
    for line in cfg.splitlines():          # "<path>   #Full Path to Output Directory" (pyazr set_output_dir writes the same)
        if line.rstrip().endswith("#Full Path to Output Directory"):
            outdir = line[:line.index("#")].strip()
            break
    if outdir:
        od = outdir if os.path.isabs(outdir) else os.path.join(base, outdir)
        home = os.path.normpath(campaign_dir or base)
        if os.path.isabs(outdir) and not os.path.normpath(od).startswith(home + os.sep) and os.path.normpath(od) != home:
            out.append(Finding("HIGH", "output-path", f"absolute output path {outdir} outside {home}: a run here writes into another directory "
                               "(a scratch run once overwrote an archived fit's output files)"))
        elif not os.path.isdir(od):
            out.append(Finding("MEDIUM", "output-path", f"output directory {od} does not exist (AZURE2 needs it before a run)"))

    # -- pairs, levels ---------------------------------------------------------------------------------
    try:
        from rmfit.azr3 import Model as AzrModel
        model = AzrModel.from_file(azr_path)
    except Exception as err:
        out.append(Finding("HIGH", "model-parse", f"pyazr could not read the model: {err}"))
        model = None
    if model is not None:
        try:
            for row in model.mass_check():
                if row.ptype != 0 and row.light is None and row.heavy is None:
                    continue
                for name, mval, tval in (("M1", row.M1, row.M1_table), ("M2", row.M2, row.M2_table)):
                    if tval is not None and abs(mval - round(mval)) < 1e-9 and abs(mval - tval) > 1e-6:
                        out.append(Finding("MEDIUM", "pair-mass", f"pair {row.pair} {name} = {mval:g} is a mass NUMBER, not an atomic mass "
                                           f"(table {tval:.6f}; Golden rule 1)"))
                    elif tval is not None and abs(mval - tval) > 5e-4:
                        out.append(Finding("MEDIUM", "pair-mass", f"pair {row.pair} {name} = {mval!r} u, table {tval!r} u"))
                if row.dS_kev is not None and abs(row.dS_kev) > 1.0:
                    out.append(Finding("HIGH" if abs(row.dS_kev) > 10 else "MEDIUM", "pair-separation",
                                       f"pair {row.pair}: S {row.S:.6f} MeV vs table {row.S_table:.6f} ({row.dS_kev:+.1f} keV)"))
        except Exception as err:
            out.append(Finding("INFO", "pair-mass", f"mass check skipped ({err})"))
        zero_free, nlev, nchan = [], 0, 0
        byjpi = defaultdict(list)
        for lv in model.levels:
            if not lv.active:
                continue
            nlev += 1
            byjpi[lv.jpi].append((float(lv.energy), lv.level_id))
            for c in lv.channels:
                if not c.active:
                    continue
                nchan += 1
                if (not c.channel_fixed) and float(c.gamma) == 0.0:
                    zero_free.append(f"{lv.jpi}@{float(lv.energy):.4f} p{c.pair} L{c.L} S{c.S}")
        if zero_free:
            out.append(Finding("MEDIUM", "zero-width", f"{len(zero_free)} of {nchan} active channels are free with width exactly 0: "
                               "a width of 0 is FIXED whatever its flag (Golden rule 3)", "; ".join(zero_free[:8]) + (" ..." if len(zero_free) > 8 else "")))
        for jpi, es in byjpi.items():
            es = sorted(es)
            for a, b in zip(es, es[1:]):
                if abs(b[0] - a[0]) < 1e-3:
                    if a[1] == b[1] or a[0] == b[0]:
                        out.append(Finding("HIGH", "level-merge", f"two active {jpi} levels read as one ({a[0]:.5f} ID {a[1]}, {b[0]:.5f} ID {b[1]})"))
                    else:
                        out.append(Finding("MEDIUM", "level-close", f"two active {jpi} levels within 1 keV ({a[0]:.5f} ID {a[1]}, {b[0]:.5f} ID {b[1]}): "
                                           "two levels since AZURE2 dev e239a12 (older builds merge them); a gap far below the widths is a merge candidate for the evaluator"))
        out.append(Finding("INFO", "model", f"{nlev} active levels, {nchan} active channels, {len(pair_lines(model))} pairs"))

    # -- segments and data ------------------------------------------------------------------------------
    try:
        from pyazr.datasets import SegmentSet, TestSegmentSet
        segs = SegmentSet.from_file(azr_path)
    except Exception as err:
        out.append(Finding("HIGH", "segments-parse", f"segment table unreadable: {err}"))
        segs = []
    tests = []
    try:
        tests = TestSegmentSet.from_file(azr_path)
    except Exception:
        pass
    thresholds = _lab_thresholds(model) if model is not None else {}
    byfile = defaultdict(list)
    for s in segs:
        byfile[s.data_file].append(s)
        if not s.active:
            continue
        cov = bool(s.covariance)        # a covariance file may carry the norm prior (t0 method): LOW, for the evaluator to confirm
        if abs(s.norm - 1.0) > 1e-9:
            out.append(Finding("LOW" if cov else "MEDIUM", "norm-nominal", f"seg {s.key} {os.path.basename(s.data_file)}: nominal norm {s.norm} != 1 "
                               + ("(covariance file present: the prior may live there, t0 method -- confirm)" if cov else
                                  "(a fitted value used as the nominal biases the penalty)")))
        if s.vary_norm and s.norm_error == 0:
            out.append(Finding("LOW" if cov else "MEDIUM", "norm-unpenalized", f"seg {s.key} {os.path.basename(s.data_file)}: free norm with no "
                               "systematic error" + (" (covariance file present: cov_free_norm_ok?)" if cov else "")))
        if 0 < s.norm_error < 1.0:
            out.append(Finding("HIGH", "norm-error-units", f"seg {s.key} {os.path.basename(s.data_file)}: normError {s.norm_error} -- AZURE2 reads PERCENT; "
                               "a fraction gives a penalty 100x too tight"))
        if s.vary_shift and s.energy_shift_error == 0:
            out.append(Finding("MEDIUM", "shift-unpenalized", f"seg {s.key} {os.path.basename(s.data_file)}: free energy shift with no error"))
        if s.observable == "differential" and (s.angle_min < -1 or s.angle_max > 180.5):
            out.append(Finding("LOW", "frames", f"seg {s.key}: angle window [{s.angle_min}, {s.angle_max}] outside 0-180"))
    active_files = sorted({s.data_file for s in segs if s.active})
    keysets = {}
    for f in active_files:
        path = f if os.path.isabs(f) else os.path.join(base, f)
        name = os.path.basename(f)
        if not os.path.exists(path):
            out.append(Finding("HIGH", "data-file", f"{name}: file not found ({path})"))
            continue
        d, probs = read_points(path)
        if probs:
            out.append(Finding("HIGH", "data-file", f"{name}: {len(probs)} unparsable lines", "; ".join(probs[:4])))
        if len(d) == 0:
            out.append(Finding("HIGH", "data-file", f"{name}: no data rows"))
            continue
        bad = ~np.isfinite(d).all(1)
        if bad.any():
            out.append(Finding("HIGH", "data-file", f"{name}: {int(bad.sum())} non-finite rows"))
            d = d[~bad]
        E, A, V, U = d.T
        segs_f = [s for s in byfile[f] if s.active]
        is_ay = all(s.observable == "analyzing-power" for s in segs_f) and segs_f
        used = np.zeros(len(d), bool)
        for s in segs_f:
            mwin = in_window(d, s)
            used |= mwin
            thr = thresholds.get((s.entrance_key, s.exit_key))
            if thr and mwin.any() and (E[mwin] < thr - 1e-6).any():
                out.append(Finding("MEDIUM", "data-below-threshold", f"seg {s.key} {name}: {int((E[mwin] < thr).sum())} rows below the "
                                   f"{s.entrance_key}->{s.exit_key} lab threshold {thr:.4f} MeV"))
        if not is_ay and (V <= 0).any():
            sig = (V <= 0) & (np.abs(V) > 3 * U)        # significantly negative: a sign or units problem, not a near-zero point
            if sig.any():
                out.append(Finding("MEDIUM", "data-values", f"{name}: {int(sig.sum())} rows significantly negative (|value| > 3 error)", _rows(d, sig)))
            out.append(Finding("LOW", "data-values", f"{name}: {int((V <= 0).sum())} rows with value <= 0 (near-zero / background-subtracted?)", _rows(d, V <= 0)))
        if (U <= 0).any():
            out.append(Finding("MEDIUM", "data-values", f"{name}: {int((U <= 0).sum())} rows with error <= 0", _rows(d, U <= 0)))
        if not is_ay and (U > np.abs(V)).any():
            out.append(Finding("LOW", "data-values", f"{name}: {int((U > np.abs(V)).sum())} rows with error > |value|", _rows(d, U > np.abs(V))))
        key = np.round(d[:, :2], 6)
        _, cnt = np.unique(key, axis=0, return_counts=True)
        if (cnt > 1).any():
            out.append(Finding("LOW", "data-duplicates", f"{name}: {int((cnt > 1).sum())} duplicate (E, angle) keys ({int((cnt[cnt > 1] - 1).sum())} extra rows)"))
        if used.sum() < len(d):
            out.append(Finding("INFO", "data-rows-unused", f"{name}: {int(len(d) - used.sum())} of {len(d)} rows outside every active window"))
        # error pattern: a FLAG column read as percent errors; implausibly small errors
        # (a 0.02 % error typo once gave a few rows the weight of thousands of points)
        if not is_ay:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.abs(U / np.where(V == 0, np.nan, V)) * 100
            r = r[np.isfinite(r)]
            tiny = int((r < 0.1).sum()) if len(r) else 0
            if tiny:
                med = float(np.median(r))
                eq_pts = float((med / max(r.min(), 1e-12)) ** 2)      # weight of the worst row in typical points
                sev = "HIGH" if (eq_pts >= 100 or tiny >= 50 or tiny > 0.1 * len(r)) else "MEDIUM"
                out.append(Finding(sev, "data-error-tiny", f"{name}: {tiny} of {len(r)} rows have an error below 0.1 % of the value "
                                   f"(min {r.min():.3g} %, median {med:.3g} %: the worst row weighs like {eq_pts:,.0f} typical points; "
                                   "seen in an archived data set)"))
            if len(r) >= 8:
                distinct = np.unique(np.round(r, 1))
                integers = np.all(np.abs(distinct - np.round(distinct)) < 0.06) and distinct.max() <= 10 and distinct.min() >= 0.5
                if 3 <= len(distinct) <= 5 and integers:
                    out.append(Finding("HIGH", "data-error-pattern", f"{name}: error/value is {sorted(np.round(distinct).astype(int).tolist())} % "
                                       "(integer percent) in every row -- an EXFOR FLAG or quality code read as the uncertainty?"))
                elif len(distinct) == 2 and integers:
                    out.append(Finding("LOW", "data-error-pattern", f"{name}: errors are {sorted(np.round(distinct).astype(int).tolist())} % of the value in every row"))
                elif len(distinct) == 1 and len(r) >= 20:
                    out.append(Finding("INFO", "data-error-pattern", f"{name}: a uniform {distinct[0]:.3g} % error on every row"))
        # spikes
        if not is_ay:
            hits = _spikes(d, used)
            if hits:
                out.append(Finding("LOW", "data-spikes", f"{name}: {len(hits)} isolated spikes/dips (> 1.8x or < 0.55x the neighbours, > 5 sigma)",
                                   "; ".join(f"E {h[0]:.5g} ang {h[1]:g} ratio {h[2]:.2f} ({h[3]:.0f} sigma)" for h in hits[:4])))
        keysets[name] = {(round(r_[0], 4), float(f"{r_[2]:.5g}"), float(f"{r_[3]:.4g}")) for r_ in d}
    names = sorted(keysets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            n = len(keysets[a] & keysets[b])
            if n >= 3:
                out.append(Finding("MEDIUM", "data-crossdup", f"{a} & {b}: {n} identical (E, value, error) rows -- the same data twice?"))

    # -- target integration on test segments ----------------------------------------------------------------
    try:
        from .model import target_int_entries
        tkeys = {t.key for t in tests if t.active}
        for active, keys, line in target_int_entries(text):
            if active and tkeys & set(keys):
                out.append(Finding("HIGH", "targetint-test", f"active <targetInt> line names test segment(s) {sorted(tkeys & set(keys))}: "
                                   "a mode-3 extrapolation here is FOLDED (Golden rule 7); run it on a copy with the line inactive"))
    except Exception as err:
        out.append(Finding("INFO", "targetint-test", f"not checked ({err})"))

    # -- legacy parameter files, stale caches ------------------------------------------------------------------
    for pf in sorted(glob.glob(os.path.join(base, "*.sav")) + glob.glob(os.path.join(base, "output", "param.*"))):
        try:
            with open(pf) as fh:
                head = [next(fh, "") for _ in range(5)]
        except OSError:
            continue
        if any(LEGACY_RE.match(l.split()[0]) for l in head if l.split()):
            out.append(Finding("HIGH", "legacy-sav", f"{os.path.relpath(pf, base)} uses pre-2025-08 parameter names (j=.._la=..): the current "
                               "AZURE2 ignores them silently -- convert (azure2-legacy-params skill)"))
    od = os.path.join(base, "output")
    for cache in glob.glob(os.path.join(od, "intEC.*")):
        if os.path.getmtime(cache) < os.path.getmtime(azr_path):
            out.append(Finding("LOW", "intec-stale", f"{os.path.relpath(cache, base)} is older than the .azr: delete it if segments changed (Golden rule 8)"))

    # -- target energies (beam vs effective) -------------------------------------------------------------------
    chk = os.path.expanduser("~/.claude/skills/azure2-target-energies/check_target_energies.py")
    if run_target_energies and os.path.exists(chk) and "<targetInt>" in text and _has_active_target_int(text):
        import subprocess, sys
        try:
            r = subprocess.run([sys.executable, chk, azr_path], capture_output=True, text=True, timeout=300, cwd=base)
            tail = (r.stdout + r.stderr).strip().splitlines()[-3:]
            sev = "INFO" if r.returncode == 0 else "MEDIUM"
            out.append(Finding(sev, "target-energies", f"check_target_energies exit {r.returncode}", " | ".join(tail)))
        except Exception as err:
            out.append(Finding("INFO", "target-energies", f"not run ({err})"))
    return out


def _has_active_target_int(text):
    for line in _block(text, "targetInt").splitlines():
        t = line.split()
        if t and t[0] == "1":
            return True
    return False


def _rows(d, mask, n=3):
    idx = np.where(mask)[0][:n]
    return "; ".join("[" + ", ".join(f"{v:.5g}" for v in d[i]) + "]" for i in idx)


def _spikes(d, used):
    hits = []
    for a in np.unique(d[:, 1]):
        idx = np.where((d[:, 1] == a) & used)[0]
        idx = idx[np.argsort(d[idx, 0], kind="stable")]
        for k in range(1, len(idx) - 1):
            p, i, n = idx[k - 1], idx[k], idx[k + 1]
            vp, v, vn = d[p, 2], d[i, 2], d[n, 2]
            nb = 0.5 * (vp + vn)
            if nb <= 0 or abs(vp - vn) / nb > 0.25 or d[i, 3] <= 0:
                continue
            r = v / nb
            if (r > 1.8 or r < 0.55) and abs(v - nb) / d[i, 3] > 5:
                hits.append((d[i, 0], a, r, (v - nb) / d[i, 3]))
    return hits


def _lab_thresholds(model):
    """``{(entrance, exit): lab threshold of the entrance projectile (MeV)}`` for particle exits, 0 when open at 0."""
    out = {}
    try:
        pairs = pair_lines(model)
    except Exception:
        return out
    for ek, ec in pairs.items():
        if ec.ptype != 0:
            continue
        try:
            k = (ec.M1 + ec.M2) / ec.M2          # E_lab = k * E_cm
        except Exception:
            continue
        for xk, xc in pairs.items():
            if xc.ptype != 0:
                continue
            # E_cm threshold of exit pair relative to the entrance pair: (S_entrance + Ex_entrance) is Ex of the compound
            # nucleus at zero entrance energy; exit opens at Ex = S_exit + e_exit
            try:
                ex_open = xc.sep_energy + xc.excitation
                ex_zero = ec.sep_energy + ec.excitation
            except Exception:
                continue
            out[(ek, xk)] = max(0.0, (ex_open - ex_zero) * k)
            out[(ek, -1)] = 0.0
    return out


def write_audit(findings: List[Finding], path: str, azr_path: str) -> dict:
    header = [f"Model `{azr_path}`.  Codes and meanings: `rmfit/audit.py` docstring.  Findings are candidates until checked "
              "against the files; data-treatment decisions are the evaluator's.", ""]
    return write_findings_md(path, "Pre-flight audit (rmfit audit)", findings, header)


def main(argv=None):
    import argparse, json
    ap = argparse.ArgumentParser(prog="rmfit audit")
    ap.add_argument("azr")
    ap.add_argument("--campaign", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--no-target-energies", action="store_true")
    ap.add_argument("--kind", default=None, choices=["input", "data", "literature", "all"],
                    help="write audits/<kind>/PROMPT.md for a background subagent (next to the .azr, or under --campaign)")
    a = ap.parse_args(argv)
    cfg = None
    if a.campaign and os.path.exists(os.path.join(a.campaign, "campaign.json")):
        cfg = json.load(open(os.path.join(a.campaign, "campaign.json")))
    if a.kind:
        from . import audit_subagent as asub
        machine = audit_model(a.azr, campaign_dir=a.campaign, campaign_cfg=cfg, run_target_energies=False)
        spec = None
        sp = os.path.join(a.campaign or os.path.dirname(os.path.abspath(a.azr)), "spec.json")
        if os.path.exists(sp):
            from .spec import FitSpec
            spec = FitSpec.load(sp)
        for kind in (asub.KINDS if a.kind == "all" else (a.kind,)):
            r = asub.write_subagent_audit(kind, a.azr, a.campaign, spec, machine)
            print(f"{kind}: prompt {r['prompt']}  report {r['report']}")
        return 0
    fs = audit_model(a.azr, campaign_dir=a.campaign, campaign_cfg=cfg, run_target_energies=not a.no_target_energies)
    if a.json:
        print(json.dumps([f.__dict__ for f in fs], indent=1))
    else:
        for f in sorted(fs, key=lambda f: {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "INFO": 3}.get(f.severity, 9)):
            print(f.line())
    if a.out:
        counts = write_audit(fs, a.out, a.azr)
        print("wrote", a.out, counts)
    return 1 if any(f.severity == "HIGH" for f in fs) else 0


if __name__ == "__main__":
    raise SystemExit(main())
