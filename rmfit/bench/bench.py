#!/usr/bin/env python3
"""Benchmark suite (DESIGN.md Phase 3.1).

A case (``cases/<name>.json``) names a REFERENCE campaign (a finished, human-judged fit re-measured on the pinned engine; with
``reference.polish_rounds`` > 0 the job first polishes it to convergence, for a reference that is only a start file) and
RUNS (campaigns the drive produced from a degraded start of that reference).  ``score`` reads the ledgers only (no engine):

    objective          incumbent objective of the run and of the reference
    data chi2          sum of the per-data-set chi2 of the incumbent (ledger segment scores); penalty = objective - data chi2
    sign distance      keys.sign_distance between the two incumbents (same structure only)
    wait decisions     DECISIONS.md entries with default 'wait' (the ones that stop a run); human turns = answered by a person
    node hours         drive_state.json wall_hours (one node), with the slot count
    new review flags   HIGH/MEDIUM findings of the run's REVIEW.md that the reference's REVIEW.md does not have

and applies the recovery yardstick adopted 2026-11-09:  PASS = data chi2 <= reference's, objective within 0.5 %, zero wait
decisions, no new HIGH/MEDIUM review flag; EXACT (secondary) = objective within 1e-4.  Rows go to ``results.tsv``.

``run`` writes the SGE job for one run from the case's ``reference`` and the run's ``start`` and ``spec`` entries (reference
init if missing, ``rmfit perturb``, spec copied, ``rmfit drive`` to DONE / WAITING / STOP), modelled on the development
benchmarks' job scripts; ``--submit`` queues it (needs a Kerberos ticket: see the crc-cluster-jobs skill)."""
import argparse, json, os, re, subprocess, sys, time
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results.tsv")
COLUMNS = ["date", "case", "tier", "run", "label", "status", "reference_objective", "objective", "objective_rel_pct", "data_chi2",
           "data_chi2_ref", "penalty", "penalty_ref", "sign_distance", "wait_decisions", "human_turns", "node_hours", "slots",
           "new_flags", "pass_data", "pass_objective", "pass_decisions", "pass_flags", "verdict", "exact_1e-4", "note"]
YARD = dict(objective_rel_max=0.005, exact_rel=1e-4)

_FLAG_RE = re.compile(r"^- \*\*(HIGH|MEDIUM|LOW|INFO)\*\* `([^`]+)`: ([^:]*)")


def load_case(path: str) -> dict:
    case = json.load(open(path))
    case["_path"] = os.path.abspath(path)
    base = case.get("base_dir") or os.path.dirname(case["_path"])
    if not os.path.isabs(base):
        base = os.path.normpath(os.path.join(os.path.dirname(case["_path"]), base))
    case["_base"] = base
    return case


def review_flags(campaign_dir: str) -> Optional[set]:
    """(severity, code, subject) of every HIGH/MEDIUM finding in <dir>/REVIEW.md; None when there is no review."""
    path = os.path.join(campaign_dir, "REVIEW.md")
    if not os.path.exists(path):
        return None
    out = set()
    for line in open(path):
        m = _FLAG_RE.match(line.strip())
        if m and m.group(1) in ("HIGH", "MEDIUM"):
            out.add((m.group(1), m.group(2), m.group(3).strip()))
    return out


def new_review_flags(run: set, ref: set) -> list:
    """The run's flags the reference does not carry: same severity and code, and a subject that is not the same
    (rmfit.review.same_flag_subject: a level that moved within its own width is the same level; 2026-11-21)."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    from rmfit.review import same_flag_subject
    return sorted(f for f in run if not any(f[0] == g[0] and f[1] == g[1] and same_flag_subject(f[2], g[2]) for g in ref))


def _campaign(d):
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    from rmfit.campaign import Campaign
    return Campaign(d, enforce_engine=False)      # ledgers only; a review (--review) runs on whatever engine is imported (noted)


def incumbent_split(c) -> dict:
    cid, x, obj, keys = c.incumbent()
    tab = c.dataset_table(cid)
    data = sum(float(d["chi2"]) for d in tab)
    return dict(cid=cid, x=x, objective=float(obj), keys=keys, data_chi2=data, penalty=float(obj) - data, table=tab)


def decisions_summary(run_dir: str, pre_answered=(), inherited_review=False) -> dict:
    """Wait decisions and human answers of a run.  ``pre_answered`` (case file, per run): decision slugs that are part of the
    case definition rather than of the fit -- e.g. `audit-high` on a reference whose archived inputs carry known findings
    (e.g. normError-as-fraction segments shared by reference and start); the supervisor answers them once and
    the yardstick does not count them.  They are listed in the row's note.  ``inherited_review`` (2026-11-10):
    every HIGH flag of the run's REVIEW.md is in the reference's own REVIEW.md, so a `review-high` decision asks about the
    reference's structure, not the run's (the evaluator's "no new review flag" criterion); it is excluded like a pre-answered one."""
    from rmfit.decisions import DecisionQueue as Decisions
    path = os.path.join(run_dir, "DECISIONS.md")
    if not os.path.exists(path):
        return dict(wait=0, human=0, n=0, pre_answered=[])
    ds = Decisions(path).all()
    pre = [d for d in ds if any(d["slug"].startswith(x) for x in pre_answered)
           or (inherited_review and d["slug"] == "review-high")]
    ds = [d for d in ds if d not in pre]
    wait = sum(1 for d in ds if (d.get("default_if_unanswered") or "wait").strip() == "wait")
    human = sum(1 for d in ds if d.get("answer") and (d.get("answered_by") or "").strip() not in ("", "default", "rmfit drive", "test"))
    return dict(wait=wait, human=human, n=len(ds), pre_answered=[f"{d.id} {d['slug']}" for d in pre])


def score_run(case: dict, run_name: str, review_missing: bool = False) -> dict:
    ref_dir = os.path.join(case["_base"], case["reference"]["campaign"])
    run = case["runs"][run_name]
    run_dir = os.path.join(case["_base"], run["dir"])
    r, p = _campaign(ref_dir), _campaign(run_dir)
    R, P = incumbent_split(r), incumbent_split(p)
    row: Dict[str, object] = dict(date=time.strftime("%Y-%m-%d"), case=case["name"], tier=case.get("tier", ""), run=run_name,
                                  label=run.get("label", ""), reference_objective=round(R["objective"], 3), objective=round(P["objective"], 3),
                                  data_chi2=round(P["data_chi2"], 1), data_chi2_ref=round(R["data_chi2"], 1),
                                  penalty=round(P["penalty"], 1), penalty_ref=round(R["penalty"], 1))
    rel = P["objective"] / R["objective"] - 1.0
    row["objective_rel_pct"] = round(100 * rel, 3)
    # sign distance (same structure only)
    try:
        from rmfit.keys import take, sign_distance
        if [str(k) for k in R["keys"]] == [str(k) for k in P["keys"]]:
            dist, _ = sign_distance(P["x"], R["x"], P["keys"], p.ledger.fixed_of(p.structure_id))
            row["sign_distance"] = int(dist)
        else:
            xr = take(R["x"], R["keys"], P["keys"], P["x"])
            dist, _ = sign_distance(P["x"], xr, P["keys"], p.ledger.fixed_of(p.structure_id))
            row["sign_distance"] = f"{int(dist)} (key sets differ)"
    except Exception as err:  # noqa: BLE001
        row["sign_distance"] = f"n/a ({err})"
    st_path = os.path.join(run_dir, "drive_state.json")
    st = json.load(open(st_path)) if os.path.exists(st_path) else {}
    row["status"] = f"{st.get('status', '?')}@{st.get('state', '?')} cycle {st.get('cycle', '?')}" if st else "no drive_state.json"
    row["node_hours"] = round(float(st.get("wall_hours", 0.0)), 2) if st else ""
    row["slots"] = p.cfg.get("nshards", "")
    # review flags
    stale = []
    for c, I in ((r, R), (p, P)):
        rv = os.path.join(c.dir, "REVIEW.md")
        have = _review_cid(rv)
        if (not os.path.exists(rv) or (have is not None and have != int(I["cid"]))):
            if review_missing:
                from rmfit.review import review_campaign
                from rmfit import enginepin
                cur, pin = enginepin.fingerprint(), c.cfg.get("engine") or {}
                review_campaign(c, grid_check=False, figures=0)
                if pin.get("so_md5") and pin["so_md5"] != cur.get("so_md5"):
                    stale.append(f"{os.path.basename(c.dir)} reviewed on engine {enginepin.short(cur)} (pinned {enginepin.short(pin)})")
            elif os.path.exists(rv):
                stale.append(f"{os.path.basename(c.dir)} REVIEW.md is of candidate {have}, incumbent {I['cid']}")
    fr, fp = review_flags(ref_dir), review_flags(run_dir)
    fresh = new_review_flags(fp, fr) if fr is not None and fp is not None else []
    inherited = fr is not None and fp is not None and not [f for f in fresh if f[0] == "HIGH"] and bool({f for f in fp if f[0] == "HIGH"})
    dec = decisions_summary(run_dir, tuple(run.get("pre_answered") or ()), inherited_review=inherited)
    row["wait_decisions"], row["human_turns"] = dec["wait"], dec["human"]
    pre_note = ((" pre-answered by the case definition" + (" / review-high on a HIGH flag the reference's REVIEW.md carries" if inherited else "")
                 + ": " + ", ".join(dec["pre_answered"])) if dec.get("pre_answered") else "")
    if stale:
        row["note"] = "; ".join(stale)
    if fr is None or fp is None:
        row["new_flags"] = "n/a (no REVIEW.md in " + ", ".join(n for n, f in (("reference", fr), ("run", fp)) if f is None) + ")"
        pass_flags = None
    else:
        new = fresh
        row["new_flags"] = "; ".join(f"{s} {c} {n}" for s, c, n in new) if new else "none"
        pass_flags = not new
    row["pass_data"] = P["data_chi2"] <= R["data_chi2"] * (1 + 1e-9)
    row["pass_objective"] = rel <= YARD["objective_rel_max"]
    row["pass_decisions"] = dec["wait"] == 0
    row["pass_flags"] = "n/a" if pass_flags is None else pass_flags
    ok = row["pass_data"] and row["pass_objective"] and row["pass_decisions"] and (pass_flags is not False)
    row["verdict"] = "PASS" if ok else "FAIL"
    row["exact_1e-4"] = abs(rel) <= YARD["exact_rel"]
    if rel < -0.01:
        # a run that beats the reference by more than 1 % says the REFERENCE is not a converged fit on this engine / objective
        # (13N tier 0, 2026-11-09: the fixture's .azr holds start values, its fit lives in param.par): the row is not a recovery score
        # -- unless the run carries levels the reference lacks (a benchmark, 2026-11-25: every deleted level recovered
        # plus two extra levels worth about the whole difference); then the gain is structure, not convergence
        extra, missing = level_difference(os.path.join(run_dir, "REVIEW.md"), os.path.join(ref_dir, "REVIEW.md"))
        if extra is None or not (extra or missing):
            row["verdict"] = "REFERENCE-SUSPECT"
            stale.append(f"run beats the reference by {-100 * rel:.1f} %: the reference is not a converged fit; rebuild it")
        else:
            row["verdict"] = "EXTRA-LEVELS"
            stale.append(f"run beats the reference by {-100 * rel:.1f} % with a different level set (levels inside the data range): "
                         f"extra {', '.join(extra) or 'none'}; missing {', '.join(missing) or 'none'} -- not a recovery score at equal "
                         f"structure; judge the extra levels (review flags) before rebuilding the reference")
    row["note"] = "; ".join(t for t in (run.get("note", ""), row.get("note", ""), "; ".join(stale[1:] if row.get("note") else stale), pre_note.strip()) if t)
    return row


def review_levels(path: str) -> Optional[List[tuple]]:
    """(J^pi, Ex MeV, Gamma text) rows of a REVIEW.md "Levels inside the data range" table, or None without one."""
    if not os.path.exists(path):
        return None
    rows, inside = [], False
    for line in open(path):
        if line.startswith("## "):
            inside = line.strip().startswith("## Levels inside the data range")
            continue
        m = re.match(r"\|\s*(\d+/?\d*[+-])\s*\|\s*([0-9.]+)\s*\|\s*([^|]*?)\s*\|", line) if inside else None
        if m and m.group(3) not in ("?", ""):     # "?" = no width: a level switched off (PRUNE zeroes, does not delete)
            rows.append((m.group(1), float(m.group(2)), m.group(3)))
    return rows if rows else None


def level_difference(run_review: str, ref_review: str, tol_mev: float = 0.03):
    """Levels of the run with no reference level of the same J^pi within tol_mev (extra) and the converse (missing), as
    "J^pi E (Gamma)" strings; (None, None) when either REVIEW.md has no level table."""
    run, ref = review_levels(run_review), review_levels(ref_review)
    if run is None or ref is None:
        return None, None

    def unmatched(a, b):
        left, out = list(b), []
        for j, e, g in a:
            hit = next((x for x in left if x[0] == j and abs(x[1] - e) <= tol_mev), None)
            if hit is None:
                out.append(f"{j} {e:.4f} ({g})")
            else:
                left.remove(hit)
        return out
    return unmatched(run, ref), unmatched(ref, run)


def _review_cid(path: str) -> Optional[int]:
    """The candidate a REVIEW.md was written for (its header line), or None."""
    if not os.path.exists(path):
        return None
    for line in open(path):
        m = re.search(r"incumbent candidate (\d+)", line)
        if m:
            return int(m.group(1))
    return None


def append_rows(rows: List[dict]):
    new = not os.path.exists(RESULTS)
    with open(RESULTS, "a") as fh:
        if new:
            fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ") for c in COLUMNS) + "\n")


def print_rows(rows: List[dict]):
    for row in rows:
        print(f"{row['case']} / {row['run']} ({row['label']}): {row['verdict']}  objective {row['objective']} vs {row['reference_objective']} "
              f"({row['objective_rel_pct']:+} %), data chi2 {row['data_chi2']} vs {row['data_chi2_ref']}, penalty {row['penalty']} vs {row['penalty_ref']}, "
              f"sign distance {row['sign_distance']}, wait decisions {row['wait_decisions']}, human turns {row['human_turns']}, "
              f"{row['node_hours']} node-hours on {row['slots']} slots, new flags: {row['new_flags']}; exact {row['exact_1e-4']}; {row['status']}")


JOB = """#!/bin/csh
{mail_line}
#$ -m a
{queue_line}
{pe_line}
#$ -N {jobname}
#$ -cwd
#$ -j y
#$ -o {log}
# rmfit bench (DESIGN.md Phase 3.1): case {case} run {run} -- {label}
# reference campaign {refdir} from {refazr}; start = `rmfit perturb` ({flips} flips, seed {seed}, jitter {jitter}, ekick {ekick},
# polished {pnfev}); spec {spec}; then `rmfit drive` to DONE / WAITING / STOP.  Score afterwards with
#   python3 -m rmfit.bench score {casefile} {run}
{setup_lines}
setenv AZURE2_ROOT {azure2_root}
{engine_guard}setenv PYTHONPATH {rmatrix}:${{AZURE2_ROOT}}
setenv OMP_NUM_THREADS 1
setenv PYTHONUNBUFFERED 1
cd {base}
echo "host `hostname` `date` job $JOB_ID"
python3 -c "import rmfit, rmfit.drive, pyazr; print('imports', pyazr.__file__, rmfit.__file__)"
if ($status != 0) exit 1
set RM = "python3 -m rmfit.cli"
{chisq}if (! -e {refdir}/campaign.json) then
  echo "`date` reference init"
  $RM {refdir} init --base {refazr} {savarg}{reffixed} --nshards {nshards} {flagsarg} {policyarg} {searcharg} --note "bench {case} reference"
  if ($status != 0) then
    echo "`date` reference init FAILED"; exit 1
  endif
endif
{polish}if (! -e {rundir}/campaign.json) then
{startblock}endif
if (! -e {rundir}/spec.json) cp {spec} {rundir}/spec.json
{searchpatch}
echo "`date` drive"
$RM {rundir} drive --max-minutes 0
echo "`date` drive exit $status"
$RM {rundir} drive --dry-run
echo "`date` all done"
"""


# the same job for a Slurm site (site.scheduler() == "slurm"; rmfit/site.py): bash, header and setup from the site file;
# also the plain-shell job ("shell"), whose header is a lock, a log redirection and the job marker instead of #SBATCH lines
JOB_SLURM = """#!/bin/bash
{slurm_header}
# rmfit bench: case {case} run {run} -- {label}
# reference campaign {refdir} from {refazr}; start = `rmfit perturb` ({flips} flips, seed {seed}, jitter {jitter}, ekick {ekick},
# polished {pnfev}); spec {spec}; then `rmfit drive` to DONE / WAITING / STOP.  Score afterwards with
#   python3 -m rmfit.bench score {casefile} {run}
{setup_lines}
export AZURE2_ROOT={azure2_root}
{engine_guard}export PYTHONPATH={rmatrix}:${{AZURE2_ROOT}}
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
cd {base} || exit 1
echo "host $(hostname) $(date) job ${{SLURM_JOB_ID:-none}}"
python3 -c "import rmfit, rmfit.drive, pyazr; print('imports', pyazr.__file__, rmfit.__file__)" || exit 1
RM="python3 -m rmfit.cli"
{chisq}if [ ! -e {refdir}/campaign.json ]; then
  echo "$(date) reference init"
  $RM {refdir} init --base {refazr} {savarg}{reffixed} --nshards {nshards} {flagsarg} {policyarg} {searcharg} --note "bench {case} reference"
  if [ $? -ne 0 ]; then
    echo "$(date) reference init FAILED"; exit 1
  fi
fi
{polish}if [ ! -e {rundir}/campaign.json ]; then
{startblock}fi
[ -e {rundir}/spec.json ] || cp {spec} {rundir}/spec.json
{searchpatch}
echo "$(date) drive"
$RM {rundir} drive --max-minutes 0
echo "$(date) drive exit $?"
$RM {rundir} drive --dry-run
echo "$(date) all done"
"""


def write_job(case: dict, run_name: str) -> str:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    from rmfit import site
    ref, run = case["reference"], case["runs"][run_name]
    base = case["_base"]
    start = dict(dict(flips=8, seed=0, jitter=0.0, ekick=0.0, polish_nfev=300), **(run.get("start") or {}))

    def jarg(flag, key):
        v = ref.get(key)
        return f"{flag} @{v}" if v else ""
    spec = run.get("spec") or ""
    if not spec:
        raise SystemExit(f"run {run_name} has no 'spec' (path of the drive spec, relative to the case's base_dir)")
    def chisq_copy(side, azr):
        # the archived chiSquared.out of the fit the .azr + .sav came from, copied into the .azr's output directory so
        # that init's reproduce check runs (2026-11-10: a reference initialised at more than twice
        # its archived objective -- the .sav's fixed widths had been ignored -- and nothing warned because the file was not there)
        src = side.get("chisq_out")
        if not src:
            return ""
        from ..review import output_dir_of
        dst = os.path.join(output_dir_of(os.path.join(base, azr)), "chiSquared.out")
        return f"mkdir -p {os.path.dirname(dst)}\ncp -p {os.path.join(base, src)} {dst}\n"
    chisq = chisq_copy(ref, ref["azr"]) + (chisq_copy(start, start["azr"]) if start.get("azr") else "")
    sched = site.scheduler()
    slurm = sched != "sge"          # bash syntax (Slurm and plain shell); the name predates the shell kind
    # shell syntax of the generated snippets: csh for the SGE job, bash for the Slurm and plain-shell ones
    DATE = "$(date)" if slurm else "`date`"
    IF_MISSING = (lambda f: f"if [ ! -e {f} ]; then") if slurm else (lambda f: f"if (! -e {f}) then")
    END_IF = "fi" if slurm else "endif"
    LOOP = (lambda n: f"for i in $(seq 1 {n}); do") if slurm else (lambda n: f"foreach i (`seq 1 {n}`)")
    END_LOOP = "done" if slurm else "end"
    FAIL_CHECK = ((lambda what: f"  if [ $? -ne 0 ]; then\n    echo \"{DATE} {what} FAILED\"; exit 1\n  fi\n") if slurm else
                  (lambda what: f"  if ($status != 0) then\n    echo \"{DATE} {what} FAILED\"; exit 1\n  endif\n"))
    polish = ""
    if int(ref.get("polish_rounds") or 0) > 0:
        # a reference that is not a converged fit (a start file, a regression fixture) is polished to convergence first;
        # the marker file keeps a resubmitted job from polishing it again
        polish = (f"{IF_MISSING(ref['campaign'] + '/.bench_polished')}\n"
                  f"  echo \"{DATE} reference polish ({ref['polish_rounds']} rounds of {ref.get('polish_nfev', 300)})\"\n"
                  f"  {LOOP(int(ref['polish_rounds']))}\n"
                  f"    $RM {ref['campaign']} polish --nfev {int(ref.get('polish_nfev') or 300)} --tighten\n"
                  f"  {END_LOOP}\n"
                  f"  touch {ref['campaign']}/.bench_polished\n"
                  f"{END_IF}\n")
    if start.get("azr"):
        # a degraded start that is a FILE (the archived state before the human's key decision), not a random perturbation
        # of the reference: the run campaign is initialised on it with the reference's engine settings
        startblock = (f"  echo \"{DATE} start from {start['azr']}\"\n"
                      f"  $RM {run['dir']} init --base {start['azr']} {('--sav ' + start['sav']) if start.get('sav') else ''}"
                      f"{(' --sav-fixed ' + start['sav_fixed']) if start.get('sav_fixed') else ''} "
                      f"--nshards {int(ref.get('nshards') or 8)} {jarg('--flags', 'flags')} {jarg('--policy', 'policy')} {jarg('--search', 'search')} "
                      f"--note \"bench {case['name']} start {start['azr']}\"\n"
                      f"{FAIL_CHECK('start init')}")
    else:
        startblock = (f"  echo \"{DATE} perturb\"\n"
                      f"  $RM {ref['campaign']} perturb --out {run['dir']} --flips {start['flips']} --seed {start['seed']} --jitter {start['jitter']} "
                      f"--ekick {start['ekick']} --polish-nfev {start['polish_nfev']}\n"
                      f"{FAIL_CHECK('perturb')}")
    searchpatch = ""
    if run.get("search"):
        searchpatch = (f"python3 -c \"import json; c=json.load(open('{run['dir']}/campaign.json')); c['search']=json.load(open('{run['search']}')); "
                       f"json.dump(c, open('{run['dir']}/campaign.json','w'), indent=2); print('search ->', c['search'])\"")
    rmatrix = os.path.dirname(os.path.dirname(HERE))
    slots = int(run.get("slots") or ref.get("nshards") or 8)
    jobname = f"bench_{case['name']}_{run_name}"[:40]
    if sched == "shell":
        header = (f"# plain shell (no scheduler; rmfit/site.py \"scheduler\": \"shell\"): start it where {slots} cores are free with\n"
                  f"#   nohup bash {os.path.join(base, 'run_crc_bench_' + run_name)} >/dev/null 2>&1 &   (or: bench run --submit)\n"
                  + site.shell_lock_lines(os.path.join(base, f".rmfit_shell_bench_{run_name}.lock")) + "\n"
                  f"exec >> {os.path.join(base, 'bench_' + run_name + '.log')} 2>&1\n"
                  "export RMFIT_SHELL_JOB=$$                # this script is the job here (the drive refuses engine states outside one)")
        extra = dict(slurm_header=header, setup_lines=site.setup_lines())
    elif slurm:
        extra = dict(slurm_header=site.slurm_header(jobname, slots, f"bench_{run_name}.log", "a"), setup_lines=site.setup_lines())
    else:
        extra = dict(mail_line=site.mail_line(), setup_lines=site.setup_lines(), queue_line=site.queue_line(),
                     pe_line=site.pe_line(slots))
    # a frozen engine of AZURE2_engines/ (freeze_engine.sh, MANIFEST): verify it before any calculation (2026-11-26,
    # skill azure2-frozen-engines "build guard"); other engines (pyazr 2 trees) have no manifest and no guard
    _root = ref.get("azure2_root") or os.path.expanduser("~/AZURE2")
    _fz = os.path.join(os.path.dirname(os.path.abspath(_root)), "freeze_engine.sh")
    engine_guard = ""
    if os.path.isfile(os.path.join(_root, "MANIFEST")) and os.path.isfile(_fz) and os.path.basename(_root).startswith("AZURE2_frozen_"):
        engine_guard = (f"bash {_fz} --check {os.path.basename(_root)[len('AZURE2_frozen_'):]} || exit 1   "
                        "# frozen engine intact (MANIFEST md5s)\n")
        if os.path.isfile(os.path.join(_root, "coulomb_table.bin")):    # the engine's own table, as its use_engine does
            ct = os.path.join(_root, "coulomb_table.bin")
            engine_guard += (f"export AZURE2_COULOMB_TABLE={ct}\n" if (slurm or site.scheduler() == "shell")
                             else f"setenv AZURE2_COULOMB_TABLE {ct}\n")
    text = (JOB_SLURM if slurm else JOB).format(**extra, slots=slots, jobname=jobname, engine_guard=engine_guard,
                      log=f"bench_{run_name}.log", case=case["name"], run=run_name, label=run.get("label", ""), refdir=ref["campaign"],
                      refazr=ref["azr"], flips=start["flips"], seed=start["seed"], jitter=start["jitter"], ekick=start["ekick"],
                      pnfev=start["polish_nfev"], spec=spec, casefile=case["_path"], azure2_root=ref.get("azure2_root") or os.path.expanduser("~/AZURE2"),
                      rmatrix=rmatrix, base=base, savarg=(f"--sav {ref['sav']}" if ref.get("sav") else ""), nshards=int(ref.get("nshards") or 8),
                      reffixed=(f" --sav-fixed {ref['sav_fixed']}" if ref.get("sav_fixed") else ""),
                      flagsarg=jarg("--flags", "flags"), policyarg=jarg("--policy", "policy"), searcharg=jarg("--search", "search"),
                      rundir=run["dir"], searchpatch=searchpatch, polish=polish, startblock=startblock, chisq=chisq)
    path = os.path.join(base, f"run_crc_bench_{run_name}")
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python3 -m rmfit.bench", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("score", help="score finished runs of a case against its reference; append rows to results.tsv")
    p.add_argument("case")
    p.add_argument("runs", nargs="*", help="run names (default: every run of the case that has a ledger)")
    p.add_argument("--review", action="store_true", help="write a REVIEW.md (engine, ~1 min per campaign) where one is missing")
    p.add_argument("--no-append", action="store_true", help="print only; do not touch results.tsv")
    p = sub.add_parser("run", help="write the SGE job for one run (reference init, perturb, spec, drive)")
    p.add_argument("case")
    p.add_argument("run")
    p.add_argument("--submit", action="store_true")
    p = sub.add_parser("list", help="cases and their runs")
    p.add_argument("case", nargs="?")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        files = [a.case] if a.case else sorted(os.path.join(HERE, "cases", f) for f in os.listdir(os.path.join(HERE, "cases")) if f.endswith(".json"))
        for f in files:
            c = load_case(f)
            print(f"{c['name']} (tier {c.get('tier')}, {c.get('reaction')}): reference {c['reference']['campaign']}; runs: "
                  + ", ".join(f"{k} [{v.get('label', '')}]" for k, v in c["runs"].items()))
        return
    case = load_case(a.case)
    if a.cmd == "run":
        path = write_job(case, a.run)
        print("job script:", path)
        if a.submit:
            from rmfit import site
            if site.scheduler() == "shell":       # no scheduler: start it detached here; it logs to bench_<run>.log itself
                pr = subprocess.Popen(["nohup", "bash", path], cwd=case["_base"], stdin=subprocess.DEVNULL,
                                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
                print(f"started in the background, pid {pr.pid}")
            else:
                out = subprocess.run([site.submit_cmd(), os.path.basename(path)], cwd=case["_base"], capture_output=True, text=True)
                print(out.stdout.strip(), out.stderr.strip())
        return
    runs = a.runs or [k for k, v in case["runs"].items() if os.path.exists(os.path.join(case["_base"], v["dir"], "ledger.sqlite"))]
    rows = [score_run(case, r, review_missing=a.review) for r in runs]
    print_rows(rows)
    if not a.no_append:
        append_rows(rows)
        print("appended to", RESULTS)


if __name__ == "__main__":
    main()
