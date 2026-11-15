"""``PYTHONPATH=<directory containing rmfit> python3 -m rmfit.cli <campaign_dir> <command> [options]``

Commands (all idempotent and resumable from the ledger):

  init --base <azr> [--sav <sav>] [--nshards N] [--flags JSON] [--policy JSON] [--search JSON] [--objective JSON] [--cov-reference NPZ]   register the base model, verify the baseline (a relvector covariance's model is fixed at the seed's, cov_reference.npz), write run_crc_rmfit
  status [--json]                                               STATUS.md / a compact JSON summary
  polish [--nfev N] [--shifts-free] [--refreeze-every N]        warm-start continuation of the incumbent
  round [--jpi 5/2+ 3/2- ...] [--rounds N] [--top-k K]          sign round(s); touch <dir>/STOP to end the loop after the current round
  crossover --donor CID [--donor-dir DIR]                       basin crossover: the donor's sign pattern per channel/group/all, all polished
  recall [--levels N] [--min-width KEV] [--nfev N] [--pinned-nfev N] [--local-release MEV] [--local-nfev N]
                                                                recall test (2026-11-09): the N broadest known levels removed one at a
                                                                time and searched for with no J^pi -- validates level-add on this model
  population [--restarts K] [--flips 6,8,10,12] [--nfev N] [--passes P] [--jitter J] [--seed S]
                                                                restarts + crossover (2026-11-08): K sign-perturbed restarts of the incumbent
                                                                polished to convergence, crossover between them; the basin-crossing stage
  levelscan [--ex E ...] [--top K] [--gamma keV] [--refit N]    residual-driven level search: template screen of every J^pi at the residual peaks
  relocate --level 9/2+#1 --ex 8.4654 [--window keV]           move an existing level and re-polish with its energy confined near --ex
  verify <cid>                                                  re-score a ledger candidate on a fresh single session
  export [--path file.azr]                                      bake the incumbent for the GUI
  log "text"                                                    append a note to the readme
  engine                                                        the campaign's pinned engine vs the one this process imports (2026-11-07)
  spec check|render|derive [--spec spec.json]                   fit specification vs the .azr / campaign.json; SPEC.md; derived init config
  audit [--azr file] [--no-target-energies]                     pre-flight checks of the model and data -> AUDIT.md (no engine session)
  audit --kind input|data|literature|all                        write audits/<kind>/PROMPT.md (+ FINDINGS.md skeleton) for a background subagent
  audit --findings | --status KIND N confirmed|rejected|decision  the subagents' FINDINGS.md: unverified rows / mark one verified
  review [--grid-check] [--figures N] [--top K] [--close-fraction f]  post-fit "GUI eye" on the incumbent -> REVIEW.md + figures/review_*.png (2026-11-07);
                                                                engine-free on any finished CLI run: python3 -m rmfit.review <file.azr> [--output DIR]
  drive [--base azr --sav sav] [--max-minutes M] [--dry-run] [--local] [--state S] [--write-job]
                                                                the closed loop (rmfit/drive.py, 2026-11-08): build phase (audit,
                                                                reproduce, spec, anchor norms, seeds, pairs, candidates, init) then
                                                                tighten / signs / adds / dataset-off / review / export per spec.json;
                                                                one leg per call, state in drive_state.json, DRIVE.md for the supervisor
  decisions                                                     pending entries of DECISIONS.md
  decide <id> "<choice>" [--by who] [--note ...]                answer one

Engine pin (2026-11-07): ``init`` records the imported engine (AZURE2_ROOT, _azure2 .so md5, commit) in
campaign.json; every evaluating command refuses to run on another engine unless RMFIT_ALLOW_ENGINE_CHANGE=1.
status/report/log/engine only warn.
"""
READ_ONLY_CMDS = {"status", "report", "log", "engine", "spec", "audit", "decisions", "decide"}

import argparse
import json
import os
import sys


def _json_arg(text):
    """A JSON option value, or ``@path`` to read it from a file (2026-09-29: csh brace-expands a
    JSON literal passed through a shell variable -- "Missing '}'." -- so long policies go in a file)."""
    if isinstance(text, str) and text.startswith("@"):
        with open(text[1:]) as fh:
            return json.load(fh)
    return json.loads(text)


def in_sge_job_env():
    """Inside a job: SGE (JOB_ID with a task, PE or slot variable), Slurm (SLURM_JOB_ID) or a plain-shell script of rmfit's
    (RMFIT_SHELL_JOB, site scheduler "shell").  Name historical."""
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("RMFIT_SHELL_JOB"):
        return True
    return bool(os.environ.get("JOB_ID")) and bool(os.environ.get("SGE_TASK_ID") or os.environ.get("PE_HOSTFILE") or os.environ.get("NSLOTS"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rmfit")
    ap.add_argument("campaign")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init")
    p.add_argument("--base", required=True)
    p.add_argument("--sav")
    p.add_argument("--nshards", type=int, default=24)
    p.add_argument("--sav-fixed", choices=("ignore", "apply"), default="ignore",
                   help="when the .sav's FIXED level parameters differ from the .azr's (2026-11-10): ignore = the .azr's are used and the "
                        ".sav seeds the free values only (default, a script-built base); apply = the pair is a CLI fit whose param.sav "
                        "holds the fixed values it ran with, so <stem>_sav.azr with the file applied becomes the base.  Either way the "
                        "disagreement and the objective both ways are logged")
    p.add_argument("--note", default="")
    p.add_argument("--policy", default="{}", help="JSON for BoundsPolicy fields")
    p.add_argument("--search", default="{}", help="JSON for SearchConfig fields")
    p.add_argument("--flags", default=None, help="JSON engine flags (default: the archive flags; CLI-only flags vary per reaction)")
    p.add_argument("--objective", default="{}", help='JSON ObjectiveSpec, e.g. {"dataset_weight": "reduced", "loss": "sivia"}; @file reads it from a file (also --policy/--search/--flags)')
    p.add_argument("--cov-reference", default=None,
                   help="relvector covariance: fix its model at this .npz (another campaign's cov_reference.npz) "
                        "instead of the seed's own (default; 2026-10-25)")
    p.add_argument("--spec", default=None, help="spec.json (rmfit/spec.py): flags/policy/search/objective/nshards derived from it "
                                                "unless given explicitly; copied into the campaign; HIGH disagreements with the .azr refuse init")
    p = sub.add_parser("status")
    p.add_argument("--json", action="store_true")
    p = sub.add_parser("polish")
    p.add_argument("--nfev", type=int, default=200)
    p.add_argument("--candidate", type=int, default=None, help="polish this ledger candidate instead of the best verified one (2026-10-05)")
    p.add_argument("--shifts-free", action="store_true")
    p.add_argument("--tighten", action="store_true", help="walk out-of-cap amplitudes into the theta^2 box in stages")
    p.add_argument("--refreeze-every", type=int, default=0,
                   help="relvector covariance WITHOUT a fixed reference only: re-freeze its model every <= N "
                        "evaluations (default 0 = never, since 2026-10-25: the iteration oscillated on a large model)")
    p = sub.add_parser("round")
    p.add_argument("--jpi", nargs="*")
    p.add_argument("--local-release", type=float, default=None,
                   help="MeV: screen + stage-1 polish free every level within +-this of the flipped one, all background "
                        "widths and the norms (2026-10-02); default: the campaign's search config (0 = J^pi group only)")
    p.add_argument("--rounds", type=int, default=1)
    p.add_argument("--top-k", type=int)
    p.add_argument("--note", default="")
    p = sub.add_parser("crossover")
    p.add_argument("--donor", type=int, required=True, help="donor candidate id")
    p.add_argument("--donor-dir", help="another campaign directory holding the donor (same base structure)")
    p.add_argument("--note", default="")
    p = sub.add_parser("levelscan")
    p.add_argument("--ex", type=float, nargs="*", help="candidate excitation energies (MeV); default: from residual runs")
    p.add_argument("--top", type=int, default=3)
    p.add_argument("--gamma", default="10", help="template total width (keV), or 'auto' = each cluster's own residual-run "
                                                 "width clamped into --gamma-bounds (2026-09-24)")
    p.add_argument("--gamma-bounds", default="10,1000", help="keV lo,hi clamp for --gamma auto")
    p.add_argument("--seg-match", default=None, help="regex on data-file names: only these segments' residuals make candidates (2026-09-24)")
    p.add_argument("--rank", default="chi2", choices=["chi2", "density"],
                   help="order residual clusters by chi2 (old) or chi2 per MeV of structure (2026-09-24)")
    p.add_argument("--ex-range", default=None, help="lo,hi MeV: drop residual candidates outside before the --top cut")
    p.add_argument("--scales", default="1", help="comma-separated multipliers of --gamma to try")
    p.add_argument("--jpi", nargs="*")
    p.add_argument("--gn-top", type=int, default=20)
    p.add_argument("--max-signs", type=int, default=32)
    p.add_argument("--refit", type=int, default=0, help="refit the best N candidates (run_level_add)")
    p.add_argument("--differential-window", type=float, default=250.0,
                   help="keV half-width: re-rank candidates on the DIFFERENTIAL segments near the candidate "
                        "energy (0 = rank on the total objective, the pre-redesign behaviour)")
    p.add_argument("--grid", default="5,20",
                   help="resonanceWidthMultiplier,pointsPerWidth for the SCREEN variant only; default 5,20 is "
                        "6.7x cheaper than the converged 20,50 for 0.03%%; 'none' keeps the model's setting")
    p.add_argument("--ex-offsets", default="-0.04,0,0.04", help="comma-separated energy offsets (MeV) tried around each candidate")
    p.add_argument("--stage-scales", default="1.2,5",
                   help="two-stage screen: every shape at these --gamma multipliers and offset 0 first; "
                        "'none' evaluates every shape at every scale and offset (tens of thousands of patterns, hours on a large model)")
    p.add_argument("--refine-top", type=int, default=3, help="shapes per J^pi carried into the all-scales x all-offsets stage")
    p.add_argument("--merge-kev", type=float, default=60.0,
                   help="residual candidates closer than this to a stronger one are not screened separately (0 = screen all)")
    p.add_argument("--note", default="")
    p = sub.add_parser("relocate")
    p.add_argument("--level", required=True, help="level label, e.g. 9/2+#1")
    p.add_argument("--ex", type=float, required=True, help="new excitation energy (MeV)")
    p.add_argument("--window", type=float, default=5.0, help="half window (keV) kept around --ex")
    p.add_argument("--nfev", type=int, default=100)
    p.add_argument("--note", default="")
    p = sub.add_parser("verify")
    p.add_argument("cid", type=int)
    p = sub.add_parser("export")
    p.add_argument("--path")
    p = sub.add_parser("report")
    p = sub.add_parser("structure")
    p.add_argument("what", choices=["add", "remove", "dataset-off"])
    p.add_argument("--candidate", type=int, help="index into candidates.json")
    p.add_argument("--file", help="data-file substring for dataset-off")
    p.add_argument("--jpi", nargs="*", help="J^pi groups to scan (default: the candidate's, or all)")
    p.add_argument("--nfev", type=int, default=150)
    p.add_argument("--local-release", type=float, default=0.0,
                   help="add: MeV half-width of a local-release polish (all J^pi near the candidate + background + norms) "
                        "between the pinned and the global stage; the control gets the same (2026-10-02)")
    p.add_argument("--local-nfev", type=int, default=80)
    p.add_argument("--pinned-nfev", type=int, default=60)
    p.add_argument("--top", type=int, default=3,
                   help="remove: how many of the cheapest-to-remove levels get an honest refit")
    p = sub.add_parser("engine", help="show the pinned engine and the imported one")
    p = sub.add_parser("spec", help="fit specification: check / render / derive")
    p.add_argument("what", choices=["check", "render", "derive"])
    p.add_argument("--spec", default=None, help="path (default <campaign>/spec.json)")
    p.add_argument("--azr", default=None, help="model to check (default the campaign's base model)")
    p = sub.add_parser("audit", help="pre-flight checks -> AUDIT.md; --kind writes a subagent prompt; --findings lists their reports")
    p.add_argument("--azr", default=None, help="model (default the campaign's base model)")
    p.add_argument("--no-target-energies", action="store_true")
    p.add_argument("--json", action="store_true")
    p.add_argument("--kind", default=None, choices=["input", "data", "literature", "all"],
                   help="write audits/<kind>/PROMPT.md (+ empty FINDINGS.md) for a background subagent instead of running the machine checks")
    p.add_argument("--findings", action="store_true", help="summarize audits/*/FINDINGS.md (unverified rows first)")
    p.add_argument("--status", nargs=3, metavar=("KIND", "N", "STATUS"), default=None,
                   help="mark finding N of audits/KIND/FINDINGS.md as confirmed | rejected | decision | candidate")
    p.add_argument("--note", default="", help="with --status: short note written next to the status")
    p = sub.add_parser("review", help="post-fit review of the incumbent (per-point summary, norm direction/trends, residual runs, "
                                      "average-vs-structure, zigzag, level support) -> REVIEW.md + figures/review_*.png")
    p.add_argument("--grid-check", action="store_true", help="re-score the incumbent on a converged target-integration grid (20/50)")
    p.add_argument("--figures", type=int, default=6, help="figures for the N worst data sets (0 = none)")
    p.add_argument("--top", type=int, default=20, help="rows in the top-points / runs tables")
    p.add_argument("--close-fraction", type=float, default=None,
                   help="flag same-J^pi pairs with a gap below this fraction of the narrower width (default 0.10; 0 = off)")
    p = sub.add_parser("uq", help="uncertainty quantification of the incumbent (rmfit.uq): export -> saturation -> Laplace -> NUTS -> gates -> "
                                  "linear vs sampled bands -> UQ.md; resumable (uq/uq_state.json); settings from spec.json 'uq' or --uq JSON")
    p.add_argument("--uq", default="{}", help="JSON overriding the spec's uq section (chains, warmup, draws, bands, ...)")
    p.add_argument("--nprocs", type=int, help="chain processes (default = chains)")
    p.add_argument("--threads", type=int, default=1, help="OMP threads per chain process")
    p.add_argument("--stop-after", choices=("EXPORT", "SETUP", "SAMPLE", "GATES", "BANDS", "REPORT"))
    p.add_argument("--reset", action="store_true", help="start the stage over (the old uq_state.json is kept with a timestamp)")
    p = sub.add_parser("drive", help="the closed loop: run states from drive_state.json until DONE / WAITING / STOP or --max-minutes")
    p.add_argument("--base", default=None, help="first call of the build phase: the base .azr (then <dir>/spec.json must exist)")
    p.add_argument("--sav", default=None, help="its parameter file (.sav)")
    p.add_argument("--max-minutes", type=float, default=0.0, help="wall cap of this leg (0 = none); the job script resubmits on CONTINUE")
    p.add_argument("--dry-run", action="store_true", help="print the plan from the current state; nothing runs")
    p.add_argument("--local", action="store_true", help="allow engine states outside an SGE job (fixture tests only)")
    p.add_argument("--state", default=None, help="jump to this state before running (e.g. after fixing an input: SCALE_ANCHOR)")
    p.add_argument("--reset", action="store_true", help="start the state machine over (the old drive_state.json is kept with a timestamp)")
    p.add_argument("--write-job", action="store_true", help="(re)write run_crc_drive and exit")
    p = sub.add_parser("decisions", help="pending entries of DECISIONS.md")
    p.add_argument("--all", action="store_true")
    p = sub.add_parser("decide", help="answer a DECISIONS.md entry")
    p.add_argument("id")
    p.add_argument("choice")
    p.add_argument("--by", default="")
    p.add_argument("--note", default="")
    p = sub.add_parser("log")
    p.add_argument("text")
    p = sub.add_parser("recall", help="recall test: the N broadest known levels removed one at a time and searched for with no J^pi given")
    p.add_argument("--levels", type=int, default=2)
    p.add_argument("--level", action="append", default=None, help="explicit level label(s) such as 3/2-#1 instead of the broadest")
    p.add_argument("--min-width", type=float, default=10.0, help="keV; only levels at least this wide are removed")
    p.add_argument("--nfev", type=int, default=150)
    p.add_argument("--pinned-nfev", type=int, default=60)
    p.add_argument("--local-release", type=float, default=0.0)
    p.add_argument("--local-nfev", type=int, default=80)
    p.add_argument("--note", default="")
    p = sub.add_parser("population", help="restarts + crossover: K sign-perturbed restarts of the incumbent polished to convergence, crossover between them")
    p.add_argument("--restarts", type=int, default=4)
    p.add_argument("--flips", default="6,8,10,12", help="flip counts drawn per restart (comma list) or one integer")
    p.add_argument("--per-group", type=int, default=3)
    p.add_argument("--nfev", type=int, default=600)
    p.add_argument("--rounds-max", type=int, default=8)
    p.add_argument("--passes", type=int, default=2)
    p.add_argument("--jitter", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--note", default="")
    p = sub.add_parser("perturb", help="recovery benchmark: a sibling campaign from a sign-perturbed incumbent")
    p.add_argument("--out", required=True)
    p.add_argument("--flips", type=int, default=8)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--jitter", type=float, default=0.3, help="relative width jitter (0 = signs only)")
    p.add_argument("--ekick", type=float, default=0.02, help="energy kick in MeV (0 = none)")
    p.add_argument("--polish-nfev", type=int, default=60)
    a = ap.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from rmfit.campaign import Campaign

    if a.cmd == "init":
        c = Campaign.init(a.campaign, a.base, sav=a.sav, nshards=a.nshards, note=a.note, sav_fixed=a.sav_fixed,
                          policy=_json_arg(a.policy), search=_json_arg(a.search),
                          flags=(_json_arg(a.flags) if a.flags else None), objective=_json_arg(a.objective),
                          cov_reference=a.cov_reference, spec=a.spec)
        print(open(os.path.join(c.dir, "STATUS.md")).read())
        print("job script:", c.write_job_script())
        return
    if a.cmd == "drive":
        # the build phase has no campaign.json yet; the driver opens the Campaign itself once INIT has run
        from rmfit.drive import main_cli
        sys.exit(main_cli(a))
    c = Campaign(a.campaign, enforce_engine=a.cmd not in READ_ONLY_CMDS)
    try:
        if a.cmd == "spec":
            from rmfit.spec import FitSpec, write_findings_md
            sp = a.spec or os.path.join(c.dir, "spec.json")
            fs = FitSpec.load(sp)
            if a.what == "render":
                out = os.path.join(c.dir, "SPEC.md")
                open(out, "w").write(fs.render_md())
                print("wrote", out)
            elif a.what == "derive":
                print(json.dumps(fs.derive_campaign_cfg(), indent=1))
            else:
                findings = fs.check(a.azr or c.base_azr, campaign_cfg=c.cfg)
                out = os.path.join(c.dir, "SPEC_CHECK.md")
                counts = write_findings_md(out, "Spec check (rmfit spec check)", findings,
                                           [f"spec `{sp}` vs model `{a.azr or c.base_azr}` and campaign.json", ""])
                for f in findings:
                    print(f.line())
                print("wrote", out, counts)
        elif a.cmd == "audit" and (a.kind or a.findings or a.status):
            from rmfit import audit_subagent as asub
            if a.status:
                kind, n, status = a.status
                path = os.path.join(c.dir, "audits", kind, "FINDINGS.md")
                ok = asub.set_status(path, int(n), status, a.note)
                print(("marked" if ok else "NOT FOUND"), kind, n, status, path)
                if ok:
                    c.note(f"audit {kind} finding {n}: {status}" + (f" ({a.note})" if a.note else ""))
            elif a.findings:
                for sm in asub.findings_summary(c.dir):
                    print(f"{sm['kind']}: {sm['n']} findings  {sm['by_severity']}  status {sm['by_status']}  ({sm['path']})")
                    for r in sm["unverified"]:
                        print(f"   #{r['n']} {r['severity']} candidate  {r['where']}: {r['what'][:110]}")
                if not asub.findings_summary(c.dir):
                    print("no audits/*/FINDINGS.md under", c.dir)
            else:
                from rmfit.audit import audit_model
                azr = a.azr or c.base_azr
                machine = audit_model(azr, campaign_dir=c.dir, campaign_cfg=c.cfg, run_target_energies=False)
                kinds = asub.KINDS if a.kind == "all" else (a.kind,)
                for kind in kinds:
                    r = asub.write_subagent_audit(kind, azr, c.dir, c.spec, machine)
                    print(f"{kind}: prompt {r['prompt']}  report {r['report']}")
                print("Hand each PROMPT.md to a background subagent (Agent tool, a cheaper model); it writes FINDINGS.md; "
                      "verify every row against the files, then `rmfit <dir> audit --status KIND N confirmed|rejected|decision`.")
                c.note(f"audit --kind {a.kind}: subagent prompt(s) written under audits/ ({', '.join(kinds)})")
        elif a.cmd == "audit":
            from rmfit.audit import audit_model, write_audit
            azr = a.azr or c.base_azr
            fs = audit_model(azr, campaign_dir=c.dir, campaign_cfg=c.cfg, run_target_energies=not a.no_target_energies)
            if a.json:
                print(json.dumps([f.__dict__ for f in fs], indent=1))
            else:
                for f in sorted(fs, key=lambda f: {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "INFO": 3}.get(f.severity, 9)):
                    print(f.line())
            counts = write_audit(fs, os.path.join(c.dir, "AUDIT.md"), azr)
            print("wrote", os.path.join(c.dir, "AUDIT.md"), counts)
            c.note(f"audit: {', '.join(f'{k} {v}' for k, v in counts.items())} (AUDIT.md)")
        elif a.cmd == "review":
            from rmfit.review import review_campaign
            r = review_campaign(c, grid_check=a.grid_check, figures=a.figures, top=a.top, close_fraction=a.close_fraction)
            for f in r["findings"]:
                print(f.line())
            print("wrote", r["path"], r["counts"])
        elif a.cmd == "decisions":
            items = c.decisions.all() if a.all else c.decisions.pending()
            for d in items:
                print(f"{d.id} [{d['status']}] {d['when']} {d['slug']}: {d.get('question', '')}")
                if d.get("options"):
                    print(f"    options: {d['options']}   default: {d.get('default_if_unanswered', 'wait')}")
                if d.get("answer"):
                    print(f"    answer: {d['answer']} ({d.get('answered_by', '')}, {d.get('answered', '')})")
            if not items:
                print("no pending decisions")
        elif a.cmd == "decide":
            d = c.decisions.answer(a.id, a.choice, by=a.by, note=a.note)
            c.note(f"decision {d.id} ({d['slug']}): {a.choice}" + (f" -- {a.by}" if a.by else "") + (f"; {a.note}" if a.note else ""))
            print(f"{d.id} answered: {a.choice}")
        elif a.cmd == "engine":
            from rmfit import enginepin
            rec = c.cfg.get("engine")
            print("pinned:  ", enginepin.short(rec) if rec else "(none: campaign made before 2026-11-07)")
            print("imported:", enginepin.short(c.engine))
            print("match" if (not rec or rec.get("so_md5") == c.engine.get("so_md5")) else "MISMATCH")
        elif a.cmd == "status":
            st = c.ledger.status()
            if a.json:
                print(json.dumps(st, default=str))
            else:
                c.ledger.write_status_md()
                print(open(os.path.join(c.dir, "STATUS.md")).read())
        elif a.cmd == "polish":
            c.run_polish(max_nfev=a.nfev, shifts_free=(True if a.shifts_free else None), tighten=a.tighten, cid=a.candidate,
                         refreeze_every=a.refreeze_every)
        elif a.cmd == "round":
            if a.top_k:
                c.search.top_k = a.top_k
            if a.local_release is not None:
                c.search.local_release = a.local_release
            stop_file = os.path.join(c.dir, "STOP")
            if os.path.exists(stop_file):
                os.remove(stop_file)
            for i in range(a.rounds):
                st = c.run_sign_round(jpi=a.jpi, note=a.note)
                if st.get("accepted", 0) == 0 and i > 0:
                    c.log("no accepted flips: stopping the round loop")
                    break
                if os.path.exists(stop_file):
                    c.log("STOP file found: stopping the round loop after this round")
                    os.remove(stop_file)
                    break
        elif a.cmd == "crossover":
            c.run_crossover(a.donor, donor_dir=a.donor_dir, note=a.note)
        elif a.cmd == "recall":
            st = c.run_recall(n=a.levels, min_width_kev=a.min_width, levels=a.level, global_nfev=a.nfev, pinned_nfev=a.pinned_nfev,
                              local_release=a.local_release, local_nfev=a.local_nfev, note=a.note)
            print(json.dumps(st, default=str, indent=1))
        elif a.cmd == "population":
            fl = [int(v) for v in str(a.flips).split(",")]
            st = c.run_population(restarts=a.restarts, n_flips=(fl[0] if len(fl) == 1 else tuple(fl)), per_group=a.per_group,
                                  nfev=a.nfev, rounds_max=a.rounds_max, passes=a.passes, jitter=a.jitter, seed=a.seed, note=a.note)
            print(json.dumps(st, default=str, indent=1))
        elif a.cmd == "levelscan":
            grid = None if str(a.grid).lower() == "none" else [float(v) for v in a.grid.split(",")]
            gamma_kev = None if str(a.gamma).lower() == "auto" else float(a.gamma)
            c.run_levelscan(ex=a.ex, top=a.top, gamma_kev=gamma_kev, scales=[float(v) for v in a.scales.split(",")],
                            gamma_bounds_kev=[float(v) for v in a.gamma_bounds.split(",")],
                            seg_match=a.seg_match, rank=a.rank,
                            ex_range=(None if not a.ex_range else [float(v) for v in a.ex_range.split(",")]),
                            jpis=a.jpi, gn_top=a.gn_top, max_signs=a.max_signs, refit=a.refit,
                            ex_offsets=[float(v) for v in a.ex_offsets.split(",")],
                            differential_window_kev=a.differential_window, grid=grid,
                            stage_scales=(None if str(a.stage_scales).lower() == "none"
                                          else [float(v) for v in a.stage_scales.split(",")]),
                            refine_top=a.refine_top, merge_kev=a.merge_kev, note=a.note)
        elif a.cmd == "relocate":
            c.run_relocate(a.level, a.ex, half_window_kev=a.window, max_nfev=a.nfev, note=a.note)
        elif a.cmd == "verify":
            x, keys = c.ledger.load(a.cid)
            cid, sc = c.verify_candidate(x, keys, move=dict(kind="reverify", of=a.cid), parent=a.cid)
            print(f"candidate {a.cid} re-verified as {cid}: {sc.objective:.3f}")
        elif a.cmd == "export":
            c.export(path=a.path)
        elif a.cmd == "uq":
            from rmfit import uq as _uq
            su = dict((getattr(c, "spec", None).get("uq") if getattr(c, "spec", None) is not None else None) or {})
            su.update(json.load(open(a.uq)) if a.uq.endswith(".json") and os.path.exists(a.uq) else json.loads(a.uq))
            if not in_sge_job_env() and not os.environ.get("RMFIT_UQ_LOCAL"):
                from rmfit import site as _site
                print(f"uq needs the engine for hours: {_site.submit_hint('run_crc_rmfit', 'uq ...')} (RMFIT_UQ_LOCAL=1 for a fixture test)")
                return 2
            st = _uq.run(campaign=c, spec_uq=su, log=c.log, nprocs=a.nprocs, threads=a.threads, stop_after=a.stop_after, reset=a.reset)
            g = st["stages"].get("GATES", {})
            print("UQ_STATUS:", "PASS" if g.get("passed") else ("FAIL" if g else "INCOMPLETE"), "->", os.path.join(c.dir, "UQ.md"))
        elif a.cmd == "report":
            print(open(c.report()).read())
        elif a.cmd == "structure" and a.what == "add":
            from rmfit.structure import LevelCandidate
            bank = json.load(open(os.path.join(c.dir, "candidates.json")))["candidates"]
            spec = dict(bank[a.candidate])
            cand = LevelCandidate(energy=spec["energy"], jpi=spec.get("jpi"), channels=spec.get("channels", []),
                                  source=spec.get("source", ""), note=spec.get("note", ""),
                                  energy_fixed=spec.get("energy_fixed", True))
            st = c.run_level_add(cand, jpis=a.jpi, global_nfev=a.nfev, pinned_nfev=a.pinned_nfev,
                                 local_release=a.local_release, local_nfev=a.local_nfev)
            print(json.dumps(st, default=str, indent=1))
        elif a.cmd == "structure" and a.what == "remove":
            st = c.run_remove(top=a.top, nfev=a.nfev)
            print(json.dumps(st, default=str, indent=1))
        elif a.cmd == "structure" and a.what == "dataset-off":
            st = c.run_dataset_off(a.file, nfev=a.nfev)
            print(json.dumps(st, default=str, indent=1))
        elif a.cmd == "log":
            c.note(a.text)
        elif a.cmd == "perturb":
            st = c.perturb(a.out, n_flips=a.flips, seed=a.seed, jitter=a.jitter, e_kick=a.ekick,
                           polish_nfev=a.polish_nfev)
            print(json.dumps(st, default=str, indent=1))
    finally:
        c.close()
        c.ledger.write_status_md()


if __name__ == "__main__":
    main()
