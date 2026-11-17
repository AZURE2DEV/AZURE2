"""``rmfit audit --kind input|data|literature``: the three idle-time reviews of the rmfit-campaign skill (section 3b)
as templated, self-contained subagent prompts, plus the parser of the FINDINGS.md they write.

rmfit cannot start a Claude subagent itself; it writes the prompt the driving session hands to one (Agent tool, a
cheaper model, background), with everything the subagent needs and nothing it should re-report: the model's segment
table, the active data files, the machine pre-flight findings (AUDIT.md), the spec's data-set roles and standing
decisions, and the head of the campaign readme.  The subagent writes ``audits/<kind>/FINDINGS.md`` in a fixed table
whose ``Status`` column the driving session fills after verifying each row (candidate -> confirmed | rejected |
decision); ``rmfit <campaign> audit --findings`` lists what is still unverified.  Later ``rmfit drive`` (DESIGN.md,
Phase 2) emits these prompts from its AUDIT state and waits for the reports.

Why (DESIGN.md, finding 2.5): input-data errors were found by hand audits, never by the
fit; the audits work (42 ``*_audit`` directories, 16-minute turnaround on a large model) but each prompt was written by hand
and each report verified by hand.  Worked example and format source: an archived campaign.
"""
import os
import re
import string
import time
from typing import List, Optional

from .spec import Finding, pair_lines
from .backend import PYAZR3

KINDS = ("input", "data", "literature")
STATUSES = ("candidate", "confirmed", "rejected", "decision")
TEMPLATES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
PYTHONPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + ":" + os.environ.get("AZURE2_ROOT", os.path.expanduser("~/AZURE2"))


def _reaction_name(azr_path):
    """``<A><El>+<proj>`` from the archive path convention (``<A><El>+<proj>/...``), else the model's file stem."""
    for part in os.path.abspath(azr_path).split(os.sep):
        if re.match(r"^\d+[A-Z][a-z]?\+[A-Za-z0-9]+$", part):
            return part
    return os.path.splitext(os.path.basename(azr_path))[0]


def _pairs_text(model):
    rows = []
    try:
        cn = model.nucleus_of_pair(1)
    except Exception:
        cn = "?"
    for key, c in pair_lines(model).items():
        if c.ptype == 10 or getattr(c, "is_photon", False):
            rows.append(f"- pair {key}: gamma + {cn}" + (f" (Ex {c.excitation:.4f})" if getattr(c, 'excitation', None) else ""))
            continue
        try:
            light = f"Z={c.Z1} M={c.M1:g} J={c.J1:g}"
            heavy = f"Z={c.Z2} M={c.M2:g} J={c.J2:g}"
            thr = c.threshold if PYAZR3 else c.entrance_sep_energy - c.sep_energy + c.excitation      # Ex of the pair's threshold
            rows.append(f"- pair {key}: ({light}) + ({heavy}), separation energy {c.sep_energy:.4f} MeV, Ex(heavy) = {c.excitation:.4f}, "
                        f"threshold at Ex = {thr:.4f} MeV, a = {c.channel_radius:g} fm")
        except Exception as err:
            rows.append(f"- pair {key}: {err}")
    return f"Compound nucleus {cn}.\n" + "\n".join(rows)


def _segment_table(segs, base):
    """The data-set table of a subagent prompt (pyazr 3 / pyazr 2 versions below)."""
    return (_segment_table_p3(segs, base) if PYAZR3 else _segment_table_p2(segs, base))

def _segment_table_p2(segs, base):
    head = ("| key | active | in->out | observable | file | E window (lab) | angle | norm (nominal, vary, err %) | shift | cov | group |\n"
            "|---|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for s in segs:
        ang = f"{s.angle_min:g}-{s.angle_max:g}" if s.observable == "differential" else "-"
        norm = f"{s.norm:g}, {'free' if s.vary_norm else 'fixed'}, {s.norm_error:g}"
        shift = f"{s.energy_shift:g}, {'free' if s.vary_shift else 'fixed'}, {s.energy_shift_error:g}" if (s.vary_shift or s.energy_shift) else "-"
        obs = s.observable + (" SUM" if s.total else "") + (" composite" if s.composite else "")
        rows.append(f"| {s.key} | {'1' if s.active else '0'} | {s.entrance_key}->{s.exit_key} | {obs} | {s.data_file} | "
                    f"{s.energy_min:g}-{s.energy_max:g} | {ang} | {norm} | {shift} | {s.covariance or '-'} | {s.norm_group or '-'} |")
    return head + "\n" + "\n".join(rows)

def _segment_table_p3(segs, base):
    head = ("| key | name | active | reaction (pairs in->out) | quantity | file | E window (lab) | angle | norm (as written) | "
            "shift (as written) | cov | shared norm |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for s in segs:
        ang = f"{s.angle_min:g}-{s.angle_max:g}" if s.differential else "-"
        obs = s.observable + (f" ({s.operation} of {len(s.terms)} reactions)" if s.composite else "")
        rows.append(f"| {s.key} | {s.name} | {'1' if s.active else '0'} | {s.tree} ({s.entrance_key}->{s.exit_key}) | {obs} | "
                    f"{s.data_file} | {s.energy_min:g}-{s.energy_max:g} | {ang} | {s.norm_text or '-'} | {s.shift_text or '-'} | "
                    f"{s.covariance or '-'} | {s.norm_group or '-'} |")
    return head + "\n" + "\n".join(rows)


def _data_files(segs, base):
    files = sorted({s.data_file for s in segs if s.active})
    rows = []
    for f in files:
        p = f if os.path.isabs(f) else os.path.join(base, f)
        keys = [str(s.key) for s in segs if s.active and s.data_file == f]
        if os.path.exists(p):
            n = sum(1 for line in open(p, errors="replace") if line.strip() and not line.lstrip().startswith("#"))
            rows.append(f"- `{p}`: {n} rows; segments {', '.join(keys)}")
        else:
            rows.append(f"- `{p}`: MISSING; segments {', '.join(keys)}")
    return "\n".join(rows) or "(no active segment)"


def _known(campaign_dir, spec, machine_findings: Optional[List[Finding]]):
    out = []
    if machine_findings:
        fs = [f for f in machine_findings if f.severity in ("HIGH", "MEDIUM")]
        if fs:
            out.append("Machine pre-flight (`rmfit audit`, AUDIT.md) already reports:")
            out += ["  " + f.line()[2:] for f in fs]
    if spec is not None:
        roles = []
        for ds in spec.data_sets:
            bits = [f"role {ds.get('role', '?')}"]
            if ds.get("norm"):
                bits.append("norm " + ", ".join(f"{k} {v}" for k, v in ds["norm"].items()))
            if ds.get("excluded_rows"):
                bits.append(f"excluded rows {ds['excluded_rows']}")
            if ds.get("window"):
                bits.append(f"window {ds['window']}")
            if ds.get("note"):
                bits.append(str(ds["note"]))
            roles.append(f"  - {ds.get('name') or ds.get('file')}: " + "; ".join(bits))
        if roles:
            out.append("Data-set roles and treatments recorded in spec.json (decided by the evaluator):")
            out += roles
        sd = spec.get("standing_decisions") or []
        if sd:
            out.append("Standing decisions (spec.json):")
            out += [f"  - {d.get('date', '')} {d.get('who', '')}: {d.get('text', '')}" for d in sd]
    dec = os.path.join(campaign_dir or "", "DECISIONS.md")
    if campaign_dir and os.path.exists(dec):
        try:
            from .decisions import DecisionQueue
            ans = [d for d in DecisionQueue(dec).all() if d["status"] == "answered"]
            if ans:
                out.append("Answered decisions (DECISIONS.md):")
                out += [f"  - {d.id} {d['slug']}: {d.get('answer', '')} ({d.get('answered_by', '')})" for d in ans]
        except Exception:
            pass
    return "\n".join(out) if out else "(nothing recorded yet)"


def _readme_head(campaign_dir, azr_path, n=14):
    for d in (campaign_dir, os.path.dirname(azr_path)):
        if not d:
            continue
        for name in ("readme", "readme.txt", "README.md"):
            p = os.path.join(d, name)
            if os.path.exists(p):
                lines = open(p, errors="replace").read().splitlines()[:n]
                return f"`{p}`:\n\n" + "\n".join("    " + ln for ln in lines)
    return "(no readme found)"


def subagent_prompt(kind: str, azr_path: str, campaign_dir: Optional[str] = None, spec=None,
                    machine_findings: Optional[List[Finding]] = None, out_dir: Optional[str] = None) -> str:
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    from .backend import AzrModel, SegmentSet
    azr_path = os.path.abspath(azr_path)
    base = os.path.dirname(azr_path)
    segs = SegmentSet.from_file(azr_path)
    model = AzrModel.from_file(azr_path)
    out_dir = out_dir or os.path.join(campaign_dir or base, "audits", kind)
    data_dir = os.path.dirname(os.path.join(base, next((s.data_file for s in segs if s.active), "data/x")))
    fields = dict(
        kind=kind, date=time.strftime("%Y-%m-%d"), azr=azr_path, fit_dir=base, data_dir=os.path.normpath(data_dir),
        campaign_dir=campaign_dir or base, out_dir=out_dir, report=os.path.join(out_dir, "FINDINGS.md"),
        reaction=_reaction_name(azr_path), pythonpath=PYTHONPATH,
        segment_table=_segment_table(segs, base), data_files=_data_files(segs, base), pairs=_pairs_text(model),
        known=_known(campaign_dir, spec, machine_findings), readme_head=_readme_head(campaign_dir, azr_path))
    common = string.Template(open(os.path.join(TEMPLATES, "audit_common.md")).read()).safe_substitute(fields)
    fields["common"] = common
    return string.Template(open(os.path.join(TEMPLATES, f"audit_{kind}.md")).read()).safe_substitute(fields)


def write_subagent_audit(kind: str, azr_path: str, campaign_dir: Optional[str] = None, spec=None,
                         machine_findings: Optional[List[Finding]] = None) -> dict:
    """Write ``<campaign>/audits/<kind>/PROMPT.md`` and an empty ``FINDINGS.md`` skeleton; return their paths.
    An existing FINDINGS.md with rows is never overwritten (the prompt is)."""
    base = os.path.dirname(os.path.abspath(azr_path))
    out_dir = os.path.join(campaign_dir or base, "audits", kind)
    os.makedirs(out_dir, exist_ok=True)
    text = subagent_prompt(kind, azr_path, campaign_dir, spec, machine_findings, out_dir)
    prompt = os.path.join(out_dir, "PROMPT.md")
    open(prompt, "w").write(text)
    report = os.path.join(out_dir, "FINDINGS.md")
    if not (os.path.exists(report) and parse_findings(report)):
        open(report, "w").write(f"# {kind} audit: (not yet written; the subagent replaces this file)\n\n"
                                "## Findings\n\n| # | Sev. | Status | Where | What | Evidence | Why it matters | Cheapest check |\n"
                                "|---|------|--------|-------|------|----------|----------------|----------------|\n")
    return dict(kind=kind, prompt=prompt, report=report, out_dir=out_dir)


_SEV = ("HIGH", "MEDIUM", "LOW", "INFO")


def parse_findings(path: str) -> List[dict]:
    """Rows of the findings table of a FINDINGS.md: ``{n, severity, status, where, what, evidence, why, check}``.
    Tolerates bold severities, missing Status column (-> candidate) and extra columns."""
    rows = []
    if not os.path.exists(path):
        return rows
    for line in open(path, errors="replace"):
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4 or not re.match(r"^\d+$", cells[0]):
            continue
        sev = cells[1].strip("*` ").upper()
        if sev not in _SEV:
            continue
        status = (cells[2].strip("*` ").lower().split() or [""])[0]      # "confirmed (note)" -> confirmed
        if status in STATUSES:
            rest = cells[3:]
        else:
            status, rest = "candidate", cells[2:]
        rest += [""] * (5 - len(rest))
        rows.append(dict(n=int(cells[0]), severity=sev, status=status, where=rest[0], what=rest[1],
                         evidence=rest[2], why=rest[3], check=rest[4]))
    return rows


def set_status(path: str, n: int, status: str, note: str = "") -> bool:
    """Change the Status cell of finding ``n`` in a FINDINGS.md (the driving session's verification mark)."""
    if status not in STATUSES:
        raise ValueError(f"status must be one of {STATUSES}")
    lines = open(path, errors="replace").read().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if not line.startswith("|"):
            continue
        cells = line.rstrip("\n").split("|")
        if len(cells) < 5 or cells[1].strip() != str(n):
            continue
        sev = cells[2].strip().strip("*` ").upper()
        if sev not in _SEV:
            continue
        cur = (cells[3].strip().strip("*` ").lower().split() or [""])[0]
        if cur in STATUSES:
            cells[3] = f" {status}" + (f" ({note})" if note else "") + " "
        else:
            cells.insert(3, f" {status}" + (f" ({note})" if note else "") + " ")
        lines[i] = "|".join(cells) + "\n"
        open(path, "w").write("".join(lines))
        return True
    return False


def findings_summary(campaign_dir: str) -> List[dict]:
    """One dict per audit kind present under ``<campaign>/audits/``: counts by severity and status."""
    out = []
    for kind in KINDS:
        p = os.path.join(campaign_dir, "audits", kind, "FINDINGS.md")
        if not os.path.exists(p):
            continue
        rows = parse_findings(p)
        out.append(dict(kind=kind, path=p, n=len(rows),
                        by_severity={s: sum(r["severity"] == s for r in rows) for s in _SEV},
                        by_status={s: sum(r["status"] == s for r in rows) for s in STATUSES},
                        unverified=[r for r in rows if r["status"] == "candidate"]))
    return out
