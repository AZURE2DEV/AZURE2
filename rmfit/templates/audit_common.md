## Ground rules (every audit)

- READ-ONLY.  Do not edit any file under `$fit_dir`, `$data_dir` or `$campaign_dir`; do not submit cluster jobs; do not
  run AZURE2 or pyazr fits.  Scripts you write go into `$out_dir/` (name them `audit_*.py`) and read the files only.
- Every finding carries its evidence: file and 1-based line number, segment key (the `<segmentsData>` lines counted
  from 1), the numbers, and how you got them.  Nothing from memory: a number you did not read from a file or a cited
  source is not a finding.
- Findings are CANDIDATES for the driving session, which verifies each one against the files before telling the
  evaluator; data-treatment decisions (drop a row, rescale a set, change a norm) belong to the evaluator.  Do not
  recommend them as done deals; state what the evidence shows and the cheapest check.
- Do not re-report what is already known (section "Already known" below): the machine pre-flight findings, the
  recorded decisions in the spec and the campaign readme.
- Severity: HIGH = a wrong number enters the objective with a large weight (a mis-typed error bar, a wrong column, a
  fraction where percent is read, a wrong frame); MEDIUM = the fit's freedom or scale is affected (free unpenalized
  norm, nominal norm != 1, dead channel, data outside the window); LOW = cosmetic or small; INFO = for the record.
- Write ONE report, `$report`, in the format below, and nothing else outside `$out_dir/`.  Finish with a one-paragraph
  summary at the top (counts by severity, the single most important finding).
- If the harness refuses the report write ("Subagents should return findings as text"), return the COMPLETE report as
  your reply, verbatim in the format below; the driving session saves it to `$report` unchanged.  Do not work around
  the refusal with a shell command.

## Report format (`$report`)

```
# $kind audit: <model or data or system>, $date

<one paragraph: what was read, by what method (scripts in this directory), counts by severity>

## Already known (not re-reported)
<one line>

## Findings

| # | Sev. | Status | Where | What | Evidence | Why it matters | Cheapest check |
|---|------|--------|-------|------|----------|----------------|----------------|
| 1 | HIGH | candidate | file:line / seg key | one sentence | numbers, lines, how obtained | one sentence | one sentence |

## Checked and clean
<the checks that found nothing, one line each -- a clean check is a result>

## To confirm with the evaluator
<intent questions, one line each, phrased as a question with the two readings>
```

`Status` is always `candidate` when you write it; the driving session changes it to `confirmed`, `rejected` or
`decision` (evaluator's call) after verifying.  `rmfit <campaign> audit --findings` reads this table.
