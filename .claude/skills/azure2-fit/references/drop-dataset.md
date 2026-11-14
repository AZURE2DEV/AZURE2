# Dropping a data set from a fit, and bringing it back

Removing one or more data sets from a fit (leave-one-out, held-out data, uncertainty runs, a set judged unusable)
must change exactly that set's contribution to the objective and nothing else. This page gives the rule
(deactivate, never delete), where the change must be recorded, how to make it, and the check that proves it.

## The rule: deactivate, never delete

Set the data set's `active` field to `no` and leave the record where it is:

```python
p = pyazr.Project.open("model.azr")
p.set("s04_thick", active="no")
p.save("model_no_s04.azr")              # in the same directory: the data paths stay relative
```

A deactivated record keeps its name, its parameters' names (`s04_thick.norm`) and its place in every table, so
an existing parameter file stays valid (its lines for the inactive set match nothing and are ignored) and
re-activation is a one-word change. A record that is deleted takes its own parameters with it; a group
(`normgroup`, `shiftgroup`) whose first member is deleted hands its shared value to the next member.

## The drop must reach the model being run

A fit directory accumulates variant projects (uncertainty model, extrapolation model, campaign seed). Deactivating
a data set in one of them does nothing for the others, and a variant built later from the master inherits the
master's field. Before a long run, list which projects have it active:

```bash
for f in *.azr */*.azr; do python3 -c "import pyazr, sys
print(sys.argv[1], pyazr.Project.open(sys.argv[1]).describe(sys.argv[2]).splitlines()[0])" "$f" s04_thick; done
```

Decide where the drop is recorded as the fit's state: a property of the analysis ("this set is not used") goes in
the master project and its readme; a property of one run goes in that run's project only, and the readme says the
master still has it on.

## Verify: the difference is exactly the dropped data set

Evaluate both projects at the same values (the fitted project written by `save_fit` carries them):

```python
with pyazr.Session("model.azr") as a, pyazr.Session("model_no_s04.azr") as b:
    ea, eb = a.evaluate(a.objective().x0), b.evaluate(b.objective().x0)
    d = ea.dataset("s04_thick")
    print(ea.objective - eb.objective, d.chi2 + d.norm_chi2 + d.shift_chi2)
```

Pass: the two numbers agree to rounding, and no remaining data set's `chi2` or `norm` moves. If the dropped set
shares a data covariance with sets that stay, the remaining sets' χ² changes too (the whitening joined them): the
check is then on the sum over the joined group, and the covariance file must still list the dropped set's data
file (`azure2-uncertainty`, data covariances). External-capture caches (`output/intEC.*`) depend on the active
capture sets: delete them when a capture set changes.

## The refit afterwards

The first objective the refit reports must equal the verified one; if not, it did not start where you think. Use
the same options (`<config>` keys of the project, `docs/source/reference/command_line.rst`) as the reference run.

## Record it

In the fit's readme: which data set, why it was dropped, in which projects it is off (and in which deliberately
still on), and the verification numbers.

## Bringing a data set back

Set `active` to `yes` in the same place. Its parameters return with the values in the record (or the parameter
file). Verify in reverse: the objective rises by that set's own contribution at those values.
