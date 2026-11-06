# Golden outputs of the present engine for T1 and T3

These files are what the engine of `feature/upgrade` c1cc92c7 (recorded in `ENGINE`) writes for the test projects
`../T1` and `../T3`. Until a phase deliberately changes a result, every later build must reproduce them; the ~96
older test folders stay as oracles alongside them (`docs/upgrade/ARCHITECTURE.md` §10). They are a regression
reference, not a truth: the truth of T1 and T3 is their generator parameters, checked by `../recover.py`.

## What is stored

Per project, three runs from a copy of the project directory, `OMP_NUM_THREADS=2`, default CLI flags (Brune basis,
exact external-capture operator) plus `--no-coulomb-table`, so the exact Coulomb functions are used whether or not
the universal table exists on the machine:

| directory | menu | files |
|---|---|---|
| `calc/` | 1, calculate segments from data at the truth | `chiSquared.out`, `AZUREOut_aa=*_R=*.out` (E_c.m., E_x, θ_c.m., σ or A_y or δ, S factor, data in c.m.), `parameters.out`, `param.par`, `scatteringLengths.out` |
| `extrap/` | 3, the `<segmentsTest>` lines | `AZUREOut_aa=*_R=*.extrap` |
| `rate/` | 5, the rate of the pair in `<lastRun>` (T1 ¹²C(p,γ₀), T3 ⁷Li(p,α)α) at T₉ = 0.01–3 | `reactionrates.out` (cm³ mol⁻¹ s⁻¹) |

## Comparing a build

```bash
python3 tests/upgrade/golden/golden.py --binary /path/to/AZURE2          # or AZURE2_BIN / AZURE2_BINARY
python3 tests/upgrade/golden/golden.py --rtol 2e-5 T1                    # another compiler or platform
```

Every non-numeric token must be equal and every number must satisfy |a − b| ≤ r_tol · max(|a|, |b|) (numbers both
below 10⁻²⁰⁰ count as equal). The default r_tol = 10⁻⁶ is ten times the last printed digit of the `AZUREOut` files
(11 significant digits) and of `reactionrates.out` (6 digits), so on one machine a pass means the results are
unchanged to print precision. Across GCC, Clang and MinGW differences near 10⁻⁵ are known (root `CLAUDE.md`,
"Determinism"); use r_tol = 2 × 10⁻⁵ there. Two runs of the same build agreed to 0 (bit for bit) on 2026-11-06.
The comparison takes about 35 s, most of it the T1 rate.

## Regenerating

Only for an understood, intended change of the physics, in its own commit stating the old and new values and the
independent evidence (`tests/CLAUDE.md`, "Re-baselining"):

```bash
python3 tests/upgrade/golden/golden.py --write --engine-commit <hash of the engine source>
```
