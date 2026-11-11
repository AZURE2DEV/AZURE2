# ratio_frame

A lab-frame RATIO composite is the ratio of the two terms' **lab** cross sections, each term taken to
the lab with its **own** pair's Jacobian.

A lab-frame differential RATIO segment (`isDiff` 1, Advanced Mode "Ratio") compares the data with
base / denominator in the lab.  Each term is evaluated in the c.m. frame of its own pair, σ_lab =
σ_cm / J, so the model is (σ_cm,base / σ_cm,den) (J_den / J_base).  Up to dev 1b6c346
`Segment::CalculateTheoreticalCrossSection` computed both factors with the **base** segment's entrance
and exit pair and the base exit's type (audit 2026-10-25, item 15).  That was right for a denominator of
the same pair, and for another particle pair at a general angle without `relativistic=1` (the
non-relativistic Jacobian uses only the point's own lab and c.m. angles), but wrong for a denominator to
another exit pair with `relativistic=1` (the Jacobian comes from the pair kinematics) or at exactly 0 or
180 deg (the mass formula), and for a photon over a particle or a particle over a photon (the wrong kind
of Jacobian for one of the two).  `Segment::RatioTermFrameFactor` now gives each term its own factor,
and 1 for a ratio that is not a lab-frame differential one (c.m., angle-integrated), as before.

A second, related defect: in mode 3 the lab output (`kinflag.dat` = 1) of a RATIO test segment was
divided by J_base, because `WorkingCopy::MakePoints` converted the test points' cross sections while
`Segment::Fill` skips that for a (dimensionless) ratio data segment.  `MakePoints` skips it too now.

## Models (`make_models.py` (deleted; see "How the models were made"))

- `particle`, `particle_rel` (`relativistic=1`): the `tests/sum_lab_frame` 13C + α model (exits
  n + 16O(6.05, 0+) key 3 and n + 16O(6.13, 3−) key 5, elastic key 2), lab 5.95–6.30 MeV, test
  segments T1 α → n0+ (0–180°), T2 α → n3−, T3 α → α at 60°, T4 α → n3− at 30°, T5 α → n3− at 0°,
  T6/T7 the c.m. versions of T1/T2, T8/T9 angle-integrated; the ratios R1 = T1/T2 (same angles),
  R2 = T1/(α,α)(60°), R3 = T1/(α,n3−)(30°), R4 = T1/(α,n3−)(0°), R5 = T1/(α,n0+)(60°), R6 = c.m.
  T6/T7, R7 = angle-integrated T8/T9.
- `gamma`: the `tests/gamma_lab_frame` 12C(α,γ) external-capture model (E1 + E2 interference), lab
  2–6 MeV: G1 α → γ (10–170°), G2 α → γ at 90°, G3 α → α at 60°, G4 α → α (10–170°); RG1 = G1/G2
  (W(θ)/W(90°)), RG2 = (α,γ)/(α,α)(60°), RG3 = (α,α)/(α,γ)(90°).
- `ratio_data`, `ratio_free`: mode 1 with `relativistic=1`, the RATIO data segment n0+/(α,α)(60°)
  reading `data/ratio.dat`, which `check.sh` writes from the plain segments' lab values; `ratio_free`
  has the seven widths free for the gradient check.

- `fixed_angle`: a term at a fixed angle.  In c.m. data (isDiff 5) the term is evaluated at the fixed angle itself
  (FC = n0+ / (α,α)(60° c.m.)); in lab data (isDiff 1) its 60° lab is converted with the term's pair to
  θ_cm = 60° + asin((m1/m2) sin 60°) = 75.4606°, a number `make_models.py` (deleted; see "How the models were made") computes and the model holds as a plain
  c.m. segment, so FL = n0+_lab / ((α,α)_cm(75.46°) / J) with J = sin²60° / sin²θ_cm cos(θ_cm − 60°) (check (g)).

## Checks (`check.sh path/to/AZURE2`)

- (a), (b): every particle ratio equals the ratio of the plain segments' lab values (c.m. /
  angle-integrated values for R6 / R7) to 1e-9, non-relativistic and relativistic; and the test is
  not vacuous (R1 differs from the c.m. ratio, `relativistic=1` changes R2).
- (c): the three gamma ratios equal the ratio of the plain lab values.
- (d): the c.m. output (`kinflag.dat` = 0) of a lab RATIO test segment prints the same ratio as the
  lab output.
- (e): the RATIO data segment against the plain segments' lab ratio has chi2 ~ 0.
- (f): its analytic gradient (`--use-gradient`, `AZURE_GRAD_CHECK=exit`) equals finite differences.

dev 1b6c346 (the c.m. output, which is not scaled): (a) R1 3.2 % at 0° and 6.5 % at 180°, R4 3.2 %;
(b) R1 6.6 %, R2 3.0 %, R3 2.8 %, R4 3.2 %; (c) RG2 31 %, RG3 1e-4; (d) every lab RATIO output 17–20 %
low, so (a)–(c) fail at the 80 % level in the lab output; (e) chi2 143 over 48 points.  R2, R3, R5, R6,
R7 and RG1 were right already in the c.m. output (R2/R3 only without `relativistic=1`).

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 2640fd8f:tests/ratio_frame/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for tests/ratio_frame (audit item 15, 2026-10-27): a lab-frame RATIO composite is the ratio of the two
terms' LAB cross sections, each term taken to the lab with its OWN pair's Jacobian.

Two reactions, all segments test segments (mode 3; the ratio is read off the lab output, kinflag.dat = 1, and
compared with the ratio of the plain segments' lab values):

  particle      13C + alpha in 17O (the tests/sum_lab_frame model: exits n + 16O(6.05, 0+) key 3 and
                n + 16O(6.13, 3-) key 5, elastic key 2; two levels, no external capture), lab 5.95-6.30 MeV:
                T1 a -> n0+ (0-180 deg), T2 a -> n3- (0-180), T3 a -> a at 60 deg, T4 a -> n3- at 30 deg,
                T5 a -> n3- at 0 deg, T6/T7 the c.m. (isDiff 5) versions of T1/T2, T8/T9 angle-integrated n0+/n3-;
                R1 = T1 / T2 (same angles), R2 = T1 / (a,a)(60), R3 = T1 / (a,n3-)(30), R4 = T1 / (a,n3-)(0),
                R5 = T1 / (a,n0+)(60) (the same pair), R6 = c.m. ratio T6 / T7, R7 = angle-integrated T8 / T9
  particle_rel  the same with relativistic=1
  gamma         12C(a,g)16O external capture (the tests/gamma_lab_frame model: alpha key 1, gamma key 2; E1 + E2
                interference), lab 2-6 MeV: G1 a -> g (10-170 deg), G2 a -> g at 90 deg, G3 a -> a at 60 deg,
                G4 a -> a (10-170); RG1 = G1 / G2 (W(theta)/W(90)), RG2 = G1 / (a,a)(60), RG3 = G4 / (a,g)(90)
  ratio_data    mode 1, relativistic=1: the RATIO data segment n0+ / (a,a)(60 deg) on the particle grid, reading
                data/ratio.dat (written by check.sh from the plain segments' lab values: chi2 must be ~0)
  ratio_free    ratio_data with the seven widths free, for the AZURE_GRAD_CHECK gradient test

Run from this directory: python3 make_models.py  (writes models/).
```
