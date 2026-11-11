# relativistic -- relativistic two-body kinematics (`relativistic=1` in `<config>`)

`include/Kinematics.h`: with the key, every channel's wave number and Sommerfeld parameter come from the
invariant mass s, k^2 = (s - m^2)(s - (m1 - m2)^2)/(4s), eta = Z1 Z2 alpha (s - m1^2 - m2^2)/(2 k sqrt(s)),
a point's c.m. kinetic energy from s = m^2 + 2 m2 E_lab, and particle lab angles and cross sections from the
exact Lorentz boost (the LANL/EDA prescription, ENDF KRL=1, LA-UR-19-31084). Level energies stay excitation
energies; `relativistic.out` lists them in EDA's E(s) convention.

`make_models.py` (deleted; see "How the models were made") writes the models (committed), each non-relativistic and with `relativistic=1`:

| model | what |
|---|---|
| `na_hs` | n + 4He, s wave only, hard sphere: sigma = 4 pi sin^2(ka)/k^2 |
| `pa_hs` | p + 4He, s wave only, hard sphere + Coulomb at 90/150 deg |
| `he5` | 5He, d + t and n + 4He with the 3/2+ level; t(d,n)a and 4He(n,d)t at equal excitation energy |
| `he5_ang`, `pa_ang` | lab-frame differential data of value 1 (the c.m. value is dOmega_lab/dOmega_cm) |

`check.sh path/to/AZURE2 [reference AZURE2]` runs them (`--no-coulomb-table`) and `verify.py` checks, against
independent formulas: (a) the c.m. energies of the lab grid and the hard-sphere cross section to 1e-9;
(b) |f_C + f_N|^2 from mpmath Coulomb functions at the relativistic eta and rho (skipped without mpmath);
(c) detailed balance g k^2 sigma between the two partitions to 1e-10; (d) the numbers in relativistic.out;
(e) the relativistic c.m. angles and Jacobians against an explicit four-vector boost (2e-7); (g) the hybrid
Coulomb functions are refused; and, given a reference binary, (f) every non-relativistic model's output is
byte-identical to it. About 20 s.

Sizes seen here: at E_lab 30 MeV n + 4He has E_cm 61 keV lower and k 0.28 % larger; the c.m. cross section of
p + 4He elastic at 20 MeV moves by up to 0.7 % and its c.m. angle by 0.19 deg.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show ec4b2faa:tests/relativistic/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models of tests/relativistic (run once; the .azr files are committed).

na_hs   n + 4He, s wave only (includeAllJPi=1, maxOrbitalMomentum=0) and one far 1/2+ level of zero width:
        pure hard-sphere scattering, sigma = 4 pi sin^2(k a) / k^2, a test of k and pi/k^2.
pa_hs   p + 4He, s wave only, the same: differential elastic at 90 and 150 deg, a test of eta and k through
        the Coulomb functions and phases.
he5     5He: d + t and n + 4He (Q = 17.589 MeV) with the 3/2+ level near the d+t threshold; t(d,n)a and its
        inverse 4He(n,d)t on test grids of equal invariant mass, for detailed balance.
Each is written twice: <name>.azr (non-relativistic) and <name>_rel.azr (relativistic=1).
```
