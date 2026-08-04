# 13N — the reference evaluation

The `12C + p` system: capture, elastic scattering, and the vector analyzing
power. Small enough to run in seconds, but it exercises most of the physics —
particle and photon channels, external capture, target integration, energy
shifts, angular distributions, and polarization.

`run_tests.sh` runs a plain calculation and compares the total and per-segment
chi-squared against `expected/chiSquared.out`.

## Segments

| # | data | observable |
|---|---|---|
| 1–2 | Artemov, Seagrave | angle-integrated capture |
| 3–5 | Meyer, at 84.3°, 114.5°, 144.1° | differential |
| 6–8 | Skowronski (LUNA HPGe, LUNA BGO, FELS) | angle-integrated |
| 9–10 | Kettner, at 0° and 55° | differential |
| 11–16 | Baumann, at six energies | **analyzing power** |

Segments 1–10 are the capture and scattering evaluation. Segments 11–16 were
added to cover the analyzing power (observable 7) against measured data; they
are the same compound nucleus and the same particle pairs, which is why they
live here rather than in an evaluation of their own.

## The analyzing-power data

R. Baumann, G. Keil, N. Kniest, E. Pfaff, M. Preiss, M. Skill and G. Clausnitzer,
*The analyzing power for elastic 12C(p,p)12C scattering below 2.1 MeV*,
Nucl. Phys. **A542** (1992) 53.

The paper publishes no table of A_y — its results are contour plots (figs. 3
and 4) plus six angular distributions in fig. 1, at E_p = 1.618, 1.658, 1.708,
1.738, 1.758 and 1.779 MeV. `data/baumann_ay.dat` is fig. 1 **digitised**: the
page rendered at 600 dpi, the axes calibrated on the printed tick labels, and
the drawn curve traced column by column. `digitize_fig1.py` regenerates it from
the PDF, which is not redistributed here.

Two things this data is not:

* What is traced is the **fitted curve** of the paper's phase-shift analysis,
  which the plotted points lie on. These are the published analysis sampled at
  10-degree intervals, not the raw measurements.
* The quoted uncertainty of 0.04 is **digitisation error, not the
  experiment's**. Baumann quote statistical errors below 0.004 — ten times
  smaller. 0.04 is roughly the drawn line width (12 px against 466 px per unit
  of A_y) plus the axis calibration residual.

Points where the automatic trace could not follow a steep section were dropped
rather than interpolated, so 71 of a possible 78 survive. One segment per energy
keeps each angular distribution separately plottable, as in the paper's figure.

Angles are centre-of-mass, as printed on the figure axis, which is what
observable 7 expects; energies are laboratory.

## Reproducing the analyzing-power fit

The evaluation as shipped is fitted to everything, so its levels sit at the
values the capture and scattering data want. Fitting A_y *alone* is a separate
exercise and a useful check that the polarization observable carries real
information:

```python
import re, pathlib
s = pathlib.Path("13N.azr").read_text()

# every level fixed except the two resonances this energy range sees
def levels(b):
    out = []
    for line in b.splitlines():
        f = line.split()
        if len(f) >= 12:
            free = f[2] in ("3.5032", "3.5453")
            f[3]  = "0" if free else "1"
            f[10] = "0" if (free and f[5] == "1") else "1"
            line = "  ".join(f)
        out.append(line)
    return "\n".join(out)
s = re.sub(r'(<levels>\n)(.*?)(</levels>)',
           lambda m: m.group(1) + levels(m.group(2)) + "\n" + m.group(3), s, flags=re.S)

# only the analyzing-power segments active
def segs(b):
    out = []
    for line in b.splitlines():
        if not line.strip():
            continue
        f = line.split()
        f[0] = "1" if f[7] == "7" else "0"
        out.append("  ".join(f))
    return "\n".join(out)
s = re.sub(r'(<segmentsData>\n)(.*?)(</segmentsData>)',
           lambda m: m.group(1) + segs(m.group(2)) + "\n" + m.group(3), s, flags=re.S)
pathlib.Path("_ay_only.azr").write_text(s)
```

Then `printf '2\n\nn\n\n' | AZURE2 --no-gui --no-readline _ay_only.azr`.

Chi-squared goes from 236.8 to 85.2 — 1.20 per point, 1.27 per degree of
freedom with four free parameters:

| parameter | fitted | Baumann table 1 | start |
|---|---|---|---|
| 3/2⁻ E_x | 3.4954 MeV | 3.499 MeV | 3.5032 MeV |
| 3/2⁻ Γ_p | 50.5 keV | 57 keV | 55.2 keV |
| 5/2⁺ E_x | 3.5430 MeV | 3.546 MeV | 3.5453 MeV |
| 5/2⁺ Γ_p | 53.6 keV | 50 keV | 49.0 keV |

Excitation energies agree to 3–4 keV and widths to about 10%, which is what
data this loose can support: the 0.04 uncertainty is an order of magnitude
looser than the measurement, so the fit cannot be expected to recover the
published parameters more tightly. Restarting the fit from the minimum returns
85.226 rather than drifting.

## Two traps this evaluation illustrates

**`<targetInt>` assigns target effects by segment index**, and `<segmentsTest>`
shares that numbering with `<segmentsData>`. The analyzing-power segments were
appended after the ten existing ones precisely so they would not inherit the
gas-target integration declared for segments 1 and 2. That matters more here
than it would for a cross section: A_y averaged over a target is weighted by
the cross section, and since Rutherford scattering diverges at low energy where
A_y is essentially zero, a thick target drives the average to about 1e-6.

**`<parameterSettings>` is also keyed by segment number.** This file previously
carried entries for segments 11 through 18, left over from a configuration that
had more segments than it does now. They were harmless only until a segment 11
existed — at which point a new data set silently inherited another segment's
nuisance prior and a freed normalization. They have been removed, and AZURE2
now drops settings entries numbered beyond the segments that actually exist.
