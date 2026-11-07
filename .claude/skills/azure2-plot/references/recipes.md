# Figure recipes

Each recipe computes what AZURE2 compares and leaves the styling to the prc-figures skill. All quantities are c.m.
unless stated. `x` is the free rwa vector of the fit (`m.load_sav("output/param.sav")`). Check every curve against the
fit column of `AZUREOut_*.out` (SKILL.md, "Rules for every figure").

## Excitation function with pulls (data mode)

The pull of point i is r_i = (c_i − d_i n)/(e_i n): c the calculation, d and e the c.m. data and uncertainty, n the
segment's norm. Σ r_i² over the segment is its χ² (without a data covariance; with a `covariance <file>` the χ² uses the
full matrix and the diagonal pulls are only indicative).

```python
import numpy as np
from pyazr import azure2

with azure2("project.azr", cwd=HERE) as m:
    x     = m.load_sav("output/param.sav")
    calc  = m.calculate_rwa(x)            # c per segment (b or b/sr)
    ecm   = m.calculate_energies(x)       # E_c.m. (MeV), including the fitted energy shift
    norms = m.segment_norms(x)            # n per segment at x
    i = 3                                 # calculated-segment index (data mode: m.datasets[i])
    c, d, e, n = (np.asarray(calc[i]), np.asarray(m.cross[i]),
                  np.asarray(m.cross_err[i]), norms[i])
    pull = (c - d * n) / (e * n)
    # upper panel: errorbar(ecm[i], d*n, e*n) and the curve; lower panel: pull against ecm[i]
```

- For a smooth curve between data points, evaluate an extrapolation grid on the same pair and angle (azure2-calculate,
  Extrapolation) with `<targetInt>` matching the data segment's: folded if the data are folded.
- Several segments of one data set: concatenate after sorting by energy; label by data set, not by segment.
- On an E_x axis use `m.calculate_excitation_energy(x)` (includes the entrance pair's E*).

## Angular distribution

- Data and fit at the measured angles: `m.calculate_angles(x)[i]` (θ_c.m., deg) with `calc[i]` and the data as above,
  for each segment of a fixed energy. The `.azr` holds lab angles; the arrays are c.m.
- A smooth curve: a `differential-cm` test segment per angle grid point, or Legendre coefficients from
  `pyazr.angular_distribution(azr_file, energies, entrance, exit, order, params=x)`, which returns
  `(energies_cm, a_k)` with W(θ) = Σ_k a_k P_k(cos θ); pass all energies in one call (it starts an engine per call).
- Analyzing power: columns hold A_y (dimensionless, signed); draw it on a linear axis.

## S factor with band and decomposition

```python
from pyazr import extrapolation_bands
bands = extrapolation_bands("project.azr", keys=[k], quantity="sfactor")   # k: 1-based <segmentsTest> number
b = bands[k]
lo, hi = b.interval(1.0)                      # 1σ band (MeV b)
# curve: b.energy (E_c.m., MeV) against b.value; band: fill_between(b.energy, lo, hi)
```

- Data S: `m.sfactor[i]`, `m.sfactor_err[i]` (data mode), multiplied by the segment norm as for σ.
- The band is linear propagation of `output/covariance.dat` at `param.sav`; its meaning and alternatives
  (profile likelihood, samples): azure2-uncertainty.
- Decomposition curves (one level, its interference, external capture): build them as cross sections
  (azure2-calculate, "Decomposing a cross section") and convert with `conv = S_full/σ_full` on the same grid.
- S(0): read it off the curve at the lowest grid energy, or fit a quadratic through the three lowest points
  (`pyazr/examples/sfactor_extrapolation.py`); state which, and keep the lowest point well below any resonance.
- Use E_c.m. in keV or MeV on a linear axis and S in keV b on a log axis when it spans decades.

## Rate ratio

The ratio of a rate to a reference rate (a compilation, or a previous evaluation) on one T9 grid:

1. Compute N_A⟨σv⟩ with mode 5 (`reactionrates.out`: T9 in GK, rate in cm³ mol⁻¹ s⁻¹) or with
   `pyazr/examples/reaction_rate.py` per sample.
2. Interpolate the reference to the same T9 values in log T9 and log rate; never interpolate across a resonance
   threshold with a linear rule.
3. Plot rate/reference against T9 on a log axis, with the band from the rate samples or from propagation
   (azure2-uncertainty) divided by the same reference.
4. If the ratio has structure, plot `cross_section_grid.dat` (E_c.m., σ, S at every node of the rate integral) first: a
   missing narrow level or a grid gap shows there.

## Check before saving

- The curve agrees with column 4 of `AZUREOut_*.out` at the data energies (relative difference at the 10⁻⁵ level for
  the same build and flags).
- No isolated spikes at channel thresholds; grid spacing below a quarter of the narrowest width shown.
- Axis labels carry frame and unit; the caption says folded or bare, parameter source and band definition.
