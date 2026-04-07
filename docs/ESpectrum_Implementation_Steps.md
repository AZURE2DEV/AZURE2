# ESpectrum Implementation Steps

Progress tracker for the ESpectrum feature. Work through phases in order; each phase should build and run cleanly before starting the next.

Status legend: `[ ]` not started · `[~]` in progress · `[x]` done

---

## Phase 1 — Scaffolding (no physics, just compiles)

- [ ] **1.1** `include/Config.h` — add `USE_SPECTRUM = (1<<19)` to the `ParameterFlags` enum
- [ ] **1.2** `src/Config.cpp` — parse `useSpectrum` keyword inside `<config>` block; set/clear `USE_SPECTRUM` bit
- [ ] **1.3** `include/ESpectrum.h` — write skeleton header: constructor, `InitializeGrid()`, `Calculate()`, `InterpolateTotal()`, `InterpolateDifferential()`, `ConvolveGaussian()`, `IntegrateTarget()`, `IsReady()`, `GetEntranceKey()`, `GetExitKey()`; all private data members declared
- [ ] **1.4** `src/ESpectrum.cpp` — stub implementations returning 0 / false; `IsReady()` always returns `false` for now
- [ ] **1.5** `CMakeLists.txt` — add `src/ESpectrum.cpp` to the source list
- [ ] **1.6** `include/EData.h` — add `#include "ESpectrum.h"`, add `spectra_` map, declare `InitializeSpectra()`, `CalculateSpectra()`, `GetSpectrum()`
- [ ] **1.7** `src/EData.cpp` — stub implementations for the three new methods (empty bodies)
- [ ] **1.8** `include/ESegment.h` — add `ESpectrum* spectrum_ = nullptr` private member; add `SetSpectrum()` / `GetSpectrum()` public methods
- [ ] **1.9** `src/AZURECalc.cpp` — add the `CalculateSpectra()` call guarded by `USE_SPECTRUM` flag (calls the stub, no effect yet)
- [ ] **1.10** Build and verify: project compiles and runs identically to before with `useSpectrum 0`

---

## Phase 2 — Grid Initialization

- [ ] **2.1** `src/EData.cpp` — implement `InitializeSpectra()`:
  - Walk `segments_`, collect unique `(entranceKey_, exitKey_)` pairs
  - For each pair, compute union energy range across all ESegments sharing that pair
  - Determine max beam sigma and max target depth from TargetEffect objects on those segments
  - Apply range padding (±5% baseline; +3σ_max for convolution; +depth for target integration)
  - Construct `ESpectrum` for each pair, call `InitializeGrid()`; store in `spectra_`
  - Call `SetSpectrum()` on each ESegment to wire the pointer
- [ ] **2.2** `src/EData.cpp` — call `InitializeSpectra()` from `EData::Initialize()` when `USE_SPECTRUM` is set
- [ ] **2.3** `src/ESpectrum.cpp` — implement `InitializeGrid()`:
  - Use `AdaptiveIntegrationGrid` + `GridConfig` to build energy grid from `eMin` to `eMax`
  - Store energies in `gridEnergies_`
  - Allocate `gridTotal_` and `gridLegendre_` with correct sizes (use `Config::maxLOrder`)
  - Set `isReady_ = false` (grid exists but cross sections not yet computed)
- [ ] **2.4** Build and verify: `InitializeSpectra()` runs without crash; correct number of ESpectrum objects created (one per unique pair)

---

## Phase 3 — Cross-Section Calculation on Grid

- [ ] **3.1** `src/ESpectrum.cpp` — implement `Calculate()` using the internal ESegment approach:
  - Create a temporary `ESegment` as an extrapolation-type segment for this pair
  - Populate it with `EPoint` objects at each energy in `gridEnergies_` (angle = 0 for total; full angular set for Legendre)
  - Call the standard `EData::Initialize()` sub-path on the temporary segment (energy-dependent values, Legendre, Coulomb)
  - For each grid point call `EPoint::Calculate(cnuc, cfg)`
  - Store results: `gridTotal_[i]` = angle-integrated cross section; `gridLegendre_[i][l]` = A_L coefficients
  - Set `isReady_ = true`
- [ ] **3.2** `src/EData.cpp` — implement `CalculateSpectra()`: iterate `spectra_` map, call `spectrum.Calculate(cnuc, cfg)` on each
- [ ] **3.3** Build and verify: with `useSpectrum 1`, spectra are computed each iteration; `IsReady()` returns `true`; no change to output (spectrum path not yet used by ESegment)

---

## Phase 4 — Interpolation Queries

- [ ] **4.1** `src/ESpectrum.cpp` — implement `InterpolateTotal()`:
  - Linear (or cubic spline) interpolation of `gridTotal_` at requested energy
  - Return 0 if energy is outside grid range
- [ ] **4.2** `src/ESpectrum.cpp` — implement `InterpolateDifferential()`:
  - Interpolate each `gridLegendre_[l]` coefficient to requested energy
  - Reconstruct σ(E, θ) = Σ_L A_L(E) · P_L(cos θ)
- [ ] **4.3** `src/ESegment.cpp` — add spectrum-path branch in `CalculateTheoreticalCrossSection()`:
  - Guard: `if (spectrum_ && spectrum_->IsReady() && !isAdvanced_)`
  - For points without target effects: call `InterpolateTotal()` or `InterpolateDifferential()` depending on segment type; apply conversion factor
- [ ] **4.4** Build and verify: with `useSpectrum 1`, non-target-effect data points use interpolation from spectrum; results should be close (not identical) to direct calculation; cross-check a known test case

---

## Phase 5 — FFT Gaussian Beam Convolution

- [ ] **5.1** `src/ESpectrum.cpp` — implement `BuildUniformGrid()` and `InterpolateAdaptiveToUniform()`:
  - Determine step `Δ = min(σ_min/10, baseStep)`
  - Interpolate `gridTotal_` and each `gridLegendre_[l]` onto equidistant `uniformEnergies_`
- [ ] **5.2** `src/ESpectrum.cpp` — implement self-contained real FFT utility (Cooley-Tukey, private to ESpectrum.cpp):
  - `FFTForward(std::vector<double>&)` — in-place real FFT
  - `FFTInverse(std::vector<double>&)` — in-place inverse FFT
- [ ] **5.3** `src/ESpectrum.cpp` — implement `ConvolveGaussian()`:
  - Select total or Legendre channel based on `isDifferential`
  - Zero-pad uniform signal to next power of 2
  - FFT → multiply by Gaussian kernel `exp(-2π²σ²k²)` → IFFT
  - Cubic spline sample at each energy in `queryEnergies`
  - If `isDifferential`, convolve each A_L separately, then reconstruct with P_L(cos θ)
- [ ] **5.4** `src/ESegment.cpp` — add convolution-only path in `CalculateViaSpectrum()`:
  - For points with `isConvolution_` TargetEffect and no target integration: call `ConvolveGaussian()`
- [ ] **5.5** Build and verify: enable `useSpectrum 1` on a dataset with Gaussian beam convolution; compare output to quadrature result; agreement should be within interpolation tolerance

---

## Phase 6 — Through-Target Integration

- [ ] **6.1** `src/ESpectrum.cpp` — implement `IntegrateTarget()`:
  - Use stopping power equation from `TargetEffect` to convert depth to energy steps
  - Walk energy grid from query energy downward (energy loss direction)
  - Weight each grid cross section by the Gaussian convolution factor (if convolution also active)
  - Sum with quadrature weights
- [ ] **6.2** `src/ESegment.cpp` — add target-integration path in `CalculateViaSpectrum()`:
  - For `IsTargetIntegration()` effects, call `IntegrateTarget()`
  - For combined convolution + target integration, call target integration first, then convolution
- [ ] **6.3** Build and verify: enable `useSpectrum 1` on a dataset with target integration; compare to existing quadrature result

---

## Phase 7 — GUI Integration

- [ ] **7.1** `gui/include/EditOptionsDialog.h` — add `QCheckBox* useSpectrumCheck` member
- [ ] **7.2** `gui/src/EditOptionsDialog.cpp` — add checkbox to Runtime Options dialog layout:
  - Label: `"Use adaptive spectrum grid for convolution (ESpectrum)"`
  - Initialize checked state from `configure->paramMask & Config::USE_SPECTRUM`
- [ ] **7.3** `gui/src/AZURESetup.cpp` — in `editOptions()`, add flag mapping for `useSpectrumCheck` → `Config::USE_SPECTRUM` (same pattern as existing checkboxes)
- [ ] **7.4** Build GUI and verify: checkbox appears in Runtime Options; toggling it sets/clears the flag correctly; calculation behaves as expected with flag on/off

---

## Phase 8 — Calculation Timing

- [ ] **8.1** `gui/src/AZUREMainThread.cpp` — add `#include <chrono>` and `#include <iomanip>`
- [ ] **8.2** `gui/src/AZUREMainThread.cpp` — in `AZUREMainThreadWorker::run()`, wrap `azureMain_()`:
  ```cpp
  auto t0 = std::chrono::steady_clock::now();
  stream_ << "Calculation started.\n";
  azureMain_();
  auto t1   = std::chrono::steady_clock::now();
  double dt = std::chrono::duration<double>(t1 - t0).count();
  if (dt < 60.0)
      stream_ << "Calculation finished in " << std::fixed
              << std::setprecision(2) << dt << " s.\n";
  else
      stream_ << "Calculation finished in "
              << static_cast<int>(dt / 60) << " min "
              << std::setprecision(2) << std::fmod(dt, 60.0) << " s.\n";
  ```
- [ ] **8.3** Build and verify: GUI log shows start and finish messages with elapsed time for all calculation types (fit, output, error analysis)

---

## Phase 9 — Validation and Cleanup

- [ ] **9.1** Run full test suite with `useSpectrum 0` — all results must be bit-identical to pre-feature baseline
- [ ] **9.2** Run with `useSpectrum 1` on a dataset with no target effects — verify cross sections match direct calculation within 0.1%
- [ ] **9.3** Run with `useSpectrum 1` on a dataset with Gaussian beam convolution — verify agreement with quadrature within tolerance; benchmark speedup
- [ ] **9.4** Run with `useSpectrum 1` on a dataset with through-target integration — verify agreement with quadrature within tolerance
- [ ] **9.5** Run multi-threaded fit with `useSpectrum 1` — verify no race conditions, results reproducible across runs
- [ ] **9.6** Verify ESpectrum is rebuilt correctly after parameter changes (Minuit iterations produce consistent chi-squared progression)
- [ ] **9.7** Check memory: no leaks when EData is cloned into the object pool and returned
- [ ] **9.8** Update `CLAUDE.md` if any new configuration keywords or architecture notes are needed

---

## Notes

- Phases 1–3 are the minimum needed to have the feature "on" without breaking anything.
- Phases 4–6 replace the existing quadrature paths progressively; each can be validated independently.
- Phase 7 and Phase 8 are independent of each other and can be done in either order after Phase 6.
- If angular distribution support turns out to be complex, Phase 4 can be split: implement `InterpolateTotal()` only first, gate `InterpolateDifferential()` behind a secondary check, and add it in a sub-step.
