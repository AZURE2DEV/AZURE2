# ESpectrum Implementation Plan

## Overview

ESpectrum is an **internal pre-computed cross-section cache** on an adaptive energy grid, one instance per unique (entrance, exit) particle pair. It acts as an "internal extrapolation segment" — like the existing `<segmentsTest>` extrapolation blocks, but auto-generated, hidden from output, and designed to serve FFT-based convolutions.

When `USE_SPECTRUM` is active, every Minuit iteration recomputes each ESpectrum first, and ESegments with target effects use FFT convolution against the spectrum instead of spawning per-point quadrature sub-points. ESegments without target effects can use the spectrum for fast cross-section interpolation.

---

## Architecture

```
EData
 ├── segments_[]           (existing data segments)
 ├── componentSegments_[]  (existing)
 └── spectra_              map<(entranceKey,exitKey), ESpectrum>  ← NEW
       └── ESpectrum
             ├── adaptiveGrid_[]     grid points on adaptive energy grid
             ├── uniformGrid_[]      equidistant grid for FFT
             └── fftConvolved_[]     pre-convolved result per σ_beam

ESegment
 ├── points_[]             (existing)
 └── spectrum_*            pointer into EData::spectra_  ← NEW (nullable)
```

---

## Design Decisions

### Angular Distribution Handling

ESpectrum must support both total (angle-integrated) and differential (angle-specific) cross sections. When an EPoint queries the spectrum it provides its type:

- **Total / angle-integrated:** ESpectrum stores σ_total(E) on the grid. The query returns the interpolated total cross section.
- **Differential:** ESpectrum stores Legendre expansion coefficients A_L(E) on the grid. The query reconstructs σ(E, θ) = Σ_L A_L(E) · P_L(cos θ) using the point's CM angle. Each A_L(E) can be FFT-convolved independently, then summed with P_L(cos θ) at the data point angle.

A single ESpectrum stores both representations simultaneously (they share the same matrix inversion at each grid point), so no duplication of physics computation.

### Target Effects

ESpectrum handles all experimental effects:
- **Gaussian beam convolution** — FFT approach (primary motivation).
- **Through-target integration (stopping power)** — numerical integration using the same adaptive grid, weighted by the stopping power equation from the segment's TargetEffect.
- **Combined convolution + target integration** — applied in sequence: stopping-power integration first, then Gaussian convolution.

For target integration the spectrum grid range is automatically extended by the maximum target depth (see Grid Range section below).

### SUM/RATIO Segments

Initial implementation: only base (non-advanced) segments use ESpectrum. Advanced SUM/RATIO segments fall back to the existing point-by-point quadrature path. This can be extended later by computing one ESpectrum per component pair and combining at query time.

### Grid Energy Range

The grid covers the union of energy ranges across **all** ESegments sharing the same (entrance, exit) pair, with the following padding:

- **No target effects:** extend ±5% beyond the min/max data energy (to avoid edge artefacts in interpolation).
- **Beam convolution active (σ > 0):** extend the range by ±3σ_max (widest beam width among all segments with this pair).
- **Target integration active:** extend the low-energy edge by the maximum target depth converted to energy loss.

Only one ESpectrum is ever created per (entrance, exit) pair regardless of how many ESegments share that pair.

### Scope

- **API mode:** not included in this implementation. API users call `EPoint::Calculate()` directly.
- **CLI mode and GUI mode:** fully supported.

---

## New Files

### `include/ESpectrum.h` / `src/ESpectrum.cpp`

```cpp
class ESpectrum {
public:
  ESpectrum(int entranceKey, int exitKey);

  // Phase 1: build adaptive energy grid (called once at initialization)
  void InitializeGrid(CNuc* cnuc, const Configure& cfg,
                      double eMin, double eMax,
                      double beamSigmaMax, double targetDepthMax);

  // Phase 2: recompute cross sections on grid (called every Minuit iteration)
  void Calculate(CNuc* cnuc, const Configure& cfg);

  // Query for total cross section at given energy
  double InterpolateTotal(double energy) const;

  // Query for differential cross section at given energy and CM angle
  double InterpolateDifferential(double energy, double cmAngle) const;

  // FFT convolution: returns convolved σ at each queryEnergy
  // isDifferential selects total vs angular mode; cmAngle used if differential
  std::vector<double> ConvolveGaussian(
      const std::vector<double>& queryEnergies,
      bool isDifferential,
      double cmAngle,
      double sigma) const;

  // Target integration at a single energy using stopping power
  double IntegrateTarget(double energy, bool isDifferential,
                         double cmAngle, const TargetEffect& effect) const;

  bool IsReady()       const { return isReady_; }
  int  GetEntranceKey() const { return entranceKey_; }
  int  GetExitKey()     const { return exitKey_; }

private:
  int  entranceKey_, exitKey_;
  bool isReady_ = false;

  // Adaptive grid
  std::vector<double>              gridEnergies_;   // energy at each node
  std::vector<double>              gridTotal_;      // σ_total(E)
  std::vector<std::vector<double>> gridLegendre_;   // A_L(E) for each L

  // Uniform grid (built once from adaptive grid, for FFT)
  std::vector<double> uniformEnergies_;
  std::vector<double> uniformTotal_;
  std::vector<std::vector<double>> uniformLegendre_;

  int maxLOrder_ = 0;

  void BuildUniformGrid();
  void InterpolateAdaptiveToUniform();
  std::vector<double> FFTConvolve(const std::vector<double>& signal,
                                  double sigma, double step) const;
};
```

### Internal Physics Approach (Option A)

ESpectrum internally creates a temporary `ESegment` populated with `EPoint` objects at each grid energy, then calls the standard initialization and calculation path on it. This reuses 100% of existing physics code (Coulomb functions, penetrabilities, Legendre polynomials, matrix inversion) with no duplication. It is identical to how `<segmentsTest>` extrapolation segments work, constructed programmatically rather than from the config file.

---

## Modified Files

### 1. `include/Config.h` — new flag

```cpp
USE_SPECTRUM = (1<<19)   // enable ESpectrum adaptive grid + FFT convolution
```

### 2. `src/Config.cpp` — parse option

Add parsing for `useSpectrum` inside the existing `<config>` block parser, analogous to `useExternalCapture`.

### 3. `include/EData.h` — spectrum map and API

```cpp
#include "ESpectrum.h"
#include <map>
#include <utility>

// New members inside EData:
std::map<std::pair<int,int>, ESpectrum> spectra_;

void InitializeSpectra(CNuc* cnuc, const Configure& cfg);
void CalculateSpectra(CNuc* cnuc, const Configure& cfg);
ESpectrum* GetSpectrum(int entranceKey, int exitKey);
```

### 4. `src/EData.cpp` — implement spectrum management

**`InitializeSpectra()`** — called from `EData::Initialize()` when `USE_SPECTRUM` is set:
1. Walk all `segments_` and collect unique `(entranceKey_, exitKey_)` pairs.
2. For each unique pair, determine the union energy range and the maximum beam sigma / target depth across all ESegments sharing that pair.
3. Construct one `ESpectrum` and call `InitializeGrid()`.
4. Store in `spectra_` map.
5. For each `ESegment`, call `segment.SetSpectrum(GetSpectrum(...))` to wire the pointer.

**`CalculateSpectra()`** — called every Minuit iteration when `USE_SPECTRUM` is set:
1. For each `ESpectrum` in `spectra_`, call `spectrum.Calculate(cnuc, cfg)`.

### 5. `include/ESegment.h` — spectrum pointer

```cpp
private:
  ESpectrum* spectrum_ = nullptr;   // non-owning; set by EData::InitializeSpectra

public:
  void SetSpectrum(ESpectrum* s) { spectrum_ = s; }
  ESpectrum* GetSpectrum() const  { return spectrum_; }
```

### 6. `src/ESegment.cpp` — use spectrum in `CalculateTheoreticalCrossSection()`

The existing function remains fully intact. Add a new early-return path when a spectrum is available:

```cpp
double ESegment::CalculateTheoreticalCrossSection(int pointIndex,
    CNuc* cnuc, const Configure& cfg)
{
  EPoint* point = GetPoint(pointIndex);

  // Spectrum-based path (only for non-advanced base segments)
  if (spectrum_ && spectrum_->IsReady() && !isAdvanced_) {
    return CalculateViaSpectrum(pointIndex, cnuc, cfg);
  }

  // ... existing logic unchanged ...
}

double ESegment::CalculateViaSpectrum(int pointIndex,
    CNuc* cnuc, const Configure& cfg)
{
  EPoint* point      = GetPoint(pointIndex);
  bool    isDiff     = IsDifferential();
  double  angle      = point->GetCMAngle();
  double  energy     = point->GetCMEnergy();
  double  conversion = point->GetCrossSectionConversionFactor(cnuc, cfg);

  if (point->HasTargetEffect()) {
    const TargetEffect& effect = *point->GetTargetEffect();
    double sigma = effect.IsConvolution() ? effect.GetSigma() : 0.0;

    if (effect.IsTargetIntegration()) {
      double integrated = spectrum_->IntegrateTarget(energy, isDiff, angle, effect);
      if (sigma > 0.0) {
        // convolve the already-integrated result
        std::vector<double> conv = spectrum_->ConvolveGaussian(
            {energy}, isDiff, angle, sigma);
        return conv[0] * conversion;
      }
      return integrated * conversion;
    }

    if (sigma > 0.0) {
      std::vector<double> conv = spectrum_->ConvolveGaussian(
          {energy}, isDiff, angle, sigma);
      return conv[0] * conversion;
    }
  }

  // No target effect: plain interpolation
  double xsec = isDiff
      ? spectrum_->InterpolateDifferential(energy, angle)
      : spectrum_->InterpolateTotal(energy);
  return xsec * conversion;
}
```

### 7. `src/AZURECalc.cpp` — trigger spectrum calculation

```cpp
// After existing parameter fills:
localCompound->FillCompoundFromParams(p);
localData->FillNormsFromParams(p);
localData->FillEnergyShiftsFromParams(p);

// NEW — recompute all spectra with updated parameters
if (configure.paramMask & Config::USE_SPECTRUM)
  localData->CalculateSpectra(localCompound, configure);

// Existing segment loop unchanged
```

Because ESpectrum is owned by EData, cloning EData (for the thread-safe object pool) automatically clones the spectra. Each thread's `localData` fills its own spectra independently — no locking required.

### 8. `CMakeLists.txt`

Add `src/ESpectrum.cpp` to the source list alongside the other `src/*.cpp` files.

---

## Adaptive Grid Configuration

Reuse the existing `AdaptiveIntegrationGrid` + `GridConfig`:

```cpp
GridConfig spectrumGridConfig;
spectrumGridConfig.maxPoints                = 2000;   // configurable
spectrumGridConfig.baseEnergyStep           = 0.005;  // MeV
spectrumGridConfig.resonanceWidthMultiplier = 10.0;   // narrower than target integration
spectrumGridConfig.pointsPerWidth           = 30.0;
spectrumGridConfig.entranceKey              = entranceKey_;
```

---

## FFT Convolution Strategy

For Gaussian beam convolution with width σ:

1. Adaptive grid → interpolate to a uniform grid in E with step `Δ = min(σ/10, baseStep)`.
2. Zero-pad to next power of 2.
3. FFT the cross section array.
4. Multiply by the Gaussian kernel in frequency space: `G(k) = exp(-2π²σ²k²)`.
5. Inverse FFT.
6. Sample at data-point energies with cubic spline interpolation.

For the FFT backend, a self-contained real Cooley-Tukey FFT (~100 lines) is implemented as a private utility inside `ESpectrum.cpp` to avoid adding an `fftw` dependency.

---

## Calculation Timing

A timing wrapper is added to report how long each full calculation run takes. The output goes to the GUI log via the existing `std::ostream` logging chain (`TextEditBuffer` → `FilteredTextEdit`).

### `gui/src/AZUREMainThread.cpp` — wrap `azureMain_()` call

```cpp
#include <chrono>

void AZUREMainThreadWorker::run() {
  // ... existing cache init ...

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

  emit done();
}
```

---

## GUI: Runtime Options Checkbox

### `gui/include/EditOptionsDialog.h`

Add alongside the existing checkboxes:

```cpp
QCheckBox* useSpectrumCheck;
```

### `gui/src/EditOptionsDialog.cpp`

Add to the constructor and layout:

```cpp
useSpectrumCheck = new QCheckBox(
    "Use adaptive spectrum grid for convolution (ESpectrum)");
useSpectrumCheck->setChecked(
    configure->paramMask & Config::USE_SPECTRUM);
// add to the options layout
```

### `gui/src/AZURESetup.cpp` — `editOptions()`

Add the flag mapping alongside the existing ones:

```cpp
if (aDialog.useSpectrumCheck->isChecked())
  configure_.paramMask |= Config::USE_SPECTRUM;
else
  configure_.paramMask &= ~Config::USE_SPECTRUM;
```

---

## File Change Summary

| File | Type | What changes |
|---|---|---|
| `include/ESpectrum.h` | **NEW** | Full class declaration |
| `src/ESpectrum.cpp` | **NEW** | Grid init, calculate, interpolate, FFT convolve, target integrate |
| `include/Config.h` | modify | Add `USE_SPECTRUM = (1<<19)` |
| `src/Config.cpp` | modify | Parse `useSpectrum` keyword |
| `include/EData.h` | modify | Add `spectra_` map + 3 new methods |
| `src/EData.cpp` | modify | Implement `InitializeSpectra()`, `CalculateSpectra()`, `GetSpectrum()` |
| `include/ESegment.h` | modify | Add `spectrum_*` pointer + getter/setter |
| `src/ESegment.cpp` | modify | Spectrum-path branch + `CalculateViaSpectrum()` |
| `src/AZURECalc.cpp` | modify | One `CalculateSpectra()` call before segment loop |
| `CMakeLists.txt` | modify | Add `src/ESpectrum.cpp` to sources |
| `gui/include/EditOptionsDialog.h` | modify | Add `useSpectrumCheck` member |
| `gui/src/EditOptionsDialog.cpp` | modify | Add checkbox to Runtime Options dialog |
| `gui/src/AZURESetup.cpp` | modify | Map checkbox to `Config::USE_SPECTRUM` in `editOptions()` |
| `gui/src/AZUREMainThread.cpp` | modify | Add `std::chrono` timing around `azureMain_()` |
