#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "ReactionRate.h"
#include "ECTable.h"
#include "AdaptiveIntegrationGrid.h"
#include <iomanip>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <time.h>

// Energy bounds (c.m., MeV) of the rate integral and of the support of sigma(E).
// The integrand sqrt(E) exp(-E/kT) of an s-wave neutron entrance (sigma ~ 1/v) carries
// about 1e-6 of the rate below 0.1 eV at T9 = 0.01; a charged entrance underflows to
// sigma = 0 long before.  The floor was 1 keV, which cost a neutron rate 49 % at
// T9 = 0.01 and 2.8 % at T9 = 0.1 (2026-10-25 audit).
static const double kRateEnergyLow = 1.0e-7;
static const double kRateEnergyHigh = 50.0;

struct gsl_reactionrate_params {
  gsl_reactionrate_params(const Config &config) :
    configure(config) {};
  const Config &configure;
  double temperature;
  CNuc *compound;
  int entranceKey;
  int exitKey;
};

// The quadratures run over u = ln E: sigma(E) E^2 exp(-E/kT) du.  Over the 8.7 decades
// of the energy range that integrand is smooth and bounded, where sigma E exp(-E/kT) dE
// with its 1/v or Gamow rise made gsl_integration_qagp stop 3e-4 short of its 1e-6
// request (n + 12C at T9 = 0.1, the same on the spline and the direct path).
double gsl_reactionrate_integrand(double u, void *p) {
  struct gsl_reactionrate_params *params = (struct gsl_reactionrate_params *)p;
  CNuc *compound = params->compound;
  const Config &configure = params->configure;
  double temperature = params->temperature;
  int entranceKey = params->entranceKey;
  int exitKey = params->exitKey;
  const double x = exp(u);

  double crossSection;
  if (x <= kRateEnergyHigh && x >= kRateEnergyLow) {
    EPoint *point = new EPoint(55.0, x, entranceKey, exitKey, false, false, false, 0.0, 0, 0);
    point->Initialize(compound, configure);
    point->Calculate(compound, configure);
    crossSection = point->GetFitCrossSection();
    delete point;
    if (!std::isfinite(crossSection) || crossSection < 0.0) crossSection = 0.0;
  } else
    crossSection = 0.0;

  return crossSection * x * x * exp(-x / temperature / boltzConst);
}

// Breakpoints in energy -> strictly increasing ln E for qagp.
static std::vector<double> LogBreakpoints(const std::vector<double> &pts) {
  std::vector<double> u;
  u.reserve(pts.size());
  for (size_t i = 0; i < pts.size(); ++i)
    if (pts[i] > 0.0 && (u.empty() || log(pts[i]) > u.back())) u.push_back(log(pts[i]));
  return u;
}

static const double kRateEpsRel = 1e-7;

/*!
 * The resonances of the parameters being evaluated, as observed energies and widths.
 * ReactionRate::ReactionRate fills the compound from the parameter vector (FitE /
 * FitGamma); the input values a level still carries (the .azr values after
 * TransformIn, formal ones under --no-brune) are not those parameters.  As
 * EData::ReplaceSubPointGrids does for the target-effect sub-grids, a private clone is
 * transformed out and refilled with the physical values, and the grid reads those.
 * Before, the clusters and breakpoints came from the input values, so a narrow level
 * that a parameter file had moved by more than its cluster (64 widths) was never
 * resolved and dropped out of the rate.
 */
static std::vector<AdaptiveIntegrationGrid::ResonanceInfo>
RateResonances(CNuc *compound, const Config &configure, int entranceKey) {
  AdaptiveIntegrationGrid::GridConfig gc;
  gc.entranceKey = entranceKey;
  try {
    std::unique_ptr<CNuc> g(compound->Clone());
    for (int j = 1; j <= compound->NumJGroups(); j++)
      for (int la = 1; la <= compound->GetJGroup(j)->NumLevels(); la++) {
        ALevel *from = compound->GetJGroup(j)->GetLevel(la);
        ALevel *to = g->GetJGroup(j)->GetLevel(la);
        to->SetFitE(from->GetFitE());
        for (int ch = 1; ch <= compound->GetJGroup(j)->NumChannels(); ch++)
          to->SetFitGamma(ch, from->GetFitGamma(ch));
      }
    g->TransformOut(configure);
    g->FillCompoundFromParamsPhysical(g->GetTransformParams(configure));
    gc.inputWidthsArePhysical = true;
    return AdaptiveIntegrationGrid(gc).GetResonances(kRateEnergyHigh, kRateEnergyLow, g.get());
  } catch (const std::exception &e) {
    configure.outStream << "\tWarning: could not transform the parameters to observed values for the "
                           "rate grid (" << e.what() << "); using the input level values." << std::endl;
    gc.inputWidthsArePhysical = false;
    return AdaptiveIntegrationGrid(gc).GetResonances(kRateEnergyHigh, kRateEnergyLow, compound);
  }
}

/*!
 * Builds a sorted list of integration breakpoints for gsl_integration_qagp.
 * The interior points are placed at each resonance energy (and a few points
 * across its width) so the adaptive integrator is guaranteed to resolve narrow
 * resonances instead of stepping over them.  Returns at least {E_LOW, E_HIGH}.
 */
std::vector<double> gsl_reactionrate_breakpoints(CNuc *compound, const Config &configure, int entranceKey) {
  std::vector<double> pts;
  pts.push_back(kRateEnergyLow);
  pts.push_back(kRateEnergyHigh);

  std::vector<AdaptiveIntegrationGrid::ResonanceInfo> resonances =
      RateResonances(compound, configure, entranceKey);

  // For each resonance, bracket the peak tightly so a Gauss-Kronrod panel is
  // forced right onto it, regardless of interference-driven shape changes.
  // Out to 1024 widths: a Lorentzian still holds 0.06 % of its area beyond that, and the
  // adaptive rule needs the scale of the tails spelled out.
  const double widthMults[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0};
  for (size_t i = 0; i < resonances.size(); ++i) {
    double e0 = resonances[i].energy;
    double w = resonances[i].totalWidth;
    for (size_t m = 0; m < sizeof(widthMults) / sizeof(widthMults[0]); ++m) {
      double lo = e0 - widthMults[m] * w;
      double hi = e0 + widthMults[m] * w;
      if (lo > kRateEnergyLow && lo < kRateEnergyHigh) pts.push_back(lo);
      if (hi > kRateEnergyLow && hi < kRateEnergyHigh) pts.push_back(hi);
    }
  }

  std::sort(pts.begin(), pts.end());
  // gsl_integration_qagp requires strictly increasing breakpoints.
  std::vector<double> unique;
  unique.push_back(pts[0]);
  for (size_t i = 1; i < pts.size(); ++i)
    if (pts[i] > unique.back() * (1.0 + 1e-12)) unique.push_back(pts[i]);
  return unique;
}

double gsl_reactionrate_integration(double temperature, CNuc *compound, const Config &configure,
                                    int entranceKey, int exitKey,
                                    const std::vector<double> &breakpoints) {
  struct gsl_reactionrate_params params(configure);
  params.temperature = temperature;
  params.compound = compound;
  params.entranceKey = entranceKey;
  params.exitKey = exitKey;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);

  gsl_function F;
  F.function = &gsl_reactionrate_integrand;
  F.params = &params;

  double result, error;

  // qagp subdivides at the supplied resonance breakpoints, so narrow resonances
  // near threshold are always sampled (qagiu could step right over them).
  std::vector<double> pts = breakpoints;
  if (pts.size() < 2) {
    pts.clear();
    pts.push_back(kRateEnergyLow);
    pts.push_back(kRateEnergyHigh);
  }
  // A benign GSL_EROUND ("roundoff") status is expected for this steep,
  // wide-dynamic-range integrand; the returned value is still accurate, so we
  // disable the throwing error handler around the call rather than abort.
  std::vector<double> u = LogBreakpoints(pts);
  gsl_error_handler_t *oldHandler = gsl_set_error_handler_off();
  gsl_integration_qagp(&F, u.data(), u.size(), 0.0, kRateEpsRel, 2000, w, &result, &error);
  gsl_set_error_handler(oldHandler);

  double rate = 1e-24 * avagadroNum * lightSpeedInCmPerS *
      sqrt(8.0 / pi / compound->GetPair(compound->GetPairNumFromKey(entranceKey))->GetRedMass() / uconv) /
      pow(boltzConst * temperature, 1.5) * result;

  gsl_integration_workspace_free(w);

  return rate;
}

// Tabulated path: compute sigma(E) once on a resonance-aware grid, fit a Steffen
// spline to ln sigma against ln E, and integrate every temperature against it.
// Selected by Config::useAdaptiveGrid; useAdaptiveGrid=false uses the per-temperature
// path.  The spline was of sigma itself against E: below the barrier sigma changes by a
// factor 2-3 between the 48-per-decade backbone nodes, and a monotone spline of convex
// data lies below it, so every non-resonant charged-particle rate came out low (1.15 %
// for 12C(a,g) at T9 = 0.1, 0.77 % for 12C(p,g) at T9 = 0.01).  In ln-ln a 1/v cross
// section is a straight line and a Gamow tail is nearly one; the error is ~1e-5.

namespace {

// 25-cell progress bar; '\r'-prefixed so the GUI updates the line in place.
void PrintRateProgress(std::ostream &os, int done, int total) {
  double frac = (total > 0) ? (double)done / (double)total : 1.0;
  if (frac > 1.0) frac = 1.0;
  int filled = (int)(frac * 25.0);
  std::string bar = "[";
  for (int j = 0; j < 25; ++j) bar += (j < filled) ? '*' : ' ';
  bar += "] ";
  os << "\r\t" << bar << (int)(frac * 100.0 + 0.5) << '%';
  os.flush();
}

// Energy grid for sigma(E): a log backbone plus a resonance cluster per level.
// Returns strictly increasing, unique energies (as gsl_spline_init needs).
std::vector<double> BuildRateCrossSectionGrid(CNuc *compound, int entranceKey,
                                              const Config &configure) {
  const double eLow = kRateEnergyLow;
  const double eHigh = kRateEnergyHigh;
  std::vector<double> pts;

  // Log backbone.
  const int pointsPerDecade = 48;
  double lgLo = std::log10(eLow), lgHi = std::log10(eHigh);
  int nLog = std::max(2, (int)std::ceil((lgHi - lgLo) * pointsPerDecade));
  for (int i = 0; i <= nLog; ++i)
    pts.push_back(std::pow(10.0, lgLo + (lgHi - lgLo) * (double)i / (double)nLog));

  std::vector<AdaptiveIntegrationGrid::ResonanceInfo> res =
      RateResonances(compound, configure, entranceKey);

  // sigma(E) has width Gamma_particle (non-RMC) or Gamma_total (RMC).
  bool rmc = (configure.paramMask & Config::USE_RMC_FORMALISM);

  // Symmetric cluster per resonance: the peak e0 itself, plus geometric offsets
  // on both sides from the core (~w/32) to the wings (~64w).  The ln-sigma spline of
  // a Lorentzian is limited by its curvature at the peak, so the core is dense.
  const int nOff = 40;
  const double coreFrac = 1.0 / 32.0;
  const double wingFactor = 64.0;
  for (size_t r = 0; r < res.size(); ++r) {
    double e0 = res[r].energy;
    double w = rmc ? res[r].totalWidth : res[r].particleWidth;
    if (!(w > 0.0)) w = res[r].totalWidth;
    if (!(w > 0.0)) continue;
    if (e0 > eLow && e0 < eHigh) pts.push_back(e0);
    double dmin = coreFrac * w, dmax = wingFactor * w;
    for (int k = 0; k <= nOff; ++k) {
      double d = dmin * std::pow(dmax / dmin, (double)k / (double)nOff);
      if (e0 - d > eLow && e0 - d < eHigh) pts.push_back(e0 - d);
      if (e0 + d > eLow && e0 + d < eHigh) pts.push_back(e0 + d);
    }
  }

  pts.push_back(eLow);
  pts.push_back(eHigh);
  std::sort(pts.begin(), pts.end());

  std::vector<double> uniq;
  uniq.reserve(pts.size());
  uniq.push_back(pts[0]);
  for (size_t i = 1; i < pts.size(); ++i)
    if (pts[i] > uniq.back() * (1.0 + 1e-9)) uniq.push_back(pts[i]);
  return uniq;
}

// Fit cross section at one CM energy (same EPoint path as the legacy integrand).
// Returns NaN if the level matrix is singular here (e.g. exactly on a pole with
// vanishing penetrability) or the result is not finite; the caller drops such points.
double RateCrossSectionAt(double energy, CNuc *compound, const Config &configure,
                          int entranceKey, int exitKey) {
  if (!(energy >= kRateEnergyLow && energy <= kRateEnergyHigh)) return 0.0;
  try {
    EPoint point(55.0, energy, entranceKey, exitKey, false, false, false, 0.0, 0, 0);
    point.Initialize(compound, configure);
    point.Calculate(compound, configure);
    double sigma = point.GetFitCrossSection();
    return std::isfinite(sigma) ? sigma : std::numeric_limits<double>::quiet_NaN();
  } catch (const std::exception &) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

// The spline holds ln sigma against ln E over the grid points where sigma > 0; sigma is
// zero outside [eMin, eMax] (an underflowed Coulomb tail, or below a threshold).
struct SplineRateParams {
  gsl_spline *spline;
  gsl_interp_accel *accel;
  double eMin;
  double eMax;
  double temperature;
};

// sigma(E) E^2 exp(-E/kT) at u = ln E, sigma from the spline (zero outside [eMin,eMax]).
double spline_rate_integrand(double u, void *p) {
  SplineRateParams *q = (SplineRateParams *)p;
  const double x = exp(u);
  if (x < q->eMin || x > q->eMax) return 0.0;
  double sigma = exp(gsl_spline_eval(q->spline, u, q->accel));
  return sigma * x * x * exp(-x / q->temperature / boltzConst);
}

// Integrate the spline rate integrand at one temperature (own accel/workspace,
// so it is safe to run concurrently over temperatures on a shared spline).  The
// knots of the spline are the breakpoints: between two knots the integrand is a cubic
// in ln sigma times E^2 exp(-E/kT), which the Gauss-Kronrod rule integrates to rounding,
// so the quadrature adds nothing to the spline's own error.  (Breakpoints at a
// resonance's +-8 widths only, as before, left qagp 5e-4 short on a 22 eV level.)
double IntegrateSplineRate(gsl_spline *spline, double eMin, double eMax,
                           double temperature, const std::vector<double> &lnKnots) {
  gsl_interp_accel *accel = gsl_interp_accel_alloc();
  SplineRateParams params{spline, accel, eMin, eMax, temperature};

  gsl_function F;
  F.function = &spline_rate_integrand;
  F.params = &params;

  std::vector<double> u = lnKnots;   // strictly increasing, from eMin to eMax
  const size_t limit = 4 * u.size() + 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);
  double result = 0.0, error = 0.0;
  gsl_error_handler_t *oldHandler = gsl_set_error_handler_off();
  gsl_integration_qagp(&F, u.data(), u.size(), 0.0, kRateEpsRel, limit, w, &result, &error);
  gsl_set_error_handler(oldHandler);
  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(accel);
  return result;
}

// Per-temperature path (Config::useAdaptiveGrid=false, or small-grid fallback).
void RunAdaptiveRateCalculation(const std::vector<double> &temps,
                                CNuc *compound, const Config &configure,
                                int entranceKey, int exitKey,
                                std::vector<RateData> &outRates) {
  std::ostream &os = configure.outStream;
  int nT = (int)temps.size();
  outRates.assign(nT, RateData(0.0, 0.0));

  os << "\t[                         ] 0%";
  os.flush();
  std::vector<double> breakpoints = gsl_reactionrate_breakpoints(compound, configure, entranceKey);

  std::atomic<int> done(0);
  time_t startTime = time(NULL);
#pragma omp parallel
  {
    CNuc *localCompound = compound->Clone();
    const Config localConfigure = configure;
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nT; ++i) {
      double rate = gsl_reactionrate_integration(temps[i], localCompound, localConfigure,
                                                 entranceKey, exitKey, breakpoints);
      outRates[i] = RateData(temps[i], rate);
      int d = ++done;
#pragma omp critical(rate_progress)
      {
        if (difftime(time(NULL), startTime) > 0.25 || d == nT) {
          startTime = time(NULL);
          PrintRateProgress(os, d, nT);
        }
      }
    }
    delete localCompound;
  }
  os << "\r\t[*************************] 100%" << std::endl;
}

// Tabulate sigma(E) once, spline it, integrate every temperature against it.
void RunTabulatedRateCalculation(const std::vector<double> &temps,
                                 CNuc *compound, const Config &configure,
                                 int entranceKey, int exitKey,
                                 std::vector<RateData> &outRates) {
  std::ostream &os = configure.outStream;
  int nT = (int)temps.size();
  outRates.assign(nT, RateData(0.0, 0.0));

  // ---- Phase 1: temperature-independent cross section on the energy grid ----
  std::vector<double> energies = BuildRateCrossSectionGrid(compound, entranceKey, configure);
  int nGrid = (int)energies.size();
  std::vector<double> sigma(nGrid, 0.0);

  os << "\tCalculating cross section on adaptive energy grid ("
     << nGrid << " points)..." << std::endl;
  os << "\t[                         ] 0%";
  os.flush();

  std::atomic<int> doneGrid(0);
  time_t startTime = time(NULL);
#pragma omp parallel
  {
    CNuc *localCompound = compound->Clone();
    const Config localConfigure = configure;
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nGrid; ++i) {
      sigma[i] = RateCrossSectionAt(energies[i], localCompound, localConfigure,
                                    entranceKey, exitKey);
      int d = ++doneGrid;
#pragma omp critical(rate_progress)
      {
        if (difftime(time(NULL), startTime) > 0.25 || d == nGrid) {
          startTime = time(NULL);
          PrintRateProgress(os, d, nGrid);
        }
      }
    }
    delete localCompound;
  }
  os << "\r\t[*************************] 100%" << std::endl;

  // Drop points with a singular level matrix (NaN); clamp negatives to zero.
  {
    std::vector<double> gE, gS;
    gE.reserve(nGrid);
    gS.reserve(nGrid);
    int nDropped = 0;
    std::streamsize oldPrec = os.precision();
    os.precision(10);
    for (int i = 0; i < nGrid; ++i) {
      if (std::isnan(sigma[i])) {
        os << "\tWARNING: singular level matrix at E_cm = " << energies[i]
           << " MeV; dropping this grid point." << std::endl;
        ++nDropped;
        continue;
      }
      gE.push_back(energies[i]);
      gS.push_back(sigma[i] < 0.0 ? 0.0 : sigma[i]);
    }
    os.precision(oldPrec);
    if (nDropped)
      os << "\tDropped " << nDropped << " singular grid point(s) total." << std::endl;
    energies.swap(gE);
    sigma.swap(gS);
    nGrid = (int)energies.size();
  }

  auto *entrancePair = compound->GetPair(compound->GetPairNumFromKey(entranceKey));
  double redMass = entrancePair->GetRedMass();

  // Dump the grid (E_cm, sigma, S-factor) for cross-checking; S = sigma*E*exp(2*pi*eta).
  {
    double z1 = entrancePair->GetZ(1);
    double z2 = entrancePair->GetZ(2);
    std::string gridFile = configure.outputdir + "cross_section_grid.dat";
    std::ofstream gout(gridFile.c_str());
    if (gout) {
      gout << std::setw(26) << "E_cm[MeV]"
           << std::setw(26) << "CrossSection"
           << std::setw(26) << "S-factor" << std::endl;
      gout << std::scientific << std::setprecision(16);
      for (int i = 0; i < nGrid; ++i) {
        double e = energies[i];
        double sfactor = 0.0;
        if (e > 0.0 && sigma[i] > 0.0) {
          // In logs: exp(2 pi eta) overflows below a few keV.  Where S still does not fit a
          // double the Coulomb tail has long underflowed and the number means nothing: 0.
          double eta = sqrt(uconv / 2.0) * fstruc * z1 * z2 * sqrt(redMass / e);
          double lnS = log(sigma[i]) + log(e) + 2.0 * pi * eta;
          sfactor = (lnS < 700.0) ? exp(lnS) : 0.0;
        }
        gout << std::setw(26) << e
             << std::setw(26) << sigma[i]
             << std::setw(26) << sfactor << std::endl;
      }
      gout.close();
      os << "\tWrote cross section grid (" << nGrid << " points) to "
         << gridFile << std::endl;
    } else {
      os << "\tWarning: could not write " << gridFile << std::endl;
    }
    if (configure.useECTable && ECTable::Lookups() + ECTable::DirectCalls() > 0)
      os << "\t" << ECTable::Report() << std::endl;
  }

  // ln sigma against ln E over the points with sigma > 0.  The rest (an underflowed
  // Coulomb tail at the low end, a closed exit below its threshold) is zero in the
  // integrand, and the spline does not cross them.
  std::vector<double> lnE, lnSigma;
  lnE.reserve(nGrid);
  lnSigma.reserve(nGrid);
  for (int i = 0; i < nGrid; ++i)
    if (sigma[i] > 0.0) {
      lnE.push_back(std::log(energies[i]));
      lnSigma.push_back(std::log(sigma[i]));
    }
  int nSpline = (int)lnE.size();

  // Steffen needs at least a handful of points; fall back if the grid collapsed.
  if (nSpline < 3) {
    os << "\tGrid too small for spline; using per-temperature integration."
       << std::endl;
    RunAdaptiveRateCalculation(temps, compound, configure, entranceKey, exitKey, outRates);
    return;
  }

  gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, nSpline);
  gsl_spline_init(spline, lnE.data(), lnSigma.data(), nSpline);
  double eMin = std::exp(lnE.front());
  double eMax = std::exp(lnE.back());

  // ---- Phase 2: integrate each temperature against the fixed spline ----
  os << "\tIntegrating reaction rate over " << nT << " temperatures..."
     << std::endl;
  os << "\t[                         ] 0%";
  os.flush();

  std::atomic<int> doneT(0);
  time_t startT2 = time(NULL);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nT; ++i) {
    double T = temps[i];
    double integral = IntegrateSplineRate(spline, eMin, eMax, T, lnE);
    double rate = 1e-24 * avagadroNum * lightSpeedInCmPerS *
        sqrt(8.0 / pi / redMass / uconv) / pow(boltzConst * T, 1.5) * integral;
    outRates[i] = RateData(T, rate);
    int d = ++doneT;
#pragma omp critical(rate_progress)
    {
      if (difftime(time(NULL), startT2) > 0.25 || d == nT) {
        startT2 = time(NULL);
        PrintRateProgress(os, d, nT);
      }
    }
  }
  os << "\r\t[*************************] 100%" << std::endl;

  gsl_spline_free(spline);
}

/*! Dispatches to the tabulated or legacy path based on Config::useAdaptiveGrid. */
void DispatchRateCalculation(const std::vector<double> &temps,
                             CNuc *compound, const Config &configure,
                             int entranceKey, int exitKey,
                             std::vector<RateData> &outRates) {
  if (configure.useAdaptiveGrid)
    RunTabulatedRateCalculation(temps, compound, configure, entranceKey, exitKey, outRates);
  else
    RunAdaptiveRateCalculation(temps, compound, configure, entranceKey, exitKey, outRates);
}

}  // namespace

ReactionRate::ReactionRate(CNuc *compound, const vector_r &params,
                           const Config &configure, int entranceKey, int exitKey) :
  configure_(configure) {
  compound_ = compound;
  compound_->FillCompoundFromParams(params);
  entrance_key_ = entranceKey;
  exit_key_ = exitKey;
}

/*!
 * Calculates the astrophysical reaction rates over a range of stellar temperatures.
 */

void ReactionRate::CalculateRates() {
  int numSteps = (configure().rateParams.tempStep != 0.) ? int((configure().rateParams.maxTemp - configure().rateParams.minTemp) / configure().rateParams.tempStep) + 1 : 1;
  std::vector<double> temps(numSteps);
  for (int i = 0; i < numSteps; ++i)
    temps[i] = configure().rateParams.minTemp + i * configure().rateParams.tempStep;

  DispatchRateCalculation(temps, compound(), configure(), entranceKey(), exitKey(), rates_);
  std::sort(rates_.begin(), rates_.end());
}

/*!
 * Calculates the astrophysical reaction rates at temperatures from file.
 */

void ReactionRate::CalculateFileRates() {
  std::ifstream inFile(configure().rateParams.temperatureFile.c_str());
  if (!inFile) {
    configure().outStream << "Couldn't open temperature file." << std::endl;
    return;
  }

  std::vector<double> temps;
  while (!inFile.eof()) {
    std::string line;
    getline(inFile, line);
    if (!inFile.eof()) {
      double temp = 0.;
      std::istringstream stm;
      stm.str(line);
      if (stm >> temp) temps.push_back(temp);
    }
  }
  inFile.close();

  // Preserve the file ordering of the temperatures in the output.
  DispatchRateCalculation(temps, compound(), configure(), entranceKey(), exitKey(), rates_);
}

void ReactionRate::WriteRates() {
  std::string outputfile = configure().outputdir + "reactionrates.out";
  std::ofstream out;
  out.open(outputfile.c_str());
  if (out) {
    out << std::setw(20) << "T9" << std::setw(20) << "Rate" << std::endl;
    for (int i = 0; i < rates_.size(); i++) {
      out << std::setw(20) << rates_[i].temperature << std::setw(20) << rates_[i].rate << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure().outStream << "Could not write reaction rate file." << std::endl;
}
