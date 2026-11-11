#include "InputFiles.h"
#include "TargetEffect.h"
#include "AdaptiveIntegrationGrid.h"
#include <sstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <cmath>

/*!
 * Builds the effect from its spec (experiment/'s SetupSpecs.h).  A spec that is not active gives an inactive effect
 * with the defaults; the blocks' consistency (density fit, Doppler fit) is checked here and reported through
 * HasDensityFitError / HasDopplerError.
 */

TargetEffect::TargetEffect(const azure2::experiment::TargetEffectSpec &line, const Config &configure)
    : targetName_(line.targetName) {
  // Straggling off; a coefficient of 0 means Bohr straggling when it is on (0.04 keV^1/2 before 2026-10-27)
  isStraggling_ = false;
  stragglingCoefficient_ = 0.0;

  // Initialize adaptive grid parameters to defaults
  // Same defaults as AdaptiveIntegrationGrid::GridConfig (see there).
  {
    AdaptiveIntegrationGrid::GridConfig defaults;
    resonanceWidthMultiplier_ = defaults.resonanceWidthMultiplier;
    pointsPerWidth_ = defaults.pointsPerWidth;
  }

  // By default the effect applies to every point of its segments, with hard
  // edges and no automatic per-point decision.
  transitionWidth_ = 0.0;
  autoTolerance_ = 0.0;

  // An inactive spec (the old reader's line that ended after the parameter count) keeps every default.
  {
    const int numIntegrationPoints = line.numIntegrationPoints;
    const double sigma = line.sigma;
    const double density = line.density;
    const int numParameters = (int)line.stoppingPowerParameters.size();
    const vector_r &parameters = line.stoppingPowerParameters;
    const int numConvCoefficients = (int)line.convolutionCoefficientValues.size();
    const vector_r &convCoefficients = line.convolutionCoefficientValues;

    isQCoefficients_ = line.qCoefficients;
    qCoefficients_ = line.qCoefficientValues;
    isConvCoefficients_ = line.convolutionCoefficients;
    convCoefficients_ = convCoefficients;

    if (line.hasStragglingFlag) {
      isStraggling_ = line.straggling;
      if (line.hasStragglingCoefficient) stragglingCoefficient_ = line.stragglingCoefficient;
    }
    // The quadrature of the sub-point grid, with or without a target (an old line gives it only after its straggling
    // flag; a format-3 `quadrature` field of a beam-only data set has no straggling)
    if (line.hasResonanceWidthMultiplier) {
      resonanceWidthMultiplier_ = line.resonanceWidthMultiplier;
      explicitResonanceWidthMultiplier_ = true;
    }
    if (line.hasPointsPerWidth) pointsPerWidth_ = line.pointsPerWidth;

    if (line.hasRanges) {
      if (line.autoWindows) autoWindows_ = true;
      if (line.hasAutoWindowTolerance) autoWindowTolerance_ = line.autoWindowTolerance;
      if (line.hasAutoWindowMargin) autoWindowMargin_ = line.autoWindowMargin;
      ranges_ = line.ranges;
      if (line.transitionWidth > 0.) transitionWidth_ = line.transitionWidth;
      if (line.autoTolerance > 0.) autoTolerance_ = line.autoTolerance;
    }

    for (const azure2::experiment::TargetEffectSpec::ProfileComponent &p : line.beamProfileComponents) {
      BeamProfileComponent c;
      c.xi = p.xi;
      c.omega = p.omega;
      c.alpha = p.alpha;
      c.weight = p.weight;
      beamProfile_.push_back(c);
    }
    beamTpcSigma_ = line.beamTpcSigma;
    beamTruncation_ = line.beamTruncation;
    beamPhotodissociation_ = line.beamDetailedBalance;
    isBeamProfile_ = line.beamProfile;
    sigmaLab_ = line.sigmaLab;
    screeningNominal_ = line.screeningUe;
    screeningSigma_ = line.screeningSigma;
    screeningError_ = line.screeningMalformed;
    if (line.screening) {
      hasScreening_ = true;
      screeningVary_ = line.screeningVary;
      screeningUe_ = screeningNominal_;
    }
    const std::string &densitySigmaToken = line.densitySigmaText;
    const bool densityFitRead = line.densityFit;
    const double densityMaxRead = line.densityMax;
    densityVary_ = line.densityVary;
    densityFitError_ = !line.densityReadError.empty();
    densityFitErrorText_ = line.densityReadError;
    densityGroup_ = line.densityGroup;
    dopplerTemperature_ = line.dopplerTemperature;
    dopplerError_ = !line.dopplerReadError.empty();
    dopplerErrorText_ = line.dopplerReadError;
    isDoppler_ = line.doppler;
    const std::string &dopplerVaryToken = line.dopplerVaryText;
    const std::string &dopplerSigmaToken = line.dopplerSigmaText;
    const std::string &dopplerMaxToken = line.dopplerMaxText;
    const bool dopplerFitRead = line.dopplerFit;
    udrFile_ = line.udrFile;
    udrFlightPath_ = line.udrFlightPath;
    udrBurstFwhm_ = line.udrBurstFwhm;
    udrChannelWidth_ = line.udrChannelWidth;
    udrCentred_ = line.udrCentred;
    isUdr_ = line.udr;

    const std::string &stoppingPowerEq = line.stoppingPowerEquation;
    const std::string &convolutionEq = line.convolutionEquation;

    isActive_ = line.active;
    datasetKeys_ = line.datasetKeys;
    numIntegrationPoints_ = numIntegrationPoints;
    isConvolution_ = line.convolution;
    sigma_ = sigma;
    isTargetIntegration_ = line.targetIntegration;
    density_ = density;
    // density 0 (or negative) with target integration = an infinitely thick target (IsThickTarget)
    isThickTarget_ = isTargetIntegration_ && !(density_ > 0.0);
    isLayerYield_ = line.layerYield && isTargetIntegration_ && !isThickTarget_;
    if (densityFitRead || (!densityGroup_.empty() && !densityFitError_)) {
      if (!isTargetIntegration_ || isThickTarget_) {
        densityFitError_ = true;
        densityFitErrorText_ = "densityfit / densitygroup need target integration with a density > 0 (a thin target)";
      } else if (isBeamProfile_ || isUdr_) {
        densityFitError_ = true;
        densityFitErrorText_ = "densityfit / densitygroup cannot be combined with a beam profile or a UDR kernel";
      }
    }
    if (densityFitRead && !densityFitError_) {
      const bool percent = densitySigmaToken.back() == '%';
      const std::string number = percent ? densitySigmaToken.substr(0, densitySigmaToken.size() - 1) : densitySigmaToken;
      char *end = nullptr;
      const double s = std::strtod(number.c_str(), &end);
      if (number.empty() || end == number.c_str() || *end != '\0' || !(s >= 0.0)) {
        densityFitError_ = true;
        densityFitErrorText_ = "the densityfit sigma \"" + densitySigmaToken + "\" is not a non-negative number or percentage";
      } else {
        hasDensityFit_ = true;
        densityNominal_ = density_;
        densitySigma_ = percent ? s / 100.0 * density_ : s;
        // the grids span the largest density the fit is expected to reach (see HasDensityFit)
        layoutDensity_ = (densitySigma_ > 0.0) ? std::max(density_ + 3.0 * densitySigma_, 1.1 * density_) : 1.5 * density_;
        if (densityMaxRead > 0.0) {
          if (!(densityMaxRead > density_)) {
            densityFitError_ = true;
            densityFitErrorText_ = "the densityfit maximum must be above the line's density";
          }
          layoutDensity_ = densityMaxRead;
        }
      }
    }
    if (dopplerFitRead && !dopplerError_) {
      auto number = [](const std::string &tok, double &v) {
        char *end = nullptr;
        v = std::strtod(tok.c_str(), &end);
        return !tok.empty() && end != tok.c_str() && *end == '\0' && std::isfinite(v);
      };
      const std::string form = "a doppler block is 'doppler T' or 'doppler T vary sigma [Tmax]': ";
      double vary = -1.0, sigma = -1.0, tmax = 0.0;
      if (!number(dopplerVaryToken, vary) || !(vary == 0.0 || vary == 1.0)) {
        dopplerError_ = true;
        dopplerErrorText_ = form + "vary \"" + dopplerVaryToken + "\" is not 0 or 1";
      } else if (!number(dopplerSigmaToken, sigma) || !(sigma >= 0.0)) {
        dopplerError_ = true;
        dopplerErrorText_ = form + "sigma \"" + dopplerSigmaToken + "\" is not a non-negative number (K)";
      } else if (!dopplerMaxToken.empty() && (!number(dopplerMaxToken, tmax) || !(tmax > dopplerTemperature_))) {
        dopplerError_ = true;
        dopplerErrorText_ = form + "Tmax \"" + dopplerMaxToken + "\" must be a temperature above T (K)";
      } else {
        dopplerFit_ = true;
        dopplerVary_ = (vary == 1.0);
        dopplerNominal_ = dopplerTemperature_;
        dopplerSigma_ = sigma;
        // the grids span the largest temperature the fit is expected to reach (see HasDopplerFit)
        dopplerLayout_ = (sigma > 0.0) ? std::max(dopplerNominal_ + 3.0 * sigma, 1.21 * dopplerNominal_)
                                       : 2.0 * dopplerNominal_;
        if (tmax > 0.0) dopplerLayout_ = tmax;
        dopplerTmax_ = tmax > 0.0 ? tmax : 0.0;
      }
    }
    if (dopplerError_) isDoppler_ = false;
    if (isTargetIntegration_) {
      stoppingPowerEq_.Initialize(stoppingPowerEq, numParameters, configure);
      for (int i = 0; i < numParameters; i++) {
        stoppingPowerEq_.SetParameter(i, parameters[i], configure);
      }
    }
    if (isConvCoefficients_) {
      convolutionEq_.Initialize(convolutionEq, numConvCoefficients, configure);
      for (int i = 0; i < numConvCoefficients; i++) {
        convolutionEq_.SetParameter(i, convCoefficients[i], configure);
      }
    }
  }
}

/*!
 * Returns true if the target effect was marked as active in the target effects
 * input file, otherwise returns false.
 */

bool TargetEffect::IsActive() const {
  return isActive_;
}

/*!
 * Returns true if the target effect contains Gaussian beam convolution,
 * otherwise returns false.
 */

bool TargetEffect::IsConvolution() const {
  return isConvolution_;
}

/*!
 * Returns true if the target effect contains target integration, otherwise
 * returns false.
 */

bool TargetEffect::IsTargetIntegration() const {
  return isTargetIntegration_;
}

/*!
 * Returns true if the target effect is a beam-profile kernel.
 */

bool TargetEffect::IsBeamProfile() const {
  return isBeamProfile_;
}

/*!
 * Returns true if the effect integrates the observable over sub-points
 * (Gaussian convolution, target integration, energy-dependent convolution or
 * a beam-profile kernel).
 */

bool TargetEffect::IsSubPointEffect() const {
  return isConvolution_ || isTargetIntegration_ || isConvCoefficients_ || isBeamProfile_ || isUdr_ || isDoppler_;
}


double TargetEffect::DopplerUAt(double cmFactor, double temperature) const {
  if (!isDoppler_ || !(cmFactor > 0.0) || !(cmFactor < 1.0)) return 0.0;
  const double kT = 8.617333262e-11 * temperature;     // MeV
  return std::sqrt(kT * (1.0 - cmFactor) / cmFactor);  // m/M = (1 - c)/c for c = M/(m+M)
}

void TargetEffect::DopplerWindow(double eCM, double cmFactor, double &lo, double &hi) const {
  lo = hi = eCM;
  const double u = this->DopplerUAt(cmFactor, this->GetDopplerLayoutTemperature());
  if (!(u > 0.0) || !(eCM > 0.0)) return;
  // the kernel is a Gaussian in velocity, exp(-(v-w)^2/u^2): standard deviation u/sqrt(2)
  const double v = std::sqrt(eCM / cmFactor), dv = convolutionRange * u / std::sqrt(2.0);
  hi = cmFactor * (v + dv) * (v + dv);
  lo = (v > dv) ? cmFactor * (v - dv) * (v - dv) : 0.0;
}

bool TargetEffect::IsUdr() const { return isUdr_; }
const std::string &TargetEffect::GetUdrFile() const { return udrFile_; }
double TargetEffect::GetUdrFlightPath() const { return udrFlightPath_; }
double TargetEffect::GetUdrBurstFwhm() const { return udrBurstFwhm_; }
double TargetEffect::GetUdrChannelWidth() const { return udrChannelWidth_; }
bool TargetEffect::IsUdrCentred() const { return udrCentred_; }

/*!
 * Reads the SAMMY UDR file: free text up to a line containing "-----", then
 * blocks of an energy line (eV; anything after the number is ignored)
 * followed by (time microseconds, density) pairs, blocks separated by blank
 * lines.  Only the blocks whose energy lies within a decade of the requested
 * range, plus the two that bracket it, are kept: the n_TOF files tabulate
 * 600 energies over twelve decades with 9000 points each.
 */

bool TargetEffect::LoadUdrTable(double eMinLab, double eMaxLab, std::ostream &log, bool cached) {
  if (!isUdr_) return false;
  const double keepMin = eMinLab / 10.0, keepMax = eMaxLab * 10.0;
  if (udrTable_ && udrTable_->file == udrFile_ && !udrTable_->energy.empty() &&
      udrTable_->loadedMin <= keepMin && udrTable_->loadedMax >= keepMax)
    return true;
  std::unique_ptr<std::istream> inStream = InputFiles::Open(udrFile_, cached);
  std::istream &in = *inStream;
  if (!in) {
    log << "ERROR: cannot open the user-defined resolution file " << udrFile_ << std::endl;
    return false;
  }
  std::shared_ptr<UdrTable> table = std::make_shared<UdrTable>();
  table->file = udrFile_;
  table->loadedMin = udrTable_ ? std::min(udrTable_->loadedMin, keepMin) : keepMin;
  table->loadedMax = udrTable_ ? std::max(udrTable_->loadedMax, keepMax) : keepMax;
  std::string line;
  bool pastHeader = false;
  // Keep every block inside the range, and the last one below and the first
  // one above it, so that interpolation at the ends has both brackets.
  std::vector<double> pendingTau, pendingR;
  double pendingEnergy = -1.0, belowEnergy = -1.0;
  std::vector<double> belowTau, belowR;
  bool haveAbove = false;
  auto flush = [&]() {
    if (pendingEnergy < 0.0 || pendingTau.size() < 2) { pendingTau.clear(); pendingR.clear(); pendingEnergy = -1.0; return; }
    if (pendingEnergy < table->loadedMin) {
      belowEnergy = pendingEnergy; belowTau.swap(pendingTau); belowR.swap(pendingR);
    } else if (pendingEnergy <= table->loadedMax || !haveAbove) {
      if (belowEnergy >= 0.0 && table->energy.empty()) {
        table->energy.push_back(belowEnergy); table->tau.push_back(belowTau); table->r.push_back(belowR);
        belowEnergy = -1.0;
      }
      if (pendingEnergy > table->loadedMax) haveAbove = true;
      table->energy.push_back(pendingEnergy);
      table->tau.push_back(std::vector<double>()); table->tau.back().swap(pendingTau);
      table->r.push_back(std::vector<double>()); table->r.back().swap(pendingR);
    }
    pendingTau.clear(); pendingR.clear(); pendingEnergy = -1.0;
  };
  while (std::getline(in, line)) {
    if (!pastHeader) {
      if (line.find("-----") != std::string::npos) pastHeader = true;
      continue;
    }
    size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) { flush(); if (haveAbove) break; continue; }
    if (line[first] == '#') continue;
    const char *s = line.c_str() + first;
    char *end = nullptr;
    double a = std::strtod(s, &end);
    if (end == s) continue;
    while (*end == ' ' || *end == '\t') end++;
    if (*end == '\0' || *end == '#' || *end == '\r' || *end == '!') {
      flush();
      if (haveAbove) break;
      pendingEnergy = a * 1.0e-6;  // eV -> MeV
      continue;
    }
    char *end2 = nullptr;
    double b = std::strtod(end, &end2);
    if (end2 == end) continue;
    if (pendingEnergy >= 0.0) { pendingTau.push_back(a); pendingR.push_back(b); }
  }
  flush();
  if (table->energy.empty()) {
    log << "ERROR: no resolution functions between " << keepMin << " and " << keepMax
        << " MeV in " << udrFile_ << std::endl;
    return false;
  }
  log << "Read " << table->energy.size() << " user-defined resolution functions from "
      << udrFile_ << " (" << table->energy.front() << " - " << table->energy.back() << " MeV)" << std::endl;
  udrTable_ = table;
  return true;
}

/*!
 * Linear interpolation of the two bracketing tabulated functions (both the
 * delay grid and the density, as SAMMY does), then centring, burst/channel
 * convolution and normalisation.
 */

std::shared_ptr<const UdrKernel> TargetEffect::BuildUdrKernel(double eLab, double mass) const {
  std::shared_ptr<UdrKernel> k = std::make_shared<UdrKernel>();
  if (!udrTable_ || udrTable_->energy.empty()) return k;
  const UdrTable &t = *udrTable_;
  size_t hi = 0;
  while (hi < t.energy.size() && t.energy[hi] < eLab) hi++;
  if (hi == 0 || hi == t.energy.size()) {
    size_t i = (hi == 0) ? 0 : t.energy.size() - 1;
    k->tau = t.tau[i];
    k->r = t.r[i];
  } else {
    size_t lo = hi - 1;
    double f = (eLab - t.energy[lo]) / (t.energy[hi] - t.energy[lo]);
    size_t n = std::min(t.tau[lo].size(), t.tau[hi].size());
    k->tau.resize(n);
    k->r.resize(n);
    for (size_t j = 0; j < n; j++) {
      k->tau[j] = (1.0 - f) * t.tau[lo][j] + f * t.tau[hi][j];
      k->r[j] = (1.0 - f) * t.r[lo][j] + f * t.r[hi][j];
    }
  }
  const size_t n = k->tau.size();
  if (n < 2) return k;
  for (size_t j = 0; j < n; j++) if (k->r[j] < 0.0) k->r[j] = 0.0;
  // Centroid of the tabulated function (trapezoid on the piecewise-linear
  // density); SAMMY realigns the function so that its centroid is at zero,
  // the mean delay being part of the nominal flight path.
  if (udrCentred_) {
    double norm = 0.0, first = 0.0;
    for (size_t j = 0; j + 1 < n; j++) {
      double dt = k->tau[j + 1] - k->tau[j];
      norm += 0.5 * (k->r[j] + k->r[j + 1]) * dt;
      first += dt * (k->r[j] * (2.0 * k->tau[j] + k->tau[j + 1]) + k->r[j + 1] * (k->tau[j] + 2.0 * k->tau[j + 1])) / 6.0;
    }
    if (norm > 0.0) {
      double centroid = first / norm;
      for (size_t j = 0; j < n; j++) k->tau[j] -= centroid;
    }
  }
  // Burst (Gaussian, FWHM) and channel width (rectangular) in time, both in
  // nanoseconds on the line; their convolution has the closed form
  // [erf((x + c/2)/(s sqrt2)) - erf((x - c/2)/(s sqrt2))] / (2c).
  const double sigma = udrBurstFwhm_ * 1.0e-3 / 2.3548200450309493;
  // A negative channel width means |n| bins per decade of energy, the
  // logarithmic binning of n_TOF data: a bin then spans t ln(10)/(2n) in
  // time at this point's flight time t (SAMMY manual eq. III C3 a.14).
  double width = udrChannelWidth_ * 1.0e-3;
  if (udrChannelWidth_ < 0.0)
    width = UdrKinematics::TimeOfFlight(eLab, mass, udrFlightPath_) * std::log(10.0) / (2.0 * std::fabs(udrChannelWidth_));
  if (sigma > 0.0 || width > 0.0) {
    double dt = (k->tau.back() - k->tau.front()) / (n - 1);
    double half = 5.0 * sigma + 0.5 * width;
    int m = static_cast<int>(std::ceil(half / dt));
    std::vector<double> h(2 * m + 1);
    double hsum = 0.0;
    for (int i = -m; i <= m; i++) {
      double x = i * dt, v;
      if (sigma > 0.0 && width > 0.0)
        v = 0.5 * (std::erf((x + 0.5 * width) / (sigma * std::sqrt(2.0))) - std::erf((x - 0.5 * width) / (sigma * std::sqrt(2.0)))) / width;
      else if (sigma > 0.0)
        v = std::exp(-0.5 * x * x / (sigma * sigma)) / (sigma * std::sqrt(2.0 * pi));
      else
        v = (std::fabs(x) <= 0.5 * width) ? 1.0 / width : 0.0;
      h[i + m] = v;
      hsum += v;
    }
    for (double &v : h) v /= hsum;  // discrete normalisation: the response sums to one
    std::vector<double> tau2(n + 2 * m), r2(n + 2 * m, 0.0);
    for (size_t j = 0; j < tau2.size(); j++) tau2[j] = k->tau.front() + (static_cast<double>(j) - m) * dt;
    for (size_t j = 0; j < n; j++) {
      if (k->r[j] == 0.0) continue;
      for (int i = -m; i <= m; i++) r2[j + m + i] += k->r[j] * h[i + m];
    }
    k->tau.swap(tau2);
    k->r.swap(r2);
  }
  // A delay cannot approach the flight time itself: t0 - tau is the true
  // flight time, and as it goes to zero the true energy diverges.  Drop the
  // part of the tabulation that would put the true energy above 100 times
  // the nominal one (no real resolution function has weight there; a badly
  // scaled or uncentred table can), so the sub-point window stays finite.
  {
    const double t0 = UdrKinematics::TimeOfFlight(eLab, mass, udrFlightPath_);
    const double tauCap = t0 - UdrKinematics::TimeOfFlight(100.0 * eLab, mass, udrFlightPath_);
    for (size_t j = 0; j < k->tau.size(); j++) if (k->tau[j] > tauCap) k->r[j] = 0.0;
  }
  // Support and normalisation.
  double peak = 0.0;
  for (double v : k->r) peak = std::max(peak, v);
  size_t a = 0, b = k->r.size() - 1;
  while (a < b && k->r[a] <= 1.0e-9 * peak) a++;
  while (b > a && k->r[b] <= 1.0e-9 * peak) b--;
  k->tauMin = k->tau[a];
  k->tauMax = k->tau[b];
  double norm = 0.0;
  for (size_t j = 0; j + 1 < k->tau.size(); j++) norm += 0.5 * (k->r[j] + k->r[j + 1]) * (k->tau[j + 1] - k->tau[j]);
  if (norm > 0.0) for (double &v : k->r) v /= norm;
  return k;
}

/*!
 * Beam profile weight at an energy: the weighted sum of the skewed-Gaussian
 * components,  G(E|xi,omega,alpha) = exp(-((E-xi)/omega)^2/2) /(omega sqrt(2 pi))
 * * (1 + erf(alpha (E-xi)/(omega sqrt 2))),  each optionally zeroed outside
 * its mean +- nCut standard deviations.
 */

double TargetEffect::BeamProfileWeight(double energy) const {
  double weight = 0.0;
  for (std::vector<BeamProfileComponent>::const_iterator c = beamProfile_.begin(); c != beamProfile_.end(); c++) {
    if (!(c->omega > 0.0)) continue;
    double z = (energy - c->xi) / c->omega;
    if (beamTruncation_ > 0.0) {
      double delta = c->alpha / std::sqrt(1.0 + c->alpha * c->alpha);
      double mean = c->xi + c->omega * delta * std::sqrt(2.0 / pi);
      double sd = c->omega * std::sqrt(1.0 - 2.0 * delta * delta / pi);
      if (std::fabs(energy - mean) > beamTruncation_ * sd) continue;
    }
    weight += c->weight * std::exp(-0.5 * z * z) / (c->omega * std::sqrt(2.0 * pi)) * (1.0 + std::erf(c->alpha * z / std::sqrt(2.0)));
  }
  return weight;
}

/*!
 * Energy interval that contains all of the beam profile: the union of
 * xi +- 5 omega over the components (tighter if the profile is truncated).
 */

void TargetEffect::BeamProfileSupport(double &low, double &high) const {
  low = 1.0e300;
  high = -1.0e300;
  for (std::vector<BeamProfileComponent>::const_iterator c = beamProfile_.begin(); c != beamProfile_.end(); c++) {
    double lo = c->xi - 5.0 * c->omega;
    double hi = c->xi + 5.0 * c->omega;
    if (beamTruncation_ > 0.0) {
      double delta = c->alpha / std::sqrt(1.0 + c->alpha * c->alpha);
      double mean = c->xi + c->omega * delta * std::sqrt(2.0 / pi);
      double sd = c->omega * std::sqrt(1.0 - 2.0 * delta * delta / pi);
      lo = std::max(lo, mean - beamTruncation_ * sd);
      hi = std::min(hi, mean + beamTruncation_ * sd);
    }
    low = std::min(low, lo);
    high = std::max(high, hi);
  }
}

double TargetEffect::GetBeamTpcSigma() const {
  return beamTpcSigma_;
}

double TargetEffect::GetBeamTruncation() const {
  return beamTruncation_;
}

bool TargetEffect::IsBeamPhotodissociation() const {
  return beamPhotodissociation_;
}

/*!
 * The c.m. beam-spread sigma at c.m. energy eCM for a point whose entrance pair has the c.m./lab factor
 * cmFactor.  The sigma and a `sigmalab` expression are kept in the lab frame of the input file and
 * converted here, per point: every segment listed on a <targetInt> line, and every component of an
 * advanced (SUM/RATIO or natural-target) segment, points to the same TargetEffect, and up to 2026-10-27
 * the sigma was converted once with the factor of the first segment on the line (and before 2026-10-01
 * once per segment), wrong for the other pairs by the ratio of their mass ratios.
 */

double TargetEffect::SigmaCM(double eCM, const Config &configure, double cmFactor) const {
  if (isConvCoefficients_ && convolutionEq_.GetEquation() != "") {
    if (sigmaLab_) return cmFactor * convolutionEq_.Evaluate(configure, eCM / cmFactor);
    return convolutionEq_.Evaluate(configure, eCM);
  }
  return cmFactor * sigma_;
}

double TargetEffect::GaussianKernel(double e, double centroid, double sigma) {
  return pow(2. * pi, -0.5) / sigma * exp(-pow(e - centroid, 2.0) / 2.0 / pow(sigma, 2.0));
}

double TargetEffect::StragglingSigma(double deltaE) const {
  if (!isStraggling_ || deltaE <= 0.0) return 0.0;
  return stragglingCoefficient_ * std::sqrt(deltaE * 1000.0) / 1000.0;
}

/*!
 * Bohr's energy straggling after traversing an areal density N x (atoms/cm^2) of the target:
 * Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x (N. Bohr, Phil. Mag. 30 (1915) 581), a lab width converted to the
 * c.m. with the entrance pair's mass ratio.  e^2 = alpha hbar c in MeV fm, 1 fm = 1e-13 cm.  Z2 is the
 * charge of the target nucleus: for a compound target whose stopping cross section is per molecule or per
 * active atom (see the user guide) the electrons of the other constituents are not counted, and an
 * explicit coefficient should be entered instead.
 */

double TargetEffect::BohrStragglingSigma(double arealDensity, double cmFactor, int z1, int z2) const {
  if (!this->IsBohrStraggling() || !(arealDensity > 0.0) || z1 == 0 || z2 <= 0) return 0.0;
  const double e2 = fstruc * hbarc * 1.0e-13;  // MeV cm
  const double omega2 = 4.0 * pi * z1 * z1 * z2 * e2 * e2 * arealDensity;  // lab MeV^2
  return cmFactor * std::sqrt(omega2);
}

/*!
 * The sub-point grid of a point reaches forwardDepth above and backwardDepth below its energy. A
 * Gaussian beam spread needs convolutionRange sigma on each side; target integration needs the target
 * thickness below, and with a beam spread the window is padded by a further factor 5 (kept from the
 * original code, which sized the grid that way for the nested integral); straggling adds its own
 * convolutionRange sigma at the back of the target. One function for what used to be written out
 * eight times in WorkingCopy.cpp.
 */

void TargetEffect::GridExtent(double eCM, double thicknessCM, const Config &configure, double cmFactor, int z1, int z2, double &forward, double &backward) const {
  forward = 0.0;
  backward = 0.0;
  double beam = this->IsGaussian() ? convolutionRange * this->SigmaCM(eCM, configure, cmFactor) : 0.0;
  if (this->IsTargetIntegration()) {
    backward = thicknessCM + (this->IsGaussian() ? 5.0 * beam : 0.0);
    forward = this->IsGaussian() ? 5.0 * beam : 0.0;
    if (this->IsBohrStraggling()) {
      // the whole target: its areal density (the largest a fitted density is laid out for), or for a thick
      // target the thickness over the surface stopping power
      double areal = this->GetLayoutDensity();
      if (isThickTarget_) {
        const double eps = cmFactor * stoppingPowerEq_.Evaluate(configure, eCM / cmFactor);
        areal = (eps > 0.0) ? thicknessCM / eps : 0.0;
      }
      backward += convolutionRange * this->BohrStragglingSigma(areal, cmFactor, z1, z2);
    } else {
      backward += convolutionRange * this->StragglingSigma(thicknessCM);
    }
  } else if (this->IsGaussian()) {
    backward = beam;
    forward = beam;
  }
}

const std::vector<TargetEffect::BeamProfileComponent> &TargetEffect::GetBeamProfile() const {
  return beamProfile_;
}

/*!
 * Returns true if the target effect contains attenuation coefficients, otherwise
 * returns false.
 */

bool TargetEffect::IsQCoefficients() const {
  return isQCoefficients_;
}

/*!
 * Returns true if the target effect contains convolution coefficients, otherwise
 * returns false.
 */

bool TargetEffect::IsConvCoefficients() const {
  return isConvCoefficients_;
}

/*!
 * Returns the number of sub-points specified for the target effect in
 * the input file.
 */

int TargetEffect::NumSubPoints() const {
  return numIntegrationPoints_;
}

/*!
 * Returns the number of attenuation coefficients for the target effect in
 * the input file.
 */

int TargetEffect::NumQCoefficients() const {
  return qCoefficients_.size();
}

/*!
 * Returns the number of convolution coefficients for the target effect in
 * the input file.
 */

int TargetEffect::NumConvCoefficients() const {
  return convCoefficients_.size();
}

/*!
 * Returns the sigma of the Guassian for beam convolution.
 */

double TargetEffect::GetSigma() const {
  return sigma_;
}

/*!
 * Returns the density of the target in atoms/cm^2.  Only needed for
 * target integration, not Gaussian beam convolution.
 */

double TargetEffect::GetDensity() const {
  return density_;
}

/*!
 * The energy loss Delta (lab MeV) of the beam across the target at lab energy E: the solution of
 *   int_{E - Delta}^{E} dE' / epsilon(E') = N,
 * the areal density N (atoms/cm^2) being reached when the beam has crossed the whole target.  Up to
 * 2026-10-27 this was the first-order epsilon(E) N, which with a stopping cross section that rises through
 * the target (every ion above the Bragg peak) misplaces the back edge by 2-8 % for targets of 5-20 % energy
 * loss and, with 1/epsilon in the yield integrand, under-weights a flat cross section by the same amount
 * (audit 2026-10-25, tier 3 entry 2).  Newton on A(Delta) = int - N with A' = 1/epsilon(E - Delta), the
 * integral by Simpson's rule on 16 intervals, from the first-order estimate (exact for a constant epsilon).
 * A non-positive epsilon(E) is returned as epsilon(E) N, i.e. <= 0, so the caller's error check fires as
 * before; a target the beam cannot cross returns E itself (the caller floors the grid at 1 keV).
 */

double TargetEffect::TargetThickness(double energy, const Config &configure) const {
  return this->TargetThickness(energy, configure, this->GetDensity());
}

double TargetEffect::TargetThickness(double energy, const Config &configure, double density) const {
  const double eps0 = stoppingPowerEq_.Evaluate(configure, energy);
  if (!(eps0 > 0.0) || !(density > 0.0) || !(energy > 0.0)) return eps0 * density;
  auto invEps = [&](double e) {
    const double eps = stoppingPowerEq_.Evaluate(configure, e);
    return (eps > 0.0) ? 1.0 / eps : 0.0;
  };
  double delta = eps0 * density;
  for (int iter = 0; iter < 30; iter++) {
    if (!(delta < energy)) return energy;
    const int n = 16;
    const double h = delta / n;
    double sum = invEps(energy) + invEps(energy - delta);
    for (int k = 1; k < n; k++) sum += ((k % 2) ? 4.0 : 2.0) * invEps(energy - k * h);
    const double residual = sum * h / 3.0 - density;
    const double slope = invEps(energy - delta);
    if (!(slope > 0.0)) break;
    double step = residual / slope;
    // never jump past the beam energy or below a tenth of the current estimate
    if (delta - step > energy) step = delta - energy;
    if (delta - step < 0.1 * delta) step = 0.9 * delta;
    delta -= step;
    if (std::fabs(step) <= 1.0e-7 * delta) break;
  }
  return delta;
}

double TargetEffect::MeanEnergyLoss(double energy, const Config &configure, double density) const {
  if (isThickTarget_) return 0.0;
  const double delta = this->TargetThickness(energy, configure, density);
  if (!(delta > 0.0) || !(delta <= energy)) return 0.0;
  auto invEps = [&](double e) {
    const double eps = stoppingPowerEq_.Evaluate(configure, e);
    return (eps > 0.0) ? 1.0 / eps : 0.0;
  };
  // Simpson on 32 panels: int (E - e)/epsilon de / int 1/epsilon de over [E - Delta, E]
  const int n = 32;
  const double h = delta / n;
  double num = 0.0, den = 0.0;
  for (int k = 0; k <= n; k++) {
    const double w = (k == 0 || k == n) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
    const double depth = k * h;
    const double ie = invEps(energy - depth);
    num += w * depth * ie;
    den += w * ie;
  }
  return (den > 0.0) ? num / den : 0.5 * delta;
}

/*!
 * Makes this line a member of a density group whose first line is \p leader: its density, prior and layout
 * density become the leader's (the group has one parameter, held by the leader).
 */

void TargetEffect::JoinDensityGroup(const TargetEffect &leader) {
  hasDensityFit_ = true;
  densityVary_ = leader.densityVary_;
  density_ = leader.density_;
  densityNominal_ = leader.densityNominal_;
  densitySigma_ = leader.densitySigma_;
  layoutDensity_ = leader.layoutDensity_;
}

/*!
 * Returns the attenuation coefficients for the given order specified in by the target effect.
 */

double TargetEffect::GetQCoefficient(int order) const {
  return (qCoefficients_.size() > order) ? qCoefficients_[order] : 1.;
}

/*!
 * Returns the convolution coefficients for the given order specified in by the target effect.
 */

double TargetEffect::GetConvCoefficient(int order) const {
  return (convCoefficients_.size() > order) ? convCoefficients_[order] : 1.;
}

/*!
 * Sets the convolution sigma to a new value.
 */

void TargetEffect::SetSigma(double sigma) {
  sigma_ = sigma;
}

/*!
 * Sets the target density to a new value.
 */

void TargetEffect::SetDensity(double density) {
  density_ = density;
}

/*!
 * Sets the number of sub-points for the TargetEffect object.
 */

void TargetEffect::SetNumSubPoints(int numPoints) {
  numIntegrationPoints_ = numPoints;
}

/*!
 * Parses and returns a vector of integers corresponding to the
 * segment list specified as a string.  The segments list contains the
 * segments for which the target effect is applicable.
 */

std::vector<int> TargetEffect::GetSegmentsList() const {
  return datasetKeys_;
}

/*!
 * Returns the Equation object corresponding to the parametrized stopping
 * cross section.
 */

const Equation *TargetEffect::GetStoppingPowerEq() const { return &stoppingPowerEq_; }

/*!
 * Returns the Equation object corresponding to the parametrized convolution.
 */

Equation *TargetEffect::GetConvolutionEq() {
  Equation *tempEquation;
  tempEquation = &convolutionEq_;
  return tempEquation;
}

/*!
 * Returns true if the target effect includes energy straggling,
 * otherwise returns false.
 */

bool TargetEffect::IsStraggling() const {
  return isStraggling_;
}

/*!
 * Returns the straggling coefficient c_s (keV^1/2, c.m.): sigma_cm = c_s sqrt(dE_cm [keV]); 0 = Bohr straggling.
 */

double TargetEffect::GetStragglingCoefficient() const {
  return stragglingCoefficient_;
}

/*!
 * Returns the resonance width multiplier for the adaptive integration grid.
 * The half-extent, in widths, of the uniform core of each resonance's lattice
 * (geometric tails continue beyond it; see AdaptiveIntegrationGrid).
 */

double TargetEffect::GetResonanceWidthMultiplier() const {
  // A value on the targetInt line always wins.  Left out, a pure Gaussian
  // convolution gets the smaller core it needs (see
  // AdaptiveIntegrationGrid::convolutionCoreWidths); every other effect keeps
  // the general default.
  if (!explicitResonanceWidthMultiplier_ && (isConvolution_ || isConvCoefficients_) &&
      !isTargetIntegration_ && !isBeamProfile_)
    return AdaptiveIntegrationGrid::convolutionCoreWidths;
  return resonanceWidthMultiplier_;
}

/*!
 * Returns the number of integration points per resonance width for the adaptive grid.
 */

double TargetEffect::GetPointsPerWidth() const {
  return pointsPerWidth_;
}

/*!
 * Returns the optional lab-energy windows the effect is restricted to.
 * An empty vector means the effect covers its segments completely.
 */

const std::vector<std::pair<double, double> > &TargetEffect::GetRanges() const {
  return ranges_;
}

/*!
 * Returns the width (MeV, lab energy) of the smooth blend applied at each
 * range edge.  Zero means the effect switches on and off abruptly.
 */

double TargetEffect::GetTransitionWidth() const {
  return transitionWidth_;
}

/*!
 * Returns the relative tolerance of the automatic per-point application.
 * Zero disables the automatic decision: the effect is always integrated.
 */

double TargetEffect::GetAutoTolerance() const {
  return autoTolerance_;
}

/*!
 * Blend weight in [0,1] at the given lab energy.  With no ranges declared the
 * weight is always one.  With ranges and a zero transition width the weight is
 * a hard 1 inside any window and 0 outside; a positive transition width turns
 * each edge into a smoothstep ramp centred on the edge, so the modelled
 * observable passes continuously between the convolved and the bare curve.
 */

double TargetEffect::BlendWeight(double labEnergy) const {
  if (ranges_.empty()) return 1.0;
  return BlendWeight(labEnergy, ranges_);
}

double TargetEffect::BlendWeight(double labEnergy, const std::vector<std::pair<double, double> > &windows) const {
  if (windows.empty()) return autoWindows_ ? 0.0 : 1.0;
  double weight = 0.0;
  for (std::vector<std::pair<double, double> >::const_iterator r = windows.begin(); r != windows.end(); r++) {
    double w;
    if (transitionWidth_ <= 0.) {
      w = (labEnergy >= r->first && labEnergy <= r->second) ? 1.0 : 0.0;
    } else {
      double up = (labEnergy - (r->first - 0.5 * transitionWidth_)) / transitionWidth_;
      double down = ((r->second + 0.5 * transitionWidth_) - labEnergy) / transitionWidth_;
      up = std::min(1.0, std::max(0.0, up));
      down = std::min(1.0, std::max(0.0, down));
      up = up * up * (3.0 - 2.0 * up);
      down = down * down * (3.0 - 2.0 * down);
      w = up * down;
    }
    if (w > weight) weight = w;
  }
  return weight;
}
