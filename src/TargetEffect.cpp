#include "TargetEffect.h"
#include <sstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include "Straggling.h"
#include <cmath>

/*!
 * Constructor reads directly from an std::ifstream pointing to the target
 * effect input file.  If a valid target effect is read, a TargetEffect object is
 * created.
 */

TargetEffect::TargetEffect(std::istream &stream, const Config &configure) {
  // Initialize straggling to safe defaults
  isStraggling_ = false;
  stragglingCoefficient_ = 0.04;

  // Initialize adaptive grid parameters to defaults
  // See AdaptiveIntegrationGrid::GridConfig for why this is 20 rather than 5.
  resonanceWidthMultiplier_ = 20.0;
  pointsPerWidth_ = 50.0;

  // By default the effect applies to every point of its segments, with hard
  // edges and no automatic per-point decision.
  transitionWidth_ = 0.0;
  autoTolerance_ = 0.0;

  int isActive;
  std::string segmentList;
  int numIntegrationPoints;
  int isConvolution;
  double sigma;
  int isTargetIntegration;
  double density;
  std::string stoppingPowerEq;
  int numParameters;
  vector_r parameters;
  int isQCoefficients;
  int numQCoefficients;
  vector_r qCoefficients;
  int isConvCoefficients;
  int numConvCoefficients;
  vector_r convCoefficients;
  std::string convolutionEq;

  stream >> isActive >> segmentList >> numIntegrationPoints >> isConvolution >> sigma >> isTargetIntegration >> density >> stoppingPowerEq >> numParameters;
  if (!stream.eof()) {
    for (int i = 0; i < numParameters; i++) {
      double tempParameter;
      stream >> tempParameter;
      parameters.push_back(tempParameter);
    }

    stream >> isQCoefficients >> numQCoefficients;
    for (int i = 0; i < numQCoefficients; i++) {
      double tempQCoefficient;
      stream >> tempQCoefficient;
      qCoefficients.push_back(tempQCoefficient);
    }
    isQCoefficients_ = (isQCoefficients == 1) ? true : false;
    qCoefficients_ = qCoefficients;

    stream >> isConvCoefficients >> convolutionEq >> numConvCoefficients;
    for (int i = 0; i < numConvCoefficients; i++) {
      double tempConvCoefficient;
      stream >> tempConvCoefficient;
      convCoefficients.push_back(tempConvCoefficient);
    }
    isConvCoefficients_ = (isConvCoefficients == 1) ? true : false;
    convCoefficients_ = convCoefficients;

    // Read straggling flag and coefficient (optional for backward compatibility)
    // Skip whitespace and check if there's a digit (not a '<' which would be </targetInt>)
    if (stream.good()) stream >> std::ws;  // Skip whitespace
    if (stream.good() && stream.peek() != '<' && std::isdigit(stream.peek())) {
      int isStraggling = 0;
      stream >> isStraggling;
      isStraggling_ = (isStraggling == 1);

      // Try to read coefficient if available
      if (stream.good()) stream >> std::ws;
      if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
        double stragglingCoeff = 0.04;
        stream >> stragglingCoeff;
        stragglingCoefficient_ = stragglingCoeff;

        // Try to read adaptive grid params (optional for backward compatibility)
        if (stream.good()) stream >> std::ws;
        if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
          double rwm;
          stream >> rwm;
          if (stream) {
            resonanceWidthMultiplier_ = rwm;
            if (stream.good()) stream >> std::ws;
            if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.' || stream.peek() == '-')) {
              double ppw;
              stream >> ppw;
              if (stream) pointsPerWidth_ = ppw;
            }
          }
        }
      }
    }

    // Optional restriction of the effect to lab-energy windows: a quoted
    // token "lo1-hi1,lo2-hi2" (empty "" means unrestricted), then optionally
    // the blend width at each edge (MeV) and the automatic-application
    // tolerance.  Old files stop before this token; old readers never look
    // past the fields they know, so the format stays compatible both ways.
    // Every probe is guarded on good(): the caller treats failbit as a
    // malformed line, and probing past the end of an exactly-consumed line
    // must not look like an error.
    if (stream.good()) {
      stream >> std::ws;
      if (stream.good() && stream.peek() == '"') {
        std::string rangesToken;
        stream >> rangesToken;
        size_t q = 0;
        while ((q = rangesToken.find('\"')) != std::string::npos) rangesToken.erase(q, 1);
        std::istringstream rs(rangesToken);
        std::string item;
        while (std::getline(rs, item, ',')) {
          size_t dash = item.find('-', 1);  // energies are positive; skip a leading sign
          if (dash == std::string::npos) continue;
          double lo = atof(item.substr(0, dash).c_str());
          double hi = atof(item.substr(dash + 1).c_str());
          if (hi > lo) ranges_.push_back(std::pair<double, double>(lo, hi));
        }
        if (stream.good()) {
          stream >> std::ws;
          if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.')) {
            double tw;
            stream >> tw;
            if (!stream.fail() && tw > 0.) transitionWidth_ = tw;
            if (stream.good()) {
              stream >> std::ws;
              if (stream.good() && stream.peek() != '<' && (std::isdigit(stream.peek()) || stream.peek() == '.')) {
                double tol;
                stream >> tol;
                if (!stream.fail() && tol > 0.) autoTolerance_ = tol;
              }
            }
          }
        }
      }
    }
    // Optional keyword blocks, introduced by a word so that older readers
    // (which only ever probe for digits or a quote) stop in front of them:
    //   beamprofile N {xi omega alpha weight}xN tpcSigma nCut dbFlag
    //   udr "<file>" flightPath_m burstFWHM_ns channel_ns centred
    while (stream.good()) {
      stream >> std::ws;
      if (!(stream.good() && std::isalpha(stream.peek()))) break;
      std::string keyword;
      stream >> keyword;
      if (keyword == "beamprofile") {
        int numComponents = 0;
        stream >> numComponents;
        for (int i = 0; i < numComponents && stream.good(); i++) {
          BeamProfileComponent c;
          stream >> c.xi >> c.omega >> c.alpha >> c.weight;
          beamProfile_.push_back(c);
        }
        int dbFlag = 0;
        stream >> beamTpcSigma_ >> beamTruncation_ >> dbFlag;
        beamPhotodissociation_ = (dbFlag == 1);
        isBeamProfile_ = !beamProfile_.empty() && !stream.fail();
        if (beamProfile_.empty()) stream.setstate(std::ios_base::failbit);
      } else if (keyword == "udr") {
        stream >> std::ws;
        std::string file;
        if (stream.peek() == '"') {
          stream.get();
          std::getline(stream, file, '"');
        } else {
          stream >> file;
        }
        int centred = 1;
        stream >> udrFlightPath_ >> udrBurstFwhm_ >> udrChannelWidth_ >> centred;
        udrFile_ = file;
        udrCentred_ = (centred != 0);
        isUdr_ = !udrFile_.empty() && udrFlightPath_ > 0.0 && !stream.fail();
        if (!isUdr_) stream.setstate(std::ios_base::failbit);
      } else {
        stream.setstate(std::ios_base::failbit);
      }
    }
    // A line consumed exactly to its end during the optional probes is not a
    // parse error; genuinely malformed fields fail without reaching eof.
    if (stream.fail() && stream.eof()) stream.clear(std::ios_base::eofbit);

    size_t found = 0;
    while (found != std::string::npos) {
      found = segmentList.find('\"');
      if (found != std::string::npos) segmentList.erase(found, 1);
    }
    found = 0;
    while (found != std::string::npos) {
      found = stoppingPowerEq.find('\"');
      if (found != std::string::npos) stoppingPowerEq.erase(found, 1);
    }
    found = 0;
    while (found != std::string::npos) {
      found = convolutionEq.find('\"');
      if (found != std::string::npos) convolutionEq.erase(found, 1);
    }
    if (isActive == 1)
      isActive_ = true;
    else
      isActive_ = false;
    segmentsList_ = segmentList;
    numIntegrationPoints_ = numIntegrationPoints;
    if (isConvolution == 1)
      isConvolution_ = true;
    else
      isConvolution_ = false;
    sigma_ = sigma;
    if (isTargetIntegration == 1)
      isTargetIntegration_ = true;
    else
      isTargetIntegration_ = false;
    density_ = density;
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
  return isConvolution_ || isTargetIntegration_ || isConvCoefficients_ || isBeamProfile_ || isUdr_;
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

bool TargetEffect::LoadUdrTable(double eMinLab, double eMaxLab, std::ostream &log) {
  if (!isUdr_) return false;
  const double keepMin = eMinLab / 10.0, keepMax = eMaxLab * 10.0;
  if (udrTable_ && udrTable_->file == udrFile_ && !udrTable_->energy.empty() &&
      udrTable_->loadedMin <= keepMin && udrTable_->loadedMax >= keepMax)
    return true;
  std::ifstream in(udrFile_.c_str());
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
 * Converts the beam-profile energies from the frame of the input file (lab)
 * to the centre of mass, exactly once, however many segments share the effect.
 */

void TargetEffect::ConvertBeamProfileToCM(double factor) {
  if (beamProfileConverted_) return;
  for (std::vector<BeamProfileComponent>::iterator c = beamProfile_.begin(); c != beamProfile_.end(); c++) {
    c->xi *= factor;
    c->omega *= factor;
  }
  beamTpcSigma_ *= factor;
  beamProfileConverted_ = true;
}

/*!
 * Converts the fixed convolution sigma from the frame of the input file (lab)
 * to the centre of mass, exactly once, however many segments share the effect.
 * Every segment listed on a <targetInt> line, and every component of an
 * advanced (SUM/RATIO) segment, points to the same TargetEffect; converting it
 * once per segment shrank the kernel by the mass ratio for each extra one.
 */

void TargetEffect::ConvertSigmaToCM(double factor) {
  if (sigmaConverted_) return;
  sigma_ *= factor;
  sigmaConverted_ = true;
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
 * Calculates the Target thickness from the stopping cross section and
 * the target density as a function of energy.
 */

double TargetEffect::TargetThickness(double energy, const Config &configure) {
  return this->GetStoppingPowerEq()->Evaluate(configure, energy) * this->GetDensity();
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
  std::vector<int> tempList;
  int i = 0;
  int lastSegNum = 0;
  bool inclusive = false;
  while (i < segmentsList_.length()) {
    if (segmentsList_[i] >= '0' && segmentsList_[i] <= '9') {
      std::string tempString;
      while (segmentsList_[i] != ',' && segmentsList_[i] != '-' &&
             i < segmentsList_.length()) {
        tempString += segmentsList_[i];
        i++;
      }
      std::istringstream stm;
      stm.str(tempString);
      int tempSegNum;
      stm >> tempSegNum;
      if (inclusive == true)
        for (int j = lastSegNum + 1; j <= tempSegNum; j++)
          tempList.push_back(j);
      else
        tempList.push_back(tempSegNum);
      lastSegNum = tempSegNum;
    }
    if (segmentsList_[i] == '-')
      inclusive = true;
    else
      inclusive = false;
    i++;
  }
  return tempList;
}

/*!
 * Returns the Equation object corresponding to the parametrized stopping
 * cross section.
 */

Equation *TargetEffect::GetStoppingPowerEq() {
  Equation *tempEquation;
  tempEquation = &stoppingPowerEq_;
  return tempEquation;
}

/*!
 * Returns the Equation object corresponding to the parametrized convolution.
 */

Equation *TargetEffect::GetConvolutionEq() {
  Equation *tempEquation;
  tempEquation = &convolutionEq_;
  return tempEquation;
}

/*!
 * Returns the multiplicative convolution factor for evaluation of the integrand
 * of a target effect.
 */

double TargetEffect::GetConvolutionFactor(double energy, double centroid) const {
  double sigma = this->GetSigma();
  return pow(2. * pi, -0.5) / sigma * exp(-pow(energy - centroid, 2.0) / 2.0 / pow(sigma, 2.0));
}

/*!
 * Calculates the sigma of the Gaussian beam convolution as a function of energy.
 */

double TargetEffect::CalculateSigma(double energy, const Config &configure) {
  double sigma = 0;
  // Check if convolution equation is defined
  if (this->GetConvolutionEq()->GetEquation() != "") {
    sigma = this->GetConvolutionEq()->Evaluate(configure, energy);
  } else
    sigma = this->GetSigma();
  return sigma;
}

/*!
 * Returns the multiplicative convolution factor for evaluation of the integrand
 * of a target effect with energy dependent sigma.
 */

double TargetEffect::CalculateConvolutionFactor(double energy, double centroid, const Config &configure) {
  double sigma = this->CalculateSigma(energy, configure);
  return pow(2. * pi, -0.5) / sigma * exp(-pow(energy - centroid, 2.0) / 2.0 / pow(sigma, 2.0));
}

/*!
 * Returns true if the target effect includes energy straggling,
 * otherwise returns false.
 */

bool TargetEffect::IsStraggling() const {
  return isStraggling_;
}

/*!
 * Returns the straggling coefficient for the target effect.
 * Typical value is 0.04 keV^0.5 per keV^0.5 of energy loss.
 */

double TargetEffect::GetStragglingCoefficient() const {
  return stragglingCoefficient_;
}

/*!
 * Returns the resonance width multiplier for the adaptive integration grid.
 * This controls how many total widths Γ are covered on each side of a resonance.
 */

double TargetEffect::GetResonanceWidthMultiplier() const {
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
  double weight = 0.0;
  for (std::vector<std::pair<double, double> >::const_iterator r = ranges_.begin(); r != ranges_.end(); r++) {
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
