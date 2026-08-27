#include "TargetEffect.h"
#include <sstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include "Straggling.h"

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
  resonanceWidthMultiplier_ = 5.0;
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
