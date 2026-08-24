#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>
#include <vector>
#include "Constants.h"
#include "Equation.h"
#include "Straggling.h"

/// An AZURE target effect entry

/*!
 * Experimential effects including gaussian beam convolution, target
 * integration, and a combination of the two are grouped under the TargetEffect
 * class. An object is created corresponding to each corresponding entry in
 * AZURESetup2.
 */

class TargetEffect {
 public:
  /// Read one target effect from the target-effects input file.
  TargetEffect(std::istream &, const Config &);
  /// Was this effect marked active in the input file?
  bool IsActive() const;
  /// Does it include Gaussian beam convolution?
  bool IsConvolution() const;
  /// Does it integrate over the target thickness?
  bool IsTargetIntegration() const;
  /// Does it carry angular attenuation coefficients?
  bool IsQCoefficients() const;
  /// Does it carry convolution coefficients?
  bool IsConvCoefficients() const;
  /// Sub-points the integral is sampled on.
  int NumSubPoints() const;
  /// Number of attenuation coefficients.
  int NumQCoefficients() const;
  /// Number of convolution coefficients.
  int NumConvCoefficients() const;
  /// Gaussian beam-spread sigma.
  double GetSigma() const;
  /// Beam sigma at a given energy, where it is energy dependent.
  double CalculateSigma(double, const Config &);
  /// Target density in atoms/cm^2. Target integration only.
  double GetDensity() const;
  /// Energy loss across the target at a given energy, from the stopping cross section and the density.
  double TargetThickness(double, const Config &);
  /// Weight of one sub-point in the convolution integrand.
  double GetConvolutionFactor(double, double) const;
  /// As above with an energy-dependent sigma.
  double CalculateConvolutionFactor(double, double, const Config &);
  /// Attenuation coefficient of the given order.
  double GetQCoefficient(int) const;
  /// Convolution coefficient of the given order.
  double GetConvCoefficient(int) const;
  void SetSigma(double);
  void SetDensity(double);
  void SetNumSubPoints(int);
  std::vector<int> GetSegmentsList() const;
  Equation *GetStoppingPowerEq();
  Equation *GetConvolutionEq();

  // ERYA straggling integration methods
  bool IsStraggling() const;
  double GetStragglingCoefficient() const;

  // Adaptive integration grid parameters
  double GetResonanceWidthMultiplier() const;
  double GetPointsPerWidth() const;

  /// The multiple of sigma above and below centroid energy to use as integration range
  static constexpr double convolutionRange = 3.;
  /*Wider window for THM (HOES) segments: their extremely narrow resonances make
  the smeared curve in the valleys a pure Gaussian tail of a neighboring spike,
  so a +-3 sigma cutoff loses up to ~25% there; +-5 sigma keeps the truncation
  error below ~1e-3.*/
  static constexpr double thmConvolutionRange = 5.;
  /*Numerical floor (in MeV) kept between a target-effect sub-point energy and the
  channel threshold, since Coulomb functions and penetrabilities are singular at
  E=0. For ordinary segments this bounds the channel CM energy directly; for THM
  (HOES) segments, which legitimately probe negative channel CM energies, it
  instead bounds the compound-system energy (CMEnergy+SepE+ExE) so the floor is
  applied relative to the physical threshold rather than an arbitrary channel origin.*/
  static constexpr double minIntegrationEnergy = 0.001;

 private:
  bool isConvolution_;
  bool isTargetIntegration_;
  bool isActive_;
  bool isQCoefficients_;
  bool isConvCoefficients_;
  int numIntegrationPoints_;
  double sigma_;
  double density_;
  Equation stoppingPowerEq_;
  std::string segmentsList_;
  vector_r qCoefficients_;
  vector_r convCoefficients_;
  Equation convolutionEq_;

  // ERYA straggling integration variables
  bool isStraggling_;
  double stragglingCoefficient_;

  // Adaptive integration grid parameters
  double resonanceWidthMultiplier_;
  double pointsPerWidth_;
};

#endif
