#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>			       
#include <vector>
#include "Constants.h"
#include "Equation.h"
#include "Straggling.h"

///An AZURE target effect entry

/*!
 * Experimential effects including gaussian beam convolution, target 
 * integration, and a combination of the two are grouped under the TargetEffect
 * class. An object is created corresponding to each corresponding entry in 
 * AZURESetup2.
 */

class TargetEffect {
 public:
  TargetEffect(std::istream &, const Config&);
  bool IsActive() const;
  bool IsConvolution() const;
  bool IsTargetIntegration() const;
  bool IsQCoefficients() const;
  bool IsConvCoefficients() const;
  int NumSubPoints() const;
  int NumQCoefficients() const;
  int NumConvCoefficients() const;
  double GetSigma() const;
  double CalculateSigma(double,const Config&);
  double GetDensity() const;
  double TargetThickness(double,const Config&);
  double GetConvolutionFactor(double, double) const;
  double CalculateConvolutionFactor(double, double, const Config&);
  double GetQCoefficient(int) const;
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

  ///The multiple of sigma above and below centroid energy to use as integration range
  static constexpr double convolutionRange=3.;
  /*Wider window for THM (HOES) segments: their extremely narrow resonances make
  the smeared curve in the valleys a pure Gaussian tail of a neighboring spike,
  so a +-3 sigma cutoff loses up to ~25% there; +-5 sigma keeps the truncation
  error below ~1e-3.*/
  static constexpr double thmConvolutionRange=5.;
  /*Numerical floor (in MeV) kept between a target-effect sub-point energy and the
  channel threshold, since Coulomb functions and penetrabilities are singular at
  E=0. For ordinary segments this bounds the channel CM energy directly; for THM
  (HOES) segments, which legitimately probe negative channel CM energies, it
  instead bounds the compound-system energy (CMEnergy+SepE+ExE) so the floor is
  applied relative to the physical threshold rather than an arbitrary channel origin.*/
  static constexpr double minIntegrationEnergy=0.001;
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
