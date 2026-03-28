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
