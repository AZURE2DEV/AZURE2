#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <iosfwd>
#include "Constants.h"
#include "Equation.h"
#include "Straggling.h"
#include "UdrKernel.h"

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

  /// Optional lab-energy windows the effect is restricted to; empty = whole segment.
  const std::vector<std::pair<double, double> > &GetRanges() const;
  /// Width (MeV, lab) of the smooth blend at each range edge; 0 = hard edges.
  double GetTransitionWidth() const;
  /// Relative tolerance for automatic per-point application; 0 = always apply.
  double GetAutoTolerance() const;
  /// Blend weight in [0,1] at the given lab energy: 1 = effect fully applied,
  /// 0 = point untouched.  Always 1 when no ranges are given.
  double BlendWeight(double labEnergy) const;

  /// The multiple of sigma above and below centroid energy to use as integration range
  static constexpr double convolutionRange = 3.;

  /// One skewed-Gaussian component of a beam energy profile.
  struct BeamProfileComponent {
    double xi;      ///< location parameter
    double omega;   ///< scale parameter
    double alpha;   ///< shape (skewness) parameter
    double weight;  ///< relative weight (luminosity) of the component
  };
  /// Beam-profile kernel: an absolute (not point-centred) beam energy profile,
  /// a Gaussian detector-resolution window per point and, optionally, the
  /// detailed-balance weight of an inverse photodissociation measurement.
  /// Written as the trailing block `beamprofile N {xi omega alpha w}xN
  /// tpcSigma nCut dbFlag` of a targetInt line, energies in the lab frame.
  bool IsBeamProfile() const;
  /// Any effect whose yield is an integral over sub-points.
  bool IsSubPointEffect() const;
  /// Beam profile weight (sum of the skewed Gaussians) at an energy.
  double BeamProfileWeight(double energy) const;
  /// Energy range outside which the beam profile is negligible.
  void BeamProfileSupport(double &low, double &high) const;
  /// Sigma of the Gaussian detector (TPC) energy resolution.
  double GetBeamTpcSigma() const;
  /// Truncate each component outside mean +- this many standard deviations; 0 = none.
  double GetBeamTruncation() const;
  /// Weight the integrand by the detailed-balance factor of the inverse reaction.
  bool IsBeamPhotodissociation() const;
  /// Scale every beam-profile energy (xi, omega, tpc sigma) once, lab -> c.m.
  void ConvertBeamProfileToCM(double factor);
  /// Scale the fixed convolution sigma once, lab -> c.m.
  void ConvertSigmaToCM(double factor);
  const std::vector<BeamProfileComponent> &GetBeamProfile() const;

  /// User-defined numerical resolution function (SAMMY's UDR): a file of
  /// time-of-flight delay distributions tabulated at a set of energies, in
  /// SAMMY's UDR format, convolved with a Gaussian burst and a rectangular
  /// channel width.  Written as the trailing block
  /// `udr "<file>" flightPath_m burstFWHM_ns channel_ns centred` of a
  /// targetInt line; a negative channel value is |n| bins per decade of
  /// energy (logarithmic binning), i.e. a width of t ln(10)/(2n) in time.
  bool IsUdr() const;
  const std::string &GetUdrFile() const;
  double GetUdrFlightPath() const;
  double GetUdrBurstFwhm() const;
  double GetUdrChannelWidth() const;
  bool IsUdrCentred() const;
  /// Read the tabulated functions whose energies lie within a decade of
  /// [eMin, eMax] (lab MeV), plus the bracketing ones; merges with what is
  /// already loaded.  Returns false if the file cannot be read.
  bool LoadUdrTable(double eMinLab, double eMaxLab, std::ostream &log);
  /// The kernel for a point at lab energy eLab (MeV) for a projectile of rest
  /// energy mass (MeV): the table interpolated in energy, centred if asked,
  /// convolved with the burst and channel width, normalised.
  std::shared_ptr<const UdrKernel> BuildUdrKernel(double eLab, double mass) const;

 private:
  struct UdrTable {
    std::string file;
    std::vector<double> energy;                ///< lab MeV, ascending
    std::vector<std::vector<double> > tau;     ///< microseconds, per energy
    std::vector<std::vector<double> > r;       ///< per microsecond, per energy
    double loadedMin = 0.0, loadedMax = 0.0;   ///< lab MeV range the table covers
  };
  bool isUdr_ = false;
  std::string udrFile_;
  double udrFlightPath_ = 0.0;   ///< metres
  double udrBurstFwhm_ = 0.0;    ///< nanoseconds
  double udrChannelWidth_ = 0.0; ///< nanoseconds
  bool udrCentred_ = true;
  std::shared_ptr<UdrTable> udrTable_;  ///< shared between copies of the effect

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

  // Optional restriction of the effect to lab-energy windows, with a smooth
  // blend of the given width at each edge, and/or automatic per-point
  // application controlled by a relative tolerance.
  std::vector<std::pair<double, double> > ranges_;
  double transitionWidth_;
  double autoTolerance_;

  // Beam-profile kernel (see IsBeamProfile).
  bool isBeamProfile_ = false;
  std::vector<BeamProfileComponent> beamProfile_;
  double beamTpcSigma_ = 0.0;
  double beamTruncation_ = 0.0;
  bool beamPhotodissociation_ = false;
  bool beamProfileConverted_ = false;
  bool sigmaConverted_ = false;
};

#endif
