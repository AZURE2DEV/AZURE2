#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <iosfwd>
#include "Constants.h"
#include "Equation.h"
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
  /// Is there a Gaussian beam spread at all (a constant sigma or an expression)?
  bool IsGaussian() const { return isConvolution_ || isConvCoefficients_; }
  /// Is an expression sigma(E) written in the lab frame (E in lab MeV, sigma in lab MeV, the `sigmalab`
  /// keyword), as the constant sigma always is? Without the keyword an expression is the pre-2026-10
  /// convention: evaluated at the c.m. energy, its value the c.m. width.
  bool IsSigmaLab() const { return sigmaLab_; }
  /*!
   * Electron screening (`screening Ue vary sigma`, Ue and sigma in eV): every point and sub-point of the
   * effect's segments gets its bare cross section times exp(pi eta U_e / E).  U_e is a fit parameter
   * (EData::FillMnParams) when vary is 1, with a Gaussian prior of width sigma about the value read when
   * sigma > 0.  The line need not carry convolution or target integration.
   */
  bool HasScreening() const { return hasScreening_; }
  /// A `screening` keyword without its three numbers (or a negative sigma).
  bool HasScreeningError() const { return screeningError_; }
  /// Current U_e in eV (the value read, or the fit parameter's).
  double GetScreeningUe() const { return screeningUe_; }
  void SetScreeningUe(double ueEv) { screeningUe_ = ueEv; }
  /// U_e as read from the line, the centre of the prior, eV.
  double GetScreeningNominal() const { return screeningNominal_; }
  bool IsScreeningVaried() const { return screeningVary_; }
  /// Prior width in eV, 0 = no prior.
  double GetScreeningSigma() const { return screeningSigma_; }
  /// ((U_e - nominal) / sigma)^2, or 0 without a prior.
  double ScreeningPenalty() const {
    if (!hasScreening_ || screeningSigma_ <= 0.0) return 0.0;
    const double r = (screeningUe_ - screeningNominal_) / screeningSigma_;
    return r * r;
  }
  /// The c.m. beam-spread sigma at c.m. energy \p eCM: the one place every kernel takes its width from.
  /// A constant sigma was converted lab -> c.m. by ConvertSigmaToCM; a lab expression is evaluated at
  /// eCM / F and its value scaled by F, the same factor; a legacy expression is used as it stands.
  double SigmaCM(double eCM, const Config &) const;
  /// The normalised Gaussian kernel exp(-(e - centroid)^2 / 2 sigma^2) / (sigma sqrt(2 pi)).
  static double GaussianKernel(double e, double centroid, double sigma);
  /// Straggling sigma (c.m. MeV) after an energy loss of \p deltaE (c.m. MeV): c_s * sqrt(dE[keV]) keV.  The
  /// coefficient c_s (keV^1/2) is a CENTRE-OF-MASS quantity: the width it gives is the c.m. width of the
  /// energy distribution after a c.m. energy loss, with no frame factor (documented as such 2026-10-27; a
  /// coefficient describing the beam in the lab has to be multiplied by sqrt(m2/(m1+m2)) first).
  double StragglingSigma(double deltaE) const;
  /// Straggling on with a coefficient of 0 (the default) means Bohr straggling, computed from the charges
  /// and the stopping power instead of a coefficient: Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x, with N x the
  /// areal density (atoms/cm^2) traversed (int dE / epsilon), Z1 the beam and Z2 the target nucleus of the
  /// entrance pair, and the c.m. width m2/(m1+m2) Omega_lab.  (Before 2026-10-27 the default coefficient
  /// was 0.04 keV^1/2, several times below Bohr for light ions in light targets: essentially no straggling.)
  bool IsBohrStraggling() const { return isStraggling_ && !(stragglingCoefficient_ > 0.0); }
  /// The Bohr straggling sigma (c.m. MeV) after traversing \p arealDensity atoms/cm^2 of the target.
  double BohrStragglingSigma(double arealDensity) const;
  /// The entrance pair's charges (beam Z1, target Z2) for Bohr straggling; the first call wins, like
  /// ConvertSigmaToCM, since every segment on the line shares the effect.
  void SetBohrCharges(int z1, int z2);
  /// How far the sub-point grid of a point at c.m. energy \p eCM must reach above (\p forward) and below
  /// (\p backward) it for this effect; \p thicknessCM is the point's target thickness (c.m. MeV, 0 if none).
  void GridExtent(double eCM, double thicknessCM, const Config &, double &forward, double &backward) const;
  /// Target density in atoms/cm^2. Target integration only.
  double GetDensity() const;
  /// A target integration with a density of 0 is an infinitely THICK target: the yield per incident
  /// particle, int_0^{E0} sigma(E) / epsilon(E) dE, not divided by a density; every point integrates
  /// down to 1 keV.  (Before 2026-10-27 this was done by overwriting the effect's density with 1e24
  /// at the first zero-thickness point, which left that point with a thickness of 0 and a yield of
  /// exactly 0 under a beam spread or straggling, and also fired for a thin target whose stopping
  /// cross section came out <= 0; audit item 17.)
  bool IsThickTarget() const { return isThickTarget_; }
  /// Energy loss (lab MeV) across the target at lab energy E: Delta with int_{E-Delta}^{E} dE'/epsilon = N, the
  /// exact range integral (first-order epsilon(E) N before 2026-10-27, 2-8 % off for 5-20 % targets).
  double TargetThickness(double, const Config &);
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

  // Energy straggling (see StragglingSigma / IsBohrStraggling)
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
  static constexpr double convolutionRange = 5.;

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
  bool isThickTarget_ = false;   ///< target integration with density 0 (IsThickTarget)
  Equation stoppingPowerEq_;
  std::string segmentsList_;
  vector_r qCoefficients_;
  vector_r convCoefficients_;
  Equation convolutionEq_;
  bool sigmaLab_ = false;        ///< the expression is lab in, lab out (`sigmalab` keyword)
  bool hasScreening_ = false;     ///< `screening` keyword present
  bool screeningError_ = false;   ///< malformed `screening` block
  double screeningUe_ = 0.0;      ///< eV
  double screeningNominal_ = 0.0; ///< eV, as read
  bool screeningVary_ = false;
  double screeningSigma_ = 0.0;   ///< eV, prior width (0 = none)
  double cmFactor_ = 1.0;        ///< c.m. / lab energy of the entrance pair, set by ConvertSigmaToCM

  // Energy straggling: sigma = c_s sqrt(dE_cm) (c.m. coefficient), or Bohr's formula when c_s = 0.
  bool isStraggling_;
  double stragglingCoefficient_;
  int bohrZ1_ = 0;               ///< beam charge for Bohr straggling (SetBohrCharges)
  int bohrZ2_ = 0;               ///< target charge for Bohr straggling
  bool bohrChargesSet_ = false;

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
