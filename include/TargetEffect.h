#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
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
  /*!
   * Fitted target density (`densityfit vary sigma` on a target-integration line, sigma in atoms/cm^2 or as a
   * percentage of the line's density, "2.3%"): the areal density is a fit parameter `density_<n>` (EData::
   * FillMnParams), with a Gaussian prior ((N - N0)/sigma)^2 about the line's value when sigma > 0.  It moves only
   * the target thickness, dE with int_{E-dE}^{E} dE'/epsilon = N, i.e. the depth window and the straggling with
   * depth: the yield is the integral over the window divided by the same N, so the depth-averaged cross section
   * carries no factor of N.  `densitygroup <label>` shares one parameter between lines (one physical target seen
   * by several detectors); the group's first line holds the parameter and its densityfit block.
   *
   * The sub-point grids are laid out once for the largest density the fit is expected to reach
   * (GetLayoutDensity: N0 + 3 sigma, at least 1.1 N0; 1.5 N0 without a prior) and the depth integral is cut at
   * E0 - dE(N), so the yield moves smoothly with N and the grids stay put (EPoint::IntegrateTargetEffect).
   */
  bool HasDensityFit() const { return hasDensityFit_; }
  /// A malformed `densityfit` / `densitygroup` block, or one on a line that does not integrate a thin target.
  bool HasDensityFitError() const { return densityFitError_; }
  const std::string &GetDensityFitErrorText() const { return densityFitErrorText_; }
  bool IsDensityVaried() const { return densityVary_; }
  /// The line's density, the centre of the prior (atoms/cm^2).
  double GetDensityNominal() const { return densityNominal_; }
  /// Prior width in atoms/cm^2 (a percentage already converted), 0 = no prior.
  double GetDensitySigma() const { return densitySigma_; }
  /// The `densitygroup` label, empty without one.
  const std::string &GetDensityGroup() const { return densityGroup_; }
  /// ((N - N0)/sigma)^2, or 0 without a prior.
  double DensityPenalty() const {
    if (!hasDensityFit_ || densitySigma_ <= 0.0) return 0.0;
    const double r = (density_ - densityNominal_) / densitySigma_;
    return r * r;
  }
  /// The density the sub-point grids span (the line's density when it is not fitted).
  double GetLayoutDensity() const { return (layoutDensity_ > 0.0) ? layoutDensity_ : density_; }
  void SetLayoutDensity(double n) { layoutDensity_ = n; }
  /// Join a density group: take its first line's density, prior and layout density.
  void JoinDensityGroup(const TargetEffect &leader);
  /// The c.m. beam-spread sigma at c.m. energy \p eCM for a point whose entrance pair has the c.m./lab
  /// factor \p cmFactor: the one place every kernel takes its width from.  A constant sigma (lab, as read)
  /// is scaled by the factor; a lab expression is evaluated at eCM / F and its value scaled by F; a legacy
  /// expression is used as it stands.
  double SigmaCM(double eCM, const Config &, double cmFactor) const;
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
  /// The Bohr straggling sigma (c.m. MeV) after traversing \p arealDensity atoms/cm^2 of the target, for a
  /// point whose entrance pair has the c.m./lab factor \p cmFactor and charges \p z1 (beam), \p z2 (target).
  double BohrStragglingSigma(double arealDensity, double cmFactor, int z1, int z2) const;
  /// How far the sub-point grid of a point at c.m. energy \p eCM must reach above (\p forward) and below
  /// (\p backward) it for this effect; \p thicknessCM is the point's target thickness (c.m. MeV, 0 if none),
  /// \p cmFactor and \p z1, \p z2 the point's entrance-pair factor and charges (Bohr straggling).
  void GridExtent(double eCM, double thicknessCM, const Config &, double cmFactor, int z1, int z2, double &forward, double &backward) const;
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
  /// The same for an areal density \p density instead of the current one.
  double TargetThickness(double, const Config &, double density);
  /// Depth-averaged energy loss (lab MeV) in the target at lab energy E: E - E_bar with
  /// E_bar = int E' dE'/epsilon / int dE'/epsilon over [E - Delta, E] (Delta = TargetThickness).  The cross section
  /// at E_bar equals the target average for a cross section linear across the target; the remainder is about
  /// sigma'' Delta^2 / 24.  It is the single energy a thin-target point is evaluated at outside the windows of an
  /// energy-restricted target integration (EPoint fallback, 2026-11-05).  0 for a thick target.
  double MeanEnergyLoss(double, const Config &, double density);
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
  /// The same for a given window list: the line's own ranges together with the automatic windows of one
  /// segment (EData, IsAutoWindows).  An empty list is "no window": 1 without "auto", 0 with it (no level
  /// asked for a window, so every point takes the effect's fallback).
  double BlendWeight(double labEnergy, const std::vector<std::pair<double, double> > &windows) const;
  /// "auto" in the windows token (2026-11-06): each segment of a thin-target integration gets windows placed
  /// from the level scheme, around every resonance narrow enough for the sigma(E_bar) fallback to miss it by
  /// more than GetAutoWindowTolerance() at its peak (Delta^2/(3 Gamma^2) for a Lorentzian of total width
  /// Gamma under a target of thickness Delta), reaching GetAutoWindowMargin() x max(Gamma, Delta) beyond the
  /// points whose target holds the resonance.  "auto:<tol>" and "auto:<tol>:<margin>" set them; explicit
  /// ranges in the same token are added.  Older engines read no range from it and integrate everywhere.
  bool IsAutoWindows() const { return autoWindows_; }
  double GetAutoWindowTolerance() const { return autoWindowTolerance_; }
  /// The margin: as given ("auto:<tol>:<margin>"), else max(3, 1/(2 sqrt(tol))) -- a window edge at d from a Lorentzian
  /// leaves the E_bar fallback a relative tail error of about Delta^2/(4 d^2), so d >= Delta/(2 sqrt(tol)) keeps it within
  /// the tolerance (5 at the default 0.01; 13C+alpha ND2021r 2026-11-07: margin 3 left 120 points > 1 % just outside edges).
  double GetAutoWindowMargin() const {
    return autoWindowMargin_ > 0.0 ? autoWindowMargin_ : std::max(3.0, 0.5 / std::sqrt(autoWindowTolerance_));
  }
  /// Any window at all, explicit or automatic.
  bool HasWindows() const { return !ranges_.empty() || autoWindows_; }

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

  /*!
   * Free-gas Doppler broadening (SAMMY's FGM, manual section "Free-Gas Model of Doppler Broadening"):
   * `doppler <T>` on the line, T the effective temperature of the target atoms in kelvin.  The cross section
   * of every sub-point is first broadened,
   *   sigma_D(E) = 1/(u sqrt(pi) v^2) int_0^inf [exp(-(v-w)^2/u^2) - exp(-(v+w)^2/u^2)] w^2 sigma(w^2) dw,
   * v = sqrt(E), w = sqrt(E'), E the lab-equivalent energy (E_cm over the point's c.m./lab factor), u^2 = kT m/M
   * with m/M the beam/target mass ratio of the entrance pair, and only then are the transmission/yield transform
   * (ESegment::TransformedValue) and the line's resolution kernel (Gaussian, UDR, target integration, beam
   * profile) applied: SAMMY's order.  A line with only `doppler` evaluates sigma_D at the point.
   */
  bool IsDoppler() const { return isDoppler_; }
  /// A malformed `doppler` block (missing or non-positive temperature, or a malformed fit form).
  bool HasDopplerError() const { return dopplerError_; }
  /// Why the `doppler` block was refused (HasDopplerError).
  const std::string &GetDopplerErrorText() const { return dopplerErrorText_; }
  /// Effective temperature, K (the current value: a fitted temperature follows the parameter vector).
  double GetDopplerTemperature() const { return dopplerTemperature_; }
  /// Set the current temperature (EData::FillDopplerFromParams); the nominal and layout temperatures stay.
  void SetDopplerTemperature(double t) { dopplerTemperature_ = t; }
  /*!
   * Fitted temperature (`doppler T vary sigma [Tmax]`): T is the fit parameter `doppler_<n>` (n = the effect
   * number, as screening_<n>), with a Gaussian prior ((T - T0)/sigma)^2 about the line's T0 when sigma > 0.  The
   * sub-point grids and the integration windows are laid out once for the layout temperature
   * (GetDopplerLayoutTemperature: Tmax, else T0 + 3 sigma but at least 1.21 T0, or 2 T0 without a prior) and the
   * Gauss-Legendre pieces for T0 (EPoint::DopplerWeights), so the broadened cross section is smooth in T and
   * neither moves during a fit.
   */
  bool HasDopplerFit() const { return dopplerFit_; }
  bool IsDopplerVaried() const { return dopplerVary_; }
  /// The line's temperature, the centre of the prior and the temperature the pieces are laid out for (K).
  double GetDopplerNominal() const { return dopplerFit_ ? dopplerNominal_ : dopplerTemperature_; }
  /// Prior width in K, 0 = no prior.
  double GetDopplerSigma() const { return dopplerSigma_; }
  /// The temperature the sub-point grids and the integration windows span (the line's T when not fitted).
  double GetDopplerLayoutTemperature() const { return (dopplerLayout_ > 0.0) ? dopplerLayout_ : dopplerTemperature_; }
  /// The lowest temperature a fit may reach: T0/4 (the pieces, u(T0)/4 wide, stay within u(T)/2).
  double GetDopplerMinimum() const { return 0.25 * this->GetDopplerNominal(); }
  /// ((T - T0)/sigma)^2, or 0 without a prior.
  double DopplerPenalty() const {
    if (!dopplerFit_ || dopplerSigma_ <= 0.0) return 0.0;
    const double r = (dopplerTemperature_ - dopplerNominal_) / dopplerSigma_;
    return r * r;
  }
  /// u = sqrt(kT m/M) in sqrt(MeV) for a pair of c.m./lab factor cmFactor = M/(m+M), at the current T.
  double DopplerU(double cmFactor) const;
  /// The same at temperature \p temperature (K).
  double DopplerUAt(double cmFactor, double temperature) const;
  /// The c.m. energy window [lo, hi] beyond which the Doppler kernel at c.m. energy eCM is below
  /// exp(-convolutionRange^2 / 2), i.e. +-convolutionRange standard deviations in velocity (lo >= 0), at the
  /// layout temperature (GetDopplerLayoutTemperature): the grids are padded and the integrals cut with it.
  void DopplerWindow(double eCM, double cmFactor, double &lo, double &hi) const;
  /// Does the line also carry a kernel applied after the Doppler step (Gaussian, UDR, target integration,
  /// beam profile)?  Without one the value is sigma_D at the point.
  bool HasResolutionStage() const { return IsGaussian() || isTargetIntegration_ || isUdr_ || isBeamProfile_; }
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
  // Every energy the line carries (the Gaussian sigma, a `sigmalab` expression, the beam-profile
  // components and TPC sigma, the stopping cross section) is kept in the LAB frame of the input file and
  // converted at the point: each kernel takes the point's own c.m./lab factor (EPoint::SetTargetKinematics).
  // Up to 2026-10-27 the sigma and the profile were converted once, with the factor of the first segment
  // on the line, so a line shared by entrance pairs of different mass ratio (natural-target components,
  // segments of several pairs) used the wrong factor for the others (audit 2026-10-25, tier 3 entry 3).
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
  /// cached: read the file through InputFiles (Config::cacheInputFiles).
  bool LoadUdrTable(double eMinLab, double eMaxLab, std::ostream &log, bool cached = false);
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
  bool isDoppler_ = false;
  bool dopplerError_ = false;
  std::string dopplerErrorText_;
  double dopplerTemperature_ = 0.0;  ///< K, current
  bool dopplerFit_ = false;          ///< `doppler T vary sigma [Tmax]`
  bool dopplerVary_ = false;
  double dopplerNominal_ = 0.0;      ///< K, as read
  double dopplerSigma_ = 0.0;        ///< K, prior width (0 = none)
  double dopplerLayout_ = 0.0;       ///< K, the temperature the grids span; 0 = the current temperature
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
  bool hasDensityFit_ = false;    ///< `densityfit` block present (or joined through a density group)
  bool densityFitError_ = false;
  std::string densityFitErrorText_;
  bool densityVary_ = false;
  double densityNominal_ = 0.0;   ///< atoms/cm^2, as read
  double densitySigma_ = 0.0;     ///< atoms/cm^2, prior width (0 = none)
  std::string densityGroup_;
  double layoutDensity_ = 0.0;    ///< density the sub-point grids span; 0 = the current density

  // Energy straggling: sigma = c_s sqrt(dE_cm) (c.m. coefficient), or Bohr's formula when c_s = 0.
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
  bool autoWindows_ = false;            ///< "auto" in the windows token (IsAutoWindows)
  double autoWindowTolerance_ = 0.01;   ///< peak error of the fallback a resonance may have outside a window
  double autoWindowMargin_ = 0.0;       ///< window reach beyond the target, in units of max(Gamma, Delta); 0 = from the tolerance

  // Beam-profile kernel (see IsBeamProfile).
  bool isBeamProfile_ = false;
  std::vector<BeamProfileComponent> beamProfile_;
  double beamTpcSigma_ = 0.0;
  double beamTruncation_ = 0.0;
  bool beamPhotodissociation_ = false;
};

#endif
