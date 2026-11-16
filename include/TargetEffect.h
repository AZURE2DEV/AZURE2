#ifndef TARGETEFFECT_H
#define TARGETEFFECT_H

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <iosfwd>
#include "Constants.h"
#include "Equation.h"
#include "TargetStopping.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/SetupSpecs.h"
#include "azure2/inference/GaussianPrior.h"

/// An AZURE target effect entry

/*!
 * Experimential effects including gaussian beam convolution, target
 * integration, and a combination of the two are grouped under the TargetEffect
 * class. An object is created corresponding to each corresponding entry in
 * AZURESetup2.
 */

class TargetEffect {
 public:
  /// Build one effect from its spec (the data sets it applies to: GetSegmentsList).
  TargetEffect(const azure2::experiment::TargetEffectSpec &, const Config &);
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
   * (WorkingCopy::FillMnParams) when vary is 1, with a Gaussian prior of width sigma about the value read when
   * sigma > 0.  The line need not carry convolution or target integration.
   */
  bool HasScreening() const { return hasScreening_; }
  /*!
   * Its place among the active target effects (1-based; 0 until read).  A
   * line listing several segments gives each segment its own copy of the
   * effect (WorkingCopy::ReadTargetEffectsFile), so effect numbers and lines differ;
   * a screening parameter belongs to the line, shared by its copies.
   */
  int GetLineNumber() const { return lineNumber_; }
  void SetLineNumber(int n) { lineNumber_ = n; }
  /// A `screening` keyword without its three numbers (or a negative sigma).
  bool HasScreeningError() const { return screeningError_; }
  /// U_e in eV as the line gives it, the start of a fit (the present value is a working copy's, TargetEffectValues).
  double GetScreeningUe() const { return screeningUe_; }
  /// The centre of the U_e prior, eV: as written, else the value read (docs/upgrade/parameter_kinds.md section 4).
  double GetScreeningNominal() const { return screeningNominal_; }
  bool IsScreeningVaried() const { return screeningVary_; }
  /// The prior covariance replaces the screening, density or Doppler prior (docs/upgrade/objective.md D3).
  void ReplaceScreeningPrior() { screeningPriorReplaced_ = true; }
  bool IsScreeningPriorReplaced() const { return screeningPriorReplaced_; }
  void ReplaceDensityPrior() { densityPriorReplaced_ = true; }
  bool IsDensityPriorReplaced() const { return densityPriorReplaced_; }
  void ReplaceDopplerPrior() { dopplerPriorReplaced_ = true; }
  bool IsDopplerPriorReplaced() const { return dopplerPriorReplaced_; }
  /// Prior width in eV, 0 = no prior.
  double GetScreeningSigma() const { return screeningSigma_; }
  /// ((U_e - nominal) / sigma)^2 at the present U_e (eV, a working copy's TargetEffectValues), or 0 without a prior or
  /// when U_e is fixed (a fixed parameter pays no prior, docs/upgrade/objective.md D5).
  double ScreeningPenalty(double ue) const {
    if (!hasScreening_ || !screeningVary_ || screeningPriorReplaced_ || screeningSigma_ <= 0.0) return 0.0;
    return azure2::inference::GaussianPriorChi2(ue, screeningNominal_, screeningSigma_);
  }
  /*!
   * Fitted target density (`densityfit vary sigma` on a target-integration line, sigma in atoms/cm^2 or as a
   * percentage of the line's density, "2.3%"): the areal density is a fit parameter `density_<n>` (WorkingCopy::
   * FillMnParams), with a Gaussian prior ((N - N0)/sigma)^2 about the line's value when sigma > 0.  It moves only
   * the target thickness, dE with int_{E-dE}^{E} dE'/epsilon = N, i.e. the depth window and the straggling with
   * depth: the yield is the integral over the window divided by the same N, so the depth-averaged cross section
   * carries no factor of N.  `densitygroup <label>` shares one parameter between lines (one physical target seen
   * by several detectors); the group's first line holds the parameter and its densityfit block.
   *
   * The sub-point grids are laid out once for the largest density the fit is expected to reach
   * (GetLayoutDensity: N0 + 3 sigma, at least 1.1 N0; 1.5 N0 without a prior; N0 the layout base) and the depth integral is cut at
   * E0 - dE(N), so the yield moves smoothly with N and the grids stay put (Unit::IntegrateTargetEffect).
   */
  bool HasDensityFit() const { return hasDensityFit_; }
  /// A malformed `densityfit` / `densitygroup` block, or one on a line that does not integrate a thin target.
  bool HasDensityFitError() const { return densityFitError_; }
  const std::string &GetDensityFitErrorText() const { return densityFitErrorText_; }
  bool IsDensityVaried() const { return densityVary_; }
  /// What the layout density and the lower limit N/100 are laid out from (atoms/cm^2): the prior's centre, or the
  /// line's density (azure2/experiment/LayoutBase.h, parameter_kinds.md section 4).
  double GetDensityLayoutBase() const { return densityLayoutBase_; }
  /// The centre of the density prior (atoms/cm^2): as written, else the line's density (parameter_kinds.md section 4).
  double GetDensityCentre() const { return densityCentre_; }
  /// The name of the target record, which names the effect's parameters ("" for an effect without a target).
  const std::string &GetTargetName() const { return targetName_; }
  /// Prior width in atoms/cm^2 (a percentage already converted), 0 = no prior.
  double GetDensitySigma() const { return densitySigma_; }
  /// The `densitygroup` label, empty without one.
  const std::string &GetDensityGroup() const { return densityGroup_; }
  /// ((N - N0)/sigma)^2 at the present density N (atoms/cm^2), or 0 without a prior or when N is fixed (objective.md D5).
  double DensityPenalty(double density) const {
    if (!hasDensityFit_ || !densityVary_ || densityPriorReplaced_ || densitySigma_ <= 0.0) return 0.0;
    return azure2::inference::GaussianPriorChi2(density, densityCentre_, densitySigma_);
  }
  /// The density the sub-point grids span (the line's density when it is not fitted).
  double GetLayoutDensity() const { return (layoutDensity_ > 0.0) ? layoutDensity_ : density_; }
  void SetLayoutDensity(double n) { layoutDensity_ = n; }
  /// Join a density group: take its first line's density, prior and layout density.
  void JoinDensityGroup(const TargetEffect &leader);
  /// The c.m. beam-spread sigma at c.m. energy \p eCM for a point whose entrance pair has the c.m./lab
  /// factor \p cmFactor: the one place every kernel takes its width from (KernelParametersFor).  A constant sigma (lab, as read)
  /// is scaled by the factor; a lab expression is evaluated at eCM / F and its value scaled by F; a legacy
  /// expression is used as it stands.
  double SigmaCM(double eCM, const Config &, double cmFactor) const;
  /*!
   * The kernels of this effect (azure2/experiment/Kernels.h): the free-gas Doppler broadening when the line has one,
   * then the resolution kernel chosen when the effect was built (target integration, else a Gaussian beam spread,
   * else UDR, else a beam profile, or the table kernel the spec names).  Straggling (a c.m. coefficient c_s in
   * keV^1/2, or Bohr's formula when it is on with c_s = 0) and the beam spread are parameters of the target kernel.
   */
  const azure2::experiment::EffectKernels &Kernels() const { return kernels_; }
  /// The parameters the kernels read, with the formulas (sigma, stopping) evaluated through this effect and \p
  /// configure: valid while both live.
  azure2::experiment::KernelParameters KernelParametersFor(const Config &configure) const;
  /// Two resolution kernels on one effect (a target integration, a beam profile, a UDR or a beam spread outside a
  /// target, two of them): refused at set-up, as one of them would be ignored.
  bool HasKernelError() const { return !kernelErrorText_.empty(); }
  const std::string &GetKernelErrorText() const { return kernelErrorText_; }
  /// Target density in atoms/cm^2 as the line gives it (the present value is a working copy's, TargetEffectValues).
  /// Target integration only.
  double GetDensity() const;
  /// A target integration with a density of 0 is an infinitely THICK target: the yield per incident
  /// particle, int_0^{E0} sigma(E) / epsilon(E) dE, not divided by a density; every point integrates
  /// down to 1 keV.  (Before 2026-10-27 this was done by overwriting the effect's density with 1e24
  /// at the first zero-thickness point, which left that point with a thickness of 0 and a yield of
  /// exactly 0 under a beam spread or straggling, and also fired for a thin target whose stopping
  /// cross section came out <= 0; audit item 17.)
  bool IsThickTarget() const { return isThickTarget_; }
  /*!
   * A yield on a layer of finite thickness (format 3: quantity yield on a target of kind layer; docs/upgrade/
   * energy_grid.md section 5): reactions per incident particle, Y(E0) = 1e-24 int_{E_out}^{E0} <sigma>/epsilon dE
   * with int_{E_out}^{E0} dE/epsilon = N (the working copy's density, fitted or not), E_out found on the point's own
   * sub-points with the yield's own quadrature (LayerExit.h).  Not divided by N, unlike a thin target's
   * target-averaged cross section.  The grid reaches LayerGridMargin() beyond the layer at its layout density.
   */
  bool IsLayerYield() const { return isLayerYield_; }
  /// The factor on the layout thickness (TargetThickness) the sub-point grid of a layer yield spans, so that E_out of
  /// the largest density a fit may reach, found by the grid's own quadrature, lies inside the grid: the two integrals
  /// of 1/epsilon differ by the quadrature's error, measured up to 8e-4 relative in tests/layer_yield (10 nodes across
  /// a layer of half the beam energy, epsilon ~ 1/E).  A grid that still ends inside the layer gives NaN.
  static constexpr double LayerGridMargin() { return 1.02; }
  /// Energy loss (lab MeV) across the target at lab energy E: Delta with int_{E-Delta}^{E} dE'/(s epsilon) = N, the
  /// exact range integral, at the line's density and stopping scale.
  double TargetThickness(double, const Config &) const;
  /// The same for an areal density \p density (atoms/cm^2) and a stopping scale \p scale instead of the line's.
  double TargetThickness(double, const Config &, double density, double scale) const;
  /// The c.m. energy loss (MeV) below a point at lab energy \p eLab, c.m. energy \p eCM (factor \p cmFactor) that the
  /// nodes of its target integration span: an infinitely thick target down to 1 keV, otherwise the thickness at the
  /// layout density and stopping scale, a layer's widened by LayerGridMargin().  The set-up lays the support out with
  /// it at the nominal energy, an energy shift moves it to the true one (ShiftedFold.h): one formula for both.
  double SupportThickness(double eLab, double eCM, double cmFactor, const Config &configure) const;
  /// The c.m. energy loss (MeV) the target kernel reads at a point (KernelPoint::targetThickness) at areal density
  /// \p density and stopping scale \p scale: an infinitely thick target down to 1 keV, otherwise the target thickness.
  double KernelThickness(double eLab, double eCM, double cmFactor, const Config &configure, double density,
                         double scale) const;
  /// Attenuation coefficient of the given order.
  double GetQCoefficient(int) const;
  /// Convolution coefficient of the given order.
  double GetConvCoefficient(int) const;
  void SetSigma(double);
  void SetDensity(double);
  void SetNumSubPoints(int);
  /// The keys of the data sets the effect applies to.
  std::vector<int> GetSegmentsList() const;
  /*!
   * The stopping the effect's target integration uses, lab MeV cm^2 per active nucleus, at lab energy \p eLab (MeV) and
   * stopping scale \p scale: s x epsilon_active(E), epsilon from the target's source (formula | evaluation | table,
   * TargetStopping).  The one function every caller asks; the scale is applied here and nowhere else.
   */
  double StoppingAt(double eLab, const Config &configure, double scale) const {
    return scale * stopping_.PerActiveNucleus(eLab, configure);
  }
  const TargetStopping &GetStopping() const { return stopping_; }
  /*!
   * The stopping-scale nuisance s (energy_grid.md section 6.2, decision D6): one parameter per (material, projectile),
   * shared by every effect with the same GetStoppingScaleKey(); made by WorkingCopy::SetUpStoppingScaleParams.  An
   * evaluation or a table always has one (fixed at 1 unless the record says otherwise), a formula only when its
   * record gives `scale`.  Its prior ((s - s0)/sigma)^2 counts when it is fitted; a fitted evaluation's scale without a
   * written prior takes the library's stated uncertainty for the projectile class (TargetStopping::
   * StatedRelUncertainty).
   */
  bool HasStoppingScale() const { return hasStoppingScale_; }
  /// The scale as written, the start value and the prior centre (1 by default).
  double GetStoppingScaleNominal() const { return stoppingScaleNominal_; }
  bool IsStoppingScaleVaried() const { return stoppingScaleVary_; }
  /// Prior width (absolute; s is dimensionless), 0 = none.
  double GetStoppingScaleSigma() const { return stoppingScaleSigma_; }
  /// The centre of the stopping scale's prior: as written, else the start scale (parameter_kinds.md section 4).
  double GetStoppingScaleCentre() const { return stoppingScaleCentre_; }
  /// True when the width is the library's default, not written (it is printed at set-up).
  bool IsStoppingScaleSigmaDefault() const { return stoppingScaleSigmaDefault_; }
  const std::string &GetStoppingScaleKey() const { return stoppingScaleKey_; }
  const std::string &GetStoppingScaleName() const { return stoppingScaleName_; }
  /// The labels the name is made of: the material ("TiN", or the target's name for a formula or table), the projectile.
  const std::string &GetStoppingScaleMaterial() const { return stoppingScaleMaterial_; }
  const std::string &GetStoppingScaleProjectile() const { return stoppingScaleProjectile_; }
  /// The largest scale the sub-point grids span (as the layout density): s0 + 3 sigma, at least 1.1 s0, when fitted
  /// with a prior; 1.5 s0 fitted without; s0 when fixed.
  double GetLayoutStoppingScale() const;
  /// What the layout scale and the lower limit are laid out from: the prior's centre, or the scale written
  /// (azure2/experiment/LayoutBase.h, parameter_kinds.md section 4).
  double GetStoppingScaleLayoutBase() const { return stoppingScaleLayoutBase_; }
  /// The smallest scale a fit may take: s0 / 10, s0 the layout base.
  double GetStoppingScaleMinimum() const { return 0.1 * stoppingScaleLayoutBase_; }
  /// A table kernel's parameter values as written (engine units, KernelParameters::values), and those a fit may vary
  /// (ApparatusParameters.h, row `kernel`).
  const std::vector<double> &GetKernelValues() const { return kernelValues_; }
  const std::vector<azure2::experiment::ApparatusParameterSpec> &GetKernelParameterFits() const { return kernelFits_; }
  /// The value a fitted kernel parameter's grid spans (its upper limit) and what that and its lower limit are laid out
  /// from (the prior's centre or the value, azure2/experiment/LayoutBase.h), by slot.
  double GetKernelLayoutValue(int slot) const { return kernelLayoutValues_.at(slot); }
  double GetKernelLayoutBase(int slot) const { return kernelLayoutBases_.at(slot); }
  Equation *GetConvolutionEq();

  // Energy straggling (a parameter of the target kernel: azure2::experiment::StragglingSigma, BohrStragglingSigma)
  bool IsStraggling() const;
  double GetStragglingCoefficient() const;

  // Adaptive integration grid parameters
  double GetResonanceWidthMultiplier() const;
  double GetPointsPerWidth() const;

  /// The multiple of sigma above and below centroid energy to use as integration range
  /// (5 since dev 7d24c88: at 3 the 0.27 % of kernel weight beyond the window was dropped and the rest rescaled).
  static constexpr double convolutionRange = azure2::experiment::kConvolutionRange;
  /*Window for THM (HOES) segments: their extremely narrow resonances make
  the smeared curve in the valleys a pure Gaussian tail of a neighboring spike,
  so a +-3 sigma cutoff loses up to ~25% there; +-5 sigma keeps the truncation
  error below ~1e-3.  Equal to convolutionRange since dev 7d24c88; kept as its
  own name so the THM window is stated where it is used.*/
  static constexpr double thmConvolutionRange = azure2::experiment::kThmConvolutionRange;
  /*Numerical floor (in MeV) kept between a target-effect sub-point energy and the
  channel threshold, since Coulomb functions and penetrabilities are singular at
  E=0. For ordinary segments this bounds the channel CM energy directly; for THM
  (HOES) segments, which legitimately probe negative channel CM energies, it
  instead bounds the compound-system energy (CMEnergy+SepE+ExE) so the floor is
  applied relative to the physical threshold rather than an arbitrary channel origin.*/
  static constexpr double minIntegrationEnergy = 0.001;

  /// One skewed-Gaussian component of a beam energy profile.
  using BeamProfileComponent = azure2::experiment::BeamProfileComponent;
  /// Beam-profile kernel: an absolute (not point-centred) beam energy profile,
  /// a Gaussian detector-resolution window per point and, optionally, the
  /// detailed-balance weight of an inverse photodissociation measurement.
  /// Written as the trailing block `beamprofile N {xi omega alpha w}xN
  /// tpcSigma nCut dbFlag` of the old line (a beam's profile field in format 3), energies in the lab frame.
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
   * (TransformedValue in src/ModelPass.cpp) and the line's resolution kernel (Gaussian, UDR, target integration, beam
   * profile) applied: SAMMY's order.  A line with only `doppler` evaluates sigma_D at the point.
   */
  bool IsDoppler() const { return isDoppler_; }
  /// A malformed `doppler` block (missing or non-positive temperature, or a malformed fit form).
  bool HasDopplerError() const { return dopplerError_; }
  /// Why the `doppler` block was refused (HasDopplerError).
  const std::string &GetDopplerErrorText() const { return dopplerErrorText_; }
  /// Effective temperature, K, as the line gives it (the present value is a working copy's, TargetEffectValues).
  double GetDopplerTemperature() const { return dopplerTemperature_; }
  /*!
   * Fitted temperature (`doppler T vary sigma [Tmax]`): T is the fit parameter `doppler_<n>` (n = the effect
   * number, as screening_<n>), with a Gaussian prior ((T - T0)/sigma)^2 about the line's T0 when sigma > 0.  The
   * sub-point grids and the integration windows are laid out once for the layout temperature
   * (GetDopplerLayoutTemperature: Tmax, else T0 + 3 sigma but at least 1.21 T0, or 2 T0 without a prior) and the
   * Gauss-Legendre pieces for T0 (Unit::DopplerWeights), so the broadened cross section is smooth in T and
   * neither moves during a fit.
   */
  bool HasDopplerFit() const { return dopplerFit_; }
  bool IsDopplerVaried() const { return dopplerVary_; }
  /// The temperature the pieces, the layout temperature and the lower limit are laid out for (K): a fitted line's
  /// prior centre or its temperature (azure2/experiment/LayoutBase.h), else the line's temperature.
  double GetDopplerLayoutBase() const { return dopplerFit_ ? dopplerLayoutBase_ : dopplerTemperature_; }
  /// Prior width in K, 0 = no prior.
  /// The Tmax field of "doppler T vary sigma Tmax" as written (K), 0 without one.
  double GetDopplerTmax() const { return dopplerTmax_; }
  double GetDopplerSigma() const { return dopplerSigma_; }
  /// The centre of the temperature prior (K): as written, else the line's temperature (parameter_kinds.md section 4).
  double GetDopplerCentre() const { return dopplerCentre_; }
  /// The temperature the sub-point grids and the integration windows span (the line's T when not fitted).
  double GetDopplerLayoutTemperature() const { return (dopplerLayout_ > 0.0) ? dopplerLayout_ : dopplerTemperature_; }
  /// The lowest temperature a fit may reach: T0/4 (the pieces, u(T0)/4 wide, stay within u(T)/2).
  double GetDopplerMinimum() const { return 0.25 * this->GetDopplerLayoutBase(); }
  /// ((T - T0)/sigma)^2 at the present temperature T (K), or 0 without a prior or when T is fixed (objective.md D5).
  double DopplerPenalty(double temperature) const {
    if (!dopplerFit_ || !dopplerVary_ || dopplerPriorReplaced_ || dopplerSigma_ <= 0.0) return 0.0;
    return azure2::inference::GaussianPriorChi2(temperature, dopplerCentre_, dopplerSigma_);
  }
  /// Does the line also carry a kernel applied after the Doppler step (Gaussian, UDR, target integration,
  /// beam profile)?  Without one the value is sigma_D at the point.
  bool HasResolutionStage() const { return kernels_.resolution != nullptr; }
  /// Sigma of the Gaussian detector (TPC) energy resolution.
  double GetBeamTpcSigma() const;
  /// Truncate each component outside mean +- this many standard deviations; 0 = none.
  double GetBeamTruncation() const;
  /// Weight the integrand by the detailed-balance factor of the inverse reaction.
  bool IsBeamPhotodissociation() const;
  // Every energy the line carries (the Gaussian sigma, a `sigmalab` expression, the beam-profile
  // components and TPC sigma, the stopping cross section) is kept in the LAB frame of the input file and
  // converted at the point: each kernel takes the point's own c.m./lab factor (Unit::SetTargetKinematics).
  // Up to 2026-10-27 the sigma and the profile were converted once, with the factor of the first segment
  // on the line, so a line shared by entrance pairs of different mass ratio (natural-target components,
  // segments of several pairs) used the wrong factor for the others (audit 2026-10-25, tier 3 entry 3).
  const std::vector<BeamProfileComponent> &GetBeamProfile() const;

  /// User-defined numerical resolution function (SAMMY's UDR): a file of
  /// time-of-flight delay distributions tabulated at a set of energies, in
  /// SAMMY's UDR format, convolved with a Gaussian burst and a rectangular
  /// channel width.  Written as the trailing block
  /// `udr "<file>" flightPath_m burstFWHM_ns channel_ns centred` of a
  /// the old line (a detector's udr field in format 3); a negative channel value is |n| bins per decade of
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
  /// convolved with the burst and channel width, normalised.  Null, with
  /// `refusal` saying why, when the burst and channel span more steps of the
  /// table than kMaxUdrHalfNodes on each side.
  std::shared_ptr<const azure2::experiment::UdrKernel> BuildUdrKernel(double eLab, double mass, std::string &refusal) const;
  /// The largest number of table steps the burst-and-channel response may span on each side of zero (2^24: the
  /// convolution then holds three arrays of 2^25 numbers, 0.8 GB).  More means a table far finer than the resolution
  /// (or a flight time far beyond any measurement); refused rather than allocated.
  static constexpr double kMaxUdrHalfNodes = 16777216.0;

 private:
  void SelectKernels(const azure2::experiment::TargetEffectSpec &, const Config &);
  std::string targetName_;
  azure2::experiment::EffectKernels kernels_;  ///< Kernels
  std::vector<double> kernelValues_;           ///< a table kernel's parameter values (KernelParameters::values)
  std::vector<double> kernelLayoutValues_;     ///< the values its grid is laid out for (KernelParameters::layoutValues)
  std::vector<double> kernelLayoutBases_;      ///< what those and the lower limits are laid out from (GetKernelLayoutBase)
  std::vector<azure2::experiment::ApparatusParameterSpec> kernelFits_;  ///< GetKernelParameterFits
  std::string kernelErrorText_;                ///< HasKernelError
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
  double dopplerLayoutBase_ = 0.0;   ///< K, GetDopplerLayoutBase
  double dopplerCentre_ = 0.0;       ///< K, the prior's centre: as written, else the value read
  double dopplerSigma_ = 0.0;        ///< K, prior width (0 = none)
  double dopplerTmax_ = 0.0;         ///< K, the Tmax field of the fit form as written; 0 = none
  double dopplerLayout_ = 0.0;       ///< K, the temperature the grids span; 0 = the current temperature
  bool isUdr_ = false;
  std::string udrFile_;
  std::string udrOpenPath_;  // udrFile_ to open (Config::OpenPath); messages print udrFile_
  double udrFlightPath_ = 0.0;   ///< metres
  double udrBurstFwhm_ = 0.0;    ///< nanoseconds
  double udrChannelWidth_ = 0.0; ///< nanoseconds
  bool udrCentred_ = true;
  std::shared_ptr<UdrTable> udrTable_;  ///< shared between copies of the effect

  bool isConvolution_ = false;
  bool isTargetIntegration_ = false;
  bool isActive_ = false;  ///< a line ending inside its first nine fields is read as inactive
  bool isQCoefficients_ = false;
  bool isConvCoefficients_ = false;
  int numIntegrationPoints_ = 0;
  double sigma_ = 0.0;
  double density_ = 0.0;
  bool isThickTarget_ = false;   ///< target integration with density 0 (IsThickTarget)
  bool isLayerYield_ = false;    ///< yield per incident particle on a layer (IsLayerYield)
  TargetStopping stopping_;  ///< formula | evaluation | table (StoppingAt)
  bool hasStoppingScale_ = false;
  double stoppingScaleNominal_ = 1.0;
  double stoppingScaleLayoutBase_ = 1.0;  ///< GetStoppingScaleLayoutBase
  bool stoppingScaleVary_ = false;
  double stoppingScaleSigma_ = 0.0;
  double stoppingScaleCentre_ = 1.0;  ///< the prior's centre: as written, else the start scale
  bool stoppingScaleSigmaDefault_ = false;
  std::string stoppingScaleKey_, stoppingScaleName_, stoppingScaleMaterial_, stoppingScaleProjectile_;
  std::vector<int> datasetKeys_;  ///< the data sets it applies to (DatasetSpec::key)
  vector_r qCoefficients_;
  vector_r convCoefficients_;
  Equation convolutionEq_;
  bool sigmaLab_ = false;        ///< the expression is lab in, lab out (`sigmalab` keyword)
  bool hasScreening_ = false;     ///< `screening` keyword present
  bool screeningError_ = false;   ///< malformed `screening` block
  int lineNumber_ = 0;
  double screeningUe_ = 0.0;      ///< eV
  double screeningNominal_ = 0.0; ///< eV, the prior's centre: as written, else the value read
  bool screeningVary_ = false;
  bool screeningPriorReplaced_ = false, densityPriorReplaced_ = false, dopplerPriorReplaced_ = false;
  double screeningSigma_ = 0.0;   ///< eV, prior width (0 = none)
  bool hasDensityFit_ = false;    ///< `densityfit` block present (or joined through a density group)
  bool densityFitError_ = false;
  std::string densityFitErrorText_;
  bool densityVary_ = false;
  double densityLayoutBase_ = 0.0;  ///< atoms/cm^2, GetDensityLayoutBase
  double densityCentre_ = 0.0;    ///< atoms/cm^2, the prior's centre: as written, else the value read
  double densitySigma_ = 0.0;     ///< atoms/cm^2, prior width (0 = none)
  std::string densityGroup_;
  double layoutDensity_ = 0.0;    ///< density the sub-point grids span; 0 = the current density

  // Energy straggling: sigma = c_s sqrt(dE_cm) (c.m. coefficient), or Bohr's formula when c_s = 0.
  bool isStraggling_ = false;
  double stragglingCoefficient_ = 0.0;

  // Adaptive integration grid parameters
  double resonanceWidthMultiplier_ = 0.0;  ///< set from AdaptiveIntegrationGrid::GridConfig in the constructor
  bool explicitResonanceWidthMultiplier_ = false;
  double pointsPerWidth_ = 0.0;            ///< likewise

  // Beam-profile kernel (see IsBeamProfile).
  bool isBeamProfile_ = false;
  std::vector<BeamProfileComponent> beamProfile_;
  double beamTpcSigma_ = 0.0;
  double beamTruncation_ = 0.0;
  bool beamPhotodissociation_ = false;
};

#endif
