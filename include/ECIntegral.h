#ifndef ECINTEGRAL_H
#define ECINTEGRAL_H

#include "CoulFunc.h"
#include "WhitFunc.h"
#include "Config.h"
#include <memory>

class EffectiveCharge;

/// A function class to calculate external capture integrals

/*!
 * The ECIntegral function class calculates external capture integrals for
 * both positive and negative energy channels.  The results are returned as
 * an ECIntResult structure.
 */

class ECIntegral {
 public:
  /// Public only so the file-local FW/GW memo helper in ECIntegral.cpp can name
  /// the type; the params_ instance itself stays private.
  typedef struct Params {
    EffectiveCharge *effectiveCharge;
    std::unique_ptr<CoulFunc> coulFunc;
    std::unique_ptr<WhitFunc> whitFunc;
    int liValue;
    int lfValue;
    int multLValue;
    double pairEnergy;
    double bindingEnergy;
    bool useLongWavelengthApprox;
    /// Reference decades subtracted from the Whittaker-Whittaker integrand, so
    /// that a sub-threshold channel (where each factor is ~1e-500) integrates
    /// as O(1) numbers instead of underflowing to zero.  See WhitFunc::Scaled.
    double log10RefIn;
    double log10RefOut;
    bool useScaledWhittaker;
  } Params;
  /*!
   * The ECIntegral object is created with reference to a PPair object.
   * The PPair object is used to create new instances of the CoulFunc and
   * WhitFunc objects.
   */
  ECIntegral(PPair *pPair, const Config &configure) {
    params_.coulFunc = std::make_unique<CoulFunc>(pPair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
    // Disable the global Coulomb cache ONLY for the reaction rate.  There the
    // energies never recur, so the radial-integral lookups (keyed on radius,
    // the integration variable) only miss and grow the cache without bound.
    // In normal fit/extrapolation calculations the energies recur at fixed
    // channel radii and the cache is a large speed-up, so it stays enabled.
    if (configure.paramMask & Config::CALCULATE_REACTION_RATE)
      params_.coulFunc->SetUseGlobalCache(false);
    params_.whitFunc = std::make_unique<WhitFunc>(pPair);
    params_.useLongWavelengthApprox = !!(configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX);
    pair_ = pPair;
    configure_ = &configure;
  };
  /*!
   * The CoulFunc and WhitFunc objects are automatically destroyed with the object.
   */
  ~ECIntegral() = default;
  complex operator()(int, int, double, double, double, double, int, char, double, double, bool, bool = false);

 private:
  void ResetIntegrals() {
    FW_ = 0.;
    GW_ = 0.;
  };
  void Integrate(double);
  static double FWIntegrand(double, void *);
  static double GWIntegrand(double, void *);
  static double WWIntegrand(double, void *);
  /// Clears the per-x FW/GW memo (see ECIntegral.cpp); called before the
  /// FW/GW pair so the two integrands share their expensive evaluations.
  static void ClearFGMemo();
  CoulFunc *coulfunction() const { return params_.coulFunc.get(); };
  WhitFunc *whitfunction() const { return params_.whitFunc.get(); };
  PPair *pair() const { return pair_; };
  double FW() const { return FW_; };
  double GW() const { return GW_; };
  Params params_;
  PPair *pair_;
  const Config *configure_;
  double FW_;
  double GW_;

  /*!
   * Calculate integral directly
   */
  complex CalculateDirect(int liValue, int lfValue, double siValue, double sfValue,
                          double jInitial, double jFinal, int multL, char radType,
                          double energy, double bindingEnergy, bool isChannelCapture);
};


#endif
