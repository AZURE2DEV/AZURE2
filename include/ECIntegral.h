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
    /// The finite-wavelength (Friar-Fallieros) charge function multiplies the
    /// integrand: electric terms without the long-wavelength approximation.
    /// Magnetic terms keep their long-wavelength operator in either mode.
    bool useChargeFunction;
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
  ECIntegral(PPair *pPair, double channelRadius, const Config &configure) {
    chanrad_ = channelRadius;
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
    params_.useChargeFunction = false;  // set per term in CalculateDirect
    pair_ = pPair;
    configure_ = &configure;
  };
  /*!
   * The CoulFunc and WhitFunc objects are automatically destroyed with the object.
   */
  ~ECIntegral() = default;
  /// External-capture radial integral for one pathway: entrance and exit l, the multipolarity and radiation type, and the energies.
  complex operator()(int, int, double, double, double, double, int, char, double, double, bool, bool = false);
  /// Relative tolerance of the adaptive radial quadrature (default 1e-4); ECTable
  /// computes its few nodes much tighter than a direct integral at every energy.
  void SetRelativeTolerance(double epsrel) { epsrel_ = epsrel; }
  /// Nonzero if a quadrature of the last call reported an error (GSL status),
  /// which happens only with the GSL error handler off.
  int Status() const { return status_; }

 private:
  /// Clear the cached FW, GW and WW integrals.
  void ResetIntegrals() {
    FW_ = 0.;
    GW_ = 0.;
    FWc_ = complex(0.0, 0.0);
    GWc_ = complex(0.0, 0.0);
  };
  /// Evaluate the three radial integrals at one energy.
  void Integrate(double);
  /// GSL integrand: regular Coulomb function against the Whittaker function.
  static double FWIntegrand(double, void *);
  /// GSL integrand: irregular Coulomb function against the Whittaker function.
  static double GWIntegrand(double, void *);
  /// GSL integrand: Whittaker against Whittaker, for a bound entrance channel.
  static double WWIntegrand(double, void *);
  /// Clears the per-x FW/GW memo (see ECIntegral.cpp); called before the
  /// FW/GW pair so the two integrands share their expensive evaluations.
  static void ClearFGMemo();
  /// Coulomb functions of this pair.
  CoulFunc *coulfunction() const { return params_.coulFunc.get(); };
  /// Whittaker functions of this pair.
  WhitFunc *whitfunction() const { return params_.whitFunc.get(); };
  /// The pair this integral belongs to.
  PPair *pair() const { return pair_; };
  /// The channel radius the integral starts at (fm): the transition's channels' radius, which since
  /// 2026-11-10 can differ from the pair's (<channelRadii>).
  double chanRad() const { return chanrad_; };
  double FW() const { return FW_; };
  double GW() const { return GW_; };
  /*!
   * External-capture radial integrals into an UNBOUND final state (E_f > 0 above the pair's threshold),
   * treated as a narrow resonance: the final state's exterior is its outgoing wave, normalized at the channel
   * radius, u_f(r) = O_lf(k_f r) / O_lf(k_f a) with O = G + iF (the resonance ANC is the amplitude of that
   * wave; Mukhamedzhanov et al., PRC 96, 024623 (2017), eq. 16).  Sets FWc_ = int F_li u_f w dr and
   * GWc_ = int G_li u_f w dr from a to infinity, w = r^L (times the charge function).  The products of two
   * continuum functions do not decay; the integral from a to R is done by Gauss-Legendre panels and the one
   * beyond R by Abel-regularized integration by parts (three terms of the asymptotic series), which is what a
   * convergence factor exp(-eps r), eps -> 0, gives.  Returns false if R would have to exceed 50000 fm.
   */
  bool IntegrateUnbound(double chanrad, double finalEnergy);
  Params params_;
  PPair *pair_;
  double chanrad_;
  const Config *configure_;
  double FW_;
  double GW_;
  /// The FW/GW integrals of an unbound final state (complex: its exterior is an outgoing wave).
  complex FWc_ = complex(0.0, 0.0);
  complex GWc_ = complex(0.0, 0.0);
  /// Coulomb functions for the unbound-final-state integrals, without the global cache: their many radial
  /// nodes at each energy would only grow it.
  std::unique_ptr<CoulFunc> unboundCoul_;
  double epsrel_ = 1e-4;
  /// Above this Sommerfeld parameter an open channel is treated through its threshold limit (CalculateDirect).
  static constexpr double kMaxOpenChannelEta = 100.0;
  int status_ = 0;

  /*!
   * Calculate integral directly
   */
  complex CalculateDirect(int liValue, int lfValue, double siValue, double sfValue,
                          double jInitial, double jFinal, int multL, char radType,
                          double energy, double bindingEnergy, bool isChannelCapture);
};


#endif
