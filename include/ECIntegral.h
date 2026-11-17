#ifndef ECINTEGRAL_H
#define ECINTEGRAL_H

#include "CoulFunc.h"
#include "WhitFunc.h"
#include "Config.h"
#include <gsl/gsl_integration.h>

#include <map>
#include <memory>
#include <vector>

class EffectiveCharge;


/// The reaction-rate mode's fixed radial grid of one final channel (ECIntegral.cpp, getRateGrid): its key and nodes.
struct ECRateGridKey {
  int z1, z2;
  double redmass;
  int lf;
  int twoJf;                       ///< the final channel's 2J key (PPair::SpinOrbitTwoJ)
  std::vector<double> potential;  ///< physics::PotentialKey of the pair's potential (the final state is W~ of it)
  double bindingEnergy;
  double chanrad;
  bool operator<(const ECRateGridKey &o) const {
    if (z1 != o.z1) return z1 < o.z1;
    if (z2 != o.z2) return z2 < o.z2;
    if (redmass != o.redmass) return redmass < o.redmass;
    if (lf != o.lf) return lf < o.lf;
    if (twoJf != o.twoJf) return twoJf < o.twoJf;
    if (potential != o.potential) return potential < o.potential;
    if (bindingEnergy != o.bindingEnergy) return bindingEnergy < o.bindingEnergy;
    return chanrad < o.chanrad;
  }
};
struct ECRateGridEntry {
  std::vector<double> x;        // radial nodes
  std::vector<double> W;        // effective weights (Laguerre weight * e^{+b(x-a)})
  std::vector<double> whitOut;  // scaled final Whittaker at each node (ref subtracted)
  double refOut;                // log10 reference decades of the final Whittaker
  bool ok;
};

/*!
 * What the external-capture integrals of one thread keep between calls: the GSL workspace, the entrance Coulomb
 * waves at the quadrature's radii (CoulRadialMemo), the final-state Whittaker values (WhitScaledMemo) and the rate
 * mode's radial grids.  Each value is keyed by everything it depends on, so a stored value is the value a
 * computation gives.  Owned by the run (RunScratch.h: one per thread of the run), never process-wide.
 */
struct ECScratch {
  ECScratch();
  ~ECScratch();
  ECScratch(const ECScratch &) = delete;
  ECScratch &operator=(const ECScratch &) = delete;
  gsl_integration_workspace *workspace;  ///< 1000 intervals
  CoulRadialMemo radial;
  WhitScaledMemo whittaker;
  std::map<ECRateGridKey, ECRateGridEntry> rateGrids;
};


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
    int twoJi;  ///< the entrance channel's 2J key (PPair::SpinOrbitTwoJ of J_i)
    int twoJf;  ///< the final channel's 2J key (PPair::SpinOrbitTwoJ of J_f)
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
    /// F, G and the scaled final W at the radii of one FW/GW pair, which share them (cleared per Integrate).
    struct FG {
      double F, G, whit;
    };
    std::map<double, FG> fgMemo;
  } Params;
  /*!
   * The ECIntegral object is created with reference to a PPair object.
   * The PPair object is used to create new instances of the CoulFunc and
   * WhitFunc objects.
   */
  /// \p radius: the radius a_c (fm) of the final-state channel, the lower limit of the integral.
  ECIntegral(const PPair *pPair, const Config &configure, double radius);
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
  };
  /// Evaluate the three radial integrals at one energy.
  void Integrate(double);
  /// GSL integrand: regular Coulomb function against the Whittaker function.
  static double FWIntegrand(double, void *);
  /// GSL integrand: irregular Coulomb function against the Whittaker function.
  static double GWIntegrand(double, void *);
  /// GSL integrand: Whittaker against Whittaker, for a bound entrance channel.
  static double WWIntegrand(double, void *);
  /// Coulomb functions of this pair.
  CoulFunc *coulfunction() const { return params_.coulFunc.get(); };
  /// Whittaker functions of this pair.
  WhitFunc *whitfunction() const { return params_.whitFunc.get(); };
  /// The pair this integral belongs to.
  const PPair *pair() const { return pair_; };
  double FW() const { return FW_; };
  double GW() const { return GW_; };
  Params params_;
  /// The run's scratch of this thread (ECScratch).
  ECScratch &scratch_;
  double radius_;  ///< the final-state channel's radius a_c (fm)
  const PPair *pair_;
  const Config *configure_;
  double FW_;
  double GW_;
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
