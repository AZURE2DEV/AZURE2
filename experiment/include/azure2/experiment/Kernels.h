// The apparatus kernels of a data set's effect (docs/upgrade/energy_grid.md section 3.3): target integration, Gaussian
// beam spread, beam profile, UDR, free-gas Doppler broadening and the others of Kernels(), each a plain function of
// the point's nodes, the point and the effect's parameters, listed in one table (src/Kernels.cpp).
#ifndef AZURE2_EXPERIMENT_KERNELS_H
#define AZURE2_EXPERIMENT_KERNELS_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/LinearForm.h"
#include "azure2/experiment/UdrKernel.h"

namespace azure2 {
namespace experiment {

/*!
 * The nodes a point's kernel folds: node j = 1 .. n in descending c.m. energy (MeV), energy[j - 1], and, for a
 * kernel that integrates over a target (Kernel::usesTarget), the stopping cross section at each node, stopping[j - 1]
 * (c.m. MeV cm^2 per active nucleus: the lab stopping times the point's c.m./lab factor; dE_cm/eps_cm = dE_lab/eps_lab).
 * The arrays belong to the caller.
 */
struct KernelNodes {
  int n = 0;
  const double *energy = nullptr;
  const double *stopping = nullptr;  ///< null for a kernel that does not use the target
  double Energy(int j) const { return energy[j - 1]; }
  double Stopping(int j) const { return stopping[j - 1]; }
};

/// What a kernel reads of its point (energies c.m. MeV unless stated).
struct KernelPoint {
  double energy = 0.0;          ///< E0, the point's c.m. energy (after an energy shift)
  double labEnergy = 0.0;       ///< its lab energy (MeV)
  double nominalLabEnergy = 0.0;  ///< its lab energy before any shift (MeV): a beam-profile window stays put
  double cmFactor = 1.0;        ///< c.m./lab energy factor of the entrance pair
  int z1 = 0, z2 = 0;           ///< beam and target charges of the entrance pair (Bohr straggling)
  double targetThickness = 0.0;  ///< energy loss across the target at the working copy's density (target kernels)
  double density = 0.0;          ///< the working copy's areal density, atoms/cm^2 (target kernels)
  /// Thin target outside or at the edge of its windows (decided by the engine's set-up and automatic tolerance):
  /// the weight of the integration against sigma(E_bar), E_bar = E0 - fallbackDepth, and whether only sigma(E_bar)
  /// is taken (fallbackOnly: the point's nodes are the fallback nodes; fallbackNow: the full grid, read at E_bar).
  double blendWeight = 1.0;
  bool fallbackOnly = false;
  bool fallbackNow = false;
  double fallbackDepth = 0.0;
  double binLow = 0.0, binHigh = 0.0;  ///< the point's reconstructed-energy window (beam profile); NaN without one
  double photoQ = 0.0;           ///< Q_gamma of the inverse photodissociation (MeV, beam profile with detailed balance)
  double photoMass = 0.0;        ///< compound-nucleus mass (MeV), 0 = no detailed-balance weight
  const UdrKernel *udr = nullptr;  ///< the point's UDR kernel (built at set-up from the effect's table)
  double udrMass = 0.0;          ///< the projectile's rest energy (MeV)
  double temperature = 0.0;      ///< the working copy's Doppler temperature (K)
};

/// One component of a beam energy profile: a skewed Gaussian in the lab frame (MeV).
struct BeamProfileComponent {
  double xi;      ///< location parameter
  double omega;   ///< scale parameter
  double alpha;   ///< shape (skewness) parameter
  double weight;  ///< relative weight (luminosity) of the component
};

/*!
 * A function of energy the effect's owner evaluates (a formula of the input): the Gaussian beam spread sigma (c.m.
 * MeV) at a c.m. energy for a pair of c.m./lab factor f, or the stopping cross section (lab MeV cm^2) at a lab energy
 * (f unused).  `owner` and `context` are the owner's (the engine passes its effect and its configuration).
 */
struct EnergyFunction {
  double (*at)(const void *owner, const void *context, double energy, double cmFactor) = nullptr;
  const void *owner = nullptr;
  const void *context = nullptr;
  double operator()(double energy, double cmFactor) const { return at(owner, context, energy, cmFactor); }
};

/// The parameters of an effect as its kernels read them: what the data set's target, beam and detector records
/// hold (energies lab, as read; each kernel converts them with the point's factor).
struct KernelParameters {
  // Gaussian beam spread (constant or a formula), also inside the target kernel.
  bool gaussian = false;
  EnergyFunction sigma;  ///< c.m. sigma at a c.m. energy
  // Target integration.
  bool thickTarget = false;  ///< density 0: the yield per incident particle down to 1 keV
  bool layerYield = false;   ///< the yield of a layer of finite thickness (reactions per incident particle)
  bool densityFit = false;   ///< the depth integral is cut at the back of the target E0 - T(N)
  bool straggling = false;
  double stragglingCoefficient = 0.0;  ///< c.m. keV^1/2; 0 with straggling on = Bohr straggling
  double layoutDensity = 0.0;          ///< the density the grids span (atoms/cm^2)
  EnergyFunction stopping;             ///< lab stopping cross section at a lab energy (the support of a thick target)
  bool BohrStraggling() const { return straggling && !(stragglingCoefficient > 0.0); }
  // Beam profile.
  const std::vector<BeamProfileComponent> *profile = nullptr;
  double tpcSigma = 0.0;    ///< lab sigma of the detector's energy resolution in the profile's window
  double truncation = 0.0;  ///< components cut at mean +- truncation sd; 0 = none
  bool detailedBalance = false;
  // UDR.
  double flightPath = 0.0;  ///< m
  // Free-gas Doppler broadening.
  double dopplerNominal = 0.0;  ///< K: the temperature the Gauss-Legendre pieces are laid out for
  double dopplerLayout = 0.0;   ///< K: the temperature the grid and the windows span
  /// The values of a table kernel's own parameters (Kernel::parameters, in that order, in their engine units).
  const double *values = nullptr;
  int numValues = 0;
};

/// The c.m. energy interval a point's nodes must span for a kernel (KernelSupport): what the set-up reads.
struct SupportPoint {
  double energy = 0.0;      ///< E0, c.m. MeV
  double labEnergy = 0.0;   ///< lab MeV
  double cmFactor = 1.0;
  int z1 = 0, z2 = 0;
  double targetThickness = 0.0;  ///< c.m. MeV at the layout density (the grid's target), 0 without a target
  double minimumEnergy = 0.0;    ///< the lowest node energy allowed (c.m. MeV)
  bool thm = false;              ///< a THM (HOES) point: Gaussian folds reach thmRange sigma, profile tails too
  double binLow = 0.0, binHigh = 0.0;  ///< reconstructed-energy window, NaN without one
  const UdrKernel *udr = nullptr;
  double udrMass = 0.0;
};

/// One parameter of a kernel: its format-3 name, its dimension (Units.h) and engine unit, whether a fit may vary it.
struct KernelParameter {
  std::string name;
  std::string dimension;
  std::string unit;
  bool fittable = false;
  std::string doc;
};

template <class V>
using FoldFunction = V (*)(const KernelParameters &, const KernelPoint &, const KernelNodes &, const V *values);
/// d(fold)/d(parameter k of the table row) for node values `values` (double only); null when the row has none.
using ParameterDerivative = double (*)(const KernelParameters &, const KernelPoint &, const KernelNodes &,
                                       const double *values, int parameter);
/// The c.m. interval [low, high] the point's nodes must span.
using SupportFunction = void (*)(const KernelParameters &, const SupportPoint &, double &high, double &low);
/// A broadening of the node values: out[r] = the broadened value at c.m. energy energies[r], r < m (one call per point).
template <class V>
using BroadenFunction = void (*)(const KernelParameters &, const KernelPoint &, const KernelNodes &, const V *values,
                                 const double *energies, int m, V *out);
/// Weights w (node, weight) of a broadening at c.m. energy e: broadened(e) = sum w.second values[w.first - 1], and,
/// if dwdT is given, their derivatives in the temperature on the same nodes.
using BroadeningWeightsFunction = void (*)(const KernelParameters &, const KernelPoint &, const KernelNodes &, double e,
                                           std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT);

/*!
 * One row of the table of kernels.  A resolution kernel folds the node values into the point's value (fold on
 * double, weights on LinearForm: the same template, so the gradient reads the forward arithmetic's own weights); a
 * broadening kernel (free-gas Doppler) maps the node values to broadened values at an energy before the resolution
 * kernel, SAMMY's order.
 */
struct Kernel {
  std::string name;
  enum class Stage { kBroadening, kResolution } stage = Stage::kResolution;
  /// What selects it in a format-3 project: "<record>.<field> [form]".
  std::string keyword;
  /// A table kernel's form is generated from `parameters` (Schema.cpp) and its values read into
  /// KernelParameters::values (BuildEffects.cpp); the older kernels have hand-written forms and named members.
  bool tableForm = false;
  std::vector<KernelParameter> parameters;
  bool usesTarget = false;  ///< reads the nodes' stopping and the point's target thickness
  FoldFunction<double> fold = nullptr;
  FoldFunction<LinearForm> weights = nullptr;
  ParameterDerivative derivative = nullptr;
  SupportFunction support = nullptr;
  // A broadening kernel: the broadened values (on double and LinearForm), their derivative in the temperature, and
  // the weights themselves (the transform of a transmission or yield applies them as a linear map).
  BroadenFunction<double> broaden = nullptr;
  BroadenFunction<LinearForm> broadenForm = nullptr;
  BroadenFunction<double> broadenDT = nullptr;
  BroadeningWeightsFunction broadeningWeights = nullptr;
  std::string doc;
};

/// The table of kernels (src/Kernels.cpp), in a fixed order.
const std::vector<Kernel> &Kernels();
/// The row named `name`, or null.
const Kernel *KernelNamed(const std::string &name);
/// The table kernel (Kernel::tableForm) whose keyword is "<record>.<field> <form>", or null.
const Kernel *KernelOfForm(const std::string &record, const std::string &field, const std::string &form);

/// The kernels of one effect: an optional broadening, then an optional resolution kernel (null: the broadened value
/// at the point).
struct EffectKernels {
  const Kernel *broadening = nullptr;
  const Kernel *resolution = nullptr;
};

/// Whether the broadening of an effect is applied in the fold or was applied by the caller already (the transform
/// of a transmission or yield broadens its terms itself before the transform).
enum class Broadening { kApply, kAlreadyApplied };

/*!
 * The effect's value at the point from its node values: the broadening at every node (or at the point, without a
 * resolution kernel), then the resolution kernel.  V = double gives the value, V = LinearForm (values the unit forms
 * of the nodes) its exact weights.  The one dispatch of the fold: one indirect call per point.
 */
double Fold(const EffectKernels &, Broadening, const KernelParameters &, const KernelPoint &, const KernelNodes &,
            const double *values);
LinearForm FoldForm(const EffectKernels &, Broadening, const KernelParameters &, const KernelPoint &,
                    const KernelNodes &, const LinearForm *values);
/// d(Fold)/dT for the node values `values`, T the broadening's temperature: the resolution kernel applied to the
/// broadening's derivative.  0 without a broadening or when it was applied by the caller.
double FoldTemperatureDerivative(const EffectKernels &, Broadening, const KernelParameters &, const KernelPoint &,
                                 const KernelNodes &, const double *values);
/// The broadening's weights at energy e (empty without a broadening).
void BroadeningWeights(const EffectKernels &, const KernelParameters &, const KernelPoint &, const KernelNodes &,
                       double e, std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT);
/// The c.m. interval [low, high] the nodes of a point must span for the effect: the resolution kernel's support (the
/// point itself without one), widened by the broadening's.
void EffectSupport(const EffectKernels &, const KernelParameters &, const SupportPoint &, double &high, double &low);

// Formulas the kernels share with the set-up (one implementation each).
/// The normalised Gaussian exp(-(e - centroid)^2 / 2 sigma^2) / (sigma sqrt(2 pi)).
double GaussianDensity(double e, double centroid, double sigma);
/// Straggling sigma (c.m. MeV) after a c.m. energy loss deltaE: c_s sqrt(dE[keV]) keV, c_s a c.m. coefficient.
double StragglingSigma(const KernelParameters &, double deltaE);
/// Bohr straggling (c.m. MeV) after an areal density (atoms/cm^2): Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x, times the
/// c.m./lab factor.  0 unless Bohr straggling is on.
double BohrStragglingSigma(const KernelParameters &, double arealDensity, double cmFactor, int z1, int z2);
/// The beam profile (sum of the skewed Gaussians) at a lab energy.
double BeamProfileWeight(const KernelParameters &, double labEnergy);
/// The lab interval outside which the beam profile is negligible.
void BeamProfileSupport(const KernelParameters &, double &low, double &high);
/// The Doppler width u (sqrt(MeV)) at temperature T (K) for a pair of c.m./lab factor f: u^2 = kT m/M.
double DopplerU(double cmFactor, double temperature);
/// The c.m. window around eCM beyond which the Doppler kernel at the layout temperature is negligible.
void DopplerWindow(const KernelParameters &, double eCM, double cmFactor, double &lo, double &hi);

/// The multiple of sigma a Gaussian fold reaches (and a THM (HOES) fold: thmRange).
constexpr double kConvolutionRange = 5.;
constexpr double kThmConvolutionRange = 5.;

}  // namespace experiment
}  // namespace azure2

#endif
