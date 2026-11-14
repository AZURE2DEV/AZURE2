#ifndef AZURE2_SITEOBSERVABLE_H
#define AZURE2_SITEOBSERVABLE_H

#include <cstdint>

#include "Config.h"
#include "Constants.h"

class CompoundSystem;
class Site;
struct TermSettings;

/*!
 * What a site's result record holds: the observables of a reaction at the site's energy (docs/upgrade/energy_grid.md
 * section 3), one kind per site, chosen by its term (ObservableOf).  A kind is one row of kSiteObservables: the slots
 * its record needs in the angular part (Results.h), and how a reader's angular row contracts them into the reader's
 * value (Site::ValueFor).  The record is written by GenMatrixFunc::CalculateBareCrossSection (the formula at E_n, one
 * branch per kind) and read back in reverse by AMatrixFunc::PointTBar (its adjoint).  A kind without an angle serves
 * every reader whatever its angle: the sharing key of the mesh has no angle for it (MeshSlots.cpp, MeshGroupKey).
 *
 * A new observable: a kind here, its row (slots, contraction) in SiteObservables.cpp, its formula in GenMatrixFunc,
 * its adjoint in AMatrixFunc::PointTBar.
 *
 * THM, the half-off-shell cross section (docs/upgrade/THM.md section 1; energy_grid.md section 3.2, row "THM
 * sigma_HOES").  What a site of a THM term holds at E_n, and what its reader brings:
 *
 *   sigma_HOES(E_n) = W(E_n) R(E_n) sum_J (2J+1) sum_c' K(E_n) 2 P_c'(E_f) sum_k w_k(E_n) sum_(s,l)
 *                     | sum_lambda,mu gamma_lambda,c' [N_C] A_lambda,mu(E_n) W_mu,(s,l)(p_k(E_n)) |^2
 *
 * Every factor is a function of the node energy, the reaction (the entrance and exit pairs) and the THM experiment:
 * A from the level matrix at E_n; P_c' and the exit phases from the node's channel rows; K(E) the <thm> kinematics
 * option; the vertex W of the experiment's model (plane wave with j_l, rho j_l' and the Coulomb term at p(E_n),
 * spectatorEnergy per entrance pair; full form factor; DW tables evaluated at E_n); the spectator-momentum window's
 * nodes p_k and weights w_k at E_n (ThmSpectatorWindow::Nodes: a function of E_n and the experiment, the fallback to
 * the nearest data energy is the experiment's list); the line shape N_C; the angular window theta= of the exit pair
 * (an average over the experiment's interval, ThmAngleWindow::Mean); the coherent background; the experiment's
 * distortion R(E_n) and the data set's weight table W(E_n), multiplied at the node.  Nothing reads the data point: not
 * its angle (the spectator angles of spectatorAngle=, spectatorAngles=, theta= are the experiment's), not its row.
 * So a THM site is a scalar without an angle: one number per node, the record's sigma, read by every reader as it is
 * (no contraction, one value per slot).  It is the kind kHoes, a row of kSiteObservables of its own so that the energy
 * pass and the gradient choose the HOES formula and its adjoint by the kind of the site, not by asking IsTHM(); the
 * reader's row is the identity, and the factors that belong to the point (the Gaussian resolution, a beam profile)
 * are its kernel's (Kernels.h).  W and R stay on the site: they are functions of E_n and the data set's or the
 * experiment's settings, and multiplying them at the node is today's arithmetic, in today's order.
 */
enum class SiteObservable : std::uint8_t {
  /// A scalar at the site's energy, without an angle (sigma, a phase shift, the total cross section): the record's
  /// sigma, read as it is by every reader.
  kValue,
  /// B_L (b/sr) of dsigma/dOmega_c.m. = sum_L B_L P_L(cos theta), L = 0 .. K-1, and A_0 = sigma_int / 4 pi of an
  /// angular distribution; no angle.  Every reaction whose exit pair is not its entrance pair, and the angular
  /// distributions (elastic ones: the nuclear amplitude's part).
  kLegendre,
  /// Elastic scattering: the nuclear B_L^NN as kLegendre, the factors of the Coulomb term and the Coulomb-nuclear
  /// coefficients D_l, X_l (complex) of dsigma/dOmega_c.m. = c_0 (|f_d|^2 + |f_x|^2) + c_x Re(f_d f_x*) +
  /// sum_L B_L P_L + Re sum_l (f_d D_l + f_x X_l) P_l, with f_d, f_x the reader's Coulomb amplitudes at the site's
  /// energy (Site::CoulombAt: f_d the (symmetrized) f_C(theta), f_x = f_C(pi - theta) for identical particles with
  /// spin, else 0); no angle.
  kElastic,
  /// The vector analysing power A_y and the differential cross section that weights it in a fold: the record of the
  /// cross section (kElastic's for elastic scattering, else kLegendre's) and the coefficients of A_y, without an angle
  /// (ay_record).  The reader's row adds the angle; its contraction gives the cross section and A_y.
  kAnalyzingPower,
  /// The half-off-shell cross section of THM data (above): a scalar at the site's energy, without an angle, the record's
  /// sigma read as it is.  Written by the HOES formula (THMMatrixFunc.h, HoesCrossSection) on the level matrix of the
  /// node, differentiated by its adjoint (THMMatrixFunc::THMPointAdjoint).
  kHoes,
};

/// The layout of the angular part of a kElastic record, packed after its counts: r[kNumB] = n B_L, r[kNumL] = m D_l
/// and X_l, c_0, c_x, then B_0 .. B_{n-1}, D_0 .. D_{m-1} (re, im), X_0 .. X_{m-1} (re, im).
namespace elastic_record {
constexpr int kNumB = 0;     ///< n
constexpr int kNumL = 1;     ///< m
constexpr int kCoulomb = 2;  ///< c_0 = geom / 100 (geom = pi/k^2 in fm^2: c_0 |f|^2 in b/sr)
constexpr int kCross = 3;    ///< c_x = 2 c_0 (-1)^(2j) / (2j + 1), 0 without spin-dependent exchange
constexpr int kFirst = 4;    ///< B_0
/// The slots after the count for at most `orders` B_L and D_l, X_l.
constexpr int Slots(int orders) { return kFirst - 1 + 5 * orders; }
}  // namespace elastic_record

/*!
 * The layout of the angular part r of a kAnalyzingPower record.  r[kWeight] = W: r + W is the angular part of the
 * cross section's record (kElastic or kLegendre, its count first); r[kForm] the form of the coefficients:
 * r[kElastic] = 1 where the cross section's record is kElastic's (elastic scattering), else 0.
 *   kParticle (Seyler 1969, Polarization::AmplitudeMatrix): r[3..6] = 2 j1, 2 j2, 2 j1', 2 j2' (entrance, exit; j1 the
 *     projectile's; elastic: the Coulomb amplitude on the diagonal), r[7] = n_a amplitudes,
 *     r[8] = n_l exit orbital momenta, from r[9] the complex C(a, l') = sum_{J,l} i sqrt(2l+1) (s l v 0|J v)
 *     (s' l' v' mu|J v) T identical, (re, im) at 9 + 2 (a n_l + l'): M_a = sum_l' C(a, l') Y_l'^mu - f_C(s) delta.
 *   kCapture (Seyler and Weller 1979, CaptureAyTerm): r[3] = K, from r[4] the K coefficients d_k = sum a_k Re(T1 T2*)
 *     of the denominator (with P_k), then the K coefficients n_k = sum b_k Re(i T1 T2*) of the numerator (with P_k^1):
 *     A_y = sum_k n_k P_k^1 / sum_k d_k P_k.
 */
namespace ay_record {
constexpr int kWeight = 0;
constexpr int kForm = 1;
constexpr int kParticle = 0, kCapture = 1;
constexpr int kElastic = 2;  ///< 1: the cross section's record is kElastic's (and the Coulomb amplitude is on M's diagonal)
constexpr int kSpins = 3, kNumAmplitudes = 7, kNumLp = 8, kCoefficients = 9;
constexpr int kNumOrders = 3, kOrders = 4;
}  // namespace ay_record

/// What a record must hold room for: the reactions' Legendre orders (one more than the largest) and the largest
/// number of channel-spin amplitudes of a particle reaction (WorkingCopy::RecordCapacity).
struct RecordCapacity {
  int orders = Config::maxLOrder + 1;
  int amplitudes = 0;
};

/// The capacity the records of the reactions of `compound` need: one more than their largest Legendre order (the
/// pathway table's exits; at most Config::maxLOrder + 1) and the amplitudes (s v, s' v') of a reaction between
/// particle pairs, ((2 j1 + 1)(2 j2 + 1))^2 for the largest pair.
RecordCapacity RecordCapacityOf(const CompoundSystem &compound);

/// What a reader brings to a slot: its angular row (P_L with its Q coefficients, Site::Row) and, for elastic
/// scattering, its Coulomb amplitudes at the slot's energy (Site::CoulombAt).
struct AngularRow {
  const double *p = nullptr;
  int length = 0;  ///< P_0 .. P_{length-1}
  complex fd = 0.0, fx = 0.0;
  /// An analysing power: Y_l^mu (Polarization::HarmonicTable, half-width muMax) and P_k^1(cos theta), k < length.
  const complex *y = nullptr;
  int muMax = 0;
  const double *p1 = nullptr;
  /// The reader's c.m. angle (rad), for the A_y factors of a reader without them (a point that is its own site).
  double theta = 0.0;
};

/// One row of the table of kinds.
struct SiteObservableRow {
  SiteObservable kind;
  const char *name;
  /// Slots of the angular part of the record after its count (Results::RecordSize) for a reaction whose Legendre
  /// orders and amplitudes are at most those of `capacity` (WorkingCopy::Capacity), 0: none.
  int (*angularSlots)(const RecordCapacity &);
  /// The value of the angular part `r` of a record under the angular row `row`, and, for a kind with two (A_y), the
  /// second into `second` (not null then); null for kValue, whose value is the record's sigma.  (A reader of one
  /// Legendre coefficient reads it without a row: Site::ValueFor.)
  double (*contract)(const double *r, const AngularRow &row, double *second);
  /// The values per slot a reader folds: 1, or 2 -- a ratio: value 0 the weight w (the cross section), value 1 the
  /// ratio's element r (A_y), folded as (K.(w r))/(K.w) (Unit::FoldsRatio, RatioPartsAt; its adjoint in AZUREGrad.cpp).
  int values;
  /// Multiply every observable of the angular part `r` (r[0] its count) by f (electron screening).
  void (*scale)(double *r, double f);
  /// NaN in the angular part `r` of `slots` slots after the count, so that every reader of it is NaN (a failed site).
  void (*poison)(double *r, int slots);
};

/// The rows, indexed by the kind.
extern const SiteObservableRow kSiteObservables[5];

/// Whether a kind's readers bring an angular row to its record (SiteObservableRow::contract): B_L, the elastic
/// amplitudes, A_y.  A scalar kind (kValue, kHoes) is read as it is.
inline bool ReadThroughRow(SiteObservable kind) { return kSiteObservables[(int)kind].contract != nullptr; }

/// The kind of the sites of a term.
SiteObservable ObservableOf(const TermSettings &settings);

#endif
