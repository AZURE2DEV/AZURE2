#ifndef AZURE2_SITEOBSERVABLE_H
#define AZURE2_SITEOBSERVABLE_H

#include <cstdint>

#include "Config.h"
#include "Constants.h"

class CNuc;
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
 */
enum class SiteObservable : std::uint8_t {
  /// The finished value at the site's own energy and angle (sigma, a phase shift, the total cross section, A_y, THM):
  /// the reader's angle is the site's, so a site of this kind is keyed by the angle.
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
RecordCapacity RecordCapacityOf(const CNuc &compound);

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
extern const SiteObservableRow kSiteObservables[4];

/// The kind of the sites of a term.
SiteObservable ObservableOf(const TermSettings &settings);

#endif
