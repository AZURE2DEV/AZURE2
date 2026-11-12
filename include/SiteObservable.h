#ifndef AZURE2_SITEOBSERVABLE_H
#define AZURE2_SITEOBSERVABLE_H

#include <cstdint>

#include "Config.h"

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
};

/// The layout of the angular part of a kElastic record (offsets from its first slot, which holds K).
namespace elastic_record {
constexpr int kOrders = Config::maxLOrder + 1;  ///< B_L and D_l, X_l slots each (l <= maxLOrder)
constexpr int kCoulomb = 1 + kOrders;           ///< c_0 = geom / 100 (geom = pi/k^2 in fm^2: c_0 |f|^2 in b/sr)
constexpr int kCross = kCoulomb + 1;            ///< c_x = 2 c_0 (-1)^(2j) / (2j + 1), 0 without spin-dependent exchange
constexpr int kNumL = kCoulomb + 2;             ///< n: D_l, X_l for l = 0 .. n-1
constexpr int kDirect = kCoulomb + 3;           ///< D_l: re, im
constexpr int kExchange = kDirect + 2 * kOrders;///< X_l: re, im
constexpr int kSlots = kExchange + 2 * kOrders - 1;  ///< the angular capacity (slots after the count)
}  // namespace elastic_record

/// One row of the table of kinds.
struct SiteObservableRow {
  SiteObservable kind;
  const char *name;
  /// Slots of the angular part of the record (Results::RecordSize), 0: none.
  int angularSlots;
  /// The reader's value from the angular part `r` (r[0] the count of B_L) at its slot `slot` (Site::RowP); null for
  /// kValue, whose value is the record's sigma.
  double (*contract)(const double *r, const Site &reader, int slot);
};

/// The rows, indexed by the kind.
extern const SiteObservableRow kSiteObservables[3];

/// The kind of the sites of a term.
SiteObservable ObservableOf(const TermSettings &settings);

#endif
