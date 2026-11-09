#ifndef AZURE2_RESULTS_H
#define AZURE2_RESULTS_H

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace azure2 {

/*!
 * What one evaluation of the model produces (docs/upgrade/data_model.md section 2.1, decision M2 (c)): the values
 * every calculation unit (a data or test point, and each of its sub-points) computes for the current level
 * parameters and nuisances.  A unit owns one record, a fixed range of `slot` at an offset assigned when the units
 * are laid out (EData::LayoutResults); nothing is looked up or allocated while the points are calculated.
 *
 * A record holds, in this order:
 *   kSigma    the unit's observable: sigma (b), dsigma/dOmega_c.m. (b/sr), a phase shift (deg), A_y (1), a
 *             Legendre coefficient a_k (1) or A_k (b/sr), the THM cross section; after the target fold for a
 *             point that carries one; the composed (SUM, RATIO, sample, transform) value once the segment's
 *             combination has run
 *   kSigmaE1  the E1 part of an angle-integrated capture cross section (b)
 *   kSigmaE2  the E2 part (b)
 *   kAy       the analysing power A_y (1), kept beside the cross section, which weights its target average
 *   kA0       A_0 = sigma_int / 4 pi (b/sr) of an angular-distribution unit
 *   then, for a unit with an angular-distribution capacity K > 0: the number of coefficients written, then K values
 *   a_0 .. a_{K-1}; then, for the point of a sample (natural) segment with N > 0 nuclei: the number of terms written,
 *   then N abundance-weighted terms, one per compound nucleus.
 *
 * The struct is plain (no engine headers), so it can move to the experiment library (Phase 4, S6/S7).
 */
struct Results {
  enum Quantity : std::size_t { kSigma = 0, kSigmaE1 = 1, kSigmaE2 = 2, kAy = 3, kA0 = 4, kNumFixed = 5 };

  std::vector<double> slot;

  /// Length of a record with angular capacity nAng and nNat natural terms.
  static std::size_t RecordSize(std::size_t nAng, std::size_t nNat) {
    return kNumFixed + (nAng ? 1 + nAng : 0) + (nNat ? 1 + nNat : 0);
  }
  /// Offset, inside a record, of the angular-distribution count (the values follow it).
  static std::size_t AngularOffset() { return kNumFixed; }
  /// Offset, inside a record, of the natural-term count (the terms follow it).
  static std::size_t NaturalOffset(std::size_t nAng) { return kNumFixed + (nAng ? 1 + nAng : 0); }

  /// Every slot quiet NaN: a slot read before the evaluation writes it then shows as NaN.
  void Poison() { std::fill(slot.begin(), slot.end(), std::numeric_limits<double>::quiet_NaN()); }

  /// AZURE_POISON_RESULTS=1: reused result storage is poisoned on entry to every evaluation.
  static bool PoisonRequested() {
    static const bool on = [] {
      const char *v = std::getenv("AZURE_POISON_RESULTS");
      return v && *v && std::strcmp(v, "0") != 0;
    }();
    return on;
  }
};

/// Copy a record laid out with (fromAng, fromNat) slots into one with (toAng, toNat): the fixed quantities, then
/// as much of the angular and natural parts as both hold (the counts included).
inline void CopyRecord(const double *from, std::size_t fromAng, std::size_t fromNat, double *to, std::size_t toAng,
                       std::size_t toNat) {
  std::copy(from, from + Results::kNumFixed, to);
  if (fromAng && toAng) {
    const std::size_t n = std::min(fromAng, toAng);
    std::copy(from + Results::AngularOffset(), from + Results::AngularOffset() + 1 + n, to + Results::AngularOffset());
    if (to[Results::AngularOffset()] > (double)toAng) to[Results::AngularOffset()] = (double)toAng;
  }
  if (fromNat && toNat) {
    const std::size_t n = std::min(fromNat, toNat);
    const double *f = from + Results::NaturalOffset(fromAng);
    double *t = to + Results::NaturalOffset(toAng);
    std::copy(f, f + 1 + n, t);
    if (t[0] > (double)toNat) t[0] = (double)toNat;
  }
}

}  // namespace azure2

#endif
