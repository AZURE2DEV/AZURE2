#ifndef AZURE2_ROWSET_H
#define AZURE2_ROWSET_H

#include <cstddef>

#include "Constants.h"
#include "SiteObservable.h"

namespace azure2 {

/*!
 * The angular rows of the readers of one block of angle-free slots (SlotGroup.h) that see the slots under one angle
 * (docs/upgrade/energy_grid.md section 3.2): what the kind's table row says a reader needs at each slot -- P_L(cos
 * theta) with the Q coefficients of the readers' target effect, and for elastic scattering the Coulomb amplitudes
 * f_d, f_x at the slot's energy.  Built once for all readers with the same angular identity (WorkingCopy::
 * BuildRowSets: the block, the kind, the c.m. and lab angle bits, the Q-coefficient effect, the entrance pair); a
 * reader holds the set and the position of its first slot in it (Site::Row, Site::CoulombAt).  Immutable once built,
 * shared by the copies of a working copy.
 *
 * The values a reader folds, the set's rows contracted with the slots' records, are those of the evaluation: they
 * live in the evaluation's buffer (Results::contracted, from `valueOffset`), written once per evaluation and set by
 * WorkingCopy::ContractRowSets between the energy pass and the folds (Unit::Calculate for a point calculated on its
 * own), and read by every reader (Site::SlotValue).  A set with valueOffset -1 is a reader's own (a point whose angle
 * moved with an energy shift: Unit::CalcNodeRows); its readers contract on each read.
 */
struct RowSet {
  SiteObservable kind = SiteObservable::kValue;
  /// The slots covered, [first, first + count) of WorkingCopy::SlotAt; a reader's slots are a range of them.
  int first = 0;
  int count = 0;
  /// Doubles per row (P_0 .. P_maxL).
  std::size_t rowLength = 0;
  /// One row per slot where the row depends on the energy (Site::RowDependsOnEnergy), else one row for all.
  bool rowPerSlot = false;
  vector_r rows;
  /// Elastic scattering: f_d, f_x per slot (Site::CoulombAmplitudesAt at the slot's energy), else empty.
  vector_c coulomb;
  /// An analysing power: Y_l^mu at the readers' angle (Polarization::HarmonicTable, l <= Config::maxLOrder, half-width
  /// muMax) and P_k^1(cos theta), k < rowLength; else empty.
  vector_c harmonics;
  int muMax = 0;
  vector_r p1;
  /// The values a reader folds per slot (the kind's, SiteObservableRow::values): 1, or 2 (cross section, A_y).
  int valuesPerSlot = 1;
  /// Where the contracted values of slot `first` lie in Results::contracted (valuesPerSlot per slot), -1: not cached.
  long valueOffset = -1;

  /// The row at slot k of the set (0-based).
  const double *Row(int k) const { return rows.data() + (rowPerSlot ? (std::size_t)k * rowLength : 0); }
  void CoulombAt(int k, complex &direct, complex &exchange) const {
    direct = coulomb[2 * (std::size_t)k];
    exchange = coulomb[2 * (std::size_t)k + 1];
  }
  /// Everything a reader brings to slot k (SiteObservable.h; the angle itself is the reader's).
  AngularRow RowAt(int k) const {
    AngularRow row;
    row.p = Row(k);
    row.length = (int)rowLength;
    if (!coulomb.empty()) CoulombAt(k, row.fd, row.fx);
    if (!harmonics.empty()) row.y = harmonics.data(), row.muMax = muMax;
    if (!p1.empty()) row.p1 = p1.data();
    return row;
  }
  /// The index in Results::contracted of value `i` of slot k.
  std::size_t ValueIndex(int k, int i) const { return (std::size_t)valueOffset + (std::size_t)k * valuesPerSlot + i; }
};

}  // namespace azure2

#endif
