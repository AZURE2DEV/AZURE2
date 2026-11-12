#ifndef AZURE2_ROWSET_H
#define AZURE2_ROWSET_H

#include <cstddef>

#include "Constants.h"
#include "SiteObservable.h"
#include "azure2/experiment/Frames.h"

namespace azure2 {

/*!
 * The angular rows of the readers of one block of angle-free slots (SlotGroup.h) that see the slots under one angle
 * (docs/upgrade/energy_grid.md section 3.2): what the kind's table row says a reader needs at each slot -- P_L(cos
 * theta) with the Q coefficients of the readers' target effect, for elastic scattering the Coulomb amplitudes f_d, f_x
 * at the slot's energy, for an analysing power Y_l'^mu and P_k^1 at the angle.  Built once for all readers with the
 * same angular identity (WorkingCopy::BuildRowSets: the block, the kind, the entrance pair, the Q-coefficient effect
 * and the angle as the data give it -- the c.m. angle bits, or for data given at a lab angle the lab angle bits); a
 * reader holds the set and the position of its first slot in it (Site::Row, Site::CoulombAt).  Immutable once built,
 * shared by the copies of a working copy.
 *
 * Data given at a lab angle (experiment::AngleIsLab: lab-frame differential data with a particle or photon exit) are
 * folded in the lab (energy_grid.md section 4.3, decision D5): at slot n the readers' c.m. angle is theta_cm(E_n;
 * theta_lab) and the value is divided by J(E_n; theta_lab) = dOmega_lab/dOmega_cm, both from experiment::LabAngleAt
 * at the slot's lab energy, so that a reader's fold is the lab fold sum_n w_n O_n(theta_cm,n) / J_n; the reader
 * multiplies it by its own Jacobian (Unit::LabFoldFactor), which its comparison divides again.
 *
 * One rule for what is held per slot: a factor that depends on the slot's energy is per slot.  f_C (and J) always do;
 * the angular factors (the angle, P_L, Y_l'^mu, P_k^1) do when `perSlot` is set -- the Q transform of a particle exit
 * makes P_L depend on the energy (Site::RowDependsOnEnergy), a lab angle makes the c.m. angle do -- and are one for
 * all slots otherwise.  The accessors take the slot and apply the stride.
 *
 * The values a reader folds, the set's rows contracted with the slots' records (divided by J for a lab angle), are
 * those of the evaluation: they live in the evaluation's buffer (Results::contracted, from `valueOffset`), written once
 * per evaluation and set by WorkingCopy::ContractRowSets between the energy pass and the folds (Unit::Calculate for a
 * point calculated on its own), and read by every reader (Site::SlotValue).  A set with valueOffset -1 is a reader's
 * own (a point whose slots move with an energy shift, off the shared mesh: Unit::CalcNodeRows); its readers contract
 * on each read.
 */
struct RowSet {
  SiteObservable kind = SiteObservable::kValue;
  /// The slots covered, [first, first + count) of WorkingCopy::SlotAt; a reader's slots are a range of them.
  int first = 0;
  int count = 0;
  /// The angular factors per slot (see above), else one for all slots.
  bool perSlot = false;
  /// The readers' c.m. angle (deg).
  vector_r angles;
  /// Doubles per row (P_0 .. P_maxL).
  std::size_t rowLength = 0;
  vector_r rows;
  /// Elastic scattering: f_d, f_x per slot (Site::CoulombAmplitudesAt at the slot's energy and angle), else empty.
  vector_c coulomb;
  /// An analysing power: Y_l^mu at the angle (Polarization::HarmonicTable, l <= Config::maxLOrder, half-width muMax,
  /// harmonicsLength values) and P_k^1(cos theta), k < rowLength; else empty.
  vector_c harmonics;
  std::size_t harmonicsLength = 0;
  int muMax = 0;
  vector_r p1;
  /// Data given at a lab angle: the readers' conversion (their term's, experiment::LabAngleAt), their lab angle
  /// (deg) and dOmega_lab/dOmega_cm per slot; `jacobian` is empty otherwise (J = 1).
  bool labFrame = false;
  experiment::FrameRule frame;
  double thetaLab = 0.0;
  vector_r jacobian;
  /// The values a reader folds per slot (the kind's, SiteObservableRow::values): 1, or 2 (cross section, A_y).
  int valuesPerSlot = 1;
  /// Where the contracted values of slot `first` lie in Results::contracted (valuesPerSlot per slot), -1: not cached.
  long valueOffset = -1;

  /// The position of slot k's angular factors (0 when they are one for all).
  std::size_t At(int k) const { return perSlot ? (std::size_t)k : 0; }
  /// The c.m. angle (deg) at slot k (0-based).
  double AngleAt(int k) const { return angles[At(k)]; }
  /// The row at slot k of the set (0-based).
  const double *Row(int k) const { return rows.data() + At(k) * rowLength; }
  void CoulombAt(int k, complex &direct, complex &exchange) const {
    direct = coulomb[2 * (std::size_t)k];
    exchange = coulomb[2 * (std::size_t)k + 1];
  }
  /// dOmega_lab/dOmega_cm at slot k: the value a reader folds is the contraction divided by it.
  double JacobianAt(int k) const { return jacobian.empty() ? 1.0 : jacobian[k]; }
  /// Everything a reader brings to slot k (SiteObservable.h).
  AngularRow RowAt(int k) const {
    AngularRow row;
    row.p = Row(k);
    row.length = (int)rowLength;
    if (!coulomb.empty()) CoulombAt(k, row.fd, row.fx);
    if (!harmonics.empty()) row.y = harmonics.data() + At(k) * harmonicsLength, row.muMax = muMax;
    if (!p1.empty()) row.p1 = p1.data() + At(k) * rowLength;
    row.theta = AngleAt(k) * pi / 180.0;
    return row;
  }
  /// The index in Results::contracted of value `i` of slot k.
  std::size_t ValueIndex(int k, int i) const { return (std::size_t)valueOffset + (std::size_t)k * valuesPerSlot + i; }
};

}  // namespace azure2

#endif
