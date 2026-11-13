// Re-energising the points of a segment after an energy shift: the shifted energies of the points (a folded point
// keeps the nodes of the shared mesh and its kernel moves, ShiftedFold.h), the recalculation of their tables, and the channel tables of the last shift states that a shift returning to an earlier
// value reuses.
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "CNuc.h"
#include "ChannelTable.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "ShiftedFold.h"
#include "TargetEffect.h"

namespace {

/*!
 * Applies an energy shift to a point energy.
 *
 * The only thing guarded here is that the result stays strictly positive: the
 * lab->CM conversions and the Coulomb functions are undefined at or below zero.
 * The floor is a fixed fraction of the point's own energy, so it can only ever
 * engage for a shift large enough to push the point through threshold -- for
 * every physically meaningful shift the mapping is the exact identity plus the
 * shift, and therefore continuous in it.
 *
 * This used to be an absolute 0.01 MeV floor, which silently snapped *every*
 * data point below 10 keV onto the same energy the moment a segment's shift
 * became non-zero.  Since the update is skipped entirely while the shift is
 * still zero, that turned chi^2 into a step function of the shift for any data
 * set reaching below 10 keV (3H+d, 3He+d, ...): an infinitesimal shift moved
 * chi^2 by thousands, the finite-difference gradient was meaningless, and the
 * fit either froze at its starting point or ran away.
 */
double ShiftedEnergy(double originalEnergy, double shift) {
  static const double kMinEnergyFraction = 1.e-6;
  double shifted = originalEnergy + shift;
  double floorEnergy = originalEnergy * kMinEnergyFraction;
  return (shifted < floorEnergy) ? floorEnergy : shifted;
}

}  // namespace

/*!
 * Updates the energies of all points in this segment based on the current energy shift.
 */

void Segment::UpdatePointEnergiesWithShift(WorkingCopy &data, CNuc *theCNuc, const Config *configure) const {
  // No lock: each thread re-energises the points of its own working copy.
  for (int t = 0; t < NumTerms(); t++)
    for (int i = 0; i < NumPoints(); i++) UpdatePointEnergyWithShift(data, t, i, theCNuc, configure);
}

std::vector<double> Segment::ShiftState(const SegmentState &segment) {
  std::vector<double> state{segment.energyShift, segment.energyShiftSqrt};
  state.insert(state.end(), segment.shiftFnValues.begin(), segment.shiftFnValues.end());
  return state;
}

void Segment::BeginShiftRefill(WorkingCopy &data, int t, CNuc *theCNuc, const Config &configure) const {
  TermState &u = data.State(*this).terms[t];
  // Tables of another system (a session that changed a radius or the Coulomb route) are dropped.
  if (!u.shiftTables.empty() && !u.shiftTables.front().table->SameSystem(theCNuc, configure)) u.shiftTables.clear();
  u.shiftRows.assign(NumPoints(), {});
}

void Segment::CommitShiftTable(WorkingCopy &data, int t, CNuc *theCNuc, const Config &configure) const {
  SegmentState &segment = data.State(*this);
  TermState &u = segment.terms[t];
  if ((int)u.shiftRows.size() != NumPoints()) return;  // no refill since BeginShiftRefill
  std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>> all;
  for (auto &v : u.shiftRows) all.insert(all.end(), v.begin(), v.end());
  u.shiftRows.clear();
  if (all.empty()) return;
  const std::vector<double> state = ShiftState(segment);
  for (auto it = u.shiftTables.begin(); it != u.shiftTables.end(); ++it)
    if (it->state == state) {
      u.shiftTables.erase(it);
      break;
    }
  u.shiftTables.insert(u.shiftTables.begin(), TermState::ShiftTable{state, ChannelTable::FromRows(theCNuc, configure, std::move(all))});
  if (u.shiftTables.size() > kShiftTables) u.shiftTables.resize(kShiftTables);
}

ChannelRowSource Segment::ShiftRowSource(TermState &u, int i) {
  ChannelRowSource rows;
  for (const TermState::ShiftTable &table : u.shiftTables) rows.tables.push_back(table.table.get());
  if (i >= 0 && i < (int)u.shiftRows.size()) rows.used = &u.shiftRows[i];
  return rows;
}

/*!
 * Applies the segment's energy shift to point i (0-based) and its sub-points
 * and recalculates their energy-dependent values.  Writes only that point, so
 * different points may be updated concurrently.
 */
void Segment::UpdatePointEnergyWithShift(WorkingCopy &data, int t, int i, CNuc *theCNuc, const Config *configure) const {
  // An ordinary segment's energies are positive and must stay so: points at
  // E <= 0 are left alone and a shift is floored above zero (ShiftedEnergy).
  // Half-off-shell data have points below the entrance threshold, and all of
  // them shift (ShiftIsFloored): a point below threshold that never moved tore
  // its fold.
  const bool floored = this->ShiftIsFloored();
  SegmentState &state = data.State(*this);
  {
    Unit *point = &TermPoints(data, t)[i];
    if (point && (!floored || point->GetOriginalEnergy() > 0)) {
      double originalEnergy = point->GetOriginalEnergy();
      double shiftedEnergy = floored ? ShiftedEnergy(originalEnergy, TotalEnergyShift(state, originalEnergy))
                                     : originalEnergy + TotalEnergyShift(state, originalEnergy);

      // Set the shifted energy
      point->SetLabEnergy(shiftedEnergy);

      // Recalculate energy dependent values
      if (theCNuc && configure) {
        PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->TermEntranceKey(t)));
        PPair *exitPair = theCNuc->GetPairFromKey(this->TermExitKey(t));  // nullptr for total capture (-1)

        // The same conversion as Segment::Fill (ratios are dimensionless: no cross-section conversion; a leaf's
        // value is always converted)
        azure2::experiment::AngleLimit lim;
        if (!this->ConvertRowToCM(data, t, i, entrancePair, exitPair, t > 0 || !this->IsDimensionlessRatio(),
                                  BeyondRange::Convert, lim)) {
          // An energy shift moved the point past the kinematic maximum lab angle mid-fit: the c.m. angle keeps
          // its value as loaded and the run says so (it cannot stop here; audit 2026-10-25, item 31).
          configure->outStream << "WARNING: segment #" << this->TermKey(t) << ": the energy shift moved the point at lab angle "
                               << point->GetLabAngle() << " deg to E_lab = " << point->GetLabEnergy()
                               << " MeV, beyond its kinematic maximum " << (lim.maxSet ? lim.max : 180.0)
                               << " deg; its c.m. angle is kept at its value as loaded." << std::endl;
        }

        if (this->IsDifferential() || this->IsCMDifferential()) {
          point->ClearLegendrePolynomials();
          // Get TargetEffect with Q-coefficients if applicable
          const TargetEffect *effect = NULL;
          if (this->IsTargetEffect()) {
            const TargetEffect *te = data.Data().GetTargetEffect(this->TermTargetEffectNum(t));
            if (te && te->IsQCoefficients()) {
              effect = te;
            }
          }
          // the point alone: its row follows its angle; its slots are nodes of the shared mesh, which stay
          point->OwnSite().CalcLegendreP(&data, configure->maxLOrder, theCNuc, effect);
        }
        const ChannelRowSource rows = ShiftRowSource(state.terms[t], i);
        point->OwnSite().RecalcEDependentValues(&data, theCNuc, *configure, &rows);
        // The nodes stay and the kernel moves with the point (ShiftedFold.h): its row set does not depend on the shift
        // (the c.m. angle of a point given at a lab angle is the row set's at each node, a function of the node's
        // energy and the lab angle: RowSet.h).  Its target is taken at the shifted, true energy (the rule of
        // ShiftedFold.h): the window of a layer, the thick target's reach, the thickness its kernel reads.
        if (point->IsTargetEffect()) azure2::MoveTargetToShift(*point, data, *configure);
      }
    }
  }
}
