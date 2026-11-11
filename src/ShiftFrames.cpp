// Re-energising the points of a segment after an energy shift: the shifted energies of the points and their
// sub-points, the recalculation of their tables, and the channel tables of the last shift states that a shift
// returning to an earlier value reuses.
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "CNuc.h"
#include "ChannelTable.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"
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
 * The earlier point of this segment that point i (0-based) is mapped onto
 * (WorkingCopy::MapData), or NULL -- also for a point mapped onto a point of another
 * segment, which this segment's shift step re-energises on its own.
 */
Unit *Segment::ShiftHost(WorkingCopy &data, int t, int i) const {
  const UnitRange<Unit> points = TermPoints(data, t);
  if (i < 0 || i >= points.size()) return nullptr;
  Unit *point = &points[i];
  if (!point->IsMapped()) return nullptr;
  EnergyMap map = point->GetMap();
  if (map.point < 1 || map.point > i) return nullptr;
  Unit *candidate = &points[map.point - 1];
  if (candidate->NumSubPoints() != point->NumSubPoints()) return nullptr;
  const azure2::MappedPoints mapped = candidate->PlanMapped(&data);
  for (int m = 0, n = mapped.size(); m < n; m++)
    if (mapped[m] == point) return candidate;
  return nullptr;
}

/*!
 * Re-energises point i and then every point mapped onto it, on the calling
 * thread.  A mapped point copies its energy-dependent values from its host
 * (UpdatePointEnergyWithShift), so the host must be finished first: the
 * parallel re-energise loops (WorkingCopy::FillEnergyShiftsFromParams, the shift
 * columns of AZURECalc::ResidualJacobian) hand out these groups, never single
 * mapped points.  A point with a host returns at once.
 */
void Segment::UpdatePointGroupEnergyWithShift(WorkingCopy &data, int t, int i, CNuc *theCNuc, const Config *configure) const {
  if (ShiftHost(data, t, i)) return;
  UpdatePointEnergyWithShift(data, t, i, theCNuc, configure);
  const UnitRange<Unit> points = TermPoints(data, t);
  Unit *point = &points[i];
  // A point mapped from another segment (WorkingCopy::MapData) is not this
  // segment's to move: it is re-energised with its own segment, on its own.
  std::less<const Unit *> before;
  const Unit *first = &points.front(), *last = &points.back();
  const azure2::MappedPoints pointMapped = point->PlanMapped(&data);
  for (int m = 0, n = pointMapped.size(); m < n; m++) {
    const Unit *mapped = pointMapped[m];
    if (before(mapped, first) || before(last, mapped)) continue;
    const int j = (int)(mapped - first);
    if (j > i && j < points.size() && ShiftHost(data, t, j) == point)
      UpdatePointEnergyWithShift(data, t, j, theCNuc, configure);
  }
}

/*!
 * Applies the segment's energy shift to point i (0-based) and its sub-points
 * and recalculates their energy-dependent values.  Writes only that point, so
 * different points may be updated concurrently -- except that a mapped point
 * reads its host (ShiftHost), so it must run after its host and not
 * concurrently with it: UpdatePointGroupEnergyWithShift.
 */
void Segment::UpdatePointEnergyWithShift(WorkingCopy &data, int t, int i, CNuc *theCNuc, const Config *configure) const {
  // An ordinary segment's energies are positive and must stay so: points at
  // E <= 0 are left alone and a shift is floored above zero (ShiftedEnergy).
  // A THM (HOES) segment has points and folding sub-points below the entrance
  // threshold, and all of them shift: leaving the sub-points at E <= 0 behind
  // (as the guard did) tore the fold of every point within its reach of
  // threshold, and a point below threshold never moved at all.
  const bool thm = this->IsTHM();
  SegmentState &state = data.State(*this);
  {
    Unit *point = &TermPoints(data, t)[i];
    // A point mapped onto an earlier point of this segment (WorkingCopy::MapData)
    // is at the same energy as its host, which has already been moved (index
    // order in UpdatePointEnergiesWithShift; host first within a group in
    // UpdatePointGroupEnergyWithShift): its energy-dependent values (and its
    // sub-points') are copied from the host instead of recomputed.  Only the
    // angle-dependent Coulomb amplitude is its own.  CopyEDependentValuesFrom
    // refuses (and the values are recomputed) if the host is not at this
    // energy.
    Unit *host = ShiftHost(data, t, i);
    if (point && (thm || point->GetOriginalEnergy() > 0)) {
      double originalEnergy = point->GetOriginalEnergy();
      double shiftedEnergy = thm ? originalEnergy + TotalEnergyShift(state, originalEnergy)
                                 : ShiftedEnergy(originalEnergy, TotalEnergyShift(state, originalEnergy));

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
          point->CalcLegendreP(&data, configure->maxLOrder, theCNuc, effect);
        }
        // The sub-points are moved and recalculated below, one by one.
        const ChannelRowSource rows = ShiftRowSource(state.terms[t], i);
        if (host && point->CopyEDependentValuesFrom(*host))
          point->CalcCoulombAmplitude(&data, theCNuc, false);
        else
          point->RecalcEDependentValues(&data, theCNuc, *configure, false, &rows);

        // IMPORTANT: Also apply energy shift to all subpoints (used in convolution/target integration)
        for (int j = 1; j <= point->NumSubPoints(); j++) {
          Unit *subPoint = point->GetSubPoint(data, j);
          if (subPoint && (thm || subPoint->GetOriginalEnergy() > 0)) {
            double subOriginalEnergy = subPoint->GetOriginalEnergy();
            // Sub-points hold c.m. energies and are the depth / resolution samples
            // of their parent point: they move RIGIDLY with it, by the parent's
            // own lab shift converted to the c.m. frame.  A shift is a calibration
            // of the measured beam energy; the target thickness, the beam profile
            // and the depth grid are physical and keep their width.  Before
            // 2026-10-26 each sub-point was moved by the shift at its own energy,
            // which for a sqrt(E) shift or a shift function stretched the window by
            // (1 + ds/dE) -- a yield factor the fit could trade against b -- and put
            // the resonance-anchored lattice (WorkingCopy::ReplaceSubPointGrids, which
            // offsets by the parent's shift) off a narrow resonance by s' x depth
            // (audit item 11, tests/shift_sqrt_subpoints).  A constant shift was
            // never affected.
            // A beta-delayed segment (entrance ptype 20) carries the light decay
            // particle's lab energy; c.m. = lab*(M1+M2)/M2 of the EXIT pair
            // (Unit::ConvertDecayEnergy), not the entrance-pair ratio.
            const double cmFactor = CMEnergyFactor(entrancePair, exitPair);
            // A THM segment's shift is not floored (its energies may lie below threshold), as for the point.
            double energyShiftCM = cmFactor * TotalEnergyShift(state, originalEnergy);
            double subShiftedEnergy = thm ? subOriginalEnergy + energyShiftCM
                                          : ShiftedEnergy(subOriginalEnergy, energyShiftCM);

            // Set the shifted energy for subpoint
            // Must set both CM and Lab energy since CalcLegendreP uses GetLabEnergy() for Q-coefficients
            subPoint->SetCMEnergy(subShiftedEnergy);
            subPoint->SetLabEnergy(subShiftedEnergy);

            // Recalculate energy dependent values for subpoint
            if (entrancePair->GetPType() == 20) {
              // The c.m. energy was set above. ConvertDecayEnergy would read the lab
              // slot (now holding the c.m. value) and multiply by (M1+M2)/M2 again,
              // pushing every convolution sub-point of a beta-delayed spectrum up by
              // 4/3 and corrupting the folded yield; only the excitation energy needs updating.
              subPoint->ConvertExcitationEnergy(exitPair);
            } else if (this->IsCMDifferential()) {
              // subPoint->ConvertLabEnergy(entrancePair);
            } else if (!this->IsCMDifferential()) {
              // subPoint->ConvertLabEnergy(entrancePair);
            }

            // A sub-point is built at its parent's c.m. angle (WorkingCopy::SetUpTargetEffects, Unit(angle, energy,
            // parent) puts that angle in both slots) and carries no data of its own, so the lab-to-c.m. angle and
            // cross-section conversions of the parent (above) must not be applied to it.  Before 2026-10-26 they
            // were, which took the parent's c.m. angle for a lab angle and moved every sub-point by the full
            // lab-c.m. difference (2-5 deg for p + light nuclei, 18 deg for alpha + 12C): a lab-frame differential
            // segment with target/resolution integration and any active shift (nonzero, varied, shiftgroup,
            // sqrtshift, shiftfn) fitted at the wrong angles, and its shift derivatives were wrong with it
            // (16O(p,p) with a 0.1 keV convolution: chi2 13,484 -> 23,097 at a 1e-7 MeV shift).  The parent's
            // c.m. angle at the shifted energy is the right value; tests/shift_subpoint_angle.
            if (this->IsDifferential() && !this->IsCMDifferential()) {
              subPoint->SetLabAngle(point->GetCMAngle());
              subPoint->SetCMAngle(point->GetCMAngle());
            }

            if (this->IsDifferential() || this->IsCMDifferential()) {
              subPoint->ClearLegendrePolynomials();
              // Get TargetEffect with Q-coefficients if applicable
              const TargetEffect *subEffect = NULL;
              if (this->IsTargetEffect()) {
                const TargetEffect *te = data.Data().GetTargetEffect(this->TermTargetEffectNum(t));
                if (te && te->IsQCoefficients()) {
                  subEffect = te;
                }
              }
              subPoint->CalcLegendreP(&data, configure->maxLOrder, theCNuc, subEffect);
            }
            if (host && subPoint->CopyEDependentValuesFrom(*host->GetSubPoint(data, j)))
              subPoint->CalcCoulombAmplitude(&data, theCNuc);
            else
              subPoint->RecalcEDependentValues(&data, theCNuc, *configure, true, &rows);
          }
        }
      }
    }
  }
}
