#ifndef AZURE2_SLOTGROUP_H
#define AZURE2_SLOTGROUP_H

#include <cstdint>
#include <vector>

#include "Site.h"

/*!
 * The slots of folded points, shared (docs/upgrade/energy_grid.md section 2.6).
 *
 * A folded point is a row of weights over the cross sections at the nodes of its range: y_i = sum_j K_ij sigma_j.
 * Where several points of one data set need the cross section at the same node with everything else that enters it
 * equal, it is calculated once.  Such points form a sharing group; the group holds one Site per node of the union of
 * its points' ranges, in node order (descending energy), as one block of WorkingCopy::SlotAt; each point's slots are a
 * contiguous range of that block (Unit::SubPointFirst, NumSubPoints), and the ranges of two points of a group overlap.
 *
 * The key of a group (SlotGroupKey), from the audit of Site's members (what a site holds and what it depends on):
 *   - node energy: the c.m., lab, compound and nominal energies, the channel rows (L_o, sqrt(P), phases), the
 *     geometrical and S-factor conversion factors, the screening k, the external-capture amplitudes -- one per node;
 *   - the term: what the energy pass reads of its settings (the observable, the pairs, UPOS L / I_c / delta, the E1/E2
 *     component, THM; TermSettings.h, SameObservable), the frame shift E_x - E_cm of the pair, the external-capture
 *     pathways -- in the key; not the segment: two data sets of one reaction on one mesh share their sites;
 *   - angle: not in the key.  A differential site on a mesh holds its observables without an angle (SiteObservable.h)
 *     and its readers bring their rows (RowSet.h); a scalar site (sigma, a phase shift, the total cross section) has
 *     none; THM, keyed by angle, is off the mesh (a group of its own).  A leaf's slots of a kind with an angle (THM)
 *     take their point's angles converted with the leaf's pairs after the grouping (WorkingCopy::InitializeLeaves); a
 *     point whose converted angles differ from its group's writer's gets slots of its own then
 *     (WorkingCopy::SeparateSlotsByAngle);
 *   - the point: no Legendre order (a coefficient reader picks its order from the site's coefficients, Site::ValueFor);
 *     the row (Site::Row is -1 on a slot, nothing reads it there); the target effect (stopping table, screening U_e and its
 *     parameter, the Q coefficients of the Legendre rows) and the c.m./lab factor of the stopping -- in the key; the
 *     target thickness and its lab energy are the point's, never set on a slot;
 *   - a result: the record (the cross sections, A_y, A_0, coefficients) is the site's, written once by the energy
 *     pass; a fold reads it and keeps any other combination of it (E1/E2, A_y sigma, a_k A_0) in arrays of its own
 *     (FoldValues), never in the record.
 *
 * Who writes a slot: the first point, in the order of the data, whose range holds it (WorkingCopy::WritesSlot): the
 * set-up's values (Legendre rows, channel rows, Coulomb and external-capture amplitudes, screening, stopping), the
 * energy pass (CalculationPlan::TaskUnits gives a plan task the slots its point writes), the result layout (one record
 * per site).  Every point of the group reads it: the fold, the model pass, the gradient (a site's cotangent is the sum
 * of its readers', summed in their plan order: UnitCotangent::reader).  A point off the shared mesh (THM, a uniform
 * grid with useAdaptiveGrid off) is a group of its own.  A point with an energy shift reads its planned range (fixed
 * for the stage), so who writes a slot does not depend on the shift (ShiftedFold.h, question 4).
 */
namespace azure2 {

/// One point of a block: its unit in the working copy (WorkingCopy::UnitAt) and its slots, [first, first + count) of
/// the block's sites.
struct SlotMember {
  int unit = -1;
  int first = 0;
  int count = 0;
};

/// The sites of one sharing group in node order and the points that read them.
struct SlotBlock {
  std::vector<Site> sites;
  std::vector<SlotMember> members;
};

/// What points must have in common to share their slots (see above), with the mesh of their nucleus.
struct SlotGroupKey {
  /// The observable the sites compute (SameObservable): two terms of the same reaction share, whatever their segment.
  const TermSettings *settings = nullptr;
  int mesh = -1;
  int targetEffect = 0;  ///< the line of the target (TargetEffect::GetLineNumber), shared by its copies
  std::uint64_t conversion = 0, frameShift = 0, screeningUe = 0;
  int screeningParam = 0;
  bool operator==(const SlotGroupKey &o) const {
    return SameObservable(*settings, *o.settings) && mesh == o.mesh && targetEffect == o.targetEffect &&
           conversion == o.conversion &&
           frameShift == o.frameShift && screeningUe == o.screeningUe && screeningParam == o.screeningParam;
  }
};

}  // namespace azure2

#endif
