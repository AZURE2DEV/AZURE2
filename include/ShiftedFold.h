#ifndef AZURE2_SHIFTEDFOLD_H
#define AZURE2_SHIFTEDFOLD_H

#include <functional>
#include <vector>

#include "MeshPlanner.h"
#include "Unit.h"

class CNuc;
class Config;
class Segment;
class WorkingCopy;
struct SegmentState;
namespace azure2 {
struct Results;
}

/*!
 * Folded data with an energy shift, on the fixed nodes of the shared mesh (docs/upgrade/energy_grid.md section 2.5).
 *
 * A point of a data set with an energy shift nu (a constant a, b sqrt(E), a shift function, a follower of a shift
 * leader; Segment::TotalEnergyShift) is measured at the nominal lab energy E_lab and modelled at E_lab + s(nu; E_lab).
 * Its c.m. shift is Delta_i(nu) = E_cm,i(nu) - E_cm,i(0).  With a fold, the nodes stay and the kernel moves:
 *
 *   m_i(nu) = sum_n w_n(Delta_i) sigma(E_n),   the kernel K_i(E; nu) = K_i(E - Delta_i; 0) on the fixed nodes E_n,
 *
 * which is the fold of the engine before (the point's sub-points moved rigidly with it) in exact arithmetic: the
 * numbers move by the quadrature difference only, chi^2 is smooth in nu, and a change of nu recomputes no node.  The
 * five questions of the design, answered:
 *
 * 1. What the planner is given (PlannedShiftSupport).  The point's unshifted support [a_i, b_i] (the kernels'
 *    EffectSupport at the nominal point), moved by a centre c_i and widened by a reach R_i on each side:
 *    [a_i + c_i - R_i, b_i + c_i + R_i]; its low edge never goes below the integration floor (it is a hard edge there).
 *    R_i = 3 sum_k |Delta_i(nu_0 + sigma_k e_k) - Delta_i(nu_0)| over the shift parameters k that a fit varies (the
 *    leader's for a shift-group follower, the owner's for a shift-function group), sigma_k the parameter's prior
 *    width; a varied parameter without a prior adds the support's own width b_i - a_i (ShiftReach).  R_i = 0 when no
 *    shift parameter of the point varies: a fixed shift moves the support once and adds no node.  The planner's nodes
 *    are functions of the anchors, the supports and the settings (MeshPlanner.h, "Determinism"): a wider support adds
 *    nodes and moves none of the narrower one's except where an interval's "kept" decision changes; the cost in nodes
 *    on the shift projects is measured in energy_grid.md section 2.5.  The reach comes from what the project says
 *    (the priors), no new field; a new kind of shift enters only through Delta_i.
 *
 * 2. When nu leaves the reach.  The centre is a function of nu alone, c_i = Delta_i(nu_0) + 2 R_i
 *    round((Delta_i(nu) - Delta_i(nu_0)) / (2 R_i)) (R_i = 0: c_i = Delta_i(nu)) (ShiftCentre), so the moved support
 *    always lies inside the planned one.  When some c_i changes, the hand-over plans the nucleus again with the copy's
 *    supports and re-places the sharing groups whose nodes changed, through the path of a moved narrow level
 *    (WorkingCopy::ReplaceSubPointGrids, decision D1): no second code path, a function of (theta, nu) only, so the
 *    same at any thread count; counted in run.log.  chi^2 steps by the quadrature difference there, as at a moved
 *    level.  Not a refusal: a fit that walks far is not an error.
 *
 * 3. Which sites a point reads at nu (ShiftedFoldNodes).  The point's slots are its planned range (fixed for the
 *    stage); per evaluation the fold reads the tightest range of them enclosing the moved support [a_i + Delta_i,
 *    b_i + Delta_i] (at least two nodes) and the kernel integrates over the moved support alone (KernelNodes::
 *    supportLow, supportHigh, the clip of section 2.4).  The block is not rewritten.  What still depends on nu, and
 *    where it lives:
 *    - the stopping the target kernel reads.  The fold before evaluated each sub-point's stopping at its unshifted
 *      energy (the rigid move kept it): the depth profile, the target thickness T and 1/epsilon are the point's as
 *      measured, translated with the beam energy.  The gather gives the kernel epsilon(E_n - Delta_i) at each node
 *      (TranslatedStopping), from the target's stopping and the copy's stopping scale: the point's, not the site's,
 *      which keeps epsilon(E_n) for every unshifted reader.  At Delta_i = 0 it is the site's value bit for bit.
 *    - the Jacobian of a point given at a lab angle, which the point converts at its own shifted energy
 *      (Segment::ConvertRowToCM) and multiplies its lab fold by (Unit::LabFoldFactor).  Its c.m. angle at a node is
 *      not the point's: it is theta_cm(E_n; theta_lab), a function of the node's energy and the lab angle, held by the
 *      row set shared by every reader of that lab angle (RowSet.h, energy_grid.md section 4.3), so the shift does not
 *      move it.
 *    - nothing else: the node's energies, channel rows, Q coefficients (a function of the node energy), screening and
 *      external-capture amplitudes are the node's.  The kernel's own nu dependence (its centre E0 = the point's
 *      shifted c.m. energy; a beam-spread formula sigma(E) at the shifted depth energies) is in the kernel's input,
 *      as before.  The kernels know nothing of shifts: they get nodes, values, stopping, a support and a centre.
 *
 * 4. Who writes a slot.  A shifted point's slot range is its planned range, fixed for the stage, so the block (the
 *    union of the planned ranges), the writer of each slot (WorkingCopy::WritesSlot: first reader in data order), the
 *    plan's tasks and the energy pass's task list do not depend on nu.  Every site of the block is computed at every
 *    evaluation: the margin's nodes are the price (the director's proposal).
 *
 * 5. The gradient.  The column of a shift parameter of a folded segment is the central difference of the fold on
 *    the fixed nodes (RefoldShiftedSegments): at nu +- h the points move (their centre, support, translated stopping),
 *    and each point is folded again from the node values already computed at nu with the same fold (FoldValues: one
 *    implementation of each fold); no node is recomputed and no level matrix solved.  The column does not plan again,
 *    so its step keeps every moved support inside the planned nodes (InsideShiftStep: the step halved; failing that,
 *    the caller's path that plans for each side): a reach below the step, or a shift near the edge of its window, made
 *    ShiftedFoldNodes throw inside the parallel refold before.  Kinds of shift enter through
 *    Delta_i(nu) only.  An analytic dK/dE0 is not within reach of this branch: it would need every kernel's energy
 *    arithmetic (its centre in sigma(E0), the clip limits, the Gaussian, the depth window) on a dual type, the kernels'
 *    folds beyond their centre; open for every kernel.  A segment given at a lab angle is no exception: its node
 *    values and its row sets do not depend on the shift, only its points' own Jacobians.
 *
 * A support narrower than MeshZeroWidth() (a layer losing less than 1 eV) gets a single node from the planner, its
 * upper edge, and the fold of a single node is the kernels' thin limit (Kernel::thinLimit).  The same holds for a
 * shifted point whose planned support [a_i + c_i - R_i, b_i + c_i + R_i] is narrower than MeshZeroWidth(): a fixed
 * shift (R_i = 0, the node at the moved support's upper edge b_i + Delta_i: the unshifted thin limit at the shifted
 * energy), or a varied shift without a prior on a layer below about 1/7 eV (R_i = 3 (b_i - a_i)).  Its fold is the thin
 * limit at the node b_i + c_i + R_i whatever nu: the moved support lies inside the planned one, so the node is at most
 * the planned support's width (below 1 eV) above the moved support, and the fold is off by at most sigma' x 1 eV
 * (ShiftedFoldNodes refuses a moved support that reaches farther than MeshZeroWidth() below the node).  Within one
 * reach the fold does not change with nu; at a re-plan (nearly every step of a fit, the reach being a fraction of an
 * eV) it steps by sigma' x 2 R_i, so chi^2(nu) is a staircase of that size, which a numerical second derivative of
 * MIGRAD sees (tests/thin_layer_limit; the follow-ups' report has the numbers).  A shifted point with a wider planned
 * support has more nodes, and its fold reads at least two of them.
 *
 * THM data with a shift are this case too.  A THM shift is not floored at zero energy (a HOES excitation function
 * reaches below the entrance threshold, THM.md): Delta_i is the point's own c.m. shift, which the shift step leaves
 * unfloored for such a data set, and the floor of the moved support is the data kind's (spec.floorEnergy, from
 * Segment::IntegrationFloor: E_x = minIntegrationEnergy of the compound system), a hard edge of the planned support as
 * a thick target's is.  The reach widens the support; the kernel moves on fixed nodes; nothing else is THM's.
 *
 * Not here: a point with a uniform grid of its own (useAdaptiveGrid off), which moves rigidly with its shift
 * (Segment::MoveOwnGridSlots) and reads a row set of its own (Unit::CalcNodeRows), for every kind of data; unfolded
 * shifted points (an exact evaluation energy that moves with nu, section 2.5).  The uniform grid is the last user of
 * the rigid move.  Proposal (not decided here: useAdaptiveGrid is a user-visible key of <config>): its place on the
 * mesh is a setting of the supports, MeshSettings::uniform -- the support's nodes the uniform lattice of `points`
 * intervals from its high edge, no anchors -- planned like any other support, so the point reads fixed nodes and the
 * rigid move, CalcNodeRows and the overflow nodes of moved sub-points go; the alternative is the key's removal with a
 * note of the converter (the adaptive mesh is never coarser than the uniform grid of the same `points`, section 2.4).
 * Keeping the own grid as it is costs MoveOwnGridSlots, Unit::CalcNodeRows and RowSet::valueOffset -1 (a set of a
 * reader's own), the own-grid branch of SetUpTargetEffects, and the overflow placement of moved slots.
 */
namespace azure2 {

/// Whether a data set's energy shift is active at the set-up (a varied, non-zero, group, sqrt or function shift).
bool HasEnergyShift(const Segment &segment);

/// Whether a folded point is on the shared mesh with an energy shift: its data set's shift was active at the set-up,
/// or its support has moved since (a shift set later).  Its fold then reads ShiftedFoldNodes.
bool IsShiftedOnMesh(const Unit &point, const WorkingCopy &data);

/// The c.m. shift Delta_i (MeV) of a point of `segment` at nominal lab energy `labEnergy` in the shift state `state`:
/// cmFactor x TotalEnergyShift (the conversion of Segment::CMEnergyFactor).
double PointShift(const Segment &segment, const SegmentState &state, double labEnergy, double cmFactor);

/// The reach R_i (c.m. MeV) of a point of `segment` (question 1): 3 x the sum over the varied shift parameters of
/// the move of its c.m. shift by one prior width, the support's width for a varied parameter without a prior; 0 when
/// no shift parameter of the point varies.  `segments` are the data set's (a follower reads its leader's prior).
double ShiftReach(const std::vector<Segment> &segments, const Segment &segment, const SegmentState &start, double labEnergy,
                  double cmFactor, double supportWidth);

/// The centre c_i (c.m. MeV) of the planned support at c.m. shift `delta` (question 2): a function of delta alone.
double ShiftCentre(double delta, double startDelta, double reach);

/// The support handed to the planner (E_x) for a shifted point of grid record `spec` at centre `centre` (question 1),
/// from the nominal support `nominal` (E_x).
MeshSupport PlannedShiftSupport(const MeshSupport &nominal, const Unit::SubGridSpec &spec, double centre);

/// The nodes a shifted point's fold reads (question 3): slots [first, first + count) of its range (0-based), the
/// moved support [low, high] (c.m. MeV) and its c.m. shift.
struct ShiftedNodes {
  int first = 0, count = 0;
  double low = 0.0, high = 0.0;
  double delta = 0.0;
};
/// For a point on the shared mesh whose data set has an energy shift (Unit::SubGridSpec::shifted); throws if the moved
/// support leaves the point's planned nodes (a planning defect: question 2 keeps it inside).
ShiftedNodes ShiftedFoldNodes(const Unit &point, const WorkingCopy &data);

/// Whether the moved support of every point of `segments` on the shared mesh with an energy shift lies inside its
/// planned nodes, the condition of ShiftedFoldNodes (question 5: a shift column moves the points without a re-plan).
bool ShiftedSupportsInside(const WorkingCopy &data, const std::vector<const Segment *> &segments);
/// The step of a central difference of a shift parameter at x0 that keeps the points of the data sets
/// `segmentNumbers` (1-based) inside their planned nodes on both sides (question 5): h, halved until move(x0 +- step)
/// does; 0 when 40 halvings do not (the caller then takes a path that plans again).  move(v) sets the parameter to v
/// and moves the points; it is left at x0.
double InsideShiftStep(const WorkingCopy &data, const std::vector<int> &segmentNumbers, double x0, double h,
                       const std::function<void(double)> &move);

/// The stopping the target kernel reads at slots [first, first + count) of a shifted point: the target's stopping at
/// each node's energy minus the point's c.m. shift (question 3), held at the nominal support's lower edge for an outer
/// node below it, c.m. MeV cm^2 per active nucleus.
void TranslatedStopping(const Unit &point, const WorkingCopy &data, const Config &configure, const ShiftedNodes &nodes,
                        std::vector<double> &stopping);

/// Whether the node values of `segment`'s folded points do not depend on its shift (question 5): every point of every
/// term folded on the shared mesh.
bool FoldOnlyShift(const WorkingCopy &data, const Segment &segment);
/// The central difference of a shift on the fixed nodes (question 5): if FoldOnlyShift holds for every one of
/// `segments` (their points moved to a new shift), each of their folded points is folded again from the node values in
/// `results`, and true is returned (their model values are then read with LeafValues::Ready); false, with nothing
/// done, otherwise (the caller calculates them).
bool RefoldShiftedSegments(WorkingCopy &data, const std::vector<const Segment *> &segments, const Config &configure,
                           Results &results);

/// The chi^2 of `segment` alone at the shift its points hold now, its data compared at the norm `norm`, with its shift
/// prior: the model of every point (refolded on the fixed nodes where FoldOnlyShift holds, else calculated), the
/// residuals, the penalty.  The local central difference of a shift column (AZURECalc::Gradient).
double SegmentLocalChi2(const Segment &segment, CNuc *compound, const Config &config, WorkingCopy &data,
                        Results &results, double norm);
/// The norm each of `segments` (1-based) is compared at in a local shift difference: the segment's own, or for a
/// profiled norm (THM, Segment::IsProfiledNorm) its profile n* at the models of the present shift, held fixed while the
/// shift moves -- at the profile d chi^2/d n* = 0, so the partial derivative at fixed n* is the total one (the envelope
/// theorem).  Calculates the profiled segments' models at the present parameters.
std::vector<double> HeldNorms(WorkingCopy &data, const std::vector<int> &segments, CNuc *compound, const Config &config,
                              Results &results);

}  // namespace azure2

#endif  // AZURE2_SHIFTEDFOLD_H
