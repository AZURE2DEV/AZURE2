#ifndef AZURE2_SHIFTEDFOLD_H
#define AZURE2_SHIFTEDFOLD_H

#include <vector>

#include "MeshPlanner.h"
#include "Unit.h"

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
 *    - the c.m. angle of a lab-frame differential point, which the point converts at its own shifted energy
 *      (Segment::ConvertRowToCM).  The sites have no angle (SiteObservable.h): the point's row set carries it (its
 *      Legendre rows, Coulomb amplitudes and A_y factors, a set of the point's own once the angle has moved:
 *      Unit::CalcNodeRows, RowSet.h).
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
 *    implementation of each fold); no node is recomputed and no level matrix solved.  Kinds of shift enter through
 *    Delta_i(nu) only.  An analytic dK/dE0 is not within reach of this branch: it would need every kernel's energy
 *    arithmetic (its centre in sigma(E0), the clip limits, the Gaussian, the depth window) on a dual type, the kernels'
 *    folds beyond their centre; open for every kernel.  A lab-frame differential segment is no exception: its node
 *    values have no angle, the moved angle is its readers' rows'.
 *
 * A support narrower than MeshZeroWidth() (a layer losing less than 1 eV) gets a single node from the planner, and
 * the fold of a single node is the kernels' thin limit (Kernel::thinLimit).  A shifted point with a reach has a
 * planned support of width 2 R_i, hence more nodes, and its fold reads at least two of them; a fixed shift of such a
 * point keeps the single node and the thin limit.
 *
 * Not here: THM (off the mesh, its sub-points still move rigidly: Segment::MoveOwnGridSlots, branch 10); a point
 * with a uniform grid of its own (useAdaptiveGrid off), which moves rigidly as well; unfolded shifted points (an exact
 * evaluation energy that moves with nu, section 2.5).
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

/// The stopping the target kernel reads at slots [first, first + count) of a shifted point: the target's stopping at
/// each node's energy minus the point's c.m. shift (question 3), held at the nominal support's lower edge for an outer
/// node below it, c.m. MeV cm^2 per active nucleus.
void TranslatedStopping(const Unit &point, const WorkingCopy &data, const Config &configure, const ShiftedNodes &nodes,
                        std::vector<double> &stopping);

/// Whether the node values of `segment`'s folded points do not depend on its shift (question 5): every point of every
/// term folded on the shared mesh, and not a lab-frame differential segment.
bool FoldOnlyShift(const WorkingCopy &data, const Segment &segment);
/// The central difference of a shift on the fixed nodes (question 5): if FoldOnlyShift holds for every one of
/// `segments` (their points moved to a new shift), each of their folded points is folded again from the node values in
/// `results`, and true is returned (their model values are then read with LeafValues::Ready); false, with nothing
/// done, otherwise (the caller calculates them).
bool RefoldShiftedSegments(WorkingCopy &data, const std::vector<const Segment *> &segments, const Config &configure,
                           Results &results);

}  // namespace azure2

#endif  // AZURE2_SHIFTEDFOLD_H
