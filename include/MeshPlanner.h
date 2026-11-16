#ifndef AZURE2_MESHPLANNER_H
#define AZURE2_MESHPLANNER_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {

/*!
 * The planner of the shared energy mesh of one compound nucleus (docs/upgrade/energy_grid.md §2.2 "Then", §2.8): a
 * function of values only.  In: the quantized anchors of the nucleus's levels, one support per folded point (the
 * interval of E_x where its kernel is non-zero, with its data set's numerical settings), the tail constants of
 * AdaptiveIntegrationGrid.  Out: the nodes (E_x, MeV, ascending, distinct) and, per support, the range of nodes its
 * kernel's quadrature uses.  Nothing is wired into the engine yet.
 *
 * Today each folded point has a private grid (AdaptiveIntegrationGrid::GenerateGrid): a uniform smooth lattice of
 * step S = window / ceil(window / baseEnergyStep) starting at the window's edge, and for every anchor finer than S its
 * own lattice.  Its local step at E, the "step of today" against which the mesh is measured, is
 *     D(E) = min(S, min over anchors r of l_r(E)),
 *     l_r(E) = Gamma_r / p                         for |E - E_r| <= m Gamma_r,
 *            = max(Gamma_r / p, t |E - E_r|)       beyond,            t = tailRatio / p * sqrt(max(c, m) / c),
 * with p points per width, m the core half-extent in widths, c = convolutionCoreWidths (Gamma_r: the anchor's particle
 * width).  The planner places nodes so that the mesh interval containing E is never longer than D(E) of any support
 * covering E (the bound of §2.4 is then never worse for any point), from lattices that are functions of the anchors
 * alone, so that every support asking for a node gets the same bits.  The rule for each kind of node:
 *
 *  - Smooth nodes (the ladder).  A support's smooth step S is quantized downward to a power of two,
 *    L = 2^floor(log2 S) MeV (exact: frexp), and its smooth nodes are n * L for the integers n from floor(low / L) to
 *    ceil(high / L): multiples of L from E_x = 0, exact in floating point (n < 2^53), the same bits whoever asks.  The
 *    ladders are nested (2^-j is a multiple of 2^-(j+1)), so where supports of different steps overlap the finer one
 *    simply has more nodes and every change of step falls on a node of both.  Cost: up to 2x the smooth nodes of a
 *    support (1.44x for steps spread evenly in log; 1.15-1.80 on the folded projects of the tree, §2.8).
 *  - Core nodes.  For an anchor (E_r, Gamma_r) and the settings (p, m) of a class of supports: E_r + k * pitch and
 *    E_r - k * pitch, pitch = Gamma_r / p, k = 0 ... nCore with nCore = ceil(m p - 1e-9) (today's floor(m p + 1e-9)
 *    when m p is an integer; one node more otherwise, so that the core reaches m widths).
 *  - Tail nodes.  E_r + o_k and E_r - o_k with o_0 = nCore * pitch and o_{k+1} = o_k + max(pitch, t * o_k), the same
 *    recurrence as today's; a function of the anchor and the class alone, never of a window.  The recurrence runs while
 *    the step is below twice the largest ladder step L of the class's supports (today's tails run up to the window's
 *    smooth step S, and S < 2 L).
 *  - Overlaps and transitions.  Every lattice (one per anchor and class, one per ladder step) is a list of intervals,
 *    each with a length lambda (pitch, the tail step max(pitch, t o_k), or L).  An interval of a lattice X is kept, both
 *    its ends becoming nodes, if X is active somewhere in it (a support of X's class, or of X's ladder step, covers a
 *    stretch of it) and no single other lattice Y dominates it: Y active over all of the part of it inside the
 *    supports and finer over the whole interval than lambda (two anchors' lattices, two ladders) or than lambda / 2
 *    (a ladder against an anchor's lattice, either way round).  The half mirrors today's grid, which keeps its smooth
 *    nodes wherever no lattice is finer than half the smooth step and a tail up to the smooth step, so that in the band
 *    between it has both; without it the mesh had up to a third fewer nodes than a point's own grid inside some
 *    supports of natural_neon.  So inside a level's reach its finer lattice wins; where two lattices cross, both keep
 *    the interval that straddles the crossing (a few nodes more than one lattice would need); the coarser one resumes
 *    at its first interval that the finer one no longer dominates, at one of its own nodes.
 *
 * Why this keeps "never coarser than today": take E in a support and a lattice F whose step attains the minimum h(E)
 * of all lattices active at E.  F's interval containing E has length lambda <= l_F(E) = h(E) (a tail interval is as
 * long as the step at its inner end); a Y dominating it would be finer than h(E) at E, so the interval is kept, and
 * the mesh interval containing E is no longer (the half in the rule only makes domination harder).  Since L <= S and the
 * anchor lattices extend at least as far as today's, h(E) <= D(E).  A soft edge of a support adds no node: the
 * quadrature clips the interval containing it, whose outer node lies at most one interval outside the support.  A
 * hard edge (MeshSupport::lowHard, highHard: a physical floor) is a node, and the support's range ends at it.
 *
 * Determinism: each node is computed by the formula of its kind from the anchor or n alone (one multiplication and one
 * addition, or the tail recurrence from o_0); which intervals are kept depends on the set of supports and anchors, not
 * on their order; the nodes are the sorted distinct values.  Nothing is accumulated across supports.
 *
 * Locality: a lattice can change the kept intervals only inside its own reach (the span of its generated nodes), so
 * nodes move only within an anchor's reach plus one interval of the longest kept length, at most twice the largest
 * ladder step (MeshChangedRegion).
 */

/// A quantized level anchor, in E_x (MeV): AdaptiveIntegrationGrid::Quantize applied in the frame of E_x.
struct MeshAnchor {
  double energy = 0.0;         ///< E_x of the anchor, MeV
  double totalWidth = 0.0;     ///< Gamma (particle + radiative), MeV; not used for the lattice
  double particleWidth = 0.0;  ///< Gamma_particle, MeV: the lattice's width (as today)
  int id = 0;                  ///< 1000 * J group + level, for the caller
};

/// The numerical settings of a data set that shape its lattices (format 3 `quadrature`).
struct MeshSettings {
  double pointsPerWidth = 50.0;  ///< p: core nodes per width; also sets the tail's growth
  double coreWidths = 20.0;      ///< m: half-extent of the uniform core, in widths
};

/// A support narrower than this (MeV) is a single exact node at its high edge, as today's grid of such a window is.
constexpr double MeshZeroWidth() { return 1.0e-6; }

/// The same for a time-of-flight window whose high edge is at E_x = high: 4096 spacings of the doubles there
/// (7e-12 MeV at E_x = 8 MeV).  Such a window is a fixed fraction of the point's energy (twice the relative width of
/// the delay distribution), so at a thermal neutron point it is far below MeshZeroWidth() and still wide against the
/// structure of the cross section; with the absolute width every such point was a single node and the fold returned
/// the cross section at the window's high edge (a bias of half the channel's relative width on 1/v).  Below 4096
/// spacings a ladder cannot be laid (its step is floored at one spacing, MeshLadderStep), and the single node is then
/// within 3e-4 of a 25 meV point's energy.
inline double MeshTimeOfFlightZeroWidth(double high) {
  const double x = std::fabs(high);
  return 4096.0 * (std::nextafter(x, std::numeric_limits<double>::infinity()) - x);
}

/// The support of one folded point's kernel, in E_x (MeV), with its data set's settings.
struct MeshSupport {
  double low = 0.0, high = 0.0;  ///< the kernel is non-zero on [low, high]
  double smoothStep = 0.0;       ///< the step asked for away from levels, MeV (today window / number of sub-points)
  MeshSettings settings;
  /// A hard edge lies at a physical floor (a thick target's lower limit, the entrance threshold, a THM floor): it is a
  /// node of the mesh and the support's range starts (ends) at it, so no node beyond the floor is evaluated for this
  /// support.  A soft edge adds no node: the quadrature clips the interval containing it.
  bool lowHard = false, highHard = false;
  int dataSet = 0;  ///< index of the data set, for the caller (the plan does not depend on it)
  std::string label;        ///< "data set NAME, point I (E_lab MeV)", for a message (the plan does not depend on it)
  double kernelWidth = 0.0; ///< MeV, the width of the point's own kernel support (unshifted), for a message
  /// A support narrower than this (MeV) is a single exact node at its high edge: MeshZeroWidth(), or
  /// MeshTimeOfFlightZeroWidth(high) for a time-of-flight window.
  double zeroWidth = MeshZeroWidth();
};

/// The most nodes the plan of one compound nucleus may have (MeshPlanBudgetError beyond).  The largest plan of the
/// tree's projects is natural_neon's 21,070 nodes (docs/upgrade/energy_grid.md, the table of section 2.8); a node of
/// the point pass takes a few kB (a 1e14 atoms/cm^2 layer with a 3 keV shift prior asked for about 1e7 nodes and took
/// 28 GB before it was stopped), so 1e6 nodes is about 50 times the largest and a few GB.
constexpr double kMeshNodeBudget = 1.0e6;

/// A plan that would exceed kMeshNodeBudget, refused before its nodes are made (what() names the support that asks for
/// the most nodes and what to do).
struct MeshPlanBudgetError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

/// The constants of the geometric tails (AdaptiveIntegrationGrid::tailRatio, convolutionCoreWidths).
struct MeshTailConstants {
  double tailRatio = 0.0;
  double convolutionCoreWidths = 0.0;
};

/// The nodes a support's quadrature uses: nodes[first] <= low and nodes[last] >= high (the outer nodes of the
/// intervals that contain the edges; a hard edge is itself the end node), or first == last for a support narrower
/// than MeshZeroWidth().
struct MeshNodeRange {
  std::size_t first = 0, last = 0;
};

struct MeshPlan {
  std::vector<double> nodes;            ///< E_x, MeV, ascending and distinct
  std::vector<MeshNodeRange> ranges;    ///< per support, in the order given
};

/// The smooth step a support's ladder uses: today's step S = window / ceil(window / smoothStep) quantized down to a
/// power of two (MeV); 0 for a support narrower than MeshZeroWidth().
double MeshLadderStep(const MeshSupport &support);

/// The tail growth t of a class of settings (step per distance from the anchor).
double MeshTailFactor(const MeshSettings &settings, const MeshTailConstants &tails);

/// Plans the mesh of one compound nucleus.  Anchors with no positive particle width are ignored.  Throws
/// MeshPlanBudgetError when the supports' ladders alone ask for more than kMeshNodeBudget nodes.
MeshPlan PlanMesh(const std::vector<MeshAnchor> &anchors, const std::vector<MeshSupport> &supports,
                  const MeshTailConstants &tails);

/// The stretches of E_x (merged, ascending) outside which PlanMesh(before, supports) and PlanMesh(after, supports)
/// have the same nodes, bit for bit: per anchor present in one list and not the other, its reach for every class of
/// settings, widened by the longest interval the plan can keep (twice the largest ladder step).  What the hand-over
/// re-places when quantized anchors change (decision D1).
std::vector<std::pair<double, double>> MeshChangedRegion(const std::vector<MeshAnchor> &before,
                                                         const std::vector<MeshAnchor> &after,
                                                         const std::vector<MeshSupport> &supports,
                                                         const MeshTailConstants &tails);

}  // namespace azure2

#endif  // AZURE2_MESHPLANNER_H
