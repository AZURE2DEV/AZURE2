#ifndef AZURE2_PHYSICS_TRANSFER_OVERLAP_H
#define AZURE2_PHYSICS_TRANSFER_OVERLAP_H

#include <vector>

#include "azure2/physics/transfer/BoundState.h"
#include "azure2/physics/transfer/Continuum.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The overlap functions of a transfer (docs/upgrade/transfer.md sections 6 and 15).  The kernel of Eq. (6.4) needs
 * of an overlap only its radial function u(r) and the potential V(r) of x relative to the core on a grid; what the
 * function is (a bound state, a bin, a scattering state at one energy) changes its unit and so the unit of the
 * cross section, nothing else.  The kinds are the rows of OverlapKinds(); a new kind (a Gamow state, a
 * pseudo-state) is one row and one RadialOverlapOf.
 */

/// u (and the potential that binds or scatters x) on r_i = i h: all the kernel reads.
struct RadialOverlap {
  int l = 0, two_j = 0;
  double step_fm = 0.0;
  std::vector<double> u;              ///< u(r_i), in the unit of the kind
  std::vector<double> potential_MeV;  ///< V(r_i), Coulomb included (the post- or prior-form binding interaction)
  std::vector<double> coulomb_MeV;    ///< its Coulomb part
  double r_min_fm = 0.0;              ///< u = 0 below it (a cut analytic tail): radial integrals start at that node
};
/// The first node of the radial integrals of an overlap: the node of r_min (0 without a cut).
int FirstNode(const RadialOverlap &overlap);

enum class OverlapKind { BoundState, Bin, ScatteringState };

/// One row: the kind, its name, the unit of u and of the cross section to it.
struct OverlapKindRow {
  OverlapKind kind;
  const char *name;
  const char *u_unit;
  const char *cross_section_unit;
};
/// The table of kinds (constant).
const std::vector<OverlapKindRow> &OverlapKinds();
const OverlapKindRow &RowOf(OverlapKind kind);

/// The radial overlap of a bound state (fm^-1/2), of a bin (fm^-1/2, Eq. (14.3)) and of a scattering state in the
/// delta(E - E') normalisation of Eq. (14.2) (MeV^-1/2 fm^-1/2: the transfer to it is d2sigma/dOmega dE).
RadialOverlap RadialOverlapOf(const BoundState &state);
RadialOverlap RadialOverlapOf(const ContinuumBin &bin);
RadialOverlap RadialOverlapOf(const ScatteringState &state, double mu_u);

/// u set to 0 beyond r_cut (fm), or up to it: the interior and exterior parts of section 16 (V is kept).
enum class RadialPart { Interior, Exterior };
RadialOverlap Truncated(const RadialOverlap &overlap, double r_cut_fm, RadialPart part);

/// One overlap of the reaction: kind, radial function, signed spectroscopic amplitude, spins.
struct Overlap {
  OverlapKind kind = OverlapKind::BoundState;
  RadialOverlap radial;
  double amplitude = 1.0;      ///< A (signed)
  int two_J_core = 0;          ///< 2J of the core (s for the projectile, A for the residual)
  int two_J_composite = 0;     ///< 2J of the composite (a, F; of the x + A state for a bin or a scattering state)
};

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
