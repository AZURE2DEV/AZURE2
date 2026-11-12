#ifndef AZURE2_PHYSICS_TRANSFER_INTERACTION_H
#define AZURE2_PHYSICS_TRANSFER_INTERACTION_H

#include <complex>
#include <vector>

#include "azure2/physics/transfer/OpticalPotential.h"
#include "azure2/physics/transfer/Overlap.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The interaction of the transfer A(a, s)F, a = s + x, F = A + x, as a sum of terms (docs/upgrade/transfer.md
 * section 17).  With H = T + V_sx + V_xA + V_sA the post and prior forms are
 *
 *   V_post  = V_sx(r_sx) + [V_sA(r_sA) - U_sF(r_sF)],   V_prior = V_xA(r_xA) + [V_sA(r_sA) - U_aA(r_aA)],   (17.1)
 *
 * the brackets the remnants.  The binding terms are the potentials of the overlaps (their tables, Coulomb included
 * or not), the core-core and optical terms are potentials of Eq. (2.1) given here (complex: FRESCO's IP2 = -1).
 * Each kind is one row of InteractionKinds(): which coordinate it depends on and its sign.  A new term (a tensor
 * force, a spin-orbit core-core term) is one row, one case in TermValue and a test.
 */
enum class InteractionKind {
  ProjectileBinding,         ///< V_sx(r_sx): the projectile overlap's potential (post form)
  ProjectileBindingNuclear,  ///< V_sx - V_C (FRESCO with a KRPOT that cancels the Coulomb term)
  ResidualBinding,           ///< V_xA(r_xA): the residual overlap's potential (prior form)
  CoreCore,                  ///< + V_sA(r_sA), r_sA = r_sx + r_xA (the remnants)
  ExitOptical,               ///< - U_sF(r_sF) (the post-form remnant)
  EntranceOptical            ///< - U_aA(r_aA) (the prior-form remnant)
};

/// The coordinate a term depends on (Eq. (5.1): r_aA = r + beta rho, r_sF = alpha r + rho, r_sA = r + rho).
enum class InteractionCoordinate { SX, XA, SA, SF, AA };

struct InteractionKindRow {
  InteractionKind kind;
  const char *name;
  InteractionCoordinate coordinate;
  double sign;  ///< +1 or -1 in Eq. (17.1)
};
const std::vector<InteractionKindRow> &InteractionKinds();
const InteractionKindRow &RowOf(InteractionKind kind);

/// One term: its kind and, for the core-core and optical terms, the potential and the charge product of the pair.
struct InteractionTerm {
  InteractionKind kind = InteractionKind::ProjectileBinding;
  OpticalPotential potential;  ///< CoreCore, ExitOptical, EntranceOptical
  int z1z2 = 0;                ///< the Coulomb charge product of that pair
};

/// The post form without remnant (part 1's default), with remnant, and the prior forms.
std::vector<InteractionTerm> PostForm();
std::vector<InteractionTerm> PostFormWithRemnant(const OpticalPotential &core_core, int z_s_z_A,
                                                 const OpticalPotential &exit, int z_s_z_F);
std::vector<InteractionTerm> PriorForm();
std::vector<InteractionTerm> PriorFormWithRemnant(const OpticalPotential &core_core, int z_s_z_A,
                                                  const OpticalPotential &entrance, int z_a_z_A);

/// The value (MeV) of a term that is not a binding term at the distances (fm) of one kernel node.
struct NodeDistances {
  double sx = 0.0, xa = 0.0, sa = 0.0, sf = 0.0, aa = 0.0;
};
std::complex<double> TermValue(const InteractionTerm &term, const NodeDistances &d);

/// V(r) u(r) of a binding term on the overlap's grid (MeV fm^-1/2), the Coulomb part dropped for the Nuclear kind.
std::vector<double> BindingSource(const InteractionTerm &term, const RadialOverlap &overlap);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
