#include "azure2/physics/transfer/Interaction.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace azure2 {
namespace physics {
namespace transfer {

const std::vector<InteractionKindRow> &InteractionKinds() {
  static const std::vector<InteractionKindRow> rows = {
      {InteractionKind::ProjectileBinding, "V_sx", InteractionCoordinate::SX, 1.0},
      {InteractionKind::ProjectileBindingNuclear, "V_sx - V_C", InteractionCoordinate::SX, 1.0},
      {InteractionKind::ResidualBinding, "V_xA", InteractionCoordinate::XA, 1.0},
      {InteractionKind::CoreCore, "V_sA", InteractionCoordinate::SA, 1.0},
      {InteractionKind::ExitOptical, "-U_sF", InteractionCoordinate::SF, -1.0},
      {InteractionKind::EntranceOptical, "-U_aA", InteractionCoordinate::AA, -1.0},
  };
  return rows;
}

const InteractionKindRow &RowOf(InteractionKind kind) {
  for (const InteractionKindRow &row : InteractionKinds())
    if (row.kind == kind) return row;
  throw std::invalid_argument("interaction: unknown kind");
}

std::vector<InteractionTerm> PostForm() { return {InteractionTerm{InteractionKind::ProjectileBinding, {}, 0}}; }

std::vector<InteractionTerm> PostFormWithRemnant(const OpticalPotential &core_core, int z_sA,
                                                 const OpticalPotential &exit, int z_sF) {
  return {InteractionTerm{InteractionKind::ProjectileBinding, {}, 0},
          InteractionTerm{InteractionKind::CoreCore, core_core, z_sA},
          InteractionTerm{InteractionKind::ExitOptical, exit, z_sF}};
}

std::vector<InteractionTerm> PriorForm() { return {InteractionTerm{InteractionKind::ResidualBinding, {}, 0}}; }

std::vector<InteractionTerm> PriorFormWithRemnant(const OpticalPotential &core_core, int z_sA,
                                                  const OpticalPotential &entrance, int z_aA) {
  return {InteractionTerm{InteractionKind::ResidualBinding, {}, 0},
          InteractionTerm{InteractionKind::CoreCore, core_core, z_sA},
          InteractionTerm{InteractionKind::EntranceOptical, entrance, z_aA}};
}

std::complex<double> TermValue(const InteractionTerm &t, const NodeDistances &d) {
  const InteractionKindRow &row = RowOf(t.kind);
  double r = 0.0;
  switch (row.coordinate) {
    case InteractionCoordinate::SA: r = d.sa; break;
    case InteractionCoordinate::SF: r = d.sf; break;
    case InteractionCoordinate::AA: r = d.aa; break;
    default: throw std::invalid_argument(std::string("interaction: ") + row.name + " is a binding term");
  }
  return row.sign * PotentialAt(t.potential, t.z1z2, 0.0, std::max(r, 1e-6));
}

std::vector<double> BindingSource(const InteractionTerm &t, const RadialOverlap &o) {
  if (t.kind != InteractionKind::ProjectileBinding && t.kind != InteractionKind::ProjectileBindingNuclear &&
      t.kind != InteractionKind::ResidualBinding)
    throw std::invalid_argument("interaction: not a binding term");
  if (o.potential_MeV.size() != o.u.size()) throw std::invalid_argument("interaction: the overlap has no potential");
  const bool nuclear = t.kind == InteractionKind::ProjectileBindingNuclear;
  std::vector<double> vu(o.u.size());
  for (std::size_t k = 0; k < vu.size(); ++k)
    vu[k] = (o.potential_MeV[k] - (nuclear ? o.coulomb_MeV[k] : 0.0)) * o.u[k];
  return vu;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
