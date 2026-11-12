#include "azure2/physics/transfer/Overlap.h"

#include <stdexcept>

namespace azure2 {
namespace physics {
namespace transfer {

const std::vector<OverlapKindRow> &OverlapKinds() {
  static const std::vector<OverlapKindRow> rows = {
      {OverlapKind::BoundState, "bound state", "fm^-1/2", "mb/sr"},
      {OverlapKind::Bin, "bin", "fm^-1/2", "mb/sr"},
      {OverlapKind::ScatteringState, "scattering state, delta(E - E')", "MeV^-1/2 fm^-1/2", "mb/(sr MeV)"},
  };
  return rows;
}

const OverlapKindRow &RowOf(OverlapKind kind) {
  for (const OverlapKindRow &row : OverlapKinds())
    if (row.kind == kind) return row;
  throw std::invalid_argument("overlap: unknown kind");
}

RadialOverlap RadialOverlapOf(const BoundState &s) {
  return {s.l, s.two_j, s.step_fm, s.u, s.potential_MeV, s.coulomb_MeV};
}

RadialOverlap RadialOverlapOf(const ContinuumBin &b) {
  return {b.l, b.two_j, b.step_fm, b.u, b.potential_MeV, b.coulomb_MeV};
}

RadialOverlap RadialOverlapOf(const ScatteringState &s, double mu_u) {
  RadialOverlap o{s.l, s.two_j, s.step_fm, s.u, s.potential_MeV, s.coulomb_MeV};
  const double c = NormalisationFactor(s, mu_u, ContinuumNormalisation::DeltaE);
  for (double &x : o.u) x *= c;
  return o;
}

RadialOverlap Truncated(const RadialOverlap &o, double r_cut, RadialPart part) {
  RadialOverlap t = o;
  for (std::size_t i = 0; i < t.u.size(); ++i) {
    const double r = i * t.step_fm;
    const bool keep = part == RadialPart::Interior ? r <= r_cut + 1e-9 : r > r_cut + 1e-9;
    if (!keep) t.u[i] = 0.0;
  }
  return t;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
