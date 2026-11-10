#ifndef AZURE2_EXPERIMENT_UNITS_H
#define AZURE2_EXPERIMENT_UNITS_H

#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * The unit table of format 3 (docs/upgrade/project_format3.md sections 3.3, 7.2): every unit token a project may
 * write, its dimension, and its factor to the dimension's base unit, which is the engine's unit of that dimension
 * (MeV, deg, b, b/sr, MeV*b, atoms/cm2, MeV*cm2/atom, K, fm, ...).  A value v written in unit u is v * u.factor in
 * the base unit.  A quantity whose engine unit is not its dimension's base (widths in eV) converts with
 * Factor(u) / Factor(engine unit).  There is no default unit: a dimensioned number always carries its token.
 *
 * Per-charge units (1/uC, counts/uC, ...) have the base 1/C; a yield per incident particle Y becomes Y / (q e) per C,
 * i.e. Y * ParticlesPerUnitCharge(u, q) per unit, with e from physics/Constants.h.
 */
struct Unit {
  const char *token;
  const char *dimension;
  double factor;  ///< to the dimension's base unit
  std::vector<std::string> aliases = {};  ///< accepted and written back as `token` (µb -> ub); none by default
};

struct Dimension {
  const char *name;
  const char *base;  ///< the engine unit (factor 1)
  const char *doc;
};

const std::vector<Unit> &Units();
const std::vector<Dimension> &Dimensions();

/// The unit of a token or alias, or nullptr.
const Unit *FindUnit(const std::string &token);
/// A unit whose token or alias equals `token` ignoring case (for a syntax-only fix), or nullptr.
const Unit *FindUnitIgnoringCase(const std::string &token);
/// The units of a dimension, in table order.
std::vector<std::string> UnitsOfDimension(const std::string &dimension);
const Dimension *FindDimension(const std::string &name);

/// Beam particles per unit of charge of a per-charge unit (1/uC: 1e-6 C / (q e)); q >= 1 is the charge state.
double ParticlesPerUnitCharge(const Unit &unit, int chargeState);

}  // namespace experiment
}  // namespace azure2

#endif
