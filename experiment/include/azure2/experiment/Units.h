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

/*!
 * v x factor, rounded as a decimal unit should be: where the factor is (to 1e-12) an integer N, v N; where its inverse
 * is, v / N; else v x factor.  Dividing by the exact integer 1000 gives the correctly rounded quotient, so a value
 * written as 12.5 mb becomes the same double as 0.0125 written in b (multiplying by the double nearest 1e-3 can be one
 * unit in the last place off).  factor 1 returns v itself.  The one place a unit factor is applied to a number.
 */
double Rescale(double value, double factor);
/// The factor that takes a value in `from` to the unit `to` of the same dimension (from.factor / to.factor).
double FactorBetween(const Unit &from, const Unit &to);
/// A value written in `unit` in the base unit of its dimension (the engine's): Rescale(value, unit->factor); the value
/// itself when there is no unit.
double ToBase(double value, const Unit *unit);
/// A value written in `from` in the unit token `to` (a row of the table, same dimension; throws otherwise).
double ToUnit(double value, const Unit *from, const std::string &to);

/// Beam particles per unit of charge of a per-charge unit (1/uC: 1e-6 C / (q e)); q >= 1 is the charge state.
double ParticlesPerUnitCharge(const Unit &unit, int chargeState);

}  // namespace experiment
}  // namespace azure2

#endif
