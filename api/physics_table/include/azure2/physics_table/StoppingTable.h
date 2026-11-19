#ifndef AZURE2_PHYSICS_TABLE_STOPPING_TABLE_H
#define AZURE2_PHYSICS_TABLE_STOPPING_TABLE_H

// The table of stopping functions (docs/upgrade/web_ui.md section 11.5, docs/upgrade/physics_tab.md section 6.2): one
// row per function of stopping/ (SRIM-2013's stopping, energy loss, range and straggling), of the physics table's
// row types (PhysicsTable.h), each with SRIM's citation.  Plain C++: no pybind11, no Qt, no engine header.  The binding
// (api/src/pybind_stopping.cpp: pyazr.stopping.table(), pyazr.kinds()["stopping"]) and the front ends' facade
// (api/front: azure2::front::StoppingFunctions, EvaluateStopping) read the same rows here.
//
// Every row takes the material (its formula, one text per call, then the numeric phase and compound correction) and
// the projectile (z1, mass_u) first (Common()), then its own arguments.  Units as stopping/'s Material.h: energies
// MeV (lab, of the projectile), stopping cross sections eV/(1e15 atoms/cm^2) per atom of the material, areal densities
// atoms/cm^2 of the material, straggling widths MeV.

#include <cstddef>
#include <string>
#include <vector>

#include "azure2/physics_table/PhysicsTable.h"
#include "azure2/stopping/Material.h"

namespace azure2 {
namespace physics_table {

/// What a row's adapter is given besides the element's numbers: the material and the projectile of the element.
struct StoppingCall {
  const azure2::stopping::Material &material;
  const azure2::stopping::Projectile &projectile;
};

/// The adapter of a stopping row: the element's numeric arguments after the common ones (a[0 .. n)), its fields
/// o[0 .. n_outs); returns a status (PhysicsTable.h).
using StoppingAdapter = int (*)(const StoppingCall &c, const double *a, double *o);

/// A row: the common arguments (the formula, then Common()) come first in every row, then `args`; returns; adapter.
struct StoppingRow {
  std::string name, summary, doc, library, axis;
  std::vector<ArgSpec> args;  ///< the row's own numeric arguments, after the common ones
  std::vector<OutSpec> outs;
  StoppingAdapter call;
};

/// The phase choices: "default" is SRIM's phase of the material.
const std::vector<std::string> &StoppingPhases();
/// The numeric arguments every row takes first: phase, compound_correction, z1, mass_u.
constexpr std::size_t kStoppingCommon = 4;
/// The text argument every row takes before them: the material's formula.
ArgSpec StoppingFormula();
std::vector<ArgSpec> StoppingCommon();

/// The rows, in the table's order (built once; immutable).
const std::vector<StoppingRow> &StoppingTable();
/// The row named `name`, or nullptr.
const StoppingRow *FindStoppingRow(const std::string &name);
/// Every numeric argument of a row in order: the common ones, then the row's own.
std::vector<ArgSpec> StoppingNumericArgs(const StoppingRow &row);

/// The material of a formula, its phase ("" = SRIM's default, "solid", "gas") and SRIM's compound correction; throws
/// what stopping/ throws for a formula it refuses.
azure2::stopping::Material StoppingMaterial(const std::string &formula, const std::string &phase,
                                            double compoundCorrection);
/// The projectile of charge z1 and atomic mass mass_u (u); std::invalid_argument for a mass <= 0.
azure2::stopping::Projectile StoppingProjectile(int z1, double massU);

/// Row `row` for the material `formula` on n elements: columns[k] points at the n values of numeric argument k
/// (StoppingNumericArgs order; a choice as its index); values[j * n + i] receives return j at element i (NaN where
/// status[i] != 0).  The material and projectile are made once per distinct (phase, compound correction, z1, mass) in
/// a run of elements.  Status 0, k > 0 argument k outside its domain, kJointly the arguments together outside it (a
/// layer at least the range), kFailed failed.  The formula must have been accepted (StoppingMaterial) by the caller.
void EvaluateStoppingRow(const StoppingRow &row, const std::string &formula, const std::vector<const double *> &columns,
                         std::size_t n, double *values, int *status);

}  // namespace physics_table
}  // namespace azure2

#endif
