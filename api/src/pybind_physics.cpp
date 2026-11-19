// The table of physics functions of the _azure2 module (docs/upgrade/pyazr3.md section 2.7): the rows are plain C++,
// api/physics_table (azure2/physics_table/PhysicsTable.h, BuildTable() in its src/PhysicsTable.cpp), which the front
// ends' facade (api/front) reads too.  Here: one template (Evaluate) applies a row element-wise to equal-length arrays
// and returns the values and a status per element, and the rows' description as Python data; pyazr/physics.py builds
// one Python function per row from physics_functions() (signature, docstring, broadcasting, DomainError).
//
// The rows' description as Python data (arguments, returned fields) is pybind_rows.h, shared with the stopping table
// of pybind_stopping.cpp.
//
// A new function is one row of the table and a test naming it (tests/consistency check 41).
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Elements.h"
#include "azure2/physics_table/PhysicsTable.h"
#include "pybind_rows.h"
#include "pybind_engine_call.h"

namespace py = pybind11;
namespace phys = azure2::physics;
using namespace azure2_rows;
using azure2::physics_table::PhysicsRow;

namespace {

const PhysicsRow &Row(const std::string &name) {
  const PhysicsRow *row = azure2::physics_table::FindPhysicsRow(name);
  if (!row) throw std::invalid_argument("no physics function named '" + name + "' (pyazr.kinds()['physics'] lists them)");
  return *row;
}


// ---------------------------------------------------------------- the template

using Column = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// The row `name` on equal-length 1-d arrays, one per argument in the row's order: (values of shape (n_outs, n),
/// status of shape (n,)), NaN where the status is not 0.
py::tuple Evaluate(const std::string &name, const std::vector<Column> &columns) {
  const PhysicsRow &row = Row(name);
  const std::size_t na = row.args.size(), no = row.outs.size();
  if (columns.size() != na)
    throw std::invalid_argument(name + ": " + std::to_string(na) + " argument arrays expected, " +
                                std::to_string(columns.size()) + " given");
  const py::ssize_t n = na ? columns[0].size() : 0;
  std::vector<const double *> in(na);
  for (std::size_t k = 0; k < na; k++) {
    if (columns[k].ndim() != 1 || columns[k].size() != n)
      throw std::invalid_argument(name + ": the argument arrays must be 1-d of one length (broadcast them first)");
    in[k] = columns[k].data();
  }
  py::array_t<double> values({static_cast<py::ssize_t>(no), n});
  py::array_t<int> status(n);
  double *v = values.mutable_data();
  int *st = status.mutable_data();
  {
    EngineCallScope engine;
    std::vector<double> a(na), o(no);
    for (py::ssize_t i = 0; i < n; i++) {
      for (std::size_t k = 0; k < na; k++) a[k] = in[k][i];
      const int s = azure2::physics_table::EvaluateElement(row, a.data(), o.data());
      for (std::size_t j = 0; j < no; j++) v[j * n + i] = s == kOk ? o[j] : kNaN;
      st[i] = s;
    }
  }
  return py::make_tuple(values, status);
}

py::list Describe() {
  py::list rows;
  for (const PhysicsRow &r : azure2::physics_table::PhysicsTable()) {
    py::list args, outs;
    for (const ArgSpec &s : r.args) args.append(DescribeArgument(s));
    for (const OutSpec &s : r.outs) outs.append(DescribeReturn(s));
    py::dict row;
    row["name"] = r.name;
    row["summary"] = r.summary;
    row["doc"] = r.doc;
    row["library"] = r.library;
    row["axis"] = r.axis;
    row["group"] = r.group;
    row["arguments"] = args;
    row["returns"] = outs;
    rows.append(row);
  }
  return rows;
}

}  // namespace

void BindStopping(py::module_ &parent);  // pybind_stopping.cpp: stopping/, the submodule _azure2.stopping
void BindLevelDoors(py::module_ &m);     // pybind_level_doors.cpp: one block's formal <-> observed levels (pyazr.endf)
void BindCollisionMatrix(py::module_ &m);  // pybind_collision_matrix.cpp: one block's A, U and T at one energy
void BindOptical(py::module_ &m);          // pybind_optical.cpp: elastic scattering on an optical potential

void BindPhysics(py::module_ &m) {
  BindStopping(m);
  BindLevelDoors(m);
  BindCollisionMatrix(m);
  BindOptical(m);
  m.def("physics_functions", &Describe,
        "The table of physics functions, one dict per row: name, summary, doc, library (the function called), axis "
        "(the argument a plot puts on its horizontal axis), group (its section of the table), arguments (name, unit, "
        "doc, kind real / integer / flag / choice, domain in words, lower and upper (None: unbounded), lower_open, "
        "upper_open, choices, default or None, example (a typical value) or None, sweep {lower, upper, spacing "
        "linear / log} (the range of a first plot over it) or None) and returns (name, unit, doc); plain JSON data.  pyazr.physics builds one Python function per row from it.");
  m.def("physics_evaluate", &Evaluate, PrecisionGuard(), py::arg("name"), py::arg("arguments"),
        "Row `name` element-wise on equal-length 1-d float arrays, one per argument in the row's order (a choice as "
        "its index, a flag as 0 or 1): (values, status), values of shape (n_returns, n), NaN where status != 0; "
        "status 0 computed, k > 0 argument k (1-based) outside its domain, -2 the arguments together outside the "
        "domain, -1 not converged or failed.");
  m.def(
      "element_symbols",
      []() {
        std::vector<std::string> out;
        for (int z = 0; z <= phys::kMaxElementZ; z++) out.push_back(phys::ElementSymbol(z));
        return out;
      },
      "Element symbols by charge number, index 0 (the neutron, \"n\") to 118 (physics/Elements.h, the one table).");
}
