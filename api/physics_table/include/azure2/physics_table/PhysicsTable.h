#ifndef AZURE2_PHYSICS_TABLE_PHYSICS_TABLE_H
#define AZURE2_PHYSICS_TABLE_PHYSICS_TABLE_H

// The table of physics functions (docs/upgrade/pyazr3.md section 2.7, docs/upgrade/physics_tab.md section 2): one row
// per function of physics/, physics/transfer or experiment/ that a front end evaluates element-wise on arrays, with
// its name, the arguments (unit, domain, example, sweep), the returned fields (unit), the docstring and an adapter
// that calls the library function on one element.  Plain C++: no pybind11, no Qt, no engine header.  The binding
// (api/src/pybind_physics.cpp: pyazr.physics, pyazr.kinds()["physics"]) and the front ends' facade (api/front:
// azure2::front::PhysicsFunctions, EvaluatePhysics) read the same rows here.
//
// The row types (an argument, a returned field, the domain check of one element, the domain in words) are shared
// with the stopping table of the binding (api/src/pybind_stopping.cpp, through pybind_rows.h).
//
// A new function is one row of BuildTable() (src/PhysicsTable.cpp) and a test naming it (tests/consistency check 41).
// No physics in the table: an adapter converts arguments (doubles to integers, masses to a ChannelKinematics) and
// calls the library.

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace physics_table {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// Status of one element: 0 in the domain and computed; k > 0 argument k (1-based, in the order of the row's numeric
/// arguments) outside its domain; kJointly the arguments together outside the domain (each in its own); kFailed not
/// converged or failed.
constexpr int kOk = 0;
constexpr int kFailed = -1;
constexpr int kJointly = -2;

/// real, integer, flag (0 or 1), choice (the index of one of `choices`; Python passes the name): numbers, one per
/// element.  text: one string per call, not broadcast (a material's formula); its example is `example_text`.
enum class ArgKind { real, integer, flag, choice, text };

/// The spacing of a sweep: the points of a first plot over an argument.
enum class Spacing { linear, log };

/// One argument of a row: name, unit, meaning, kind and domain [lo, hi] (open ends where said); a default makes it
/// optional in Python.  Every numeric argument must be finite; an integer a whole number; a flag 0 or 1; a choice the
/// index of one of `choices`.  For the front ends' Physics pages: an example (a typical value, the form's start value
/// where there is no default) and a sweep (the range of a first plot over this argument; an argument with a sweep may
/// be a plot's axis).  Neither changes a computed number.
struct ArgSpec {
  std::string name, unit, doc;
  ArgKind kind = ArgKind::real;
  double lo = -kInf, hi = kInf;
  bool lo_open = false, hi_open = false;
  bool has_default = false;
  double default_value = 0.0;
  std::vector<std::string> choices;
  bool has_example = false;
  double example_value = 0.0;
  std::string example_text;
  bool has_sweep = false;
  double sweep_lo = 0.0, sweep_hi = 0.0;
  Spacing sweep_spacing = Spacing::linear;

  ArgSpec gt(double v) const { ArgSpec a = *this; a.lo = v; a.lo_open = true; return a; }
  ArgSpec ge(double v) const { ArgSpec a = *this; a.lo = v; a.lo_open = false; return a; }
  ArgSpec lt(double v) const { ArgSpec a = *this; a.hi = v; a.hi_open = true; return a; }
  ArgSpec le(double v) const { ArgSpec a = *this; a.hi = v; a.hi_open = false; return a; }
  ArgSpec otherwise(double v) const { ArgSpec a = *this; a.has_default = true; a.default_value = v; return a; }
  ArgSpec example(double v) const { ArgSpec a = *this; a.has_example = true; a.example_value = v; return a; }
  ArgSpec sweep(double lo, double hi, Spacing spacing) const {
    ArgSpec a = *this;
    a.has_sweep = true, a.sweep_lo = lo, a.sweep_hi = hi, a.sweep_spacing = spacing;
    return a;
  }
};

inline ArgSpec Real(const std::string &name, const std::string &unit, const std::string &doc) {
  ArgSpec a;
  a.name = name, a.unit = unit, a.doc = doc;
  return a;
}
inline ArgSpec Integer(const std::string &name, const std::string &doc) {
  ArgSpec a = Real(name, "", doc);
  a.kind = ArgKind::integer;
  return a;
}
inline ArgSpec Flag(const std::string &name, const std::string &doc, bool otherwise) {
  ArgSpec a = Real(name, "", doc).ge(0).le(1).otherwise(otherwise ? 1.0 : 0.0);
  a.kind = ArgKind::flag;
  return a;
}
inline ArgSpec Choice(const std::string &name, const std::string &doc, std::vector<std::string> choices,
                      int otherwise) {
  ArgSpec a = Real(name, "", doc).ge(0).le(static_cast<double>(choices.size()) - 1).otherwise(otherwise);
  a.kind = ArgKind::choice;
  a.choices = std::move(choices);
  return a;
}
/// A text argument with its example (one value per call).
inline ArgSpec Text(const std::string &name, const std::string &doc, const std::string &example) {
  ArgSpec a = Real(name, "", doc);
  a.kind = ArgKind::text;
  a.example_text = example;
  return a;
}

/// One returned field: name, unit, meaning.
struct OutSpec {
  std::string name, unit, doc;
};

/// 0 if every numeric argument a[k] is in the domain of args[k], else k + 1 of the first that is not.
inline int CheckDomain(const std::vector<ArgSpec> &args, const double *a) {
  for (std::size_t k = 0; k < args.size(); k++) {
    const ArgSpec &s = args[k];
    const double x = a[k];
    bool ok = std::isfinite(x);
    if (ok && s.kind != ArgKind::real) ok = x == std::floor(x);
    if (ok) ok = s.lo_open ? x > s.lo : x >= s.lo;
    if (ok) ok = s.hi_open ? x < s.hi : x <= s.hi;
    if (!ok) return static_cast<int>(k) + 1;
  }
  return kOk;
}

inline const char *KindName(ArgKind k) {
  switch (k) {
    case ArgKind::integer: return "integer";
    case ArgKind::flag: return "flag";
    case ArgKind::choice: return "choice";
    case ArgKind::text: return "text";
    default: return "real";
  }
}

inline std::string Number(double x) {
  if (x == std::floor(x) && std::fabs(x) < 1e15) return std::to_string(static_cast<long long>(x));
  std::string s = std::to_string(x);
  while (!s.empty() && s.back() == '0') s.pop_back();
  return s;
}

/// The domain of an argument in words: "finite > 0", "integer in [0, 10000]", "0 or 1", "one of ...", "text".
inline std::string DomainText(const ArgSpec &s) {
  if (s.kind == ArgKind::flag) return "0 or 1 (False or True)";
  if (s.kind == ArgKind::text) return "text";
  if (s.kind == ArgKind::choice) {
    std::string out = "one of";
    for (std::size_t i = 0; i < s.choices.size(); i++) out += (i ? ", " : " ") + s.choices[i];
    return out;
  }
  std::string out = s.kind == ArgKind::integer ? "integer" : "finite";
  const bool lo = std::isfinite(s.lo), hi = std::isfinite(s.hi);
  if (lo && hi)
    out += std::string(" in ") + (s.lo_open ? "(" : "[") + Number(s.lo) + ", " + Number(s.hi) +
           (s.hi_open ? ")" : "]");
  else if (lo)
    out += (s.lo_open ? " > " : " >= ") + Number(s.lo);
  else if (hi)
    out += (s.hi_open ? " < " : " <= ") + Number(s.hi);
  return out;
}

/// The adapter of a row: reads the row's arguments a[0 .. n_args), writes its fields o[0 .. n_outs), returns a status.
using Adapter = int (*)(const double *a, double *o);

/// One function of the table.
struct PhysicsRow {
  std::string name;
  std::string summary;   ///< one line: what it computes, with units
  std::string doc;       ///< the formula, conventions, validity
  std::string library;   ///< the library function called (header: function)
  std::string axis;      ///< the argument a plot puts on its horizontal axis
  std::vector<ArgSpec> args;
  std::vector<OutSpec> outs;
  Adapter call;
  std::string group = "";  ///< the row's section of the table (set by BuildTable from its section lines)
};

/// The rows, in the table's order (built once; immutable).
const std::vector<PhysicsRow> &PhysicsTable();

/// The row named `name`, or nullptr.
const PhysicsRow *FindPhysicsRow(const std::string &name);

/// The row at one element: the domain check (CheckDomain), then the adapter.  `a` holds the row's args.size()
/// arguments (a choice as its index, a flag as 0 or 1), `o` receives outs.size() fields, NaN where the status is not 0.
/// Returns the status (kOk, k > 0, kJointly, kFailed, or a row's own positive code as its doc says).
inline int EvaluateElement(const PhysicsRow &row, const double *a, double *o) {
  int s = CheckDomain(row.args, a);
  if (s == kOk) s = row.call(a, o);
  if (s != kOk)
    for (std::size_t j = 0; j < row.outs.size(); j++) o[j] = kNaN;
  return s;
}

}  // namespace physics_table
}  // namespace azure2

#endif
