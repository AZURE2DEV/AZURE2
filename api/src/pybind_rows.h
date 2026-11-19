// The binding's side of the function tables: the physics table (rows in api/physics_table, read by
// pybind_physics.cpp) and the stopping table (pybind_stopping.cpp).  The row types -- an argument (name, unit, meaning,
// kind, domain, default, example, sweep), a returned field (name, unit, meaning), the domain check of one element and
// the domain in words -- are plain C++ (azure2/physics_table/PhysicsTable.h, shared with the front ends' facade); here
// is their description as plain JSON data (the rows of pyazr.kinds()["physics"] and ["stopping"], the web UI's Physics
// tab: docs/upgrade/web_ui.md section 6.6).  No physics here.
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "azure2/physics_table/PhysicsTable.h"

namespace azure2_rows {

namespace py = pybind11;

using azure2::physics_table::ArgKind;
using azure2::physics_table::ArgSpec;
using azure2::physics_table::CheckDomain;
using azure2::physics_table::Choice;
using azure2::physics_table::DomainText;
using azure2::physics_table::Flag;
using azure2::physics_table::Integer;
using azure2::physics_table::KindName;
using azure2::physics_table::kFailed;
using azure2::physics_table::kInf;
using azure2::physics_table::kJointly;
using azure2::physics_table::kNaN;
using azure2::physics_table::kOk;
using azure2::physics_table::Number;
using azure2::physics_table::OutSpec;
using azure2::physics_table::Real;
using azure2::physics_table::Spacing;
using azure2::physics_table::Text;

/// An argument as plain data: name, unit, doc, kind, domain in words, lower and upper (None: unbounded), lower_open,
/// upper_open, choices, default or None, example (a number, a text for a text argument) or None, sweep {lower,
/// upper, spacing} or None.
inline py::dict DescribeArgument(const ArgSpec &s) {
  py::dict d;
  d["name"] = s.name;
  d["unit"] = s.unit;
  d["doc"] = s.doc;
  d["kind"] = KindName(s.kind);
  d["domain"] = DomainText(s);
  d["lower"] = std::isfinite(s.lo) ? py::object(py::float_(s.lo)) : py::object(py::none());  // None: unbounded
  d["upper"] = std::isfinite(s.hi) ? py::object(py::float_(s.hi)) : py::object(py::none());
  d["lower_open"] = s.lo_open;
  d["upper_open"] = s.hi_open;
  d["choices"] = s.choices;
  if (!s.has_default)
    d["default"] = py::none();
  else if (s.kind == ArgKind::choice)
    d["default"] = s.choices.at(static_cast<std::size_t>(s.default_value));
  else if (s.kind == ArgKind::flag)
    d["default"] = s.default_value != 0.0;
  else
    d["default"] = s.default_value;
  if (s.kind == ArgKind::text)
    d["example"] = s.example_text;
  else
    d["example"] = s.has_example ? py::object(py::float_(s.example_value)) : py::object(py::none());
  if (s.has_sweep) {
    py::dict w;
    w["lower"] = s.sweep_lo;
    w["upper"] = s.sweep_hi;
    w["spacing"] = s.sweep_spacing == Spacing::log ? "log" : "linear";
    d["sweep"] = w;
  } else {
    d["sweep"] = py::none();
  }
  return d;
}

/// A returned field as plain data: name, unit, doc.
inline py::dict DescribeReturn(const OutSpec &s) {
  py::dict d;
  d["name"] = s.name;
  d["unit"] = s.unit;
  d["doc"] = s.doc;
  return d;
}

}  // namespace azure2_rows
