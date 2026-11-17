// The array conversions of the binding's session (pybind_api.cpp) and of its THM reports (pybind_thm_reports.cpp):
// a std::vector<double> to a numpy array and back, each a copy.  Functions only; no arithmetic.
#ifndef AZURE2_API_PYBIND_ARRAYS_H
#define AZURE2_API_PYBIND_ARRAYS_H

#include <cstring>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Constants.h"

// Convert a std::vector<double> to a new float64 numpy array (a copy, so the
// array owns its memory and outlives the C++ object).
inline pybind11::array_t<double> to_array(const vector_r &v) {
  return pybind11::array_t<double>(static_cast<pybind11::ssize_t>(v.size()), v.data());
}

// The mirror image: a flat copy of any float64-coercible array.  The engine
// takes vector_r&, so it needs its own buffer either way.  Bindings that
// release the GIL for the whole call (py::call_guard) take their array by
// const reference: a by-value py::array is reference-counted inside the
// GIL-free call (copy in, destroy out), and a temporary converted from a
// Python list then crashed (segfault on calculate_chi2_rwa([...]), October
// 2026).  The others convert their argument with the GIL held and release it
// only around the engine call (an inner block: EngineCallScope, so that
// the GIL is taken back before the result is converted); they may take the
// array by value.
inline vector_r to_vector(const pybind11::array_t<double, pybind11::array::forcecast> &a) {
  vector_r v(a.size());
  if (a.size() > 0) std::memcpy(v.data(), a.data(), v.size() * sizeof(double));
  return v;
}

#endif
