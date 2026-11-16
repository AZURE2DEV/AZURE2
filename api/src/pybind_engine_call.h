#ifndef AZURE2_API_PYBIND_ENGINE_CALL_H
#define AZURE2_API_PYBIND_ENGINE_CALL_H

#include <pybind11/pybind11.h>

#include "azure2/support/FloatingPointPrecision.h"

/*!
 * The binding's entry into the engine: the GIL released for the length of the call, and the calling thread's x87
 * precision at 64 bits (azure2/support/FloatingPointPrecision.h), so that a call computes the same bits in a Python
 * thread as in the main thread on every system (on Windows a Python thread starts with 53 bits).  Every place of
 * api/src that releases the GIL takes this, as a local (`EngineCallScope engine;`) or as
 * `py::call_guard<EngineCallScope>()`, never pybind11::gil_scoped_release alone [repo_consistency check 60].
 * The precision is set before the GIL is released and put back after it is taken again.
 */
struct EngineCallScope {
  azure2::support::ExtendedPrecisionScope precision;
  pybind11::gil_scoped_release release;
};

#endif
