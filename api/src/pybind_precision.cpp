// The x87 control words of the threads that run the engine: a diagnostic and a test hook
// (azure2/support/FloatingPointPrecision.h; tests/pyazr/x87_precision_test.py).
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "azure2/support/FloatingPointPrecision.h"
#include "pybind_engine_call.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

namespace {

// What the threads of a call hold: the calling thread's word as the call found it and inside the engine's scope, and
// the word each thread of one OpenMP region found at the region's start (as the runtime made or left it) and holds
// inside the region's scope.  All 0 where there is no x87 (arm64).
py::dict ControlWords() {
  py::dict d;
  d["calling"] = azure2::support::X87ControlWord();
  std::vector<unsigned> inCall, found, inRegion;
  {
    EngineCallScope engine;
    inCall.push_back(azure2::support::X87ControlWord());
#ifdef _OPENMP
    const int n = omp_get_max_threads();
#else
    const int n = 1;
#endif
    found.assign(n, 0);
    inRegion.assign(n, 0);
#pragma omp parallel num_threads(n)
    {
      const azure2::support::ExtendedPrecisionScope fpPrecision;
#ifdef _OPENMP
      const int t = omp_get_thread_num();
#else
      const int t = 0;
#endif
      found[t] = fpPrecision.Found();
      inRegion[t] = azure2::support::X87ControlWord();
    }
  }
  d["in_call"] = inCall[0];
  d["omp_found"] = found;
  d["omp_in_region"] = inRegion;
  d["after_call"] = azure2::support::X87ControlWord();
  return d;
}

}  // namespace

void BindPrecision(py::module_ &m) {
  m.def("_x87_control_words", &ControlWords,
        "The x87 control words (0x37F: 64-bit precision; 0x27F: 53-bit) of the calling thread before, during and "
        "after a call into the engine, and of each thread of an OpenMP region at its start and inside the region's "
        "scope; all 0 without an x87.  A diagnostic: azure2/support/FloatingPointPrecision.h.");
  m.def(
      "_set_x87_control_word",
      [](unsigned word) {
        azure2::support::SetX87ControlWord(static_cast<unsigned short>(word));
        return azure2::support::X87ControlWord();
      },
      py::arg("word"),
      "A test hook: load `word` into the calling thread's x87 control word (it stays after the call) and return the "
      "word now in force; 0 and no effect without an x87.  tests/pyazr/x87_precision_test.py lowers a thread's "
      "precision to 53 bits with it, as a new thread has it on Windows.");
}
