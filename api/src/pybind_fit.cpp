// The engine's fitters from Python (docs/upgrade/fit.md; pyazr/session_fit.py): the fitter table (fitters()), a stop
// request another thread can set (StopRequest), and the session's fit (FitWithEngineDict, Session.fit_with_engine).
// The fit runs without the GIL; the progress callable is called with the GIL taken, from the engine's threads.
#include <algorithm>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <atomic>
#include <cmath>
#include <memory>
#include <mutex>

#include <Minuit2/FunctionMinimum.h>

#include "AZUREAPI.h"
#include "Fitters.h"
#include "pybind_engine_call.h"

namespace py = pybind11;

namespace {

// The stop request a caller owns: set from any thread, read by the engine (FitControl::stop).
struct StopRequest {
  std::atomic<bool> flag{false};
};

py::array_t<double> ArrayOf(const std::vector<double> &v) {
  py::array_t<double> a(v.size());
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

}  // namespace

py::dict FitWithEngineDict(AZUREAPI &api, const std::string &method, bool minos, const std::vector<double> &x0,
                           unsigned int maxCalls, const py::object &progress, const py::object &stopObject) {
  StopRequest *stop = stopObject.is_none() ? nullptr : stopObject.cast<StopRequest *>();
  AZUREAPI::EngineFitRequest request;
  request.method = method;
  request.minos = minos;
  request.x0 = x0;
  request.maxCalls = maxCalls;
  StopRequest own;
  StopRequest *request_stop = stop ? stop : &own;
  request.stop = &request_stop->flag;
  const std::vector<bool> fixed = api.params_fixed();
  // A Python error in the progress callable stops the fit and is raised when the fit has ended.
  std::mutex errorMutex;
  std::atomic<bool> failed{false};
  std::unique_ptr<py::error_already_set> error;
  if (!progress.is_none()) {
    request.progress = [&](const FitProgress &record) {
      if (failed.load()) return;
      py::gil_scoped_acquire gil;
      try {
        py::dict d;
        d["fitter"] = record.fitter;
        d["stage"] = record.stage;
        d["kind"] = record.kind == FitProgress::Kind::kEvaluation ? "evaluation" : "iteration";
        d["evaluation"] = record.evaluation;
        d["iteration"] = record.iteration;
        d["calls"] = record.calls;
        d["objective"] = record.objective;
        d["best"] = record.best;
        d["measure"] = record.measure;
        d["error_analysis"] = record.errorAnalysis;
        if (record.x && record.x->size() == fixed.size()) {
          std::vector<double> free;
          for (size_t i = 0; i < fixed.size(); i++)
            if (!fixed[i]) free.push_back((*record.x)[i]);
          d["x"] = ArrayOf(free);
        } else {
          d["x"] = py::none();
        }
        progress(d);
      } catch (py::error_already_set &e) {
        std::lock_guard<std::mutex> lock(errorMutex);
        if (!error) error.reset(new py::error_already_set(std::move(e)));
        failed.store(true);
        request_stop->flag.store(true);
      }
    };
  }
  AZUREAPI::EngineFit fit;
  {
    EngineCallScope engine;
    fit = api.FitWithEngine(request);
  }
  if (error) throw *error;
  const FitResult &r = fit.result;
  py::dict d;
  d["method"] = r.fitter;
  d["x_full"] = ArrayOf(r.x);
  d["objective"] = fit.objective;
  d["fitter_objective"] = r.objective;
  d["status"] = r.status;
  d["converged"] = r.converged;
  d["edm"] = r.edm;
  d["edm_goal"] = r.edmGoal;
  d["calls"] = r.calls;
  d["free_indices"] = r.freeIndices;
  d["covariance"] = ArrayOf(r.covariance);
  d["errors_full"] = ArrayOf(r.errors);
  d["covariance_from_hesse"] = r.covarianceFromHesse;
  d["fell_back_to_migrad"] = r.fellBackToMigrad;
  py::list m;
  for (const auto &interval : r.minos) m.append(py::make_tuple(interval.first, interval.second));
  d["minos"] = m;
  const auto side = [](const MinosSide &s) {
    py::dict v;
    v["valid"] = s.valid;
    v["at_limit"] = s.atLimit;
    v["call_limit"] = s.callLimit;
    v["new_minimum"] = s.newMinimum;
    return v;
  };
  py::list mv;
  for (const MinosValidity &v : r.minosValidity) {
    py::dict sides;
    sides["lower"] = side(v.lower);
    sides["upper"] = side(v.upper);
    mv.append(sides);
  }
  d["minos_validity"] = mv;
  d["outcome"] = FitOutcome(r);
  d["notices"] = r.notices;
  d["log"] = fit.log;
  return d;
}

void BindFit(py::module_ &m) {
  py::class_<StopRequest>(m, "StopRequest",
                          "A request to stop an engine fit (Session.fit_with_engine), set from any thread.")
      .def(py::init<>())
      .def("set", [](StopRequest &s) { s.flag.store(true); }, "Ask the fit to stop.")
      .def("clear", [](StopRequest &s) { s.flag.store(false); }, "Withdraw the request (before the next fit).")
      .def("is_set", [](const StopRequest &s) { return s.flag.load(); }, "Whether the request is set.");
  m.def(
      "fitters",
      []() {
        py::list rows;
        for (const FitterRow &row : Fitters()) {
          py::dict d;
          d["name"] = row.name;
          d["option"] = row.option;
          d["gradient"] = row.gradient;
          d["residual_rows"] = row.residualRows;
          d["bounds"] = row.bounds;
          d["covariance"] = row.covariance;
          d["minos"] = row.minos;
          d["call_limit"] = row.callLimit;
          d["fallback"] = row.fallback;
          d["progress_label"] = row.progressLabel;
          d["doc"] = row.doc;
          rows.append(d);
        }
        return rows;
      },
      "The engine's fitters (the rows of Fitters()), in the program's precedence.");
}
