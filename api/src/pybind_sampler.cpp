// The engine's ensemble sampler from Python (docs/upgrade/pyazr3.md section 12; pyazr/session_sampler.py):
// Session.sample_with_engine.  The run is without the GIL; the progress callable and the stop request are looked at
// after every step, with the GIL taken.
#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "AZUREAPI.h"
#ifdef USE_MCMC
#include "AZUREAPISampler.h"
#endif

namespace py = pybind11;

#ifndef USE_MCMC
py::dict SampleWithEngineDict(AZUREAPI &, const py::kwargs &) {
  throw std::invalid_argument("this build of the engine has no ensemble sampler (cmake -DUSE_MCMC=ON); "
                              "pyazr.describe()[\"engine\"][\"mcmc\"] says whether a build has it.");
}
#else
namespace {

py::array_t<double> ArrayOf(const std::vector<double> &v, std::vector<py::ssize_t> shape) {
  py::array_t<double> a(std::move(shape));
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

}  // namespace

/// The session's sampler run.  `options`, each required (pyazr/session_sampler.py gives the defaults): walkers,
/// steps, sampler_threads (whole numbers), basis ("rwa" or "physical"), parameter_spread (percent), energy_spread
/// (keV), overwrite_samples (bool), x0 (the free start vector; empty: the session's values), progress (a callable or
/// None), stop (an object with is_set() or None).
py::dict SampleWithEngineDict(AZUREAPI &api, const py::kwargs &options) {
  const auto given = [&options](const char *name) {
    if (!options.contains(name)) throw py::type_error(std::string("sample_with_engine: missing keyword ") + name);
    return py::reinterpret_borrow<py::object>(options[name]);
  };
  const py::object progress = given("progress"), stop = given("stop");
  EngineSampleRequest request;
  // std::invalid_argument: the command line's refusal of an answer
  request.options = SamplerRequest(given("walkers").cast<int>(), given("steps").cast<int>(),
                                   given("sampler_threads").cast<int>(), given("basis").cast<std::string>(),
                                   given("parameter_spread").cast<double>(), given("energy_spread").cast<double>(),
                                   given("overwrite_samples").cast<bool>());
  request.x0 = given("x0").cast<std::vector<double>>();
  // A Python error in the progress callable (or in stop.is_set) ends the run and is raised when it has ended.
  std::unique_ptr<py::error_already_set> error;
  bool stopRequested = false;
  int lastStep = -1;
  if (!progress.is_none() || !stop.is_none()) {
    request.hooks.afterStep = [&](int step, int total, double logProbability, double logLikelihood, double logPrior) {
      if (error || stopRequested) return true;
      if (step == lastStep) return false;  // one record per step
      lastStep = step;
      py::gil_scoped_acquire gil;
      try {
        if (!progress.is_none()) {
          py::dict d;
          d["step"] = step;
          d["steps"] = total;
          d["log_prob"] = logProbability;
          d["log_likelihood"] = logLikelihood;
          d["log_prior"] = logPrior;
          progress(d);
        }
        if (!stop.is_none() && stop.attr("is_set")().cast<bool>()) stopRequested = true;
      } catch (py::error_already_set &e) {
        error.reset(new py::error_already_set(std::move(e)));
      }
      return error || stopRequested;
    };
  }
  EngineSample sample;
  {
    py::gil_scoped_release release;
    sample = SampleWithEngine(api, request);
  }
  if (error) throw *error;
  py::dict d;
  d["outcome"] = sample.outcome == SamplerOutcome::kSampled   ? "sampled"
                 : sample.outcome == SamplerOutcome::kRefused ? "refused"
                                                              : "not_sampled";
  d["stop_requested"] = stopRequested;
  d["steps"] = sample.steps;
  d["walkers"] = sample.walkers;
  d["parameters"] = sample.parameters;
  d["samples_file"] = sample.samplesFile;
  d["walkers_file"] = sample.walkersFile;
  d["output_dir"] = sample.outputDirectory;
  d["chain"] = ArrayOf(sample.chain, {sample.steps, sample.walkers, sample.parameters});
  d["log_prob"] = ArrayOf(sample.logProbability, {sample.steps, sample.walkers});
  d["log_likelihood"] = ArrayOf(sample.logLikelihood, {sample.steps, sample.walkers});
  d["log_prior"] = ArrayOf(sample.logPrior, {sample.steps, sample.walkers});
  return d;
}
#endif
