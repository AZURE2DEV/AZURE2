// The analytic gradient over the energy tasks (src/GradEnergyTasks.cpp) on tests/upgrade/T2 (folded and unfolded
// points, SUM and RATIO components, natural targets, shifts, a data covariance), at its start parameters:
//   1. the scalar gradient (AccumulateEGammaGradient), the residual Jacobian (ComputeResidualJacobian) and the model
//      gradients (ComputeModelGradients) are the same, bit for bit, at 1 and at 3 threads: the sums are formed in an
//      order fixed by the tasks (64 interleaved chunks; rows task by task), not by the threads;
//   2. with every parameter fixed (no packed column) the residual Jacobian still gives the residuals, equal to those
//      of a plain evaluation, and an empty Jacobian.  Shown to fail (2026-11-12): before the fix the row pass divided
//      by the column count and wrote past its buffer (a least-squares fit of thick_target_resonance, which has no
//      free parameter, crashed).
//
//   energy_task_gradient_test <T2.azr>   (run in a directory holding T2's data/, cov_diff.txt and prior_cov.txt)
#include <cstdio>
#include <cstring>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SetUpForTests.h"
#include "AZUREGrad.h"
#include "ProjectSetup.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ModelPass.h"
#include "Results.h"
#include "Segment.h"
#include "WorkingCopy.h"

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

bool SameBits(const std::vector<double> &a, const std::vector<double> &b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

void SetThreads(int n) {
#ifdef _OPENMP
  omp_set_num_threads(n);
#else
  (void)n;
#endif
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: energy_task_gradient_test <T2.azr>\n");
    return 2;
  }
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = argv[1];
  WorkingCopy data;
  AZUREParams params;
  CompoundSetUp setUp;
  if (cfg.ReadConfigFile() < 0 || (setUp = test_setup::WithData(cfg, data)).refusal != CompoundSetUp::Refusal::kNone) {
    check("T2 loads", false, sink.str().substr(0, 400));
    return 1;
  }
  const CompoundSystem &compound = *setUp.compound;
  LevelState &state = setUp.state;  // the levels, as the run holds them
  LevelInputs &inputs = setUp.inputs;
  compound.FillMnParams(params, state, inputs, &cfg);
  data.FillMnParams(params);
  if (data.Initialize(&compound, cfg) == -1) {
    check("T2 initializes", false, sink.str().substr(0, 400));
    return 1;
  }
  const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const vector_r p = mp.Params();
  std::vector<bool> fixed(p.size());
  for (unsigned int i = 0; i < p.size(); i++) fixed[i] = mp.Parameter(i).IsFixed();

  azure2::Results results;
  compound.FillCompoundFromParams(p, state);
  data.FillNormsFromParams(p);
  data.FillEnergyShiftsFromParams(p, results, &compound, &state, &inputs, &cfg);
  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (cfg.paramMask & Config::USE_BRUNE_FORMALISM) {
    compound.CalcShiftFunctions(cfg, state, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
    shiftDeriv = BuildShiftDerivTable(&compound, state, cfg);
    sdp = &shiftDeriv;
  }
  const ParamIndexMap pmap = BuildParamIndexMap(&compound, &data, fixed);

  // 1. Thread independence.
  auto scalar = [&](int threads) {
    SetThreads(threads);
    GradAccum accum;
    accum.Init(&compound, state);
    FitBarFn fitBar = [](const Segment *, int, int, double model) { return model; };  // d(sum m^2 / 2)/dm
    const bool ok = AccumulateEGammaGradient(&compound, state, &data, results, cfg, pmap, sdp, fitBar, accum);
    vector_r g(pmap.NumFull(), 0.0);
    if (ok) accum.Scatter(pmap, g);
    return std::make_pair(ok, g);
  };
  auto jacobian = [&](int threads, const ParamIndexMap &map, vector_r &r) {
    SetThreads(threads);
    vector_r jac;
    int nCols = 0;
    const bool ok = ComputeResidualJacobian(&compound, state, &data, results, cfg, map, sdp, r, jac, nCols, true);
    return std::make_pair(ok, jac);
  };
  auto models = [&](int threads) {
    SetThreads(threads);
    std::map<Unit *, vector_r> rows;
    const bool ok = ComputeModelGradients(&compound, state, &data, results, cfg, pmap, sdp, rows, true);
    vector_r flat;
    for (int s = 1; s <= data.Data().NumSegments(); s++)
      for (int i = 1; i <= data.Data().GetSegment(s)->NumPoints(); i++) {
        auto it = rows.find(data.Data().GetSegment(s)->GetPoint(data, i));
        if (it != rows.end()) flat.insert(flat.end(), it->second.begin(), it->second.end());
      }
    return std::make_pair(ok, flat);
  };
  const auto s1 = scalar(1), s3 = scalar(3);
  check("scalar gradient: analytic at 1 and 3 threads", s1.first && s3.first);
  check("scalar gradient: the same bits at 1 and 3 threads (" + std::to_string(s1.second.size()) + " entries)",
        SameBits(s1.second, s3.second));
  vector_r r1, r3;
  const auto j1 = jacobian(1, pmap, r1), j3 = jacobian(3, pmap, r3);
  check("residual Jacobian: analytic at 1 and 3 threads", j1.first && j3.first);
  check("residual Jacobian: the same bits at 1 and 3 threads (" + std::to_string(j1.second.size()) + " entries)",
        SameBits(j1.second, j3.second) && SameBits(r1, r3));
  const auto m1 = models(1), m3 = models(3);
  check("model gradients: the same bits at 1 and 3 threads (" + std::to_string(m1.second.size()) + " entries)",
        m1.first && m3.first && SameBits(m1.second, m3.second));

  // 2. No free parameter.
  const ParamIndexMap none = BuildParamIndexMap(&compound, &data, std::vector<bool>(p.size(), true));
  vector_r r0;
  const auto j0 = jacobian(3, none, r0);
  check("no free parameter: the residual Jacobian runs, with no column", j0.first && j0.second.empty() && none.NumPacked() == 0);
  check("no free parameter: the residuals are those of the free-parameter pass (" + std::to_string(r0.size()) + " rows)",
        SameBits(r0, r1));

  std::printf("energy_task_gradient: %s (%d failure%s)\n", failures ? "FAIL" : "all checks passed", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
