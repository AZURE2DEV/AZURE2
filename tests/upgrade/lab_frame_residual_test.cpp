// The residual of a datum is formed in the frame the datum was given in (docs/upgrade/energy_grid.md section 4,
// include/ModelPass.h DataResidual): r = (f/J - n y)/(n dy) with f the c.m. model, J = dOmega_lab/dOmega_cm of the row
// at its present energy shift, and y, dy as read.  Model: tests/upgrade/T2 (lab-differential segments with and
// without a fitted shift, c.m. and angle-integrated data, SUM and RATIO components, a data covariance, THM), evaluated
// at its start parameters and with every energy shift moved by 3 keV.  For every row of every segment:
//   1. the pull equals the c.m. comparison (f - n J y)/(n J dy) to rounding: |difference| <= 8 eps (|f/J| + |n y|)/(n dy);
//   2. where J is 1 (c.m. data, angle-integrated data, ratios) the pull equals (f - y n)/(dy n) bit for bit;
//   3. the c.m. value and error shown by the writers are J y and J dy, bit for bit;
// and T2 has rows of both kinds, with J moving with the shift.  Shown to fail (2026-11-11) with DataResidual
// multiplying the model by J instead of dividing (check 1).
//
//   lab_frame_residual_test <T2.azr>   (run in a directory holding T2's data/, cov_diff.txt and prior_cov.txt)
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "SetUpForTests.h"
#include "AZUREParams.h"
#include "ProjectSetup.h"
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

bool SameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

struct Counts {
  int lab = 0, unit = 0;               // rows with J != 1 and J == 1
  int bad1 = 0, bad2 = 0, bad3 = 0;    // rows failing checks 1, 2, 3
  double worst = 0.0;                  // largest |difference| / (8 eps scale) of check 1
  std::string first;
};

// Every row of every segment of `data` at the model values of `results`.
Counts Rows(const WorkingCopy &data, const azure2::Results &results) {
  Counts c;
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *seg = data.Data().GetSegment(s);
    if (!seg) continue;
    const SegmentState &state = data.State(*seg);
    const double n = state.norm;
    for (int i = 0; i < seg->NumPoints(); i++) {
      const double f = seg->GetPoints(data)[i].GetModelValue(results);
      const double j = Segment::Jacobian(state, 0, i), y = seg->LabValue(i), dy = seg->LabError(i);
      if (dy == 0.0 || n == 0.0) continue;
      const double pull = azure2::PointResidual(data, *seg, i, results, n).Pull();
      const double yCM = y * j, dyCM = dy * j;
      const double cm = (f - yCM * n) / (dyCM * n);
      const double scale = 8.0 * DBL_EPSILON * (std::fabs(f / j) + std::fabs(n * y)) / std::fabs(n * dy);
      const double d = std::fabs(pull - cm);
      if (scale > 0.0) c.worst = std::max(c.worst, d / scale);
      if (!(d <= scale)) {
        if (!c.bad1++) {
          std::ostringstream o;
          o.precision(17);
          o << "segment " << s << " row " << i << ": pull " << pull << ", c.m. form " << cm << ", J " << j;
          c.first = o.str();
        }
      }
      if (j == 1.0) {
        c.unit++;
        if (!SameBits(pull, (f - y * n) / (dy * n))) c.bad2++;
      } else {
        c.lab++;
      }
      if (!SameBits(seg->CMValue(state, i), yCM) || !SameBits(seg->CMError(state, i), dyCM)) c.bad3++;
    }
  }
  return c;
}

void Report(const std::string &at, const Counts &c) {
  std::printf("%s: %d rows with J != 1, %d with J = 1; largest |pull - c.m. form| %.3g of the 8-eps bound\n",
              at.c_str(), c.lab, c.unit, c.worst);
  check(at + ": 1. pull = (f - n J y)/(n J dy) to 8 eps", c.bad1 == 0, std::to_string(c.bad1) + " rows; " + c.first);
  check(at + ": 2. J = 1: pull = (f - y n)/(dy n) bit for bit", c.bad2 == 0, std::to_string(c.bad2) + " rows");
  check(at + ": 3. c.m. value and error = J y, J dy bit for bit", c.bad3 == 0, std::to_string(c.bad3) + " rows");
  check(at + ": rows of both kinds", c.lab > 0 && c.unit > 0);
}

std::vector<double> Jacobians(const WorkingCopy &data) {
  std::vector<double> js;
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *seg = data.Data().GetSegment(s);
    if (!seg) continue;
    for (int i = 0; i < seg->NumPoints(); i++) js.push_back(Segment::Jacobian(data.State(*seg), 0, i));
  }
  return js;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: lab_frame_residual_test <T2.azr>\n");
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
  const vector_r p0 = mp.Params();
  vector_r p1 = p0;
  int shifts = 0;
  for (unsigned int i = 0; i < p0.size(); i++)
    if (params.Record(i).kind == ParamKind::EnergyShift) p1[i] = p0[i] + 0.003, shifts++;
  check("T2 has energy-shift parameters (" + std::to_string(shifts) + ")", shifts > 0);

  azure2::Results results;
  auto evaluate = [&](const vector_r &p) {
    compound.FillCompoundFromParams(p, state);
    data.FillNormsFromParams(p);
    data.FillEnergyShiftsFromParams(p, results, &compound, &state, &inputs, &cfg);
    if (cfg.paramMask & Config::USE_BRUNE_FORMALISM) compound.CalcShiftFunctions(cfg, state, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
    azure2::EvaluateModel(data, &compound, state, cfg, results, azure2::LeafValues::Calculate);
  };
  evaluate(p0);
  const std::vector<double> j0 = Jacobians(data);
  Report("start", Rows(data, results));
  evaluate(p1);
  Report("shifts + 3 keV", Rows(data, results));
  const std::vector<double> j1 = Jacobians(data);
  int moved = 0;
  for (size_t k = 0; k < j0.size() && k < j1.size(); k++) moved += !SameBits(j0[k], j1[k]);
  check("the shift moves the Jacobian of " + std::to_string(moved) + " rows", moved > 0);

  std::printf("lab_frame_residual: %s (%d failure%s)\n", failures ? "FAIL" : "all checks passed", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
