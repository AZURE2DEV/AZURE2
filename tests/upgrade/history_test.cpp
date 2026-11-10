// A working copy's state is a function of the parameters alone (docs/upgrade/data_model.md S1, section 1.3): a copy
// that evaluated other parameters before gives the bits of a fresh copy.  A fit evaluates pooled copies of the data
// set (AZURECalc::GetPooledData), each carrying what its earlier evaluations left -- applied energy shifts, sub-point
// grids, the mapping, the angle of a point a shift moved -- and which copy serves which evaluation depends on the
// scheduling.  Anything a copy keeps that is not a function of (theta, nu) makes the fit depend on the thread count.
//
// For each project, two parameter vectors: A, the project's start, and B, every level energy moved by 2 keV, every
// width by 1 % (fixed ones too: the copies are filled from the whole vector), every free norm by 2 % and every free
// energy-shift parameter by 1 keV (sqrt(E) coefficient and shift-function terms by 1e-3 of their unit).  Then
//   fresh A, fresh B           a copy of the master evaluated once;
//   pooled B -> A              a copy evaluated at B, then at A: its per-segment chi^2, total and model values equal
//                              fresh A's bit for bit;
//   pooled A -> B -> A -> B    the same copy back and forth: equal to fresh B's at the end.
// Also: a second load of the project in the same process gives its component segments the same keys, all distinct.
// An evaluation is AZURECalc::Chi2On's: the hand-over, the parallel point pass, the model value of every point, the
// objective (Objective.h) with the per-segment chi^2.
//
//   history_test <dir/project.azr> ...   (each project is read in its own directory: data paths are relative)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "Objective.h"
#include "ParameterKind.h"
#include "Results.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

// Bit for bit (a NaN equals a NaN with the same bits).
bool SameBits(const std::vector<double> &a, const std::vector<double> &b, std::string &why) {
  if (a.size() != b.size()) {
    why = std::to_string(a.size()) + " and " + std::to_string(b.size()) + " values";
    return false;
  }
  for (size_t i = 0; i < a.size(); i++)
    if (std::memcmp(&a[i], &b[i], sizeof(double)) != 0) {
      std::ostringstream s;
      s.precision(17);
      s << "first difference at " << i << ": " << a[i] << " and " << b[i];
      why = s.str();
      return false;
    }
  return true;
}

struct Evaluation {
  std::vector<double> segmentChi2, models;  // segmentChi2[0] = total
};

struct Copy {
  std::unique_ptr<CNuc> compound;
  std::unique_ptr<EData> data;
  azure2::Results results;
};

// One evaluation, as AZURECalc::Chi2On: hand-over, point pass, model values, objective.
Evaluation Evaluate(Copy &c, const vector_r &p, const Config &cfg) {
  CNuc *lc = c.compound.get();
  EData *ld = c.data.get();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, c.results, ld, lc, &cfg);
  if (cfg.paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(cfg, false);
  std::vector<ESegment *> segments;
  for (int i = 1; i <= ld->NumSegments(); i++) segments.push_back(ld->GetSegment(i));
  const bool pre = EData::CalculatePointsInParallel(segments, lc, cfg, c.results) >= 0;
  Evaluation e;
  for (ESegment *segment : segments)
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      const double m = segment->CalculateTheoreticalCrossSection(pid, lc, cfg, ld, c.results, pre);
      segment->GetPoint(pid + 1)->SetFitCrossSection(c.results, m);
      e.models.push_back(m);
    }
  ObjectiveContext context{lc, ld, c.results, p, cfg};
  context.segmentChi2 = &e.segmentChi2;
  const double total = Objective::Engine().Chi2(context);
  e.segmentChi2[0] = total;
  return e;
}

void Compare(const std::string &what, const Evaluation &got, const Evaluation &want) {
  std::string why;
  check(what + ": per-segment chi^2 and total", SameBits(got.segmentChi2, want.segmentChi2, why), why);
  check(what + ": model values", SameBits(got.models, want.models, why), why);
}

int RunProject(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  const std::string dir = slash == std::string::npos ? "." : path.substr(0, slash);
  const std::string file = slash == std::string::npos ? path : path.substr(slash + 1);
  if (chdir(dir.c_str()) != 0) {
    check(path + " directory", false);
    return 1;
  }
  std::printf("%s\n", path.c_str());
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file;
  g_config = &cfg;
  CNuc compound;
  EData data;
  AZUREParams params;
  if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1 || data.Fill(cfg, &compound) == -1) {
    check("loads", false, sink.str().substr(0, 400));
    return 1;
  }
  compound.Initialize(cfg);
  compound.FillMnParams(params, &cfg);
  data.FillMnParams(params);
  if (data.Initialize(&compound, cfg) == -1) {
    check("initializes", false, sink.str().substr(0, 400));
    return 1;
  }
  // The component segments' keys depend on the project alone: a second load in this process gives the same keys,
  // and no two components share one (H10: one counter per data set, not one per process and per kind).
  {
    std::ostringstream sink2;
    Config cfg2(sink2);
    cfg2.configfile = file;
    CNuc compound2;
    EData data2;
    AZUREParams params2;
    const bool loaded = cfg2.ReadConfigFile() >= 0 && compound2.Fill(cfg2) != -1 && data2.Fill(cfg2, &compound2) != -1;
    if (loaded) {
      compound2.Initialize(cfg2);
      compound2.FillMnParams(params2, &cfg2);
      data2.FillMnParams(params2);
    }
    check("a second load in the same process", loaded && data2.Initialize(&compound2, cfg2) != -1);
    g_config = &cfg;
    auto keys = [](EData &d) {
      std::vector<int> k;
      for (int i = 1; i <= d.NumSegments(); i++)
        for (const ESegment *c : d.GetSegment(i)->GetComponentSegments()) k.push_back(c->GetSegmentKey());
      return k;
    };
    const std::vector<int> k1 = keys(data), k2 = keys(data2);
    std::vector<int> sorted = k1;
    std::sort(sorted.begin(), sorted.end());
    check(std::to_string(k1.size()) + " component keys: the same in both loads, all distinct",
          k1 == k2 && std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end(),
          k1.empty() || k2.empty() ? "" : "first keys " + std::to_string(k1[0]) + " and " + std::to_string(k2[0]));
  }
  const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const vector_r pA = mp.Params();
  vector_r pB = pA;
  int moved = 0;
  for (unsigned int i = 0; i < pA.size(); i++) {
    const ParamKind kind = params.Kind((int)i);
    const bool free = !mp.Parameter(i).IsFixed();
    double v = pA[i];
    if (kind == ParamKind::LevelEnergy) v += 0.002;
    else if (kind == ParamKind::Gamma) v *= 1.01;
    else if (!free) continue;  // a fixed shift or norm keeps its value (see the report: a fixed shift is not undone)
    else if (kind == ParamKind::EnergyShift) v += 0.001;
    else if (kind == ParamKind::EnergyShiftSqrt || kind == ParamKind::EnergyShiftFn) v += 0.001;
    else if (kind == ParamKind::Norm) v *= 1.02;
    if (v != pA[i]) pB[i] = v, moved++;
  }
  check("B moves " + std::to_string(moved) + " parameters", moved > 0);

  auto copy = [&]() {
    Copy c;
    c.compound.reset(compound.Clone());
    c.data.reset(data.Clone());
    return c;
  };
  Copy freshA = copy(), freshB = copy();
  const Evaluation eA = Evaluate(freshA, pA, cfg), eB = Evaluate(freshB, pB, cfg);
  std::string why;
  std::printf("  %zu segments, %zu points; chi2 A %.17g, B %.17g\n", eA.segmentChi2.size() - 1, eA.models.size(),
              eA.segmentChi2[0], eB.segmentChi2[0]);
  check("A and B give different chi^2", eA.segmentChi2[0] != eB.segmentChi2[0]);

  Copy pooled = copy();
  Evaluate(pooled, pB, cfg);
  Compare("pooled B -> A = fresh A", Evaluate(pooled, pA, cfg), eA);
  Copy back = copy();
  Evaluate(back, pA, cfg);
  Evaluate(back, pB, cfg);
  Evaluate(back, pA, cfg);
  Compare("pooled A -> B -> A -> B = fresh B", Evaluate(back, pB, cfg), eB);
  return 0;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: history_test <dir/project.azr> ...\n");
    return 2;
  }
  char start[4096];
  if (!getcwd(start, sizeof(start))) return 2;
  for (int k = 1; k < argc; k++) {
    if (chdir(start) != 0) return 2;
    RunProject(argv[k]);
  }
  std::printf("history: %s (%d failure%s)\n", failures ? "FAIL" : "all checks passed", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
