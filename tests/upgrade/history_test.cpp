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
//   pooled A -> B -> A -> B    the same copy back and forth: equal to fresh B's at the end;
//   pooled C -> A -> C         C = A with every fixed energy-shift parameter moved by 1 keV: equal to fresh A's, then
//                              to fresh C's (a shift that becomes inactive again must be undone).
//   copy of a copy             a copy made from a copy evaluated at B, evaluated at A; its source then evaluated at D
//                              (B with the free shifts 1 keV further): the copy's buffer does not move, and the copy
//                              at A again equals fresh A (no state shared with the source).
// Also: a second load of the project in the same process gives its leaves the same keys, all distinct; a copy shares
// its source's data set (DataSet.h: one object, never written after the set-up), and every evaluation and hand-over of
// the copies above leaves a fingerprint of that data set as it was after the set-up.
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
#include "DataSet.h"
#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "ModelPass.h"
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

// A fingerprint of what a data set holds that an evaluation could read: the settings, start values, priors and terms of
// every segment, its measurement's rows, the target effects' lines and priors, the parameter layout, the THM
// experiments and covariance groups (FNV-1a over the values' bits).
struct Fingerprint {
  std::uint64_t h = 1469598103934665603ULL;
  void Bytes(const void *p, size_t n) {
    const unsigned char *c = static_cast<const unsigned char *>(p);
    for (size_t i = 0; i < n; i++) h = (h ^ c[i]) * 1099511628211ULL;
  }
  void Add(double v) { Bytes(&v, sizeof v); }
  void Add(int v) { Bytes(&v, sizeof v); }
  void Add(bool v) { Add(v ? 1 : 0); }
  void Add(const std::string &v) { Bytes(v.data(), v.size()), Add((int)v.size()); }
  void Add(const std::vector<double> &v) {
    for (double x : v) Add(x);
    Add((int)v.size());
  }
};

std::uint64_t DataSetFingerprint(const DataSet &d) {
  Fingerprint f;
  f.Add(d.NumSegments());
  f.Add(d.GetNormParamOffset()), f.Add(d.GetEnergyShiftParamOffset()), f.Add(d.GetEnergyShiftSqrtParamOffset());
  f.Add(d.GetEnergyShiftFnParamOffset()), f.Add(d.ScreeningParamOffset()), f.Add(d.DensityParamOffset());
  f.Add(d.DopplerParamOffset()), f.Add(d.ThicknessParamOffset()), f.Add(d.ThmCoherentParamOffset());
  for (const DataSet::LeafRef &leaf : d.LeafRefs()) f.Add(leaf.first), f.Add(leaf.second);
  for (const ESegment &s : d.GetSegments()) {
    f.Add(s.GetSegmentKey()), f.Add(s.Pos()), f.Add(s.NumTerms()), f.Add(s.NumPoints());
    f.Add(s.StartNorm()), f.Add(s.GetNominalNorm()), f.Add(s.GetNormError()), f.Add(s.IsVaryNorm());
    f.Add(s.StartEnergyShift()), f.Add(s.GetNominalEnergyShift()), f.Add(s.GetEnergyShiftError());
    f.Add(s.StartEnergyShiftSqrt()), f.Add(s.GetEnergyShiftSqrtError()), f.Add(s.StartThickness());
    f.Add(s.GetNormLeaderPos()), f.Add(s.GetShiftLeaderPos()), f.Add(s.GetThicknessLeaderPos());
    f.Add(s.IsNormPriorReplaced()), f.Add(s.IsShiftPriorReplaced()), f.Add(s.IsThicknessPriorReplaced());
    f.Add(s.HasShiftFn()), f.Add(s.NumShiftFnParams()), f.Add(s.GetShiftFnParamStart());
    f.Add(s.GetNominalShiftFnValues()), f.Add(s.GetShiftFnErrors()), f.Add(s.GetCovarianceGroup());
    for (int t = 0; t < s.NumTerms(); t++) {
      f.Add(s.TermEntranceKey(t)), f.Add(s.TermExitKey(t)), f.Add(s.TermKey(t)), f.Add(s.TermTargetEffectNum(t));
      f.Add(s.LeafScaling(t));
      const TermSettings &ts = *s.TermSettingsOf(t);
      f.Add(ts.entranceKey), f.Add(ts.exitKey), f.Add(ts.differential), f.Add(ts.total), f.Add(ts.angularDist);
      f.Add(ts.maxAngDistOrder), f.Add(ts.crossSectionComponent);
    }
    if (s.HasMeasurement()) {
      const azure2::experiment::MeasurementRows &r = *s.GetMeasurement().rows;
      f.Add(r.energy), f.Add(r.angle), f.Add(r.value), f.Add(r.error);
    }
  }
  for (int n = 1; n <= d.NumTargetEffects(); n++) {
    const TargetEffect &te = *d.GetTargetEffect(n);
    f.Add(te.GetScreeningUe()), f.Add(te.GetDensity()), f.Add(te.GetDopplerTemperature());
    f.Add(te.GetDensityNominal()), f.Add(te.GetDensitySigma()), f.Add(te.GetScreeningSigma());
    f.Add(te.GetLayoutDensity()), f.Add(te.GetDopplerLayoutTemperature()), f.Add(te.IsScreeningPriorReplaced());
  }
  for (int k = 0; k < d.NumThmCoherentParams(); k++) f.Add(d.GetThmCoherentParam(k).value);
  for (int g = 0; g < d.NumThmGroups(); g++) f.Add(d.GetThmGroup(g).name), f.Add((int)d.GetThmGroup(g).segments.size());
  f.Add(d.NumCovarianceGroups());
  return f.h;
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
  ld->FillEnergyShiftsFromParams(p, c.results, lc, &cfg);
  if (cfg.paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(cfg, false);
  std::vector<const ESegment *> segments;
  for (int i = 1; i <= ld->Data().NumSegments(); i++) segments.push_back(ld->Data().GetSegment(i));
  const bool pre = ld->CalculatePointsInParallel(segments, lc, cfg, c.results) >= 0;
  Evaluation e;
  azure2::EvaluateModel(*ld, lc, cfg, c.results, pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
  e.models = c.results.model;
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
  // The leaves' keys depend on the project alone: a second load in this process gives the same keys,
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
      for (int i = 1; i <= d.Data().NumSegments(); i++)
        for (int t = 1; t < d.Data().GetSegment(i)->NumTerms(); t++) k.push_back(d.Data().GetSegment(i)->TermKey(t));
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
  // The data set is shared by the copies and only read by their evaluations.
  const std::uint64_t setUp = DataSetFingerprint(data.Data());
  check("the set-up froze the data set", data.IsFrozen());
  Copy freshA = copy(), freshB = copy();
  check("a copy shares its source's data set (one object)",
        freshA.data->DataSetPointer().get() == data.DataSetPointer().get() &&
            freshB.data->DataSetPointer().get() == data.DataSetPointer().get());
  const Evaluation eA = Evaluate(freshA, pA, cfg), eB = Evaluate(freshB, pB, cfg);
  check("an evaluation at A and one at B leave the shared data set as the set-up made it",
        DataSetFingerprint(freshB.data->Data()) == setUp);
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

  // A copy of a working copy shares no mutable state with its source (EData::Clone is a plain copy): the source
  // evaluated at other parameters after the copy was made moves neither the copy's buffer nor its next evaluation.
  // D is B with every free energy shift moved 1 keV further, so that the source re-energises its points.
  {
    vector_r pD = pB;
    for (unsigned int i = 0; i < pA.size(); i++)
      if (!mp.Parameter(i).IsFixed() && params.Kind((int)i) == ParamKind::EnergyShift) pD[i] += 0.001;
    Copy source = copy();
    Evaluate(source, pB, cfg);
    Copy second;
    second.compound.reset(source.compound->Clone());
    second.data.reset(source.data->Clone());
    Compare("a copy of a copy (evaluated at B) at A = fresh A", Evaluate(second, pA, cfg), eA);
    const std::vector<double> slots = second.results.slot, models = second.results.model;
    Evaluate(source, pD, cfg);
    check("the source evaluated at D leaves the copy's buffer", SameBits(second.results.slot, slots, why) &&
                                                                     SameBits(second.results.model, models, why),
          why);
    Compare("then the copy at A again = fresh A", Evaluate(second, pA, cfg), eA);
  }

  // C: A with every fixed energy-shift parameter (constant and sqrt(E) coefficient) moved by 1 keV (1e-3).  A fixed
  // parameter is still a component of the vector, and pyazr can change it between evaluations: a copy that evaluated
  // C and then A -- where most fixed shifts are 0 again -- must give fresh A's bits (the points must not keep C's
  // shift because the segment's shift is inactive at A).
  vector_r pC = pA;
  int movedFixed = 0;
  for (unsigned int i = 0; i < pA.size(); i++) {
    const ParamKind kind = params.Kind((int)i);
    if (mp.Parameter(i).IsFixed() && (kind == ParamKind::EnergyShift || kind == ParamKind::EnergyShiftSqrt))
      pC[i] += 0.001, movedFixed++;
  }
  if (movedFixed > 0) {
    std::printf("  C moves %d fixed energy-shift parameters\n", movedFixed);
    Copy freshC = copy();
    const Evaluation eC = Evaluate(freshC, pC, cfg);
    check("A and C give different chi^2", eA.segmentChi2[0] != eC.segmentChi2[0]);
    Copy fixedShift = copy();
    Evaluate(fixedShift, pC, cfg);
    Compare("pooled C -> A = fresh A (fixed shifts moved and back)", Evaluate(fixedShift, pA, cfg), eA);
    Compare("pooled C -> A -> C = fresh C", Evaluate(fixedShift, pC, cfg), eC);
  }
  check("every evaluation and hand-over above leaves the shared data set as the set-up made it",
        DataSetFingerprint(data.Data()) == setUp && DataSetFingerprint(freshA.data->Data()) == setUp);
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
