// The results of an evaluation belong to its caller (include/Results.h, docs/upgrade/data_model.md S3): a working copy
// of the data set -- and every point copied with it -- writes only into the buffer its own caller hands it, never into
// the buffer of the data set it was copied from or of another copy.  Model: tests/upgrade/T2 (sub-point folds whose
// grids follow the level energies, energy shifts, mapped points, SUM and RATIO components, a data covariance, THM).
// Two parameter vectors: p0, the project's start, and p1, every free level energy moved by 2 keV and every free width
// by 1 %, which re-places sub-point grids, so the copies' layouts are not all the same.
//
//   reference  RA, RB: a fresh copy of the master evaluated alone at p0 and at p1;
//   1. X, a copy of the master, evaluated at p0; Y, a copy of X made after that (its points carry X's offsets),
//      evaluated at p1 into its own buffer: X's buffer still equals RA bit for bit, Y's equals RB;
//   2. X at p1 and Y at p0 at the same time, one per thread: X's buffer equals RB, Y's equals RA;
//   3. X evaluated again at p0 into the buffer it last used at p1: equals RA (a buffer laid out for another layout
//      is laid out again, its values moved, never read at stale offsets);
//   4. every model value (the per-point loop of the objective) equal to the reference's, bit for bit.
// The test sets AZURE_POISON_RESULTS=1, so every buffer is NaN on entry to each evaluation and holds exactly what that
// evaluation wrote: a slot no evaluation writes is NaN everywhere, and buffers can be compared whole.
// Shown to fail (2026-11-10) with EPoint::Rec writing into the first buffer it ever saw instead of the one it is
// given -- a buffer shared behind the caller's back, as a point's pointer into its data set's buffer was before S3 (b).
//
//   results_owner_test <T2.azr>   (run in a directory holding T2's data/, cov_diff.txt and prior_cov.txt)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "ModelPass.h"
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

// One evaluation, as AZURECalc::Chi2On does it: the hand-over into `results`, then every point's model value
// (Results::model).
std::vector<double> Evaluate(CNuc *lc, EData *ld, azure2::Results &results, const vector_r &p, const Config &cfg) {
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, results, ld, lc, &cfg);
  if (cfg.paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(cfg, false);
  azure2::EvaluateModel(*ld, lc, cfg, results, azure2::LeafValues::Calculate);
  return results.model;
}

struct Copy {
  std::unique_ptr<CNuc> compound;
  std::unique_ptr<EData> data;
  azure2::Results results;
  std::vector<double> models;
};

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: results_owner_test <T2.azr>\n");
    return 2;
  }
  setenv("AZURE_POISON_RESULTS", "1", 1);  // read once, at the first evaluation
  check("AZURE_POISON_RESULTS on", azure2::Results::PoisonRequested());
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = argv[1];
  g_config = &cfg;
  CNuc compound;
  EData data;
  AZUREParams params;
  if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1 || data.Fill(cfg, &compound) == -1) {
    check("T2 loads", false, sink.str().substr(0, 400));
    return 1;
  }
  compound.Initialize(cfg);
  compound.FillMnParams(params, &cfg);
  data.FillMnParams(params);
  if (data.Initialize(&compound, cfg) == -1) {
    check("T2 initializes", false, sink.str().substr(0, 400));
    return 1;
  }
  const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const vector_r p0 = mp.Params();
  vector_r p1 = p0;
  int moved = 0;
  for (unsigned int i = 0; i < p0.size(); i++) {
    if (mp.Parameter(i).IsFixed()) continue;
    const std::string &name = mp.GetName(i);
    if (name.rfind("energy_", 0) == 0) p1[i] = p0[i] + 0.002, moved++;
    else if (name.rfind("width_", 0) == 0) p1[i] = p0[i] * 1.01, moved++;
  }
  check("p1 moves " + std::to_string(moved) + " level parameters", moved > 0);

  auto copyOf = [](CNuc *c, EData *d) {
    Copy k;
    k.compound.reset(c->Clone());
    k.data.reset(d->Clone());
    return k;
  };
  std::string why;
  Copy refA = copyOf(&compound, &data), refB = copyOf(&compound, &data);
  refA.models = Evaluate(refA.compound.get(), refA.data.get(), refA.results, p0, cfg);
  refB.models = Evaluate(refB.compound.get(), refB.data.get(), refB.results, p1, cfg);
  std::printf("reference: %zu model values, buffers of %zu (p0) and %zu (p1) slots\n", refA.models.size(),
              refA.results.slot.size(), refB.results.slot.size());
  check("p0 and p1 give different model values", !SameBits(refA.models, refB.models, why));

  // 1. Y copied from X after X's evaluation, then evaluated at other parameters.
  Copy x = copyOf(&compound, &data);
  x.models = Evaluate(x.compound.get(), x.data.get(), x.results, p0, cfg);
  Copy y = copyOf(x.compound.get(), x.data.get());
  y.models = Evaluate(y.compound.get(), y.data.get(), y.results, p1, cfg);
  check("1. X at p0 = reference at p0 (buffer)", SameBits(x.results.slot, refA.results.slot, why), why);
  check("1. Y, copied from X, at p1 = reference at p1 (buffer)", SameBits(y.results.slot, refB.results.slot, why), why);
  check("1. X's buffer untouched by Y's evaluation (= reference at p0)",
        SameBits(x.results.slot, refA.results.slot, why), why);

  // 2. X and Y at the same time, one per thread.
#ifdef _OPENMP
  omp_set_num_threads(2);
#endif
#pragma omp parallel sections
  {
#pragma omp section
    x.models = Evaluate(x.compound.get(), x.data.get(), x.results, p1, cfg);
#pragma omp section
    y.models = Evaluate(y.compound.get(), y.data.get(), y.results, p0, cfg);
  }
  check("2. X at p1, concurrently with Y = reference at p1 (buffer)", SameBits(x.results.slot, refB.results.slot, why),
        why);
  check("2. Y at p0, concurrently with X = reference at p0 (buffer)", SameBits(y.results.slot, refA.results.slot, why),
        why);
  check("2. model values: X = reference at p1, Y = reference at p0",
        SameBits(x.models, refB.models, why) && SameBits(y.models, refA.models, why), why);

  // 3. X back at p0, into the buffer it used at p1.
  x.models = Evaluate(x.compound.get(), x.data.get(), x.results, p0, cfg);
  check("3. X at p0 again, in its buffer of p1 = reference at p0 (buffer)",
        SameBits(x.results.slot, refA.results.slot, why), why);
  check("4. model values of X at p0 = reference", SameBits(x.models, refA.models, why), why);

  std::printf("results_owner: %s (%d failure%s)\n", failures ? "FAIL" : "all checks passed", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
