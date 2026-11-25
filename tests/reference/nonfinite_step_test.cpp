/*!
 * A non-finite parameter vector is a failed evaluation, never an abort (2026-11-20).
 *
 * A fit step that went to NaN (an LM step of a fit with a potential background) reached the channel functions with
 * eta = rho = NaN; the Coulomb functions throw there ("Coulomb functions failed at l = 1, eta = nan, rho = nan"), the
 * hybrid pair's fallback to the point-Coulomb functions threw again, and the program ended with an uncaught
 * std::runtime_error.  A pair without a potential ended the same way.  Now AZURECalc refuses such a vector before
 * evaluating it (AZURECalc::FiniteVector): operator() returns the chi-squared ceiling the minimizers already get for a
 * non-finite chi-squared, Chi2Value NaN (LM: not a decrease), ResidualJacobian and ResidualsOnly false (the fitters'
 * failure path); and a parameter file with a non-finite value is refused when it is read.
 *
 * Models (synthetic, nothing of a user's project): tests/potential_background/models/fit_level.azr (3He + 4He, the
 * Gaussian potential VD with the s-wave interior background, a free 1/2+ level on top, phase-shift data) and the same
 * project without its <potentials> section (a pair without a potential).  Checked for each, with NaN, +inf and -inf
 * in the level energy and in its width:
 *   - operator() = 1e100, no exception; Chi2Value is NaN; ResidualJacobian and ResidualsOnly return false;
 *   - at the start vector all four evaluate (finite, true) and Chi2Value = operator();
 *   - a parameter file with "nan" or "inf" as a value is refused (ReadUserParameters false), one with a finite
 *     value is applied.
 * Without the guard every non-finite case throws (the test then ends with "exception").
 *
 * Run from the directory CMake prepares (data/ of the model, output/, checks/):
 *   nonfinite_step_test <path to fit_level.azr>   (ctest: nonfinite_step)
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../upgrade/SetUpForTests.h"
#include "AZURECalc.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "WorkingCopy.h"

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

/// The project without its <potentials> section.
bool WithoutPotential(const std::string &source, const std::string &target) {
  std::ifstream in(source.c_str());
  if (!in) return false;
  std::ofstream out(target.c_str());
  std::string line;
  bool skip = false, removed = false;
  while (std::getline(in, line)) {
    if (line == "<potentials>") skip = removed = true;
    if (!skip) out << line << "\n";
    if (line == "</potentials>") skip = false;
  }
  return removed;
}

struct Model {
  Config cfg;
  std::unique_ptr<const CompoundSystem> owner;
  const CompoundSystem *compound = nullptr;
  LevelState state;
  LevelInputs inputs;
  FitterStates fitterStates;
  WorkingCopy *data = nullptr;
  AZUREParams params;
  explicit Model(const std::string &azr) : cfg(std::cout) { cfg.configfile = azr; }
  ~Model() { delete data; }
  bool Load() {
    if (cfg.ReadConfigFile() < 0) return false;
    data = new WorkingCopy();
    CompoundSetUp setUp = test_setup::WithData(cfg, *data);
    if (setUp.refusal != CompoundSetUp::Refusal::kNone) return false;
    owner = std::move(setUp.compound);
    compound = owner.get();
    state = std::move(setUp.state);
    inputs = std::move(setUp.inputs);
    compound->FillMnParams(params, state, inputs, &cfg);
    data->FillMnParams(params);
    return data->Initialize(compound, cfg) != -1;
  }
};

void Variant(const std::string &label, const std::string &azr) {
  std::printf("== %s\n", label.c_str());
  Model m(azr);
  if (!m.Load()) {
    check(label + ": project loads", false, azr);
    return;
  }
  ROOT::Minuit2::MnUserParameters &mp = m.params.GetMinuitParams();
  const int n = (int)mp.Params().size();
  vector_r p(n);
  std::vector<bool> fixed(n);
  std::vector<int> freeIdx;
  for (int i = 0; i < n; i++) {
    p[i] = mp.Value(i);
    fixed[i] = mp.Parameter(i).IsFixed();
    if (!fixed[i]) freeIdx.push_back(i);
  }
  check(label + ": the level's energy and width are free", freeIdx.size() == 2, std::to_string(freeIdx.size()));
  if (freeIdx.size() != 2) return;

  azure2::Results results;
  AZURECalc calc(m.data, results, m.compound, m.state, m.inputs, m.fitterStates, m.cfg);
  calc.SetErrorDef(1.0);
  calc.SetFixedMask(fixed);
  m.data->SetFit(true);

  try {
    const double full = calc(p), value = calc.Chi2Value(p);
    vector_r r, jac, ro;
    std::vector<int> p2f;
    const bool rj = calc.ResidualJacobian(p, r, jac, p2f), rs = calc.ResidualsOnly(p, ro);
    check(label + ": the start vector evaluates (chi2 " + std::to_string(full) + ")",
          std::isfinite(full) && full < 1e99 && std::fabs(value - full) <= 1e-10 * full && rj && rs);
  } catch (const std::exception &e) {
    check(label + ": the start vector evaluates", false, std::string("exception: ") + e.what());
  }

  const double bad[] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
  for (int k : freeIdx)
    for (double b : bad) {
      vector_r q = p;
      q[k] = b;
      std::ostringstream what;
      what << label << ": " << mp.GetName(k) << " = " << b;
      try {
        const double full = calc(q), value = calc.Chi2Value(q);
        vector_r r, jac, ro;
        std::vector<int> p2f;
        const bool rj = calc.ResidualJacobian(q, r, jac, p2f), rs = calc.ResidualsOnly(q, ro);
        check(what.str() + ": a failed evaluation (ceiling, NaN, false, false)",
              full == 1e100 && std::isnan(value) && !rj && !rs,
              std::to_string(full) + " " + std::to_string(value) + " " + std::to_string(rj) + " " + std::to_string(rs));
      } catch (const std::exception &e) {
        check(what.str() + ": a failed evaluation", false, std::string("exception: ") + e.what());
      }
    }

  // a parameter file: "nan" and "inf" refused, a finite value applied
  for (const char *value : {"nan", "inf", "12.5"}) {
    const std::string file = label + "_start.par";
    {
      std::ofstream out(file.c_str());
      out << "#parametrization 1.0 0.0\n" << mp.GetName(freeIdx[0]) << " " << value << " 1.0\n";
    }
    std::ostringstream messages;
    std::streambuf *saved = std::cout.rdbuf(messages.rdbuf());  // ReadUserParameters says why on std::cout
    const bool applied = m.params.ReadUserParameters(file);
    std::cout.rdbuf(saved);
    const bool finite = std::string(value) == "12.5";
    check(label + ": a parameter file with the value " + value + (finite ? " is applied" : " is refused"),
          applied == finite && (finite || messages.str().find("not a finite number") != std::string::npos),
          messages.str());
  }
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: nonfinite_step_test <fit_level.azr>\n");
    return 2;
  }
  Variant("potential background", argv[1]);
  if (WithoutPotential(argv[1], "no_potential.azr"))
    Variant("no potential", "no_potential.azr");
  else
    check("the project without its potential is written", false);
  std::printf(failures ? "nonfinite_step: %d FAILED\n" : "nonfinite_step: all passed\n", failures);
  return failures ? 1 : 0;
}
