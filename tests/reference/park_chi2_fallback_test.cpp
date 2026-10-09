/*!
 * The Park J > 0 penalty in AZURECalc::Chi2Value, and so in the finite-
 * difference fallback of AZURECalc::Gradient.
 *
 * Under --use-park the objective MIGRAD minimizes (AZURECalc::operator()) is
 * the data chi-squared plus sum_lambda (J_lambda / 1e-3)^2 over the levels with
 * J < 0.  The analytic gradient adds that penalty's derivative
 * (AddParkPenaltyGradient), but the side-effect-free Chi2Value -- which the
 * gradient differentiates where the analytic adjoint does not apply (all of it
 * when the adjoint bails), and which LM and the GSL trust region take as their
 * cost -- left the penalty out.  Beyond the wall the fallback gradient was
 * then the derivative of a different function than the one minimized, by
 * ~1e6 per unit of a reduced width here.
 *
 * Model: tests/identical_pp_res with the 0+ level's width raised to 51 MeV,
 * past the bound (J = -0.54, penalty 2.9e5, as tests/park_formalism section 4);
 * free: that width and the 2- level's first width (J > 0).  Two variants:
 *   A  as is -- the analytic adjoint applies;
 *   B  a 0.1 keV Gaussian convolution on the analyzing-power segment: sub-
 *      points under A_y are not handled by the adjoint, which bails, and
 *      Gradient() is finite differences of Chi2Value throughout.
 * Checked, at the starting point:
 *   - AccumulateEGammaGradient succeeds for A and bails for B (the fallback is
 *     really taken);
 *   - Chi2Value(p) = operator()(p), the penalty included (rel. 1e-10);
 *   - Gradient(p) = central differences of operator() (rel. 1e-5), A and B;
 *   - A and B agree on the wall width's derivative (rel. 1e-3: the two
 *     models differ by the convolution only).
 *
 * Run from the directory CMake prepares (data/, output/, checks/):
 *   park_chi2_fallback_test <path to identical_pp_res.azr>   (ctest: park_chi2_fallback)
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AZURECalc.h"
#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"

Config *g_config = nullptr;

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

std::string num(double x) {
  std::ostringstream s;
  s.precision(10);
  s << x;
  return s.str();
}

// The project with every energy and width fixed but two widths, the 0+ one
// raised past Park's bound; with `convolve`, a target-effect line on segment 4.
bool writeProject(const std::string &source, const std::string &target, bool convolve) {
  std::ifstream in(source.c_str());
  if (!in) return false;
  std::ofstream out(target.c_str());
  std::string line;
  bool inLevels = false, wallDone = false;
  int levelLine = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
    if (line == "<levels>") { inLevels = true; out << line << "\n"; continue; }
    if (line == "</levels>") inLevels = false;
    if (inLevels) {
      std::istringstream fields(line);
      std::vector<std::string> f;
      std::string token;
      while (fields >> token) f.push_back(token);
      if (f.size() > 12) {
        levelLine++;
        f[3] = "1";   // energy fixed
        f[10] = "1";  // width fixed
        if (!wallDone && f[11] == "300000.0") {  // 0+: past the bound, free
          f[11] = "51000000";
          f[10] = "0";
          wallDone = true;
        } else if (levelLine == 4) {             // 2-, first channel: free
          f[10] = "0";
        }
        for (size_t i = 0; i < f.size(); i++) out << (i ? " " : "") << f[i];
        out << "\n";
        continue;
      }
    }
    out << line << "\n";
    if (line == "<targetInt>" && convolve) out << "1 \"4\" 10 1 0.0001 0 0 \"\" 0 0 0 0 \"\" 0\n";
  }
  return wallDone;
}

struct Model {
  Config cfg;
  CNuc *compound = nullptr;
  EData *data = nullptr;
  AZUREParams params;
  explicit Model(const std::string &azr) : cfg(std::cout) {
    cfg.configfile = azr;
    cfg.paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  }
  ~Model() {
    delete compound;
    delete data;
  }
  bool Load() {
    g_config = &cfg;
    if (cfg.ReadConfigFile() < 0) return false;
    compound = new CNuc();
    data = new EData();
    if (compound->Fill(cfg) == -1 || data->Fill(cfg, compound) == -1) return false;
    compound->Initialize(cfg);
    compound->FillMnParams(params.GetMinuitParams(), &cfg);
    data->FillMnParams(params.GetMinuitParams());
    return data->Initialize(compound, cfg) != -1;
  }
};

// Whether the analytic energy/width adjoint applies at p (AZURECalc::Gradient
// falls back to finite differences of Chi2Value where it does not).
bool AdjointApplies(Model &m, const vector_r &p) {
  CNuc *lc = m.compound->Clone();
  EData *ld = m.data->Clone();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &m.cfg);
  lc->CalcShiftFunctions(m.cfg);
  vector_matrix_r sd = BuildShiftDerivTable(lc, m.cfg);
  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, std::vector<bool>());
  GradAccum accum;
  accum.Init(lc);
  FitBarFn one = [](ESegment *, int, int, double) { return 1.0; };
  const bool ok = AccumulateEGammaGradient(lc, ld, m.cfg, pmap, &sd, one, accum);
  delete lc;
  delete ld;
  return ok;
}

// Runs the checks on one variant; returns the gradient of the wall width.
double Variant(const std::string &label, const std::string &azr, bool expectAdjoint) {
  std::printf("== %s\n", label.c_str());
  Model m(azr);
  if (!m.Load()) {
    check(label + ": project loads", false, azr);
    return NAN;
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
  check(label + ": two free parameters (the 0+ width first)", freeIdx.size() == 2 && freeIdx[0] == 1,
        std::to_string(freeIdx.size()));
  if (freeIdx.size() != 2) return NAN;

  const bool adjoint = AdjointApplies(m, p);
  check(label + (expectAdjoint ? ": analytic adjoint applies" : ": analytic adjoint bails (fallback path)"),
        adjoint == expectAdjoint);

  AZURECalc calc(m.data, m.compound, m.cfg);
  calc.SetErrorDef(1.0);
  calc.SetFixedMask(fixed);
  m.data->SetFit(true);

  const double full = calc(p);
  const double value = calc.Chi2Value(p);
  m.compound->FillCompoundFromParams(p);
  m.compound->CalcShiftFunctions(m.cfg);
  const double penalty = m.compound->ParkNormPenalty();
  check(label + ": starting point beyond the wall (Park penalty " + num(penalty) + ")", penalty > 1.0e3);
  check(label + ": Chi2Value = operator() (" + num(value) + " vs " + num(full) + ")",
        std::fabs(value - full) <= 1e-10 * std::fabs(full), "difference " + num(full - value));

  std::vector<double> g = calc.Gradient(p);
  double wallGrad = NAN;
  for (int k = 0; k < 2; k++) {
    const int idx = freeIdx[k];
    const double h = 1.0e-6 * (std::fabs(p[idx]) + 1.0);
    vector_r pp = p, pm = p;
    pp[idx] += h;
    pm[idx] -= h;
    const double fd = (calc(pp) - calc(pm)) / (2.0 * h);
    const double rel = std::fabs(g[idx] - fd) / (std::fabs(fd) + 1.0);
    check(label + ": d/d" + mp.GetName(idx) + " = differences of the objective (" + num(g[idx]) +
              " vs " + num(fd) + ")",
          rel < 1e-5, "rel " + num(rel));
    if (k == 0) wallGrad = g[idx];
  }
  return wallGrad;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s path/to/identical_pp_res.azr\n", argv[0]);
    return 2;
  }
  if (!writeProject(argv[1], "park_a.azr", false) || !writeProject(argv[1], "park_b.azr", true)) {
    std::fprintf(stderr, "cannot read or adapt %s\n", argv[1]);
    return 2;
  }
  const double gA = Variant("A: analytic gradient", "park_a.azr", true);
  const double gB = Variant("B: finite-difference fallback", "park_b.azr", false);
  check("A and B agree on the wall width's derivative (" + num(gA) + " vs " + num(gB) + ")",
        std::fabs(gA - gB) <= 1e-3 * std::fabs(gA));
  std::printf("%s\n", failures ? "FAIL" : "PASS");
  return failures ? 1 : 0;
}
