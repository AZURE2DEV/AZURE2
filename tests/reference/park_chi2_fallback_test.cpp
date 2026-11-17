/*!
 * The Park J > 0 penalty in AZURECalc::Chi2Value, and so in the finite-
 * difference fallback of AZURECalc::Gradient, on a project whose gradient takes
 * that fallback.
 *
 * Under --use-park the objective MIGRAD minimizes (AZURECalc::operator()) is
 * the data chi-squared plus sum_lambda (J_lambda / 1e-3)^2 over the levels with
 * J < 0.  The analytic gradient adds that penalty's derivative
 * (ParkTerm's AddGradient), but the side-effect-free Chi2Value -- which the
 * gradient differentiates where the analytic adjoint does not apply (all of it
 * when the adjoint bails), and which LM and the GSL trust region take as their
 * cost -- left the penalty out.  Beyond the wall the fallback gradient was
 * then the derivative of a different function than the one minimized, by
 * ~1e6 per unit of a reduced width here.
 *
 * Model: tests/identical_pp_res with the 0+ level's width raised to 51 MeV,
 * past the bound (J = -0.54, penalty 2.9e5, as tests/park_formalism section 4);
 * free: that width and the 2- level's first width (J > 0).  Two variants:
 *   A  as is;
 *   B  a 0.1 keV Gaussian convolution on the analyzing-power segment (a fold
 *      of a ratio: the adjoint goes through the fold's weights).
 * The analytic adjoint applies to both (a convolved A_y made it bail until
 * dev's target-integrated A_y adjoint, 2e726db; the automatic-application
 * tolerance that made B bail after that is deleted: the shared mesh folds
 * every point): no construct of a p+p project is left that the adjoint does
 * not cover.  Variant C takes the fallback:
 *   C  tests/13N_capture_ay (12C(p,g)13N, a capture analyzing power) with a
 *      2 keV Gaussian beam spread on the analyzing power: a folded capture A_y,
 *      which the adjoint refuses (GradTargetEffectAdjoint, src/AZUREGrad.cpp:
 *      the attenuation folded into P_k has no counterpart for the P_k^1 of the
 *      numerator), so Gradient() differentiates Chi2Value by central
 *      differences; the 1/2- background level's proton width raised past the
 *      bound and free, the 1/2+ resonance's proton width free.
 * Variant D reaches another gate of the adjoint that input reaches (AMatrixFunc::PointTBar), without Park's
 * formalism (the fallback's gradient against central differences of the objective):
 *   D  variant A's project with --use-rmc (Reich-Moore capture, --no-brune): the adjoint does not cover it, for
 *      every point.
 * The gates of PointTBar that no project reaches are assertions there (the reason at each).
 *   E  variant A with the wall level's J set to NaN: the value, the rows and the gradient all skip it (the gradient
 *      tested J >= 0 until 2026-11-16 and gave NaN; NanVariant).
 * Checked, at the starting point:
 *   - AccumulateEGammaGradient succeeds for A and B and bails for C;
 *   - Chi2Value(p) = operator()(p), the penalty included (rel. 1e-10);
 *   - Gradient(p) = central differences of operator() (rel. 1e-5), A, B and C
 *     (for C the fallback's own differences of Chi2Value: without the penalty
 *     in Chi2Value they differ by the penalty's derivative);
 *   - A and B agree on the wall width's derivative (rel. 1e-3: the two
 *     models differ by the convolution only).
 *
 * Run from the directory CMake prepares (data/ with both projects' data, output/, checks/):
 *   park_chi2_fallback_test <path to identical_pp_res.azr> <path to 13N_capture_ay.azr>   (ctest: park_chi2_fallback)
 */
#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../upgrade/SetUpForTests.h"
#include "AZURECalc.h"
#include "Objective.h"
#include "ProjectSetup.h"
#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "CaptureAyRequests.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"

namespace {

int failures = 0;
// The proton width (eV) of 13N_capture_ay's 1/2- background level in variant C: past Park's bound.
const char *kCaptureWall = "1e+09";

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
// raised past Park's bound; with `convolve`, a beam spread on data set 4 (ay_60).
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
    // The records the converter writes for the old <targetInt> line 1 "4" 10 1 0.0001 ...
    if (line == "dataset ay_60" && convolve) out << "    beam           spread_ay_60\n    quadrature     points 10\n";
    if (line == "<targets>" && convolve) out << "beam spread_ay_60\n    spread         0.0001 MeV lab\n\n";
  }
  return wallDone;
}

// tests/13N_capture_ay with every energy and width fixed but two widths, the 1/2- background level's proton width
// (20 MeV, 5e6 eV) raised past Park's bound, and a beam spread on the capture analyzing power (data set capture_ay).
bool writeCaptureProject(const std::string &source, const std::string &target) {
  std::ifstream in(source.c_str());
  if (!in) return false;
  std::ofstream out(target.c_str());
  std::string line;
  bool inLevels = false, wallDone = false, secondDone = false, inAy = false;
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
        f[3] = "1";   // energy fixed
        f[10] = "1";  // width fixed
        if (!wallDone && f[0] == "0.5" && f[1] == "-1" && f[2] == "20" && f[11] == "5e+06") {
          f[11] = kCaptureWall;  // past the bound, free
          f[10] = "0";
          wallDone = true;
        } else if (!secondDone && f[2] == "2.3682" && f[11] == "34000") {  // the 1/2+ resonance's proton width: free
          f[10] = "0";
          secondDone = true;
        }
        for (size_t i = 0; i < f.size(); i++) out << (i ? " " : "") << f[i];
        out << "\n";
        continue;
      }
    }
    if (line == "dataset capture_ay") inAy = true;
    out << line << "\n";
    if (inAy && line.find("    norm ") == 0) {
      out << "    beam           spread_ay\n    quadrature     points 10\n";
      inAy = false;
    }
    if (line == "<targets>") out << "beam spread_ay\n    spread         0.002 MeV lab\n\n";
  }
  return wallDone && secondDone;
}

struct Model {
  Config cfg;
  std::unique_ptr<const CompoundSystem> owner;  // made by the set-up (test_setup::WithData)
  const CompoundSystem *compound = nullptr;
  LevelState state;  // the model's levels (the run's, as AZUREMain holds them)
  LevelInputs inputs;
  FitterStates fitterStates;
  WorkingCopy *data = nullptr;
  AZUREParams params;
  // `externalCapture`: the project has external capture (as the CLI's checkExternalCapture sets it); `mask`: the
  // formalism flags (Park's by default).
  Model(const std::string &azr, bool externalCapture,
        unsigned mask = Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM)
      : cfg(std::cout) {
    cfg.configfile = azr;
    cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
    cfg.paramMask |= mask;
    if (externalCapture) cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
  }
  ~Model() {
    delete data;
  }
  bool Load() {
    if (cfg.ReadConfigFile() < 0) return false;
    data = new WorkingCopy();
    CompoundSetUp setUp = test_setup::WithData(cfg, *data, test_setup::CaptureAy::kRequested);
    if (setUp.refusal != CompoundSetUp::Refusal::kNone) {
      return false;
    }
    owner = std::move(setUp.compound);
    compound = owner.get();
    state = std::move(setUp.state);
    inputs = std::move(setUp.inputs);
    compound->FillMnParams(params, state, inputs, &cfg);
    data->FillMnParams(params);
    return data->Initialize(compound, cfg) != -1;
  }
};

// Whether the analytic energy/width adjoint applies at p (AZURECalc::Gradient
// falls back to finite differences of Chi2Value where it does not).
bool AdjointApplies(Model &m, const vector_r &p) {
  const CompoundSystem *lc = m.compound;
  LevelState ls(m.state);  // a copy, as a fit's evaluation takes one
  WorkingCopy *ld = m.data->Clone();
  lc->FillCompoundFromParams(p, ls);
  ld->FillNormsFromParams(p);
  azure2::Results results;
  ld->FillEnergyShiftsFromParams(p, results, lc, &ls, &m.inputs, &m.cfg);
  lc->CalcShiftFunctions(m.cfg, ls);
  vector_matrix_r sd = BuildShiftDerivTable(lc, ls, m.cfg);
  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, std::vector<bool>());
  GradAccum accum;
  accum.Init(lc, ls);
  FitBarFn one = [](const Segment *, int, int, double) { return 1.0; };
  const bool ok = AccumulateEGammaGradient(lc, ls, ld, results, m.cfg, pmap, &sd, one, accum);
  delete ld;
  return ok;
}

// Runs the checks on one variant; returns the gradient of the wall width (the free parameter `wallIndex`, the first
// free one).
double Variant(const std::string &label, const std::string &azr, bool expectAdjoint, int wallIndex,
               bool externalCapture) {
  std::printf("== %s\n", label.c_str());
  Model m(azr, externalCapture);
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
  check(label + ": two free parameters (the wall's width first)", freeIdx.size() == 2 && freeIdx[0] == wallIndex,
        std::to_string(freeIdx.size()) + (freeIdx.empty() ? "" : ", first " + std::to_string(freeIdx[0])));
  if (freeIdx.size() != 2) return NAN;

  const bool adjoint = AdjointApplies(m, p);
  check(label + (expectAdjoint ? ": analytic adjoint applies" : ": analytic adjoint bails (fallback path)"),
        adjoint == expectAdjoint);

  azure2::Results results;
  AZURECalc calc(m.data, results, m.compound, m.state, m.inputs, m.fitterStates, m.cfg);
  calc.SetErrorDef(1.0);
  calc.SetFixedMask(fixed);
  m.data->SetFit(true);

  const double full = calc(p);
  const double value = calc.Chi2Value(p);
  m.compound->FillCompoundFromParams(p, m.state);
  m.compound->CalcShiftFunctions(m.cfg, m.state);
  double penalty = 0.0;  // the objective's Park term alone
  for (const auto &term : Objective::Engine().Terms())
    if (term->Name() == "Park J > 0") term->Add(ObjectiveContext{m.compound, m.state, m.inputs, m.data, results, p, m.cfg}, penalty);
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

// A variant without Park's formalism whose adjoint bails at a gate of PointTBar: the gradient of every free parameter
// against central differences of the objective.
void FallbackVariant(const std::string &label, const std::string &azr, unsigned mask) {
  std::printf("== %s\n", label.c_str());
  Model m(azr, false, mask);
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
  check(label + ": free parameters (" + std::to_string(freeIdx.size()) + ")", !freeIdx.empty());
  check(label + ": analytic adjoint bails (fallback path)", !AdjointApplies(m, p));
  azure2::Results results;
  AZURECalc calc(m.data, results, m.compound, m.state, m.inputs, m.fitterStates, m.cfg);
  calc.SetErrorDef(1.0);
  calc.SetFixedMask(fixed);
  m.data->SetFit(true);
  const double full = calc(p);
  check(label + ": objective finite (" + num(full) + ")", std::isfinite(full) && full < 1e99);
  std::vector<double> g = calc.Gradient(p);
  for (int idx : freeIdx) {
    const double h = 1.0e-6 * (std::fabs(p[idx]) + 1.0);
    vector_r pp = p, pm = p;
    pp[idx] += h;
    pm[idx] -= h;
    const double fd = (calc(pp) - calc(pm)) / (2.0 * h);
    const double rel = std::fabs(g[idx] - fd) / (std::fabs(fd) + 1.0);
    check(label + ": d/d" + mp.GetName(idx) + " = differences of the objective (" + num(g[idx]) + " vs " + num(fd) +
              ")",
          rel < 1e-5, "rel " + num(rel));
  }
}

// A level whose overlap J is not a number (a NaN width or shift-function derivative) is not on the wall: the value
// (ParkTerm), the rows (its AppendRows) and the gradient (its AddGradient) skip it alike.  Until 2026-11-16 the
// gradient tested J >= 0, which a NaN fails, so it added 2 J/s^2 (...) = NaN to the level's energy and widths while
// the value and the rows added nothing.  On variant A: the wall level's J set to NaN after the shift functions.
void NanVariant(const std::string &azr) {
  const std::string label = "E: a level with J = NaN";
  std::printf("== %s\n", label.c_str());
  Model m(azr, false);
  if (!m.Load()) {
    check(label + ": project loads", false, azr);
    return;
  }
  const vector_r p = m.params.GetMinuitParams().Params();
  m.compound->FillCompoundFromParams(p, m.state);
  m.compound->CalcShiftFunctions(m.cfg, m.state);
  const ParamIndexMap pmap = BuildParamIndexMap(m.compound, m.data, std::vector<bool>());
  auto gradient = [&]() {
    GradAccum accum;
    accum.Init(m.compound, m.state);
    azure2::Results results;
    for (const auto &term : Objective::Engine().Terms())
      if (term->Name() == "Park J > 0")
        term->AddGradient(ObjectiveContext{m.compound, m.state, m.inputs, m.data, results, p, m.cfg}, pmap, nullptr, accum);
    vector_r g(pmap.NumFull(), 0.0);
    accum.Scatter(pmap, g);
    return g;
  };
  auto wallValue = [&]() {
    azure2::Results results;
    double v = 0.0;
    for (const auto &term : Objective::Engine().Terms())
      if (term->Name() == "Park J > 0") term->Add(ObjectiveContext{m.compound, m.state, m.inputs, m.data, results, p, m.cfg}, v);
    return v;
  };
  // The level beyond the wall (J < 0): there must be one, or the check below tests nothing.
  int wallJ = 0, wallLa = 0;
  for (int j = 1; j <= m.compound->NumJGroups() && !wallJ; j++)
    for (int la = 1; la <= m.compound->GetJGroup(j)->NumLevels(); la++)
      if (m.state.ParkNorm(j, la) < 0.0) {
        wallJ = j;
        wallLa = la;
        break;
      }
  check(label + ": a level beyond the wall at the start", wallJ > 0);
  if (!wallJ) return;
  double before = 0.0;
  for (double x : gradient()) before = std::max(before, std::fabs(x));
  check(label + ": with J < 0 the wall has a gradient (largest " + num(before) + ")", before > 0.0 && std::isfinite(before));
  m.state.SetParkNorm(wallJ, wallLa, NAN);
  const vector_r g = gradient();
  bool zero = true;
  for (double x : g) zero = zero && x == 0.0;
  check(label + ": no gradient from it (every component 0, none NaN)", zero);
  check(label + ": no value from it (" + num(wallValue()) + ")", wallValue() == 0.0);
  vector_r rows;
  {
    azure2::Results results;
    Objective::Engine().AppendModelRows(ObjectiveContext{m.compound, m.state, m.inputs, m.data, results, p, m.cfg}, &pmap, 0,
                                        rows, nullptr);
  }
  bool rowsZero = !rows.empty();
  for (double r : rows) rowsZero = rowsZero && r == 0.0;
  check(label + ": its row is 0 (" + std::to_string(rows.size()) + " rows)", rowsZero);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s path/to/identical_pp_res.azr path/to/13N_capture_ay.azr\n", argv[0]);
    return 2;
  }
  if (!writeProject(argv[1], "park_a.azr", false) || !writeProject(argv[1], "park_b.azr", true) ||
      !writeCaptureProject(argv[2], "park_c.azr")) {
    std::fprintf(stderr, "cannot read or adapt %s or %s\n", argv[1], argv[2]);
    return 2;
  }
  const double gA = Variant("A: analytic gradient", "park_a.azr", true, 1, false);
  const double gB = Variant("B: a convolved analyzing power", "park_b.azr", true, 1, false);
  Variant("C: a folded capture analyzing power (the fallback)", "park_c.azr", false, 4, true);
  FallbackVariant("D: Reich-Moore capture (--use-rmc)", "park_a.azr", Config::USE_RMC_FORMALISM);
  NanVariant("park_a.azr");
  check("A and B agree on the wall width's derivative (" + num(gA) + " vs " + num(gB) + ")",
        std::fabs(gA - gB) <= 1e-3 * std::fabs(gA));
  std::printf("%s\n", failures ? "FAIL" : "PASS");
  return failures ? 1 : 0;
}
