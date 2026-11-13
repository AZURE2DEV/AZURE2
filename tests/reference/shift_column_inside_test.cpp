/*!
 * The difference step of a shift column stays inside the planned nodes (ShiftedFold.h, question 5).
 *
 * The column of a varied energy shift (the scalar gradient's local shift path in AZURECalc::Gradient, and the residual
 * Jacobian's in FillShiftColumns) moves the points of its data sets by +-h without planning the mesh again.  A point
 * whose shift lay within h of the edge of its reach left its planned nodes, and ShiftedFoldNodes threw inside the
 * parallel refold: the process ended (std::logic_error, terminate).  A varied shift without a prior on a thin layer
 * has a reach of 3 support widths, below h ~ 1 eV up to about 3e14 atoms/cm^2, so a fit with --use-gradient or
 * --use-lm ended at its first gradient; on thicker layers it ended whenever a shift came within h of a window's edge.
 * Now the step is halved until both sides stay inside (InsideShiftStep), or the column takes the path that plans
 * again.
 *
 * Model: tests/thin_layer_limit/models/thin_limit.azr (a 15N layer, every level fixed) with "shift offset 0.015 MeV fit"
 * on its two folded data sets (a yield and a target-averaged cross section, data 0 with errors near the model, so
 * that the residuals are of order 1), at 3e14, 1e15 and 1e16 atoms/cm^2.  The shifts are set to 0.015 MeV + 0.97 of
 * the smallest reach of the points (in lab MeV), i.e. 3% of the reach from the edge of the window the hand-over keeps,
 * closer than h.  Checked at that point, for each density:
 *   1. the gradient (AZURECalc::Gradient) is finite and equals central differences of the objective (rel. 1e-3);
 *   2. the residual Jacobian (ResidualJacobian, the --use-lm columns) equals central differences of the residuals
 *      (ResidualsOnly) (rel. 1e-3 of the column's largest element).
 * The reference differences take a step of 1e-8 MeV.  The fold on fixed nodes is smooth in the shift only piecewise
 * (the moved support's edges cross nodes), so differences of other steps agree to about 2e-4 (measured at 1e16 on the
 * target-averaged cross section): the tolerance.
 * Before the fix the test ends in std::terminate at the first density.
 *
 *   shift_column_inside_test <path to thin_limit.azr>   (ctest: shift_column_inside; run in a scratch directory)
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AZURECalc.h"
#include "ProjectSetup.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CaptureAyRequests.h"
#include "Config.h"
#include "WorkingCopy.h"

namespace {

int failures = 0;
const double kA0 = 4.5e-22, kA1 = -0.6;  // epsilon = a0 x^a1, lab MeV cm^2 per active nucleus (thin_limit.azr)
const double kStart = 0.015;             // MeV lab

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

// thin_limit.azr at areal density `n` with a varied shift on the folded data sets, their data in data/y.dat and
// data/s.dat (0 with errors near the model), the bare data set as in the file.
bool writeProject(const std::string &source, const std::string &target, const std::string &n) {
  std::ifstream in(source.c_str());
  if (!in) return false;
  std::ofstream out(target.c_str());
  std::string line;
  int norms = 0, files = 0;
  while (std::getline(in, line)) {
    if (line.find("areal_density  1e+14 ") != std::string::npos) line.replace(line.find("1e+14"), 5, n);
    if (line == "    file           data/e0.dat" && files < 2) line = files++ ? "    file           data/s.dat" : "    file           data/y.dat";
    out << line << "\n";
    if (line.find("    norm ") == 0 && norms < 2) {
      out << "    shift          offset 0.015 MeV fit\n";
      norms++;
    }
  }
  std::ofstream y("data/y.dat"), s("data/s.dat");
  for (double e : {0.6, 0.7, 0.8}) {
    y << e << " 0 0 " << 1e-24 * std::stod(n) * 3e-9 << "\n";
    s << e << " 0 0 3e-9\n";
  }
  return norms == 2 && files == 2;
}

struct Model {
  Config cfg;
  CNuc *compound = nullptr;
  WorkingCopy *data = nullptr;
  AZUREParams params;
  explicit Model(const std::string &azr) : cfg(std::cout) { cfg.configfile = azr; }
  ~Model() {
    delete compound;
    delete data;
  }
  bool Load() {
    if (cfg.ReadConfigFile() < 0) return false;
    compound = new CNuc();
    data = new WorkingCopy();
    if (compound->Fill(cfg) == -1 || data->Fill(cfg, compound) != 0) return false;
    RequestCaptureAyTerms(*compound, *data);
    RequestDataEntrances(*compound, cfg);
    compound->Initialize(cfg);
    compound->FillMnParams(params, &cfg);
    data->FillMnParams(params);
    return data->Initialize(compound, cfg) != -1;
  }
};

void Density(const std::string &source, const std::string &n) {
  std::printf("== %s atoms/cm^2\n", n.c_str());
  const std::string azr = "shift_" + n + ".azr";
  if (!writeProject(source, azr, n)) {
    check("project written", false, azr);
    return;
  }
  Model m(azr);
  if (!m.Load()) {
    check("project loads", false, azr);
    return;
  }
  ROOT::Minuit2::MnUserParameters &mp = m.params.GetMinuitParams();
  const int np = (int)mp.Params().size();
  vector_r p(np);
  std::vector<bool> fixed(np);
  std::vector<int> freeIdx;
  for (int i = 0; i < np; i++) {
    p[i] = mp.Value(i);
    fixed[i] = mp.Parameter(i).IsFixed();
    if (!fixed[i]) freeIdx.push_back(i);
  }
  check("two free shifts", freeIdx.size() == 2, std::to_string(freeIdx.size()));
  if (freeIdx.size() != 2) return;
  // 0.97 of the smallest reach (the 0.8 MeV point's: 3 x the layer grid's 1.02 T, in lab MeV)
  const double reach = 3.0 * 1.02 * kA0 * std::pow(0.8, kA1) * std::stod(n);
  for (int i : freeIdx) p[i] = kStart + 0.97 * reach;

  azure2::Results results;
  AZURECalc calc(m.data, results, m.compound, m.cfg);
  calc.SetErrorDef(1.0);
  calc.SetFixedMask(fixed);
  m.data->SetFit(true);
  const double chi = calc(p);
  check("objective finite (" + num(chi) + ")", std::isfinite(chi) && chi < 1e99);

  const std::vector<double> g = calc.Gradient(p);
  const double h = 1.0e-8;
  for (int i : freeIdx) {
    vector_r pp = p, pm = p;
    pp[i] += h;
    pm[i] -= h;
    const double fd = (calc(pp) - calc(pm)) / (2.0 * h);
    const double rel = std::fabs(g[i] - fd) / std::max(std::fabs(fd), 1e-300);
    check("gradient d/d" + mp.GetName(i) + " = differences of chi^2 (" + num(g[i]) + " vs " + num(fd) + ")",
          std::isfinite(g[i]) && rel < 1e-3, "rel " + num(rel));
  }

  vector_r r, jac;
  std::vector<int> packedToFull;
  const bool ok = calc.ResidualJacobian(p, r, jac, packedToFull);
  check("residual Jacobian computed", ok && !packedToFull.empty());
  if (!ok) return;
  const size_t cols = packedToFull.size(), rows = r.size();
  for (size_t a = 0; a < cols; a++) {
    const int i = packedToFull[a];
    vector_r rp, rm, pp = p, pm = p;
    pp[i] += h;
    pm[i] -= h;
    if (!calc.ResidualsOnly(pp, rp) || !calc.ResidualsOnly(pm, rm) || rp.size() != rows || rm.size() != rows) {
      check("residuals at +-h", false);
      continue;
    }
    double worst = 0.0, scale = 0.0;
    for (size_t k = 0; k < rows; k++) scale = std::max(scale, std::fabs((rp[k] - rm[k]) / (2.0 * h)));
    for (size_t k = 0; k < rows; k++)
      worst = std::max(worst, std::fabs(jac[k * cols + a] - (rp[k] - rm[k]) / (2.0 * h)));
    check("Jacobian column " + mp.GetName(i) + " = differences of the residuals (largest " + num(scale) + ")",
          scale > 0.0 && worst <= 1e-3 * scale, "worst " + num(worst));
  }
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s path/to/thin_limit.azr\n", argv[0]);
    return 2;
  }
  for (const char *n : {"3.162278e+14", "1e+15", "1e+16"}) Density(argv[1], n);
  std::printf("%s\n", failures ? "FAIL" : "PASS");
  return failures ? 1 : 0;
}
