// physics/coulomb at eta = 0 (docs/upgrade/coulomb.md section 10): F, F', G, G' of the default algorithm (coulsm)
// for a neutral pair against mpmath at 40 digits (data/neutral_coulomb.csv from data/make_neutral_coulomb.py:
// F = rho j_l, G = -rho y_l by mpmath's besselj and bessely of order l + 1/2), l = 0-8 and rho = 1e-6-5e4 (small rho
// for each l, rho ~ l, large rho).  Every function within 1e-13 relative (none is near a zero or subnormal on this
// grid), and the Wronskian F'G - FG' = 1 to 1e-13.
//
// What it detects: F and F' from an upward recurrence in l, which is unstable for F at rho < l (the routine's own
// eta = 0 branch, wrong by 1.6e-10 at l = 1 and by 2.9e13 at l = 4 for rho = 7.5e-4: this test fails on it); a
// change of G and G', which climb in their stable direction; a sign or normalisation (the Wronskian).
// Usage: physics_coulomb_neutral_test <path to neutral_coulomb.csv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/coulomb/Coulomb.h"

using namespace azure2::physics::coulomb;

namespace {

struct Row {
  int l;
  double rho;
  double v[4];  // F, F', G, G'
};

std::vector<Row> Read(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    Row r{};
    std::string tok;
    ss >> r.l >> tok;
    r.rho = std::strtod(tok.c_str(), nullptr);
    for (double &x : r.v) {
      ss >> tok;
      x = std::strtod(tok.c_str(), nullptr);
    }
    rows.push_back(r);
  }
  return rows;
}

const int kLMax = 8;
const double kTolerance = 1e-13;  // asked of the repair; 6.1e-15 measured (GCC 11.5, x86-64, -O3)

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s neutral_coulomb.csv\n", argv[0]);
    return 2;
  }
  const std::vector<Row> rows = Read(argv[1]);
  if (rows.size() != 369) {
    std::printf("FAIL: reference file %s incomplete (%zu rows, 369 expected)\n", argv[1], rows.size());
    return 1;
  }
  int failures = 0;
  // Largest relative error per l of F, F', G, G' and |F'G - FG' - 1|, and where the largest of F, F' lies.
  double err[kLMax + 1][5] = {};
  double rho_at[kLMax + 1] = {};
  for (const Row &r : rows) {
    const CoulombWaves w = coulomb(r.l, 0.0, r.rho);
    if (w.status != Status::ok) {
      std::printf("FAIL: l=%d rho=%g: status not ok\n", r.l, r.rho);
      ++failures;
      continue;
    }
    const double got[4] = {w.F, w.Fp, w.G, w.Gp};
    double e[5];
    for (int i = 0; i < 4; ++i) e[i] = std::fabs(got[i] - r.v[i]) / std::fabs(r.v[i]);
    e[4] = std::fabs(w.Fp * w.G - w.F * w.Gp - 1.0);
    for (int i = 0; i < 5; ++i) {
      if (std::isnan(e[i])) e[i] = INFINITY;
      if (i < 2 && e[i] > err[r.l][0] && e[i] > err[r.l][1]) rho_at[r.l] = r.rho;
      err[r.l][i] = std::fmax(err[r.l][i], e[i]);
    }
    if (!(e[0] <= kTolerance && e[1] <= kTolerance && e[2] <= kTolerance && e[3] <= kTolerance &&
          e[4] <= kTolerance)) {
      std::printf("FAIL: l=%d rho=%g: F %.1e, F' %.1e, G %.1e, G' %.1e, F'G - FG' - 1 %.1e\n", r.l, r.rho, e[0], e[1],
                  e[2], e[3], e[4]);
      ++failures;
    }
  }
  std::printf("| l | F | F' | G | G' | F'G - FG' - 1 | rho of the largest error of F, F' |\n"
              "|---|---|---|---|---|---|---|\n");
  for (int l = 0; l <= kLMax; ++l)
    std::printf("| %d | %.1e | %.1e | %.1e | %.1e | %.1e | %g |\n", l, err[l][0], err[l][1], err[l][2], err[l][3],
                err[l][4], rho_at[l]);
  if (failures) std::printf("%d failure(s) of 369 points (tolerance %.0e)\n", failures, kTolerance);
  else std::printf("physics/coulomb at eta = 0: all checks passed\n");
  return failures ? 1 : 0;
}
