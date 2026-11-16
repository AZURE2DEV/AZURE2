// inference_park_wall: Park's wall J > 0 (T6) on arrays (azure2/inference/ParkWall.h).  One level of two particle
// channels whose S'_c(E) = a_c + b_c E is linear (S''_c = b_c), so that J(gamma, E) = 1 - sum_c gamma_c^2 S'_c(E) is
// known in closed form: J by hand; the value, the gradient's form and the row's form against central differences in
// gamma and E; their Gauss-Newton product against the gradient; J >= 0 and a NaN J give no value, no gradient and a
// zero row.
#include "azure2/inference/ParkWall.h"

#include <cmath>
#include <limits>

#include "check.h"

using namespace azure2::inference;

namespace {
const double a[2] = {0.5, 0.4};    // S'_c at E = 0, MeV^-1
const double b[2] = {0.1, -0.05};  // S''_c, MeV^-2
const double kScale = 1.0e-3;

// The level at amplitudes g (MeV^1/2) and energy E (MeV), as the engine's gather fills it.
ParkLevel LevelAt(const double g[2], double E) {
  ParkLevel l;
  l.group = 2;
  l.level = 3;
  l.channels = {1, 3};
  l.J = 1.0;
  for (int c = 0; c < 2; c++) {
    l.gamma.push_back(g[c]);
    l.dS.push_back(a[c] + b[c] * E);
    l.d2S.push_back(b[c]);
    l.J -= g[c] * g[c] * (a[c] + b[c] * E);
  }
  return l;
}
double Value(const double g[2], double E) { return ParkWallChi2({LevelAt(g, E)}, kScale); }
double Row(const double g[2], double E) {
  double r;
  std::vector<double> dg;
  double dE;
  ParkWallRow(LevelAt(g, E), kScale, r, dg, dE);
  return r;
}
bool Close(double x, double y, double rel) { return std::fabs(x - y) <= rel * (std::fabs(x) + std::fabs(y)); }
}  // namespace

int main() {
  // gamma = (1.2, 0.8), E = 2: S' = (0.7, 0.3), J = 1 - (1.44 0.7 + 0.64 0.3) = -0.2; value (0.2/1e-3)^2 = 40000.
  const double g[2] = {1.2, 0.8}, E = 2.0;
  const ParkLevel level = LevelAt(g, E);
  check("J by hand: -0.2 (1e-15)", std::fabs(level.J + 0.2) < 1e-15);
  check("value by hand: 40000 (1e-8)", std::fabs(ParkWallChi2({level}, kScale) - 40000.0) < 1e-8);
  check("value = pow(J/scale, 2.0), bit for bit", ParkWallChi2({level}, kScale) == std::pow(level.J / kScale, 2.0));

  // The gradient's form, by hand: factor 2 J/s^2 = -4e5; dgamma = -4e5 (-2 gamma S') = (672000, 192000);
  // dE = -4e5 (-(1.44 0.1 - 0.64 0.05)) = 44800.
  std::vector<double> dg;
  double dE;
  const bool on = ParkWallGradient(level, kScale, dg, dE);
  check("gradient by hand: (672000, 192000), dE 44800 (1e-9 relative)",
        on && dg.size() == 2 && Close(dg[0], 672000.0, 1e-9) && Close(dg[1], 192000.0, 1e-9) && Close(dE, 44800.0, 1e-9));
  bool central = true;
  for (int c = 0; c < 2; c++) {
    const double h = 1e-7;
    double up[2] = {g[0], g[1]}, down[2] = {g[0], g[1]};
    up[c] += h;
    down[c] -= h;
    if (!Close((Value(up, E) - Value(down, E)) / (2.0 * h), dg[c], 1e-6)) central = false;
  }
  const double hE = 1e-6;
  if (!Close((Value(g, E + hE) - Value(g, E - hE)) / (2.0 * hE), dE, 1e-6)) central = false;
  check("gradient = central differences of the value in gamma and E (1e-6 relative)", central);

  // The row's form: r = J/s = -200, dr = dJ/s against central differences of r.
  double r;
  std::vector<double> rg;
  double rE;
  const bool rowOn = ParkWallRow(level, kScale, r, rg, rE);
  check("row by hand: r = -200, dr/dgamma = (-1680, -480), dr/dE = -112 (1e-9 relative)",
        rowOn && Close(r, -200.0, 1e-12) && Close(rg[0], -1680.0, 1e-9) && Close(rg[1], -480.0, 1e-9) &&
            Close(rE, -112.0, 1e-9));
  bool rowCentral = true;
  for (int c = 0; c < 2; c++) {
    const double h = 1e-6;
    double up[2] = {g[0], g[1]}, down[2] = {g[0], g[1]};
    up[c] += h;
    down[c] -= h;
    if (!Close((Row(up, E) - Row(down, E)) / (2.0 * h), rg[c], 1e-7)) rowCentral = false;
  }
  if (!Close((Row(g, E + hE) - Row(g, E - hE)) / (2.0 * hE), rE, 1e-7)) rowCentral = false;
  check("row derivatives = central differences of r in gamma and E (1e-7 relative)", rowCentral);
  check("value = r^2 (1e-12 relative)", Close(r * r, ParkWallChi2({level}, kScale), 1e-12));
  // Gauss-Newton: 2 r dr of the row's form = the gradient's form, to rounding (the two forms round differently).
  check("2 r dr = the gradient (1e-12 relative)",
        Close(2.0 * r * rg[0], dg[0], 1e-12) && Close(2.0 * r * rg[1], dg[1], 1e-12) && Close(2.0 * r * rE, dE, 1e-12));

  // Off the wall: J = 1 - (0.25 0.7 + 0.25 0.3) = 0.75 >= 0, and a NaN J.
  const double small[2] = {0.5, 0.5};
  ParkLevel off = LevelAt(small, E);
  ParkLevel nan = level;
  nan.J = std::numeric_limits<double>::quiet_NaN();
  for (const ParkLevel *l : {&off, &nan}) {
    const char *name = l == &off ? "J = 0.75" : "J = NaN";
    std::vector<double> og;
    double oE = 1.0, orow = 1.0;
    const bool gOn = ParkWallGradient(*l, kScale, og, oE);
    check(std::string(name) + ": no value", ParkWallChi2({*l}, kScale) == 0.0);
    check(std::string(name) + ": no gradient (false, zeros)", !gOn && og.size() == 2 && og[0] == 0.0 && og[1] == 0.0 && oE == 0.0);
    const bool rOn = ParkWallRow(*l, kScale, orow, og, oE);
    check(std::string(name) + ": a zero row (false, zeros)", !rOn && orow == 0.0 && og[0] == 0.0 && og[1] == 0.0 && oE == 0.0);
  }
  check("J = 0 is not on the wall", !OnParkWall(0.0) && OnParkWall(-1e-300));
  // A list: only the level on the wall counts, in list order.
  check("a list of three levels: only the one on the wall counts",
        ParkWallChi2({off, level, nan}, kScale) == ParkWallChi2({level}, kScale));
  return result();
}
