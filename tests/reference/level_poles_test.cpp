// NarrowPolesOf (LevelPoles.h): the poles of interfering levels of one J^pi, for the energy mesh.  The mesh took a
// level's width as the sum of its partial widths, which is the width of its resonance only when the level is
// isolated: 22Ne(a,n), 1- levels at E_x 11.8847 and 11.8913 MeV with 3.06 and 93.95 keV in the n0 channel, has a
// peak 33 eV wide where the mesh assumed 3.08 keV, and a folded point was 25 % off (2026-11-25).
#include "LevelPoles.h"

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}

int main() {
  using azure2::NarrowPole;
  using azure2::PoleLevel;
  auto level = [](double energy, std::vector<double> widths, std::vector<double> signs) {
    PoleLevel l;
    l.energy = energy;
    for (std::size_t c = 0; c < widths.size(); c++) {
      l.width += widths[c];
      l.amplitude.push_back(signs[c] * std::sqrt(widths[c]));
    }
    return l;
  };

  // 1. Two levels in one channel, strongly overlapping: 4 dE^2 G1 G2 / (G1 + G2)^3 to leading order.
  {
    const double g1 = 3.059e-3, g2 = 93.952e-3, dE = 6.666e-3;
    std::vector<NarrowPole> p = azure2::NarrowPolesOf({level(1.270, {g1}, {1.0}), level(1.270 + dE, {g2}, {-1.0})});
    const double estimate = 4.0 * dE * dE * g1 * g2 / std::pow(g1 + g2, 3);
    check(p.size() == 1, "two overlapping levels in one channel: one narrow pole");
    if (p.size() == 1) {
      std::printf("    pole at %.6f MeV, width %.2f eV (leading-order estimate %.2f eV; the level's own %.0f eV)\n", p[0].energy,
                  p[0].width * 1e6, estimate * 1e6, g1 * 1e6);
      check(std::abs(p[0].width / estimate - 1.0) < 0.1, "its width is the strong-overlap estimate within 10 %");
      check(p[0].width < 0.03 * g1, "and far below the narrow level's partial width");
      check(std::abs(p[0].energy - 1.270) < 1e-3, "it lies at the narrow level, within 1 keV");
    }
    // the sign of the amplitudes does not matter for one channel
    std::vector<NarrowPole> q = azure2::NarrowPolesOf({level(1.270, {g1}, {1.0}), level(1.270 + dE, {g2}, {1.0})});
    check(q.size() == 1 && p.size() == 1 && std::abs(q[0].width - p[0].width) < 1e-12, "the same with equal signs");
  }
  // 2. The 22Ne(a,n) pair with its five channels (alpha, n0, gamma left out, n1 s = 0 and 1), eV.
  {
    PoleLevel a = level(1.270, {0.100867e-6, 3058.973e-6, 0.093026e-6, 22.698344e-6}, {1, 1, -1, -1});
    PoleLevel b = level(1.270 + 6.666e-3, {0.267304e-6, 93952.363e-6, 682.390584e-6, 2794.349e-6}, {1, -1, 1, 1});
    std::vector<NarrowPole> p = azure2::NarrowPolesOf({a, b});
    check(p.size() == 1, "the 22Ne(a,n) 1- pair: one narrow pole");
    if (p.size() == 1) {
      std::printf("    pole at %.6f MeV, width %.1f eV\n", p[0].energy, p[0].width * 1e6);
      check(p[0].width > 10e-6 && p[0].width < 150e-6, "tens of eV wide, where the level's partial widths sum to 3.08 keV");
    }
  }
  // 3. Isolated levels, levels without a common channel, and a single level: nothing.
  check(azure2::NarrowPolesOf({level(1.0, {1e-3}, {1.0}), level(1.5, {2e-3}, {1.0})}).empty(), "two levels far apart: none");
  check(azure2::NarrowPolesOf({level(1.0, {1e-3, 0.0}, {1.0, 1.0}), level(1.0005, {0.0, 2e-3}, {1.0, 1.0})}).empty(),
        "two close levels in different channels: none");
  check(azure2::NarrowPolesOf({level(1.0, {1e-3}, {1.0})}).empty(), "one level: none");
  check(azure2::NarrowPolesOf({}).empty(), "no level: none");
  // 4. Mild overlap (dE = 2 (G1 + G2) is outside the screen; dE = G1 + G2 changes the widths by less than 20 %): none.
  check(azure2::NarrowPolesOf({level(1.0, {1e-3}, {1.0}), level(1.0 + 4e-3, {1e-3}, {1.0})}).empty(), "levels two summed widths apart: none");
  {
    std::vector<NarrowPole> p = azure2::NarrowPolesOf({level(1.0, {1e-3}, {1.0}), level(1.0 + 1.9e-3, {1e-3}, {1.0})});
    check(p.empty(), "equal levels just inside the screen: the poles are not 20 % narrower, none");
  }
  // 5. Three levels: a narrow one between two broad ones in the same channel, and a far one that plays no part.
  {
    std::vector<NarrowPole> p = azure2::NarrowPolesOf({level(1.000, {50e-3}, {1.0}), level(1.004, {2e-3}, {1.0}), level(3.0, {1e-3}, {1.0})});
    check(p.size() == 1 && p[0].width < 1.6e-3 && std::abs(p[0].energy - 1.004) < 2e-3, "three levels: the pole of the narrow one in the broad one");
  }
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
