// The adaptive sub-point grid (AdaptiveIntegrationGrid) against closed forms.
//
// Every experimental effect that integrates over energy interpolates the cross
// section linearly between sub-points (EPoint.cpp GridQuadrature), so the
// integral of a resonance on a grid is the trapezoid rule on that grid.  Here the
// "cross section" is a Lorentzian of width G centred on E_R,
//     f(E) = (G/2)^2 / ((E - E_R)^2 + (G/2)^2),
// whose integral over [a, b] is (G/2) [atan((b - E_R)/(G/2)) - atan((a - E_R)/(G/2))].
// The resonance is moved across a fixed window in small steps, as a fit moves a
// level, and the grid is laid out again at each position, as the engine does
// (EData::ReplaceSubPointGrids).  Two things are checked, relative to the
// resonance's area pi G / 2:
//   - accuracy: the trapezoid integral equals the exact one;
//   - smoothness: the error does not swing as the resonance moves.  A grid whose
//     nodes near the resonance do not move with it integrates the tail
//     differently at every position, and the objective of a fit oscillates in the
//     level energy (17O+a 2026-10-22: a 0.8 eV level under a 50-point target
//     integration, chi2 +-17 with a 10 eV period, and an analytic gradient of the
//     wrong sign compared with finite differences).
//
// Usage: adaptive_grid_reference_test [-v]   (-v prints every case's numbers)

#include "AdaptiveIntegrationGrid.h"
#include "Config.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

int failures = 0;
bool verbose = false;

void check(const char *name, bool ok, const char *fmt = "", double a = 0.0, double b = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", name);
  if (!ok || verbose) {
    std::printf("  -- ");
    std::printf(fmt, a, b);
  }
  std::printf("\n");
  if (!ok) failures++;
}

double Exact(double a, double b, double eR, double g) {
  const double h = 0.5 * g;
  return h * (std::atan((b - eR) / h) - std::atan((a - eR) / h));
}

double Lorentz(double e, double eR, double g) {
  const double h = 0.5 * g;
  return h * h / ((e - eR) * (e - eR) + h * h);
}

// Trapezoid on a grid in descending order.
double Trapezoid(const std::vector<double> &grid, double eR, double g) {
  double s = 0.0;
  for (size_t k = 0; k + 1 < grid.size(); k++)
    s += 0.5 * (grid[k] - grid[k + 1]) * (Lorentz(grid[k], eR, g) + Lorentz(grid[k + 1], eR, g));
  return s;
}

struct Result {
  double maxErr = 0.0;     // max |trapezoid - exact| / area over the scan
  double swing = 0.0;      // max - min of (trapezoid - exact) / area over the scan
  double maxStepErr = 0.0; // max |change of the error between neighbouring positions| / area
  size_t maxNodes = 0;
};

// Moves a resonance of width g (MeV) across the window [lo, hi] (MeV) in `steps` steps
// over [from, to], laying out the grid at each position.
Result Scan(double lo, double hi, int numPoints, double g, double mult, double ppw, double from, double to,
            int steps) {
  AdaptiveIntegrationGrid::GridConfig gc;
  gc.maxPoints = numPoints;
  gc.baseEnergyStep = (hi - lo) / numPoints;  // as EData::ReadTargetEffectsFile sets it
  gc.resonanceWidthMultiplier = mult;
  gc.pointsPerWidth = ppw;
  AdaptiveIntegrationGrid grid(gc);
  const double area = M_PI * 0.5 * g;
  Result r;
  double emin = 1e300, emax = -1e300, prev = 0.0;
  for (int i = 0; i <= steps; i++) {
    const double eR = from + (to - from) * i / steps;
    std::vector<AdaptiveIntegrationGrid::ResonanceInfo> res(1);
    res[0].energy = eR;
    res[0].particleWidth = g;
    res[0].totalWidth = g;
    std::vector<double> nodes = grid.GenerateGrid(hi, lo, res);
    const double err = (Trapezoid(nodes, eR, g) - Exact(lo, hi, eR, g)) / area;
    r.maxErr = std::max(r.maxErr, std::fabs(err));
    emin = std::min(emin, err);
    emax = std::max(emax, err);
    if (i > 0) r.maxStepErr = std::max(r.maxStepErr, std::fabs(err - prev));
    prev = err;
    r.maxNodes = std::max(r.maxNodes, nodes.size());
  }
  r.swing = emax - emin;
  return r;
}

void Case(const char *name, double window, int numPoints, double g, double mult, double ppw, double tolErr,
          double tolSwing, double from = 0.25, double to = 0.75) {
  // the window is [10, 10 + window] MeV; the resonance moves over [from, to] of it
  const double lo = 10.0, hi = 10.0 + window;
  Result r = Scan(lo, hi, numPoints, g, mult, ppw, lo + from * window, lo + to * window, 2000);
  std::printf("%s: G = %.3g eV, window %.3g keV, %d points, mult %g, ppw %g -> max error %.2e, swing %.2e, "
              "largest jump %.2e, up to %zu nodes\n",
              name, g * 1e6, window * 1e3, numPoints, mult, ppw, r.maxErr, r.swing, r.maxStepErr, r.maxNodes);
  char label[256];
  std::snprintf(label, sizeof label, "%s: integral within %.0e of the exact area", name, tolErr);
  check(label, r.maxErr < tolErr, "max error %.3e of the area (tolerance %.0e)", r.maxErr, tolErr);
  std::snprintf(label, sizeof label, "%s: error swings by less than %.0e as the level moves", name, tolSwing);
  check(label, r.swing < tolSwing, "swing %.3e of the area (tolerance %.0e)", r.swing, tolSwing);
}

}  // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; i++)
    if (!std::strcmp(argv[i], "-v")) verbose = true;

  std::printf("Lorentzian moved across a fixed window, grid laid out again at each position\n");
  // The 17O+a case: a 0.8 eV level in a target-integration window, default settings.  Before the
  // graded tail (2026-10-22) this one was 7.9e-2 off and swung by 3.5e-2.
  Case("eV level, 10 keV window, 50 points", 10e-3, 50, 0.8e-6, 20, 50, 2e-4, 1e-4);
  Case("eV level, 500 eV window, 50 points", 0.5e-3, 50, 0.8e-6, 20, 50, 2e-4, 1e-4);
  Case("eV level, 10 keV window, 2000 points", 10e-3, 2000, 0.8e-6, 20, 50, 2e-4, 1e-4);
  // A keV-wide level (thick_target_resonance's 14N(p,g) 259.7 keV) in a 100 keV window.
  Case("keV level, 100 keV window, 20 points", 0.1, 20, 0.97e-3, 20, 50, 2e-4, 1e-4);
  // A level as wide as the coarse step: anchored, its lattice reaches past the graded zone.
  Case("level as wide as the step, 50 points", 10e-3, 50, 0.2e-3, 20, 50, 2e-4, 1e-4);
  // The level leaves the window: its nodes inside must not jump when it crosses the edge
  // (IdentifyResonances keeps a resonance whose zone reaches into the window).  The integral
  // is the partial one; the swing is of the error, so it stays small too.
  Case("eV level crossing the window edge", 10e-3, 50, 0.8e-6, 20, 50, 2e-4, 1e-4, 0.9, 1.1);
  Case("keV level crossing the window edge", 0.1, 20, 0.97e-3, 20, 50, 2e-4, 1e-4, 0.9, 1.1);

  std::printf("\n");
  if (failures) {
    std::printf("FAILED: %d check(s)\n", failures);
    return 1;
  }
  std::printf("all checks passed\n");
  return 0;
}
