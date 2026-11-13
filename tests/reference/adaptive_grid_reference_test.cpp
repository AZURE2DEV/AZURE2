// The adaptive sub-point grid (AdaptiveIntegrationGrid) against closed forms.
//
// Every experimental effect that integrates over energy interpolates the cross
// section linearly between sub-points (Unit.cpp GridQuadrature), so the
// integral of a resonance on a grid is the trapezoid rule on that grid.  Here the
// "cross section" is a Lorentzian of width G centred on E_R,
//     f(E) = (G/2)^2 / ((E - E_R)^2 + (G/2)^2),
// whose integral over [a, b] is (G/2) [atan((b - E_R)/(G/2)) - atan((a - E_R)/(G/2))].
// The resonance is moved across a fixed window in small steps, as a fit moves a
// level, and the grid is laid out again at each position, as the engine does
// (WorkingCopy::ReplaceSubPointGrids).  Two things are checked, relative to the
// resonance's area pi G / 2:
//   - accuracy: the trapezoid integral equals the exact one;
//   - smoothness: the error does not swing as the resonance moves.  A grid whose
//     nodes near the resonance do not move with it integrates the tail
//     differently at every position, and the objective of a fit oscillates in the
//     level energy (17O+a 2026-10-22: a 0.8 eV level under a 50-point target
//     integration, chi2 +-17 with a 10 eV period, and an analytic gradient of the
//     wrong sign compared with finite differences).
// For a level that crosses the window edge the error is not expected to stay
// constant: its tail-bias part leaves the window with the level (a one-way drift,
// bounded by the accuracy check), and a lattice cut by a fixed edge has an error
// term that follows the slope of the peak at the edge (CrossingCase).  There the
// smoothness check is on the variation of the error beyond its net change,
// i.e. on the back-and-forth part only (decision 8 (a), docs/upgrade/HANDOFF.md).
//
// Usage: adaptive_grid_reference_test [-v]   (-v prints every case's numbers)

#include "AdaptiveIntegrationGrid.h"
#include "Config.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

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
  double net = 0.0;        // error at the last position minus the error at the first / area
  double excess = 0.0;     // total variation of the error minus |net|: what is not a one-way drift / area
  size_t maxNodes = 0;
};

// Moves a resonance of width g (MeV) across the window [lo, hi] (MeV) in `steps` steps
// over [from, to], laying out the grid at each position.
Result Scan(double lo, double hi, int numPoints, double g, double mult, double ppw, double from, double to,
            int steps) {
  AdaptiveIntegrationGrid::GridConfig gc;
  gc.maxPoints = numPoints;
  gc.baseEnergyStep = (hi - lo) / numPoints;  // as WorkingCopy::SetUpTargetEffects sets it
  gc.resonanceWidthMultiplier = mult;
  gc.pointsPerWidth = ppw;
  AdaptiveIntegrationGrid grid(gc);
  const double area = M_PI * 0.5 * g;
  Result r;
  double emin = 1e300, emax = -1e300, prev = 0.0, first = 0.0, variation = 0.0;
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
    if (i > 0) variation += std::fabs(err - prev);
    if (i == 0) first = err;
    prev = err;
    r.maxNodes = std::max(r.maxNodes, nodes.size());
  }
  r.swing = emax - emin;
  r.net = prev - first;
  r.excess = variation - std::fabs(r.net);
  return r;
}

void Case(const char *name, double window, int numPoints, double g, double mult, double ppw, double tolErr,
          double tolSwing, double from = 0.25, double to = 0.75) {
  // the window is [10, 10 + window] MeV; the resonance moves over [from, to] of it
  const double lo = 10.0, hi = 10.0 + window;
  Result r = Scan(lo, hi, numPoints, g, mult, ppw, lo + from * window, lo + to * window, 2000);
  std::printf("%s: G = %.3g eV, window %.3g keV, %d points, mult %g, ppw %g -> max error %.2e, swing %.2e, "
              "net change %+.2e, variation beyond it %.2e, largest jump %.2e, up to %zu nodes\n",
              name, g * 1e6, window * 1e3, numPoints, mult, ppw, r.maxErr, r.swing, r.net, r.excess, r.maxStepErr,
              r.maxNodes);
  char label[256];
  std::snprintf(label, sizeof label, "%s: integral within %.0e of the exact area", name, tolErr);
  check(label, r.maxErr < tolErr, "max error %.3e of the area (tolerance %.0e)", r.maxErr, tolErr);
  std::snprintf(label, sizeof label, "%s: error swings by less than %.0e as the level moves", name, tolSwing);
  check(label, r.swing < tolSwing, "swing %.3e of the area (tolerance %.0e)", r.swing, tolSwing);
}

// Error term of a lattice of pitch h = G/ppw cut by a fixed window edge, as the
// peak sweeps across that edge, relative to the area pi G / 2.
//
// With the edge at distance D from E_R, the cell next to the edge is clipped to
// length s < h.  Moving the edge by dD changes the trapezoid-minus-exact error by
// [f(D - s) - f(D) + s f'(D)] dD / 2 = s^2 f''(D) dD / 4 + O(s^3), so over one
// cell it changes by h^3 f''/12, and over many cells by (h^2/12) [f'(D2) - f'(D1)]:
// the Euler-Maclaurin end term, with the edge as the end.  It is monotone where
// f'' keeps its sign and reverses where it does not, which for a Lorentzian is
// at |E - E_R| = G/(2 sqrt 3), where |f'| = 3 sqrt 3 / (4 G) is largest.  As the
// level crosses, f'(edge) goes 0 -> -max -> +max -> 0, a total variation of
// 4 max|f'| with no net change, so the error moves back and forth by
//     (h^2/12) * 4 * 3 sqrt 3 / (4 G) / (pi G / 2) = sqrt 3 / (2 pi ppw^2),
// 1.10e-4 at 50 points per width.  The same term exists for any trapezoid
// lattice with a node at a fixed edge: it is the O(h^2) accuracy of the core,
// spread over about one width of level motion, not an oscillation.  Measured
// with a scan of 0.00125 G steps across the edge (eV level, E_R from 0.999 to
// 1.001 of the window): 1.102e-4.  Its slope in E_R is at most
// (h^2/12) max|f''| / area = 4/(3 pi ppw^2 G) = 1.7e-4/G, against 2/(pi G) for
// the in-window area itself at the edge: a 3e-4 change of the slope, which
// cannot reverse the gradient of a fit as the 17O+a ripple did.
double EdgeTermVariation(double ppw) { return std::sqrt(3.0) / (2.0 * M_PI * ppw * ppw); }

// A level moved across the window edge.  The accuracy check is the one of Case.
// The smoothness check is on r.excess, the total variation of the error over
// the scan minus |net change|.  A one-way drift costs nothing: here the tail
// bias leaving the window with the level (the graded tails over-count the area
// by about r^2 G/(2 pi core) of it, 2e-4 for 20-width cores, measured 1.69e-4;
// AdaptiveIntegrationGrid::TailFactor), which the accuracy check bounds.  Each
// back-and-forth excursion of size d costs 2d.  A sampled scan never shows more
// excess than the continuous curve, so the tolerance holds for any scan step.
// Tolerance: the edge term above, which every lattice of this pitch has, plus
// tolOther for everything else, the same 1e-4 of the area that the in-window
// cases allow the error to swing (here as a total, so stricter).
void CrossingCase(const char *name, double window, int numPoints, double g, double mult, double ppw,
                  double tolErr, double tolOther, double from, double to) {
  const double lo = 10.0, hi = 10.0 + window;
  Result r = Scan(lo, hi, numPoints, g, mult, ppw, lo + from * window, lo + to * window, 2000);
  const double edge = EdgeTermVariation(ppw);
  const double tolExcess = edge + tolOther;
  std::printf("%s: G = %.3g eV, window %.3g keV, %d points, mult %g, ppw %g -> max error %.2e, swing %.2e, "
              "net change %+.2e, variation beyond it %.2e (edge term %.2e), largest jump %.2e, up to %zu nodes\n",
              name, g * 1e6, window * 1e3, numPoints, mult, ppw, r.maxErr, r.swing, r.net, r.excess, edge,
              r.maxStepErr, r.maxNodes);
  char label[256];
  std::snprintf(label, sizeof label, "%s: integral within %.0e of the exact area", name, tolErr);
  check(label, r.maxErr < tolErr, "max error %.3e of the area (tolerance %.0e)", r.maxErr, tolErr);
  std::snprintf(label, sizeof label, "%s: error varies beyond its net change by less than %.2e", name, tolExcess);
  check(label, r.excess < tolExcess, "variation beyond the net change %.3e of the area (tolerance %.3e)", r.excess,
        tolExcess);
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
  // (AnchorsInReach keeps a resonance whose zone reaches into the window).  The integral is
  // the partial one; the check is on the variation of its error beyond the net change
  // (CrossingCase), tolerance sqrt3/(2 pi ppw^2) + 1e-4 = 2.10e-4 for both cases.
  //  - eV level: the net change is the tail bias leaving with the level (-1.69e-4).  The
  //    variation beyond it is 2.0e-5 here because the 1 eV scan step (1.25 G) samples the
  //    edge term coarsely; with 100 times finer steps it is 1.14e-4, the edge term (1.10e-4)
  //    plus 4e-6.  Before the graded tails (dev 07581f03) it was 2.87e-1 (swing 7.9e-2).
  //  - keV level: E_R stays within 10.3 G of the edge, inside the 20-width core, so the
  //    tail bias does not leave and the net change is 4e-7.  The variation beyond it is
  //    1.79e-4 (1.83e-4 with 100 times finer steps).  Split by cell: the cells from 20 G
  //    below E_R upward (core and edge) give 1.102e-4, the edge term; the cells below the
  //    core (the graded tail and its junction with the fixed 5 keV coarse lattice) give a
  //    ripple of 6.8e-6 peak to peak with many reversals, the rest.  The in-window keV case
  //    above has the same ripple (swing 2.4e-5, within its 1e-4).
  CrossingCase("eV level crossing the window edge", 10e-3, 50, 0.8e-6, 20, 50, 2e-4, 1e-4, 0.9, 1.1);
  CrossingCase("keV level crossing the window edge", 0.1, 20, 0.97e-3, 20, 50, 2e-4, 1e-4, 0.9, 1.1);

  std::printf("\n");
  if (failures) {
    std::printf("FAILED: %d check(s)\n", failures);
    return 1;
  }
  std::printf("all checks passed\n");
  return 0;
}
