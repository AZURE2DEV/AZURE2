// transfer_bound_state (docs/upgrade/transfer.md sections 3 and 9): the bound state of two clusters, each property
// against an answer that does not come from the library's solver.
//
//  1. Square well, l = 0, no Coulomb (mu = 1 u, R = 3 fm, B = 2 MeV), 0 and 1 interior nodes: the depth against the
//     root of K cot(K R) = -kappa solved here by bisection, to 1e-8 relative at h = 0.5e-3 fm (the step on the node
//     at R takes the mean of the two sides; the error is O(h^2), and the test checks the factor 4 between h and h/2).
//  2. Normalisation: int u^2 dr = 1 by the trapezoid rule with the Euler-Maclaurin end correction (not Simpson), 1e-9;
//     u > 0 next to the origin and the sign of the ANC (-1)^n.
//  3. Tail: u(r) = C W_{-eta, l+1/2}(2 kappa r) over the last 5 fm of the grid, to 1e-6 relative, for a neutral s state
//     (W = e^{-kappa r}, written out here) and for p + 15N in p1/2 with Coulomb and spin-orbit (W from the Whittaker
//     functions of physics/coulomb, which are checked against mpmath by physics_whittaker).
//  4. The analytic tails: Yukawa C e^{-kappa r} and Whittaker C W(2 kappa r) (= C e^{-kappa r} for eta = 0, l = 0),
//     zero below r_min.
//  5. Refusals: B <= 0, an imaginary well, a matching radius outside the grid.
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Whittaker.h"
#include "azure2/physics/transfer/BoundState.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

const double kPi = std::acos(-1.0);

double TwoMuOverHbar2(double mu_u) { return 2.0 * mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm); }

// ---------------------------------------------------------------- 1. square well
double SquareWellDepth(double mu_u, double R, double B, int nodes) {
  const double k2mu = TwoMuOverHbar2(mu_u), kappa = std::sqrt(k2mu * B);
  auto f = [&](double V) {
    const double K = std::sqrt(k2mu * (V - B));
    return K * std::cos(K * R) + kappa * std::sin(K * R);
  };
  // the n-node root has K R in ((n + 1/2) pi, (n + 1) pi)
  double lo = std::pow((nodes + 0.5) * kPi / R, 2) / k2mu + B, hi = std::pow((nodes + 1.0) * kPi / R, 2) / k2mu + B;
  lo *= 1.0 + 1e-12;
  hi *= 1.0 - 1e-12;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    (f(mid) * f(lo) > 0.0 ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

BoundState SquareWell(double R, double B, int nodes, double h) {
  BoundStateRequest q;
  q.pair.mu_u = 1.0;
  q.binding_MeV = B;
  q.nodes = nodes;
  q.step_fm = h;
  q.r_max_fm = 15.0;
  q.matching_radius_fm = R;
  auto form = [R](double r) { return std::fabs(r - R) < 1e-9 ? 0.5 : (r < R ? 1.0 : 0.0); };
  return SolveBoundStateOfForm(form, [](double) { return 0.0; }, q);
}

void TestSquareWell() {
  const double R = 3.0, B = 2.0;
  for (int nodes = 0; nodes <= 1; ++nodes) {
    const double exact = SquareWellDepth(1.0, R, B, nodes);
    const BoundState fine = SquareWell(R, B, nodes, 0.0005);
    const double err = std::fabs(fine.depth_MeV - exact) / exact;
    char buf[200];
    std::snprintf(buf, sizeof buf, "square well, %d nodes: depth %.12f against %.12f (rel %.2e)", nodes,
                  fine.depth_MeV, exact, err);
    Check(err < 1e-8, buf);
    const double e1 = SquareWell(R, B, nodes, 0.004).depth_MeV - exact;
    const double e2 = SquareWell(R, B, nodes, 0.002).depth_MeV - exact;
    std::snprintf(buf, sizeof buf, "square well, %d nodes: error ratio h/(h/2) = %.3f, not 4 (O(h^2))", nodes, e1 / e2);
    Check(std::fabs(e1 / e2 - 4.0) < 0.1, buf);

    // 2. normalisation by the trapezoid rule with the end correction h^2/12 [f'(0) - f'(r_max)] (f = u^2, f'(0) = 0
    //    for l = 0 since f ~ r^2; f'(r_max) negligible)
    double sum = 0.0;
    const std::vector<double> &u = fine.u;
    for (std::size_t i = 0; i < u.size(); ++i) sum += (i == 0 || i + 1 == u.size() ? 0.5 : 1.0) * u[i] * u[i];
    sum *= fine.step_fm;
    std::snprintf(buf, sizeof buf, "square well, %d nodes: int u^2 = %.12f", nodes, sum);
    Check(std::fabs(sum - 1.0) < 1e-9, buf);
    Check(u[1] > 0.0, "u > 0 next to the origin");
    Check((nodes % 2 == 0) == (fine.anc_fm_half > 0.0), "the sign of the ANC is (-1)^n");
    // 3a. neutral s-wave tail: u = C e^{-kappa r} beyond R exactly
    double worst = 0.0;
    const std::size_t n = u.size();
    for (std::size_t i = n - 1 - static_cast<std::size_t>(std::lround(5.0 / fine.step_fm)); i < n; ++i)
      worst = std::max(worst, std::fabs(u[i] / (fine.anc_fm_half * std::exp(-fine.kappa_fm * i * fine.step_fm)) - 1.0));
    std::snprintf(buf, sizeof buf, "square well, %d nodes: tail against C e^{-kappa r}, worst %.2e", nodes, worst);
    Check(worst < 1e-6, buf);
    Check(std::fabs(fine.kappa_fm - std::sqrt(TwoMuOverHbar2(1.0) * B)) < 1e-14, "kappa");
  }
}

// ---------------------------------------------------------------- 3b. charged tail (p + 15N, p1/2)
void TestChargedTail() {
  OpticalPotential w;
  w.convention = RadiusConvention::TargetMass;
  w.target_mass_number = 15;
  w.coulomb_r0_fm = 1.25;
  w.real_volume = {50.0, 1.25, 0.65};
  w.real_spin_orbit = {6.0, 1.25, 0.65};
  BoundStateRequest q;
  q.pair.z1 = 1;
  q.pair.z2 = 7;
  q.pair.mu_u = 15.0 / 16.0;
  q.l = 1;
  q.two_j = 1;
  q.two_s = 1;
  q.binding_MeV = 12.1274;
  q.step_fm = 0.02;
  q.r_max_fm = 25.0;
  const BoundState s = SolveBoundState(w, q);
  double worst = 0.0;
  const std::size_t n = s.u.size();
  for (std::size_t i = n - 1 - static_cast<std::size_t>(std::lround(5.0 / s.step_fm)); i < n; ++i) {
    const coulomb::WhittakerScaled ws = coulomb::whittaker_scaled(1, s.eta, s.kappa_fm * i * s.step_fm);
    const double W = ws.mantissa * std::pow(10.0, static_cast<double>(ws.log10_scale));
    worst = std::max(worst, std::fabs(s.u[i] / (s.anc_fm_half * W) - 1.0));
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "p + 15N p1/2: tail against C W(2 kappa r) over the last 5 fm, worst %.2e", worst);
  Check(worst < 1e-6, buf);
  const double eta = q.pair.z1 * q.pair.z2 * kFineStructure * std::sqrt(q.pair.mu_u * kAtomicMassUnit_MeV / (2.0 * q.binding_MeV));
  Check(std::fabs(s.eta - eta) < 1e-12, "eta_b = Z1 Z2 alpha sqrt(mu c^2 / 2B)");
}

// ---------------------------------------------------------------- 4. analytic tails
void TestAnalyticTails() {
  BoundStateRequest q;
  q.pair.mu_u = 0.5;
  q.binding_MeV = 2.2246;
  q.step_fm = 0.05;
  q.r_max_fm = 20.0;
  const double kappa = std::sqrt(TwoMuOverHbar2(0.5) * q.binding_MeV);
  const BoundState y = AnalyticTail(TailShape::Yukawa, q, 0.88, 1.0);
  const BoundState w = AnalyticTail(TailShape::Whittaker, q, 0.88, 1.0);
  double worst = 0.0;
  for (std::size_t i = 0; i < y.u.size(); ++i) {
    const double r = i * q.step_fm;
    const double expect = r < 1.0 ? 0.0 : 0.88 * std::exp(-kappa * r);
    worst = std::max({worst, std::fabs(y.u[i] - expect), std::fabs(w.u[i] - expect)});
  }
  Check(worst < 1e-13, "analytic tails: Yukawa and Whittaker (eta = 0, l = 0) against 0.88 e^{-kappa r} beyond 1 fm");
}

// ---------------------------------------------------------------- 5. refusals
void TestRefusals() {
  OpticalPotential w;
  w.real_volume = {50.0, 3.0, 0.65};
  BoundStateRequest q;
  q.pair.mu_u = 1.0;
  q.binding_MeV = 0.0;
  bool threw = false;
  try {
    SolveBoundState(w, q);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "B = 0 refused");
  q.binding_MeV = 1.0;
  w.imaginary_volume = {5.0, 3.0, 0.65};
  threw = false;
  try {
    SolveBoundState(w, q);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "an imaginary well refused");
  threw = false;
  q.matching_radius_fm = 100.0;
  try {
    SolveBoundStateOfForm([](double) { return 1.0; }, [](double) { return 0.0; }, q);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "a matching radius outside the grid refused");
}

}  // namespace

int main() {
  TestSquareWell();
  TestChargedTail();
  TestAnalyticTails();
  TestRefusals();
  if (failures) {
    std::printf("transfer_bound_state: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_bound_state: all checks passed\n");
  return 0;
}
