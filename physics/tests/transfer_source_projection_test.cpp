// transfer_source_projection (docs/upgrade/transfer.md section 25): the projection S_lm(a) of the source of the
// distorted-wave THM vertex, against answers that do not come from it.
//
//  1. The numerics it is built of: SixPointLagrange reproduces a polynomial of degree 5 (1e-12) and its derivative
//     (1e-11);
//     SphericalHarmonicsWithDerivative equals SphericalHarmonicsAtPhiZero at theta >= 0 (1e-14), Y(-theta, 0) =
//     (-1)^m Y(theta, 0), and its derivative a central difference (1e-8).
//  2. The closed form: plane waves in both channels and a Yukawa source phi = e^{-kappa u}/u, phi~(q) =
//     4 pi/(q^2 + kappa^2), S(r) = phi~(q) e^{i p.r}, p = k_i - alpha k_f, q = k_f - beta k_i, so
//     S_lm(a) = 4 pi i^l j_l(pa) Y*_lm(p^) phi~(q) and S_lm'(a) = 4 pi i^l p j_l'(pa) Y*_lm(p^) phi~(q) (j_l from the
//     Coulomb functions at eta = 0), l <= 3, three directions: 1e-8 of the largest |S_lm|; SourceMomentumAmplitude
//     against 4 pi/(q^2 + kappa^2): 1e-9.
//  3. EntranceVertex built from the projection (w = r S_lm on the level's grid, w'(a) = S_lm(a) + a S_lm'(a)) against
//     PlaneWaveEntranceVertex times 4 pi phi~(q) i^l Y*_lm(p^) (transfer.md (25.3)): volume, surface and interior
//     terms, 1e-8 of the largest.
//  4. The engine's numbers before the projection moved here (its own quadrature; printed by a scratch driver from the
//     base build, tests/reference/thm_dw_vertex_test.cpp's point-Coulomb 19F(d,n) cases at 55 MeV): the Gram entries
//     of (S_lm, a S_lm')/4 pi phi~ summed over m, l = 0, 1, 3, to 1e-9 of the largest entry.
//  5. The order of the sum: the blocks summed in reverse order give ProjectSource's numbers bit for bit.
//  6. Refinement of case 4: twice the nodes in u and in cos(theta_u), partial waves + 20 and the wave tables' steps
//     halved change the Gram entries by less than 1e-7 of the largest (measured 1.2e-8 to 3.7e-8).
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/BoundState.h"
#include "azure2/physics/transfer/RadialNumerics.h"
#include "azure2/physics/transfer/SourceProjection.h"
#include "azure2/physics/transfer/ThmVertex.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);
int failures = 0;
bool verbose = false;

void Check(const std::string &what, double distance, double tolerance) {
  const bool ok = distance <= tolerance;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-60s %.3e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), distance, tolerance);
}

cplx IPower(int l) {
  const cplx p[4] = {cplx(1, 0), cplx(0, 1), cplx(-1, 0), cplx(0, -1)};
  return p[((l % 4) + 4) % 4];
}

// j_l(x) and j_l'(x) from the Coulomb functions at eta = 0: F = x j_l, F' = d(x j_l)/dx.
void Bessel(int l, double x, double &j, double &dj) {
  const coulomb::CoulombWaves c = CoulombForMatching(l, 0.0, x);
  j = c.F / x;
  dj = (c.Fp - j) / x;
}

void TestNumerics() {
  // A polynomial of degree 5 and its derivative.
  auto f = [](double r) {
    return 0.3 - 1.1 * r + 0.7 * r * r - 0.05 * r * r * r + 0.02 * std::pow(r, 4) - 0.001 * std::pow(r, 5);
  };
  auto df = [](double r) { return -1.1 + 1.4 * r - 0.15 * r * r + 0.08 * std::pow(r, 3) - 0.005 * std::pow(r, 4); };
  const double h = 0.1;
  const int n = 60;
  std::vector<double> g(n);
  for (int i = 0; i < n; i++) g[i] = f(i * h);
  double worst = 0.0, worst_d = 0.0;
  for (double r : {0.0, 0.013, 0.25, 1.77, 3.333, 5.88, 5.9}) {
    const SixPointStencil st = SixPointLagrange(r, h, n);
    double v = 0.0, dv = 0.0;
    for (int i = 0; i < 6; i++) {
      v += st.value[i] * g[st.first + i];
      dv += st.derivative[i] * g[st.first + i];
    }
    worst = std::max(worst, std::fabs(v - f(r)));
    worst_d = std::max(worst_d, std::fabs(dv - df(r)));
  }
  Check("six-point stencil: value of a quintic", worst, 1e-12);
  Check("six-point stencil: derivative of a quintic", worst_d, 1e-11);

  const int l_max = 30, m_max = 4;
  double same = 0.0, mirror = 0.0, deriv = 0.0;
  for (double th : {0.0, 0.3, 1.2, 2.0, 3.0, kPi}) {
    std::vector<double> y, dy, ym, dym, yp, dyp, yq, dyq;
    SphericalHarmonicsWithDerivative(l_max, m_max, th, y, dy);
    const std::vector<double> x = SphericalHarmonicsAtPhiZero(l_max, m_max, std::cos(th));
    SphericalHarmonicsWithDerivative(l_max, m_max, -th, ym, dym);
    const double e = 1e-6;
    SphericalHarmonicsWithDerivative(l_max, m_max, th + e, yp, dyp);
    SphericalHarmonicsWithDerivative(l_max, m_max, th - e, yq, dyq);
    for (int l = 0; l <= l_max; l++)
      for (int m = 0; m <= std::min(l, m_max); m++) {
        const int i = l * (m_max + 1) + m;
        same = std::max(same, std::fabs(y[i] - x[i]));
        mirror = std::max(mirror, std::fabs(ym[i] - (m % 2 ? -1.0 : 1.0) * y[i]));
        deriv = std::max(deriv, std::fabs(dy[i] - (yp[i] - yq[i]) / (2 * e)) / std::max(1.0, std::fabs(dy[i])));
      }
  }
  Check("Y_lm(theta, 0): the same as SphericalHarmonicsAtPhiZero", same, 1e-12);
  Check("Y_lm(-theta, 0) = (-1)^m Y_lm(theta, 0)", mirror, 1e-14);
  Check("dY_lm/dtheta against a central difference", deriv, 1e-7);
}

// d + 19F -> n + 20Ne-like masses (u) for the plane-wave cases.
const double kMd = 2.0135532134, kMF = 18.9934652, kMn = 1.0086649159, kMp = 1.0072764675;

WaveChannel Plane(double mu_u, double energy) {
  WaveChannel c;
  c.plane_wave = true;
  c.pair.mu_u = mu_u;
  c.energy_MeV = energy;
  return c;
}

struct PlaneCase {
  SourceSetup setup;
  WaveChannel exit;
  double k_i = 0.0, k_f = 0.0;
};

PlaneCase MakePlaneCase(double radius, double kappa, const std::vector<int> &ls) {
  PlaneCase c;
  const double mu_i = kMd * kMF / (kMd + kMF), mu_f = kMn * (kMF + kMp) / (kMn + kMF + kMp);
  c.setup.entrance = Plane(mu_i, 7.0);
  c.exit = Plane(mu_f, 4.0);
  c.k_i = rho(c.setup.entrance.pair, 1.0, 7.0);
  c.k_f = rho(c.exit.pair, 1.0, 4.0);
  c.setup.overlap.phi = [kappa](double u) { return u > 0.0 ? std::exp(-kappa * u) / u : 0.0; };
  c.setup.overlap.kappa_fm = kappa;
  c.setup.entrance_step_fm = c.setup.exit_step_fm = 0.02;
  c.setup.exit_k_max_fm = c.k_f;
  c.setup.alpha = kMF / (kMF + kMp);
  c.setup.beta = kMn / kMd;
  c.setup.radius_fm = radius;
  c.setup.l_values = ls;
  return c;
}

// p^ (polar angle, azimuth pi or 0) and |p|, |q| for k_i along z, k_f at cos_fi in the xz plane (x > 0).
void PlaneKinematics(const PlaneCase &c, double cos_fi, double &p, double &theta_p, int &azimuth_sign, double &q) {
  const double s = std::sqrt(1.0 - cos_fi * cos_fi);
  const double px = -c.setup.alpha * c.k_f * s, pz = c.k_i - c.setup.alpha * c.k_f * cos_fi;
  p = std::hypot(px, pz);
  theta_p = std::atan2(std::fabs(px), pz);
  azimuth_sign = px < 0.0 ? -1 : 1;  // Y_lm(theta, pi) = (-1)^m Y_lm(theta, 0)
  q = std::hypot(c.k_f * s, c.k_f * cos_fi - c.setup.beta * c.k_i);
}

void TestClosedForm() {
  const double kappa = 0.6, a = 5.0;
  const std::vector<int> ls = {0, 1, 2, 3};
  PlaneCase c = MakePlaneCase(a, kappa, ls);
  const SourceProjectionPlan plan = PrepareSourceProjection(c.setup, SourceProjectionRequest{});
  const SourceAmplitudes amp = ProjectSource(plan, {c.exit})[0];
  double worst = 0.0, worst_phi = 0.0;
  for (double cos_fi : {1.0, 0.8, -0.3}) {
    double p, theta_p, q;
    int sign;
    PlaneKinematics(c, cos_fi, p, theta_p, sign, q);
    const double phit = 4.0 * kPi / (q * q + kappa * kappa);
    worst_phi = std::max(worst_phi, std::fabs(SourceMomentumAmplitude(plan, q) / phit - 1.0));
    const SourceMultipoles s = SourceMultipolesAt(plan, amp, cos_fi);
    const std::vector<double> y = SphericalHarmonicsAtPhiZero(3, 3, std::cos(theta_p));
    double largest = 0.0, distance = 0.0;
    for (int li = 0; li < (int)ls.size(); li++) {
      const int l = ls[li];
      double j, dj;
      Bessel(l, p * a, j, dj);
      for (int m = 0; m <= l; m++) {
        const double ylm = y[l * 4 + m] * (sign < 0 && m % 2 ? -1.0 : 1.0);
        const cplx want = 4.0 * kPi * IPower(l) * j * ylm * phit;
        const cplx want_d = 4.0 * kPi * IPower(l) * p * dj * ylm * phit;
        largest = std::max({largest, std::abs(want), a * std::abs(want_d)});
        distance = std::max({distance, std::abs(s.value[li][m] - want), a * std::abs(s.derivative[li][m] - want_d)});
      }
    }
    worst = std::max(worst, distance / largest);
  }
  Check("plane waves, Yukawa source: S_lm(a), a S_lm'(a) against the closed form", worst, 1e-8);
  Check("phi~(q) on the plan's nodes against 4 pi/(q^2 + kappa^2)", worst_phi, 1e-9);
}

void TestEntranceVertex() {
  // A level function on [0, a] (u = sin(K r), V = -30 MeV), and the wave w = r S_lm(r) at the grid's radii.
  const double kappa = 0.9, a = 4.0;
  const int nodes = 21, l = 2, m = 1;
  const double cos_fi = 0.7;
  LevelFunction level;
  level.step_fm = a / (nodes - 1);
  level.mu_MeV = kMp * kMF / (kMp + kMF) * kAtomicMassUnit_MeV;
  level.energy_MeV = 1.3;
  for (int i = 0; i < nodes; i++) {
    level.u.push_back(std::sin(0.83 * i * level.step_fm));
    level.potential_MeV.push_back(-30.0);
  }
  const double B = BoundaryValue(level);
  std::vector<cplx> w(nodes, 0.0);
  cplx dw;
  double p = 0.0, theta_p = 0.0, q = 0.0;
  int sign = 1;
  double phit = 0.0;
  for (int i = 1; i < nodes; i++) {
    PlaneCase c = MakePlaneCase(i * level.step_fm, kappa, {l});
    const SourceProjectionPlan plan = PrepareSourceProjection(c.setup, SourceProjectionRequest{});
    const SourceMultipoles s = SourceMultipolesAt(plan, ProjectSource(plan, {c.exit})[0], cos_fi);
    w[i] = i * level.step_fm * s.value[0][m];
    if (i == nodes - 1) {
      dw = s.value[0][m] + a * s.derivative[0][m];
      PlaneKinematics(c, cos_fi, p, theta_p, sign, q);
      phit = SourceMomentumAmplitude(plan, q);
    }
  }
  const double eps = kHbarC_MeV_fm * kHbarC_MeV_fm * p * p / (2.0 * level.mu_MeV);
  const EntranceVertexValue got = EntranceVertex(level, w, dw, eps, B);
  const EntranceVertexValue pw = PlaneWaveEntranceVertex(level, l, p);
  const std::vector<double> y = SphericalHarmonicsAtPhiZero(l, l, std::cos(theta_p));
  const cplx factor = 4.0 * kPi * phit * IPower(l) * y[l * (l + 1) + m] * (sign < 0 && m % 2 ? -1.0 : 1.0);
  const double scale =
      std::max({std::abs(factor * pw.volume), std::abs(factor * pw.surface), std::abs(factor * pw.interior)});
  Check("EntranceVertex from the projection: volume term", std::abs(got.volume - factor * pw.volume) / scale, 1e-8);
  Check("EntranceVertex from the projection: surface term", std::abs(got.surface - factor * pw.surface) / scale, 1e-8);
  Check("EntranceVertex from the projection: interior term", std::abs(got.interior - factor * pw.interior) / scale,
        1e-8);
}

// The engine's 19F(d,n) point-Coulomb cases: set-up and Gram entries at the spectatorAngle direction, as the base
// build's ThmDwVertex printed them (work/dwvertex/driver/engine_numbers.cpp of the branch that moved the projection).
struct EngineCase {
  int l;
  double mu_aa, e_aa, mu_sf, e_sf, step_aa, step_sf, k_sf_max, alpha, beta, kappa, mu_sx, bind, x_fi, q;
  double g[4];
};
const double kMuAA = 1.8205512145724445, kEAA = 5.271829850273158, kMuSF = 0.96023875457417629;
const double kAlpha = 0.94963804421629217, kBeta = 0.50093779950161421, kKappa = 0.23160616098692202;
const double kMuSx = 0.5039851067819181, kBind = 2.2245662705769096;
const EngineCase kEngineCases[] = {
    {0, kMuAA, kEAA, kMuSF, 2.7233635796962483, 0.02, 0.02, 0.35630934440347073, kAlpha, kBeta, kKappa, kMuSx, kBind,
     1.0, 0.014264322530445985,
     {0.31510015023488164, 0.024096617091680119, 0.0073599880518974824, 0.086825562143972024}},
    {1, kMuAA, kEAA, kMuSF, 2.8345635796962485, 0.02, 0.02, 0.36340772658248438, kAlpha, kBeta, kKappa, kMuSx, kBind,
     1.0, 0.021413618814143432,
     {0.039526617978248482, 0.032460304728206456, 0.035562645010152362, -0.0042830298074184739}},
    {3, kMuAA, kEAA, kMuSF, 2.8345635796962485, 0.02, 0.02, 0.36340772658248438, kAlpha, kBeta, kKappa, kMuSx, kBind,
     1.0, 0.021413618814143432,
     {0.00022478536355824329, 0.0019567757905444663, 0.00066312293350615956, -1.1069441616096671e-05}},
};

// The Gram entries (S_lm, a S_lm')/(4 pi phi~) summed over m, as the engine builds them.
void Gram(const SourceProjectionPlan &plan, const SourceAmplitudes &amp, double x_fi, double q, double g[4]) {
  const double phit = SourceMomentumAmplitude(plan, q);
  const SourceMultipoles s = SourceMultipolesAt(plan, amp, x_fi);
  const int l = plan.setup.l_values[0];
  const double a = plan.setup.radius_fm;
  double g11 = 0.0, g22 = 0.0;
  cplx g12(0.0, 0.0);
  for (int m = 0; m <= l; m++) {
    const cplx sm = s.value[0][m] / (4.0 * kPi * phit), dm = s.derivative[0][m] * (a / (4.0 * kPi * phit));
    const double mult = m ? 2.0 : 1.0;
    g11 += mult * std::norm(sm);
    g22 += mult * std::norm(dm);
    g12 += mult * std::conj(sm) * dm;
  }
  const double f = 4.0 * kPi / (2.0 * l + 1.0);
  g[0] = f * g11;
  g[1] = f * g22;
  g[2] = f * g12.real();
  g[3] = f * g12.imag();
}

SourceSetup EngineSetup(const EngineCase &c, WaveChannel &exit) {
  SourceSetup s;
  BoundStateRequest bound;
  bound.pair.z1 = 0;
  bound.pair.z2 = 1;
  bound.pair.mu_u = c.mu_sx;
  bound.binding_MeV = c.bind;
  s.overlap.phi = [bound](double u) {
    return u > 0.0 ? AnalyticTailAt(TailShape::Yukawa, bound, 1.0, 0.0, u) / u : 0.0;
  };
  s.overlap.kappa_fm = c.kappa;
  s.entrance.pair.z1 = 1;
  s.entrance.pair.z2 = 9;
  s.entrance.pair.mu_u = c.mu_aa;
  s.entrance.energy_MeV = c.e_aa;
  s.entrance_step_fm = c.step_aa;
  s.exit_step_fm = c.step_sf;
  s.exit_k_max_fm = c.k_sf_max;
  s.alpha = c.alpha;
  s.beta = c.beta;
  s.radius_fm = 5.136;
  s.l_values = {c.l};
  exit.pair.z1 = 0;
  exit.pair.z2 = 10;
  exit.pair.mu_u = c.mu_sf;
  exit.energy_MeV = c.e_sf;
  return s;
}

double GramDistance(const double *g, const double *want) {
  const double scale = std::max(std::fabs(want[0]), std::fabs(want[1]));
  double d = 0.0;
  for (int i = 0; i < 4; i++) d = std::max(d, std::fabs(g[i] - want[i]) / scale);
  return d;
}

void TestEngineNumbers() {
  for (const EngineCase &c : kEngineCases) {
    WaveChannel exit;
    const SourceSetup setup = EngineSetup(c, exit);
    const SourceProjectionPlan plan = PrepareSourceProjection(setup, SourceProjectionRequest{});
    const std::vector<SourceAmplitudes> amp = ProjectSource(plan, {exit});
    double g[4];
    Gram(plan, amp[0], c.x_fi, c.q, g);
    char label[96];
    std::snprintf(label, sizeof label, "the engine's former Gram entries, 19F(d,n) l=%d", c.l);
    Check(label, GramDistance(g, c.g), 1e-9);

    // The blocks in reverse order: the same numbers, bit for bit.
    std::vector<SourceExitWaves> exits = {SolveSourceExit(plan, exit)};
    const int blocks = SourceSumBlocks(plan, 1);
    const std::size_t size = SourceBlockSize(plan, 1);
    std::vector<cplx> sums(size * blocks, cplx(0.0, 0.0));
    for (int blk = blocks - 1; blk >= 0; blk--) SumSourceBlock(plan, exits, blk, blocks, &sums[blk * size]);
    const SourceAmplitudes again = CombineSourceBlocks(plan, exits, sums, blocks)[0];
    bool identical = again.value == amp[0].value && again.derivative == amp[0].derivative;
    std::snprintf(label, sizeof label, "blocks in reverse order bit for bit, l=%d (%d blocks)", c.l, blocks);
    Check(label, identical ? 0.0 : 1.0, 0.0);

    // Refinement: twice the nodes, partial waves + 20, the tables' steps halved.
    SourceSetup fine = setup;
    fine.entrance_step_fm *= 0.5;
    fine.exit_step_fm *= 0.5;
    SourceProjectionRequest r;
    r.nodes_per_panel *= 2;
    r.angle_nodes_per_l *= 2.0;
    r.angle_nodes_extra *= 2;
    r.partial_wave_margin += 20;
    const SourceProjectionPlan fine_plan = PrepareSourceProjection(fine, r);
    double gf[4];
    Gram(fine_plan, ProjectSource(fine_plan, {exit})[0], c.x_fi, c.q, gf);
    std::snprintf(label, sizeof label, "refined (nodes x 2, L + 20, steps / 2), l=%d", c.l);
    Check(label, GramDistance(g, gf), 1e-7);
    if (verbose) std::printf("      refined - present: %.3e\n", GramDistance(g, gf));
  }
}

}  // namespace

int main(int argc, char **argv) {
  verbose = argc > 2 && std::string(argv[2]) == "-v";
  if (argc > 1 && std::string(argv[1]) == "-v") verbose = true;
  TestNumerics();
  TestClosedForm();
  TestEntranceVertex();
  TestEngineNumbers();
  if (failures) {
    std::printf("transfer_source_projection: %d check(s) failed\n", failures);
    return 1;
  }
  std::printf("transfer_source_projection: all checks passed\n");
  return 0;
}
