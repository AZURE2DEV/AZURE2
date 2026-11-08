// physics/external_capture against hand values and identities (docs/upgrade/external_capture.md section 5).
#include "azure2/physics/Constants.h"
#include "azure2/physics/ExternalCapture.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>

using namespace azure2::physics;

namespace {
int failures = 0;
void check(bool ok, const char *what, double got = 0.0, double want = 0.0) {
  if (!ok) {
    std::printf("FAIL %s: got %.17g want %.17g\n", what, got, want);
    failures++;
  }
}
bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }
}  // namespace

int main() {
  // 1. N_f, Eq. (1): two channels, by hand.  2 mu a / hbar^2 = 2 * 3 u * 931.49410242 MeV/u * 5 fm / (197.3269804)^2.
  {
    const std::vector<BoundChannel> ch = {{5.0, 3.0, 0.4, 2.5}, {4.0, 1.5, -0.2, 0.0}};
    const double k = 2.0 * 3.0 * kAtomicMassUnit_MeV * 5.0 / (kHbarC_MeV_fm * kHbarC_MeV_fm);
    const double want = 1.0 + k * 0.16 * 2.5;
    check(near(bound_state_normalization(1.0, ch), want, 1e-15), "N_f Brune", bound_state_normalization(1.0, ch), want);
    check(near(bound_state_normalization(0.7, ch), want - 0.3, 1e-15), "N_f Park n0", bound_state_normalization(0.7, ch),
          want - 0.3);
    check(near(inverse_sqrt_normalization(1.0, ch), 1.0 / std::sqrt(want), 1e-15), "N_f^-1/2");
    bool threw = false;
    try {
      inverse_sqrt_normalization(-1.0, {});
    } catch (const std::domain_error &) {
      threw = true;
    }
    check(threw, "N_f <= 0 throws");
  }

  // 2. Selection rules.  s = 1/2 everywhere.  E1: p-wave (l=1) to s-wave J_f=1/2 allowed, l=2 to l_f=0 not; E2:
  // d-wave to s-wave allowed, p-wave to s-wave not (l_f = 0 outside |1-2|..3); M1: l = l_f with the spin rule.
  check(couples(Multipole::E, 1, 1, 0.5, 0, 0.5, 0.5), "E1 p->s");
  check(!couples(Multipole::E, 1, 2, 0.5, 0, 0.5, 0.5), "E1 d->s forbidden");
  check(couples(Multipole::E, 2, 2, 0.5, 0, 0.5, 0.5), "E2 d->s");
  check(!couples(Multipole::E, 2, 1, 0.5, 0, 0.5, 0.5), "E2 p->s forbidden");
  check(!couples(Multipole::E, 1, 1, 0.5, 0, 0.5, 1.5), "E1 p->s J_f=3/2 forbidden");
  check(!couples(Multipole::E, 1, 1, 0.5, 0, 1.5, 0.5), "E1 s != s_f forbidden");
  check(couples(Multipole::M, 1, 0, 0.5, 0, 1.5, 1.5), "M1 spin flip s->s");
  check(!couples(Multipole::E, 1, 0, 0.5, 0, 1.5, 1.5), "E1 l=l_f forbidden");

  // 3. Eq. (4) with a stub integral: two entrance channels (partition 1, l = 1 and l = 3) and two final channels
  // (partition 1 l_f = 0; partition 2 l_f = 0).  E1 to J_f = 1/2: only (l=1, partition 1) couples.
  {
    const std::vector<CaptureChannel> lev = {{1, 1, 0.5, 0.3}, {1, 3, 0.5, 0.9}};
    const std::vector<FinalChannel> fin = {{1, 0, 0.5, 2.0}, {2, 0, 0.5, 7.0}};
    int calls = 0;
    const RadialIntegral I = [&](std::size_t c, std::size_t cf) {
      calls++;
      return std::complex<double>(1.0 + c, 0.5 + cf);
    };
    const std::complex<double> g = external_reduced_width(Multipole::E, 1, 0.5, lev, fin, I);
    const std::complex<double> want = -2.0 * 0.3 * std::complex<double>(1.0, 0.5);
    check(std::abs(g - want) < 1e-15, "gamma^ext", g.real(), want.real());
    check(calls == 1, "integral called for coupled pairs only", calls, 1);
    check(external_reduced_width(Multipole::E, 3, 0.5, lev, fin, I) == std::complex<double>(0.0, 0.0), "E3 is zero");
  }

  // 4. capture_mask: diagonal T gives P; Pi^2 = Pi; and with A = T^T A~ T the masked formal bilinear equals the
  // masked Brune one (gamma~ = T gamma), for a random non-orthogonal T.
  {
    const std::size_t n = 3;
    const std::vector<double> D = {2.0, 0, 0, 0, -0.5, 0, 0, 0, 1.0};
    const std::vector<double> Pi = capture_mask(D, n, {0, 1, 0});
    const std::vector<double> P = {1, 0, 0, 0, 0, 0, 0, 0, 1};
    double worst = 0.0;
    for (std::size_t i = 0; i < n * n; i++) worst = std::max(worst, std::fabs(Pi[i] - P[i]));
    check(worst < 1e-15, "diagonal T gives P", worst, 0.0);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> T(n * n), At(n * n), g1(n), g2(n);
    for (double &x : T) x = u(rng);
    for (std::size_t i = 0; i < n; i++) T[i * n + i] += 2.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j <= i; j++) At[i * n + j] = At[j * n + i] = u(rng);
    for (std::size_t i = 0; i < n; i++) g1[i] = u(rng), g2[i] = u(rng);
    const std::vector<unsigned char> m = {1, 0, 1};
    const std::vector<double> Q = capture_mask(T, n, m);
    double idem = 0.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        double s = 0.0;
        for (std::size_t k = 0; k < n; k++) s += Q[i * n + k] * Q[k * n + j];
        idem = std::max(idem, std::fabs(s - Q[i * n + j]));
      }
    check(idem < 1e-12, "Pi^2 = Pi", idem, 0.0);
    auto mul = [&](const std::vector<double> &M, const std::vector<double> &v) {
      std::vector<double> r(n, 0.0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) r[i] += M[i * n + j] * v[j];
      return r;
    };
    auto bil = [&](const std::vector<double> &M, const std::vector<double> &a, const std::vector<double> &b) {
      const std::vector<double> Mb = mul(M, b);
      double s = 0.0;
      for (std::size_t i = 0; i < n; i++) s += a[i] * Mb[i];
      return s;
    };
    std::vector<double> A(n * n, 0.0);  // T^T A~ T
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++)
        for (std::size_t k = 0; k < n; k++)
          for (std::size_t l = 0; l < n; l++) A[i * n + j] += T[k * n + i] * At[k * n + l] * T[l * n + j];
    std::vector<double> b1 = mul(T, g1), b2 = mul(T, g2);
    for (std::size_t i = 0; i < n; i++)
      if (m[i]) b1[i] = b2[i] = 0.0;
    const double brune = bil(At, b1, b2), formal = bil(A, mul(Q, g1), mul(Q, g2));
    check(near(formal, brune, 1e-12), "masked bilinear, formal = Brune", formal, brune);
    // Without the mask the bilinear is the same in both bases; a wrong transpose fails the masked one at order 1.
    const double unmaskedB = bil(At, mul(T, g1), mul(T, g2)), unmaskedF = bil(A, g1, g2);
    check(near(unmaskedF, unmaskedB, 1e-12), "bilinear, formal = Brune", unmaskedF, unmaskedB);
  }

  if (failures == 0) std::printf("physics_external_capture_unit: all checks passed\n");
  return failures ? 1 : 0;
}
