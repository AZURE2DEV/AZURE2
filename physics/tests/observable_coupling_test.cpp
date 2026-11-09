// The 9j of AngularCoupling.h and the observable coupling coefficients of ObservableCoupling.h (docs/upgrade/
// cnuc_split.md section 9), checked against things that do not use the code under test:
//
//  1. the 9j with a zero argument against its reduction to a 6j (Edmonds 6.4.14);
//  2. the 9j against its definition as a sum over all projections of six 3j (Edmonds 6.4.1), a route that shares
//     only wigner_3j with the sum over x of three 6j, for all 2j <= 4;
//  3. the orthogonality of the 9j (Edmonds 6.4.8);
//  4. the particle interference coefficient (9.1) against the angular distribution summed in the m-scheme: for
//     pathways (l1 s J1 -> l1' s') and (l2 s J2 -> l2' s') with entrance amplitude l^ (l 0 s nu | J nu) and exit
//     amplitude (l' m' s' nu' | J nu) Y_l'm'(theta, 0), the sum over nu, nu' of the products is
//     (1/pi) sum_k c_k P_k(cos theta), for l <= 3, 2J <= 9 and two spin pairs;
//  5. the photon interference coefficient (9.2) the same way, with the photon amplitude
//     L^ (L mu jf mf | J nu) d^L_{mu lambda}(theta), times lambda for a magnetic multipole, summed over nu, mf and
//     the helicity lambda = +-1: 8 sum_k c_k P_k(cos theta), for L <= 2, E and M, 2jf <= 4.
// Y_lm, P_k and d^L are computed here from their own formulas.  A transposed 6j argument, a wrong phase or a
// missing hat fails at least one of these by more than 1e-3.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <string>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/ObservableCoupling.h"

using namespace azure2::physics;

namespace {

int failures = 0;
long checks = 0;
double worst = 0.0;

void close(const std::string &what, double got, double want, double tol) {
  checks++;
  const double d = std::fabs(got - want);
  if (d > worst) worst = d;
  if (!(d <= tol)) {
    std::printf("  FAIL  %s: got %.17g, want %.17g (diff %.3g)\n", what.c_str(), got, want, d);
    failures++;
  }
}

double Phase(int twice) { return ((std::abs(twice) / 2) % 2 == 0) ? 1.0 : -1.0; }

double Fact(int n) { return std::tgamma(n + 1.0); }

// Associated Legendre P_l^m(x), m >= 0, with the Condon-Shortley phase, by the upward recurrence in l.
double Plm(int l, int m, double x) {
  double pmm = 1.0;
  const double s = std::sqrt(1.0 - x * x);
  for (int i = 1; i <= m; i++) pmm *= -(2.0 * i - 1.0) * s;
  if (l == m) return pmm;
  double p1 = x * (2.0 * m + 1.0) * pmm;
  if (l == m + 1) return p1;
  double p0 = pmm;
  for (int ll = m + 2; ll <= l; ll++) {
    const double p2 = ((2.0 * ll - 1.0) * x * p1 - (ll + m - 1.0) * p0) / (ll - m);
    p0 = p1;
    p1 = p2;
  }
  return p1;
}

// Y_lm(theta, 0).
double Ylm(int l, int m, double theta) {
  const int am = std::abs(m);
  double v = std::sqrt((2.0 * l + 1.0) / (4.0 * M_PI) * Fact(l - am) / Fact(l + am)) * Plm(l, am, std::cos(theta));
  if (m < 0 && am % 2) v = -v;
  return v;
}

// Wigner d^j_{m'm}(beta) for integer j (Wigner's formula).
double SmallD(int j, int mp, int m, double beta) {
  double s = 0.0;
  for (int k = 0; k <= 2 * j; k++) {
    const int a1 = j + m - k, a3 = j - k - mp, a4 = k - m + mp;
    if (a1 < 0 || a3 < 0 || a4 < 0) continue;
    s += Phase(2 * (k - m + mp)) * std::pow(std::cos(beta / 2), 2 * j + m - mp - 2 * k) *
         std::pow(std::sin(beta / 2), 2 * k - m + mp) / (Fact(a1) * Fact(k) * Fact(a3) * Fact(a4));
  }
  return std::sqrt(Fact(j + mp) * Fact(j - mp) * Fact(j + m) * Fact(j - m)) * s;
}

// 9j as the sum over all projections of the six 3j of its rows and columns (Edmonds 6.4.1).
double NineJFromThreeJ(const int j[9]) {
  double s = 0.0;
  for (int m0 = -j[0]; m0 <= j[0]; m0 += 2)
    for (int m1 = -j[1]; m1 <= j[1]; m1 += 2) {
      const int m2 = -m0 - m1;
      if (std::abs(m2) > j[2]) continue;
      for (int m3 = -j[3]; m3 <= j[3]; m3 += 2)
        for (int m4 = -j[4]; m4 <= j[4]; m4 += 2) {
          const int m5 = -m3 - m4;
          if (std::abs(m5) > j[5]) continue;
          const int m6 = -m0 - m3, m7 = -m1 - m4, m8 = -m2 - m5;
          if (std::abs(m6) > j[6] || std::abs(m7) > j[7] || std::abs(m8) > j[8]) continue;
          if (m6 + m7 + m8 != 0) continue;
          s += wigner_3j(j[0], j[1], j[2], m0, m1, m2) * wigner_3j(j[3], j[4], j[5], m3, m4, m5) *
               wigner_3j(j[6], j[7], j[8], m6, m7, m8) * wigner_3j(j[0], j[3], j[6], m0, m3, m6) *
               wigner_3j(j[1], j[4], j[7], m1, m4, m7) * wigner_3j(j[2], j[5], j[8], m2, m5, m8);
        }
    }
  return s;
}

}  // namespace

int main() {
  std::printf("9j with a zero\n");
  // {a b e; c d e; f f 0} = (-1)^(b + c + e + f) / sqrt((2e+1)(2f+1)) {a b e; d c f}
  for (int a = 0; a <= 6; a++)
    for (int b = 0; b <= 6; b++)
      for (int e = 0; e <= 6; e++)
        for (int c = 0; c <= 6; c++)
          for (int d = 0; d <= 6; d++)
            for (int f = 0; f <= 6; f++) {
              const double want = wigner_6j(a, b, e, d, c, f) == 0.0
                  ? 0.0
                  : Phase(b + c + e + f) / std::sqrt((e + 1.0) * (f + 1.0)) * wigner_6j(a, b, e, d, c, f);
              close("9j zero", wigner_9j(a, b, e, c, d, e, f, f, 0), want, 1e-14);
            }
  const long zeroChecks = checks;
  const double zeroWorst = worst;

  std::printf("9j vs six-3j sum\n");
  worst = 0.0;
  long nonzero = 0;
  {
    int j[9];
    for (j[0] = 0; j[0] <= 4; j[0]++)
      for (j[1] = 0; j[1] <= 4; j[1]++)
        for (j[2] = 0; j[2] <= 4; j[2]++) {
          if (!triangle(j[0], j[1], j[2])) continue;
          for (j[3] = 0; j[3] <= 4; j[3]++)
            for (j[4] = 0; j[4] <= 4; j[4]++)
              for (j[5] = 0; j[5] <= 4; j[5]++) {
                if (!triangle(j[3], j[4], j[5])) continue;
                for (j[6] = 0; j[6] <= 4; j[6]++)
                  for (j[7] = 0; j[7] <= 4; j[7]++)
                    for (j[8] = 0; j[8] <= 4; j[8]++) {
                      if (!triangle(j[6], j[7], j[8])) continue;
                      const double got = wigner_9j(j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8]);
                      if (got != 0.0) nonzero++;
                      close("9j m-scheme", got, NineJFromThreeJ(j), 1e-14);
                    }
              }
        }
  }
  const long sumChecks = checks - zeroChecks;
  const double sumWorst = worst;

  std::printf("9j orthogonality\n");
  worst = 0.0;
  // sum_{j12, j34} (2j12+1)(2j34+1)(2j13+1)(2j24+1) {j1 j2 j12; j3 j4 j34; j13 j24 J}{j1 j2 j12; j3 j4 j34; j13' j24' J}
  //   = delta(j13 j13') delta(j24 j24')
  for (int j1 = 0; j1 <= 3; j1++)
    for (int j2 = 0; j2 <= 3; j2++)
      for (int j3 = 0; j3 <= 3; j3++)
        for (int j4 = 0; j4 <= 3; j4++)
          for (int tJ = 0; tJ <= 6; tJ++)
            for (int j13 = std::abs(j1 - j3); j13 <= j1 + j3; j13 += 2)
              for (int j24 = std::abs(j2 - j4); j24 <= j2 + j4; j24 += 2)
                for (int j13p = std::abs(j1 - j3); j13p <= j1 + j3; j13p += 2)
                  for (int j24p = std::abs(j2 - j4); j24p <= j2 + j4; j24p += 2) {
                    if (!triangle(j13, j24, tJ) || !triangle(j13p, j24p, tJ)) continue;
                    double s = 0.0;
                    for (int j12 = std::abs(j1 - j2); j12 <= j1 + j2; j12 += 2)
                      for (int j34 = std::abs(j3 - j4); j34 <= j3 + j4; j34 += 2)
                        s += (j12 + 1.0) * (j34 + 1.0) * std::sqrt((j13 + 1.0) * (j24 + 1.0) * (j13p + 1.0) * (j24p + 1.0)) *
                             wigner_9j(j1, j2, j12, j3, j4, j34, j13, j24, tJ) *
                             wigner_9j(j1, j2, j12, j3, j4, j34, j13p, j24p, tJ);
                    close("9j orthogonality", s, (j13 == j13p && j24 == j24p) ? 1.0 : 0.0, 1e-13);
                  }
  const long orthChecks = checks - zeroChecks - sumChecks;
  const double orthWorst = worst;

  const double thetas[3] = {0.4, 1.3, 2.6};

  std::printf("particle coefficient vs m-scheme\n");
  worst = 0.0;
  long particleCases = 0;
  const int spinPairs[2][4] = {{1, 3, 1, 3}, {1, 0, 2, 1}};  // 2I1, 2I2 entrance; 2I1', 2I2' exit
  for (const auto &sp : spinPairs)
    for (int ts = std::abs(sp[0] - sp[1]); ts <= sp[0] + sp[1]; ts += 2)
      for (int tsp = std::abs(sp[2] - sp[3]); tsp <= sp[2] + sp[3]; tsp += 2)
        for (int l1 = 0; l1 <= 3; l1++)
          for (int l1p = 0; l1p <= 3; l1p++)
            for (int l2 = 0; l2 <= 3; l2++)
              for (int l2p = 0; l2p <= 3; l2p++)
                for (int tJ1 = 0; tJ1 <= 9; tJ1++)
                  for (int tJ2 = 0; tJ2 <= 9; tJ2++) {
                    if (!triangle(2 * l1, ts, tJ1) || !triangle(2 * l1p, tsp, tJ1) || !triangle(2 * l2, ts, tJ2) ||
                        !triangle(2 * l2p, tsp, tJ2))
                      continue;
                    if ((l1 + l1p + l2 + l2p) % 2) continue;  // the two levels have one parity change between them
                    for (double th : thetas) {
                      double S = 0.0;
                      for (int nu = -ts; nu <= ts; nu += 2)
                        for (int nup = -tsp; nup <= tsp; nup += 2) {
                          const int tm = nu - nup;
                          if (std::abs(tm) > 2 * l1p || std::abs(tm) > 2 * l2p) continue;
                          const double a1 = std::sqrt(2. * l1 + 1.) * clebsch_gordan(2 * l1, 0, ts, nu, tJ1, nu) *
                                            clebsch_gordan(2 * l1p, tm, tsp, nup, tJ1, nu) * Ylm(l1p, tm / 2, th);
                          const double a2 = std::sqrt(2. * l2 + 1.) * clebsch_gordan(2 * l2, 0, ts, nu, tJ2, nu) *
                                            clebsch_gordan(2 * l2p, tm, tsp, nup, tJ2, nu) * Ylm(l2p, tm / 2, th);
                          S += a1 * a2;
                        }
                      double B = 0.0;
                      for (int k = 0; k <= 8; k++)
                        B += particle_interference(l1, l1p, tJ1, l2, l2p, tJ2, ts, tsp, k) * Plm(k, 0, std::cos(th));
                      particleCases++;
                      close("particle l=" + std::to_string(l1) + std::to_string(l1p) + std::to_string(l2) +
                                std::to_string(l2p) + " 2J=" + std::to_string(tJ1) + "," + std::to_string(tJ2),
                            S, B / M_PI, 1e-13);
                    }
                  }
  const double particleWorst = worst;

  std::printf("photon coefficient vs m-scheme\n");
  worst = 0.0;
  long photonCases = 0;
  for (int ta1 : {1, 2})
    for (int ta2 : {0, 1, 3})
      for (int tjf = 0; tjf <= 4; tjf++)
        for (int ts = std::abs(ta1 - ta2); ts <= ta1 + ta2; ts += 2)
          for (int l1 = 0; l1 <= 3; l1++)
            for (int l2 = 0; l2 <= 3; l2++)
              for (int L1 = 1; L1 <= 2; L1++)
                for (int L2 = 1; L2 <= 2; L2++)
                  for (int e1 = 0; e1 < 2; e1++)  // 1 electric, 0 magnetic
                    for (int e2 = 0; e2 < 2; e2++)
                      for (int tJ1 = 0; tJ1 <= 9; tJ1++)
                        for (int tJ2 = 0; tJ2 <= 9; tJ2++) {
                          if (!triangle(2 * l1, ts, tJ1) || !triangle(2 * l2, ts, tJ2) || !triangle(2 * L1, tjf, tJ1) ||
                              !triangle(2 * L2, tjf, tJ2))
                            continue;
                          // One final-state parity: (-1)^l1 (-1)^(L1 + 1 - e1) = (-1)^l2 (-1)^(L2 + 1 - e2).
                          if ((l1 + L1 + 1 - e1) % 2 != (l2 + L2 + 1 - e2) % 2) continue;
                          for (double th : thetas) {
                            double S = 0.0;
                            for (int nu = -ts; nu <= ts; nu += 2)
                              for (int mf = -tjf; mf <= tjf; mf += 2) {
                                const int tmu = nu - mf;
                                if (std::abs(tmu) > 2 * L1 || std::abs(tmu) > 2 * L2) continue;
                                for (int lam : {-1, 1}) {
                                  const double a1 = std::sqrt(2. * l1 + 1.) * std::sqrt(2. * L1 + 1.) *
                                                    clebsch_gordan(2 * l1, 0, ts, nu, tJ1, nu) *
                                                    clebsch_gordan(2 * L1, tmu, tjf, mf, tJ1, nu) *
                                                    SmallD(L1, tmu / 2, lam, th) * (e1 ? 1.0 : (double)lam);
                                  const double a2 = std::sqrt(2. * l2 + 1.) * std::sqrt(2. * L2 + 1.) *
                                                    clebsch_gordan(2 * l2, 0, ts, nu, tJ2, nu) *
                                                    clebsch_gordan(2 * L2, tmu, tjf, mf, tJ2, nu) *
                                                    SmallD(L2, tmu / 2, lam, th) * (e2 ? 1.0 : (double)lam);
                                  S += a1 * a2;
                                }
                              }
                            // The engine's selection rule: l1 + l2 + k and L1 + L2 + e1 + e2 + k even.
                            double B = 0.0;
                            for (int k = 0; k <= 6; k++) {
                              if ((l1 + l2 + k) % 2 || (L1 + L2 + e1 + e2 + k) % 2) continue;
                              B += photon_interference(l1, L1, tJ1, l2, L2, tJ2, ts, tjf, k) * Plm(k, 0, std::cos(th));
                            }
                            photonCases++;
                            close("photon l=" + std::to_string(l1) + std::to_string(l2) + " L=" + std::to_string(L1) +
                                      std::to_string(L2) + " EM=" + std::to_string(e1) + std::to_string(e2),
                                  S, 8.0 * B, 1e-12);
                          }
                        }
  const double photonWorst = worst;

  std::printf("  9j with a zero:         %ld checks, largest difference %.2g\n", zeroChecks, zeroWorst);
  std::printf("  9j vs six-3j sum:       %ld checks (%ld non-zero), largest difference %.2g\n", sumChecks, nonzero,
              sumWorst);
  std::printf("  9j orthogonality:       %ld checks, largest difference %.2g\n", orthChecks, orthWorst);
  std::printf("  particle vs m-scheme:   %ld cases, largest difference %.2g\n", particleCases, particleWorst);
  std::printf("  photon vs m-scheme:     %ld cases, largest difference %.2g\n", photonCases, photonWorst);
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
