#include "azure2/physics/AngularCoupling.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace azure2 {
namespace physics {

namespace {

// n! for n = 0 .. 2 kMaxTwiceJ + 1 (the largest argument of the sums below is j1 + j2 + j3 + 1 <= 3 jmax + 1 for
// 3j and the 6j's t + 1 <= 4 jmax + 1, both below 2 kMaxTwiceJ + 2).  In long double: 161! ~ 1e289 fits.
constexpr int kMaxFactorial = 2 * kMaxTwiceJ + 2;

const std::array<long double, kMaxFactorial + 1> &Factorials() {
  static const std::array<long double, kMaxFactorial + 1> table = [] {
    std::array<long double, kMaxFactorial + 1> f{};
    f[0] = 1.0L;
    for (int n = 1; n <= kMaxFactorial; n++) f[n] = f[n - 1] * (long double)n;
    return f;
  }();
  return table;
}

long double Fact(int n) {
  if (n < 0 || n > kMaxFactorial) throw std::out_of_range("angular coupling: factorial argument " + std::to_string(n));
  return Factorials()[n];
}

void CheckRange(int two_j) {
  if (two_j > kMaxTwiceJ)
    throw std::out_of_range("angular coupling: 2j = " + std::to_string(two_j) + " exceeds " +
                            std::to_string(kMaxTwiceJ));
}

// Delta(abc) = (a+b-c)! (a-b+c)! (-a+b+c)! / (a+b+c+1)!, arguments doubled and already triangle-checked.
long double Delta(int ta, int tb, int tc) {
  return Fact((ta + tb - tc) / 2) * Fact((ta - tb + tc) / 2) * Fact((-ta + tb + tc) / 2) / Fact((ta + tb + tc) / 2 + 1);
}

}  // namespace

bool triangle(int two_a, int two_b, int two_c) {
  if (two_a < 0 || two_b < 0 || two_c < 0) return false;
  if ((two_a + two_b + two_c) % 2 != 0) return false;
  return two_c >= std::abs(two_a - two_b) && two_c <= two_a + two_b;
}

double wigner_3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
  CheckRange(tj1);
  CheckRange(tj2);
  CheckRange(tj3);
  if (tm1 + tm2 + tm3 != 0) return 0.0;
  if (!triangle(tj1, tj2, tj3)) return 0.0;
  if (std::abs(tm1) > tj1 || std::abs(tm2) > tj2 || std::abs(tm3) > tj3) return 0.0;
  if ((tj1 + tm1) % 2 != 0 || (tj2 + tm2) % 2 != 0 || (tj3 + tm3) % 2 != 0) return 0.0;
  // Edmonds 3.7.3 (Racah):
  // (j1 j2 j3; m1 m2 m3) = (-1)^(j1 - j2 - m3) sqrt(Delta(j1 j2 j3) (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (j3+m3)!
  //   (j3-m3)!) sum_k (-1)^k / [k! (j3-j2+k+m1)! (j3-j1+k-m2)! (j1+j2-j3-k)! (j1-k-m1)! (j2-k+m2)!]
  const int a1 = (tj3 - tj2 + tm1) / 2;  // j3 - j2 + m1
  const int a2 = (tj3 - tj1 - tm2) / 2;  // j3 - j1 - m2
  const int b1 = (tj1 + tj2 - tj3) / 2;  // j1 + j2 - j3
  const int b2 = (tj1 - tm1) / 2;        // j1 - m1
  const int b3 = (tj2 + tm2) / 2;        // j2 + m2
  int kmin = 0;
  if (-a1 > kmin) kmin = -a1;
  if (-a2 > kmin) kmin = -a2;
  int kmax = b1;
  if (b2 < kmax) kmax = b2;
  if (b3 < kmax) kmax = b3;
  long double sum = 0.0L;
  for (int k = kmin; k <= kmax; k++) {
    const long double term = 1.0L / (Fact(k) * Fact(a1 + k) * Fact(a2 + k) * Fact(b1 - k) * Fact(b2 - k) * Fact(b3 - k));
    sum += (k % 2 == 0) ? term : -term;
  }
  const long double norm = Delta(tj1, tj2, tj3) * Fact((tj1 + tm1) / 2) * Fact((tj1 - tm1) / 2) * Fact((tj2 + tm2) / 2) *
      Fact((tj2 - tm2) / 2) * Fact((tj3 + tm3) / 2) * Fact((tj3 - tm3) / 2);
  const int phase = (tj1 - tj2 - tm3) / 2;  // j1 - j2 - m3 = (j1 + m1) + (j2 + m2) - 2 j2: an integer
  const long double value = std::sqrt(norm) * sum;
  return (double)((phase % 2 == 0) ? value : -value);
}

double clebsch_gordan(int tj1, int tm1, int tj2, int tm2, int tJ, int tM) {
  const double threej = wigner_3j(tj1, tj2, tJ, tm1, tm2, -tM);
  if (threej == 0.0) return 0.0;
  const int phase = (tj1 - tj2 + tM) / 2;  // j1 - j2 + M is an integer when the 3j is not zero
  const double value = std::sqrt((double)(tJ + 1)) * threej;
  return (phase % 2 == 0) ? value : -value;
}

double clebsch_gordan_integer(int j1, int m1, int j2, int m2, int J, int M) {
  if (m1 + m2 != M || J < std::abs(j1 - j2) || J > j1 + j2 || std::abs(m1) > j1 || std::abs(m2) > j2 ||
      std::abs(M) > J)
    return 0.0;
  auto lf = [](int n) { return lgammal((long double)n + 1.0L); };
  const long double pre = 0.5L * (logl(2.0L * J + 1.0L) + lf(J + j1 - j2) + lf(J - j1 + j2) + lf(j1 + j2 - J) -
                                  lf(j1 + j2 + J + 1) + lf(J + M) + lf(J - M) + lf(j1 - m1) + lf(j1 + m1) +
                                  lf(j2 - m2) + lf(j2 + m2));
  long double sum = 0.0L;
  const int kMin = std::max(0, std::max(j2 - J - m1, j1 - J + m2));
  const int kMax = std::min(j1 + j2 - J, std::min(j1 - m1, j2 + m2));
  for (int k = kMin; k <= kMax; k++) {
    const long double t = expl(pre - lf(k) - lf(j1 + j2 - J - k) - lf(j1 - m1 - k) - lf(j2 + m2 - k) -
                               lf(J - j2 + m1 + k) - lf(J - j1 - m2 + k));
    sum += (k % 2) ? -t : t;
  }
  return (double)sum;
}

double wigner_6j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
  CheckRange(tj1);
  CheckRange(tj2);
  CheckRange(tj3);
  CheckRange(tj4);
  CheckRange(tj5);
  CheckRange(tj6);
  // The four triads (j1 j2 j3), (j1 j5 j6), (j4 j2 j6), (j4 j5 j3).
  if (!triangle(tj1, tj2, tj3) || !triangle(tj1, tj5, tj6) || !triangle(tj4, tj2, tj6) || !triangle(tj4, tj5, tj3))
    return 0.0;
  // Edmonds 6.3.7 (Racah):
  // {j1 j2 j3; j4 j5 j6} = sqrt(Delta(j1j2j3) Delta(j1j5j6) Delta(j4j2j6) Delta(j4j5j3))
  //   sum_t (-1)^t (t+1)! / [(t-a1)! (t-a2)! (t-a3)! (t-a4)! (b1-t)! (b2-t)! (b3-t)!]
  // with a = the four triad sums and b = the three sums of opposite pairs.
  const int a1 = (tj1 + tj2 + tj3) / 2;
  const int a2 = (tj1 + tj5 + tj6) / 2;
  const int a3 = (tj4 + tj2 + tj6) / 2;
  const int a4 = (tj4 + tj5 + tj3) / 2;
  const int b1 = (tj1 + tj2 + tj4 + tj5) / 2;
  const int b2 = (tj2 + tj3 + tj5 + tj6) / 2;
  const int b3 = (tj3 + tj1 + tj6 + tj4) / 2;
  int tmin = a1;
  if (a2 > tmin) tmin = a2;
  if (a3 > tmin) tmin = a3;
  if (a4 > tmin) tmin = a4;
  int tmax = b1;
  if (b2 < tmax) tmax = b2;
  if (b3 < tmax) tmax = b3;
  long double sum = 0.0L;
  for (int t = tmin; t <= tmax; t++) {
    const long double term = Fact(t + 1) /
        (Fact(t - a1) * Fact(t - a2) * Fact(t - a3) * Fact(t - a4) * Fact(b1 - t) * Fact(b2 - t) * Fact(b3 - t));
    sum += (t % 2 == 0) ? term : -term;
  }
  const long double norm = Delta(tj1, tj2, tj3) * Delta(tj1, tj5, tj6) * Delta(tj4, tj2, tj6) * Delta(tj4, tj5, tj3);
  return (double)(std::sqrt(norm) * sum);
}

double racah_w(int ta, int tb, int tc, int td, int te, int tf) {
  const double sixj = wigner_6j(ta, tb, te, td, tc, tf);
  if (sixj == 0.0) return 0.0;
  // a + b + c + d is an integer when the 6j is not zero: (a b e) and (c d e) are triads.
  return (((ta + tb + tc + td) / 2) % 2 == 0) ? sixj : -sixj;
}

double wigner_9j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6, int tj7, int tj8, int tj9) {
  CheckRange(tj1);
  CheckRange(tj2);
  CheckRange(tj3);
  CheckRange(tj4);
  CheckRange(tj5);
  CheckRange(tj6);
  CheckRange(tj7);
  CheckRange(tj8);
  CheckRange(tj9);
  // Rows and columns are triads.
  if (!triangle(tj1, tj2, tj3) || !triangle(tj4, tj5, tj6) || !triangle(tj7, tj8, tj9) ||
      !triangle(tj1, tj4, tj7) || !triangle(tj2, tj5, tj8) || !triangle(tj3, tj6, tj9))
    return 0.0;
  // Edmonds 6.4.3:
  // {j1 j2 j3; j4 j5 j6; j7 j8 j9} = sum_x (-1)^(2x) (2x+1) {j1 j4 j7; j8 j9 x} {j2 j5 j8; j4 x j6} {j3 j6 j9; x j1 j2}
  // with x running over the values allowed by (j1 j9 x), (j8 j4 x) and (j6 j2 x).
  int txmin = std::abs(tj1 - tj9);
  if (std::abs(tj8 - tj4) > txmin) txmin = std::abs(tj8 - tj4);
  if (std::abs(tj6 - tj2) > txmin) txmin = std::abs(tj6 - tj2);
  int txmax = tj1 + tj9;
  if (tj8 + tj4 < txmax) txmax = tj8 + tj4;
  if (tj6 + tj2 < txmax) txmax = tj6 + tj2;
  long double sum = 0.0L;
  for (int tx = txmin; tx <= txmax; tx += 2) {
    const long double term = (long double)(tx + 1) * wigner_6j(tj1, tj4, tj7, tj8, tj9, tx) *
        wigner_6j(tj2, tj5, tj8, tj4, tx, tj6) * wigner_6j(tj3, tj6, tj9, tx, tj1, tj2);
    sum += (tx % 2 == 0) ? term : -term;
  }
  return (double)sum;
}

}  // namespace physics
}  // namespace azure2
