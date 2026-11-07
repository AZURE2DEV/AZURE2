// 3j, 6j and Clebsch-Gordan coefficients of physics/ (AngularCoupling.h) and the recoupling coefficients of
// level_parameters (recoupling_coefficient), checked against things that do not use the code under test:
//
//  1. exact values: tabulated closed forms (Edmonds Table 2; the stretched 3j; the 6j with a zero; {111;111} = 1/6,
//     {222;222} = -3/70, {1/2 1/2 1; 1/2 1/2 1} = 1/6);
//  2. orthogonality of the 3j (sum over m1, m2) and of the 6j (sum over x of (2x+1)(2j+1) {a b x; c d j}{a b x; c d j'});
//  3. the 6j against its definition as a sum of four 3j over all projections (Edmonds 6.2.4 / Racah), a route that
//     shares only wigner_3j with it;
//  4. every recoupling coefficient of the Convention record (LS heavy-first and l-first, jj light and heavy with
//     either orbital order) against the overlap of the two states built in the m-scheme from Clebsch-Gordan
//     products, for 2j1, 2j2 <= 4 and l <= 3: a wrong phase or a wrong 6j argument order fails it;
//  5. the iˡ sign (-1)^floor(l/2) and the orthogonality of each (alpha, l) block of the recoupling.
// A wrong factorial, a wrong sign in a Racah sum, a transposed 6j argument or a wrong phase of the recoupling fails
// at least one of these by more than 1e-3.
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/LevelParameters.h"

using namespace azure2::physics;

namespace {

int failures = 0;
long checks = 0;
double worst = 0.0;

void close(const std::string &what, double got, double want, double tol = 1e-14) {
  checks++;
  const double d = std::fabs(got - want);
  if (d > worst) worst = d;
  if (!(d <= tol)) {
    std::printf("  FAIL  %s: got %.17g, want %.17g (diff %.3g)\n", what.c_str(), got, want, d);
    failures++;
  }
}

double Phase(int twice) { return ((twice / 2) % 2 == 0) ? 1.0 : -1.0; }

// Product-basis vector of a three-spin coupled state.  Spins doubled: index 0 = j1, 1 = j2, 2 = l (2l).
// Tree: |((p q) x, r) J M> if xFirst, else |(r, (p q) x) J M>.
using MState = std::map<std::tuple<int, int, int>, double>;
MState Coupled(const int tj[3], int p, int q, int tx, int r, bool xFirst, int tJ, int tM) {
  MState out;
  for (int mp = -tj[p]; mp <= tj[p]; mp += 2)
    for (int mq = -tj[q]; mq <= tj[q]; mq += 2) {
      const int mx = mp + mq, mr = tM - mx;
      if (std::abs(mx) > tx || std::abs(mr) > tj[r]) continue;
      const double a = clebsch_gordan(tj[p], mp, tj[q], mq, tx, mx);
      const double c = xFirst ? clebsch_gordan(tx, mx, tj[r], mr, tJ, tM) : clebsch_gordan(tj[r], mr, tx, mx, tJ, tM);
      if (a == 0.0 || c == 0.0) continue;
      int m[3];
      m[p] = mp;
      m[q] = mq;
      m[r] = mr;
      out[std::make_tuple(m[0], m[1], m[2])] += a * c;
    }
  return out;
}

double Dot(const MState &a, const MState &b) {
  double s = 0.0;
  for (const auto &kv : a) {
    auto it = b.find(kv.first);
    if (it != b.end()) s += kv.second * it->second;
  }
  return s;
}

// 6j from four 3j (Edmonds 6.2.4 in the form of A. Messiah, Quantum Mechanics, App. C, Eq. C.36):
// {j1 j2 j3; j4 j5 j6} = sum_{all m} (-1)^(sum (j - m)) (j1 j2 j3; -m1 -m2 -m3)(j1 j5 j6; m1 -m5 m6)
//                                                       (j4 j2 j6; m4 m2 -m6)(j4 j5 j3; -m4 m5 m3)
double SixJFromThreeJ(int a, int b, int c, int d, int e, int f) {
  double s = 0.0;
  for (int m1 = -a; m1 <= a; m1 += 2)
    for (int m2 = -b; m2 <= b; m2 += 2) {
      const int m3 = -m1 - m2;  // first 3j: -m1 - m2 - m3 = 0
      if (std::abs(m3) > c) continue;
      for (int m5 = -e; m5 <= e; m5 += 2) {
        const int m6 = -m1 + m5;  // second: m1 - m5 + m6 = 0
        if (std::abs(m6) > f) continue;
        const int m4 = m6 - m2;  // third: m4 + m2 - m6 = 0
        if (std::abs(m4) > d) continue;
        if (-m4 + m5 + m3 != 0) continue;
        const double t = wigner_3j(a, b, c, -m1, -m2, -m3) * wigner_3j(a, e, f, m1, -m5, m6) *
            wigner_3j(d, b, f, m4, m2, -m6) * wigner_3j(d, e, c, -m4, m5, m3);
        const int ph = (a - m1) + (b - m2) + (c - m3) + (d - m4) + (e - m5) + (f - m6);
        s += Phase(ph) * t;
      }
    }
  return s;
}

}  // namespace

int main() {
  std::printf("exact values\n");
  close("(1 1 0; 0 0 0)", wigner_3j(2, 2, 0, 0, 0, 0), -1.0 / std::sqrt(3.0));
  close("(1 1 2; 0 0 0)", wigner_3j(2, 2, 4, 0, 0, 0), std::sqrt(2.0 / 15.0));
  close("(2 2 2; 0 0 0)", wigner_3j(4, 4, 4, 0, 0, 0), -std::sqrt(2.0 / 35.0));
  close("(1/2 1/2 1; 1/2 -1/2 0)", wigner_3j(1, 1, 2, 1, -1, 0), 1.0 / std::sqrt(6.0));
  close("(1/2 1/2 1/2 -1/2 | 1 0)", clebsch_gordan(1, 1, 1, -1, 2, 0), 1.0 / std::sqrt(2.0));
  close("(1/2 1/2 1/2 -1/2 | 0 0)", clebsch_gordan(1, 1, 1, -1, 0, 0), 1.0 / std::sqrt(2.0));
  close("(1/2 -1/2 1/2 1/2 | 0 0)", clebsch_gordan(1, -1, 1, 1, 0, 0), -1.0 / std::sqrt(2.0));
  close("(1 0 1/2 1/2 | 3/2 1/2)", clebsch_gordan(2, 0, 1, 1, 3, 1), std::sqrt(2.0 / 3.0));
  close("(1 0 1/2 1/2 | 1/2 1/2)", clebsch_gordan(2, 0, 1, 1, 1, 1), -std::sqrt(1.0 / 3.0));
  close("(1 1 1/2 -1/2 | 1/2 1/2)", clebsch_gordan(2, 2, 1, -1, 1, 1), std::sqrt(2.0 / 3.0));
  close("{1 1 1; 1 1 1}", wigner_6j(2, 2, 2, 2, 2, 2), 1.0 / 6.0);
  close("{2 2 2; 2 2 2}", wigner_6j(4, 4, 4, 4, 4, 4), -3.0 / 70.0);
  close("{1/2 1/2 1; 1/2 1/2 1}", wigner_6j(1, 1, 2, 1, 1, 2), 1.0 / 6.0);
  close("{1/2 1/2 1; 1/2 1/2 0}", wigner_6j(1, 1, 2, 1, 1, 0), 0.5);
  // Stretched 3j (Edmonds 3.7.12): (j1 j2 j1+j2; m1 m2 -M) =
  // (-1)^(j1-j2+M) [(2j1)! (2j2)! (J+M)! (J-M)! / ((2J+1)! (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)!)]^1/2
  for (int tj1 = 0; tj1 <= 8; tj1++)
    for (int tj2 = 0; tj2 <= 8; tj2++)
      for (int tm1 = -tj1; tm1 <= tj1; tm1 += 2)
        for (int tm2 = -tj2; tm2 <= tj2; tm2 += 2) {
          const int tJ = tj1 + tj2, tM = tm1 + tm2;
          auto fact = [](int n) { return std::tgamma(n + 1.0); };
          const double v = Phase(tj1 - tj2 + tM) *
              std::sqrt(fact(tj1) * fact(tj2) * fact((tJ + tM) / 2) * fact((tJ - tM) / 2) /
                        (fact(tJ + 1) * fact((tj1 + tm1) / 2) * fact((tj1 - tm1) / 2) * fact((tj2 + tm2) / 2) *
                         fact((tj2 - tm2) / 2)));
          close("stretched 3j", wigner_3j(tj1, tj2, tJ, tm1, tm2, -tM), v, 2e-15);
        }
  // 6j with a zero (Edmonds 6.3.2): {a b c; 0 c b} = (-1)^(a+b+c) / sqrt((2b+1)(2c+1)).
  for (int a = 0; a <= 8; a++)
    for (int b = 0; b <= 8; b++)
      for (int c = std::abs(a - b); c <= a + b; c += 2)
        close("{a b c; 0 c b}", wigner_6j(a, b, c, 0, c, b), Phase(a + b + c) / std::sqrt((b + 1.0) * (c + 1.0)));
  close("triangle violation is 0", wigner_3j(2, 2, 6, 0, 0, 0), 0.0);
  close("projection sum is 0", wigner_3j(2, 2, 2, 2, 0, 0), 0.0);
  const long exactChecks = checks;
  const double exactWorst = worst;

  std::printf("orthogonality\n");
  worst = 0.0;
  long before = checks;
  for (int tj1 = 0; tj1 <= 6; tj1++)
    for (int tj2 = 0; tj2 <= 6; tj2++)
      for (int tj3 = std::abs(tj1 - tj2); tj3 <= tj1 + tj2; tj3 += 2)
        for (int tj3p = std::abs(tj1 - tj2); tj3p <= tj1 + tj2; tj3p += 2)
          for (int tm3 = -std::min(tj3, tj3p); tm3 <= std::min(tj3, tj3p); tm3 += 2) {
            if ((tj3 + tm3) % 2 || (tj3p + tm3) % 2) continue;
            double s = 0.0;
            for (int tm1 = -tj1; tm1 <= tj1; tm1 += 2) {
              const int tm2 = -tm1 - tm3;
              s += wigner_3j(tj1, tj2, tj3, tm1, tm2, tm3) * wigner_3j(tj1, tj2, tj3p, tm1, tm2, tm3);
            }
            close("3j orthogonality", (tj3 + 1) * s, tj3 == tj3p ? 1.0 : 0.0, 1e-13);
          }
  for (int a = 0; a <= 6; a++)
    for (int b = 0; b <= 6; b++)
      for (int c = 0; c <= 6; c++)
        for (int d = 0; d <= 6; d++)
          for (int j = 0; j <= 8; j++)
            for (int jp = 0; jp <= 8; jp++) {
              if (!triangle(a, d, j) || !triangle(c, b, j) || !triangle(a, d, jp) || !triangle(c, b, jp)) continue;
              double s = 0.0;
              for (int x = 0; x <= 12; x++) s += (x + 1) * (j + 1) * wigner_6j(a, b, x, c, d, j) * wigner_6j(a, b, x, c, d, jp);
              close("6j orthogonality", s, j == jp ? 1.0 : 0.0, 1e-13);
            }
  const long orthChecks = checks - before;
  const double orthWorst = worst;

  std::printf("6j against four 3j\n");
  worst = 0.0;
  before = checks;
  for (int a = 0; a <= 4; a++)
    for (int b = 0; b <= 4; b++)
      for (int c = 0; c <= 4; c++)
        for (int d = 0; d <= 4; d++)
          for (int e = 0; e <= 4; e++)
            for (int f = 0; f <= 4; f++) {
              if (!triangle(a, b, c) || !triangle(a, e, f) || !triangle(d, b, f) || !triangle(d, e, c)) continue;
              close("6j vs 3j sum", wigner_6j(a, b, c, d, e, f), SixJFromThreeJ(a, b, c, d, e, f), 1e-13);
            }
  const long sixChecks = checks - before;
  const double sixWorst = worst;

  std::printf("recoupling coefficients against the m-scheme\n");
  worst = 0.0;
  before = checks;
  struct Variant {
    const char *name;
    Convention c;
  };
  std::vector<Variant> variants;
  {
    Convention c;
    variants.push_back({"LS internal", c});
    c.spin_order = SpinOrder::HeavyFirst;
    variants.push_back({"LS heavy first", c});
    c.spin_order = SpinOrder::LightFirst;
    c.orbital_order = OrbitalOrder::OrbitalFirst;
    variants.push_back({"LS l first", c});
    c.spin_order = SpinOrder::HeavyFirst;
    variants.push_back({"LS heavy, l first", c});
    for (SpinOrder so : {SpinOrder::LightFirst, SpinOrder::HeavyFirst})
      for (OrbitalOrder oo : {OrbitalOrder::SpinFirst, OrbitalOrder::OrbitalFirst}) {
        Convention j;
        j.coupling = CouplingScheme::JJ;
        j.spin_order = so;
        j.orbital_order = oo;
        variants.push_back({so == SpinOrder::LightFirst ? (oo == OrbitalOrder::SpinFirst ? "jj light (j1 l)" : "jj light (l j1)")
                                                        : (oo == OrbitalOrder::SpinFirst ? "jj heavy (j2 l)" : "jj heavy (l j2)"),
                            j});
      }
  }
  long mscheme = 0;
  for (const Variant &v : variants)
    for (int tj1 = 0; tj1 <= 4; tj1++)
      for (int tj2 = 0; tj2 <= 4; tj2++)
        for (int l = 0; l <= 3; l++)
          for (int tJ = 0; tJ <= tj1 + tj2 + 2 * l; tJ++) {
            const int tj[3] = {tj1, tj2, 2 * l};
            for (int ts = std::abs(tj1 - tj2); ts <= tj1 + tj2; ts += 2) {
              if (!triangle(ts, 2 * l, tJ)) continue;
              const MState in = Coupled(tj, 0, 1, ts, 2, true, tJ, tJ);  // |((j1 j2) s, l) J J>
              const bool jj = v.c.coupling == CouplingScheme::JJ;
              const bool light = v.c.spin_order == SpinOrder::LightFirst;
              const bool lFirst = v.c.orbital_order == OrbitalOrder::OrbitalFirst;
              // The project's x values and states.
              const int xa = light ? 0 : 1, xb = light ? 1 : 0;
              const int lo = jj ? std::abs(2 * l - tj[xa]) : std::abs(tj1 - tj2);
              const int hi = jj ? 2 * l + tj[xa] : tj1 + tj2;
              double norm = 0.0;
              for (int tx = lo; tx <= hi; tx += 2) {
                MState out;
                if (!jj) {
                  if (!triangle(tx, 2 * l, tJ)) continue;
                  // spins: (j1 j2) or (j2 j1); orbital: (s l) or (l s)
                  out = Coupled(tj, light ? 0 : 1, light ? 1 : 0, tx, 2, !lFirst, tJ, tJ);
                } else {
                  if (!triangle(tx, tj[xb], tJ)) continue;
                  // j = (ja l) or (l ja); J = (j, jb)
                  out = lFirst ? Coupled(tj, 2, xa, tx, xb, true, tJ, tJ) : Coupled(tj, xa, 2, tx, xb, true, tJ, tJ);
                }
                const double want = Dot(in, out);
                const double got = recoupling_coefficient(tj1, tj2, l, ts, tJ, v.c, tx);
                close(std::string(v.name) + " 2j1=" + std::to_string(tj1) + " 2j2=" + std::to_string(tj2) +
                          " l=" + std::to_string(l) + " 2s=" + std::to_string(ts) + " 2J=" + std::to_string(tJ) +
                          " 2x=" + std::to_string(tx),
                      got, want, 1e-13);
                norm += got * got;
                mscheme++;
              }
              close(std::string(v.name) + ": row norm", norm, 1.0, 1e-13);
            }
          }
  const long recChecks = checks - before;
  const double recWorst = worst;

  std::printf("i^l sign\n");
  {
    Convention c;
    c.orbital_phase = OrbitalPhase::IPowL;
    const double want[6] = {1, 1, -1, -1, 1, 1};
    for (int l = 0; l < 6; l++) close("i^l sign l=" + std::to_string(l), recoupling_coefficient(0, 0, l, 0, 2 * l, c, 0), want[l]);
    c.orbital_phase = OrbitalPhase::IPowMinusL;
    for (int l = 0; l < 6; l++) close("i^-l sign l=" + std::to_string(l), recoupling_coefficient(0, 0, l, 0, 2 * l, c, 0), want[l]);
  }

  std::printf("  exact values:            %ld checks, largest difference %.2g\n", exactChecks, exactWorst);
  std::printf("  orthogonality:           %ld checks, largest difference %.2g\n", orthChecks, orthWorst);
  std::printf("  6j vs four-3j sum:       %ld checks, largest difference %.2g\n", sixChecks, sixWorst);
  std::printf("  recoupling vs m-scheme:  %ld coefficients (%zu conventions), %ld checks, largest difference %.2g\n",
              mscheme, variants.size(), recChecks, recWorst);
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
