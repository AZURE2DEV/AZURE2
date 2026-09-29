/*!
 * Reference check of the fixed-angle THM observable (ThmAngular.h,
 * src/ThmAngular.cpp): the spin-summed angular distribution of the exit pair
 * for an entrance with m_l = 0 along p_xA,
 *
 *   sum_{nu nu'} |F(theta)|^2 = (1/pi) sum_L b_L P_L(cos theta),
 *
 * computed by the engine in the Blatt-Biedenharn form (GSL 3j/6j symbols) and
 * averaged over a window by Gauss-Legendre in cos theta, against
 * thm_fixed_angle_reference.py, which does the literal M-sum
 *   F_{nu nu'} = sum sqrt(2l+1) <s nu l 0|J nu><s' nu' l' m'|J nu> x Y_l'^m'(theta, 0)
 * with sympy's exact Clebsch-Gordan coefficients and mpmath's spherical
 * harmonics and quadrature.  The amplitudes x come from a toy two-level
 * R-matrix assembled as THMMatrixFunc does (two entrance l per channel spin).
 *
 *  (d) b_L, dsigma/dOmega at 0, 37, 90, 143, 180 deg and window averages
 *      against the reference, 1e-12;
 *  (b) the 0-180 window is sum (2J+1)|x|^2 / 4 pi (the angle-integrated
 *      HOES cross section over 4 pi), 1e-13;
 *  (c) theta = 0 is the m_l = 0 amplitude along the axis (nu' = nu), and a
 *      0-0.5 deg window approaches it; for a' + a exit (s' = 0) it is
 *      (1/4 pi) sum_s |sum_{J,l} sqrt((2l+1)(2J+1)) <s 0 l 0|J 0> x|^2 --
 *      coherent in l and J with Clebsch-Gordan weights, which is not the
 *      entranceL=coherent recipe (equal weights, J incoherent);
 *  (e) an identical-boson exit (7Li + p -> a + a like): no odd L, symmetric
 *      about 90 deg; the asymmetric case is not.
 *
 * Run:  tests/reference/thm_fixed_angle_test      (ctest: thm_fixed_angle)
 *       tests/reference/thm_fixed_angle_test -v   (print every value)
 */
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "AngCoeff.h"
#include "Config.h"
#include "ThmAngular.h"

// Each consumer of the engine defines this itself (see thm_coulomb_term_test).
Config *g_config = nullptr;

namespace {

bool verbose = false;
int failures = 0;

void Check(const std::string &what, double got, double want, double tol, double scale = 0.0) {
  double den = std::max(std::fabs(want), scale > 0.0 ? scale : 1.0e-300);
  double rel = std::fabs(got - want) / den;
  bool ok = rel <= tol;
  if (!ok) failures++;
  if (!ok || verbose)
    std::printf("  %s  %-40s got %.15e want %.15e rel %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got,
                want, rel, tol);
}

struct Wave {
  ThmPartialWave w;
  double re, im;
};

struct Case {
  const char *name;
  std::vector<Wave> waves;
  double sigma;
  std::vector<double> b;
  double ds[5];  // 0, 37, 90, 143, 180 deg
  double axis;
  double win[5];  // 50-70, 0-180, 0-0.5, 20-160, 110-130
};

ThmPartialWave W(double J, double s, int l, double sp, int lp) {
  ThmPartialWave w;
  w.J = J;
  w.s = s;
  w.l = l;
  w.sp = sp;
  w.lp = lp;
  return w;
}

// python3 tests/reference/thm_fixed_angle_reference.py
std::vector<Case> Cases() {
  std::vector<Case> c(2);
  c[0].name = "li7";
  c[0].waves = {{W(2, 1, 1, 0, 2), -2.65248456898447726e-01, 2.53372452773144241e-01},
                {W(2, 1, 3, 0, 2), -1.54436106945734404e-02, 1.71205812999355236e-02},
                {W(2, 2, 1, 0, 2), -2.10588932890686864e-01, 1.93468194501487278e-01},
                {W(2, 2, 3, 0, 2), 2.23473578776968267e-03, -4.24176515420260272e-03},
                {W(0, 1, 1, 0, 0), 1.74291593540974593e-01, 3.10487376964923123e-01}};
  c[0].sigma = 1.21121292424355143e+00;
  c[0].b = {0.30280323106088791, 0.0, 0.015918916479243384, 0.0, -0.046505065609345961};
  double ds0[5] = {0.0866493883666719, 0.10220810349654398, 0.088300554465831447, 0.10220810349654398,
                   0.0866493883666719};
  double win0[5] = {0.099554431645706023, 0.0963852620150753, 0.086651917065770935, 0.09677746104687714,
                    0.099554431645706023};
  std::memcpy(c[0].ds, ds0, sizeof ds0);
  std::memcpy(c[0].win, win0, sizeof win0);
  c[0].axis = 0.086649388366671908;

  c[1].name = "asym";
  c[1].waves = {{W(0.5, 0.5, 0, 0.5, 1), -2.47186883930311518e-01, -2.88088934765062032e-01},
                {W(0.5, 0.5, 1, 0.5, 0), -2.24952906706327477e-01, -1.91914163458938719e-01},
                {W(1.5, 0.5, 1, 0.5, 2), 2.40898267811578137e-01, -1.25135991325982743e-01}};
  c[1].sigma = 7.57826859625051807e-01;
  c[1].b = {0.18945671490626293, 0.063900559511654296, 0.013340198818149568, 0.0, 0.0};
  double ds1[5] = {0.084892442351276984, 0.078489747266618613, 0.058182786774830269, 0.04600096748825784,
                   0.044212082700806543};
  double win1[5] = {0.069790739107410211, 0.060305945358567431, 0.084891812578374835, 0.060057582984164955,
                    0.049759572317856066};
  std::memcpy(c[1].ds, ds1, sizeof ds1);
  std::memcpy(c[1].win, win1, sizeof win1);
  c[1].axis = 0.084892442351276986;
  return c;
}

double At(const std::vector<double> &b, double thmin, double thmax) {
  ThmAngleWindow w;
  BuildThmAngleWindow(thmin, thmax, w);
  return w.Mean(b);
}

}  // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; i++)
    if (!std::strcmp(argv[i], "-v")) verbose = true;

  for (const Case &c : Cases()) {
    std::vector<ThmPartialWave> waves;
    std::vector<complex> x;
    for (const Wave &w : c.waves) {
      waves.push_back(w.w);
      x.push_back(complex(w.re, w.im));
    }
    std::vector<double> b;
    ThmLegendreCoefficients(waves, x, b);
    std::string n = c.name;

    std::printf("(d) %s: Legendre coefficients, dsigma/dOmega and windows against the M-sum\n", c.name);
    const double scale = c.b[0];  // odd / vanishing b_L: absolute, in units of b_0
    for (size_t L = 0; L < c.b.size(); L++)
      Check(n + " b_" + std::to_string(L), L < b.size() ? b[L] : 0.0, c.b[L], 1.0e-12, scale);
    const double angles[5] = {0.0, 37.0, 90.0, 143.0, 180.0};
    for (int k = 0; k < 5; k++)
      Check(n + " dsigma/dOmega(" + std::to_string((int)angles[k]) + ")", At(b, angles[k], angles[k]), c.ds[k],
            1.0e-12);
    const double lo[5] = {50.0, 0.0, 0.0, 20.0, 110.0}, hi[5] = {70.0, 180.0, 0.5, 160.0, 130.0};
    for (int k = 0; k < 5; k++) {
      char label[64];
      std::snprintf(label, sizeof label, " window %g-%g", lo[k], hi[k]);
      Check(n + label, At(b, lo[k], hi[k]), c.win[k], 1.0e-12);
    }

    std::printf("(b) %s: the 0-180 window is the angle-integrated sigma / 4 pi\n", c.name);
    double sigma = 0.0;
    for (size_t i = 0; i < x.size(); i++) sigma += (2.0 * waves[i].J + 1.0) * std::norm(x[i]);
    Check(n + " sum (2J+1)|x|^2", sigma, c.sigma, 1.0e-14);
    Check(n + " 4 pi <.>_{0-180} / sigma", 4.0 * M_PI * At(b, 0.0, 180.0) / sigma, 1.0, 1.0e-13);
    Check(n + " 4 b_0 / sigma", 4.0 * b[0] / sigma, 1.0, 1.0e-13);

    std::printf("(c) %s: theta -> 0 is the m_l = 0 amplitude along the axis\n", c.name);
    Check(n + " dsigma/dOmega(0) vs axis sum", At(b, 0.0, 0.0), c.axis, 1.0e-12);
    double d1 = At(b, 0.0, 0.25) - At(b, 0.0, 0.0), d2 = At(b, 0.0, 0.5) - At(b, 0.0, 0.0);
    Check(n + " window 0-t converges as t^2", d2 / d1, 4.0, 1.0e-4);

    std::printf("(e) %s: symmetry about 90 deg\n", c.name);
    double f = At(b, 20.0, 40.0), bwd = At(b, 140.0, 160.0);
    if (n == "li7") {
      for (size_t L = 1; L < b.size(); L += 2) Check(n + " odd b_" + std::to_string(L), b[L], 0.0, 1.0e-15, b[0]);
      Check(n + " 20-40 == 140-160", f, bwd, 1.0e-14);
    } else {
      bool asym = std::fabs(f / bwd - 1.0) > 0.1;
      if (!asym) failures++;
      if (!asym || verbose)
        std::printf("  %s  asym 20-40 / 140-160 = %.6f (not symmetric)\n", asym ? "ok  " : "FAIL", f / bwd);
    }
  }

  std::printf("(c) a + a exit: theta = 0 is coherent in l and J with CG weights\n");
  {
    // (1/4 pi) sum_s |sum_{J,l} sqrt((2l+1)(2J+1)) <s 0 l 0|J 0> x|^2 for the
    // li7 case, from AngCoeff::ClebGord -- and it is not the entranceL=coherent
    // recipe sum_J (2J+1) sum_s |sum_l x|^2 / 4 pi.
    Case c = Cases()[0];
    double direct = 0.0, recipe = 0.0;
    for (double s : {1.0, 2.0}) {
      complex f(0.0, 0.0);
      for (const Wave &w : c.waves)
        if (w.w.s == s)
          f += std::sqrt((2.0 * w.w.l + 1.0) * (2.0 * w.w.J + 1.0)) *
               AngCoeff::ClebGord(s, w.w.l, w.w.J, 0.0, 0.0, 0.0) * complex(w.re, w.im);
      direct += std::norm(f) / (4.0 * M_PI);
      for (double J : {0.0, 2.0}) {
        complex g(0.0, 0.0);
        for (const Wave &w : c.waves)
          if (w.w.s == s && w.w.J == J) g += complex(w.re, w.im);
        recipe += (2.0 * J + 1.0) * std::norm(g) / (4.0 * M_PI);
      }
    }
    Check("li7 CG-weighted axis sum", direct, c.axis, 1.0e-12);
    bool differs = std::fabs(recipe / direct - 1.0) > 0.05;
    if (!differs) failures++;
    if (!differs || verbose)
      std::printf("  %s  entranceL=coherent recipe / theta=0 = %.6f (differ)\n", differs ? "ok  " : "FAIL",
                  recipe / direct);
  }

  if (failures) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all checks passed\n");
  return 0;
}
