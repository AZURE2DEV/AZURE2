#include "ThmAngular.h"

#include <gsl/gsl_sf_legendre.h>

#include <algorithm>
#include <cmath>

#include "AngCoeff.h"

namespace {

// Zb(l1 J1 l2 J2; s L) as CNuc::CalcAngularDists writes it.
double Zbar(int l1, double J1, int l2, double J2, double s, int L) {
  return std::sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * J1 + 1.0) * (2.0 * J2 + 1.0)) *
         AngCoeff::ClebGord(l1, l2, L, 0.0, 0.0, 0.0) * AngCoeff::Racah(l1, J1, l2, J2, s, L);
}

struct Term {
  int i, j, L;
  double c;  // includes the 2 of i != j
};

struct Table {
  std::vector<ThmPartialWave> key;
  std::vector<Term> terms;
  int maxL = 0;
  bool built = false;
};

void Build(const std::vector<ThmPartialWave> &waves, Table &t) {
  t.key = waves;
  t.terms.clear();
  t.maxL = 0;
  const int n = (int)waves.size();
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      const ThmPartialWave &a = waves[i], &b = waves[j];
      if (a.s != b.s || a.sp != b.sp) continue;  // channel spins: incoherent
      int lo = std::max(std::max(std::abs(a.l - b.l), std::abs(a.lp - b.lp)), (int)std::lround(std::fabs(a.J - b.J)));
      int hi = std::min(std::min(a.l + b.l, a.lp + b.lp), (int)std::lround(a.J + b.J));
      for (int L = lo; L <= hi; L++) {
        if ((a.l + b.l + L) % 2 || (a.lp + b.lp + L) % 2) continue;
        double sign = ((int)std::lround(a.sp - a.s)) % 2 ? -1.0 : 1.0;
        double c = sign / 4.0 * Zbar(a.l, a.J, b.l, b.J, a.s, L) * Zbar(a.lp, a.J, b.lp, b.J, a.sp, L);
        if (i != j) c *= 2.0;
        if (std::fabs(c) < 1.0e-14) continue;
        t.terms.push_back(Term{i, j, L, c});
        t.maxL = std::max(t.maxL, L);
      }
    }
  t.built = true;
}

}  // namespace

void ThmLegendreCoefficients(const std::vector<ThmPartialWave> &waves, const std::vector<complex> &x,
                             std::vector<double> &b) {
  thread_local Table table;
  if (!table.built || !(table.key == waves)) Build(waves, table);
  b.assign(table.maxL + 1, 0.0);
  for (const Term &t : table.terms) b[t.L] += t.c * std::real(x[t.i] * std::conj(x[t.j]));
}

double ThmAngleWindow::Mean(const std::vector<double> &b) const {
  double sum = 0.0;
  const int n = std::min((int)b.size(), (int)meanP.size());
  for (int L = 0; L < n; L++) sum += b[L] * meanP[L];
  return sum / M_PI;  // exact pi (Constants.h pi is 3.14159265): 4 pi <.> over 0-180 == sigma
}

void BuildThmAngleWindow(double thetaMin, double thetaMax, ThmAngleWindow &out) {
  out.thetaMin = thetaMin;
  out.thetaMax = thetaMax;
  out.meanP.assign(ThmAngleWindow::kMaxL + 1, 0.0);
  std::vector<double> p(ThmAngleWindow::kMaxL + 1);
  // cos of the ends; the exact values at 0, 90 and 180 degrees.
  auto cosDeg = [](double t) { return t == 0.0 ? 1.0 : t == 90.0 ? 0.0 : t == 180.0 ? -1.0 : std::cos(t * M_PI / 180.0); };
  const double xHi = cosDeg(thetaMin), xLo = cosDeg(thetaMax);
  if (thetaMin == thetaMax) {
    gsl_sf_legendre_Pl_array(ThmAngleWindow::kMaxL, xHi, p.data());
    out.meanP = p;
    return;
  }
  // Gauss-Legendre nodes on [xLo, xHi] by Newton's method on P_n (GSL's
  // fixed tables are good to ~1e-12 only at this order).
  const int nodes = 48;  // exact for polynomials of degree <= 95
  const double mid = 0.5 * (xHi + xLo), half = 0.5 * (xHi - xLo);
  for (int k = 0; k < (nodes + 1) / 2; k++) {
    double z = std::cos(M_PI * (k + 0.75) / (nodes + 0.5)), dp = 0.0;
    for (int it = 0; it < 100; it++) {
      double p0 = 1.0, p1 = 0.0;
      for (int j = 1; j <= nodes; j++) {
        double p2 = p1;
        p1 = p0;
        p0 = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p2) / j;
      }
      dp = nodes * (z * p0 - p1) / (z * z - 1.0);
      double dz = p0 / dp;
      z -= dz;
      if (std::fabs(dz) < 1.0e-16) break;
    }
    const double wk = 2.0 / ((1.0 - z * z) * dp * dp) * half;
    for (double xk : {mid + half * z, mid - half * z}) {
      gsl_sf_legendre_Pl_array(ThmAngleWindow::kMaxL, xk, p.data());
      for (int L = 0; L <= ThmAngleWindow::kMaxL; L++) out.meanP[L] += wk * p[L];
    }
  }
  for (int L = 0; L <= ThmAngleWindow::kMaxL; L++) out.meanP[L] /= (xHi - xLo);
  out.meanP[0] = 1.0;
}
