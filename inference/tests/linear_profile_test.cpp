// inference_linear_profile: the closed-form profiles of the objective (azure2/inference/LinearProfile.h), the first
// unit test of SolveThmProfile.  The line s m + a_0 + a_1 E + a_2 E^2 recovered by hand with 0, 1, 2, 3 background
// terms; chi^2 of noisy data equal to the minimum of a scan; a singular normal matrix (with and without the model
// column); s <= 0; the profiled norm's closed form against a scan, its degenerate branch by hand, and against the
// profile with no background; the reported norm and background with their covariance; ThmProfileDerivative against
// central differences of the re-solved profile.
#include "azure2/inference/LinearProfile.h"

#include <cmath>

#include "check.h"

using namespace azure2::inference;

namespace {
const int kPoints = 9;
// A model m_i(p) = p_0 exp(-p_1 E_i) (b) on E_i = 0.5 .. 2.5 MeV, and its derivative columns d m_i/d p (n x 2).
std::vector<double> Energies() {
  std::vector<double> E;
  for (int i = 0; i < kPoints; i++) E.push_back(0.5 + 0.25 * i);
  return E;
}
std::vector<double> Model(const double p[2], const std::vector<double> &E) {
  std::vector<double> m;
  for (double x : E) m.push_back(p[0] * std::exp(-p[1] * x));
  return m;
}
std::vector<double> ModelColumns(const double p[2], const std::vector<double> &E) {
  std::vector<double> Jm;
  for (double x : E) {
    Jm.push_back(std::exp(-p[1] * x));
    Jm.push_back(-x * p[0] * std::exp(-p[1] * x));
  }
  return Jm;
}
// A fixed pattern of deviations, so that the data are not on the curve (chi^2 > 0).
double Noise(int i) { return 0.03 * std::sin(1.7 * i + 0.3); }
bool Close(double x, double y, double tol) { return std::fabs(x - y) <= tol; }
double Chi2At(const std::vector<double> &m, const std::vector<double> &d, const std::vector<double> &e,
              const std::vector<double> &E, double s, const double a[3], int terms) {
  double chi2 = 0.0;
  for (size_t i = 0; i < m.size(); i++) {
    if (e[i] == 0.0) continue;
    double f = s * m[i];
    for (int k = 0; k < terms; k++) f += a[k] * std::pow(E[i], k);
    chi2 += (f - d[i]) * (f - d[i]) / (e[i] * e[i]);
  }
  return chi2;
}
}  // namespace

int main() {
  const std::vector<double> E = Energies();
  const double p[2] = {2.0, 0.8};
  const std::vector<double> m = Model(p, E);
  std::vector<double> e(kPoints, 0.05);
  e[4] = 0.0;  // a point without weight

  // 1. The line by hand: d = s m + a_0 + a_1 E + a_2 E^2 exactly, s = 0.8 (n = 1.25), a = (0.3, -0.1, 0.02).
  const double s0 = 0.8, a0[3] = {0.3, -0.1, 0.02};
  for (int terms = 0; terms <= 3; terms++) {
    std::vector<double> d;
    for (int i = 0; i < kPoints; i++) {
      double f = s0 * m[i];
      for (int k = 0; k < terms && k < 3; k++) f += a0[k] * std::pow(E[i], k);
      d.push_back(i == 4 ? 99.0 : f);  // the weightless point is ignored
    }
    const ThmProfile prof = SolveThmProfile(m, d, e, E, terms < 3 ? terms : 3);
    bool ok = prof.scaleProfiled && Close(prof.s, s0, 1e-9) && Close(prof.Norm(), 1.25, 1e-9) && prof.chi2 < 1e-15 &&
              prof.points == kPoints - 1 && prof.status == "profiled";
    for (int k = 0; k < terms; k++) ok = ok && Close(prof.a[k], a0[k], 1e-8);
    check("the line by hand, " + std::to_string(terms) + " background terms: s = 0.8, a, chi2 = 0", ok);
  }

  // 2. Noisy data, a linear background: chi^2 = the minimum of a scan in (s, a_0, a_1) around the solution, and the
  //    scan finds no lower value.
  std::vector<double> d;
  for (int i = 0; i < kPoints; i++) d.push_back(s0 * m[i] + a0[0] + a0[1] * E[i] + Noise(i));
  const ThmProfile lin = SolveThmProfile(m, d, e, E, 2);
  const double aLin[3] = {lin.a[0], lin.a[1], 0.0};
  check("chi2 = the chi2 of the line at the solution (1e-12)", Close(lin.chi2, Chi2At(m, d, e, E, lin.s, aLin, 2), 1e-12));
  double scanMin = 1e300;
  for (int i = -10; i <= 10; i++)
    for (int j = -10; j <= 10; j++)
      for (int k = -10; k <= 10; k++) {
        const double a[3] = {lin.a[0] + 1e-3 * j, lin.a[1] + 1e-3 * k, 0.0};
        scanMin = std::min(scanMin, Chi2At(m, d, e, E, lin.s + 1e-3 * i, a, 2));
      }
  check("no point of a scan in (s, a0, a1) lies below chi2", scanMin >= lin.chi2 - 1e-12);
  check("the scan's minimum is the profile's chi2 (the centre of the scan)", Close(scanMin, lin.chi2, 1e-12));

  // 3. A singular normal matrix.  (a) All points at one energy: the columns 1 and E are collinear, with and without
  //    the model: s = 1, no background.  (b) A constant model with a constant background: the model column is
  //    collinear with 1, the background alone is regular: s = 1, a_0 profiled on d - m.
  {
    const std::vector<double> E1(kPoints, 1.0);
    const ThmProfile sing = SolveThmProfile(m, d, e, E1, 2);
    check("one energy, 2 terms: s = 1, nothing profiled",
          !sing.scaleProfiled && !sing.backgroundProfiled && sing.s == 1.0 && sing.a[0] == 0.0 &&
              sing.status.find("collinear") != std::string::npos);
    const std::vector<double> mc(kPoints, 0.4);
    std::vector<double> dc;
    for (int i = 0; i < kPoints; i++) dc.push_back(0.4 + 0.25 + Noise(i));
    const ThmProfile bg = SolveThmProfile(mc, dc, e, E, 1);
    double mean = 0.0;
    int n = 0;
    for (int i = 0; i < kPoints; i++)
      if (e[i] != 0.0) mean += dc[i] - mc[i], n++;
    mean /= n;
    check("constant model, 1 term: s = 1, a0 = the mean of d - m (1e-12)",
          !bg.scaleProfiled && bg.backgroundProfiled && bg.s == 1.0 && Close(bg.a[0], mean, 1e-12) &&
              bg.status.find("background profiled") != std::string::npos);
  }

  // 4. s <= 0: data falling where the model rises.  With no background (S_md <= 0) and with one (c_0 <= 0).
  {
    std::vector<double> neg;
    for (int i = 0; i < kPoints; i++) neg.push_back(-m[i] + Noise(i));
    const ThmProfile p0 = SolveThmProfile(m, neg, e, E, 0);
    check("s <= 0, no background: s = 1, not profiled", !p0.scaleProfiled && p0.s == 1.0 && p0.cov.empty());
    const ThmProfile p1 = SolveThmProfile(m, neg, e, E, 1);
    check("s <= 0, one term: s = 1, the background profiled alone", !p1.scaleProfiled && p1.s == 1.0 && p1.backgroundProfiled && p1.cov.size() == 1);
    const double a1[3] = {p1.a[0], 0.0, 0.0};
    check("  its chi2 is that of 1 m + a0 (1e-12)", Close(p1.chi2, Chi2At(m, neg, e, E, 1.0, a1, 1), 1e-12));
  }

  // 5. The profiled norm (T1a): n* = S_mm/S_md against a scan of n; equal to 1/s of the profile without background;
  //    the degenerate branch by hand.
  {
    std::vector<double> dn;
    for (int i = 0; i < kPoints; i++) dn.push_back(m[i] / 1.25 + Noise(i));
    double norm = 0.0;
    const double chi2 = ProfiledNormChi2(m, dn, e, norm);
    auto chi2AtNorm = [&](double n) {
      double c = 0.0;
      for (int i = 0; i < kPoints; i++)
        if (e[i] != 0.0) c += std::pow((m[i] - n * dn[i]) / (n * e[i]), 2.0);
      return c;
    };
    double best = 1e300, bestN = 0.0;
    for (int k = -2000; k <= 2000; k++) {
      const double n = norm + 1e-6 * k;
      if (chi2AtNorm(n) < best) best = chi2AtNorm(n), bestN = n;
    }
    check("profiled norm: chi2 = the scan's minimum (1e-10), n* its place (1e-6)",
          Close(chi2, best, 1e-10) && Close(norm, bestN, 1e-6) && Close(chi2, chi2AtNorm(norm), 1e-10));
    const ThmProfile p0 = SolveThmProfile(m, dn, e, E, 0);
    check("profiled norm = 1/s of the profile with no background, chi2 equal (1e-12)",
          Close(norm, p0.Norm(), 1e-12) && Close(chi2, p0.chi2, 1e-10));
    // Degenerate: m = (1, 2, 5), d = (-3, 1, 7), e = (0.5, 1, 0): S_md = -3/0.25 + 2 = -10 <= 0; the third point has
    // no weight.
    const std::vector<double> m2 = {1.0, 2.0, 5.0}, d2 = {-3.0, 1.0, 7.0}, e2 = {0.5, 1.0, 0.0};
    double n2 = 0.0;
    const double c2 = ProfiledNormChi2(m2, d2, e2, n2);
    check("degenerate: n = 1, chi2 = ((1+3)/0.5)^2 + ((2-1)/1)^2 = 65", n2 == 1.0 && c2 == 65.0);
    double nEmpty = 0.0;
    check("no points: degenerate, n = 1, chi2 = 0", ProfiledNormChi2({}, {}, {}, nEmpty) == 0.0 && nEmpty == 1.0);
  }

  // 6. Reported: (n, b_k) = (1/s, a_k/s); the covariance of n is cov(s)/s^4.
  {
    double value[4], cov[16];
    lin.Reported(value, cov);
    check("reported n = 1/s, b = a/s, b2 = 0 (1e-14)",
          Close(value[0], 1.0 / lin.s, 1e-14) && Close(value[1], lin.a[0] / lin.s, 1e-14) &&
              Close(value[2], lin.a[1] / lin.s, 1e-14) && value[3] == 0.0);
    check("reported var(n) = var(s)/s^4 (1e-12 relative)",
          Close(cov[0], lin.cov[0] / std::pow(lin.s, 4), 1e-12 * cov[0]) && cov[15] == 0.0);
  }

  // 7. ThmProfileDerivative against central differences: J = dr_i/dp with the profile re-solved at every p, and
  //    G = d f_i/dp / s(p0), f the curve on the data scale; with a quadratic background and with none.
  for (int terms : {0, 2}) {
    const ThmProfile prof = SolveThmProfile(m, d, e, E, terms);
    std::vector<double> Jr, G;
    ThmProfileDerivative(prof, m, d, e, E, ModelColumns(p, E), 2, Jr, G);
    bool okJ = true, okG = true;
    for (int c = 0; c < 2; c++) {
      const double h = 1e-6 * p[c];
      double up[2] = {p[0], p[1]}, down[2] = {p[0], p[1]};
      up[c] += h;
      down[c] -= h;
      const std::vector<double> mu = Model(up, E), md = Model(down, E);
      const ThmProfile pu = SolveThmProfile(mu, d, e, E, terms), pd = SolveThmProfile(md, d, e, E, terms);
      for (int i = 0; i < kPoints; i++) {
        const double fdJ = (pu.Residual(mu[i], d[i], e[i], E[i]) - pd.Residual(md[i], d[i], e[i], E[i])) / (2.0 * h);
        if (!Close(Jr[i * 2 + c], fdJ, 1e-6 * (1.0 + std::fabs(fdJ)))) okJ = false;
        const double fu = pu.s * mu[i] + pu.s * pu.Background(E[i]), fd = pd.s * md[i] + pd.s * pd.Background(E[i]);
        const double fdG = (fu - fd) / (2.0 * h) / prof.s;
        if (!Close(G[i * 2 + c], fdG, 1e-6 * (1.0 + std::fabs(fdG)))) okG = false;
      }
    }
    check("ThmProfileDerivative J = central differences, " + std::to_string(terms) + " terms (1e-6)", okJ);
    check("ThmProfileDerivative G = central differences, " + std::to_string(terms) + " terms (1e-6)", okG);
    check("  the weightless point has J = 0", Jr[4 * 2] == 0.0 && Jr[4 * 2 + 1] == 0.0);
  }
  return result();
}
