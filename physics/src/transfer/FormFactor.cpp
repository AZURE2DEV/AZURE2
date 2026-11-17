#include "azure2/physics/transfer/FormFactor.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace physics {
namespace transfer {
namespace form_factor {

namespace {
const double kPi = std::acos(-1.0);
}  // namespace


double Well::H2M() const { return kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * mu); }

double Well::Potential(double V0, double r) const {
  double v = -V0 / (1.0 + std::exp((r - R) / diff));
  if (zze2 != 0.0) v += r < rc ? zze2 * (3.0 - (r / rc) * (r / rc)) / (2.0 * rc) : zze2 / r;
  return v;
}

double SphericalJ(int l, double x) {
  if (x <= 0.0) return l == 0 ? 1.0 : 0.0;
  if (x < std::max(2.5, l + 1.0)) {
    // x^l/(2l+1)!! sum_k (-x^2/2)^k / (k! (2l+3)(2l+5)...(2l+2k+1)).
    double front = 1.0;
    for (int k = 1; k <= l; k++) front *= x / (2.0 * k + 1.0);
    double term = 1.0, sum = 1.0;
    const double y = -0.5 * x * x;
    for (int k = 1; k < 200; k++) {
      term *= y / (k * (2.0 * l + 2.0 * k + 1.0));
      sum += term;
      if (std::fabs(term) < 1.0e-17 * std::fabs(sum)) break;
    }
    return front * sum;
  }
  const double s = std::sin(x), c = std::cos(x);
  double j0 = s / x;
  if (l == 0) return j0;
  double j1 = s / (x * x) - c / x;
  for (int n = 1; n < l; n++) {
    double j2 = (2.0 * n + 1.0) / x * j1 - j0;
    j0 = j1;
    j1 = j2;
  }
  return j1;
}

double SurfaceM(int l, double rho, double B) {
  const double jl = SphericalJ(l, rho);
  // j_l'(x) = (l/x) j_l - j_{l+1}; at x = 0 only l = 1 has j_l' = 1/3.
  const double djl = rho > 0.0 ? l * jl / rho - SphericalJ(l + 1, rho) : (l == 1 ? 1.0 / 3.0 : 0.0);
  return (B - 1.0) * jl - rho * djl;
}

namespace {

// Simpson weight of grid point i of n steps (n even), times h.
inline double SimpsonWeight(int i, int n, double h) {
  if (i == 0 || i == n) return h / 3.0;
  return (i % 2 ? 4.0 : 2.0) * h / 3.0;
}

// The Pruefer angle pi n + acot(a u'/u) in (n pi, (n+1) pi].
double Theta(int nodes, double u, double av) {
  double f = std::atan2(u, av);
  if (f <= 0.0) f += kPi;
  return kPi * nodes + f;
}

}  // namespace

namespace {

// The regular solution at depth V0 on r_i = i a/steps (u, unnormalized),
// with a u'(a) and the node count in (0, a].
void Integrate(const Well &w, int l, double e, double V0, std::vector<double> &uu, double &av, int &nodes) {
  const int n = w.steps;
  const double h = w.a / n, h2m = w.H2M();
  // q(r) = l(l+1)/r^2 + (V(r) - e)/(hbar^2/2mu) on the half grid r = k h/2, k >= 2.
  std::vector<double> q(2 * n + 1, 0.0);
  for (int k = 2; k <= 2 * n; k++) {
    const double r = 0.5 * k * h;
    q[k] = l * (l + 1.0) / (r * r) + (w.Potential(V0, r) - e) / h2m;
  }
  // The regular solution r^(l+1) (1 + c1 r^2 + c2 r^4) at r = h (constant q near 0).
  const double k0 = (w.Potential(V0, 0.0) - e) / h2m;
  const double c1 = k0 / (4.0 * l + 6.0), c2 = k0 * c1 / (8.0 * l + 20.0);
  double u = std::pow(h, l + 1) * (1.0 + c1 * h * h + c2 * h * h * h * h);
  double v = (l + 1.0) * std::pow(h, l) + c1 * (l + 3.0) * std::pow(h, l + 2) + c2 * (l + 5.0) * std::pow(h, l + 4);
  uu.assign(n + 1, 0.0);
  uu[1] = u;
  for (int i = 1; i < n; i++) {
    const double qa = q[2 * i], qm = q[2 * i + 1], qb = q[2 * i + 2];
    const double k1u = v, k1v = qa * u;
    const double k2u = v + 0.5 * h * k1v, k2v = qm * (u + 0.5 * h * k1u);
    const double k3u = v + 0.5 * h * k2v, k3v = qm * (u + 0.5 * h * k2u);
    const double k4u = v + h * k3v, k4v = qb * (u + h * k3u);
    u += h / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
    v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    uu[i + 1] = u;
    if (std::fabs(u) > 1.0e150) {  // a strongly repulsive well: keep the scale finite
      for (int k = 1; k <= i + 1; k++) uu[k] *= 1.0e-150;
      u *= 1.0e-150;
      v *= 1.0e-150;
    }
  }
  nodes = 0;
  for (int i = 1; i < n; i++)
    if ((uu[i] > 0.0) != (uu[i + 1] > 0.0) && uu[i + 1] != 0.0) nodes++;
  av = w.a * v;
}

}  // namespace

void SolveAtDepth(const Well &w, int l, double e, double V0, double &logDerivative, int &nodes,
                  std::vector<double> *uOut) {
  std::vector<double> local;
  std::vector<double> &u = uOut ? *uOut : local;
  double av;
  Integrate(w, l, e, V0, u, av, nodes);
  logDerivative = av / u.back();
}

LevelFunction SolveLevel(const Well &w, int l, double e, double B, int nodes) {
  LevelFunction f;
  f.l = l;
  f.e = e;
  f.B = B;
  f.nodes = nodes;
  if (!(w.mu > 0.0) || !(w.a > 0.0) || w.steps < 4 || w.steps % 2 || nodes < 0 || l < 0 || !std::isfinite(B)) {
    f.why = "invalid well or level";
    return f;
  }
  const double target = kPi * nodes + (0.5 * kPi - std::atan(B));
  std::vector<double> u;
  auto F = [&](double V0) {
    double av;
    int n;
    Integrate(w, l, e, V0, u, av, n);
    return Theta(n, u.back(), av) - target;
  };
  // theta increases with the depth: bracket, then the Illinois method.
  const double vLimit = 1.0e4;
  double lo = 0.0, hi = 0.0, flo = F(0.0), fhi = flo;
  if (flo == 0.0) {
    hi = lo;
  } else if (flo < 0.0) {
    double step = 10.0;
    while (fhi < 0.0) {
      lo = hi;
      flo = fhi;
      hi = lo + step;
      step *= 2.0;
      if (hi > vLimit) {
        std::ostringstream s;
        s << "no Woods-Saxon depth below " << vLimit << " MeV gives " << nodes << " node(s) and a u'/u = " << B;
        f.why = s.str();
        return f;
      }
      fhi = F(hi);
    }
  } else {
    double step = 10.0;
    while (flo > 0.0) {
      hi = lo;
      fhi = flo;
      lo = hi - step;
      step *= 2.0;
      if (lo < -vLimit) {
        std::ostringstream s;
        s << "no Woods-Saxon depth above " << -vLimit << " MeV gives " << nodes << " node(s) and a u'/u = " << B;
        f.why = s.str();
        return f;
      }
      flo = F(lo);
    }
  }
  double x = lo;
  if (flo != 0.0 && fhi != 0.0) {
    int side = 0;
    for (int it = 0; it < 300; it++) {
      x = (lo * fhi - hi * flo) / (fhi - flo);
      if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
      const double fx = F(x);
      if (fx == 0.0) break;
      if ((fx < 0.0) == (flo < 0.0)) {
        lo = x;
        flo = fx;
        if (side == -1) fhi *= 0.5;
        side = -1;
      } else {
        hi = x;
        fhi = fx;
        if (side == 1) flo *= 0.5;
        side = 1;
      }
      if (hi - lo < 1.0e-13 * (1.0 + std::fabs(x))) break;
    }
  } else if (fhi == 0.0) {
    x = hi;
  }
  double av;
  int n;
  Integrate(w, l, e, x, u, av, n);
  const double L = av / u.back();
  if (n != nodes || !(std::fabs(L - B) <= 1.0e-6 * (1.0 + std::fabs(B)))) {
    std::ostringstream s;
    s << "the depth search did not converge (V0 = " << x << " MeV: " << n << " node(s), a u'/u = " << L
      << ", wanted " << nodes << " and " << B << ")";
    f.why = s.str();
    return f;
  }
  const int steps = w.steps;
  const double h = w.a / steps;
  double norm = 0.0;
  for (int i = 0; i <= steps; i++) norm += SimpsonWeight(i, steps, h) * u[i] * u[i];
  const double scale = (u[steps] < 0.0 ? -1.0 : 1.0) / std::sqrt(norm);
  for (double &x2 : u) x2 *= scale;
  f.u.swap(u);
  f.ua = f.u[steps];
  f.depth = x;
  f.ok = true;
  return f;
}

double Overlap(const Well &w, const LevelFunction &f, double p) {
  const int n = w.steps;
  const double h = w.a / n;
  double s = 0.0;
  for (int i = 1; i <= n; i++) {
    const double r = i * h;
    s += SimpsonWeight(i, n, h) * f.u[i] * r * SphericalJ(f.l, p * r);
  }
  return s;
}

double PotentialOverlap(const Well &w, const LevelFunction &f, double p) {
  const int n = w.steps;
  const double h = w.a / n;
  double s = 0.0;
  for (int i = 1; i <= n; i++) {
    const double r = i * h;
    s += SimpsonWeight(i, n, h) * f.u[i] * w.Potential(f.depth, r) * r * SphericalJ(f.l, p * r);
  }
  return s;
}

double VolumeTerm(const Well &w, const LevelFunction &f, double eps) {
  const double h2m = w.H2M();
  const double p = std::sqrt(eps / h2m);
  return (f.e - eps) / h2m * Overlap(w, f, p) / f.ua;
}

int WildermuthNodes(int l, int coreCount) {
  if (l < 0 || coreCount < 0) return -1;
  int N = 0;
  while ((N + 1) * (N + 2) * (N + 3) / 3 <= coreCount) N++;
  while (N < l || (N - l) % 2) N++;
  return (N - l) / 2;
}

std::vector<double> Chebyshev::Nodes(double lo, double hi, int n) {
  std::vector<double> x(n);
  for (int k = 0; k < n; k++) x[k] = 0.5 * (lo + hi) + 0.5 * (hi - lo) * std::cos(kPi * (k + 0.5) / n);
  return x;
}

void Chebyshev::Fit(double l0, double h0, const std::vector<double> &values) {
  lo = l0;
  hi = h0;
  const int n = (int)values.size();
  c.assign(n, 0.0);
  for (int j = 0; j < n; j++) {
    double s = 0.0;
    for (int k = 0; k < n; k++) s += values[k] * std::cos(kPi * j * (k + 0.5) / n);
    c[j] = 2.0 * s / n;
  }
  if (n) c[0] *= 0.5;
}

double Chebyshev::operator()(double x) const {
  const double t = (2.0 * x - lo - hi) / (hi - lo);
  double b1 = 0.0, b2 = 0.0;
  for (int j = (int)c.size() - 1; j >= 1; j--) {
    const double b0 = 2.0 * t * b1 - b2 + c[j];
    b2 = b1;
    b1 = b0;
  }
  return t * b1 - b2 + (c.empty() ? 0.0 : c[0]);
}

}  // namespace form_factor
}  // namespace transfer
}  // namespace physics
}  // namespace azure2
