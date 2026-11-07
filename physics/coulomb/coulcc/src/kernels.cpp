// Subprograms of coulcc.f (I. J. Thompson and A. R. Barnett, Comput. Phys.
// Commun. 36 (1985) 363) other than COULCC: CF1C, CF1R, CF2, F11, F20, CF1A,
// RCF, CLOGAM/CDIGAM/LOGAM and TIDY.  The arithmetic follows the Fortran
// statement by statement (order of operations, mixed-mode promotion,
// single-precision intrinsics), so that the results are bit-identical to
// gfortran's; the comments give the Fortran labels.  DESIGN.md section 3.
#include "kernels.h"

#include <algorithm>
#include <cmath>

namespace coulcc {
namespace detail {

namespace {

// ---------------------------------------------------------------------------
// Extended precision for F11 (REAL*16 in the Fortran).
#if defined(__SIZEOF_FLOAT128__) && !defined(COULCC_FORCE_DOUBLE_DOUBLE)
using wide = __float128;  // IEEE binary128, as gfortran's REAL*16
inline wide wabs(wide x) { return x == 0 ? wide(0) : (x < 0 ? -x : x); }
inline double to_double(wide x) { return double(x); }
#else
// Double-double (106-bit significand) where __float128 is not available.
// Not bit-identical to REAL*16; its accuracy is tested in test/.
struct wide {
  double hi = 0, lo = 0;
  wide() = default;
  wide(double h) : hi(h), lo(0) {}  // NOLINT
  wide(double h, double l) : hi(h), lo(l) {}
};
inline wide two_sum(double a, double b) {
  double s = a + b, bb = s - a;
  return {s, (a - (s - bb)) + (b - bb)};
}
inline wide quick_two_sum(double a, double b) {
  double s = a + b;
  return {s, b - (s - a)};
}
inline wide operator+(wide a, wide b) {
  wide s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = quick_two_sum(s.hi, s.lo);
  s.lo += t.lo;
  return quick_two_sum(s.hi, s.lo);
}
inline wide operator-(wide a) { return {-a.hi, -a.lo}; }
inline wide operator-(wide a, wide b) { return a + (-b); }
inline wide operator*(wide a, wide b) {
  double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  return quick_two_sum(p, e);
}
inline wide operator/(wide a, wide b) {
  double q1 = a.hi / b.hi;
  wide r = a - b * wide(q1);
  double q2 = r.hi / b.hi;
  r = r - b * wide(q2);
  double q3 = r.hi / b.hi;
  return quick_two_sum(q1, q2) + wide(q3);
}
inline bool operator<(wide a, wide b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
inline wide wabs(wide x) { return x.hi < 0 || (x.hi == 0 && x.lo < 0) ? -x : wide(std::fabs(x.hi), x.lo); }
inline double to_double(wide x) { return x.hi + x.lo; }
#endif

// The Bernoulli-number coefficients of Stirling's series and the constants
// that LOGAM(ACC) sets for CLOGAM and CDIGAM; built once for ACC = ACC8 and
// read-only afterwards.
struct LogGammaTable {
  double b[16] = {};  // B(1..NT)
  int nt = 0, nx0 = 6;
  double accur = 0, pi = 0, alpi = 0, hl2p = 0;
};

LogGammaTable make_loggamma_table(double acc) {
  static const double bn[16] = {0,         +1.0,      -1.0,          +1.0,       -1.0,
                                +5.0,      -691.0,    +7.0,          -3617.0,    43867.0,
                                -174611.0, 854513.0,  -236364091.0,  +8553103.0, -23749461029.0,
                                8615841276005.0};
  static const double bd[16] = {0, 6.0, 30.0, 42.0, 30.0, 66.0, 2730.0, 6.0, 510.0, 798.0, 330.0, 138.0,
                                2730.0, 6.0, 870.0, 14322.0};
  LogGammaTable t;
  t.nx0 = 6;
  const double x0 = t.nx0 + 1.0;
  t.pi = 4.0 * std::atan(1.0);
  t.alpi = std::log(t.pi);
  t.hl2p = std::log(2.0 * t.pi) * 0.5;
  t.accur = acc;
  const int nb = 15;
  double err = 0, f21 = 0;
  int k;
  for (k = 1; k <= nb; ++k) {
    f21 = k * 2 - 1.0;
    t.b[k] = bn[k] / (bd[k] * k * 2.0 * f21);
    err = std::fabs(t.b[k]) * k * 2.0 / std::pow(x0, f21);
    if (err < acc) break;  // 120 ... GO TO 130
  }
  if (k > nb) {
    t.nx0 = int(std::pow(err / acc, 1.0 / f21) * x0);
    k = nb;
  }
  t.nt = k;  // 130
  return t;
}

const LogGammaTable& loggamma_table() {
  static const LogGammaTable table = make_loggamma_table(kAcc8);  // CALL LOGAM(ACC8)
  return table;
}

}  // namespace

// ---------------------------------------------------------------------------
fcomplex tidy(fcomplex z, double acc) {
  double x = real(z), y = imag(z);
  const double az = (std::fabs(x) + std::fabs(y)) * acc * 5;
  if (std::fabs(x) < az) x = 0.0;
  if (std::fabs(y) < az) y = 0.0;
  return {x, y};
}

// ---------------------------------------------------------------------------
fcomplex clogam(fcomplex z, unsigned& warnings) {
  const LogGammaTable& t = loggamma_table();
  const double x = real(z), tt = imag(z);
  // MX = INT(REAL(ACCUR*100 - X)): REAL of a REAL*8 is single precision.
  const int mx = int(float(t.accur * 100 - x));
  if (std::fabs(std::fabs(x) - mx) + std::fabs(tt) < t.accur * 50) {  // 60: argument is a non-positive integer
    warnings |= kWarnLogGammaPole;
    return 0.0;
  }
  double f = std::fabs(tt);
  fcomplex v{x, f};
  if (x < 0.0) v = 1.0 - v;
  fcomplex h = 0.0;
  double c = real(v);
  const int n = t.nx0 - int(c);
  if (n >= 0) {
    h = v;
    const double d = imag(v);
    double a = std::atan2(d, c);
    for (int i = 1; i <= n; ++i) {  // 10
      c = c + 1.0;
      v = fcomplex(c, d);
      h = h * v;
      a = a + std::atan2(d, c);
    }
    h = fcomplex(0.5 * std::log(real(h) * real(h) + imag(h) * imag(h)), a);  // 20
    v = v + 1.0;
  }
  const fcomplex r = 1.0 / (v * v);  // 30
  fcomplex ser = t.b[t.nt];
  for (int j = 2; j <= t.nt; ++j) {
    const int k = t.nt + 1 - j;
    ser = t.b[k] + r * ser;  // 40
  }
  fcomplex result = t.hl2p + (v - 0.5) * log(v) - v + ser / v - h;
  if (x < 0.0) {
    const double a = int(x) - 1.0;
    const double c2 = t.pi * (x - a);
    const double d = t.pi * f;
    double e = 0.0;
    f = -2.0 * d;
    if (f > kFplmin) e = std::exp(f);
    f = std::sin(c2);
    e = d + 0.5 * std::log(e * (f * f) + 0.25 * ((1.0 - e) * (1.0 - e)));
    f = std::atan2(std::cos(c2) * std::tanh(d), f) - a * t.pi;
    result = fcomplex(t.alpi) - fcomplex(e, f) - result;
  }
  if (sign(1.0, tt) < -0.5) result = conj(result);  // 50
  return result;
}

fcomplex cdigam(fcomplex z, unsigned& warnings) {
  const LogGammaTable& t = loggamma_table();
  fcomplex u = z;
  const double x = real(u);
  const double a = std::fabs(x);
  if (std::fabs(imag(u)) + std::fabs(a + int(x)) < t.accur) {  // 110
    warnings |= kWarnLogGammaPole;
    return 0.0;
  }
  if (x < 0.0) u = -u;
  fcomplex v = u;
  fcomplex h = 0.0;
  const int n = t.nx0 - int(a);
  if (n >= 0) {
    h = 1.0 / v;
    for (int i = 1; i <= n; ++i) {  // 70
      v = v + 1.0;
      h = h + 1.0 / v;
    }
    v = v + 1.0;  // 80
  }
  const fcomplex r = 1.0 / (v * v);  // 90
  fcomplex ser = t.b[t.nt] * (2 * t.nt - 1);
  for (int j = 2; j <= t.nt; ++j) {
    const int k = t.nt + 1 - j;
    ser = t.b[k] * (2 * k - 1) + r * ser;  // 100
  }
  fcomplex result = log(v) - 0.5 / v - r * ser - h;
  if (x >= 0.0) return result;
  h = t.pi * u;
  return result + 1.0 / u + t.pi * cos(h) / sin(h);
}

// ---------------------------------------------------------------------------
fcomplex cf1c(fcomplex x, fcomplex eta, fcomplex zl, double eps, fcomplex& fcl, fcomplex& tpk1, bool etane0,
              int limit, double& err, int& nfp, double acch, double fpmin, double fpmax, unsigned& warnings) {
  fcl = 1.0;
  const fcomplex xi = 1.0 / x;
  fcomplex pk = zl + 1.0;
  const double px = real(pk + limit);
  fcomplex ek, rk2, f, pk1, tk, sl, d, df;
  for (;;) {  // 10
    ek = eta / pk;
    rk2 = 1.0 + ek * ek;
    f = (ek + pk * xi) * fcl + (fcl - 1.0) * xi;
    pk1 = pk + 1.0;
    tpk1 = pk + pk1;
    tk = tpk1 * (xi + ek / pk1);
    if (!etane0 || absc(tk) > acch) break;
    // b1 = 0 (e.g. negative eta): exact fix-up, two steps at once.
    const fcomplex e1 = eta / pk1;
    fcl = rk2 / (1.0 + e1 * e1);
    sl = tpk1 * xi * (tpk1 + 2.0) * xi;
    pk = 2.0 + pk;
  }
  d = 1.0 / tk;  // 20
  df = -(fcl * rk2 * d);
  if (real(pk) > real(zl) + 2.0) fcl = -(rk2 * sl);
  fcl = fcl * d * tpk1 * xi;
  f = f + df;
  // CF1 loop on pk = k = lambda + 1.
  double rk = 1.0;
  const double small = std::sqrt(fpmin);
  for (;;) {  // 30
    pk = pk1;
    pk1 = pk1 + 1.0;
    tpk1 = pk + pk1;
    if (etane0) {
      ek = eta / pk;
      rk2 = 1.0 + ek * ek;
    }
    tk = tpk1 * (xi + ek / pk1);
    d = tk - d * rk2;
    if (!(absc(d) > acch)) {
      warnings |= kWarnCf1AccuracyLoss;
      rk = rk + 1.0;
      if (rk > 2.0) {  // 50
        warnings |= kWarnCf1NotConverged;
        err = 2.0;
        return f;
      }
    }
    d = 1.0 / d;  // 40
    fcl = fcl * d * tpk1 * xi;
    if (absc(fcl) < small) fcl = fcl / small;
    if (absc(fcl) > fpmax) fcl = fcl / fpmax;
    df = df * (d * tk - 1.0);
    f = f + df;
    if (real(pk) > px) {  // 50
      warnings |= kWarnCf1NotConverged;
      err = 2.0;
      return f;
    }
    if (!(absc(df) >= absc(f) * eps)) break;
  }
  nfp = int(real(pk - zl - 1));
  err = eps * double(std::sqrt(float(nfp)));  // EPS * SQRT(REAL(NFP)): single-precision root
  return f;
}

double cf1r(double x, double eta, double zl, double eps, double& fcl, double& tpk1, bool etane0, int limit,
            double& err, int& nfp, double acch, double fpmin, double fpmax, unsigned& warnings) {
  fcl = 1.0;
  const double xi = 1.0 / x;
  double pk = zl + 1.0;
  const double px = pk + limit;
  double ek, rk2, f, pk1, tk, sl = 0, d, df;
  for (;;) {  // 10
    ek = eta / pk;
    rk2 = 1.0 + ek * ek;
    f = (ek + pk * xi) * fcl + (fcl - 1.0) * xi;
    pk1 = pk + 1.0;
    tpk1 = pk + pk1;
    tk = tpk1 * (xi + ek / pk1);
    if (!etane0 || std::fabs(tk) > acch) break;
    const double e1 = eta / pk1;
    fcl = rk2 / (1.0 + e1 * e1);
    sl = tpk1 * xi * (tpk1 + 2.0) * xi;
    pk = 2.0 + pk;
  }
  d = 1.0 / tk;  // 20
  df = -fcl * rk2 * d;
  if (pk > zl + 2.0) fcl = -rk2 * sl;
  fcl = fcl * d * tpk1 * xi;
  f = f + df;
  double rk = 1.0;
  const double small = std::sqrt(fpmin);
  for (;;) {  // 30
    pk = pk1;
    pk1 = pk1 + 1.0;
    tpk1 = pk + pk1;
    if (etane0) {
      ek = eta / pk;
      rk2 = 1.0 + ek * ek;
    }
    tk = tpk1 * (xi + ek / pk1);
    d = tk - d * rk2;
    if (!(std::fabs(d) > acch)) {
      warnings |= kWarnCf1AccuracyLoss;
      rk = rk + 1.0;
      if (rk > 2.0) {
        warnings |= kWarnCf1NotConverged;
        err = 2.0;
        return f;
      }
    }
    d = 1.0 / d;  // 40
    fcl = fcl * d * tpk1 * xi;
    if (std::fabs(fcl) < small) fcl = fcl / small;
    if (std::fabs(fcl) > fpmax) fcl = fcl / fpmax;
    df = df * (d * tk - 1.0);
    f = f + df;
    if (pk > px) {
      warnings |= kWarnCf1NotConverged;
      err = 2.0;
      return f;
    }
    if (!(std::fabs(df) >= std::fabs(f) * eps)) break;
  }
  nfp = int(pk - zl - 1);
  err = eps * double(std::sqrt(float(nfp)));
  return f;
}

// ---------------------------------------------------------------------------
fcomplex cf2(fcomplex x, fcomplex eta, fcomplex zl, fcomplex pm, double eps, int limit, double& err, int& npq,
             double acc8, double acch, double accur, unsigned& warnings) {
  const double ta = 2.0 * limit;
  const fcomplex e2mm1 = eta * eta + zl * zl + zl;
  const fcomplex etap = eta * pm;
  const fcomplex xi = 1.0 / x;
  const fcomplex wi = 2.0 * etap;
  double rk = 0.0;
  fcomplex pq = (1.0 - eta * xi) * pm;
  fcomplex aa = -e2mm1 + etap;
  fcomplex bb = 2.0 * (x - eta + pm);
  fcomplex rl = xi * pm;
  if (absc(bb) < acch) {
    rl = rl * aa / (aa + rk + wi);
    pq = pq + rl * (bb + 2.0 * pm);
    aa = aa + 2.0 * (rk + 1.0 + wi);
    bb = bb + (2.0 + 2.0) * pm;
    rk = rk + (2.0 + 2.0);
  }
  fcomplex dd = 1.0 / bb;
  fcomplex dl = aa * dd * rl;
  do {  // 10
    pq = pq + dl;
    rk = rk + 2.0;
    aa = aa + rk + wi;
    bb = bb + 2.0 * pm;
    dd = 1.0 / (aa * dd + bb);
    dl = dl * (bb * dd - 1.0);
    err = absc(dl) / absc(pq);
  } while (err >= std::max(eps, acc8 * rk * 0.5) && rk <= ta);
  npq = int(rk / 2.0);
  pq = pq + dl;
  if (npq >= limit - 1 && err > accur) warnings |= kWarnCf2NotConverged;
  return pq;
}

// ---------------------------------------------------------------------------
fcomplex f11(fcomplex x, fcomplex eta, fcomplex zl, fcomplex p, double eps, int limit, int kind, double& err,
             int& nits, double fpmax, double acc8, double acc16, unsigned& warnings) {
  const fcomplex aa = zl + 1.0 - eta * p;
  const fcomplex bb = 2.0 * (zl + 1.0);
  const fcomplex z = 2.0 * p * x;
  const bool zllin = real(bb) <= 0.0 && abs(bb - nintc(bb)) < std::pow(acc8, 0.25);
  if (zllin && !(real(bb) + limit < 1.5)) {
    nits = -1;
    return 0.0;
  }
  if (limit <= 0) {  // 10
    err = 0.0;
    nits = 1;
    return 0.0;
  }
  double ta = 1.0, rk = 1.0;
  int i;
  if (kind <= 0 && absc(z) * absc(aa) > absc(bb) * 1.0) {
    // Real arithmetic in extended precision on the real and imaginary parts.
    wide dr = 1.0, di = 0.0, gr = 1.0, gi = 0.0;
    wide ar = real(aa), br = real(bb), fi = 0.0;
    const wide aai = imag(aa), bbi = imag(bb), zr = real(z), zi = imag(z);
    for (i = 2; i <= limit; ++i) {
      const wide fi1 = fi + wide(1.0);
      wide tr = br * fi1;
      wide ti = bbi * fi1;
      const wide den = wide(1.0) / (tr * tr + ti * ti);
      const wide ur = (ar * tr + aai * ti) * den;
      const wide ui = (aai * tr - ar * ti) * den;
      tr = ur * gr - ui * gi;
      ti = ur * gi + ui * gr;
      gr = zr * tr - zi * ti;
      gi = zr * ti + zi * tr;
      dr = dr + gr;
      di = di + gi;
      err = to_double(wabs(gr) + wabs(gi));
      if (err > fpmax) {  // 60, with F11 = 0
        nits = i;
        return 0.0;
      }
      rk = to_double(wabs(dr) + wabs(di));
      ta = std::max(ta, rk);
      if (err < rk * eps || (i >= 4 && err < acc16)) break;  // 30
      fi = fi1;
      ar = ar + wide(1.0);
      br = br + wide(1.0);
    }
    // 30: F11 = DR + CI*DI, evaluated in COMPLEX*32 and rounded to COMPLEX*16.
    const wide ci_di_re = wide(0.0) * di - wide(1.0) * wide(0.0);
    const wide ci_di_im = wide(0.0) * wide(0.0) + wide(1.0) * di;
    const fcomplex value{to_double(dr + ci_di_re), to_double(wide(0.0) + ci_di_im)};
    err = acc16 * ta / rk;
    nits = i;
    return value;
  }
  // Complex double arithmetic (KIND > 0, or no cancellation expected).
  fcomplex g = 1.0;
  fcomplex f = 1.0;
  if (kind >= 2) f = cdigam(aa, warnings) - cdigam(bb, warnings) - cdigam(g, warnings);
  fcomplex dd = f;
  for (i = 2; i <= limit; ++i) {
    const fcomplex ai = aa + (i - 2);
    const fcomplex bi = bb + (i - 2);
    const double r = i - 1.0;
    g = g * z * ai / (bi * r);
    if (kind >= 2) f = f + 1.0 / ai - 1.0 / bi - 1.0 / r;  // psi(a+r) - psi(b+r) - psi(1+r)
    const fcomplex t = g * f;
    dd = dd + t;
    err = absc(t);
    if (err > fpmax) {  // 60, with F11 = 0
      nits = i;
      return 0.0;
    }
    rk = absc(dd);
    ta = std::max(ta, rk);
    if (err < rk * eps || (err < acc8 && i >= 4)) break;  // 50
  }
  err = acc8 * ta / rk;  // 50
  nits = i;
  return dd;
}

// ---------------------------------------------------------------------------
bool rcf(const fcomplex* a, fcomplex* b, int ibeg, int& inum, fcomplex (*xx)[kJmax + 2], double eps,
         RcfState& s) {
  const int ibn = inum;  // b(ibn) is the last value set on this call
  enum class From { start, l20, l50 } from = From::start;
  if (ibeg > 4) {
    if (s.m != ibeg - 1) return false;  // 90: restart does not continue the last call
    from = From::l50;
  } else if (ibeg == 4) {
    from = From::l20;
  }
  if (from == From::start) {
    b[1] = a[1];
    if (ibn >= 2) b[2] = -(a[2] / a[1]);
    if (ibn < 3) return true;  // 10
    s.x0 = a[3] / a[2];
    xx[2][1] = b[2];
    xx[1][1] = -s.x0;
    xx[1][2] = 0.0;
    b[3] = -s.x0 - b[2];
    s.x0 = -(b[3] * a[2]);
    s.m = 3;
    s.mp12 = 2;
    s.even = true;
    if (!(ibn > 3)) return true;
    from = From::l20;
  }
  if (from == From::l20) {
    if (abs(b[3]) < eps * abs(s.x0)) {  // 80
      inum = -s.m;
      return true;
    }
    s.m = 4;
  }
  for (;;) {
    if (from != From::l50) {  // 30
      s.x1 = a[s.m];
      s.m2m1 = s.mp12;
      s.mp12 = s.m2m1 + 1;
      if (s.even) s.mp12 = s.m2m1;
      for (int k = 2; k <= s.mp12; ++k) s.x1 = s.x1 + a[s.m - k + 1] * xx[1][k - 1];  // 40
      b[s.m] = -(s.x1 / s.x0);
      if (s.m >= ibn) return true;
    }
    from = From::start;
    if (abs(b[s.m]) < eps * abs(s.x0)) {  // 50 -> 80
      inum = -s.m;
      return true;
    }
    int k = s.m2m1;
    do {  // 60
      xx[2][k] = xx[1][k] + b[s.m] * xx[2][k - 1];
      k = k - 1;
    } while (k > 1);
    xx[2][1] = xx[1][1] + b[s.m];
    for (k = 1; k <= s.m2m1; ++k) {  // 70
      s.x0 = xx[2][k];
      xx[2][k] = xx[1][k];
      xx[1][k] = s.x0;
    }
    s.x0 = s.x1;
    xx[1][s.m2m1 + 1] = 0.0;
    s.m = s.m + 1;
    s.even = !s.even;
  }
}

// ---------------------------------------------------------------------------
fcomplex f20(fcomplex aa, fcomplex bb, fcomplex z, double eps, int jmax, double& re, double fpmax, int& n,
             PadeWork& w, unsigned& warnings) {
  fcomplex* x = w.a;  // X(I,1): terms of the series
  fcomplex* c = w.b;  // X(I,2): continued-fraction coefficients
  re = 0.0;
  x[1] = 1.0;
  fcomplex sum = x[1];
  double atl = absc(x[1]);
  fcomplex f = sum, d = 1.0, df = sum;
  int j = 0;
  const double ep = eps * jmax * 10.0;
  const int ma = -nintc(aa);
  const int mb = -nintc(bb);
  const bool finite = (std::fabs(std::fabs(real(aa)) - ma) < ep && std::fabs(imag(aa)) < ep) ||
                      (std::fabs(std::fabs(real(bb)) - mb) < ep && std::fabs(imag(bb)) < ep);
  int imax = jmax;
  if (finite && ma >= 0) imax = std::min(ma + 1, imax);
  if (finite && mb >= 0) imax = std::min(mb + 1, imax);
  double at = 0;
  int i;
  for (i = 2; i <= imax; ++i) {
    x[i] = x[i - 1] * z * (aa + i - 2) * (bb + i - 2) / (i - 1);
    if (absc(x[i]) > fpmax) {  // 40
      n = 0;
      if (!finite) re = at / absc(sum);
      return sum;
    }
    at = absc(x[i]);
    if (j == 0) {
      sum = sum + x[i];
      if (at < absc(sum) * eps) {  // 20
        n = i;
        if (!finite) re = at / absc(sum);
        return sum;
      }
    }
    if (!finite) {
      if (j > 0 || at > atl || i >= jmax - 2) j = j + 1;
      if (j != 0) {
        int inum = i;
        if (!rcf(x, c, j, inum, w.xx, eps, w.rcf)) {
          warnings |= kWarnRcfRestart;
          inum = -1;
        }
        if (inum < 0) {  // 40
          n = 0;
          if (!finite) re = at / absc(sum);
          return sum;
        }
        for (int k = std::max(j, 2); k <= i; ++k) {  // 50
          d = 1.0 / (d * c[k] + 1.0);
          df = df * (d - 1.0);
          f = f + df;
          if (absc(df) < absc(f) * eps || (df == 0.0 && f == 0.0 && i >= 4)) {  // 30
            re = absc(df) / absc(f);
            n = k;
            return f;
          }
        }
        j = i;
      }
    }
    atl = at;  // 10
  }
  if (!finite) i = -jmax;
  n = i;  // 20
  if (!finite) re = at / absc(sum);
  return sum;
}

// ---------------------------------------------------------------------------
fcomplex cf1a(fcomplex rho, fcomplex eta, fcomplex xl, fcomplex psi, double eps, int nmax, int& nused,
              fcomplex& fcl, double& re, double fpmax, PadeWork& w, unsigned& warnings) {
  fcomplex* g = w.a;  // G(N): differences of the partial sums
  fcomplex* c = w.b;  // C(N): continued-fraction coefficients
  double t1 = std::sin(real(psi));
  double t2 = std::cos(real(psi));
  double atl = std::tanh(imag(psi));
  // cos(psi)/cosh(Im psi), which always has the correct sign.
  const fcomplex cosl{t2, -(t1 * atl)};
  const fcomplex tanl = fcomplex(t1, t2 * atl) / cosl;
  re = 0.0;
  const fcomplex xll1 = xl * (xl + 1.0);
  const fcomplex etasq = eta * eta;
  fcomplex sl1 = 1.0, sl = sl1, sc1 = 0.0, sc = sc1, tl1 = sc, tl = tl1;
  fcomplex tc1 = 1.0 - eta / rho, tc = tc1;
  fcl = tl + sl * tanl;
  g[1] = (tc + sc * tanl) / fcl;
  fcomplex glast = g[1];
  atl = absc(glast);
  fcomplex f = glast, d = 1.0, df = glast;
  int j = 0;
  bool failed = false;
  for (int n = 2; n <= nmax; ++n) {
    t1 = n - 1;
    t2 = 2.0 * t1 - 1.0;
    const double t3 = t1 * (t1 - 1.0);
    const fcomplex denom = 2.0 * rho * t1;
    const fcomplex c1 = (eta * t2) / denom;
    const fcomplex c2 = (etasq + xll1 - t3) / denom;
    const fcomplex sl2 = c1 * sl1 - c2 * tl1;
    const fcomplex tl2 = c1 * tl1 + c2 * sl1;
    const fcomplex sc2 = c1 * sc1 - c2 * tc1 - sl2 / rho;
    const fcomplex tc2 = c1 * tc1 + c2 * sc1 - tl2 / rho;
    sl = sl + sl2;
    tl = tl + tl2;
    sc = sc + sc2;
    tc = tc + tc2;
    sl1 = sl2;
    tl1 = tl2;
    sc1 = sc2;
    tc1 = tc2;
    fcl = tl + sl * tanl;
    if (absc(fcl) > fpmax || absc(fcl) < 1.0 / fpmax) {  // 40
      failed = true;
      break;
    }
    const fcomplex gsum = (tc + sc * tanl) / fcl;
    g[n] = gsum - glast;
    glast = gsum;
    const double at = absc(g[n]);
    if (at < absc(gsum) * eps) {  // 20
      fcl = fcl * cosl;
      re = at / absc(gsum);
      nused = n;
      return gsum;
    }
    if (j > 0 || at > atl || n >= nmax - 2) j = j + 1;
    if (j != 0) {
      int inum = n;
      if (!rcf(g, c, j, inum, w.xx, eps, w.rcf)) {
        warnings |= kWarnRcfRestart;
        inum = -1;
      }
      if (inum < 0) {  // 40
        failed = true;
        break;
      }
      for (int k = std::max(j, 2); k <= n; ++k) {  // 60
        d = 1.0 / (d * c[k] + 1.0);
        df = df * (d - 1.0);
        f = f + df;
        if (absc(df) < absc(f) * eps || (df == 0.0 && f == 0.0 && n >= 4)) {  // 30
          fcl = fcl * cosl;
          re = absc(df) / absc(f);
          nused = k;
          return f;
        }
      }
      j = n;
    }
    atl = at;  // 10
  }
  if (!failed) {  // loop exhausted: K = -NMAX, GO TO 30
    fcl = fcl * cosl;
    re = absc(df) / absc(f);
    nused = -nmax;
    return f;
  }
  fcl = 1.0;  // 40
  re = 1.0;
  nused = 0;
  return g[1];
}

}  // namespace detail
}  // namespace coulcc
