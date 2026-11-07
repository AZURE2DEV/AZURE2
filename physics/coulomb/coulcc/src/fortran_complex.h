// Complex arithmetic with the semantics of COMPLEX*16 in gfortran, and the
// statement functions of coulcc.f (ABSC, NINTC, NPINT).
//
// Fortran converts a real or integer operand of a mixed expression to complex
// and then performs a full complex operation; std::complex works on the
// components.  The results differ only in the sign of zero parts, but that
// sign selects the side of the branch cuts of log, sqrt and atan2, so it can
// change the sign of a Coulomb or Bessel function on the negative real axis.
// Division is Smith's algorithm with gfortran's branch condition.  DESIGN.md
// section 3.3.
#ifndef COULCC_FORTRAN_COMPLEX_H
#define COULCC_FORTRAN_COMPLEX_H

#include <cmath>
#include <complex>

namespace coulcc {
namespace detail {

struct fcomplex {
  double re = 0, im = 0;
  constexpr fcomplex() = default;
  constexpr fcomplex(double r) : re(r), im(0) {}  // NOLINT: Fortran promotion of a real
  constexpr fcomplex(double r, double i) : re(r), im(i) {}
  fcomplex(const std::complex<double>& z) : re(z.real()), im(z.imag()) {}  // NOLINT
  std::complex<double> std() const { return {re, im}; }
};

// Every operand is promoted to complex first: (r) -> (r, +0).
inline fcomplex operator+(fcomplex a, fcomplex b) { return {a.re + b.re, a.im + b.im}; }
inline fcomplex operator-(fcomplex a, fcomplex b) { return {a.re - b.re, a.im - b.im}; }
inline fcomplex operator-(fcomplex a) { return {-a.re, -a.im}; }
inline fcomplex operator*(fcomplex a, fcomplex b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
inline fcomplex operator/(fcomplex a, fcomplex b) {
  // gfortran (tree-complex.c, expand_complex_div_wide), Smith's algorithm.
  if (std::fabs(b.re) < std::fabs(b.im)) {
    const double ratio = b.re / b.im;
    const double div = b.re * ratio + b.im;
    return {(a.re * ratio + a.im) / div, (a.im * ratio - a.re) / div};
  }
  const double ratio = b.im / b.re;
  const double div = b.im * ratio + b.re;
  return {(a.im * ratio + a.re) / div, (a.im - a.re * ratio) / div};
}
// Real operands: exact-match overloads, so that a double argument is never
// converted to int by the integer overloads below.
inline fcomplex operator+(fcomplex a, double b) { return a + fcomplex(b); }
inline fcomplex operator-(fcomplex a, double b) { return a - fcomplex(b); }
inline fcomplex operator*(fcomplex a, double b) { return a * fcomplex(b); }
inline fcomplex operator/(fcomplex a, double b) { return a / fcomplex(b); }
inline fcomplex operator+(double a, fcomplex b) { return fcomplex(a) + b; }
inline fcomplex operator-(double a, fcomplex b) { return fcomplex(a) - b; }
inline fcomplex operator*(double a, fcomplex b) { return fcomplex(a) * b; }
inline fcomplex operator/(double a, fcomplex b) { return fcomplex(a) / b; }
inline bool operator==(fcomplex a, double b) { return a.re == b && a.im == 0.0; }
// Integer operands are promoted through double, as Fortran does.
inline fcomplex operator+(fcomplex a, int b) { return a + fcomplex(double(b)); }
inline fcomplex operator-(fcomplex a, int b) { return a - fcomplex(double(b)); }
inline fcomplex operator*(fcomplex a, int b) { return a * fcomplex(double(b)); }
inline fcomplex operator+(int a, fcomplex b) { return fcomplex(double(a)) + b; }
inline fcomplex operator-(int a, fcomplex b) { return fcomplex(double(a)) - b; }
inline fcomplex operator*(int a, fcomplex b) { return fcomplex(double(a)) * b; }
inline fcomplex operator/(fcomplex a, int b) { return a / fcomplex(double(b)); }
inline bool operator==(fcomplex a, fcomplex b) { return a.re == b.re && a.im == b.im; }

inline double real(fcomplex z) { return z.re; }
inline double imag(fcomplex z) { return z.im; }
inline fcomplex conj(fcomplex z) { return {z.re, -z.im}; }
// Elementary functions: the C library's, as gfortran calls them (cabs, cexp, clog, csqrt, csin, ccos).
inline double abs(fcomplex z) { return std::abs(z.std()); }
inline fcomplex exp(fcomplex z) { return std::exp(z.std()); }
inline fcomplex log(fcomplex z) { return std::log(z.std()); }
inline fcomplex sqrt(fcomplex z) { return std::sqrt(z.std()); }
inline fcomplex sin(fcomplex z) { return std::sin(z.std()); }
inline fcomplex cos(fcomplex z) { return std::cos(z.std()); }

const fcomplex CI{0.0, 1.0};

// Fortran SIGN(A,B) = |A| with the sign of B (gfortran honours the sign of zero).
inline double sign(double a, double b) { return std::copysign(std::fabs(a), b); }

// ABSC(W) = |Re W| + |Im W|.
inline double absc(fcomplex w) { return std::fabs(w.re) + std::fabs(w.im); }
// NINTC(W) = NINT(REAL(REAL(W))): the inner REAL of a REAL*8 gives default
// (single-precision) REAL, so the value is rounded to float before NINT.
inline int nintc(fcomplex w) { return int(std::lround(float(w.re))); }
// NPINT(W,ACCB): W is within ACCB of a non-positive integer.
inline bool npint(fcomplex w, double accb) { return absc(nintc(w) - w) < accb && w.re < 0.5; }

}  // namespace detail
}  // namespace coulcc

#endif
