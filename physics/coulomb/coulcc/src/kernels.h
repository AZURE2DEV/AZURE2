// The subprograms of coulcc.f other than COULCC itself: continued fractions,
// hypergeometric series, log-gamma and TIDY.  Each function is the Fortran
// routine of the same name with the same arguments and the same arithmetic;
// arguments that the Fortran writes are reference parameters and are written
// in exactly the cases where the Fortran writes them (some callers rely on an
// argument being left unchanged).  Messages become bits in `warnings`.
// DESIGN.md section 3.2.
#ifndef COULCC_KERNELS_H
#define COULCC_KERNELS_H

#include "fortran_complex.h"

namespace coulcc {
namespace detail {

// Machine-dependent constants of coulcc.f (DATA statements), kept unchanged.
constexpr double kAcc8 = 2e-16;   // "smallest number with 1+ACC8 .ne. 1" as set in the Fortran
constexpr double kAcc16 = 3e-33;  // the same for REAL*16
constexpr double kFpmax = 1e60;
constexpr double kFpmin = 1e-60;
constexpr double kFplmax = 140.0;
constexpr double kFplmin = -140.0;
constexpr int kLimit = 20000;  // maximum iterations of CF1, CF2 and the 1F1 series
constexpr int kJmax = 50;      // size of the Pade work arrays of CF1A and F20

// Warning bits (the Fortran WRITE statements); public names in coulcc.h.
enum WarningBit : unsigned {
  kWarnReflectionNotUsedL = 1u << 0,
  kWarnReflectionNotUsedX = 1u << 1,
  kWarnLinearIndependenceLost = 1u << 2,
  kWarnIrregularOverflow = 1u << 3,
  kWarnIrregularUnderflow = 1u << 4,
  kWarnRegularOverflow = 1u << 5,
  kWarnRegularUnderflow = 1u << 6,
  kWarnCf1aAndCf1Unstable = 1u << 7,
  kWarnF11Overflow = 1u << 8,
  kWarnPolesAndInstabilities = 1u << 9,
  kWarnRoundoff = 1u << 10,
  kWarnCf1AccuracyLoss = 1u << 11,
  kWarnCf1NotConverged = 1u << 12,
  kWarnCf2NotConverged = 1u << 13,
  kWarnLogGammaPole = 1u << 14,
  kWarnRcfRestart = 1u << 15,
  kWarnIrregularTooSmallRho = 1u << 16,
};

// State that RCF keeps between successive calls for one series: COMMON
// /RCFCM2/ X1, M2M1, MP12, EVEN, M, and RCF's local X0, which the Fortran
// needs to be SAVEd.  One RcfState per CF1A or F20 evaluation.
struct RcfState {
  fcomplex x0, x1;
  int m2m1 = 0, mp12 = 0, m = 0;
  bool even = false;
};

// Work space of CF1A and F20 (the array XRCF(JMAX,4) of COULCC): the series
// coefficients A(k), the continued-fraction coefficients B(k) and RCF's
// auxiliary vector XX(2,k), all 1-based.
struct PadeWork {
  fcomplex a[kJmax + 1];
  fcomplex b[kJmax + 1];
  fcomplex xx[3][kJmax + 2];  // xx[1][k], xx[2][k] = XX(1,k), XX(2,k)
  RcfState rcf;
};

// TIDY: parts smaller than 5*acc*(|Re z|+|Im z|) set to zero.
fcomplex tidy(fcomplex z, double acc);

// CLOGAM: ln Gamma(z); CDIGAM: psi(z) = d ln Gamma / dz.  At a non-positive
// integer both return 0 and set kWarnLogGammaPole.
fcomplex clogam(fcomplex z, unsigned& warnings);
fcomplex cdigam(fcomplex z, unsigned& warnings);

// CF1C: f = F'(zl,eta,x)/F(zl,eta,x) by the continued fraction CF1 in complex
// arithmetic.  fcl and tpk1 are written as the Fortran writes them; on
// failure err = 2 and nfp is not written.
fcomplex cf1c(fcomplex x, fcomplex eta, fcomplex zl, double eps, fcomplex& fcl, fcomplex& tpk1, bool etane0,
              int limit, double& err, int& nfp, double acch, double fpmin, double fpmax, unsigned& warnings);

// CF1R: the same in real arithmetic (real x, eta, zl).
double cf1r(double x, double eta, double zl, double eps, double& fcl, double& tpk1, bool etane0, int limit,
            double& err, int& nfp, double acch, double fpmin, double fpmax, unsigned& warnings);

// CF2: p + i omega q = H'/H for H = H^omega(eta,x) at lambda zl, pm = i omega.
fcomplex cf2(fcomplex x, fcomplex eta, fcomplex zl, fcomplex pm, double eps, int limit, double& err, int& npq,
             double acc8, double acch, double accur, unsigned& warnings);

// F11: the series 1F1(zl+1-eta p; 2(zl+1); 2 p x) to accuracy eps with at most
// limit terms; kind 0 in extended precision (REAL*16), 1 in complex double, 2
// with the digamma factor of the logarithmic solution.  When the series
// cannot start (b a non-positive integer) nits = -1 and err is not written.
fcomplex f11(fcomplex x, fcomplex eta, fcomplex zl, fcomplex p, double eps, int limit, int kind, double& err,
             int& nits, double fpmax, double acc8, double acc16, unsigned& warnings);

// F20: the series 2F0(aa, bb;; z), continued into a continued fraction by RCF
// when the terms diverge.  n > 0 on success, n <= 0 on failure.
fcomplex f20(fcomplex aa, fcomplex bb, fcomplex z, double eps, int jmax, double& re, double fpmax, int& n,
             PadeWork& work, unsigned& warnings);

// CF1A: f = F'/F from the asymptotic (2F0) expansions of H+ and H-, valid for
// Re rho > 0 and best for |rho| >> eta^2, xl.  nused < 0 or > 0 as in the
// Fortran; fcl is the asymptotic estimate of F.
fcomplex cf1a(fcomplex rho, fcomplex eta, fcomplex xl, fcomplex psi, double eps, int nmax, int& nused,
              fcomplex& fcl, double& re, double fpmax, PadeWork& work, unsigned& warnings);

// RCF: converts the series a(1..inum) into continued-fraction coefficients
// b(ibeg..inum) (Patry and Gupta).  inum becomes -m on breakdown.  Returns
// false when a restart is requested that does not continue the previous call
// (the Fortran stops there).
bool rcf(const fcomplex* a, fcomplex* b, int ibeg, int& inum, fcomplex (*xx)[kJmax + 2], double eps,
         RcfState& state);

}  // namespace detail
}  // namespace coulcc

#endif
