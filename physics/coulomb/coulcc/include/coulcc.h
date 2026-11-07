// COULCC: complex Coulomb wave functions by Steed's method.
//
// C++17 translation of COULCC by I. J. Thompson and A. R. Barnett, Comput.
// Phys. Commun. 36 (1985) 363 (coulcc.f in this directory), computing the
// same numbers bit for bit (README.md, DESIGN.md).
//
// For complex rho, eta and lambda the functions returned are the regular and
// irregular solutions of
//     u'' + [1 - 2 eta/rho - lambda(lambda+1)/rho^2] u = 0
// normalised so that, for large |rho|,
//     F -> sin(theta), G -> cos(theta), H+- = G +- i F -> exp(+-i theta),
//     theta = rho - eta ln(2 rho) - lambda pi/2 + sigma_lambda,
//     sigma_lambda = [ln Gamma(lambda+1+i eta) - ln Gamma(lambda+1-i eta)] / 2i,
// with derivatives d/drho and the Wronskian F'G - FG' = 1.  rho = k r and eta
// = Z1 Z2 e^2 mu / (hbar^2 k) are dimensionless.  With eta = 0 the same code
// gives spherical, cylindrical and modified Bessel functions.
//
// Thread-safe: no mutable global or static state.  Nothing is printed;
// failures are reported in ifail/status and warnings.
#ifndef COULCC_H
#define COULCC_H

#include <complex>
#include <vector>

namespace coulcc {

using complex = std::complex<double>;

// Which functions (the Fortran argument KFN).
enum class Functions : int {
  Coulomb = 0,             // F, G (or H+-) and the Coulomb phase shifts sigma
  CoulombNoPhase = -1,     // F, G (or H+-); sigma not computed
  SphericalBessel = 1,     // j_lambda, y_lambda (eta ignored)
  CylindricalBessel = 2,   // J_nu, Y_nu with nu = lambda (eta ignored)
  ModifiedBessel = 3,      // I_nu, K_nu (eta ignored; H+- not available)
};

// The second (irregular) function (the Fortran |MODE1|/10).
enum class Irregular : int {
  G = 0,       // G, y_lambda or Y_nu
  Hplus = 1,   // H+ = G + iF, h(1) = j + iy, or H(1) = J + iY
  Hminus = 2,  // H- = G - iF, h(2), H(2)
};

struct Options {
  Functions kind = Functions::Coulomb;
  Irregular second = Irregular::G;
  bool irregular = true;     // compute the second function (MODE 1, 2) or F only (MODE 3, 4)
  bool derivatives = true;   // compute derivatives (MODE 1, 3)
  bool scaled = false;       // MODE1 < 0: values multiplied by an exponential factor, see below
  // Scaling: with s = Im rho (Re rho for ModifiedBessel), F (j, J, I) and G
  // (y, Y) are multiplied by exp(-|s|), H+ (h(1), H(1), K) by exp(s), H- (h(2),
  // H(2)) by exp(-s).

  int mode1() const;                                // the Fortran MODE1 equivalent
  static Options from_fortran(int mode1, int kfn);  // inverse; mode1, kfn as in coulcc.f
};

// The Fortran IFAIL, by name.
enum class Status : int {
  Ok = 0,                 // all values satisfactory
  Partial = 1,            // IFAIL > 0: the last IFAIL values over- or underflowed
  ContinuedFraction = -1, // a continued fraction failed, or the arithmetic check before the last recursion
  RhoTooSmall = -2,       // |rho| too small for the irregular solution
  NotFoundAtMin = -3,     // values at lambda_min not found (over- or underflow)
  Roundoff = -4,          // roundoff errors make the results meaningless (rerr > 0.1)
  InvalidArgument = -5,   // nl < 1, or MODE1/KFN outside the documented values (not in the Fortran)
};

// Conditions the Fortran reports by printing a message.
enum Warning : unsigned {
  ReflectionNotUsedL = 1u << 0,        // lambda <= -1: reflection formula advised, accuracy reduced
  ReflectionNotUsedX = 1u << 1,        // Re rho < 0: CF2 inaccurate, reflection advised
  LinearIndependenceLost = 1u << 2,    // F and H(omega) nearly dependent (Coulomb eigenstate, or CF1 unstable)
  IrregularOverflow = 1u << 3,         // irregular solution would exceed 1e60
  IrregularUnderflow = 1u << 4,        // irregular solution would be below 1e-60
  RegularOverflow = 1u << 5,           // regular solution would exceed 1e60
  RegularUnderflow = 1u << 6,          // regular solution would be below 1e-60
  Cf1aAndCf1Unstable = 1u << 7,        // (eta, lambda)/rho too large for CF1A and CF1 unstable
  F11Overflow = 1u << 8,               // overflow in the 1F1 series
  PolesAndInstabilities = 1u << 9,     // bound-state poles and F instabilities: split the lambda range
  RoundoffWarning = 1u << 10,          // rerr > ACCUR^(1/4) = 3e-4
  Cf1AccuracyLoss = 1u << 11,          // CF1 denominator near zero
  Cf1NotConverged = 1u << 12,          // CF1 did not converge in 20000 iterations
  Cf2NotConverged = 1u << 13,          // CF2 did not converge fully
  LogGammaPole = 1u << 14,             // ln Gamma or psi requested at a non-positive integer (0 used)
  RcfRestart = 1u << 15,               // internal inconsistency of the continued-fraction conversion
  IrregularTooSmallRho = 1u << 16,     // |rho| < 1e-60 with the irregular solution requested
};

// Iteration counts and the cases used (COMMON /STEED/ of the Fortran).
struct Diagnostics {
  int nfp = 0;            // CF1 terms (< 0: CF1A terms)
  int n11 = 0;            // 1F1 series terms
  int npq[2] = {0, 0};    // CF2 terms for H+ and H-
  int n20 = 0;            // 2F0 terms
  int kase[2] = {0, 0};   // method for the irregular solution, first and reverse block (1-6)
};

// Outputs other than the arrays, for the Fortran-style call.
struct Info {
  double rerr = 0;        // estimated relative roundoff error of the results
  unsigned warnings = 0;  // bit set of Warning
  Diagnostics diag;
};

struct Result {
  // Index i corresponds to lambda = lambda_min + i.  Entries i >= n_valid are
  // NaN; vectors of functions not requested are empty.
  std::vector<complex> F, G, Fp, Gp, sigma;
  int ifail = 0;           // the Fortran IFAIL
  Status status = Status::Ok;
  int n_valid = 0;         // nl - ifail if ifail >= 0, else 0
  double rerr = 0;         // estimated relative roundoff error
  unsigned warnings = 0;   // bit set of Warning
  Diagnostics diag;
  bool ok() const { return ifail == 0; }
};

// Functions at lambda = lambda_min, ..., lambda_min + nl - 1.
Result coulcc(complex rho, complex eta, complex lambda_min, int nl, const Options& options = Options());

// Fortran-style call with caller-owned arrays of length nl and the MODE1 and
// KFN arguments of coulcc.f; returns IFAIL.  As in the Fortran, arrays not
// required by mode1/kfn are not accessed (they may be null), gc and gcp serve
// as scratch for the recurrence when the irregular solution is requested, and
// sig is written only for kfn = 0.  Allocates nothing.
int coulcc(complex xx, complex eta1, complex zlmin, int nl, complex* fc, complex* gc, complex* fcp, complex* gcp,
           complex* sig, int mode1, int kfn, Info& info);

}  // namespace coulcc

#endif
