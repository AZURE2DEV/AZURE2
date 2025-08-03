#ifndef COMPLEX_FUNCTIONS_H
#define COMPLEX_FUNCTIONS_H

#define _USE_MATH_DEFINES // for C++
#if defined(_USE_MATH_DEFINES) && !defined(_MATH_DEFINES_DEFINED)
#define _MATH_DEFINES_DEFINED

/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG2E    - log2(e)
 * M_LOG10E   - log10(e)
 * M_LN2      - ln(2)
 * M_LN10     - ln(10)
 * M_PI       - pi
 * M_PI_2     - pi/2
 * M_PI_4     - pi/4
 * M_1_PI     - 1/pi
 * M_2_PI     - 2/pi
 * M_2_SQRTPI - 2/sqrt(pi)
 * M_SQRT2    - sqrt(2)
 * M_SQRT1_2  - 1/sqrt(2)
 */

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401

#endif  /* _USE_MATH_DEFINES */

#include <complex>
#include <iostream>
#include <cstdlib>
#define SIGN(a) (((a) < 0) ? (-1) : (1))
const double precision = 1E-10,sqrt_precision = 1E-5;

// Infinite norm of a complex number.
// ----------------------------------
// It is max(|Re[z]|,|Im[z]|)

inline double inf_norm (const std::complex<double> &z)
{
  return std::max (std::abs (real (z)),std::abs (imag (z)));
}

// Test of finiteness of a complex number
// --------------------------------------
// If real or imaginary parts are finite, true is returned.
// Otherwise, false is returned

inline bool isfinite (const std::complex<double> &z)
{
  const double x = real (z), y = imag (z);

  return (std::isfinite (x) && std::isfinite (y));
}

// Usual operator overloads of complex numbers with integers
// ---------------------------------------------------------
// Recent complex libraries do not accept for example z+n or z==n with n integer, signed or unsigned.
// The operator overload is done here, by simply putting a cast on double to the integer.

inline std::complex<double> operator + (const std::complex<double> &z,const int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}

inline std::complex<double> operator + (const std::complex<double> &z,const unsigned int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const unsigned int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const unsigned int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const unsigned int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}

inline bool operator == (const std::complex<double> &z,const int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}

inline bool operator == (const std::complex<double> &z,const unsigned int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const unsigned int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}

extern std::complex<double> expm1 (const std::complex<double> &z);

extern std::complex<double> log1p (const std::complex<double> &z);

extern std::complex<double> log_Gamma (const std::complex<double> &z);

extern std::complex<double> sigma_l_calc (const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> log_Cl_eta_calc (const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> log_cut_constant_AS_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> log_cut_constant_CFa_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> log_cut_constant_CFb_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> sin_chi_calc (const std::complex<double> &l,const std::complex<double> &eta);

extern std::complex<double> exp_I_omega_chi_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta);

#endif
