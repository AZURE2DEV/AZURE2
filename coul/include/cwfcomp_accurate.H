#ifndef CWFCOMP_ACCURATE_H
#define CWFCOMP_ACCURATE_H

#include <complex>

#define _USE_MATH_DEFINES
#if defined(_USE_MATH_DEFINES) && !defined(_MATH_DEFINES_DEFINED)
#define _MATH_DEFINES_DEFINED

#ifndef M_E
#define M_E        2.71828182845904523536
#endif
#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340736
#endif
#ifndef M_LOG10E
#define M_LOG10E   0.434294481903251827651
#endif
#ifndef M_LN2
#define M_LN2      0.693147180559945309417
#endif
#ifndef M_LN10
#define M_LN10     2.30258509299404568402
#endif
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4     0.785398163397448309616
#endif
#ifndef M_1_PI
#define M_1_PI     0.318309886183790671538
#endif
#ifndef M_2_PI
#define M_2_PI     0.636619772367581343076
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif
#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2  0.707106781186547524401
#endif

#endif  /* _USE_MATH_DEFINES */

// Simplified Coulomb wave functions class using accurate_cwf implementation  
class Coulomb_wave_functions_accurate
{
public:
  // Constructor - matches original interface
  Coulomb_wave_functions_accurate(const bool is_it_normalized_c, const std::complex<double> &l_c, const std::complex<double> &eta_c);
  
  // Destructor
  ~Coulomb_wave_functions_accurate();

  // Initialization method for F and dF at a point (for compatibility)
  void F_dF_init(const std::complex<double> &z, const std::complex<double> &F, const std::complex<double> &dF);

  // Main calculation methods - match original interface
  void F_dF(const std::complex<double> &z, std::complex<double> &F, std::complex<double> &dF); 
  void G_dG(const std::complex<double> &z, std::complex<double> &G, std::complex<double> &dG);
  void H_dH(const int omega, const std::complex<double> &z, std::complex<double> &H, std::complex<double> &dH);
  void H_dH_scaled(const int omega, const std::complex<double> &z, std::complex<double> &H, std::complex<double> &dH);

  void compute(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
               std::complex<double> &G, std::complex<double> &dG) const
  {
      // Core solver using accurate_cwf implementation
      solve_cwf(rho, F, dF, G, dG);
  }

  // Public constants - match original interface
  const std::complex<double> l, eta;
  const bool is_it_normalized;

private:
  // Helper methods for the accurate_cwf implementation
  void solve_cwf(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
                 std::complex<double> &G, std::complex<double> &dG) const;
  
  // Continued fraction method for general complex parameters
  void solve_continued_fraction(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
                               std::complex<double> &G, std::complex<double> &dG) const;
  
  // Helper methods for continued fraction implementation
  void solve_small_rho_series(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
                             std::complex<double> &G, std::complex<double> &dG) const;
  void solve_large_rho_asymptotic(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
                                 std::complex<double> &G, std::complex<double> &dG) const;
  void solve_intermediate_continued_fraction(const std::complex<double> &rho, std::complex<double> &F, std::complex<double> &dF,
                                           std::complex<double> &G, std::complex<double> &dG) const;
  
  // Mathematical helper functions
  std::complex<double> compute_Coulomb_constant(const std::complex<double> &l, const std::complex<double> &eta) const;
  std::complex<double> compute_coulomb_phase_shift(const std::complex<double> &l, const std::complex<double> &eta) const;
  
  // Core continued fraction methods (from original coul library)
  void compute_H_functions(const std::complex<double> &z, std::complex<double> &H_plus, std::complex<double> &dH_plus) const;
  void compute_F_functions(const std::complex<double> &z, std::complex<double> &F, std::complex<double> &dF) const;
  void compute_H_from_logarithmic_derivative(const std::complex<double> &z, const std::complex<double> &h, 
                                           int omega, std::complex<double> &H_plus, std::complex<double> &dH_plus) const;
  void compute_F_from_logarithmic_derivative(const std::complex<double> &z, const std::complex<double> &f,
                                           std::complex<double> &F, std::complex<double> &dF) const;
  double inf_norm(const std::complex<double> &z) const;
  
  // Recursion relations for real parameters when Michel solver fails
  double get_next_el_cwf(const double &rho, const double &eta, const int &eL, 
                         const double &u_L, const double &upr_L) const;
  double get_next_el_cwf_deriv(const double &rho, const double &eta, const int &eL, 
                               const double &u_L_minus, const double &u_L) const;
  void compute_asymptotic_bessel_solution(double rho_r, double eta_r, int eL,
                                          double &F, double &dF, double &G, double &dG) const;

  // Constants for numerical stability
  static constexpr double NUMERICAL_MAX_EXPONENT = 700.0; // log(max double)
};

#endif