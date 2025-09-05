#ifndef COULOMB_WF_COMPLEX_H
#define COULOMB_WF_COMPLEX_H

#include<complex>
#include<algorithm>
#include <vector>
#include<math.h>

#include "complex_functions.H"

// Define constants locally instead of using ScaleUtils
namespace ScaleUtils {
    namespace Math {
        const double PI = 3.14159265358979323846;
    }
}
// Direct integration of the Coulomb equation
// ------------------------------------------
// One uses the Burlisch-Stoer-Henrici method, where one integrates on different meshes
// with the Henrici method, and then uses the Richardson method to get the final result by extrapolation.
// Numerical Recipes, Chap. 16.4 .

class ODE_integration
{
public:
  ODE_integration(const std::complex<double>& l_1,
                  const std::complex<double>& two_eta_1)
      : l(l_1)
      , ll_plus_one(l_1 * (l_1 + 1.0))
      , two_eta(two_eta_1)
      , precision(1E-10)
      , m_tab(8)
      , one_over_m_tab(8)
      , interpolation_term_tab(8)
  {
      for (int n = 0; n < 8; n++)
      {
          interpolation_term_tab[n].resize(8);
          for (int i = 0; i < n; i++)
          {
              interpolation_term_tab[n][i] = 1.0;

              for (int j = 0; j < n; j++)
                  if (i != j)
                      interpolation_term_tab[n][i] *= (i + 1.0) / (i - j);
          }
      }
      for (unsigned int k = 0; k < 8; k++)
          m_tab[k] = 2 * (k + 1);
      for (unsigned int k = 0; k < 8; k++)
          one_over_m_tab[k] = 1.0 / static_cast<double>(m_tab[k]);
  }

  void operator() (const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
		   const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const;

  private:
  std::complex<double> extrapolation_in_zero (const unsigned int n,const std::complex<double> *f) const;
  std::complex<double> F_r_u (const std::complex<double> &r,const std::complex<double> &u) const;
  void integration_Henrici (const unsigned int m,const std::complex<double> &h,
			    const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
			    const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const;

  const std::complex<double> l,ll_plus_one;  // angular momentum,l(l+1).
  const std::complex<double> two_eta;        // 2.eta, with eta the Sommerfeld parameter.
  double precision;//note that this should match the value given in the Coulomb_wave_functions class
  std::vector<unsigned int> m_tab;     // integers used in the extrapolation
                                       // method.
  std::vector<double> one_over_m_tab;  // doubles used in the extrapolation
                                       // method.
  std::vector<std::vector<double>> interpolation_term_tab;
};






// Class to calculate the Coulomb wave functions
// ---------------------------------------------

class Coulomb_wave_functions
{
public:
  // Constructor.
  // ------------
  // Constants are defined in the constructor,
  // plus a pointer to class ODE_integration, ODE_ptr, to integrate numerically
  // the regular Coulomb wave function.
  //
  // Variables:
  // ----------
  // is_it_normalized_c : true if one wants normalized functions, i.e. the
  // standard normalization,
  //                      false if one wants F -> F/C(l,eta) and H+/H-/G ->
  //                      H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1
  //                      and |z| small.
  // l_c : orbital angular momentum.
  // eta_c : Sommerfeld parameter.

  Coulomb_wave_functions(const bool                  is_it_normalized_c,
                         const std::complex<double>& l_c,
                         const std::complex<double>& eta_c)
      : l(l_c)
      , eta(eta_c)
      , is_it_normalized(is_it_normalized_c)
      , precision(1E-10)
      , sqrt_precision(1E-5)
      , // these two are intimately linked with prec_first_order_expansion -
        // update it as well if these change
      neg_int_omega_one(
          (rint(real(l_c + std::complex<double>(-imag(eta_c), real(eta_c))))
           == l_c + std::complex<double>(-imag(eta_c), real(eta_c)))
          && (rint(real(1 + l_c
                        + std::complex<double>(-imag(eta_c), real(eta_c))))
              <= 0.0))
      , neg_int_omega_minus_one(
            (rint(real(l_c - std::complex<double>(-imag(eta_c), real(eta_c))))
             == l_c - std::complex<double>(-imag(eta_c), real(eta_c)))
            && (rint(real(1 + l_c
                          - std::complex<double>(-imag(eta_c), real(eta_c))))
                <= 0.0))
      , sigma_l(sigma_l_calc(l_c, eta_c))
      , log_Cl_eta(log_Cl_eta_calc(l_c, eta_c))
      , Cl_eta(exp(log_Cl_eta_calc(l_c, eta_c)))
      , exp_I_chi(exp_I_omega_chi_calc(1, l_c, eta_c))
      , exp_minus_I_chi(exp_I_omega_chi_calc(-1, l_c, eta_c))
      , one_over_sin_chi(1.0 / sin_chi_calc(l_c, eta_c))
      , log_cut_constant_CFa_plus(
            log_cut_constant_CFa_calc(is_it_normalized_c, 1, l_c, eta_c))
      , log_cut_constant_CFa_minus(
            log_cut_constant_CFa_calc(is_it_normalized_c, -1, l_c, eta_c))
      , cut_constant_CFa_plus(
            exp(log_cut_constant_CFa_calc(is_it_normalized_c, 1, l_c, eta_c)))
      , cut_constant_CFa_minus(
            exp(log_cut_constant_CFa_calc(is_it_normalized_c, -1, l_c, eta_c)))
      , log_cut_constant_CFb_plus(
            log_cut_constant_CFb_calc(is_it_normalized_c, 1, l_c, eta_c))
      , log_cut_constant_CFb_minus(
            log_cut_constant_CFb_calc(is_it_normalized_c, -1, l_c, eta_c))
      , log_cut_constant_AS_plus(log_cut_constant_AS_calc(1, l_c, eta_c))
      , log_cut_constant_AS_minus(log_cut_constant_AS_calc(-1, l_c, eta_c))
      , cut_constant_CFb_plus(
            exp(log_cut_constant_CFb_calc(is_it_normalized_c, 1, l_c, eta_c)))
      , cut_constant_CFb_minus(
            exp(log_cut_constant_CFb_calc(is_it_normalized_c, -1, l_c, eta_c)))
      , log_sym_constant_arg_neg(
            (is_it_normalized_c)
                ? (-ScaleUtils::Math::PI
                   * (eta_c + (l_c + 1) * std::complex<double>(0.0, 1.0)))
                : (-ScaleUtils::Math::PI * (l_c + 1)
                   * std::complex<double>(0.0, 1.0)))
      , log_sym_constant_arg_pos(
            (is_it_normalized_c)
                ? (-ScaleUtils::Math::PI
                   * (eta_c - (l_c + 1) * std::complex<double>(0.0, 1.0)))
                : (ScaleUtils::Math::PI * (l_c + 1)
                   * std::complex<double>(0.0, 1.0)))
      , sym_constant_arg_neg(
            (is_it_normalized_c)
                ? (exp(-ScaleUtils::Math::PI
                       * (eta_c + (l_c + 1) * std::complex<double>(0.0, 1.0))))
                : (exp(-ScaleUtils::Math::PI * (l_c + 1)
                       * std::complex<double>(0.0, 1.0))))
      , sym_constant_arg_pos(
            (is_it_normalized_c)
                ? (exp(-ScaleUtils::Math::PI
                       * (eta_c - (l_c + 1) * std::complex<double>(0.0, 1.0))))
                : (exp(ScaleUtils::Math::PI * (l_c + 1)
                       * std::complex<double>(0.0, 1.0))))
      , turning_point(std::max(
            1.0,
            std::abs(eta_c)
                + sqrt(std::abs(l_c * (l_c + 1.0)) + std::abs(eta_c * eta_c))))
      , is_H_dir_int_naive(false)
      , cwf_real_ptr(0)
      , cwf_real_eta_plus_ptr(0)
      , cwf_real_eta_minus_ptr(0)
      , cwf_real_l_plus_ptr(0)
      , cwf_real_l_minus_ptr(0)
      , cwf_minus_eta_ptr(0)
      , cwf_lp_ptr(0)
      , prec_first_order_expansion(0.1 * sqrt_precision)
  {
    ODE_ptr = new class ODE_integration (l,2.0*eta);

    debut = 0.0;

    if (real (l) >= 0.0)
    {
      F_debut = 0.0;
      dF_debut = (l == 0) ? ((is_it_normalized) ? (Cl_eta) : (1.0)) : (0.0);
    }
  }

  ~Coulomb_wave_functions (void)
  {
    delete cwf_real_ptr;

    delete cwf_real_l_plus_ptr;
    delete cwf_real_l_minus_ptr;

    delete cwf_real_eta_plus_ptr;
    delete cwf_real_eta_minus_ptr;

    delete cwf_minus_eta_ptr;
    delete cwf_lp_ptr;

    delete ODE_ptr;
  }

  void F_dF_init (const std::complex<double> &z,const std::complex<double> &F,const std::complex<double> &dF);

  void F_dF(const std::complex<double>& z,
            std::complex<double>&       F,
            std::complex<double>&       dF);
  void G_dG (const std::complex<double> &z,std::complex<double> &G,std::complex<double> &dG);
  void H_dH (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH);
  void H_dH_scaled (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH);

  const std::complex<double> l,eta; // Angular momentum and Sommerfeld parameter.

  const bool is_it_normalized;
  // true if F(z) ~ C(l,eta).z^{l+1} in 0, false if F(z) ~ z^{l+1} in 0.

private:

  void asymptotic_series (const int omega,const std::complex<double> &one_over_z,std::complex<double> sum[],std::complex<double> dsum[],bool &is_it_successful);

  std::complex<double> continued_fraction_f (const std::complex<double> &z,const int omega);
  std::complex<double> continued_fraction_h (const std::complex<double> &z,const int omega);

  void F_dF_power_series (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF);

  void asymptotic_expansion_F_dF (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF,bool &is_it_successful);
  void asymptotic_expansion_H_dH_scaled (const int omega,const std::complex<double> &one_over_z,
					 std::complex<double> &H_scaled,std::complex<double> &dH_scaled,bool &is_it_successful);

  void F_dF_direct_integration (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF,bool &is_it_successful);
  void H_dH_direct_integration (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH);

  void partial_derivatives (const bool is_it_regular,const bool is_it_eta,const double x,double &d_chi_Bx,double &d_chi_dBx);
  void first_order_expansions (const bool is_it_regular,const std::complex<double> &z,std::complex<double> &B,std::complex<double> &dB);
  void H_dH_from_first_order_expansions (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH);

  void H_dH_with_F_dF_and_CF (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH);
  void H_dH_with_expansion (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH,bool &is_it_successful);

  void F_dF_with_symmetry_relations (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF);

  const double precision;
  const double sqrt_precision;

  const bool neg_int_omega_one,neg_int_omega_minus_one;
  // neg_int_omega_one : true if 1+l+i.eta is negative integer, false if not.
  // neg_int_omega_minus_one : true if 1+l-i.eta is negative integer, false if not.

  const std::complex<double> sigma_l,log_Cl_eta,Cl_eta; // sigma(l,eta), log[C(l,eta)], C(l,eta)

  const std::complex<double> exp_I_chi, exp_minus_I_chi, one_over_sin_chi;
  // exp (i.chi), exp (-i.chi), 1/sin (chi) with chi = sigma(l,eta) - sigma(-l-1,eta) - (l+1/2).Pi .
  // They are used to calculate H+/H- from F(l,eta,z) and F(-l-1,eta,z) if |Im[l]| >= 1 and |z| <= 1 .

  const std::complex<double> log_cut_constant_CFa_plus,log_cut_constant_CFa_minus;
  const std::complex<double> cut_constant_CFa_plus,cut_constant_CFa_minus;
  const std::complex<double> log_cut_constant_CFb_plus,log_cut_constant_CFb_minus;
  const std::complex<double> log_cut_constant_AS_plus,log_cut_constant_AS_minus;
  const std::complex<double> cut_constant_CFb_plus,cut_constant_CFb_minus;
  // cut constants and their logs for continued fractions (CFa and CFb) and asymptotic series (AS).
  // plus,minus is for omega = 1 or -1.
  // See functions log_cut_constant_AS_calc, log_cut_constant_CFa_calc and log_cut_constant_CFb_calc.


  const std::complex<double> log_sym_constant_arg_neg,log_sym_constant_arg_pos,sym_constant_arg_neg,sym_constant_arg_pos;
  // Multiplicative constants and their logs used in the following reflection
  // formulas : F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta-i.l)] if arg (z) > 0
  // and is_it_normalized is true, so sym_constant_arg_pos =
  // -exp[-Pi.(eta-i.l)], F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta+i.l)] if arg
  // (z) <= 0 and is_it_normalized is true, so sym_constant_arg_neg =
  // -exp[-Pi.(eta-i.l)], F(l,eta,z) = -F(l,-eta,-z).exp[i.Pi.l)] if arg (z) >
  // 0 and is_it_normalized is false, so sym_constant_arg_pos = -exp[i.Pi.l)],
  // F(l,eta,z) = -F(l,-eta,-z).exp[-i.Pi.l)] if arg (z) <= 0 and
  // is_it_normalized is false, so sym_constant_arg_neg = -exp[-i.Pi.l)].

  const double turning_point; // turning_point : max (1,||eta| + sqrt[|l(l+1)| + |eta|^2]|).

  bool is_H_dir_int_naive; // true if one integrates H+/H- forward without considering |H+/H-|, false if not. It is false except in continued_fraction_h.

  std::complex<double> debut,F_debut,dF_debut;
  // Coulomb wave functions and derivative at z = debut.
  // It is used to integrate the Coulomb wave function faster,
  // as debut is usually close to the argument of the Coulomb wave function so
  // that the integration is quicker and more stable.

  class ODE_integration *ODE_ptr;  // pointer to class ODE_integration to integrate numerically the Coulomb equation.

  class Coulomb_wave_functions *cwf_real_ptr,*cwf_real_eta_plus_ptr,*cwf_real_eta_minus_ptr,*cwf_real_l_plus_ptr,*cwf_real_l_minus_ptr;
  // pointers to classes Coulomb_wave_functions of parameters (l_r,eta_r) (one
  // has eta_r = Re[eta], eta_i = Im[eta], l_r = Re[l] and l_i = Im[l]), (l_r
  // +/- epsilon[l],eta_r) and (l_r,eta_r +/- epsilon[eta]). They are first put
  // to zero and allocated in the program when they are needed. They are used
  // for the first order expansion method when |l_i| << 1, |eta_i| << 1 and
  // |Im[z]| << Re[z] with Re[z] > 0.

  class Coulomb_wave_functions *cwf_minus_eta_ptr,*cwf_lp_ptr;
  // pointers to classes Coulomb_wave_functions of parameters (l,-eta), (lp = -l-1,eta) and (l_r +/- precision,eta).
  // They are first put to zero and allocated in the program when they are needed.
  // They are used for symmetry relations : F(l,eta,z) \propto F(l,-eta,-z) and h[omega](l,eta,z) = -h[omega](l,-eta,-z)
  // and to calculate H+/H- from F(l,eta,z) and F(lp,eta,z) if |Im[l]| >= 1 and |z| <= 1.

  const double prec_first_order_expansion; // prec_first_order_expansion : 0.1*sqrt_precision. It is the precision used for first_order_expansions.
};


#endif //end of COULOMB_WF_COMPLEX_H
