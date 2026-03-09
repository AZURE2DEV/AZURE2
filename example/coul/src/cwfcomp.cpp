#include "cwfcomp.H"

// Calculation of h(omega) = H(omega)'/H(omega) with a continued fraction.
// -----------------------------------------------------------------------
// One calculates the ratio h = H'/H with the continued fraction of the associated hypergeometric confluent function.
// One uses Lentz's method.
// One has : h = [b[0] + a[1]/b[1]+ a[2]/b[2]+ ... a[n]/b[n]+ ...].i.omega/z with :
// b[0] = z - eta, a[n] = (1 + l + i.omega + n-1).(i.omega.eta - l + n-1), b[n] = 2[z - eta] + i.omega.n .
//
// If the number of iterations reaches 100000 and |z| > 0.5, the convergence is too slow. One is probably very close to the imaginary axis.
// If l=0 and |z| <= 0.5, the direct integration is still done as it is stable.
// If 1+l+/-i.eta is negative integer, it has to be done otherwise it is too long even if |z| is not exceedingly small.
// One first considers Re[z] >= 0.
// One takes the starting point z0 = 2 + i.(Im[z] + 2.sign[Im[z]]) (0.6 + i.sign[Re[z]].0.6.sign[Im[z]] if |z| <= 0.5).
// Then, one calculates H[omega],H[omega]' at the starting point, and one integrates numerically H[omega] until z.
// h(omega)(z) is then H[omega](z)'/H[omega](z).
// If Re[z] < 0, one calculates H[-omega],H[-omega]' with l, eta -> -eta, and z -> -z.
// One uses the Coulomb wave functions class cwf_minus_eta_ptr defined with l and -eta. 
// The ratio h(omega,l,eta,z) is then equal to -h(-omega,l,-eta,-z).
// To avoid infinite loops, continued_fraction_h must not be used in this integration. So, the starting point is chosen so |H[+/-]| is very likely
// to increase in modulus. Then, one always integrates forward. Forward integration is enforced putting
// is_H_dir_int_naive to true. It is put to false again at the end of the calculation.
//
// Variables:
// ----------
// z : variable of the Coulomb wave function.
// omega : 1 for the outgoing wave function ratio H+'/H+, -1 for the incoming wave function ratio H-'/H-.
// large,small : 1E50,1E-50. They are used in the case of vanishing denominators or numerators.
// I_omega,I_omega_eta,two_I_omega,z_minus_eta,two_z_minus_eta : i.omega,i.omega.eta, 2.i.omega, z-eta, 2(z-eta).
// a,c : 1 + l + i.omega.eta, i.omega.eta - l.
// b0,Cn,Dn,an,bn,Delta_n : variables used in the Lentz method.
// n,nm1 : index of a[n] and b[n], n-1 
// bn_plus_an_Dn : bn + an.Dn. Dn = 1/[bn + an.Dn] or 1E50 if bn + an.Dn is zero.
// bn_plus_an_over_Cn : bn + an/Cn. Cn = bn + an/Dn or 1E-50 if bn + an/Dn is zero.
// hn : value of the continuous fraction during the iteration process.
// test : test of convergence of hn.
// cwf : reference to *this if Re[z] >= 0, reference to cwf_minus_eta_ptr if Re[z] < 0
//       cwf_minus_eta_ptr is allocated first if it is zero, in the case of Re[z] < 0 .
//       It is used to integrate the Coulomb wave functions with l,eta for Re[z] >= 0 or l,-eta for Re[z] < 0. 
// z00 : 2 + i.sign[Re[z]].(Im[z] + 2.sign[Im[z]]).
// z01 : 0.6 + i.0.6.sign[Re[z]].sign[Im[z]].
// abs_z : |z|
// z_start,F_start,dF_start : starting point of the direct integration, F(l,+/-eta,z_start),F'(l,+/-eta,z_start) 
// H,dH : Coulomb wave function and derivative in l,eta,omega,z if Re[z] > 0, in l,-eta,-omega,-z if Re[z] < 0.
// h : value of H(omega)'/H(omega).
// debut_cwf,F_debut_cwf,dF_debut_cwf : values stored in cwf put back in cwf at the end of direct integration as they change values.

std::complex<double> Coulomb_wave_functions::continued_fraction_h (const std::complex<double> &z,const int omega)
{ 
  const double small = 1E-50,large = 1E50,abs_z = abs (z);
  const std::complex<double> I_omega(0.0,omega),two_I_omega(0.0,2.0*omega),I_omega_eta = I_omega*eta;
  const std::complex<double> a = I_omega_eta + l + 1.0,c = I_omega_eta - l,z_minus_eta = z - eta,two_z_minus_eta = 2.0*z_minus_eta;

  std::complex<double> b0 = z_minus_eta,hn = (b0 != 0.0) ? (b0) : (1E-50), Cn = hn, Dn = 0.0;
  int n = 1;
  double test;
  do
  {
    const int nm1 = n-1;
    const std::complex<double> an = (a + nm1)*(c + nm1),bn = two_z_minus_eta + n*two_I_omega,bn_plus_an_Dn = bn + an*Dn,bn_plus_an_over_Cn = bn + an/Cn;

    Dn = (bn_plus_an_Dn != 0.0) ? (1.0/bn_plus_an_Dn) : (large);
    Cn = (bn_plus_an_over_Cn != 0.0) ? (bn_plus_an_over_Cn) : (small);

    const std::complex<double> Delta_n = Dn*Cn;
    hn *= Delta_n;
    test = inf_norm (1.0 - Delta_n);
     
    if ((n++ > 100000) && ((l == 0.0) || (abs_z > 0.5) || neg_int_omega_one || neg_int_omega_minus_one))
    {
      if ((real (z) < 0.0) && (cwf_minus_eta_ptr == 0)) cwf_minus_eta_ptr = new class Coulomb_wave_functions (is_it_normalized,l,-eta);
      class Coulomb_wave_functions &cwf = (real (z) < 0.0) ? (*cwf_minus_eta_ptr) : (*this);

      const std::complex<double> z00(2.0,SIGN (real (z))*(imag (z) + 2.0*SIGN (imag (z)))),z01(0.6,0.6*SIGN (real (z))*SIGN (imag (z)));
      const std::complex<double> z_start = (abs_z > 0.5) ? (z00) : (z01),debut_cwf = cwf.debut,F_debut_cwf = cwf.F_debut,dF_debut_cwf = cwf.dF_debut;
      std::complex<double> F_start,dF_start,H,dH;
      cwf.F_dF (z_start,F_start,dF_start);
      cwf.is_H_dir_int_naive = true, cwf.H_dH_direct_integration (SIGN (real (z))*omega,SIGN (real (z))*z,H,dH), cwf.is_H_dir_int_naive = false;
      cwf.debut = debut_cwf, cwf.F_debut = F_debut_cwf, cwf.dF_debut = dF_debut_cwf;
      return (SIGN (real (z))*dH/H);
    }
  }
  while (test > 1E-15);

  const std::complex<double> h = hn*I_omega/z;
  return h;
}










// Calculation of H(omega) and dH(omega)/dz (scaled) with asymptotic series
// ------------------------------------------------------------------------
// H[omega](z) = exp[i.omega.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S[omega](z) for Re[z] >= 0.
// H[omega]'(z) = exp[i.omega.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S[omega]'(z) + i.omega.(1 - eta/z).S[omega](z)] for Re[z] >= 0.
//
// S[omega](z) is the asymptotic series and S[omega]'(z) its derivative calculated in asymptotic_series.
// If they did not converge, one leaves the routine.
//
// The negative cut is taken into account if log [cut_constant_AS] is finite, i.e. cut_constant_AS is not exactly zero : 
//
// If Re[z] < 0.0 and omega.Im[z] < 0.0 :
// --------------------------------------
// H[omega] = H[omega][ASd] + cut_constant_AS_plus.H[-omega][ASd].
// H[omega][ASd] is given by directly by the asymptotic series.
// H[-omega][ASd] is given by the asymptotic series, which gives the good result H[-omega] in this case as one is not in its bad quadrant.
//
// The function is scaled, so one returns : H[omega](z).exp[-i.omega.[z - eta.log[2z]]] and H[omega]'(z).exp[-i.omega.[z - eta.log[2z]]] .
//
// In the case of overflows or underflows, one uses logs.
//
//
// Variables :
// -----------
// omega : 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
// one_over_z : 1/z.
// z : variable of the Coulomb wave function.
// H_scaled,dH_scaled : H[omega](z).exp(-i.omega.[z - eta.log[2z]]) and H[omega]'(z).exp(-iomega.[z - eta.log[2z]]).
// is_it_successful : true if the asymptotic expansions converged, false it not
// sum,dsum : {S[omega](z),S[-omega](z)} and {S[omega]'(z),S[-omega]'(z)}
// I_omega,two_I_omega : i.omega, 2i.omega.
// I_omega_one_minus_eta_over_z : i.omega.(1 - eta/z).
// phase_shift : i.omega.[-l.Pi/2 + sigma(l,eta)]
// exp_phase_shift : exp (i.omega.(-l.Pi/2 + sigma(l,eta))).
// log_cut_constant_AS : log_cut_constant_AS_plus if omega = 1, log_cut_constant_AS_minus if omega =-1.
//                       If log_cut_constant_AS is not finite (i.e. cut_constant_AS = 0 if it is defined), there is no branch cut to consider.
// factor : exp (log_cut_constant_AS - 2.i.omega.(z - eta.log(2z)) - phase_shift).

void Coulomb_wave_functions::asymptotic_expansion_H_dH_scaled (const int omega,const std::complex<double> &one_over_z,
							       std::complex<double> &H_scaled,std::complex<double> &dH_scaled,bool &is_it_successful)
{  
  std::complex<double> sum[2],dsum[2];

  asymptotic_series (omega,one_over_z,sum,dsum,is_it_successful);
  if (!is_it_successful) return;

  const std::complex<double> I_omega(0,omega),two_I_omega(0,2*omega),I_omega_one_minus_eta_over_z = I_omega*(1.0 - eta*one_over_z);

  const std::complex<double> phase_shift = I_omega*(sigma_l - l*M_PI_2),exp_phase_shift = exp (phase_shift);

  H_scaled = sum[0]*exp_phase_shift;
  dH_scaled = (dsum[0] + sum[0]*I_omega_one_minus_eta_over_z)*exp_phase_shift;

  const std::complex<double> log_cut_constant_AS = (omega == 1) ? (log_cut_constant_AS_plus) : (log_cut_constant_AS_minus);

  if (one_over_z != 0.0)
  {
    const std::complex<double> z = 1.0/one_over_z;
 
    if (isfinite (log_cut_constant_AS) && (real (z) < 0.0) && (SIGN (imag (z)) == -omega))
    {
      const std::complex<double> factor = exp (-two_I_omega*(z - eta*(M_LN2 + log (z))) - phase_shift + log_cut_constant_AS);

      H_scaled += sum[1]*factor;
      dH_scaled += (dsum[1] - sum[1]*I_omega_one_minus_eta_over_z)*factor;
    }
  }

  if (!is_it_normalized)
  {
    if ((Cl_eta == 0.0) || (!isfinite (Cl_eta))) 
      H_scaled = exp (log_Cl_eta + log (H_scaled)),dH_scaled = exp (log_Cl_eta + log (dH_scaled));
    else 
      H_scaled *= Cl_eta,dH_scaled *= Cl_eta;
  }
}






// Calculation of H[omega](z) and H[omega]'(z) by direct integration.
// ------------------------------------------------------------------
// To calculate H[omega](z) and H'[omega](z), one integrates numerically H[omega]''(z) = [l(l+1)/z^2 + 2.eta/z - 1].H[omega](z)
// starting from debut,H_debut = H[omega](debut) and dH_debut = H'[omega](debut).
// There is no branch cut problem as Re[z] >= 0.
// The starting point comes for the regular function from the stored values debut, F_debut and dF_debut.
// If debut = 0, one puts debut = debut_omega = z/|2z| and calculates F(debut) and F'(debut) with power series.
// Then, the starting point {debut_omega,H[omega](debut),H'[omega](debut)} is calculated 
// from {debut,F(debut),F'(debut)} and the continued fraction h[omega](debut).
// The first order expansions method is used if one is very close to the real axes of l,eta and z (Re[z] > 0).
// The step of the integration is (z - debut)/N_num, with N_num = [|z - debut|/min (0.1,10.turning_point)] + 1.
// The value of min (0.1,10/turning_point) gives a smaller step when turning_point increases, 
// as calculations become there more difficult as |H[omega]| typically varies faster in this case.
// The intermediates points are called z_aft. They go from debut to z, 
// and (debut_omega,H_debut,dH_debut) is put to {z_aft,H[omega](z_aft),H'[omega](z_aft)} at each step.
// If |H[omega]| increases along the path, the integration is stable.
// If is_H_dir_int_naive is true, one has to integrate forward, as this integration is used to calculate the continued fraction.
// If not, and if |H[omega]| decreases, 
// one integrates backwards from z_aft to debut_omega with the knowledge of h[omega](z_aft) = H'[omega](z_aft)/H[omega](z_aft).
// H'[omega](z_aft)/H[omega](z_aft) is given by the continued fraction formula.
// One then obtains by direct integration C.H[omega](debut) and C.H'[omega](debut).
// Knowing H_debut, one deduces H[omega](z) = 1/C and H'[omega](z) = f(z_aft).H[omega](z).
// Increase or decrease is known using the Taylor expansion of H[omega] near debut_omega in z up to second order.
//
// If H(z_aft) is not finite, one stops the integration.
//
// Variables
// ---------
// omega : 1 for H+,H+', -1 for H-,H-'.
// z : variable of the Coulomb wave function.
// H,dH : wave function H[omega] and derivative H'[omega] to calculate.
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// debut_omega, H_debut,dH_debut : starting point or the integration. debut_omega = debut at the beginning.
// ODE : reference to pointer ODE_ptr, which performs direct integration from one point to another.
// step_abs : length of the integration step from z to debut. It is min (0.1,10/turning_point)
// N_num : number of integrations to do from debut to z. It is |z - debut|/step_abs + 1.
// step_num : std::complex step from debut to z. It is (z - debut)/N_num .
// ll_plus_one,two_eta : l(l+1), 2.eta. They are used in the Taylor expansion test.
// z_aft : next point in which H[omega] is calculated by direct integration . It is z - i*step_num, with i from N_num-1 to 0.
// one_over_debut,log_H_debut_der,d2H_debut_over_H_debut: 1/debut, H'[omega](debut)/H[omega](debut), H[omega]''(debut)/H[omega](debut). 
//                                                        They are used in the Taylor expansion test. 
// h : continued fraction h[omega] = H[omega]'/H[omega] in z.
// H_debut_not_normed,dH_debut_not_normed : unnormed regular wave function and derivative in debut_omega,
//                                          after integration from the starting point {z_aft, 1.0, H'[omega](z_aft)/H[omega](z_aft)}.
//                                          Then, H[omega](z_aft) =  H_debut/H_debut_not_normed and H'[omega](z_aft) = h[omega].H[omega](z_aft) .

void Coulomb_wave_functions::H_dH_direct_integration (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH)
{ 
  const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);

  if (debut == 0.0) 
  {
    debut = 0.5*z/abs (z);
    F_dF_power_series (debut,F_debut,dF_debut);
  }

  std::complex<double> debut_omega = debut,H_debut,dH_debut;

  if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
      && (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
      && (!neg_int_omega_one && !neg_int_omega_minus_one))
    first_order_expansions (omega,debut_omega,H_debut,dH_debut);
  else
    H_dH_with_F_dF_and_CF (omega,debut_omega,H_debut,dH_debut);

  const class ODE_integration &ODE = *ODE_ptr;
  const double step_abs = std::min(0.1,10.0/turning_point);
  const unsigned int N_num = static_cast<unsigned int> ((abs (z-debut_omega)/step_abs) + 1);
  const std::complex<double> step_num = (z - debut_omega)/static_cast<double> (N_num),ll_plus_one = l*(l+1.0),two_eta = 2.0*eta;

  for (unsigned int i = N_num-1 ; i <= N_num ; i--)
  {
    const std::complex<double> z_aft = z - i*step_num,one_over_debut = 1.0/debut,log_H_debut_der = dH_debut/H_debut;
    const std::complex<double> d2H_debut_over_H_debut = (ll_plus_one*one_over_debut + two_eta)*one_over_debut - 1.0;
  
    if (is_H_dir_int_naive)
      ODE (debut_omega,H_debut,dH_debut,z_aft,H,dH);
    else if (abs (1.0 + step_num*(log_H_debut_der + 0.5*step_num*d2H_debut_over_H_debut)) < 1.0)
    {
      const std::complex<double> h = continued_fraction_h (z_aft,omega);      
      std::complex<double> H_debut_not_normed,dH_debut_not_normed;

      ODE (z_aft,1.0,h,debut_omega,H_debut_not_normed,dH_debut_not_normed);
      H = H_debut/H_debut_not_normed;
      dH = h*H; 
    }
    else ODE (debut_omega,H_debut,dH_debut,z_aft,H,dH);

    if (!isfinite (H) || !isfinite (dH)) std::cout<<"Numerical failure encountered in H_dH_direct_integration."<<std::endl,exit (1);

    debut_omega = z_aft,H_debut = H,dH_debut = dH;
  }
}







// Calculation of H[omega](z), H'[omega](z) with the first order expansion method.
// -------------------------------------------------------------------------------
// When imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
// one has to separate the calculation of the real and imaginary parts of H[omega] and H'[omega], as they can differ by tens of orders of magnitude.
//
// For that, one expands F(z),G(z),F'(z),G'(z) in first order in y, eta_i and l_i in first_order_expansions.
//
// Then, H[omega](z) = G(z) + i.omega.norm.F(z) and H'[omega](z) = G'(z) + i.omega.norm.F'(z),
// with norm = 1 if the wave functions are normalized and C(l,eta)^2 if not.
// One uses logs if norm underflows or overflows.
//
// Variables:
// ----------
// omega : 1 if one calculates H+(z),H+'(z), -1 if one calculates H-(z),H-'(z).
// z : variable of the Coulomb wave function.
// H,dH : wave function H[omega] and derivative H'[omega] to calculate.
// F,dF,G,dG : regular and irregular wave functions and derivatives in z, calculated with first order expansions.
// I_omega : i.omega
// norm_functions,log_norm : 1 if the wave functions are normalized and C(l,eta)^2 if not, its log.


void Coulomb_wave_functions::H_dH_from_first_order_expansions (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH)
{
  std::complex<double> F,dF,G,dG;
  first_order_expansions (true,z,F,dF);
  first_order_expansions (false,z,G,dG);

  const std::complex<double> I_omega(0,omega),norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);

  if ((norm_functions == 0.0) || (!isfinite (norm_functions))) 
  {
    const std::complex<double> log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta) : (0.0);

    H = G + I_omega*exp (log (F) + log_norm);
    dH = dG + I_omega*exp (log (dF) + log_norm);
  }
  else
  {
    H = G + I_omega*norm_functions*F;
    dH = dG + I_omega*norm_functions*dF;
  }
}













// Calculation of H(omega) and H(omega)' with F(z), F'(z) and continued fractions.
// -------------------------------------------------------------------------------
// If 1+l-i.omega.eta is negative integer, one uses H[omega] = 1/F/(f - h[omega]) as it is the only available solution
// linearly independent of F.
// If 1+l+i.omega.eta is negative integer, H[omega] is proportional to F so the values F and F' are arbitrarily chosen for H[omega] and H'[omega].
// 
// If not, one calculates h[sign(Im[z])](z), and one has f = F'(z)/F(z).
// One chooses sign(Im[z]), as h converges fastest for z in this region.
// If |f - h[sign(Im[z])]|oo >= 1, F and H[sign] are numerically linearly independent so the continued fraction is meaningful.
// Then, h[-sign(Im[z])](z) is not needed and is put to f, the worst value it can have. If not, h[-sign(Im[z])](z) is needed and calculated.
// Then, h[omega] and h[-omega] take their values from h[sign(Im[z])]| and h[-sign(Im[z])]|.
// 
// If |f - h[omega]|oo > |f - h[-omega]|oo, one uses the continued fraction h[omega].
// If Re[z] > 0 or sign(Im[z]) = omega (good quadrants), H[omega] = 1/F/(f - h[omega]) and H'[omega] = h[omega].H[omega].
// If not, one has to take the branch cut into account, as one is in the bad quadrant of H[omega] : 
// H[omega] = 1/F/(f - h[omega]) - cut_constant.F and H'[omega] = h[omega].H[omega] - cut_constant.F' .
// cut_constant is cut_constant_CFa_plus if omega = 1, and cut_constant_CFa_minus if omega = -1.
// If log [cut_constant] is not finite, it means that cut_constant is exactly zero if it is defined, so that there is no branch cut to consider.
//
// If |f - h[omega]|oo < |f - h[-omega]|oo, one uses the continued fraction h[-omega],
// calculating H[omega] and H'[omega] using H[omega] = H[-omega] + constant.F .
// constant is 2.i.omega.norm if Re[z] > 0 or sign(Im[z]) = -omega (good quadrants), and cut_constant if not.
// cut_constant is cut_constant_CFb_plus if omega = 1, cut_constant_CFb_minus if omega = -1.
// norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
//
// If cut_constant underflows or overflows, one uses logs of F,F' and log_cut_constant for the calculation.
//
// Variables:
// ----------
// omega : 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
// z : variable of the Coulomb wave function.
// F,dF: regular Coulomb wave function and derivative.
// H,dH : H+(z) and H+'(z) if omega=1, H-(z) and H-'(z) if omega=-1.
// x,y : Re[z], Im[z].
// f,two_I_omega : ratio F'(z)/F(z), 2.i.omega .
// h_sign,h_minus_sign : continuous fractions h[-SIGN[y]] and h[SIGN[y]]. See before for their calculations
// h_omega,h_minus_omega : continuous fraction h[omega] and h[-omega].
// cut_constant_CFa,log_cut_constant_CFa : cut constant for H[omega], so in its bad quadrant H[omega](z) =  1/F/(f - h[omega]) - cut_constant.F and its log.
//                                         cut_constant is cut_constant_CFa_plus if omega = 1, cut_constant_CFa_minus if omega = -1.
//                                         If the log is not finite, cut_constant has to be strictly zero if it is defined so there is no branch cut to consider.
// H_minus_omega,dH_minus_omega : H-(z) and H-'(z) if omega=1, H+(z) and H+'(z) if omega=-1 in good quadrants. 
//                                In bad quadrants, they are not equal to the previous functions, but combined with branch cut formulas, one calculates them.
// cut_constant_CFb,constant,norm_functions : cut_constant is cut_constant_CFb_plus if omega = 1, cut_constant_CFb_minus if omega = -1.
//                                            constant is 2.i.omega.norm_functions in good quadrants, cut_constant in bad quadrants.
//                                            norm_functions is 1 for normalized functions, C(l,eta)^2 for unnormalized wave functions.
//                                            One has : H[-omega](z) =  1/F/(f - h[-omega]) + constant.F in all quadrants.
// log_cut_constant_CFb,log_constant,log_norm : logs of previous values.
// log_two_I_omega : log[2.i.omega] = log[2] + i.omega.Pi/2 .

void Coulomb_wave_functions::H_dH_with_F_dF_and_CF (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH)
{  
  std::complex<double> F,dF;
  F_dF (z,F,dF);

  const double x = real (z),y = imag (z);
  const std::complex<double> f = dF/F,two_I_omega(0.0,2.0*omega);

  if (((neg_int_omega_one && (omega == -1))) || ((neg_int_omega_minus_one && (omega == 1))))
  {
    const std::complex<double> h_omega = continued_fraction_h (z,omega);
    H = 1.0/(F*(f - h_omega)),dH =  H*h_omega;
  }
  else if (neg_int_omega_one || neg_int_omega_minus_one)
    H = F,dH = dF;
  else
  {
    const std::complex<double> h_sign = continued_fraction_h (z,SIGN (y)),h_minus_sign = (abs (f - h_sign) < 1.0) ? (continued_fraction_h (z,-SIGN (y))) : (f);
    const std::complex<double> h_omega = (omega == SIGN (y)) ? (h_sign) : (h_minus_sign),h_minus_omega = (omega == SIGN (y)) ? (h_minus_sign) : (h_sign);

    if (inf_norm (f - h_omega) > inf_norm (f - h_minus_omega))
    {
      H = 1.0/(F*(f - h_omega)),dH =  H*h_omega;

      const std::complex<double> log_cut_constant_CFa = (omega == 1) ? (log_cut_constant_CFa_plus) : (log_cut_constant_CFa_minus);

      if ((isfinite (log_cut_constant_CFa)) && (x < 0.0) && (SIGN (y) == -omega))
      {
	const std::complex<double> cut_constant_CFa = (omega == 1) ? (cut_constant_CFa_plus) : (cut_constant_CFa_minus); 
      
	if ((cut_constant_CFa == 0.0) || (!isfinite (cut_constant_CFa)))
	  H -= exp (log_cut_constant_CFa + log (F)),dH -= exp (log_cut_constant_CFa + log (dF));
	else 
	  H -= cut_constant_CFa*F,dH -= cut_constant_CFa*dF;
      }
    }
    else
    {
      const std::complex<double> H_minus_omega = 1.0/(F*(f - h_minus_omega)),dH_minus_omega = H_minus_omega*h_minus_omega;
      const std::complex<double> norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);
      const std::complex<double> cut_constant_CFb = (omega == 1) ? (cut_constant_CFb_plus) : (cut_constant_CFb_minus); 
      const std::complex<double> constant = ((x < 0.0) && (SIGN (y) == omega)) ? (cut_constant_CFb) : (two_I_omega*norm_functions);

      if ((constant == 0.0) || (!isfinite (constant))) 
      {
	const std::complex<double> log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta) : (0.0),log_two_I_omega(M_LN2,omega*M_PI_2);
	const std::complex<double> log_cut_constant_CFb = (omega == 1) ? (log_cut_constant_CFb_plus) : (log_cut_constant_CFb_minus);
	const std::complex<double> log_constant = ((x < 0.0) && (SIGN (y) == omega)) ? (log_cut_constant_CFb) : (log_two_I_omega + log_norm);
	
	H = exp (log_constant + log (F)) + H_minus_omega,dH = exp (log_constant + log (dF)) + dH_minus_omega;
      }
      else H = constant*F + H_minus_omega,dH = constant*dF + dH_minus_omega;
    }
  }
}












// Calculation of H[omega],H'[omega] for std::complex l with the expansion formula.
// ---------------------------------------------------------------------------
// When 2l is non-integer, one can expand H[omega] with F[l,eta,z] and F[-l-1,eta,z].
// H[omega] = (exp[i.omega.chi].F - Fp)/sin (chi), H'[omega] = (exp[i.omega.chi].F' - Fp')/sin (chi) if wave functions are normalized,
// H[omega] = (exp[i.omega.chi].F.C(l,eta)^2/sin (chi) + Fp/(2l+1), H'[omega] = (exp[i.omega.chi].F'.C(l,eta)^2/sin (chi) + Fp'/(2l+1) if not.
// chi is sigma(l,eta) - sigma(-l-1,eta) - (l+0.5).Pi, and Fp is F(-l-1,eta).
// Fp is calculated using a class Coulomb_wave_functions with parameters -l-1 and eta.
// To avoid numerical imprecisions, sin (chi) is calculated with the stable formula -(2l+1).C(l,eta).C(-l-1,eta).
// The validity of this expansion is checked with the wronskian of F and Fp, which must be correct up to precision :
// F'.Fp - Fp'.F = sin (chi) if wave functions are normalized, Fp'.F - F'.Fp = 2l + 1 if not.
// If the wronskians are numerically correct, one does the calculation and is_it_successful is put to true.
// If not, one puts is_it_successful to false and quits the routine.
// If C(l,eta)^2 underflows or overflows, one uses logs of F,F' and C(l,eta)^2 for the calculation.
//
// Variables
// ---------
// omega : 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
// z : variable of the Coulomb wave function.
// F,dF: regular Coulomb wave function and derivative in l,eta and z.
// H,dH : H+(z) and H+'(z) if omega=1, H-(z) and H-'(z) if omega=-1.
// is_it_successful : false if the wronskian between F(l,eta,z) and F(-l-1,eta,z) is not equal to zero up to precision, true if not.
// Fp,dFp : F(-l-1,eta,z), F'(-l-1,eta,z).
// exp_I_omega_chi : exp[i.omega.chi]
// one_over_2lp1 : 1/(2l + 1)
// Cl_eta_2,exp_I_omega_chi_over_sin_chi : C(l,eta)^2, exp[i.omega.chi]/sin (chi)
// F_Cl_eta_2, dF_Cl_eta_2 : F(z).C(l,eta)^2, F'(z).C(l,eta)^2.

void Coulomb_wave_functions::H_dH_with_expansion (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH,bool &is_it_successful)
{
  if (cwf_lp_ptr == 0) cwf_lp_ptr = new class Coulomb_wave_functions (is_it_normalized,-l-1,eta);

  std::complex<double> F,dF,Fp,dFp;
  F_dF (z,F,dF);
  cwf_lp_ptr->F_dF (z,Fp,dFp);
    
  const std::complex<double> exp_I_omega_chi = (omega == 1) ? (exp_I_chi) : (exp_minus_I_chi);
    
  if (is_it_normalized)
  { 
    if (inf_norm ((dFp*F - dF*Fp)*one_over_sin_chi - 1.0) < precision)
    {
      H = (exp_I_omega_chi*F - Fp)*one_over_sin_chi;
      dH = (exp_I_omega_chi*dF - dFp)*one_over_sin_chi;
    }
    else {is_it_successful = false; return;}
  }
  else  
  {
    const std::complex<double> one_over_2lp1 = 1.0/(2*l+1);

    if (inf_norm ((dF*Fp - dFp*F)*one_over_2lp1 - 1.0) < precision)
    {
      const std::complex<double> Cl_eta_2 = Cl_eta*Cl_eta,exp_I_omega_chi_over_sin_chi = exp_I_omega_chi*one_over_sin_chi;
      const std::complex<double> F_Cl_eta_2 = ((Cl_eta_2 == 0.0) || (!isfinite (Cl_eta_2))) ? (exp (2.0*log_Cl_eta + log (F))) : (Cl_eta_2*F);
      const std::complex<double> dF_Cl_eta_2 = ((Cl_eta_2 == 0.0) || (!isfinite (Cl_eta_2))) ? (exp (2.0*log_Cl_eta + log (dF))) : (Cl_eta_2*dF);
      
      H = exp_I_omega_chi_over_sin_chi*F_Cl_eta_2 + Fp*one_over_2lp1;
      dH = exp_I_omega_chi_over_sin_chi*dF_Cl_eta_2 + dFp*one_over_2lp1;
    }
      else {is_it_successful = false; return;}
  }
  
  is_it_successful = true;
}







// Calculation of F and F' by power series.
// ----------------------------------------
// It is used only when |z| <= 0.5, to avoid numerical inaccuracies.
//
// F(z) = norm.z^(l+1).\sum a[n], n in [0:+oo[, where :
// a[0] = 1.0.
// a[1] = z.eta/(l+1).
// a[n] = (2.z.eta.a[n-1] - a[n-2].(z^2))/(n.(n+2l+1)),n >= 2.
//
// The z = 0 case is treated first. It is defined only for Re[l] > 0 or l = 0. The program aborts for other cases.
// Norm is C(l,eta) if one uses normalized wave functions, 1.0 if not.
// So, one multiplies by C(l,eta) at the end if one uses normalized functions.
// If there is overflow or underflow for C(l,eta) in this last case, one uses logs of F,F' and C(l,eta) for the calculation.
//
// Variables:
// ----------
// z : variable of the Coulomb wave function.
// F,dF : regular wave function and derivative
// z_square,z_two_eta,z_pow_l_plus_one : z^2, 2.z.eta, z^{l+1}.
// n : index of the power series term. It begins at two.
// an,an_minus_one,an_minus_two : a[n], a[n-1], a[n-2].
// The test of convergence is |(n+l-1).a[n-2]|oo + |(n+l).a[n-1]|oo, as one of the two can be zero even before convergence.

void Coulomb_wave_functions::F_dF_power_series (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF)
{
  if (z == 0.0)
  {
    if (l == 0) F = 0.0,dF = (is_it_normalized) ? (Cl_eta) : (1.0);
    else if (real (l) > 0) F = dF = 0.0;
    else std::cout<<"F(z=0) and/or F'(z=0) are undefined."<<std::endl, abort ();
  }
  else
  {
    const std::complex<double> z_square = z*z,z_two_eta = 2.0*eta*z;

    int n = 2;
    std::complex<double> an_minus_two = 1.0,an_minus_one = z*eta/(l+1.0);
 
    F = an_minus_two + an_minus_one;
    dF = (l+1.0)*an_minus_two + (l+2.0)*an_minus_one;

    while (inf_norm (an_minus_two*(n+l-1.0)) + inf_norm (an_minus_one*(n+l)) > precision)
    {
      const std::complex<double> an = (z_two_eta*an_minus_one - an_minus_two*z_square)/(n*(n + l + l + 1.0));

      F += an;
      dF += an*(n + l + 1.0);
    
      n++;
      an_minus_two = an_minus_one;
      an_minus_one = an;
    }
    
    const std::complex<double> z_pow_l_plus_one = pow (z,l+1.0);
    F *= z_pow_l_plus_one;
    dF *= z_pow_l_plus_one/z; 
    
    if (is_it_normalized)
    {
      if ((Cl_eta == 0.0) || (!isfinite (Cl_eta))) F = exp (log_Cl_eta + log (F)),dF = exp (log_Cl_eta + log (dF));
      else F *= Cl_eta,dF *= Cl_eta;
    }
  }
}







// Calculation of f = F'/F with a continued fraction.
// --------------------------------------------------
// One calculates the ratio f = F'/F with the continued fraction of the associated hypergeometric confluent function.
// One uses Lentz's method.
// One has : f = [b[0] + a[1]/b[1]+ a[2]/b[2]+ ... a[n]/b[n]+ ...]/z with :
// b[0] = l + 1 + i.omega.z, a[n] = -2.i.omega.[1 + l + i.omega.eta] + (n-1).[-2.i.omega.z], b[n] = 2l + 2 + 2.i.omega.z + n-1.
// omega is 1 or -1, and theoretically the result is the same.
// If they are not equal numerically, omega = sign[-Im [z]] gives usually the best result.
// If 1+l+i.omega.eta is a negative integer, f[omega] is finite and must be used.
//
// Variables:
// ----------
// z : variable of the Coulomb wave function.
// omega : 1 or -1. Both values should be tried to test stability.
// large,small : 1E50,1E-50. They are used in the case of vanishing denominators or numerators.
// I_omega,a,b : i.omega, 1 + l + i.omega.eta, 2l + 2.
// minus_two_I_omega_z,minus_two_I_omega_a_z,b_plus_two_I_omega_z : -2.i.omega.z, -2.i.omega.a.z, b + 2.i.omega.z
// b0,Cn,Dn,an,bn,Delta_n : variables used in the Lentz method.
// n,nm1 : index of a[n] and b[n], n-1 
// bn_plus_an_Dn : bn + an.Dn. Dn = 1/[bn + an.Dn] or large if bn + an.Dn is zero.
// bn_plus_an_over_Cn : bn + an/Cn. Cn = bn + an/Dn or small if bn + an/Dn is zero.
// fn : value of the continuous fraction during the iteration process.
// test : test of convergence of fn.
// f : value of F'(z)/F(z).

std::complex<double> Coulomb_wave_functions::continued_fraction_f (const std::complex<double> &z,const int omega)
{
  const double small = 1E-50,large = 1E50;

  const std::complex<double> I_omega(0.0,omega);
  const std::complex<double> a = I_omega*eta + l + 1.0,b = 2*l + 2;
  const std::complex<double> minus_two_I_omega_z = -2.0*I_omega*z,minus_two_I_omega_a_z = minus_two_I_omega_z*a,b_plus_two_I_omega_z = b - minus_two_I_omega_z;

  const std::complex<double> b0 = l + 1.0 + I_omega*z;
  std::complex<double> fn = (b0 != 0.0) ? (b0) : (small), Cn = fn, Dn = 0.0;  

  int n = 1;
  double test;
  do
  {
    const int nm1 = n-1;
    const std::complex<double> an = minus_two_I_omega_a_z + nm1*minus_two_I_omega_z;
    const std::complex<double> bn = b_plus_two_I_omega_z + nm1;
    
    const std::complex<double> bn_plus_an_Dn = bn + an*Dn,bn_plus_an_over_Cn = bn + an/Cn;

    Dn = (bn_plus_an_Dn != 0.0) ? (1.0/bn_plus_an_Dn) : (large);
    Cn = (bn_plus_an_over_Cn != 0.0) ? (bn_plus_an_over_Cn) : (small);

    const std::complex<double> Delta_n = Dn*Cn;
    fn *= Delta_n;
    test = inf_norm (1.0 - Delta_n);
    n++;
  }
  while (test > 1E-15);

  const std::complex<double> f = fn/z;

  return f;
}








// Calculation of F(z) and F'(z) with asymptotic series
// ----------------------------------------------------
// F(z) = [H+(z) - H-(z)]/[2.i.norm]. 
// F'(z) = [H+'(z) - H-'(z)]/[2.i.norm]. 
// In this routine, Re[z] >= 0, so there is no branch cut problem.
//
// H+(z) = exp[i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S+(z) .
// H-(z) = exp[-i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S-(z) .
//
// H+'(z) = exp[i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S+'(z) + S+(z).i.(1 - eta/z)] .
// H-'(z) = exp[-i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S-'(z) - S-(z).i.(1 - eta/z)] .
//
//
// S+ and S- and derivatives are calculated in asymptotic_series. If is_it_successful is true, the series are meaningful. If not, one leaves the routine.
// Norm is C(l,eta) if one uses normalized wave functions, 1.0 if not.
// If there is overflow or underflow for C(l,eta) in this last case, one uses logs of F,F' and C(l,eta) for the calculation.
//
// Variables :
// -----------
// z : variable of the Coulomb wave function.
// one_over_z : 1/z
// F,dF : regular wave function and derivative to calculate.
// is_it_successful : true if the calculation converged, i.e. the series are good up to precision and the wronskian of H+,H- up to precision, false if not.
// sum,dsum : asymptotic series. sum[0] = S+(z), sum[1] = S-(z), dsum[0] = S+'(z), dsum[1] = S-'(z).
// I,one_over_two_I : i,1/[2i].
// I_one_minus_eta_over_z : i*(1 - eta/z).
// exp_phase_shift_plus : exp (i*(z - eta.log(2z) - l.Pi/2 + sigma(l,eta)))).
// exp_phase_shift_minus : exp (-i*(z - eta.log(2z) - l.Pi/2 + sigma(l,eta)))).
// H_plus,dH_plus,H_minus,dH_minus : exp (+/-i*(z - eta.log(2z) - l.Pi/2 + sigma(l,eta))).S(+/-)(z) and derivatives.

void Coulomb_wave_functions::asymptotic_expansion_F_dF (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF,bool &is_it_successful)
{
  std::complex<double> sum[2],dsum[2];

  const std::complex<double> one_over_z = 1.0/z;

  asymptotic_series (1,one_over_z,sum,dsum,is_it_successful);
  if (!is_it_successful) return;

  const std::complex<double> I(0,1),one_over_two_I(0.0,-0.5),I_one_minus_eta_over_z = I*(1.0 - eta*one_over_z);

  const std::complex<double> exp_phase_shift_plus = exp (I*(z - eta*(M_LN2 + log (z)) - l*M_PI_2 + sigma_l));
  const std::complex<double> exp_phase_shift_minus = 1.0/exp_phase_shift_plus;

  const std::complex<double> H_plus = sum[0]*exp_phase_shift_plus,dH_plus = (dsum[0] + sum[0]*I_one_minus_eta_over_z)*exp_phase_shift_plus;
  const std::complex<double> H_minus = sum[1]*exp_phase_shift_minus,dH_minus = (dsum[1] - sum[1]*I_one_minus_eta_over_z)*exp_phase_shift_minus;

  F = (H_plus - H_minus)*one_over_two_I;
  dF = (dH_plus - dH_minus)*one_over_two_I;

  if (!is_it_normalized)
  {
    if ((Cl_eta == 0.0) || (!isfinite (Cl_eta))) F = exp (log (F) - log_Cl_eta),dF = exp (log (dF) - log_Cl_eta);
    else F /= Cl_eta,dF /= Cl_eta;
  }


}






// Calculation of F(z) and F'(z) by direct integration.
// ----------------------------------------------------
// To calculate F(z) and F'(z), one integrates numerically F''(z) = [l(l+1)/z^2 + 2.eta/z - 1].F(z)
// starting from debut, F_debut = F(debut) and dF_debut = F'(debut).
// One always has Re[z] >= 0.0, so there is no branch cut problem.
// If z = debut, the previous values are returned.
// The starting point come from the stored values debut, F_debut and dF_debut.
// If debut = 0, one puts debut = z/|2z| and calculates F(debut) and F'(debut) with power series.
// The step of the integration is (z - debut)/N_num, with N_num = [|z - debut|/min (0.1,10.turning_point)] + 1.
// The value of min (0.1,10/turning_point) gives a smaller step when turning_point increases, 
// as calculations become there more difficult as |F| typically varies faster in this case.
// The intermediates points are called z_aft. They go from debut to z, and (debut,F_debut,dF_debut) is put to {z_aft,F(z_aft),F'(z_aft)} at each step.
// If |F| increases along the path, the integration is stable.
// If it decreases, and if z_aft decreases in modulus or one does not integrate with constant argument (i.e. theta constant in z = |z|.exp[i.theta]) 
// for Re[l] > -1, one reintegrates F(z) from debut = z/|2z|, as integration is usually stable at constant argument for Re[l] > -1.
// Increase or decrease is known using the Taylor expansion of F near debut in z up to second order.
// If one integrates with constant argument and |F| decreases, 
// one integrates backwards from z_aft to debut with the knowledge of f(z_aft) = F'(z_aft)/F(z_aft).
// F'(z_aft)/F(z_aft) is given by the continued fraction formula.
// One then obtains by direct integration C.F(debut) and C.F'(debut).
// Knowing F_debut, one deduces F(z) = 1/C and F'(z) = f(z_aft).F(z).
// If 1+l+i.omega.eta is a negative integer, f[omega] is finite and is used.
// Otherwise, f(z_aft) is calculated with omega = 1 and -1. If they are equal up to precision, f(omega) is correct and used.
// omega is chosen so Re[-2.i.omega.z] < 0, 
// for which the anomalous convergence phenomenon of Gautschi of f is the smallest (W. Gautschi, Math. Comp. Vol. 31 p.994).
// If not, but |norm.F| < 0.1, one still has to use f[omega], as it is probably correct as F is the minimal solution, and also one has no other way to calculate F.
// If |norm.F| > 0.1 in this case, one stops the procedure and F will be calculated from H+ and H-, given by direct integration and continued fraction formulae.
// In this case, is_it_successful is put to false, and otherwise the integration worked and it is put to true.
// Norm is 1.0 if one uses normalized wave functions, C(l,eta) if not.
//
// If F(z_aft) is not finite, one stops the integration and is_it_successful is put to false.
//
// Variables
// ---------
// z : variable of the Coulomb wave function.
// F,dF : regular wave function and derivative to calculate.
// is_it_successful : false is the calculation is unstable, 
// i.e. |F| > 0.1 decreasing on the integration path and f(omega) is not equal to f(-omega) up to precision, true if not.
// ODE : reference to pointer ODE_ptr, which performs direct integration from one point to another.
// step_abs : length of the integration step from z to debut. It is min (0.1,10/turning_point)
// N_num : number of integrations to do from debut to z. It is |z - debut|/step_abs + 1.
// step_num : std::complex step from debut to z. It is (z - debut)/N_num .
// ll_plus_one,two_eta : l(l+1), 2.eta. They are used in the Taylor expansion test.
// z_aft : next point in which F is calculated by direct integration . It is z - i*step_num, with i from N_num-1 to 0.
// one_over_debut,log_F_debut_der,d2F_debut_over_F_debut: 1/debut, F'(debut)/F(debut), F''(debut)/F(debut). They are used in the Taylor expansion test. 
// ratio : debut/z. It is used to know if one integrates with constant argument or if |z_aft| increases. 
//         If not and |F| increases and Re[l] > -1, one reintegrates from z/|2z|.
// f_omega,f_minus_omega : continued fractions F'/F with omega = sign(-Im[z_aft]) and -omega.
// fp,fm : continued fractions F'/F with omega = 1 or -1 used in the case of 1+l +/- i.eta negative integer, as continued fractions are finite in this case.
// norm_functions : C(l,eta)^2 for unnormalized functions, 1 for normalized functions.
// F_debut_not_normed,dF_debut_not_normed : unnormed regular wave function and derivative in debut,
//                                          after integration from the starting point {z_aft, 1.0, F'(z_aft)/F(z_aft)}.
//                                          Then, F(z_aft) =  F_debut/F_debut_not_normed and F'(z_aft) = f[omega].F .

void Coulomb_wave_functions::F_dF_direct_integration (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF,bool &is_it_successful)
{
  if (z == debut) {F = F_debut,dF = dF_debut,is_it_successful = true; return;}
  if (debut == 0.0) debut = 0.5*z/abs (z),F_dF_power_series (debut,F_debut,dF_debut);

  const class ODE_integration &ODE = *ODE_ptr;
  const double step_abs = std::min(0.1,10.0/turning_point);
  const unsigned int N_num = static_cast<unsigned int> (rint (abs (z-debut)/step_abs) + 1); 
  const std::complex<double> step_num = (z - debut)/static_cast<double> (N_num),ll_plus_one = l*(l+1.0),two_eta = 2.0*eta;
  const std::complex<double> norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);

  for (unsigned int i = N_num-1 ; i <= N_num ; i--)
  {
    const std::complex<double> z_aft = z - i*step_num,one_over_debut = 1.0/debut,log_F_debut_der = dF_debut/F_debut;
    const std::complex<double> d2F_debut_over_F_debut = (ll_plus_one*one_over_debut + two_eta)*one_over_debut - 1.0;

    if (abs (1.0 + step_num*(log_F_debut_der + 0.5*step_num*d2F_debut_over_F_debut)) < 1.0)
    {
      const std::complex<double> ratio = debut/z; 
      if ((real (l) > -1.0) && ((std::abs (imag (ratio)) > precision) || (real (ratio) > 1.0)))
	{debut = 0.0, F_dF_direct_integration (z,F,dF,is_it_successful); return;}

      std::complex<double> F_debut_not_normed,dF_debut_not_normed;

      if (neg_int_omega_one)
      {
	const std::complex<double> fp = continued_fraction_f (z_aft,1);
	ODE (z_aft,1.0,fp,debut,F_debut_not_normed,dF_debut_not_normed);
	F = F_debut/F_debut_not_normed; 
	dF = fp*F;
      }
      else if (neg_int_omega_minus_one)
      {
	const std::complex<double> fm = continued_fraction_f (z_aft,-1);
	ODE (z_aft,1.0,fm,debut,F_debut_not_normed,dF_debut_not_normed); 
	F = F_debut/F_debut_not_normed; 
	dF = fm*F;
      }
      else
      {
	const std::complex<double> f_omega = continued_fraction_f (z_aft,SIGN(-imag(z_aft))),f_minus_omega = continued_fraction_f (z_aft,-SIGN(-imag(z_aft)));

	//// Comment the following line if you want to accept the value of f(omega) = F'/F inconditionally
	if ((abs (F*norm_functions) > 0.1) && (abs (f_minus_omega/f_omega - 1.0) > precision)) {is_it_successful = false; return;}

	ODE (z_aft,1.0,f_omega,debut,F_debut_not_normed,dF_debut_not_normed);
	F = F_debut/F_debut_not_normed; 
	dF = f_omega*F;
      }
    }
    else ODE (debut,F_debut,dF_debut,z_aft,F,dF);

    debut = z_aft,F_debut = F,dF_debut = dF;

    if (!isfinite (F) || !isfinite (dF)) std::cout<<"Numerical failure encountered in F_dF_direct_integration."<<std::endl,exit (1);
  }
  is_it_successful = true;
}



// Regular wave function and derivative from symmetry relations.
// -------------------------------------------------------------
// If |z| > 0.5 and Re[z] < 0, one calculates F(l,eta,z),F'(l,eta,z) from F(l,-eta,-z),F'(l,-eta,-z) using the formulas :
// F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta-i.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-Pi.(eta-i.l)] if arg (z) > 0 and is_it_normalized is true,
// F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta+i.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-Pi.(eta+i.l)] if arg (z) <= 0 and is_it_normalized is true,
// F(l,eta,z) = -F(l,-eta,-z).exp[i.Pi.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[i.Pi.l)] if arg (z) > 0 and is_it_normalized is false,
// F(l,eta,z) = -F(l,-eta,-z).exp[-i.Pi.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-i.Pi.l)] if arg (z) <= 0 and is_it_normalized is false.
//
// F(l,-eta,-z) is calculated using the class cwf_minus_eta_ptr defined with (l,-eta).
// The debut point of the class cwf_minus_eta_ptr is initialized with {debut,F_debut,dF_debut} and previous relations
// if cwf_minus_eta_ptr->debut and -debut are different and debut non zero.
//
// If the normalization constant underflows or overflows, one uses logs.
//
// Variables
// ---------
// z : variable of the Coulomb wave function.
// F,dF : regular wave function and derivative to calculate.
// arg_debut, arg_z : argument angles of debut and z.
// sym_constant_debut,log_sym_constant_debut : constant C so F(l,-eta,-debut) = C.F(l,eta,debut), its log.
// sym_constant,log_sym_constant : constant C so F(l,eta,z) = C.F(l,-eta,-z), its log.

void Coulomb_wave_functions::F_dF_with_symmetry_relations (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF)
{
  if (cwf_minus_eta_ptr == 0) cwf_minus_eta_ptr = new class Coulomb_wave_functions (is_it_normalized,l,-eta);

  if ((debut != 0.0) && (cwf_minus_eta_ptr->debut != -debut))
  {
    const double arg_debut = arg (debut);
    const std::complex<double> sym_constant_debut = (arg_debut <= 0.0) ? (1.0/sym_constant_arg_neg) : (1.0/sym_constant_arg_pos);

    cwf_minus_eta_ptr->debut = -debut; 

    if ((sym_constant_debut == 0.0) || (!isfinite (sym_constant_debut))) 
    {
      const std::complex<double> log_sym_constant_debut = (arg (debut) <= 0.0) ? (-log_sym_constant_arg_neg) : (-log_sym_constant_arg_pos);

      cwf_minus_eta_ptr->F_debut = exp (log_sym_constant_debut + log (F_debut)); 
      cwf_minus_eta_ptr->dF_debut = -exp (log_sym_constant_debut + log (dF_debut));
    } 
    else
    {
      cwf_minus_eta_ptr->F_debut = F_debut*sym_constant_debut; 
      cwf_minus_eta_ptr->dF_debut = -dF_debut*sym_constant_debut;
    }
  }

  const double arg_z = arg (z);
  const std::complex<double> sym_constant = (arg_z <= 0.0) ? (sym_constant_arg_neg) : (sym_constant_arg_pos);

  cwf_minus_eta_ptr->F_dF (-z,F,dF);

  if ((sym_constant == 0.0) || (!isfinite (sym_constant))) 
  {
    const std::complex<double> log_sym_constant = (arg_z <= 0.0) ? (log_sym_constant_arg_neg) : (log_sym_constant_arg_pos);

    F = exp (log_sym_constant + log (F)); 
    dF = -exp (log_sym_constant + log (dF));
  } 
  else 
  {
    F *= sym_constant; 
    dF *= -sym_constant;
  }
}





// Calculation of the asymptotic series
// ------------------------------------
// Asymptotic expansion: 
// S(+/-)(z) = 1.0+\sum_(n=1)^N a[n] with a[n+1]=a[n].[n.[n+1+/-2i.eta]+i.eta.(i.eta+/-1)-l(l+1)]/[+/-2i.(n+1)]/z, n >= 0 and a[0] = 1.
// This expansion diverges : it is only useful with the smallest term summation method.
// The test of convergence is max(|a[n]|oo,|n.a[n]/z|oo), so the largest norm of the term of series of function and derivative.
// Practically, one stops when test < precision (it worked) or when test is not finite (it failed). 
// After that, one tests the series with the wronskian of H[omega] and H[-omega]. 
//
// Variables :
// -----------
// z is the variable of the Coulomb wave function.
// omega : 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
// one_over_z : 1/z.
// sum : sum[0] is the series in H(omega), sum[1] the one in H(-omega).
// dsum : dsum[0] is the series in H'(omega), dsum[1] the one in H'(-omega).
// is_it_successful : true if the series converged and |wronskian - 2.i.omega|oo is smaller than precision, false if not.
// test : test of convergence of the asymptotic series. It is max(|a[n]|oo,|n.a[n]/z|oo).
// sign : if i=0, it is omega, if i=1, it is -omega.
// Ieta,two_I_eta_sign : i.eta, 2.i.eta if sign = 1, -2.i.eta if sign = -1.
// Ieta_Ieta_plus_sign_minus_ll_plus_one : i.eta(i.eta+1) - l(l+1) if sign=1, i.eta(i.eta-1) - l(l+1) if sign=-1.
// n,an_sign : index of the series, a[n+1]
// n_plus_one,sign_one_over_two_I_n_plus_one : n+1, 1/(2.i.(n+1)) if sign=1, -1/(2.i.(n+1)) if sign=-1.
// sum_term, dsum_term : an_sign, (n+1).an_sign/z
// two_I_omega : 2.i.omega .

void Coulomb_wave_functions::asymptotic_series (const int omega,const std::complex<double> &one_over_z,
						std::complex<double> sum[],std::complex<double> dsum[],bool &is_it_successful)
{
  sum[0] = sum[1] = 1.0;
  dsum[0] = dsum[1] = 0.0;

  double test;

  for (unsigned int i = 0 ; i <= 1 ; i++)
  {
    const int sign = (i == 0) ? (omega) : (-omega); 
    const std::complex<double> Ieta(-imag (eta),real (eta)),two_I_eta_sign = 2.0*sign*Ieta,Ieta_Ieta_plus_sign_minus_ll_plus_one = Ieta*(Ieta + sign) - l*(l+1.0);
    
    int n = 0;
    std::complex<double> an_sign = 1.0;

    do
    {
      const double n_plus_one = n + 1.0;
      const std::complex<double> sign_one_over_two_I_n_plus_one(0,-sign*0.5/n_plus_one);

      an_sign *= one_over_z*(n*(n_plus_one + two_I_eta_sign) + Ieta_Ieta_plus_sign_minus_ll_plus_one)*sign_one_over_two_I_n_plus_one;

      const std::complex<double> sum_term = an_sign,dsum_term = n_plus_one*an_sign*one_over_z;

      sum[i] += sum_term;
      dsum[i] -= dsum_term;

      test = std::max (inf_norm (sum_term),inf_norm (dsum_term));
      n++;
    }
    while ((test > precision) && (isfinite (test)));

    if (!isfinite (test)) {is_it_successful = false; return;}
  }

  const std::complex<double> two_I_omega(0.0,2.0*omega);
  is_it_successful = (inf_norm (sum[1]*dsum[0] - sum[0]*dsum[1] + two_I_omega*(1.0 - eta*one_over_z)*sum[0]*sum[1] - two_I_omega) < precision);
}









// Numerical partial derivatives according to l or eta.
// ----------------------------------------------------
// One calculates here the partial derivatives according to l or eta
// of the Coulomb wave functions F(x) or G(x) and of their derivatives with x F'(x) or G'(x),
// with l_r=Re[l] and eta_r=Re[eta].
// One considers here the parameters l_r and eta_r with the argument x.
// For this, one uses the standard formula : df/d_chi(x) = [f(x,chi+eps) - f(x,chi-eps)]/[2.eps], 
// with chi = l_r or eta_r and eps = prec_first_order_expansion*chi. f is either F, F', G or G'.
// One then calculates F or G, F' or G' with chi+eps and chi-eps.
// With them, one obtains dF/d_chi(x),dF'/d_chi(x) or dG/d_chi(x),dG'/d_chi(x), with chi = l or eta.
// All functions are real, so double values are returned, taking the real part of std::complex variables.
//
// Variables
// ---------
// is_it_regular : true if one calculates derivatives of F(z),F'(z), false if one calculates derivatives of G(z),G'(z).
// is_it_eta : true if one calculates the partial derivatives according to eta,
//             false if one calculates the partial derivatives according to l.
// cwf_plus,cwf_minus : references on class Coulomb_wave_functions with parameters l+eps,eta and l-eps,eta or l,eta+eps and l,eta-eps.
//                      If one calculates the partial derivatives according to eta, they are *cwf_real_eta_plus_ptr and *cwf_real_eta_minus_ptr.
//                      If one calculates the partial derivatives according to l, they are *cwf_real_l_plus_ptr and *cwf_real_l_minus_ptr.
// x : Re[z], z the argument of the wave function.
// d_chi_Bx,d_chi_dBx : partial derivatives of F and F' in l_r,eta_r,x according to chi = l or eta if is_it_regular is true, or G and G' if not.
// l_r,eta_r,chi_r : Re[l],Re[eta], Re[eta] is is_it_eta is true, Re[l] if not.
// chi_r_plus,chi_r_minus : chi_r+eps, chi_r-eps.
// A_plus,dA_plus,A_minus,dA_minus : F(x,chi+eps),F'(x,chi+eps),F(x,chi-eps),F'(x,chi-eps) if is_it_regular is true,
//                                   H+(x,chi+eps),H+'(x,chi+eps),H+(x,chi-eps),H+'(x,chi-eps) if is_it_regular is false.
// B_plus,B_minus,dB_plus,dB_minus : F(x,chi+eps),F'(x,chi+eps),F(x,chi-eps),F'(x,chi-eps) if is_it_regular is true,
//                                   G(x,chi+eps),G'(x,chi+eps),G(x,chi-eps),G'(x,chi-eps) if is_it_regular is false.

void Coulomb_wave_functions::partial_derivatives (const bool is_it_regular,const bool is_it_eta,const double x,double &d_chi_Bx,double &d_chi_dBx)
{
  const double l_r = real (l),eta_r = real (eta),chi_r = (is_it_eta) ? (eta_r) : (l_r);
  const double chi_r_plus = (chi_r != 0.0) ? (chi_r*(1.0 + prec_first_order_expansion)) : (prec_first_order_expansion);
  const double chi_r_minus = (chi_r != 0.0) ? (chi_r*(1.0 - prec_first_order_expansion)) : (-prec_first_order_expansion);

  if (is_it_eta)
  { 
    if (cwf_real_eta_plus_ptr == 0) cwf_real_eta_plus_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,chi_r_plus);
    if (cwf_real_eta_minus_ptr == 0) cwf_real_eta_minus_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,chi_r_minus);
  }
  else
  {
    if (cwf_real_l_plus_ptr == 0) cwf_real_l_plus_ptr = new class Coulomb_wave_functions (is_it_normalized,chi_r_plus,eta_r);
    if (cwf_real_l_minus_ptr == 0) cwf_real_l_minus_ptr = new class Coulomb_wave_functions (is_it_normalized,chi_r_minus,eta_r);
  }
 
  class Coulomb_wave_functions &cwf_plus = (is_it_eta) ? (*cwf_real_eta_plus_ptr) : (*cwf_real_l_plus_ptr);
  class Coulomb_wave_functions &cwf_minus = (is_it_eta) ? (*cwf_real_eta_minus_ptr) : (*cwf_real_l_minus_ptr);

  std::complex<double> A_plus,dA_plus,A_minus,dA_minus;

  if (is_it_regular)
  {
    cwf_plus.F_dF (x,A_plus,dA_plus);
    cwf_minus.F_dF (x,A_minus,dA_minus);
  }
  else
  {
    cwf_plus.H_dH (1,x,A_plus,dA_plus);
    cwf_minus.H_dH (1,x,A_minus,dA_minus);
  }

  const double B_plus = real (A_plus),B_minus = real (A_minus),dB_plus = real (dA_plus),dB_minus = real (dA_minus);

  d_chi_Bx = (B_plus - B_minus)/(chi_r_plus - chi_r_minus);
  d_chi_dBx = (dB_plus - dB_minus)/(chi_r_plus - chi_r_minus);
}








// Calculation of F(z),F'(z) or G(z),G'(z) with the first order expansion method.
// ------------------------------------------------------------------------------
// When imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
// one has to separate the calculation of the real and imaginary parts of (F,G)(z) and (F',G')(z), as they can differ by tens of orders of magnitude.
// Re[z] < 0 is not considered, as G(z) is std::complex with Im[z] = Im[l] = Im[eta] = 0.
// So, with z = x+i.y, eta = eta_r + i.eta_i and l = l_r+i.l_i, one calculates (F,G)[l_r,eta_r](x) and (F,G)'[l_r,eta_r](x).
//
// One considers here the parameters l_r and eta_r with the argument x, and the parameters l and eta with the argument z.
//
// One then has F(x),F'(x) from usual relations and one takes their real parts.
// One also has H+(x),H+'(x) from usual relations and G(x) = Re[H+(x)],G'(x) = Re[H+'(x)].
//
// After that, one expands (F,G)(z) and (F',G')(z) in first order in y, eta_i and l_i and one has, up to y^2,eta_i^2 and l_i^2:
//
// (F,G)(z) = (F,G)(x) + i.y.d(F,G)/dx[omega](x) + i.eta_i.d(F,G)/d[eta](x) + i.l_i.d(F,G)/dl(x).
// (F',G')[omega](z) = (F',G')[omega](x) + i.y.(F'',G'')[omega](x) + i.eta_i.d(F',G')/d[eta](x) + i.l_i.d(F',G')/dl(x).
//
// Variables:
// ----------
// is_it_regular : true if one calculates F(z),F'(z), false if one calculates G(z),G'(z).
// z : variable of the Coulomb wave function.
// B,dB : F(z),F'(z) if is_it_regular is true, G(z),G'(z) if it is false.
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// A_x,dA_x : F(x),F'(x) if is_it_regular is true, H+(x),H+'(x) if it is false.
// Bx,dBx : F(x),F'(x) if is_it_regular is true, G(x),G'(x) if it is false.
// cwf_real : class Coulomb_wave_functions of parameters l_r and eta_r.
// one_over_x,two_eta_r : 1/x, 2.eta_r.
// d_l_Bx,d_l_dBx,d_eta_Bx,d_eta_dBx : partial derivatives according to l and eta of B(x) and B'(x)
//                                     They are initialized at zero, and not calculated if they are multiplied by zero after.
// Coulomb_eq_x,d2Bx : l_r(l_r+1)/[x^2] + 2.eta_r/x - 1, B''(x) = [l_r(l_r+1)/[x^2] + 2.eta_r/x - 1].B(x),


void Coulomb_wave_functions::first_order_expansions (const bool is_it_regular,const std::complex<double> &z,std::complex<double> &B,std::complex<double> &dB)
{
  const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);
    
  if (cwf_real_ptr == 0) cwf_real_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,eta_r);
  class Coulomb_wave_functions &cwf_real = *cwf_real_ptr;
    
  std::complex<double> A_x,dA_x;
  if (is_it_regular) 
    cwf_real.F_dF (x,A_x,dA_x);
  else
    cwf_real.H_dH (1,x,A_x,dA_x);

  const double Bx = real (A_x),dBx = real (dA_x);
  
  double d_l_Bx = 0.0,d_l_dBx = 0.0,d_eta_Bx = 0.0,d_eta_dBx = 0.0;
  if (l_i != 0.0) partial_derivatives (is_it_regular,false,x,d_l_Bx,d_l_dBx);
  if (eta_i != 0.0) partial_derivatives (is_it_regular,true,x,d_eta_Bx,d_eta_dBx);
  
  const double one_over_x = 1.0/x,Coulomb_eq_x = (l_r*(l_r+1)*one_over_x + 2.0*eta_r)*one_over_x - 1.0;
  const double d2Bx = Coulomb_eq_x*Bx;

  B = std::complex<double> (Bx,y*dBx + l_i*d_l_Bx + eta_i*d_eta_Bx); 
  dB = std::complex<double> (dBx,y*d2Bx + l_i*d_l_dBx + eta_i*d_eta_dBx);
}













// Regular wave function and derivative.
// -------------------------------------
// One calculates F(z) and F'(z), so F(z) ~ z^{l+1} for z -> 0 if is_it_normalized is false, 
//                                   F(z) ~ C(l,eta).z^{l+1} for z -> 0 if is_it_normalized is true.
// If |z| <= 0.5, one uses the power series.
//
// If |z| > 0.5 and Re[z] < 0, one calculates F(z) from F[l,-eta,-z] with F_dF_with_symmetry_relations.
//
// If |z| > 0.5 and Re[z] >= 0, and 1+l+/-i.eta no negative integer, one first tries the asymptotic series formula.
// If it failed, one integrates directly the regular Coulomb wave function with F_dF_direct_integration.
// If 1+l+/-i.eta is a negative integer, this is the only available method besides power series so it is accepted.
// If 1+l+/-i.eta is no negative integer but it failed again, 
// one calculates the Coulomb wave function H[omega] with H_dH_direct_integration and omega = sign(Im[z]).
// omega is chosen so one cannot encounter the branch cut of h[omega].
// H[-omega] is calculated from H[omega] and continued fractions h[omega] and h[-omega].
// One then has F(z) = (H[omega] - H[-omega])/(2.i.omega.norm), F'(z) = (H'[omega] - H'[-omega])/(2.i.omega.norm),
// with norm = 1 if the wave functions are normalized and C(l,eta)^2 if not.
// The formula is stable as one uses this case only when |F(z)| > 0.1 .
//
// One takes only real parts if l, eta and z are real.
// At the end of the function, one puts {debut,F_debut,dF_debut} equal to {z,F,dF}.
//
// Variables:
// ----------
// z : variable of the Coulomb wave function.
// F,dF : regular Coulomb wave function in z and derivative in z.
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// is_it_successful : true if the asymptotic expansions converges, false it not (after asymptotic_expansion_F_dF).
//                    true if the direct integration worked, false it not (after F_dF_direct_integration).
// omega : sign[Im[z]]. H[omega](z) is calculated if the direct integration of F failed.
// two_I_omega, two_I_term : 2.i.omega, 2.i.omega if is_it_normalized is true, 2.i.omega.Cl_eta^2 if not.
// h_omega, h_minus_omega : log derivatives of H[omega](z) and H[-omega](z) calculated with continued fractions.
// H_omega,dH_omega,H_minus_omega,dH_minus_omega : H[omega](z), H[-omega](z) and derivatives

void Coulomb_wave_functions::F_dF (const std::complex<double> &z,std::complex<double> &F,std::complex<double> &dF)
{  
  const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);

  if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
      && (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
      && (!neg_int_omega_one && !neg_int_omega_minus_one))
    first_order_expansions (true,z,F,dF);
  else if (abs (z) <= 0.5)
    F_dF_power_series (z,F,dF);
  else 
  {
    if (real (z) < 0.0) 
      F_dF_with_symmetry_relations (z,F,dF);
    else
    {
      bool is_it_successful = false;
      if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_F_dF (z,F,dF,is_it_successful);

      if (!is_it_successful) 
      {
	F_dF_direct_integration (z,F,dF,is_it_successful);

	if (!neg_int_omega_one && !neg_int_omega_minus_one && !is_it_successful)
	{
	  const int omega = SIGN(imag (z));
	  const std::complex<double> two_I_omega(0.0,2.0*omega),two_I_term = (is_it_normalized) ? (two_I_omega) : (two_I_omega*Cl_eta*Cl_eta);
	  const std::complex<double> h_omega = continued_fraction_h (z,omega),h_minus_omega = continued_fraction_h (z,-omega),one_over_two_I_term = 1.0/two_I_term;

	  std::complex<double> H_omega,dH_omega;
	  H_dH_direct_integration (omega,z,H_omega,dH_omega);
	  
	  const std::complex<double> H_minus_omega = two_I_term/((h_omega - h_minus_omega)*H_omega),dH_minus_omega = h_minus_omega*H_minus_omega;  	  
	  F = (H_omega - H_minus_omega)*one_over_two_I_term;
	  dF = (dH_omega - dH_minus_omega)*one_over_two_I_term;
	}}}
  }

  if (!isfinite (F) || !isfinite (dF)) std::cout<<"Numerical failure encountered in F_dF."<<std::endl,exit (1);

  if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0)) F = real (F),dF = real (dF); 
  debut = z,F_debut = F,dF_debut = dF;
}









// Calculation of G(z) and G'(z).
// ------------------------------
// One calculates the irregular Coulomb wave function from H+ and F.
// If 1+l+i.omega.eta is a negative integer, G is by definition H[-omega]. 
// If not, one uses the formulas : 
// G(z) = H+(z) - i.F(z), G'(z) = H+'(z) - i.F'(z) if is_it_normalized is true,
// G(z) = H+(z) - i.Cl_eta^2.F(z), G'(z) = H+'(z) - i.Cl_eta^2.F'(z) if not.
// There is no numerical inaccuracy as G is never a minimal solution.
// One takes only real parts if l, eta and z are real.
//
// Variables :
// -----------
// z : variable of the Coulomb wave function.
// G,dG : irregular Coulomb wave functions and derivatives.
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// H_plus,dH_plus : H+(z) and H+'(z).
// F,dF : regular Coulomb wave functions and derivatives.
// I_Cl_eta_square : i.C(l,eta)^2.

void Coulomb_wave_functions::G_dG (const std::complex<double> &z,std::complex<double> &G,std::complex<double> &dG)
{
  const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);

  if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
      && (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
      && (!neg_int_omega_one && !neg_int_omega_minus_one))
    first_order_expansions (false,z,G,dG);
  else
  {
    if (neg_int_omega_one) 
      H_dH (-1,z,G,dG);
    else if (neg_int_omega_minus_one)
      H_dH (1,z,G,dG);
    else
    {
      std::complex<double> F,dF;
      F_dF (z,F,dF);

      std::complex<double> H_plus,dH_plus;
      H_dH (1,z,H_plus,dH_plus);

      const std::complex<double> I(0.0,1.0);

      if (is_it_normalized)
	G = H_plus - I*F,dG = dH_plus - I*dF;
      else
      {
	const std::complex<double> I_Cl_eta_square = I*Cl_eta*Cl_eta;
	
	if ((I_Cl_eta_square == 0.0) || (!isfinite (I_Cl_eta_square)))
	  G = H_plus - I*exp (2.0*log_Cl_eta + log (F)),dG = dH_plus - I*exp (2.0*log_Cl_eta + log (dF));
	else
	  G = H_plus - I_Cl_eta_square*F,dG = dH_plus - I_Cl_eta_square*dF;
      }
      
      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0)) G = real (G),dG = real (dG);
    }
  }
}




// Calculation of H[omega](z) and H'[omega](z). 
// --------------------------------------------
// One first tries the asymptotic expansion formula if 1+l+/-i.eta is no negative integer.
// On uses logs if the unscaling factor underflows or overflows.
// If it failed, and imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
// one calculates H[omega](z) and H'[omega](z) with the first order expansion method.
// If one is not in this case, one calculates F(z)and F'(z).
// If |Im[l]| > 1 and |z| <= 1, one tries the expansion formula with F(l,eta,z) and F(-l-1,eta,z).
// If not, or if it failed, one uses the continued fraction formula.
// If l,eta and z are real, one rewrites H[omega] as H[omega] = Re[H[omega]] + i.omega.norm.Re[F] to avoid numerical inaccuracies for Im[H[omega]].
// norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
// 
//
// Variables :
// -----------
// omega : 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
// z : variable of the Coulomb wave function.
// H,dH : H+(z) and H+'(z) if omega=1, H-(z) and H-'(z) if omega=-1.
// is_it_successful : true if the asymptotic expansions converges, false it not (in asymptotic_expansion_H_dH_scaled).
//                    true if the expansion formula with F(l,eta,z) and F(-l-1,eta,z) worked, false it not (in H_dH_with_expansion).
// H_scaled,dH_scaled : H[omega](z).exp[-i.omega.[z - eta.log[2z]]) and H'[omega](z).exp[-i.omega.z.[z - eta.log[2z]]) given by the asymptotic expansion formula.
// I_omega_z,log_unscale,unscale : i.omega,i.omega.[z - eta.log[2z]],exp[i.omega.[z - eta.log[2z]]]
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// F,dF : regular Coulomb wave function and derivative in z.
// omega_norm_functions : for l,eta,z real, omega.C(l,eta)^2 if is_it_normalized is false, omega if it is true.
//                        One indeed has H[omega] = G + I.omega.norm.F, with in this case G and omega.norm.F as its real and imaginary parts.

void Coulomb_wave_functions::H_dH (const int omega,const std::complex<double> &z,std::complex<double> &H,std::complex<double> &dH)
{
  bool is_it_successful = false;

  std::complex<double> H_scaled,dH_scaled;
  if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_H_dH_scaled (omega,1.0/z,H_scaled,dH_scaled,is_it_successful);

  if (is_it_successful)
  {
    const std::complex<double> I_omega(0,omega),log_unscale = I_omega*(z - eta*(M_LN2 + log (z))),unscale = exp (log_unscale);

    if ((unscale == 0.0) || (!isfinite (unscale)))
      H = exp (log (H_scaled) + log_unscale),dH = exp (log (dH_scaled) + log_unscale);
    else    
      H = H_scaled*unscale,dH = dH_scaled*unscale;
  }
  else
  {
    const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);
    
    if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
	&& (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
	&& (!neg_int_omega_one && !neg_int_omega_minus_one))
      H_dH_from_first_order_expansions (omega,z,H,dH);
    else
    { 
      //// Replace the following line by : if (!neg_int_omega_one && !neg_int_omega_minus_one) H_dH_with_expansion (omega,z,H,dH,is_it_successful);
      //// if you want H_dH_with_expansion to be used if |l_i| < 1 or |z| > 1.
      if (!neg_int_omega_one && !neg_int_omega_minus_one && (std::abs (l_i) >= 1.0) && (std::abs (z) <= 1.0)) H_dH_with_expansion (omega,z,H,dH,is_it_successful);

      if (!is_it_successful) H_dH_with_F_dF_and_CF (omega,z,H,dH);

      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0))
      {
	std::complex<double> F,dF;
	F_dF (z,F,dF);

	const double omega_norm_functions = (!is_it_normalized) ? (omega*real (Cl_eta)*real (Cl_eta)) : (omega);
	H = std::complex<double> (real (H),omega_norm_functions*real (F));
	dH = std::complex<double> (real (dH),omega_norm_functions*real (dF));
      }
    }
  }

  if (!isfinite (H) || !isfinite (dH)) std::cout<<"Numerical failure encountered in H_dH."<<std::endl,exit (1);
}







// Calculation of the scaled H[omega](z) and H'[omega](z). 
// -------------------------------------------------------
// They are H(omega)(z).exp[-i.omega.[z - eta.log[2z]]] and dH(omega)/dz(z).exp[-i.omega.[z - eta.log[2z]]].
// One first tries the asymptotic expansion formula if 1+l+/-i.eta is no negative integer.
// If it failed, and imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
// one calculates H[omega](z) and H'[omega](z) with the first order expansion method.
// If one is not in this case, one calculates F(z) and F'(z).
// If |Im[l]| > 1 and |z| <= 1, one tries the expansion formula with F(l,eta,z) and F(-l-1,eta,z).
// If not, or if it failed, one uses the continued fraction formula.
// If l,eta and z are real, one rewrites H[omega] as H[omega] = Re[H[omega]] + I.omega.norm.Re[F], to avoid numerical inaccuracies for Im[H[omega]].
// norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
// One uses logs if the scaling factor underflows or overflows.
//
// Variables :
// -----------
// omega : 1 if one calculates H+(z) and H+'(z) scaled, -1 if one calculates H-(z) and H-'(z) scaled.
// z : variable of the Coulomb wave function.
// H_scaled,dH_scaled : H[omega](z).exp[-i.omega.[z - eta.log[2z]]) and H'[omega](z).exp[-i.omega.[z - eta.log[2z]]).
// is_it_successful : true if the asymptotic expansions converges, false it not (in asymptotic_expansion_H_dH_scaled).
//                    true if the expansion formula with F(l,eta,z) and F(-l-1,eta,z) worked, false it not (in H_dH_with_expansion).
// x,y,l_r,l_i,eta_r,eta_i : Re[z], Im[z], Re[l], Im[l], Re[eta], Im[eta].
// F,dF : regular Coulomb wave function and derivative in z.
// H,dH : H+(z) and H+'(z) if omega=1, H-(z) and H-'(z) if omega=-1.
// omega_norm_functions : for l,eta,z real, omega.C(l,eta)^2 if is_it_normalized is false, omega if it is true.
//                        One indeed has H[omega] = G + I.omega.norm.F, with in this case G and omega.norm.F as its real and imaginary parts.
// I_omega_z,log_scale,scale : i.omega,-i.omega.[z - eta.log[2z]],exp[-i.omega.[z - eta.log[2z]]]

void Coulomb_wave_functions::H_dH_scaled (const int omega,const std::complex<double> &z,std::complex<double> &H_scaled,std::complex<double> &dH_scaled)
{ 
  bool is_it_successful = false;
  if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_H_dH_scaled (omega,1.0/z,H_scaled,dH_scaled,is_it_successful);

  if (!is_it_successful)
  {
    std::complex<double> H,dH;
   
    const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);

    if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
	&& (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
	&& (!neg_int_omega_one && !neg_int_omega_minus_one))
      H_dH_from_first_order_expansions (omega,z,H,dH);
    else 
    { 
      //// Replace the following line by : if (!neg_int_omega_one && !neg_int_omega_minus_one) H_dH_with_expansion (omega,z,H,dH,is_it_successful);
      //// if you want H_dH_with_expansion to be used if |l_i| < 1 or |z| > 1.
      if (!neg_int_omega_one && !neg_int_omega_minus_one && (std::abs (l_i) >= 1.0) && (abs (z) <= 1.0)) H_dH_with_expansion (omega,z,H,dH,is_it_successful);

      if (!is_it_successful) H_dH_with_F_dF_and_CF (omega,z,H,dH);
      
      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0))
      {  
	std::complex<double> F,dF;
	F_dF (z,F,dF);

	const double omega_norm_functions = (!is_it_normalized) ? (omega*real (Cl_eta)*real (Cl_eta)) : (omega);
	H = std::complex<double> (real (H),omega_norm_functions*real (F));
	dH = std::complex<double> (real (dH),omega_norm_functions*real (dF));
      }
    }
  
    const std::complex<double> I_omega(0,omega),log_scale = -I_omega*(z - eta*(M_LN2 + log (z))),scale = exp (log_scale);
    
    if ((scale == 0.0) || (!isfinite (scale)))
      H_scaled = exp (log (H) + log_scale),dH_scaled = exp (log (dH) + log_scale);
    else
      H_scaled = H*scale,dH_scaled = dH*scale; 
  }

  if (!isfinite (H_scaled) || !isfinite (dH_scaled)) std::cout<<"Numerical failure encountered in H_dH_scaled."<<std::endl,exit (1);
}


// Storage of initial conditions debut,F(debut),F'(debut)
// ------------------------------------------------------

void Coulomb_wave_functions::F_dF_init (const std::complex<double> &z,const std::complex<double> &F,const std::complex<double> &dF)
{
  debut = z, F_debut = F, dF_debut = dF;
}

