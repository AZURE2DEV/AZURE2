#ifndef ODE_INT_H
#define ODE_INT_H

#include "complex_functions.H"

// Direct integration of the Coulomb equation
// ------------------------------------------
// One uses the Burlisch-Stoer-Henrici method, where one integrates on different meshes
// with the Henrici method, and then uses the Richardson method to get the final result by extrapolation.
// Numerical Recipes, Chap. 16.4 .

class ODE_integration
{
public:
  ODE_integration (const std::complex<double> &l_1,
		   const std::complex<double> &two_eta_1)
    : l (l_1), ll_plus_one (l_1*(l_1+1.0)), two_eta (two_eta_1)
    {
      for (int n = 0 ; n < 8 ; n++)
	for (int i = 0 ; i < n ; i++)
	{
	  interpolation_term_tab[n][i] = 1.0;

	  for (int j = 0 ; j < n ; j++)
	    if (i != j)
	      interpolation_term_tab[n][i] *= (i+1.0)/(i-j);
	}

      for (unsigned int k = 0 ; k < 8 ; k++) m_tab[k] = 2*(k+1);
      for (unsigned int k = 0 ; k < 8 ; k++) one_over_m_tab[k] = 1.0/static_cast<double> (m_tab[k]);
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

  unsigned int m_tab[8];                                 // integers used in the extrapolation method.
  double one_over_m_tab[8],interpolation_term_tab[8][8]; // doubles used in the extrapolation method.
};


#endif
