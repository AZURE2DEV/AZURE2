#include "ode_int.H"

// Extrapolation in h=0 of a table of function values h close to h=0
// -----------------------------------------------------------------
//
// Variables:
// ----------
// n : number of points of the function f near h=0.
// T : table containing the points f[h(0)]...f[h(n-1)] close to h=0.
// f_in_zero : extrapolated value of the points f[h(0)]...f[h(n-1)] in h=0.

 
std::complex<double> ODE_integration::extrapolation_in_zero (const unsigned int n,const std::complex<double> T[]) const
{  
  std::complex<double> f_in_zero = 0.0;

  for (unsigned int i = 0 ; i < n ; i++)
    f_in_zero += interpolation_term_tab[n][i]*T[i];

  return f_in_zero;
}



// Calculation of F(z,u(z)) in u''(z) = F(z,u(z))  
// ----------------------------------------------
//
// F(z,u(z))=(l(l+1)/(z^2) + 2.eta/z - 1).u(z), 
//
// Variables:
// ----------
// z : parameter of the wave function.
// u : discretized wave function in z.
// one_over_z : 1.0/z

std::complex<double> ODE_integration::F_r_u (const std::complex<double> &z,const std::complex<double> &u) const
{
  if (l == 0) return (two_eta/z - 1.0)*u;

  const std::complex<double> one_over_z = 1.0/z;

  return ((ll_plus_one*one_over_z + two_eta)*one_over_z - 1.0)*u;
}



// Integration with discretization of u''(r)=F(r,u(r)) with the Henrici method.
// ----------------------------------------------------------------------------
//
// See Numerical Recipes for the method.
//
// Initials conditions : r0,u(r0),du/dr(r0).
// Obtained functions : r,u(r),du/dr(r).
//
// Variables:
// ----------
// m : number of intervals between r0 and r
// h : integration step (r-r0)/m .
// r0,u0,du0 : r0,u(r0),du/dr(r0).
// r,u,du : r,u(r),du/dr(r).
// h_square : h*h
// half_h = 0.5*h
// delta : value used in the Henrici method.

void ODE_integration::integration_Henrici (const unsigned int m,const std::complex<double> &h,
					   const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
					   const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const
{
  const std::complex<double> h_square = h*h,half_h = 0.5*h;

  std::complex<double> delta = h*(du0 + half_h*F_r_u (r0,u0));
  u = u0 + delta;

  for (unsigned int i = 1 ; i < m ; i++)
  {
    delta += h_square*F_r_u (r0 + i*h,u);
    u += delta;
  }
  
  du = delta/h + half_h*F_r_u (r,u);
}










// Integration of u''(r) = F(r,u(r)) with the Bulirsch-Stoer-Henrici method.
// -------------------------------------------------------------------------
//
// Initials conditions : r0,u0=u(r0),du0=du/dr(r0)
// Obtained functions : r,u=u(r),du=du/dr(r)
//
// See Numerical Recipes for the method.
//
// Variables:
// ---------- 
// r0,u0,du0 : r0,u(r0),du/dr(r0).
// r,u,du : r,u(r),du/dr(r).
// H,r_debut,r_end,u_debut,du_debut : length of an integration interval, debut and end of the integration interval, u(r_debut), u'(r_debut).
//                                    H is equal to r-r0 at the beginning and is divided by 2 each time the extrapolation fails with 16 sub-intervals 
//                                    between r_debut and r_end. If H = [r-r0]/N, with N=2,4,8,16,...,
//                                    the integration intervals are [r_debut = r0:r_end = r0+H], ..., [r_debut = r0+(N-1).H,r_end = r].
// u_end,du_end,k : tables of u(r_end) and u'(r_end) values calculated with the Henrici method with 2,4,6,...,2(k+1) sub-intervals between r_debut and r_end,
//                  with 0 <= k <= 7.
// H_over_m_tab : H/m for m=2,4,6,...,16.
// inf_norm_half_H : |H/2|oo. It is used to know if r = r_debut up to numerical accuracy, as one has |r-r_debut|oo <= |H/2|oo for this case only.
// u_extrapolated,u_extrapolated_next : values of u extrapolated from the points of the table u_end with k and k->k+1 points, k >= 2.
// test : test to know if the method worked, i.e., |u_extrapolated/u_extrapolated_next - 1|oo < precision.
// du_extrapolated_next : value of u'(r_end) extrapolated from k points of the table du_end, k >= 3.
// r_debut_plus_H : r_debut+H. r_debut+H at the end of integration is not necessarily r because of numerical cancellations.
//                             In this case, r_end must be put equal to r.

void ODE_integration::operator() (const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
				  const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const
{
  if (r == r0) {u = u0; du = du0; return;}

  std::complex<double> r_debut = r0,u_debut = u0,du_debut = du0,H = r-r0,u_end[8],du_end[8],u_extrapolated_next,du_extrapolated_next;
  double test = 1.0;

  while (test > precision)
  {
    std::complex<double> H_over_m_tab[8];
    for (unsigned int k = 0 ; k < 8 ; k++) H_over_m_tab[k] = H*one_over_m_tab[k];
    const double inf_norm_half_H = inf_norm (H_over_m_tab[0]);

    while (inf_norm (r_debut - r) > inf_norm_half_H)
    {
      const std::complex<double> r_debut_plus_H = r_debut + H, r_end = (inf_norm (r - r_debut_plus_H) > inf_norm_half_H) ? (r_debut_plus_H) : (r);

      integration_Henrici (2,H_over_m_tab[0],r_debut,u_debut,du_debut,r_end,u_end[0],du_end[0]);
      integration_Henrici (4,H_over_m_tab[1],r_debut,u_debut,du_debut,r_end,u_end[1],du_end[1]);
      std::complex<double> u_extrapolated = extrapolation_in_zero (2,u_end); 

      unsigned int k = 2; 
      do
      {
	integration_Henrici (m_tab[k],H_over_m_tab[k],r_debut,u_debut,du_debut,r_end,u_end[k],du_end[k]);
	u_extrapolated_next = extrapolation_in_zero (++k,u_end);
	test = inf_norm (u_extrapolated/u_extrapolated_next - 1.0);
	u_extrapolated = u_extrapolated_next;
      }
      while ((test > precision) && (k < 7));

      r_debut += H;
      u_debut = u_extrapolated_next;
      du_debut = du_extrapolated_next = extrapolation_in_zero (k,du_end);
    }

    H *= 0.5;
    r_debut = r0;
    u_debut = u0;
    du_debut = du0;
  }

  u = u_extrapolated_next;
  du = du_extrapolated_next;
}
