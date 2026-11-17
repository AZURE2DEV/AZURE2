#ifndef AZURE2_PHYSICS_TRANSFER_FORMFACTOR_H
#define AZURE2_PHYSICS_TRANSFER_FORMFACTOR_H

#include <string>
#include <vector>

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The full form factor of the THM entrance vertex (`vertexModel=ff` of a <thm> experiment;
 * docs/source/theory/thm_implementation.rst, "Full form factor").
 *
 * La Cognata et al., ApJ 708 (2010) 796, Eqs. 20-21 write the transfer
 * amplitude to level lambda with the form factor W_l(p) = <phi_lambda|V|j_l(pr)>,
 * phi_lambda = u_lambda/r a real single-particle level function in the nuclear
 * interior.  Inside the channel radius a it obeys (T + V) u = E_lambda u with
 * the logarithmic derivative a u'(a)/u(a) = B_lambda, so Green's theorem gives
 * exactly
 *
 *   W_l = (hbar^2/2mu) u(a) [ M_l(rho; B_lambda) + Delta_lambda(E) ],
 *   M_l = (B - 1) j_l(rho) - rho j_l'(rho)                      (surface),
 *   Delta_lambda = (2mu/hbar^2) (E_lambda - eps) / u(a) Int_0^a u r j_l(pr) dr (volume),
 *
 * with eps = hbar^2 p^2/2mu = E + B_xs + T_s and rho = p a.  The surface term
 * is AZURE2's plane-wave vertex (gamma_lambda M_l, since gamma^2 = (hbar^2/2mu a)
 * u(a)^2 up to the spectroscopic factor); vertexModel=ff replaces M_l by
 * M_l + Delta_lambda, i.e. multiplies the vertex of each level by
 * R_lambda(E) = 1 + Delta_lambda/M_l.  Delta does not depend on the
 * normalization of u.
 *
 * The level function is the regular solution in a Woods-Saxon well
 * -V0/(1 + exp((r - R)/a0)), R = r0 A^(1/3) (A the mass number of the heavier
 * nucleus of the entrance pair, the core), plus the potential of a uniformly
 * charged sphere of radius R for a charged pair; no spin-orbit term.  For a
 * given node count n_r the depth V0 is the one (it is unique) for which
 * a u'(a)/u(a) = B_lambda at E = E_lambda.
 *
 * Numerics: RK4 on r_i = i a/N (N = 1000), started at r_1 from the series
 * r^(l+1) (1 + c1 r^2 + c2 r^4); the depth from the Pruefer angle
 * theta(V0) = pi n(V0) + acot(a u'/u), continuous and increasing in V0, by
 * the Illinois method; the overlap Int u r j_l(pr) dr by Simpson's rule on the
 * same grid; as a function of p (entire in p) by a Chebyshev series on
 * [0, p_max] from 48 nodes, built once per level function.
 */
namespace form_factor {


/// The single-particle well of one entrance pair, and the radial grid.
struct Well {
  double mu = 0.0;     ///< reduced mass, MeV
  double a = 0.0;      ///< channel radius, fm
  double R = 0.0;      ///< Woods-Saxon radius, fm
  double diff = 0.65;  ///< Woods-Saxon diffuseness, fm
  double rc = 0.0;     ///< Coulomb (uniform sphere) radius, fm
  double zze2 = 0.0;   ///< Z1 Z2 e^2 (MeV fm); 0: no Coulomb
  int steps = 1000;    ///< radial steps on [0, a] (even)
  /// hbar^2/2mu (MeV fm^2).
  double H2M() const;
  /// Woods-Saxon + Coulomb at r for depth V0 (MeV).
  double Potential(double V0, double r) const;
};

/// One level function u_lambda of the well.
struct LevelFunction {
  bool ok = false;
  std::string why;         ///< when !ok
  int l = 0, nodes = 0;
  double e = 0.0;          ///< E_lambda relative to the pair's threshold (MeV)
  double B = 0.0;          ///< a u'(a)/u(a)
  double depth = 0.0;      ///< V0 (MeV)
  double ua = 0.0;         ///< u(a) > 0, with Int_0^a u^2 dr = 1
  std::vector<double> u;   ///< on r_i = i a/steps, normalized, u(a) > 0
};

/// The level function of `nodes` interior nodes, eigenvalue e and boundary B.
LevelFunction SolveLevel(const Well &w, int l, double e, double B, int nodes);
/// a u'(a)/u(a), the node count in (0, a] and the normalized u of the regular solution at depth V0.
void SolveAtDepth(const Well &w, int l, double e, double V0, double &logDerivative, int &nodes,
                  std::vector<double> *u = nullptr);
/// Int_0^a u r j_l(p r) dr (Simpson; the normalized u of `f`).
double Overlap(const Well &w, const LevelFunction &f, double p);
/// Int_0^a u V r j_l(p r) dr, V the well at f.depth: the form factor W_l computed directly.
double PotentialOverlap(const Well &w, const LevelFunction &f, double p);
/// M_l(rho; B) = (B - 1) j_l(rho) - rho j_l'(rho) with the analytic derivative.
double SurfaceM(int l, double rho, double B);
/// Delta = (E_lambda - eps)/(hbar^2/2mu) * Overlap/u(a): the volume term.
double VolumeTerm(const Well &w, const LevelFunction &f, double eps);
/// Spherical Bessel j_l(x), x >= 0 (series below max(2.5, l + 1), upward recursion above).
double SphericalJ(int l, double x);

/*!
 * The Wildermuth node count of a nucleon x on a core: x occupies the lowest
 * oscillator shell N not filled by the core's nucleons of its kind
 * (`coreCount` protons for a proton, neutrons for a neutron; shells hold
 * (N+1)(N+2) of each kind), raised to the first N >= l with N - l even;
 * n_r = (N - l)/2.  p + 18O: N = 2 (1s, 0d), 3 (1p, 0f).
 */
int WildermuthNodes(int l, int coreCount);

/// A Chebyshev series on [lo, hi].
struct Chebyshev {
  double lo = 0.0, hi = 1.0;
  std::vector<double> c;
  double operator()(double x) const;
  /// The Chebyshev nodes x_k of n terms on [lo, hi].
  static std::vector<double> Nodes(double lo, double hi, int n);
  /// The coefficients from the values at Nodes(lo, hi, n).
  void Fit(double lo, double hi, const std::vector<double> &values);
};

}  // namespace form_factor
}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
