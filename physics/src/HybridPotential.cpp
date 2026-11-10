#include "azure2/physics/HybridPotential.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace physics {

namespace {

/// e^2 (MeV fm) of the Coulomb term V_C.
constexpr double kCoulombE2_MeV_fm = 1.4399764;

}  // namespace

bool operator==(const FormFactor &a, const FormFactor &b) {
  return a.shape == b.shape && a.radius_fm == b.radius_fm && a.diffuseness_fm == b.diffuseness_fm &&
         a.width_fm == b.width_fm;
}

bool operator==(const HybridPotential &a, const HybridPotential &b) {
  return a.coulomb_radius_fm == b.coulomb_radius_fm && a.central_depth_MeV == b.central_depth_MeV &&
         a.central == b.central && a.matching_radius_fm == b.matching_radius_fm && a.step_fm == b.step_fm;
}

HybridPotential WoodsSaxonPotential(double depth_MeV, double radius_fm, double diffuseness_fm) {
  if (depth_MeV <= 0.0 || radius_fm <= 0.0 || diffuseness_fm <= 0.0)
    throw std::invalid_argument("WoodsSaxon: V0, R, a must be positive");
  HybridPotential p;
  p.central_depth_MeV = depth_MeV;
  p.central.shape = FormShape::WoodsSaxon;
  p.central.radius_fm = radius_fm;
  p.central.diffuseness_fm = diffuseness_fm;
  return p;
}

HybridPotential GaussianPotential(double depth_MeV, double width_fm) {
  if (depth_MeV <= 0.0 || width_fm <= 0.0) throw std::invalid_argument("Gaussian: V0, r0 must be positive");
  HybridPotential p;
  p.central_depth_MeV = depth_MeV;
  p.central.shape = FormShape::Gaussian;
  p.central.width_fm = width_fm;
  return p;
}

double NuclearTerm(const HybridPotential &p, double r) {
  if (p.central.shape == FormShape::Gaussian) {
    double x = r / p.central.width_fm;
    return -p.central_depth_MeV * std::exp(-x * x);
  }
  return -p.central_depth_MeV / (1.0 + std::exp((r - p.central.radius_fm) / p.central.diffuseness_fm));
}

double RadialPotential(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double r) {
  const double Rc = p.coulomb_radius_fm;
  // Coulomb: uniformly charged sphere inside R_C, point charge outside.
  double Vc;
  if (r < Rc) {
    double x = r / Rc;
    Vc = z1 * z2 * kCoulombE2_MeV_fm / (2.0 * Rc) * (3.0 - x * x);
  } else {
    Vc = z1 * z2 * kCoulombE2_MeV_fm / r;
  }
  double V_nuclear = NuclearTerm(p, r);
  double V_centrifugal = 0.0;
  if (r > 0.0) {
    double hbar2_over_2mu = (kHbarC_MeV_fm * kHbarC_MeV_fm) / (2.0 * mu_MeV);
    V_centrifugal = hbar2_over_2mu * l * (l + 1) / (r * r);
  }
  return Vc + V_nuclear + V_centrifugal;
}

std::vector<double> PotentialKey(const std::optional<HybridPotential> &p) {
  if (!p) return {0.0};
  return {1.0, (double)p->central.shape, p->central_depth_MeV, p->central.radius_fm, p->central.diffuseness_fm,
          p->central.width_fm, p->coulomb_radius_fm, p->matching_radius_fm, p->step_fm};
}

RadialValue IntegrateInward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double a_fm, double E_MeV,
                            const std::function<std::complex<double>(double r_fm)> &outer) {
  // The grid: a step close to step_fm with a on node n_a; the two nodes beyond r_m carry the outer solution.
  const int n_a = std::max(2, static_cast<int>(std::lround(a_fm / p.step_fm)));
  const double dr = a_fm / n_a;
  const int n_top = static_cast<int>(std::ceil(p.matching_radius_fm / dr));
  const int npts = std::max(n_top + 3, n_a + 3);
  // Nodes below n_a - 2 are never read (the difference below reaches n_a - 2), so neither g nor u is made there.
  const int lowest = n_a - 2;

  std::vector<std::complex<double>> u(npts);
  std::vector<double> g(npts);
  for (int i = lowest; i < npts; ++i) {
    const double r = i * dr;
    if (i >= npts - 2) u[i] = outer(r);
    // g(r) = 2 mu / hbar^2 (E - V(r)), fm^-2
    g[i] = 2.0 * mu_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm) * (E_MeV - RadialPotential(p, z1, z2, mu_MeV, l, r));
  }

  // Numerov, inward from r_m.
  const double dr2 = dr * dr;
  for (int i = npts - 3; i >= lowest; --i) {
    double denom = 1.0 + (dr2 / 12.0) * g[i];
    if (std::abs(denom) > 1e-20) {
      u[i] = ((2.0 - (5.0 * dr2 / 6.0) * g[i + 1]) * u[i + 1] - (1.0 + (dr2 / 12.0) * g[i + 2]) * u[i + 2]) / denom;
    } else {
      u[i] = u[i + 1];
    }
  }

  RadialValue v;
  v.u = u[n_a];
  v.du_dr = (u[n_a - 2] - 8.0 * u[n_a - 1] + 8.0 * u[n_a + 1] - u[n_a + 2]) / (12.0 * dr);
  return v;
}

}  // namespace physics
}  // namespace azure2
