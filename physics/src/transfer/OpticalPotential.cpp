#include "azure2/physics/transfer/OpticalPotential.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "azure2/physics/HybridPotential.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

/// The Woods-Saxon form factor f(r; R, a) of physics/HybridPotential.h (NuclearTerm of a unit-depth well), so that
/// the library has one implementation of it.
double WoodsSaxonForm(double radius_fm, double diffuseness_fm, double r_fm) {
  HybridPotential unit;
  unit.central_depth_MeV = 1.0;
  unit.central.shape = FormShape::WoodsSaxon;
  unit.central.radius_fm = radius_fm;
  unit.central.diffuseness_fm = diffuseness_fm;
  return -NuclearTerm(unit, r_fm);
}

/// V_C of a uniformly charged sphere of radius R_C (0: point charge) for the charge product z1 z2, by
/// physics/HybridPotential.h's RadialPotential with no nuclear term and l = 0.
double CoulombSphere(int z1z2, double rc_fm, double r_fm) {
  if (z1z2 == 0) return 0.0;
  HybridPotential coulomb_only;
  coulomb_only.coulomb_radius_fm = rc_fm;
  coulomb_only.central_depth_MeV = 0.0;
  coulomb_only.central.shape = FormShape::WoodsSaxon;
  coulomb_only.central.radius_fm = 1.0;
  coulomb_only.central.diffuseness_fm = 1.0;
  return RadialPotential(coulomb_only, z1z2, 1, 1.0, 0, 0.0, r_fm);
}

void CheckGeometry(const WoodsSaxonTerm &t, double R, const char *name) {
  if (t.depth_MeV != 0.0 && (!(R > 0.0) || !(t.diffuseness_fm > 0.0)))
    throw std::invalid_argument(std::string("optical potential: the ") + name +
                                " term has a depth but no positive radius and diffuseness");
}

}  // namespace

bool operator==(const WoodsSaxonTerm &a, const WoodsSaxonTerm &b) {
  return a.depth_MeV == b.depth_MeV && a.r0_fm == b.r0_fm && a.diffuseness_fm == b.diffuseness_fm;
}

bool operator==(const OpticalPotential &a, const OpticalPotential &b) {
  return a.convention == b.convention && a.target_mass_number == b.target_mass_number &&
         a.projectile_mass_number == b.projectile_mass_number && a.coulomb_r0_fm == b.coulomb_r0_fm &&
         a.real_volume == b.real_volume && a.imaginary_volume == b.imaginary_volume &&
         a.imaginary_surface == b.imaginary_surface && a.real_spin_orbit == b.real_spin_orbit &&
         a.imaginary_spin_orbit == b.imaginary_spin_orbit;
}

double RadiusOf(const OpticalPotential &p, double r0_fm) {
  switch (p.convention) {
    case RadiusConvention::Absolute:
      return r0_fm;
    case RadiusConvention::TargetMass:
      if (!(p.target_mass_number > 0.0)) throw std::invalid_argument("optical potential: A_t must be positive");
      return r0_fm * std::cbrt(p.target_mass_number);
    case RadiusConvention::BothMasses:
      if (!(p.target_mass_number > 0.0) || !(p.projectile_mass_number > 0.0))
        throw std::invalid_argument("optical potential: A_t and A_p must be positive");
      return r0_fm * (std::cbrt(p.target_mass_number) + std::cbrt(p.projectile_mass_number));
  }
  return r0_fm;
}

PotentialTerms TermsAt(const OpticalPotential &p, int z1z2, double r) {
  PotentialTerms t;
  t.coulomb_MeV = CoulombSphere(z1z2, RadiusOf(p, p.coulomb_r0_fm), r);
  auto form = [&](const WoodsSaxonTerm &term, const char *name) {
    const double R = RadiusOf(p, term.r0_fm);
    CheckGeometry(term, R, name);
    return term.depth_MeV == 0.0 ? 0.0 : WoodsSaxonForm(R, term.diffuseness_fm, r);
  };
  t.volume_MeV = {-p.real_volume.depth_MeV * form(p.real_volume, "real volume"),
                  -p.imaginary_volume.depth_MeV * form(p.imaginary_volume, "imaginary volume")};
  const double fd = form(p.imaginary_surface, "imaginary surface");
  t.surface_MeV = {0.0, -4.0 * p.imaginary_surface.depth_MeV * fd * (1.0 - fd)};
  auto so = [&](const WoodsSaxonTerm &term, const char *name) {
    const double R = RadiusOf(p, term.r0_fm);
    CheckGeometry(term, R, name);
    if (term.depth_MeV == 0.0) return 0.0;
    SpinOrbitTerm s;
    s.depth_MeV = term.depth_MeV;
    s.form.shape = FormShape::WoodsSaxon;
    s.form.radius_fm = R;
    s.form.diffuseness_fm = term.diffuseness_fm;
    return SpinOrbitRadialTerm(s, r);
  };
  t.spin_orbit_form_MeV = {so(p.real_spin_orbit, "real spin-orbit"), so(p.imaginary_spin_orbit, "imaginary spin-orbit")};
  return t;
}

std::complex<double> PotentialAt(const OpticalPotential &p, int z1z2, double two_l_dot_s, double r) {
  const PotentialTerms t = TermsAt(p, z1z2, r);
  std::complex<double> u = t.coulomb_MeV + t.volume_MeV + t.surface_MeV;
  if (two_l_dot_s != 0.0) u += t.spin_orbit_form_MeV * two_l_dot_s;
  return u;
}

std::string Describe(const OpticalPotential &p) {
  std::ostringstream s;
  s.precision(6);
  const char *conv = p.convention == RadiusConvention::Absolute     ? "R = r0"
                     : p.convention == RadiusConvention::TargetMass ? "R = r0 At^1/3"
                                                                      : "R = r0 (Ap^1/3 + At^1/3)";
  s << "radii " << conv << " (At = " << p.target_mass_number << ", Ap = " << p.projectile_mass_number << ")\n";
  s << "Coulomb sphere   Rc = " << RadiusOf(p, p.coulomb_r0_fm) << " fm\n";
  auto line = [&](const char *name, const WoodsSaxonTerm &t) {
    if (t.depth_MeV == 0.0) return;
    s << name << " depth = " << t.depth_MeV << " MeV  R = " << RadiusOf(p, t.r0_fm) << " fm  a = " << t.diffuseness_fm
      << " fm\n";
  };
  line("real volume     ", p.real_volume);
  line("imag volume     ", p.imaginary_volume);
  line("imag surface    ", p.imaginary_surface);
  line("real spin-orbit ", p.real_spin_orbit);
  line("imag spin-orbit ", p.imaginary_spin_orbit);
  return s.str();
}

OpticalPotential NoPotential() { return OpticalPotential{}; }

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
