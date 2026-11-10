#ifndef COULFUNC_H
#define COULFUNC_H

#include <cmath>
#include <complex>
#include <memory>
#include <optional>
#include <gsl/gsl_sf_coulomb.h>
#include "Kinematics.h"
#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/HybridPotential.h"

class PPair;
class CoulombTable;

/// The return structure of the CoulFunc function class.

/*!
 * The CoulWaves structure contains both the irregular and regular solutions to the
 * Coulomb equation, as well as their derivatives with respect to \f$ \rho \f$.
 */

struct CoulWaves {
  /// Regular solution
  double F;
  /// Derivative of regular solution with respect to \f$ \rho \f$
  double dF;
  /// Irregular solution
  double G;
  /// Derivative of irregular solution with respect to \f$ \rho \f$
  double dG;
};

/// F, F', G, G' at (l, eta, rho) from physics/coulomb (docs/upgrade/coulomb.md) by the engine's algorithm, coulsm;
/// throws std::runtime_error when the evaluation fails (the library returns NaN with a status, never a number).
CoulWaves CoulombWavesAt(int l, double eta, double rho);

/*!
 * exp(i phi) of the hard-sphere phase, (G - i F)/sqrt(F^2 + G^2).  Just above a
 * charged-particle threshold (eta >~ 110) G^2 overflows and the plain formula
 * returns 0 instead of a unit phase; the scaled form (hypot) is used only then,
 * so every value the plain formula gets right is unchanged to the bit.
 */
/// The engine's Coulomb waves as physics/'s (same numbers; physics/ formulas take them).
inline azure2::physics::coulomb::CoulombWaves PhysicsWaves(const CoulWaves &c) {
  return azure2::physics::coulomb::CoulombWaves{c.F, c.dF, c.G, c.dG, 0.0, azure2::physics::coulomb::Status::ok};
}
inline std::complex<double> CoulombHardSpherePhase(const CoulWaves &c) { return azure2::physics::exp_minus_i_phi(PhysicsWaves(c)); }

/// A function class to calculate Coulomb functions for positive energy channels

/*!
 * The CoulFunc function class calculates the solutions to the Coulomb equation,
 * as well as other useful quantities such as shift functions and their energy
 * derivative and penetrabilities.
 */

class CoulFunc {
 public:
  /// Build for a particle pair, with that pair's hybrid nuclear potential if it has one.
  CoulFunc(PPair *pPair, bool useGSLFunctions);
  /// Charge number of the first particle.
  int z1() const;
  /// Charge number of the second particle.
  int z2() const;
  /// Reduced mass of the pair, in u.
  double redmass() const;
  /// Orbital angular momentum of the last evaluation.
  int lLast() const;
  /// Radius of the last evaluation, fm.
  double radiusLast() const;
  /// Centre-of-mass energy of the last evaluation, MeV.
  double energyLast() const;
  /// The last functions computed; a one-deep memo for repeated calls at the same point.
  struct CoulWaves coulLast() const;
  /// Record an evaluation as the memoized one.
  void setLast(int, double, double, CoulWaves);
  /// Coulomb functions at (l, radius fm, centre-of-mass energy MeV).
  CoulWaves operator()(int, double, double);
  /// Penetrability \f$P_l\f$ at (l, radius, energy).
  double Penetrability(int, double, double);
  /// Shift function \f$S_l\f$ at (l, radius, energy), positive energies.
  double PEShift(int, double, double);
  /// Hard-sphere phase factor exp(-i phi_l) = (G - iF)/|G + iF| at (l, radius, energy):
  /// from the Coulomb table when one is in use and covers the point, else from the waves.
  std::complex<double> HardSpherePhase(int, double, double);
  /// Energy derivative \f$dS_l/dE\f$ -- the term that makes the observed-width transformation singular when it grows too large.
  double PEShift_dE(int, double, double, double step = 1.0e-6);
  /// Second energy derivative of the shift function (numerical).
  double PEShift_d2E(int, double, double);

  /// The pair's hybrid nuclear potential (PPair::GetHybridPotential), none for plain Coulomb functions.
  const std::optional<azure2::physics::HybridPotential> &HybridPotential() const { return potential_; }
  /// Whether P, S and the phase come from the universal Coulomb table.
  bool UsesTable() const { return table_ != nullptr; }

  /// The external-capture integrals ask for F, G at the radii the adaptive quadrature chooses, the same radii for
  /// every term of one entrance l and energy.  With this on, operator() keeps the waves of each (pair, l, energy)
  /// block of the calling thread, keyed by the exact radius (at most 2^20 radii): a thread-private scratch whose
  /// values are what recomputing gives (docs/upgrade/channel_table.md section 3).
  void SetRadialScratch(bool use) { radialScratch_ = use; }

 private:
  /// The waves of the radial equation with the pair's potential (spin_orbit_potential.md section 1.1): the
  /// solution equal to H+ = G + iF beyond the matching radius, integrated inward to the radius (physics/
  /// IntegrateInward); the plain Coulomb functions at a radius beyond the matching radius.
  CoulWaves computeHybrid(int l, double radius, double energy);
  /// The point-Coulomb functions by the selected routine (physics/coulomb or GSL).
  CoulWaves PlainWaves(int l, double radius, double energy) const;

  bool useGSLFunctions_;
  bool radialScratch_ = false;
  int z1_;
  int z2_;
  int lLast_;
  double redmass_;
  /// k and eta of this pair (relativistic or not, Kinematics.h).
  Kinematics::Channel kin_;
  double radiusLast_;
  double energyLast_;
  struct CoulWaves coulLast_;

  // The universal Coulomb table (charged pairs without a nuclear potential, when
  // Config::useCoulombTable is set), with a one-deep memo of its last lookup: a
  // point asks for P, S and the phase of the same (l, radius, energy) in turn.
  const CoulombTable *table_ = nullptr;
  bool TableLookup(int l, double radius, double energy);
  /// PEShift by the exact routine (never the table): one route for every point of a difference stencil.
  double PEShiftExact(int l, double radius, double energy);
  int tLastL_ = -1;
  double tLastR_ = 0.0, tLastE_ = 0.0, tP_ = 0.0, tS_ = 0.0, tPhi_ = 0.0;
  bool tOk_ = false;

  /// The pair's hybrid nuclear potential; none: plain Coulomb functions.
  std::optional<azure2::physics::HybridPotential> potential_;
};

#endif
