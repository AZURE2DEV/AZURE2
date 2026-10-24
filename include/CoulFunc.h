#ifndef COULFUNC_H
#define COULFUNC_H

#include <cmath>
#include <complex>
#include <memory>
#include <gsl/gsl_sf_coulomb.h>
#include "Kinematics.h"

class PPair;
class CoulFuncCache;
class CoulombTable;
class NuclearPotential;

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

/*!
 * exp(i phi) of the hard-sphere phase, (G - i F)/sqrt(F^2 + G^2).  Just above a
 * charged-particle threshold (eta >~ 110) G^2 overflows and the plain formula
 * returns 0 instead of a unit phase; the scaled form (hypot) is used only then,
 * so every value the plain formula gets right is unchanged to the bit.
 */
inline std::complex<double> CoulombHardSpherePhase(const CoulWaves &c) {
  const double s = std::sqrt(std::pow(c.F, 2.0) + std::pow(c.G, 2.0));
  const double re = c.G / s, im = -c.F / s;
  if (std::isfinite(re) && std::isfinite(im) && std::isfinite(s) && s > 0.0) return std::complex<double>(re, im);
  const double h = std::hypot(c.F, c.G);
  return std::complex<double>(c.G / h, -c.F / h);
}

/// A function class to calculate Coulomb functions for positive energy channels

/*!
 * The CoulFunc function class calculates the solutions to the Coulomb equation,
 * as well as other useful quantities such as shift functions and their energy
 * derivative and penetrabilities.
 */

class CoulFunc {
 public:
  /// Build for a particle pair; resolves that pair's hybrid nuclear potential.
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
  double PEShift_dE(int, double, double);
  /// Second energy derivative of the shift function (numerical).
  double PEShift_d2E(int, double, double);

  // Hybrid method support
  /// Override this instance's nuclear potential.
  void setNuclearPotential(std::shared_ptr<NuclearPotential> potential);
  /// The nuclear potential in force for this pair.
  std::shared_ptr<NuclearPotential> getNuclearPotential() const;
  /// Radius at which the Numerov solution is matched to the Coulomb functions.
  void setMatchingRadius(double r_match);
  /// Matching radius, fm.
  double getMatchingRadius() const;
  /// Step of the outward Numerov integration, fm.
  void setNumerovGridStep(double dr);
  /// Numerov step, fm.
  double getNumerovGridStep() const;
  /// Turn the hybrid model on or off for this instance.
  void setUseHybridMethod(bool useHybrid);
  /// Is the hybrid model active here?
  bool getUseHybridMethod() const;
  /// The key of the pair this object was built for, which is what selects
  /// its nuclear potential out of NuclearPotentialManager.
  int pairKey() const { return pairKey_; }
  /// The nuclear potential these waves are computed under (0 = plain Coulomb);
  /// part of the ECTable term key.
  long hybridTag() const { return hybridTag_; }
  /// Whether P, S and the phase come from the universal Coulomb table.
  bool UsesTable() const { return table_ != nullptr; }

  // The global (radius-keyed, mutex-protected) Coulomb cache only helps when the
  // same radius is queried across many energies (e.g. penetrabilities at the
  // channel radius).  For external-capture integrals the radius is the
  // integration variable, so the cache never hits and its per-call mutex lock
  // serializes the OpenMP threads.  Disable it for those instances.
  void SetUseGlobalCache(bool use) { useGlobalCache_ = use; }
  bool GetUseGlobalCache() const { return useGlobalCache_; }

 private:
  // Hybrid method calculation
  CoulWaves computeHybrid(int l, double radius, double energy);

  static double thisPEShift(double, void *);
  static double thisPEShift_dE(double, void *);
  typedef struct DEShiftParams {
    CoulFunc *coulFunc;
    int lValue;
    double radius;
  } DEShiftParams;
  DEShiftParams dEShiftParams_;
  bool useGSLFunctions_;
  bool useHybridMethod_;
  bool useGlobalCache_ = true;
  int pairKey_;
  /// Identifies the nuclear potential these waves were computed under; part of
  /// the Coulomb memo key.  0 means the plain Coulomb solution.
  long hybridTag_;
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
  int tLastL_ = -1;
  double tLastR_ = 0.0, tLastE_ = 0.0, tP_ = 0.0, tS_ = 0.0, tPhi_ = 0.0;
  bool tOk_ = false;

  // Hybrid method parameters
  std::shared_ptr<NuclearPotential> nuclearPotential_;
  double rMatch_;     // Matching radius for boundary conditions
  double drNumerov_;  // Numerov grid step
};

#endif
