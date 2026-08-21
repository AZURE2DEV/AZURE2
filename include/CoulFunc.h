#ifndef COULFUNC_H
#define COULFUNC_H

#include <memory>
#include <gsl/gsl_sf_coulomb.h>

class PPair;
class CoulFuncCache;
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
  /// Energy derivative \f$dS_l/dE\f$ -- the term that makes the observed-width transformation singular when it grows too large.
  double PEShift_dE(int, double, double);

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
  double radiusLast_;
  double energyLast_;
  struct CoulWaves coulLast_;

  // Hybrid method parameters
  std::shared_ptr<NuclearPotential> nuclearPotential_;
  double rMatch_;     // Matching radius for boundary conditions
  double drNumerov_;  // Numerov grid step
};

#endif
