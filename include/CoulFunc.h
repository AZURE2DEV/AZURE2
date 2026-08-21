#ifndef COULFUNC_H
#define COULFUNC_H

#include <memory>
#include <gsl/gsl_sf_coulomb.h>

class PPair;
class CoulFuncCache;
class NuclearPotential;

///The return structure of the CoulFunc function class.

/*! 
 * The CoulWaves structure contains both the irregular and regular solutions to the
 * Coulomb equation, as well as their derivatives with respect to \f$ \rho \f$.
 */

struct CoulWaves {
  ///Regular solution
  double F; 
  ///Derivative of regular solution with respect to \f$ \rho \f$
  double dF;
  ///Irregular solution
  double G;
  ///Derivative of irregular solution with respect to \f$ \rho \f$
  double dG;
};

///A function class to calculate Coulomb functions for positive energy channels

/*! 
 * The CoulFunc function class calculates the solutions to the Coulomb equation, 
 * as well as other useful quantities such as shift functions and their energy
 * derivative and penetrabilities. 
 */

class CoulFunc {
 public:
  CoulFunc(PPair *pPair, bool useGSLFunctions);
  int z1() const;
  int z2() const;
  double redmass() const;
  int lLast() const;
  double radiusLast() const;
  double energyLast() const;
  struct CoulWaves coulLast() const;
  void setLast(int, double, double, CoulWaves);
  CoulWaves operator()(int,double,double);
  double Penetrability(int,double,double);
  double PEShift(int,double,double);
  double PEShift_dE(int,double,double);

  // Hybrid method support
  void setNuclearPotential(std::shared_ptr<NuclearPotential> potential);
  std::shared_ptr<NuclearPotential> getNuclearPotential() const;
  void setMatchingRadius(double r_match);
  double getMatchingRadius() const;
  void setNumerovGridStep(double dr);
  double getNumerovGridStep() const;
  void setUseHybridMethod(bool useHybrid);
  bool getUseHybridMethod() const;
  ///The key of the pair this object was built for, which is what selects
  ///its nuclear potential out of NuclearPotentialManager.
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

  static double thisPEShift(double,void*);
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
  ///Identifies the nuclear potential these waves were computed under; part of
  ///the Coulomb memo key.  0 means the plain Coulomb solution.
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
  double rMatch_;      // Matching radius for boundary conditions
  double drNumerov_;   // Numerov grid step
};

#endif
