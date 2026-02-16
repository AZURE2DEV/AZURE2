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
