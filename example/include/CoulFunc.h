#ifndef COULFUNC_H
#define COULFUNC_H

class PPair;

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
  double PEShift_dE(int,double,double);;
 private:
  static double thisPEShift(double,void*);
  typedef struct DEShiftParams {
    CoulFunc *coulFunc;
    int lValue;
    double radius;
  } DEShiftParams;
  DEShiftParams dEShiftParams_;
  bool useGSLFunctions_;
  int z1_;
  int z2_;
  int lLast_;
  double redmass_;
  double radiusLast_;
  double energyLast_;
  struct CoulWaves coulLast_;
};

#endif
