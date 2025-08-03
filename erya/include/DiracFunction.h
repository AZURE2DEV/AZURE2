#ifndef DIRACFUNCTION_H
#define DIRACFUNCTION_H

#include <cmath>

// Dirac Delta Distribution
class DiracFunction
{
 private:
 double Dirac(double x);
 public:
 DiracFunction(){};
 double GetValue(double x){return this->Dirac(x);}
 double GetDiracStep(){return 1;};
 double GetDiracMinimum(){return 0;};
 double GetDiracMaximum(){return 0;};
};

#endif