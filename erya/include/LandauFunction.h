#ifndef LANDAUFUNCTION_H
#define LANDAUFUNCTION_H

#include <cmath>

// Laudau Function
class LandauFunction
{
private:
double LandauStep,LandauMinimum,LandauMaximum,LandauEnergy,LandauXi,LandauBeta,LandauK;
double L(double lambda);
public:
LandauFunction(){};
void SetLandauStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative);
double GetLandauStep(){return LandauStep;};
double GetLandauMinimum(){return LandauMinimum;};
double GetLandauMaximum(){return LandauMaximum;};
double GetValue(double AtEnergy);
};

#endif