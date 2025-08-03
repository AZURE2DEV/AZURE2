#ifndef GAUSSFUNCTION_H
#define GAUSSFUNCTION_H

#include <cmath>

// Gauss Normal Distribution
class GaussFunction
{
private:
double GaussDelta,GaussMinimum,GaussMaximum,GaussAverage,GaussStandardDesviation;
double FG(double Mean, double StandardDesviation, double x);
public:
GaussFunction(){};
void SetGaussStep(double Mean, double StandardDesviation, unsigned int numberstep, bool TrimNegative);
double GetGaussStep(){return GaussDelta;};
double GetGaussMinimum(){return GaussMinimum;};
double GetGaussMaximum(){return GaussMaximum;};
double GetValue(double Mean, double StandardDesviation, double x){return this->FG(Mean,StandardDesviation,x);};
double GetValue(double AtEnergy){return this->FG(GaussAverage,GaussStandardDesviation,AtEnergy);};
};

#endif