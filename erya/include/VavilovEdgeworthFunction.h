#ifndef VAVILOVEDGEWORTHFUNCTION_H
#define VAVILOVEDGEWORTHFUNCTION_H

#include <cmath>

// Vavilov Function by Edgeworth approximation
class VavilovEdgeworthFunction
{
private:
double EdgeworthStep,EdgeworthMinimum,EdgeworthMaximum;
double VE(double delta, double xi, double beta, double k);
double VEalpha(unsigned int n);
double VEmu(unsigned int n);
double VEStdDev, VExi, VEk, VEbeta, VEdelta, VEDEM;
public:
VavilovEdgeworthFunction(){};
void SetEdgeworthStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative);
double GetEdgeworthStep(){return EdgeworthStep;};
double GetEdgeworthMinimum(){return EdgeworthMinimum;};
double GetEdgeworthMaximum(){return EdgeworthMaximum;};
double GetValue(double delta, double xi, double beta, double k){return this->VE(delta,xi,beta,k);};
double GetValue(double AtEnergy);
};

#endif