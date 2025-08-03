#ifndef VAVILOVMOYALFUNCTION_h
#define VAVILOVMOYALFUNCTION_h

#include <cmath>

// Vavilov Function by Moyal approximation
class VavilovMoyalFunction
{
private:
double MoyalStep,MoyalMinimum,MoyalMaximum;
double VMCk, VMCbeta, VMDEM, VMxi;
double VMTable(unsigned int z, unsigned int i, unsigned int m, unsigned int n);
double VMbeta(double beta);
double VMk(double k);
double VMa(unsigned int i, double k, double beta);
double VMTwoSigma(double k, double beta, double lambda);
double VMThreeSigma(double k, double beta, double lambda);
double VMMain(double k, double beta, double lambda);
void MakeTable();
double VMT[3][13][10];
public:
VavilovMoyalFunction(){this->MakeTable();};
void SetMoyalStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative);
double GetMoyalStep(){return MoyalStep;};
double GetMoyalMinimum(){return MoyalMinimum;};
double GetMoyalMaximum(){return MoyalMaximum;};
double GetValue(double k, double beta, double lambda){return this->VMMain(k,beta,lambda);};
double GetValue(double AtEnergy);
};

#endif