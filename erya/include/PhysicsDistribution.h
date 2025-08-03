#ifndef PHYSICSDISTRIBUTION_H
#define PHYSICSDISTRIBUTION_H

#include <cmath>

#include "include/DiracFunction.h"
#include "include/GaussFunction.h"
#include "include/LandauFunction.h"
#include "include/VavilovEdgeworthFunction.h"
#include "include/VavilovAiryFunction.h"
#include "include/VavilovMoyalFunction.h"

// Handle the distributions related to the straggling, thermal Doppler and Beam Resolution
class PhysicsDistribution
{
private:
double AverageBeamResolution, ThermalDopplerEnergy, AverageDopplerEnergy;
double ThermalStep,ThermalMinimum,ThermalMaximum,StraggStep,StraggMinimum,StraggMaximum;
GaussFunction ThermalFunction,StraggGauss;
LandauFunction StraggLandau;
VavilovEdgeworthFunction StraggEdgeworth;
VavilovAiryFunction StraggAiry;
VavilovMoyalFunction StraggMoyal;
DiracFunction StraggDirac,ThermalDirac;
unsigned int PDMode;
bool IsDefined;
public:
PhysicsDistribution();
void SetSigma(double Sigma);
double GetBeamResolution(){return AverageBeamResolution;}
double GetThermalDoppler(double AtEnergy, double TargetMolarMass);
double GetThermalVariance(){return std::sqrt(AverageBeamResolution*AverageBeamResolution + AverageDopplerEnergy * AverageDopplerEnergy);};
bool CheckValidity(){return IsDefined;}
bool SetDistribution(double xi, double beta, double k, double DEM, double VEM, unsigned int Gauss, unsigned int Moyal, unsigned int Edgeworth, unsigned int Airy, unsigned int Landau, bool StrictGaussian);
double GetStraggStep(){return StraggStep;};
double GetStraggMin(){return StraggMinimum;};
double GetStraggMax(){return StraggMaximum;};
double GetThermalStep(){return ThermalStep;};
double GetThermalMin(){return ThermalMinimum;};
double GetThermalMax(){return ThermalMaximum;};
double GetValue(double StraggEnergy, double ThermalEnergy);
};

#endif