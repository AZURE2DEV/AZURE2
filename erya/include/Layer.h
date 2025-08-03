#ifndef LAYER_H
#define LAYER_H

#include <cmath>

#include "include/SRIM.h"
#include "include/PhysicsDistribution.h"
#include "include/ResonanceFunction.h"

// The atomic class to handle the layers calculations
class Layer
{
private:

SRIM srim;

int Element;
double ThicknessStep;

double Xi(double E);
double Emax(double E);
double MakeBeta(double E);

public:

Layer(int LayerElement);
~Layer( );

PhysicsDistribution Distribution;
ResonanceFunction Resonance;

void setThicknessStep(double Step){ this->ThicknessStep = Step; };

void SetSigma(double Sigma);
void SetDistribution(double xi, double beta, double k, double DEM, double VEM, unsigned int Gauss, unsigned int Moyal, unsigned int Edgeworth, unsigned int Airy, unsigned int Landau, bool StrictGaussian);

void SetAtomicMass( double AtomicMass ){Resonance.SetAtomicMass(AtomicMass);};
void SetResonance(double Energy, double Width, double Peak, double min, double max);

double GetK(double E){return this->Xi(E)/this->Emax(E);};
double GetXi(double E){return this->Xi(E);};
double GetBeta(double E){return this->MakeBeta(E);};

double EvaluateSigma(double AtEnergy);
double EvaluateZiegler(double AtEnergy);
double EvaluateIntegral(double Energy);
double EvaluateYield(double Energy, double Thickness, double Steps);

double GetGVL(double E);
double GetVVL(double E);
double GetDEML(double E);
};

#endif