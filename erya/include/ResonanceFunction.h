#ifndef RESONANCEFUNCTION_H
#define RESONANCEFUNCTION_H

#include <cmath>

// Handle the Custom or built-in Briet-Wigner ressonance distribution
class ResonanceFunction
{
private:
double BWEnergy, BWWidth, BWPeak, BWmin, BWmax;
double GBWmin, GBWmax, GBAtomicMass;
double StrenghtEnergy(double Energy);
double BrietWigner(double Energy);
public:
ResonanceFunction(){};
void SetResonance(double Energy, double Width, double Peak, double min, double max);
void SetAtomicMass(double AtomicMass){GBAtomicMass = AtomicMass; return;}
double GetValue(double Energy);
double GetDomainMinimum();
double GetDomainMaximum();
};

#endif