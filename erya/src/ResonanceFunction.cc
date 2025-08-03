#include "include/ResonanceFunction.h"

void ResonanceFunction::SetResonance(double Energy, double Width, double Peak, double min, double max)
{
 BWEnergy = Energy;
 BWWidth = Width;
 BWPeak = Peak;
 BWmin = min;
 BWmax = max;
 GBWmin = min;
 GBWmax = max;
 return;
}

//Evaluate the Briet-Wigner at a certain energy.
// f(E) = (peak * width^2 / 4 ) / (width^2 / 4 + (E-energy)^2)
double ResonanceFunction::BrietWigner(double Energy)
{
 if(GBWmin <= Energy && Energy <= GBWmax)
 {
    double BWSum = 0;
    BWSum = BWSum +  (BWPeak * BWWidth * BWWidth / 4 ) / (BWWidth * BWWidth / 4 + std::pow(Energy - BWEnergy,2));
    return BWSum;
 }
 else
 {
   return 0;
 }
}

//Evaluate the Briet-Wigner with a Resonance Strenght at a certain energy.
// f(E) = (K * strenght * width / 4 )/ (width^2 / 4 + (E-energy)^2)
// where K = 2607472.5 (A+1)/(A*E) mili-barn
double ResonanceFunction::StrenghtEnergy(double Energy)
{
 if(GBWmin <= Energy && Energy <= GBWmax)
 {
    double BWSum = 0;
    BWSum = BWSum + ((BWPeak * BWWidth * 2607472.5 * (GBAtomicMass + 1) ) / (4 * Energy * GBAtomicMass) ) / (BWWidth * BWWidth / 4 + std::pow(Energy - BWEnergy,2));
    return BWSum;
 }
 else
 {
   return 0;
 }
}

//Get the Ressonance function value
double ResonanceFunction::GetValue(double Energy)
{
  return this->StrenghtEnergy(Energy);
}

double ResonanceFunction::GetDomainMinimum()
{
  return GBWmin;
}

double ResonanceFunction::GetDomainMaximum()
{
  return GBWmax;
}