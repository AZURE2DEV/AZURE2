#include "include/GaussFunction.h"

double GaussFunction::FG(double Mean, double StandardDesviation, double x)
{
 double pi = 2*std::acos(0.0);
 double arg = (x - Mean)/StandardDesviation;
 double f = (1/std::sqrt(2 * pi * StandardDesviation * StandardDesviation)) * std::exp(-1 * arg * arg / 2);
 return f;
}

void GaussFunction::SetGaussStep(double Mean, double StandardDesviation, unsigned int numberstep, bool TrimNegative)
{
 GaussDelta = (6 * StandardDesviation)/numberstep;
 GaussMinimum = Mean - 3 * StandardDesviation;
 GaussMaximum = Mean + 3 * StandardDesviation;
 GaussAverage = Mean;
 GaussStandardDesviation = StandardDesviation;
 // Cut negative minimum values
 if(GaussMinimum < 0 && TrimNegative)
 {
  GaussMinimum = 0;
 }
 return;
}