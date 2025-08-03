#include "include/VavilovEdgeworthFunction.h"
#include "include/HermitePolynomial.h"

// Auxiliary private functions
double VavilovEdgeworthFunction::VEalpha(unsigned int n)
{
 double alpha = (std::pow(VExi,n)/std::pow(VEk,n-1)) * ( (1 / ( n - 1)) - (VEbeta * VEbeta) / n);
 return alpha;
}

double VavilovEdgeworthFunction::VEmu(unsigned int n)
{
 if (n==0)
  return 1;
 else if (n==1)
  return 0;
 else if (n==2)
  return this->VEalpha(2);
 else if (n==3)
  return this->VEalpha(3);
 else if (n==4)
  return 3 * this->VEalpha(2) * this->VEalpha(2) + this->VEalpha(4);
 else if (n==5)
  return 10 * this->VEalpha(2) * this->VEalpha(3) + this->VEalpha(5);
 else
  return 0;
}

// Evaluate the value of the Vavilov distribution by the Edgeworth approximation
double VavilovEdgeworthFunction::VE(double delta, double xi, double beta, double k)
{
 // Initialize global variables
 double pi = 2 * std::acos(0);
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 VExi = xi;
 VEk = k;
 VEbeta = beta;
 VEdelta = delta;
 VEStdDev = std::sqrt(this->VEmu(2));
 // Evaluate the function
 HermitePolynomial h;
 double gauss = (1 / std::sqrt(2 * pi * VEStdDev * VEStdDev)) * std::exp((-1 * VEdelta * VEdelta) / (2 * VEStdDev * VEStdDev));
 double edgeworth3 = (1 / 6) * (this->VEmu(3) / std::pow(VEStdDev,3)) * h.GetValue(3,VEdelta/VEStdDev);
 double edgeworth4 = (1 / 24) * (this->VEmu(4) / std::pow(VEStdDev,4) - 3) * h.GetValue(4,VEdelta/VEStdDev);
 double edgeworth5 = (1 / 120) * ((this->VEmu(5) / std::pow(VEStdDev,5)) - 10 * (this->VEmu(3) / std::pow(VEStdDev,3))) * h.GetValue(5,VEdelta/VEStdDev);
 double edgeworth6 = (1 / 72) * std::pow(this->VEmu(3) / std::pow(VEStdDev,3),2) * h.GetValue(6,VEdelta/VEStdDev);
 double edgeworth7 = (1 / 144) * (this->VEmu(3) / std::pow(VEStdDev,3)) * (this->VEmu(4) / std::pow(VEStdDev,4) - 3) * h.GetValue(7,VEdelta/VEStdDev);
 double edgeworth9 = (1 / 1296) * std::pow( this->VEmu(3) / std::pow(VEStdDev,3) , 3) * h.GetValue(9,VEdelta/VEStdDev);
 double edgeworth = 1 + edgeworth3 + edgeworth4 + edgeworth5 + edgeworth6 + edgeworth7 + edgeworth9;
 double f = gauss * edgeworth;
 return (f > 0) ? f : 0.0;
}

// Define the distribution domain
void VavilovEdgeworthFunction::SetEdgeworthStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative)
{
 // The integral domain may be evaluated from k and beta parameters
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 VExi = xi;
 VEk = k;
 VEbeta = beta;
 VEDEM = DEM;
 double lowlambda = (-0.0322 * VEbeta * VEbeta - 0.0743)*VEk + (-0.2453 * VEbeta * VEbeta + 0.0701)/std::sqrt(VEk) + (-0.5561 * VEbeta * VEbeta - 3.1579);
 double highlambda = (-0.0135 * VEbeta * VEbeta - 0.0488)*VEk + (-1.6921 * VEbeta * VEbeta + 8.3656)/std::sqrt(VEk) + (-0.7327 * VEbeta * VEbeta - 3.5226);
 double lambdastep = (highlambda - lowlambda)/numberstep;
 // Define the Function Domain
 EdgeworthStep = VExi * lambdastep;
 EdgeworthMinimum = VExi * (lowlambda + 1 - euler + VEbeta*VEbeta + std::log(VEk)) + DEM;
 EdgeworthMaximum = VExi * (highlambda + 1 - euler + VEbeta*VEbeta + std::log(VEk)) + DEM;
 // Cut energy loss greater than average
 if((EdgeworthMinimum < 0) && TrimNegative)
 {
  EdgeworthMinimum = 0;
 }
 return;
}

double VavilovEdgeworthFunction::GetValue(double AtEnergy)
{
 VEdelta = AtEnergy - VEDEM;
 return this->VE(VEdelta,VExi,VEbeta,VEk);
}