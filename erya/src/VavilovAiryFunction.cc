#include "include/VavilovAiryFunction.h"

// Solve numerically the equation Ai'(tp)/Ai(tp)+a=0, which is equivalent to find the maximum value of log(Ai(tp))+tp*a=0, where tp are the optimal value for the maximum of Vavilov-Airy Distribution
double VavilovAiryFunction::MaximumFunction(double a)
{
 // The method adopts the Golden-Search Algorithm which are adequate for f'(x)/f(x) kind of functions
 double dx = 1e-3;
 double x1,x2,xl,xu,d;
 // Initial Values
 xl = -2.338107; // First Airy Root
 xu = a*a; // Gaussian condition: tp - a*a = 0
 d = (std::sqrt(5)-1)*(xu-xl)/2;
 x1 = xl + d;
 x2 = xu - d;
 // Cycle condition
 while(std::abs(xu-xl)<dx)
 {
  if((std::log(this->Airy(x1))+a*x1) > (std::log(this->Airy(x2))+a*x2)) // Recalculate x1
  {
    xl = x2;
    x2 = x1;
    xu = xu;
    x1 = xl + (std::sqrt(5)-1)*(xu-xl)/2;
  }
  else // Recalculate x2
  {
    xl = xl;
    xu = x1;
    x1 = x2;
    x2 = xu - (std::sqrt(5)-1)*(xu-xl)/2;
  }
 }
 // Found a solution
 return (xu+xl)/2;
}

// Set the Airy Distribution domain, using the same procedure of Laudau variable on Edgeworth function
void VavilovAiryFunction::SetAiryStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative)
{
 // The integral domain may be evaluated from k and beta parameters
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 VAxi = xi;
 VAk = k;
 VAbeta = beta;
 VADEM = DEM;
 double lowlambda = (-0.0322 * VAbeta * VAbeta - 0.0743)*VAk + (-0.2453 * VAbeta * VAbeta + 0.0701)/std::sqrt(VAk) + (-0.5561 * VAbeta * VAbeta - 3.1579);
 double highlambda = (-0.0135 * VAbeta * VAbeta - 0.0488)*VAk + (-1.6921 * VAbeta * VAbeta + 8.3656)/std::sqrt(VAk) + (-0.7327 * VAbeta * VAbeta - 3.5226);
 double lambdastep = (highlambda - lowlambda)/numberstep;
 // For lower k values, the Edgeworth optimal interval may truncate the left side of distribution. To fix this, a new lambda value should be evauated from the first zero of Airy's function.
 double t0 = -2.33811;
 double eta = (xi * std::pow((1 - (2*VAbeta*VAbeta)/3),1.0/3.0))/(std::pow(2.0*VAk,2.0/3.0));
 double a = (std::pow(2.0*VAk,1.0/3.0) * (1 - (VAbeta*VAbeta)/2)) / (std::pow((1 - (2.0*VAbeta*VAbeta)/3.0),2.0/3.0));
 double deltafix = eta*(t0-a*a);
 double tp = this->MaximumFunction(a);
 double deltamax = eta*(tp-a*a);
 // Define the Function Domain
 AiryMinimum = DEM + deltafix;
 AiryMaximum = DEM - deltamax - deltafix;
 AiryStep = 1.0 / numberstep;
 // Cut negative minimum values
 if((AiryMinimum - DEM) < 0 && TrimNegative)
 {
  AiryMinimum = DEM;
 }
 return;
}

// Evaluate the Airy-Ai function between the first negative zero (-2.33811) to all positive values. Around t=5 it will use the taylor expansion, and beyond that an assymptotic exponential function

double VavilovAiryFunction::VAa(unsigned int n)
{
 if (n == 0)
  return 1 / (std::pow(3.0,2.0/3.0) * std::tgamma(2.0/3.0));
 else if (n==1)
  return -1 / (std::pow(3.0,1.0/3.0) * std::tgamma(1.0/3.0));
 else if (n==2)
  return 0;
 else
  return (this->VAa(n-3))/(1.0*n*(1.0*n-1));
}

double VavilovAiryFunction::Airy(double t)
{
 // Apply a domain cut-off
 if ((t>-2.33811) && (t<=5))
 {
   double x=0;
   for(unsigned int n=0; n<123; n++)
   {
     x = x + this->VAa(n)*std::pow(t,n*1.0);
   }
  return x;
 }
 else if (t>5)
 {
  return ((std::exp((-2.0/3.0)*std::pow(t,3.0/2.0)))/(std::sqrt(16.0*std::atan(1.0))*std::pow(t,1.0/4.0)))*(1-5/(48*std::pow(t,3.0/2.0)));
 }
 else
 {
  return 0;
 }
}

// Evaluate the Airy approximation value of Vavilov distribution
double VavilovAiryFunction::VA(double delta, double xi, double beta, double k)
{
 // Set the following model parameters
 double eta = (xi * std::pow((1 - (2*VAbeta*VAbeta)/3),1.0/3.0))/(std::pow(2.0*VAk,2.0/3.0));
 double a = (std::pow(2.0*VAk,1.0/3.0) * (1 - (VAbeta*VAbeta)/2)) / (std::pow((1 - (2.0*VAbeta*VAbeta)/3.0),2.0/3.0));
 double t = delta / eta + a*a;
 double f = this->Airy(t) * std::exp(a*t-(a*a*a)/3) / eta;
 return (f>0) ? f : 0.0;
}


// Return the distribution value in terms of energy
double VavilovAiryFunction::GetValue(double AtEnergy)
{
 VAdelta = AtEnergy - VADEM;
 return this->VA(VAdelta,VAxi,VAbeta,VAk)*VAxi;
}