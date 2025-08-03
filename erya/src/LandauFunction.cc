#include "include/LandauFunction.h"

// Evaluate the value of the Laudau function
double LandauFunction::L(double lambda)
{
// Wilkinson interpolation formula.
 if(lambda<=-3.4)
 {
  return 0;
 }
 else if (lambda > -3.4 && lambda <= -1)
 {
  double f = (1/std::sqrt(4*std::acos(0)))*std::exp((std::abs(lambda)-1)/2 - exp(std::abs(lambda)-1));
  double g = 1 + 0.01*(6.7853 + 4.884*lambda + 1.4488*std::pow(lambda,2) + 0.20802*std::pow(lambda,3) + 0.012057*std::pow(lambda,4));
  return (f*g > 0) ? f*g : 0.0;
 }
 else if (lambda >-1 && lambda <=3)
 {
  double f = 0.17885481 - 0.015464468 * lambda - 0.030040482 * std::pow(lambda,2) + 0.013781358 * std::pow(lambda,3) - 0.0011491556 * std::pow(lambda,4) - 0.0012835837 * std::pow(lambda,5) + 0.00062395162 * std::pow(lambda,6) - 0.0001262613 * std::pow(lambda,7) + 0.000010108918 * std::pow(lambda,8);
  return (f > 0) ? f : 0.0;
 }
 else if (lambda > 3 && lambda <= 150)
 {
  double x = std::log(lambda);
  double f = std::exp(-1.5669876 - 1.5811676 * x + 1.677088 * std::pow(x,2) - 1.4972908 * std::pow(x,3) + 0.57062974 * std::pow(x,4) - 0.11777036 * std::pow(x,5) + 0.01343737 * std::pow(x,6) - 0.00076315158 * std::pow(x,7) + 0.000014881108 * std::pow(x,8));
  return (f > 0) ? f : 0.0;
 }
 else if (lambda > 150)
 {
  double x = std::log(lambda);
  double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
  double g = std::pow(lambda,-2) - ( 3 - 2 * euler - 2 * x) * std::pow(lambda,-3);
  double f = 1 - 0.01 * std::exp(5.157 - 1.42 * x);
  return (g/f > 0) ? g/f : 0.0;
 }
 else
 {
  return 0; //No valid argument
 }
}

// Define the default function domain
void LandauFunction::SetLandauStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative)
{
 // Function domain in terms of lambda
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 double lambdalow = -3.4;
 double lambdahigh = 25;
 double lambdastep = (lambdahigh - lambdalow)/numberstep;
 // Function domain in terms of energy, using the Landau's variable E = xi(lambda - lambda_M) + DEM , where lambda_M = -1 + euler - beta^2 - log(k)
 LandauStep = xi * lambdastep;
 LandauMinimum = xi * (lambdalow + 1 - euler + beta*beta + std::log(k)) + DEM;
 LandauMaximum = xi * (lambdahigh + 1 - euler + beta*beta + std::log(k)) + DEM;
 LandauEnergy = DEM;
 LandauXi = xi;
 LandauBeta = beta;
 LandauK = k;
 // Cut energy loss greater than average
 if((LandauMinimum < 0) && TrimNegative)
 {
  LandauMinimum = 0;
 }
 return;
}

double LandauFunction::GetValue(double AtEnergy)
{
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 double EnergyLambda = (AtEnergy - LandauEnergy)/LandauXi - 1 + euler - LandauBeta*LandauBeta - std::log(LandauK);
 return this->L(EnergyLambda);
}