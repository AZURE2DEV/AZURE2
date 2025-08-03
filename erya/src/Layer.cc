#include "include/Layer.h"

Layer::Layer(int LayerElement)
{
 Element = LayerElement;
 ThicknessStep = 10; 
}

Layer::~Layer( )
{
}

//Evaluate the k and x factors
double Layer::MakeBeta(double E)
{
 double mp = 938272; //proton's mass in keV
 return std::sqrt(2 * mp * E + E * E) / (mp + E);
}

double Layer::Xi(double E)
{
 //Xi = 2.5507e-7 * Z * M * n / (A * b^2) [keV]
 double XiFactor = 2.5507e-7;
 double TotalCharge = Element;
 double TotalMass = srim.getElementMass(Element);
 double MolarMass = srim.getElementMass(Element);
 double Beta = this->MakeBeta(E);
 double f = (XiFactor * TotalCharge * MolarMass * this->ThicknessStep)/(TotalMass * Beta * Beta);
 return f;
}

double Layer::Emax(double E)
{
 //E = (2 * m * c^2 * b^2 * g^2)/(1 + 2*g*m/M)
 double me = 511; //electron's mass in keV
 double mp = 938272; //proton's mass in keV
 double Beta = this->MakeBeta(E);
 double Gamma = 1/std::sqrt(1-Beta*Beta);
 double f = (2 * me * Beta * Beta * Gamma * Gamma) / (1 + (2*Gamma*me)/mp);
 return f;
}

//Get the Bohr variance of the current layer
double Layer::GetGVL(double E)
{
 double BohrFactor = (8 * this->GetXi(E)) / (3);
 double Beta = this->MakeBeta(E);
 double me = 511; //electron mass in keV
 double SumIonization = srim.getElementIonization(Element) * std::log(2 * me * Beta * Beta / srim.getElementIonization(Element));
 return std::sqrt(BohrFactor * SumIonization);
}

//Get the Vavilov variance of the current layer
double Layer::GetVVL(double E)
{
 double Xi = this->GetXi(E);
 double Beta = this->GetBeta(E);
 double K = this->GetK(E);
 return std::sqrt(Xi*Xi*(1-Beta*Beta/2)/K);
}

double Layer::GetDEML(double E)
{
 return this->EvaluateZiegler(E) * this->ThicknessStep * 0.001;
}

void Layer::SetSigma(double Sigma)
{
   Distribution.SetSigma(Sigma);
}

// Set the total distribution, making the convolution of the thermal distribution with the straggling, that requires the following new parameters
void Layer::SetDistribution(double xi, double beta, double k, double DEM, double VEM, unsigned int Gauss, unsigned int Moyal, unsigned int Edgeworth, unsigned int Airy, unsigned int Landau, bool StrictGaussian)
{
   Distribution.SetDistribution(xi,beta,k,DEM,VEM,Gauss,Moyal,Edgeworth,Airy,Landau,StrictGaussian);
}

void Layer::SetResonance(double Energy, double Width, double Peak, double min, double max)
{
 Resonance.SetResonance(Energy,Width,Peak,min,max);
}

double Layer::EvaluateSigma(double Energy)
{
   if(Resonance.GetDomainMinimum() <= Energy && Energy <= Resonance.GetDomainMaximum()) return Resonance.GetValue(Energy);
   else return 0;
}

// Evaluate the Ziegler function
double Layer::EvaluateZiegler(double AtEnergy)
{
 if ( AtEnergy == 0)
 {
    return 0;
 }
 else if ( AtEnergy >= 0 && AtEnergy < 25)  // Simplification to get Ziegler(0)=0
 {
    return srim.getElementValue1( Element )*std::pow(AtEnergy,srim.getElementValue2( Element )) + srim.getElementValue3( Element )*std::pow(AtEnergy,srim.getElementValue4( Element ));
 }
 else
 {
    double StoppingLow, StoppingHigh, Stopping;
    StoppingLow = srim.getElementValue1( Element )*std::pow(AtEnergy,srim.getElementValue2( Element )) + srim.getElementValue3( Element )*std::pow(AtEnergy,srim.getElementValue4( Element ));
    StoppingHigh = (srim.getElementValue5( Element )/std::pow(AtEnergy,srim.getElementValue6( Element )))*std::log( (srim.getElementValue7( Element )/AtEnergy) + (srim.getElementValue8( Element )*AtEnergy));
    Stopping = (StoppingLow*StoppingHigh)/(StoppingHigh+StoppingLow);
    return Stopping;
 }
}

double Layer::EvaluateIntegral(double Energy)
{
 //Implements a double numerical integration using the Simpson method
 double DDT = this->Distribution.GetThermalStep();
 double DDS = this->Distribution.GetStraggStep();
 double DTmin = this->Distribution.GetThermalMin();
 double DTmax = this->Distribution.GetThermalMax();
 double DSmin = this->Distribution.GetStraggMin();
 double DSmax = this->Distribution.GetStraggMax();
 // Set the number of integration steps
 unsigned int Tsteps = std::abs(DTmax-DTmin)<1e-9 ? 0 : std::floor(std::abs((DTmax-DTmin)/(DDT)));
 unsigned int Ssteps = std::abs(DSmax-DSmin)<1e-9 ? 0 : std::floor(std::abs((DSmax-DSmin)/(DDS)));
 // Fixes the convolution integration domain
 double Tmin = DTmin + DSmin;
 double Tmax = DTmax + DSmax;
 double Smin = DTmin + DSmin;
 double Smax = DTmax + DSmax;
 double DT = std::abs(DTmax-DTmin)<1e-9 ? DDT : DDT*(Tmax - Tmin)/(DTmax-DTmin);
 double DS = std::abs(DSmax-DSmin)<1e-9 ? DDS : DDS*(Smax - Smin)/(DSmax-DSmin);
 if(Tsteps==0)
 {
  Tmin = 0;
  Tmax = 0;
 }
 if(Ssteps==0)
 {
  Smin = 0;
  Smax = 0;
 }
 // Evaluate the integral itself, including the renormalization integral.
 double DoubleSimpsonSum = 0;
 double RenormalizationSum = 0;
 for(unsigned int i=0; i<=Ssteps; i++)
 {
   for(unsigned int j=0; j<=Tsteps; j++)
   {
     double S = Smin + i * DS ;
     double T = Tmin + j * DT ;
     double LocalStoppingPower = this->EvaluateZiegler(Energy-S);
     double LocalDistribution = Distribution.GetValue(S-T,T);

     DoubleSimpsonSum = DoubleSimpsonSum + LocalDistribution * this->EvaluateSigma(Energy-S) / LocalStoppingPower;
     RenormalizationSum = RenormalizationSum + LocalDistribution;
   }
 }
 double CrossSigmaSum = (DS * DT) * DoubleSimpsonSum;
 double DistributionSum = (DS * DT) * RenormalizationSum;
 // If the renormalization itself are zero, return zero, since the first integral will also the zero.
 if(DistributionSum == 0.0)
   return 0;
 else
   return CrossSigmaSum / DistributionSum;
}

double Layer::EvaluateYield( double Energy, double Thickness, double Steps )
{
 double Yield = 0;
 double Step = Thickness / Steps;
 for(unsigned int i=0; i<=Steps; i++)
 {
   double LocalEnergy = Energy - i * Step;
   Yield = Yield + this->EvaluateIntegral(LocalEnergy);
 }
 return Yield;
}