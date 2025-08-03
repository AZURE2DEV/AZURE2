#include "include/PhysicsDistribution.h"

// Distribution main constructor
PhysicsDistribution::PhysicsDistribution( )
{
    ThermalDopplerEnergy = 0;
    AverageBeamResolution = 0;
}

void PhysicsDistribution::SetSigma(double Sigma)
{
 AverageBeamResolution = Sigma;
}

// Obtain the Thermal Doppler variance, that depends from energy and ion mass target
double PhysicsDistribution::GetThermalDoppler(double AtEnergy, double TargetMolarMass)
{
 double ProtonMolarMass = 1.007276;
 double RatioMolarMass = ProtonMolarMass / TargetMolarMass; //non-dimensional parameter
 AverageDopplerEnergy = 2.355 * std::sqrt(2 * RatioMolarMass * AtEnergy * ThermalDopplerEnergy); //keV
 return AverageDopplerEnergy; //Requires to store the last value on certain calculations
}

// Set the total distribution, making the convolution of the thermal distribution with the straggling, that requires the following new parameters
bool PhysicsDistribution::SetDistribution(double xi, double beta, double k, double DEM, double VEM, unsigned int Gauss, unsigned int Moyal, unsigned int Edgeworth, unsigned int Airy, unsigned int Landau, bool StrictGaussian)
{
 // The thermal distribution is always Gaussian, but if the variance are zero, then collapse to a Dirac's delta.
 double ThermalVariance = std::sqrt(AverageBeamResolution*AverageBeamResolution + AverageDopplerEnergy * AverageDopplerEnergy);
 unsigned int VarianceMode;
 if(std::abs(ThermalVariance)<1e-9)
 {
  ThermalDirac = DiracFunction();
  ThermalStep = ThermalDirac.GetDiracStep();
  ThermalMinimum = ThermalDirac.GetDiracMinimum();
  ThermalMaximum = ThermalDirac.GetDiracMaximum();
  VarianceMode = 0;
 }
 else
 {
  ThermalFunction = GaussFunction();
  ThermalFunction.SetGaussStep(0.0,ThermalVariance,Gauss,false);
  ThermalStep = ThermalFunction.GetGaussStep();
  ThermalMinimum = ThermalFunction.GetGaussMinimum();
  ThermalMaximum = ThermalFunction.GetGaussMaximum();
  VarianceMode = 6;
 }
 // The second distribution depends from the actual value of k-factor
  if (StrictGaussian) // Applies the Gaussian Distribution for k>0
  {
    if(k==0) //happens on first layer.
    {
    StraggDirac = DiracFunction();
    StraggStep = StraggDirac.GetDiracStep();
    StraggMaximum = StraggDirac.GetDiracMaximum();
    StraggMinimum = StraggDirac.GetDiracMinimum();
    PDMode = 0 + VarianceMode;
    return true;
    }
    else
    {
      //Gauss Distribution
     if (VEM <= 0) //It collapses to a Dirac's Delta
     {
     StraggDirac = DiracFunction();
     StraggStep = StraggDirac.GetDiracStep();
     StraggMaximum = StraggDirac.GetDiracMaximum();
     StraggMinimum = StraggDirac.GetDiracMinimum();
     PDMode = 0 + VarianceMode;
     return true;
     }
     else
     {
     StraggGauss = GaussFunction();
     StraggGauss.SetGaussStep(DEM,VEM,Gauss,false);
     StraggStep = StraggGauss.GetGaussStep();
     StraggMaximum = StraggGauss.GetGaussMaximum();
     StraggMinimum = StraggGauss.GetGaussMinimum();
     PDMode = 5 + VarianceMode;
     return true;
     }
    }
  }
  else
  {
   if(k==0) //happens on first layer.
   {
    StraggDirac = DiracFunction();
    StraggStep = StraggDirac.GetDiracStep();
    StraggMaximum = StraggDirac.GetDiracMaximum();
    StraggMinimum = StraggDirac.GetDiracMinimum();
    PDMode = 0 + VarianceMode;
    return true;
   }
   else if(k>0 && k<0.02) //Landau Distribution
   {
    StraggLandau = LandauFunction();
    StraggLandau.SetLandauStep(xi,beta,k,DEM,Landau,true);
    StraggStep = StraggLandau.GetLandauStep();
    StraggMaximum = StraggLandau.GetLandauMaximum();
    StraggMinimum = StraggLandau.GetLandauMinimum();
    PDMode = 1 + VarianceMode;
    return true;
   }
   else if(k>=0.02 && k<0.29) //Vavilov-Moyal Distribution
   {
    StraggMoyal = VavilovMoyalFunction();
    StraggMoyal.SetMoyalStep(xi,beta,k,DEM,Moyal,false);
    StraggStep = StraggMoyal.GetMoyalStep();
    StraggMaximum = StraggMoyal.GetMoyalMaximum();
    StraggMinimum = StraggMoyal.GetMoyalMinimum();
    PDMode = 2 + VarianceMode;
    return true;
   }
   else if(k>=0.29 && k<22.00) //Vavilov-Airy Distribution
   {
    StraggAiry = VavilovAiryFunction();
    StraggAiry.SetAiryStep(xi,beta,k,DEM,Airy,false);
    StraggStep = StraggAiry.GetAiryStep();
    StraggMaximum = StraggAiry.GetAiryMaximum();
    StraggMinimum = StraggAiry.GetAiryMinimum();
    PDMode = 3 + VarianceMode;
    return true;
   }
   else if(k>=22.00 && k<22.00) //Vavilov-Edgeworth Distribution
   {
    StraggEdgeworth = VavilovEdgeworthFunction();
    StraggEdgeworth.SetEdgeworthStep(xi,beta,k,DEM,Edgeworth,false);
    StraggStep = StraggEdgeworth.GetEdgeworthStep();
    StraggMaximum = StraggEdgeworth.GetEdgeworthMaximum();
    StraggMinimum = StraggEdgeworth.GetEdgeworthMinimum();
    PDMode = 4 + VarianceMode;
    return true;
   }
   else
   {
    //Gauss Distribution
    if (VEM <= 0) //It collapses to a Dirac's Delta
    {
     StraggDirac = DiracFunction();
     StraggStep = StraggDirac.GetDiracStep();
     StraggMaximum = StraggDirac.GetDiracMaximum();
     StraggMinimum = StraggDirac.GetDiracMinimum();
     PDMode = 0 + VarianceMode;
     return true;
    }
    else
    {
     StraggGauss = GaussFunction();
     StraggGauss.SetGaussStep(DEM,VEM,Gauss,false);
     StraggStep = StraggGauss.GetGaussStep();
     StraggMaximum = StraggGauss.GetGaussMaximum();
     StraggMinimum = StraggGauss.GetGaussMinimum();
     PDMode = 5 + VarianceMode;
     return true;
    }
   }
  }
  return false;
}


// Get the distribution at an certain point
double PhysicsDistribution::GetValue(double StraggEnergy, double ThermalEnergy)
{
  if(PDMode ==0)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggDirac.GetValue(StraggEnergy);
  if(PDMode ==1)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggLandau.GetValue(StraggEnergy);
  if(PDMode ==2)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggMoyal.GetValue(StraggEnergy);
  if(PDMode ==3)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggAiry.GetValue(StraggEnergy);
  if(PDMode ==4)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggEdgeworth.GetValue(StraggEnergy);
  if(PDMode ==5)
   return ThermalDirac.GetValue(ThermalEnergy) * StraggGauss.GetValue(StraggEnergy);
  if(PDMode ==6)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggDirac.GetValue(StraggEnergy);
  if(PDMode ==7)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggLandau.GetValue(StraggEnergy);
  if(PDMode ==8)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggMoyal.GetValue(StraggEnergy);
  if(PDMode ==9)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggAiry.GetValue(StraggEnergy);
  if(PDMode ==10)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggEdgeworth.GetValue(StraggEnergy);
  if(PDMode ==11)
   return ThermalFunction.GetValue(ThermalEnergy) * StraggGauss.GetValue(StraggEnergy);
  else
   return 0;
}