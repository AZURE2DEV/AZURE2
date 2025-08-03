#ifndef VAVILOVAIRYFUNCTION_H
#define VAVILOVAIRYFUNCTION_H

#include <cmath>

// Vavilov Function by Airy approximation
class VavilovAiryFunction
{
 private:
 double AiryStep,AiryMinimum,AiryMaximum;
 double VAxi, VAk, VAbeta, VAdelta, VADEM;
 double VAa(unsigned int n);
 double Airy(double t);
 double VA(double delta, double xi, double beta, double k);
 double MaximumFunction(double a);
 public:
 VavilovAiryFunction(){};
 void SetAiryStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative);
 double GetAiryStep(){return AiryStep;};
 double GetAiryMinimum(){return AiryMinimum;};
 double GetAiryMaximum(){return AiryMaximum;};
 double GetValue(double delta, double xi, double beta, double k){return this->VA(delta,xi,beta,k);};
 double GetValue(double AtEnergy);
};

#endif