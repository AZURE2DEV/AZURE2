#include "include/DiracFunction.h"

// Dirac's delta function gives one when x=0
double DiracFunction::Dirac(double x)
{
 if(std::abs(x) < 1e-9)
 {
   return 1;
 }
 else
 {
   return 0;
 }
}