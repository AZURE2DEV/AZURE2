#include "include/ChebyshevPolynomial.h"

// Evaluate the value of a Chebyshev polynomial.
double ChebyshevPolynomial::T(unsigned int Grade, double Argument)
{
 // Recurrence Chebyshev Polynomial:
 // For n>=2: T_{n}(x) = 2 * x * T_{n-1}(x) - T_{n-2}(x)
 // And n=0 : T_{0} = 1 ; n = 1 : T_{1} = x
 if(Grade == 0)
  return 1;
 else if (Grade == 1)
  return Argument;
 else
  return 2 * Argument * this->T(Grade-1,Argument) - this->T(Grade-2,Argument);
}