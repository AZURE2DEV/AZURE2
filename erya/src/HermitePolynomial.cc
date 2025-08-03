#include "include/HermitePolynomial.h"

// Evaluate the value of a Hermite polynomial.
double HermitePolynomial::H(unsigned int Grade, double Argument)
{
 // Recurrence Hermite Polynomial:
 // For n>=2: H_{n}(x) = x * H_{n-1}(x) - (n-1) * H_{n-2}(x)
 // And n=0 : H_{0} = 1 ; n = 1 : H_{1} = x
 if(Grade == 0)
  return 1;
 else if (Grade == 1)
  return Argument;
 else
  return Argument * this->H(Grade-1,Argument) - (Grade - 1 ) * this->H(Grade-2,Argument);
}