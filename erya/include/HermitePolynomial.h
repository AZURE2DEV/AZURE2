#ifndef HERMITEPOLYNOMIAL_H
#define HERMITEPOLYNOMIAL_H

// Hermite polynomial
class HermitePolynomial
{
private:
double H(unsigned int Grade, double Argument);
public:
HermitePolynomial(){};
double GetValue(unsigned int Grade, double Argument){return this->H(Grade,Argument);};
};

#endif