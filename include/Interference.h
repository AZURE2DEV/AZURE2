#ifndef INTERFERENCE_H
#define INTERFERENCE_H

#include <string>

///An AZURE \f$ l_1,l_2,l_1',l_2',J_1,J_2 \f$ combination

/*!
 * In the differential cross section formula of R-Matrix, nested inside the \f$ s,s',L \f$ sum
 * is a sum over \f$ l_1,l_2,l_1',l_2',J_1,J_2 \f$.  In the language of AZURE, these are equivalent
 * to combinations of two reaction pathways.  If the pathways are the same, the term represents the
 * actual contribution from the pathway to the cross section.  If they are different, the term 
 * represents the interference between the two.
 */

class Interference {
 public:
 Interference(int, int, double, std::string);
 Interference(int, int, double, double, std::string);
  std::string GetInterferenceType() const;
  int GetM1() const;
  int GetM2() const;
  double GetZ1Z2() const;
  double GetZ1Z2_UPOS() const;
 private:
  int m1_;
  int m2_;
  double z1z2_;
  double z1z2_upos_;
  std::string intertype_;
};

#endif
