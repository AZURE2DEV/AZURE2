#ifndef INTERFERENCE_H
#define INTERFERENCE_H

#include <string>

/// An AZURE \f$ l_1,l_2,l_1',l_2',J_1,J_2 \f$ combination

/*!
 * In the differential cross section formula of R-Matrix, nested inside the \f$ s,s',L \f$ sum
 * is a sum over \f$ l_1,l_2,l_1',l_2',J_1,J_2 \f$.  In the language of AZURE, these are equivalent
 * to combinations of two reaction pathways.  If the pathways are the same, the term represents the
 * actual contribution from the pathway to the cross section.  If they are different, the term
 * represents the interference between the two.
 */

class Interference {
 public:
  /// Build from two pathway positions, their \\f$Z_1Z_2\\f$ coefficient and the type.
  Interference(int, int, double, std::string);
  /// As above with a second coefficient, for UPOS reactions.
  Interference(int, int, double, double, std::string);
  /// "RR", "ER", "RE" or "EE" -- which vector each index refers to,
  /// MGroup (R, resonant) or ECMGroup (E, external).
  std::string GetInterferenceType() const;
  /// 1-based position of the first pathway, in the vector the type names.
  int GetM1() const;
  /// 1-based position of the second pathway.
  int GetM2() const;
  /// Angular-distribution coefficient \\f$Z_1Z_2\\f$ for this pathway pair.
  double GetZ1Z2() const;
  /// The same coefficient for UPOS reactions.
  double GetZ1Z2_UPOS() const;

 private:
  int m1_;
  int m2_;
  double z1z2_;
  double z1z2_upos_;
  std::string intertype_;
};

#endif
