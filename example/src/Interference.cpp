#include "Interference.h"

/*!
 * The pathways combination is created specifically using references to two positions
 * in the MGroup and ECMGroup vectors under the corresponding KGroup object.  Additionally,
 * the \f$ Z_1 Z_2 \f$ coefficients are passed along with the interference type.
 * The interference type is either RR, ER, RE,or EE, indicating which vector, the MGroup or ECMGroup, the 
 * stored indices refer to.
 */

Interference::Interference(int mGroupNum1, int mGroupNum2, double z1z2Coeff, std::string interferenceType) :
  m1_(mGroupNum1), m2_(mGroupNum2),z1z2_(z1z2Coeff),intertype_(interferenceType) {};

/*!
 * Overloaded to accommidate additional z1z2 coefficient for unobserved primary, observed secondar
 * reactions.
 */

Interference::Interference(int mGroupNum1, int mGroupNum2, double z1z2Coeff, double z1z2Coeff_upos, std::string interferenceType) :
  m1_(mGroupNum1), m2_(mGroupNum2), z1z2_(z1z2Coeff), z1z2_upos_(z1z2Coeff_upos), intertype_(interferenceType) {};

/*!
 * Returns the interference type.
 */

std::string Interference::GetInterferenceType() const {
  return intertype_;
}

/*!
 * Returns the position in the MGroup or ECMGroup vector of the first pathway.
 */

int Interference::GetM1() const {
  return m1_;
}

/*!
 * Returns the position in the MGroup or ECMGroup vector of the second pathway.
 */

int Interference::GetM2() const {
  return m2_;
}

/*!
 * Returns the corresponding \f$ Z_1 Z_2 \f$ coefficient.
 */

double Interference::GetZ1Z2() const {
  return z1z2_;
}

/*!
 * Returns the corresponding \f$ Z_1 Z_2 \f$ coefficient for unobserved primary observed secondary reactions.
 */

double Interference::GetZ1Z2_UPOS() const {
  return z1z2_upos_;
}

