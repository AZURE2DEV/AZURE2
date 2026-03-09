#include "MGroup.h"

/*!
 * This constructor is used to create an MGroup object with reference to positions in the JGroup and subsequent AChannel vectors.
 */

MGroup::MGroup(int jGroupNum, int channelNum, int channelPrimeNum) :
  jnum_(jGroupNum), ch_(channelNum), chp_(channelPrimeNum), statspinfactor_(0.0) {
}

/*!
 * Returns the position of the entrance channel in the AChannel vector below the corresponding JGroup object.
 */

int MGroup::GetChNum() const {
  return ch_;
}

/*!
 * Returns the position of the exit channel in the AChannel vector below the corresponding JGroup object.
 */

int MGroup::GetChpNum() const {
  return chp_;
}

/*! 
 * Returns the position of the \f$ J^\pi \f$ group in the JGroup vector.
 */

int MGroup::GetJNum() const {
  return jnum_;
}

/*!
 * Returns the statistical spin factor, \f$ g_J \f$, for the reaction pathway.
 */

double MGroup::GetStatSpinFactor() const {
  return statspinfactor_;
}

/*!
 * Sets the statistical spin factor, \f$ g_J \f$, for the reaction pathway.
 */

void MGroup::SetStatSpinFactor(double spinFactor) {
  statspinfactor_=spinFactor;
}
