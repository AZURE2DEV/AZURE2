#include "AChannel.h"
#include "NucLine.h"
#include <math.h>
#include <assert.h>


/*!
 * This constructor can be used if a channel is to be created from a specific line of the nuclear input file.
 */

AChannel::AChannel(NucLine nucLine, int pairNum) {
  l_=nucLine.l();
  s_=nucLine.s();
  pair_=pairNum;
  if(nucLine.pType()==0) {
    radtype_='P';
  } else if(nucLine.pType()==10) {
    if(nucLine.levelPi()*nucLine.pi2()==pow(-1,nucLine.l())) radtype_='E';
    else radtype_='M';
  } else if(nucLine.pType()==20) {
    if(nucLine.l()==0) radtype_='F';
    else radtype_='G';
  }
}

/*!
 * This constructor can be used if a channel is to be created directly using specified channel couplings.
 */

AChannel::AChannel(int lValue, double sValue, int pairNum, char radType) {
  l_=lValue;
  s_=sValue;
  pair_=pairNum;
  radtype_=radType;
};

/*!
 * Returns the pair number, or position of the corresponding particle pair in the PPair vector, for the channel.
 */

int AChannel::GetPairNum() const {
  return pair_;
}

/*!
 * Returns the orbital angular momentum for the channel.
 */

int AChannel::GetL() const {
  return l_;
}

/*!
 * Returns the coupled channel spin for the channel.
 */

double AChannel::GetS() const {
  return s_;
}

/*!
 * Returns the boundary condition for the channel.  The energy of the boundary condition is fixed at the first level given in the nuclear input file for a given \f$ J^\pi \f$ group.
 */

double AChannel::GetBoundaryCondition() const {
  return boundary_condition_;
}

/*!
 * Returns the radiation type for the channel. Radiation types are:
 * - P: Particle Radiation
 * - E: EL Electromagnetic Radiation
 * - M: ML Electromagnetic Radiation
 */

char AChannel::GetRadType() const {
  return radtype_;
}

/*!
 * This function is used to set the boundary condition for the channel.
 */
void AChannel::SetBoundaryCondition(double boundaryCondition) {
  boundary_condition_=boundaryCondition;
}

