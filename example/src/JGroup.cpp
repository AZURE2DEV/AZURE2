#include "JGroup.h"
#include "NucLine.h"

/*!
 * This constructor is used when a \f$ J^\pi \f$ group is created from an entry in the nuclear input file.
 */

JGroup::JGroup(NucLine nucLine):
  pi_(nucLine.levelPi()), j_(nucLine.levelJ()), isinrmatrix_(true) {};
  
/*!
 * This constructor is used when a \f$ J^\pi \f$ group is created from specified values of spin and parity.
 */

JGroup::JGroup(double j,int pi):
  pi_(pi),j_(j),isinrmatrix_(false) {};

/*!
 * Returns true if the \f$ J^\pi \f$ group is to be included in the A-/R-Matrix calculation, otherwise returns false. 
 * A \f$ J^\pi \f$ group may specify only a bound state for external capture, but may not correspond to an R-Matrix state (i.e. subthreshold state).
 */

bool JGroup::IsInRMatrix() const {
  return isinrmatrix_;
}

/*!
 * This function tests if a given level already exists in the vector of ALevel objects.  If the level exists 
 * the position of the level in the vector is returned, otherwise the function returns 0.
 */

int JGroup::IsLevel(ALevel level) {
  bool b=false;
  int c=0;
  double tol=1e-3;
  while(!b&&c<this->NumLevels())
    {
      if(this->GetLevel(c+1)->GetE()-tol<=level.GetE()&&
	 level.GetE()<=this->GetLevel(c+1)->GetE()+tol) b=true;
      c++;
    }
  if(b) return c;
  else return 0;
}

/*!
 * Returns the parity of the the \f$ J^\pi \f$ group as \f$ \pm1 \f$.
 */

int JGroup::GetPi() const {
  return pi_;
}

/*!
 * Returns the number of levels in the ALevel vector.
 */

int JGroup::NumLevels() const {
  return levels_.size();
}

/*!
 * Returns the number of channels in the AChannel vector.
 */

int JGroup::NumChannels() {
  return channels_.size();
}

/*!
 * This function tests if a given channel already exists in the vector of AChannel objects.  If the channel exists 
 * the position of the channel in the vector is returned, otherwise the function returns 0.
 */

int JGroup::IsChannel(AChannel channel) {
  bool b=false;
  int c=0;
  while(!b&&c<this->NumChannels()) {
    if(channel.GetL()==this->GetChannel(c+1)->GetL()&&
       channel.GetS()==this->GetChannel(c+1)->GetS()&&
       channel.GetPairNum()==this->GetChannel(c+1)->GetPairNum()) b=true;
    c++;
  }
  if(b) return c;
  else return 0;
}

/*!
 * Returns the spin value of the \f$ J^\pi \f$ group.
 */

double JGroup::GetJ() const {
  return j_;
}

/*!
 * Adds a new level to the vector of ALevel objects.
 */

void JGroup::AddLevel(ALevel level) {
  levels_.push_back(level);
}

/*!
 * Adds a new channel to the vector of AChannel objects.
 */

void JGroup::AddChannel(AChannel channel) {
  channels_.push_back(channel);
}

/*! 
 * Returns a pointer to a specified channel in the AChannel vector.
 */


AChannel *JGroup::GetChannel(int channelNum) {
  AChannel *b=&channels_[channelNum-1];
  return b;
}

/*! 
 * Returns a pointer to a specified level in the ALevel vector.
 */

ALevel *JGroup::GetLevel(int levelNum) {
  ALevel *b=&levels_[levelNum-1];
  return b;
}
