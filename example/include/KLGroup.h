#ifndef KLGROUP_H
#define KLGROUP_H

#include <vector>
#include "Interference.h"

///An AZURE \f$ s,s',L \f$ group

/*!
 * Differential cross sections in R-Matrix theory contains terms nested inside a sum over entrance and exit spins
 * as well as Legendre polynomial orders, \f$ L \f$.  In AZURE, an \f$ s,s' \f$ combination is given by a KGroup object.
 * It is therefore convenient to group KGroup objects with a specified polynomial orders for the calculation of differenial
 * cross sections.  The KLGroup object serves as a container class for a vector of Interference objects.
 */

class KLGroup {
 public:
  KLGroup(int,int);
  int GetK() const;
  int GetLOrder() const;
  int NumInterferences() const;
  int IsInterference(Interference);
  void AddInterference(Interference);
  Interference *GetInterference(int);
 private:
  int k_;
  int lorder_;
  std::vector<Interference> interferences_;
};

#endif
