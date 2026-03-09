#ifndef ULGROUP_H
#define ULGROUP_H

#include <vector>
#include "Interference.h"

///An AZURE \f$ s,s1',s2',L \f$ group

/*!
 * Differential cross sections for unobserved particles in R-Matrix theory contain terms nested inside a sum over entrance 
 * and exit spins as well as Legendre polynomial orders, \f$ L \f$.  In AZURE, an \f$ s,s1',s2' \f$ combination is given by a 
 * UGroup object. It is therefore convenient to group UGroup objects with a specified polynomial orders for the calculation 
 * of differenial cross sections. The ULGroup object serves as a container class for a vector of Interference objects.
 */

class ULGroup {
 public:
  ULGroup(int,int);
  int GetU() const;
  int GetLOrder() const;
  int NumInterferences() const;
  int IsInterference(Interference);
  void AddInterference(Interference);
  Interference *GetInterference(int);
 private:
  int u_;
  int lorder_;
  std::vector<Interference> interferences_;
};

#endif
