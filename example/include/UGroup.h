#ifndef UGROUP_H
#define UGROUP_H

#include "CGroup.h"

/// An AZURE \f$ s,s1',s2' \f$ group.

/*!
 * In R-Matrix formalism, the equations required to calculate the cross section usually nested 
 * inside sums over entrance and exit channel spins.  For this reason AZURE groups reaction pathways 
 * according to their entrance and exit channel spins.  Each UGroup object is a container for vectors
 * of CGroup objects.  
 */

class UGroup {
 public:
  UGroup(double, double);
  int NumCGroups() const;
  int IsCGroup(CGroup);
  double GetS() const;
  double GetS1p() const;
  double GetS2p() const;
  void AddCGroup(CGroup);
  CGroup *GetCGroup(int);
 private:
  double s_;
  double s1p_;
  double s2p_;
  std::vector<MGroup> cgroups_;
};

#endif
