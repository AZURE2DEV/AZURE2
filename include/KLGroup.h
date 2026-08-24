#ifndef KLGROUP_H
#define KLGROUP_H

#include <vector>
#include <map>
#include <tuple>
#include <string>
#include "Interference.h"

/// An AZURE \f$ s,s',L \f$ group

/*!
 * Differential cross sections in R-Matrix theory contains terms nested inside a sum over entrance and exit spins
 * as well as Legendre polynomial orders, \f$ L \f$.  In AZURE, an \f$ s,s' \f$ combination is given by a KGroup object.
 * It is therefore convenient to group KGroup objects with a specified polynomial orders for the calculation of differenial
 * cross sections.  The KLGroup object serves as a container class for a vector of Interference objects.
 */

class KLGroup {
 public:
  /// Build from a KGroup position and a Legendre order.
  KLGroup(int, int);
  /// 1-based position of the \\f$s,s'\\f$ group this belongs to.
  int GetK() const;
  /// Order of the Legendre polynomial this term multiplies.
  int GetLOrder() const;
  int NumInterferences() const;
  /// 1-based position of the pathway pair, or 0 if absent.
  int IsInterference(Interference);
  void AddInterference(Interference);
  /// Pathway pair \\p i, 1-based.
  Interference *GetInterference(int);

 private:
  int k_;
  int lorder_;
  std::vector<Interference> interferences_;
  // Index of (m1,m2,type) -> 1-based position in interferences_, so that
  // IsInterference() is O(log N) instead of a linear scan.  Without this the
  // angular-distribution setup is O(M^4) in the number of pathways and stalls
  // for external-capture channels with many EC pathways (e.g. capture to
  // excited states).
  std::map<std::tuple<int, int, std::string>, int> interferenceIndex_;
};

#endif
