#ifndef KGROUP_H
#define KGROUP_H

#include "MGroup.h"
#include "ECMGroup.h"

/// An AZURE \f$ s,s' \f$ group.

/*!
 * In R-Matrix formalism, the equations required to calculate the cross section usually nested
 * inside sums over entrance and exit channel spins.  For this reason AZURE groups reaction pathways
 * according to their entrance and exit channel spins.  Each KGroup object is a container for vectors
 * of MGroup and ECMGroup objects.
 */

class KGroup {
 public:
  /// Build from an entrance and an exit channel spin.
  KGroup(double, double);
  /// Build with a second exit channel spin, for an unobserved-primary,
  /// observed-secondary (UPOS) reaction.
  KGroup(double, double, double);
  /// Number of internal (resonant) pathways.
  int NumMGroups() const;
  /// Number of external-capture pathways.
  int NumECMGroups() const;
  /// 1-based position of the internal pathway, or 0 if absent.
  int IsMGroup(MGroup);
  /// Entrance channel spin \\f$s\\f$.
  double GetS() const;
  /// Exit channel spin \\f$s'\\f$.
  double GetSp() const;
  /// Second exit channel spin, UPOS reactions only.
  double GetSp2() const;
  void AddMGroup(MGroup);
  void AddECMGroup(ECMGroup);
  /// Internal pathway \\p i, 1-based.
  MGroup *GetMGroup(int);
  /// External-capture pathway \\p i, 1-based.
  ECMGroup *GetECMGroup(int);

 private:
  double s_;
  double sp_;
  double sp2_;
  std::vector<MGroup> mgroups_;
  std::vector<ECMGroup> ec_mgroups_;
};

#endif
