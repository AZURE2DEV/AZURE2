#ifndef JGROUP_H
#define JGROUP_H

#include "Constants.h"
#include "ALevel.h"
#include "AChannel.h"

class NucLine;

/// An AZURE \f$ J^\pi \f$ group.

/*!
 * In R-Matrix theory, levels are grouped according to their \f$ J^\pi \f$ values.  There is one R-/A-Matrix,
 * and thus one T-Matrix, for each \f$ J^\pi \f$ group. A JGroup object holds vectors of ALevel and AChannel objects.
 */

class JGroup {
 public:
  /// Build from a line of the nuclear input file.
  JGroup(NucLine);
  /// Build directly from a spin and a parity (\f$\pm1\f$).
  JGroup(double, int);
  /// Is this group part of the A-/R-Matrix calculation?
  /// False for a group that exists only to carry a bound state for external
  /// capture -- a subthreshold state has no R-Matrix level of its own.
  bool IsInRMatrix() const;
  /// 1-based position of the level in the group, or 0 if it is not there.
  int IsLevel(ALevel);
  /// Parity of the group, \f$\pm1\f$.
  int GetPi() const;
  /// Number of levels in the group.
  int NumLevels() const;
  /// Number of channels shared by every level of the group.
  int NumChannels();
  /// 1-based position of the channel in the group, or 0 if it is not there.
  int IsChannel(AChannel);
  /// Total angular momentum of the group.
  double GetJ() const;
  /// Append a level.
  void AddLevel(ALevel);
  /// Append a channel; it applies to every level in the group.
  void AddChannel(AChannel);
  /// Channel \p i, 1-based. Every level of the group shares this channel set.
  AChannel *GetChannel(int);
  /// Level \p i, 1-based, numbered within this group rather than globally.
  ALevel *GetLevel(int);

 private:
  bool isinrmatrix_;
  int pi_;
  double j_;
  std::vector<ALevel> levels_;
  std::vector<AChannel> channels_;
};

#endif
