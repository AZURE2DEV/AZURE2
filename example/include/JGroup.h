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
  JGroup(NucLine);
  JGroup(double,int);
  bool IsInRMatrix() const;
  int IsLevel(ALevel);
  int GetPi() const;
  int NumLevels() const;
  int NumChannels();
  int IsChannel(AChannel);
  double GetJ() const;
  void AddLevel(ALevel);
  void AddChannel(AChannel);
  AChannel *GetChannel(int);
  ALevel *GetLevel(int);
 private:
  bool isinrmatrix_;
  int pi_;
  double j_;
  std::vector<ALevel> levels_;
  std::vector<AChannel> channels_;
};

#endif
