#ifndef MGROUP_H
#define MGROUP_H

#include "Constants.h"

/// An AZURE internal reaction pathway.

/*!
 * An MGroup in AZURE represents a given entrance and exit channel through a \f$ J^\pi \f$ group. These can be visualized
 * as paths entering one row of the T-Matrix, and exiting through a column.
 */

class MGroup {
 public:
  /// Build from entrance channel, exit channel and J-group positions.
  MGroup(int, int, int);
  /// 1-based entrance channel, within the J-group's AChannel vector.
  int GetChNum() const;
  /// 1-based exit channel, within the same vector.
  int GetChpNum() const;
  /// 1-based J-group this pathway runs through.
  int GetJNum() const;
  /// Statistical spin factor \\f$g_J = (2J+1)/[(2I_1+1)(2I_2+1)]\\f$.
  double GetStatSpinFactor() const;
  void SetStatSpinFactor(double);

 private:
  int jnum_;
  int ch_;
  int chp_;
  double statspinfactor_;
  complex tmatrix_;
};

#endif
