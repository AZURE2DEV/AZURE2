#ifndef MGROUP_H
#define MGROUP_H

#include "Constants.h"

///An AZURE internal reaction pathway.

/*!
 * An MGroup in AZURE represents a given entrance and exit channel through a \f$ J^\pi \f$ group. These can be visualized
 * as paths entering one row of the T-Matrix, and exiting through a column.
 */

class MGroup {
 public:
  MGroup(int, int, int);
  int GetChNum() const;
  int GetChpNum() const;
  int GetJNum() const;
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
