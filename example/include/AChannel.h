#ifndef ACHANNEL_H
#define ACHANNEL_H

class NucLine;

///An AZURE channel object.

/*!
 * An R-Matrix channel for a given \f$J^\pi\f$ group represents a specfic combination of \f$ \alpha,s,l \f$ couplings.
 */

class AChannel {
 public:
  AChannel(NucLine, int);
  AChannel(int, double, int, char);
  int GetPairNum() const;
  int GetL() const;
  double GetS() const;
  double GetBoundaryCondition() const;
  char GetRadType() const;
  void SetBoundaryCondition(double);
 private:
  int l_;
  int pair_;
  double s_;
  char radtype_;
  double boundary_condition_;
};

#endif
