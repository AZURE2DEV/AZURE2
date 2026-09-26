#ifndef CHANNELFUNC_H
#define CHANNELFUNC_H

#include <memory>

class PPair;
class CoulFunc;
class ShftFunc;

/*!
 * Shift function, penetrability and dS/dE of a particle channel at a channel
 * energy e (MeV, measured from the pair's threshold, sepE + exE), continuous
 * through threshold.  Every place that evaluates them at a level energy, or at
 * any energy that may sit at or near threshold, goes through here, so the
 * threshold handling is written once:
 *
 *   e < 0             S from the Whittaker function (ShftFunc), P = 0
 *   e = 0             S = S_l(0), the exact zero-energy limit, P = 0
 *   e > 0, eta > 100  "Coulomb threshold": P < exp(-2 pi 100) ~ 1e-273 is zero
 *                     in double precision and the Coulomb-wave routines fail
 *                     (NaN, garbage or seconds per call), so P = 0 and
 *                     S = 2 S(0) - S(-e), exact to O(e^2)
 *   otherwise         the Coulomb functions (CoulFunc)
 *
 * dS/dE is the central difference AZURE2 always used (step 1e-6 MeV) of the
 * branch the stencil lies in; a stencil that crosses threshold differences the
 * continuous S above instead of calling a routine outside its domain (the
 * Coulomb functions at e <= 0, or ShftFunc, which takes |E - threshold| and so
 * mirrored the probe above threshold back below it: dS/dE came out 0 for a level
 * exactly at threshold).
 *
 * Neutral pairs never reach the Coulomb-threshold branch.  The radius is the
 * pair's channel radius, as at every call site.
 */
class ChannelFunc {
 public:
  ChannelFunc(PPair *pair, bool useGSLFunctions);
  ~ChannelFunc();

  /// e > 0 and a charged pair with Sommerfeld parameter eta above kEtaClosed.
  static bool IsCoulombThreshold(PPair *pair, double e);
  /// Closed to double precision: e <= 0, or a Coulomb threshold.  P = 0 there.
  bool IsClosed(double e) const;

  /// Shift function S_l(e), continuous through threshold.
  double Shift(int l, double e);
  /// Penetrability P_l(e); 0 wherever IsClosed(e).
  double Penetrability(int l, double e);
  /// dS_l/dE at e (MeV^-1).
  double ShiftDerivative(int l, double e);

  /// Sommerfeld parameter above which a positive-energy channel counts as closed.
  static constexpr double kEtaClosed = 100.0;
  /// Step of the central differences (MeV), as CoulFunc/ShftFunc use.
  static constexpr double kStep = 1.0e-6;

 private:
  CoulFunc &coul();
  ShftFunc &shft();
  static double shiftAdaptor(double e, void *p);

  PPair *pair_;
  bool useGSL_;
  double threshold_;
  double radius_;
  int lDeriv_ = 0;
  std::unique_ptr<CoulFunc> coul_;
  std::unique_ptr<ShftFunc> shft_;
};

#endif
