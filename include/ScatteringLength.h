#ifndef SCATTERINGLENGTH_H
#define SCATTERINGLENGTH_H

#include "Constants.h"
#include <string>
#include <vector>

class CNuc;
class LevelState;
class Config;
class PPair;
struct GradAccum;

/// s-wave scattering lengths from the level matrix at a particle pair's threshold.
/*!
 * Paneru et al., Phys. Rev. C 99, 045807 (2019), Eq. (7) (also Phys. Rev. C 109, 015802 (2024), Eq. (26)):
 * \f[ a_0 = -a\left[\frac{M_{cc}}{x^2K_1^2(x)} - \frac{2I_1(x)}{x^2K_1(x)}\right],\quad
 *     M_{cc} = \tilde\gamma_c^T\tilde A\,\tilde\gamma_c\ \mathrm{at}\ E = 0,\quad
 *     x = (8Z_1Z_2e^2\mu a/\hbar^2)^{1/2}, \f]
 * and \f$a_0 = a(1-M_{cc})\f$ for a neutral pair.  The definition is the Coulomb-modified effective-range
 * expansion \f$C_0^2k\cot\delta_0 + 2k\eta h(\eta) = -1/a_0 + \ldots\f$.  M is evaluated with the entrance
 * channel at exactly zero energy: P = 0 and \f$S_l(0) = -l - (x/2)K_{2l}(x)/K_{2l+1}(x)\f$
 * (ZeroEnergyShift); the other channels of the J-group take their values at that total energy.  An
 * open channel (a lower threshold) makes M and a0 complex.  The result does not depend on the
 * boundary condition and is the same in the standard, Brune and Park parametrizations.
 */
namespace ScatteringLength {

/// Channels within this energy (MeV) of their threshold take the exact zero-energy limit.
const double kThresholdTolerance = 1.0e-10;

/// \f$x = (8Z_1Z_2e^2\mu a/\hbar^2)^{1/2}\f$ at the pair's channel radius (0 for a neutral pair).
double CoulombX(const PPair *pair);

/// Shift function \f$S_l\f$ of the pair at zero channel energy: \f$-l - (x/2)K_{2l}(x)/K_{2l+1}(x)\f$,
/// \f$-l\f$ for a neutral pair.
double ZeroEnergyShift(const PPair *pair, int l);

/// a0 = c1 M + c0 (fm), the coefficients of the closed formula for the pair.
struct Coefficients {
  double c1;
  double c0;
};
Coefficients FormulaCoefficients(const PPair *pair);

/// One s-wave channel spin of a particle pair.
struct Channel {
  int pairKey = 0;     ///< pair key as written in the .azr
  int pairNum = 0;     ///< position in CNuc
  double s = 0.0;      ///< channel spin = J for l = 0
  int parity = 1;      ///< pi1 pi2
  int jGroup = 0;      ///< J-group holding the channel, 0 if none
  int channel = 0;     ///< channel index in that J-group, 0 if none
  bool hasLevels = false;  ///< false: no level has a width in this channel, M = 0 (hard sphere)
};

/// Every allowed s-wave channel spin of a particle pair (identical pairs: s even only).
std::vector<Channel> Channels(CNuc *compound, int pairKey);

/// Value of one channel: a0 (fm, complex if a channel is open at the threshold) and M.
struct Value {
  Channel ch;
  complex a0;
  complex M;
};

/// Evaluate one channel with the compound's current parameters (Brune and Park need
/// CNuc::CalcShiftFunctions first, as for any calculation).
Value Evaluate(CNuc *compound, const LevelState &state, const Config &configure, const Channel &ch);

/*!
 * Accumulate Re(conj(aBar) d a0/d parameters) into accum (level energies and widths), exactly.
 * Returns false when the level-matrix adjoint does not support the configuration (RMC formalism,
 * levels left out of the R-matrix); the caller then differentiates numerically.
 */
bool Adjoint(CNuc *compound, const LevelState &state, const Config &configure, const Channel &ch, complex aBar, GradAccum &accum,
             const vector_matrix_r *shiftDeriv);

/// A measured scattering length: a <modelData> datum of quantity scattering_length (ModelData.cpp).
struct Constraint {
  bool active = true;
  int pairKey = 0;
  double s = 0.0;           ///< channel spin; < 0: spin-weighted (coherent) sum over the channel spins
  double value = 0.0;       ///< Re a (fm)
  double sigma = 0.0;       ///< uncertainty of Re a (fm); 0: Re not constrained
  double imValue = 0.0;     ///< Im a (fm)
  double imSigma = 0.0;     ///< uncertainty of Im a (fm); 0: Im not constrained
  int line = 0;             ///< line number in the block, for messages
  std::vector<Channel> terms;   ///< the channel spin(s) and their weights
  std::vector<double> weights;
};

/// The channel spin(s) and weights of a constraint from its pair key and s (s < 0: the coherent combination), as the
/// <modelData> reader (ModelData.cpp) resolves them.  "" or the reason it cannot be.
std::string ResolveConstraint(CNuc *compound, Constraint &c);

/// Model value of a constraint (fm).
complex ConstraintValue(CNuc *compound, const LevelState &state, const Config &configure, const Constraint &c);

/// Residuals of one constraint, appended to r: (Re a - value)/sigma and (Im a - imValue)/imSigma for each
/// nonzero sigma; none when it is inactive.
void ConstraintResiduals(CNuc *compound, const LevelState &state, const Config &configure, const Constraint &c, std::vector<double> &r);

/*!
 * Gradient of one constraint's chi2 (row = -1) or of its row-th residual (0: the first of Re, Im that it has) into
 * accum.  Returns false where Adjoint does.
 */
bool ConstraintGradient(CNuc *compound, const LevelState &state, const Config &configure, const Constraint &c, GradAccum &accum,
                        const vector_matrix_r *shiftDeriv, int row);

/// Number of residual rows of one constraint (0, 1 or 2).
int NumResiduals(const Constraint &c);

/// Table of every particle pair's s-wave scattering lengths, as written to scatteringLengths.out (measured ones are
/// <modelData> data: ModelData::Report).
std::string Report(CNuc *compound, const LevelState &state, const Config &configure);

}  // namespace ScatteringLength

#endif
