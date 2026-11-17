#include "ScatteringLength.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Unit.h"
#include "JGroup.h"
#include "PPair.h"
#include "ParameterLabel.h"
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace ScatteringLength {

double CoulombX(const PPair *pair, double radius) {
  const double z1z2 = (double)pair->GetZ(1) * pair->GetZ(2);
  if (z1z2 <= 0.0) return 0.0;
  // x^2 = 8 Z1 Z2 e^2 mu a / hbar^2 = 8 Z1 Z2 alpha (mu c^2) a / (hbar c)
  return sqrt(8.0 * z1z2 * fstruc * pair->GetRedMass() * uconv * radius / hbarc);
}

double ZeroEnergyShift(const PPair *pair, int l, double radius) {
  const double x = CoulombX(pair, radius);
  if (x <= 0.0) return -(double)l;
  // The zero-energy irregular solution is sqrt(r) K_{2l+1}(x), x = sqrt(8 r / a_B); r d/dr = (x/2) d/dx
  // and K_nu' = -K_{nu-1} - (nu/x) K_nu give S_l(0) = -l - (x/2) K_{2l}/K_{2l+1}.  Scaled Bessel
  // functions: the ratio is unchanged and nothing overflows.
  const double kNum = gsl_sf_bessel_Kn_scaled(2 * l, x);
  const double kDen = gsl_sf_bessel_Kn_scaled(2 * l + 1, x);
  return -(double)l - 0.5 * x * kNum / kDen;
}

Coefficients FormulaCoefficients(const PPair *pair, double radius) {
  const double a = radius;
  const double x = CoulombX(pair, radius);
  if (x <= 0.0) return {-a, a};  // a0 = a (1 - M)
  // a0 = -a [ M / (x^2 K1^2) - 2 I1 / (x^2 K1) ] with K1 = e^-x K1s, I1 = e^x I1s.
  const double k1s = gsl_sf_bessel_K1_scaled(x);
  const double i1s = gsl_sf_bessel_I1_scaled(x);
  const double e2x = exp(2.0 * x);
  return {-a * e2x / (x * x * k1s * k1s), 2.0 * a * e2x * i1s / (x * x * k1s)};
}

std::vector<Channel> Channels(const CompoundSystem *compound, int pairKey) {
  std::vector<Channel> out;
  if (!compound->IsPairKey(pairKey)) return out;
  const int pairNum = compound->GetPairNumFromKey(pairKey);
  const PPair *pair = compound->GetPair(pairNum);
  if (pair->GetPType() != 0) return out;
  const double j1 = pair->GetJ(1), j2 = pair->GetJ(2);
  const int parity = pair->GetPi(1) * pair->GetPi(2);
  for (double s = fabs(j1 - j2); s <= j1 + j2 + 1.0e-6; s += 1.0) {
    // Bose/Fermi symmetry of an identical pair: l + s even, so s even for l = 0.
    if (pair->IsIdentical() && ((int)(s + 0.5)) % 2 != 0) continue;
    Channel c;
    c.pairKey = pairKey;
    c.pairNum = pairNum;
    c.s = s;
    c.parity = parity;
    c.radius = pair->GetPairRadius();
    JGroup probe(s, parity);
    probe.SetNucleus(pair->GetNucleus(), pair->GetNucleusZ(), pair->GetNucleusA());
    const int j = compound->IsJGroup(probe);
    if (j > 0) {
      const JGroup *jg = compound->GetJGroup(j);
      const int ch = jg->IsChannel(AChannel(0, s, pairNum, 'P'));
      if (ch > 0) {
        c.jGroup = j;
        c.channel = ch;
        c.radius = jg->GetChannel(ch)->GetRadius();
        c.hasLevels = jg->IsInRMatrix() && jg->NumLevels() > 0;
      }
    }
    out.push_back(c);
  }
  return out;
}

namespace {
/// Fill and invert the level matrices at the pair's threshold.
void ThresholdMatrices(const CompoundSystem *compound, const LevelState &state, const Config &configure, Unit &point, AMatrixFunc &f) {
  point.OwnSite().CalcEDependentValues(nullptr, compound, configure, nullptr);
  f.ClearMatrices();
  f.FillMatrices(&point);
  f.InvertMatrices();
}
}  // namespace

Value Evaluate(const CompoundSystem *compound, const LevelState &state, const Config &configure, const Channel &ch) {
  Value v;
  v.ch = ch;
  v.M = complex(0.0, 0.0);
  if (ch.hasLevels) {
    Unit point(0.0, 0.0, ch.pairKey, ch.pairKey, false, false, false, 0.0, 0, 0);
    AMatrixFunc f(compound, state, configure);
    ThresholdMatrices(compound, state, configure, point, f);
    v.M = f.GetUBilinear(ch.jGroup, ch.channel, ch.channel);
  }
  const Coefficients c = FormulaCoefficients(compound->GetPair(ch.pairNum), ch.radius);
  v.a0 = c.c1 * v.M + c.c0;
  return v;
}

bool Adjoint(const CompoundSystem *compound, const LevelState &state, const Config &configure, const Channel &ch, complex aBar, GradAccum &accum,
             const vector_matrix_r *shiftDeriv) {
  if (!ch.hasLevels) return true;  // a0 does not depend on any parameter
  Unit point(0.0, 0.0, ch.pairKey, ch.pairKey, false, false, false, 0.0, 0, 0);
  AMatrixFunc f(compound, state, configure);
  ThresholdMatrices(compound, state, configure, point, f);
  const Coefficients c = FormulaCoefficients(compound->GetPair(ch.pairNum), ch.radius);
  return f.BilinearAdjoint(&point, ch.jGroup, ch.channel, c.c1 * aBar, accum, shiftDeriv);
}

std::string ResolveConstraint(const CompoundSystem *compound, Constraint &c) {
  if (!compound->IsPairKey(c.pairKey)) return "no particle pair with key " + std::to_string(c.pairKey);
  const PPair *pair = compound->GetPair(compound->GetPairNumFromKey(c.pairKey));
  if (pair->GetPType() != 0) return "pair " + std::to_string(c.pairKey) + " is not a particle pair";
  std::vector<Channel> all = Channels(compound, c.pairKey);
  if (c.s < 0.0) {
    if (pair->IsIdentical()) return "the coherent combination (s = -1) is not defined for an identical pair";
    const double norm = (2.0 * pair->GetJ(1) + 1.0) * (2.0 * pair->GetJ(2) + 1.0);
    for (const Channel &ch : all) {
      c.terms.push_back(ch);
      c.weights.push_back((2.0 * ch.s + 1.0) / norm);
    }
  } else {
    for (const Channel &ch : all)
      if (fabs(ch.s - c.s) < 1.0e-6) {
        c.terms.push_back(ch);
        c.weights.push_back(1.0);
      }
    if (c.terms.empty()) {
      std::ostringstream m;
      m << "pair " << c.pairKey << " has no s-wave channel with channel spin " << c.s;
      return m.str();
    }
  }
  return "";
}

complex ConstraintValue(const CompoundSystem *compound, const LevelState &state, const Config &configure, const Constraint &c) {
  complex a(0.0, 0.0);
  for (size_t t = 0; t < c.terms.size(); t++) a += c.weights[t] * Evaluate(compound, state, configure, c.terms[t]).a0;
  return a;
}

int NumResiduals(const Constraint &c) {
  if (!c.active) return 0;
  return (c.sigma > 0.0 ? 1 : 0) + (c.imSigma > 0.0 ? 1 : 0);
}

void ConstraintResiduals(const CompoundSystem *compound, const LevelState &state, const Config &configure, const Constraint &c, std::vector<double> &r) {
  if (!c.active || (c.sigma <= 0.0 && c.imSigma <= 0.0)) return;
  complex a = ConstraintValue(compound, state, configure, c);
  if (c.sigma > 0.0) r.push_back((std::real(a) - c.value) / c.sigma);
  if (c.imSigma > 0.0) r.push_back((std::imag(a) - c.imValue) / c.imSigma);
}

bool ConstraintGradient(const CompoundSystem *compound, const LevelState &state, const Config &configure, const Constraint &c, GradAccum &accum,
                        const vector_matrix_r *shiftDeriv, int row) {
  if (!c.active || (c.sigma <= 0.0 && c.imSigma <= 0.0)) return true;
  // d(chi2) = Re(conj(aBar) da):  aBar = 2 r_re/sigma + i 2 r_im/imSigma.  One residual row:
  // aBar = 1/sigma (Re row) or i/imSigma (Im row), since Re(conj(i/s) da) = Im(da)/s.
  complex aBar(0.0, 0.0);
  if (row < 0) {
    complex a = ConstraintValue(compound, state, configure, c);
    if (c.sigma > 0.0) aBar += 2.0 * (std::real(a) - c.value) / (c.sigma * c.sigma);
    if (c.imSigma > 0.0) aBar += complex(0.0, 2.0 * (std::imag(a) - c.imValue) / (c.imSigma * c.imSigma));
  } else if (row == 0 && c.sigma > 0.0) {
    aBar = 1.0 / c.sigma;
  } else if (row < NumResiduals(c)) {
    aBar = complex(0.0, 1.0 / c.imSigma);
  } else {
    return true;
  }
  for (size_t t = 0; t < c.terms.size(); t++)
    if (!Adjoint(compound, state, configure, c.terms[t], c.weights[t] * aBar, accum, shiftDeriv)) return false;
  return true;
}

std::string Report(const CompoundSystem *compound, const LevelState &state, const Config &configure) {
  std::ostringstream out;
  out << "s-wave scattering lengths at each particle pair's threshold" << std::endl
      << "(a0 = -a [M/(x^2 K1^2) - 2 I1/(x^2 K1)], M = gamma^T A gamma at E = 0; Paneru et al., "
         "PRC 99, 045807 (2019), Eq. 7)" << std::endl;
  out << std::setw(6) << "Pair" << "  " << std::left << std::setw(14) << "Particles" << std::setw(6) << "J^pi"
      << std::right << std::setw(10) << "x" << std::setw(16) << "Re a0 (fm)" << std::setw(16) << "Im a0 (fm)"
      << "  note" << std::endl;
  for (int aa = 1; aa <= compound->NumPairs(); aa++) {
    const PPair *pair = compound->GetPair(aa);
    if (pair->GetPType() != 0) continue;
    for (const Channel &ch : Channels(compound, pair->GetPairKey())) {
      Value v = Evaluate(compound, state, configure, ch);
      out << std::setw(6) << ch.pairKey << "  " << std::left << std::setw(14) << AZURELabel::Pair(compound, aa)
          << std::setw(6) << AZURELabel::JPi(ch.s, ch.parity) << std::right << std::fixed << std::setprecision(4)
          << std::setw(10) << CoulombX(pair, ch.radius) << std::scientific << std::setprecision(6) << std::setw(16)
          << std::real(v.a0) << std::setw(16) << std::imag(v.a0);
      if (!ch.hasLevels) out << "  no level in this channel (hard sphere)";
      out << std::defaultfloat << std::endl;
    }
  }
  return out.str();
}

}  // namespace ScatteringLength
