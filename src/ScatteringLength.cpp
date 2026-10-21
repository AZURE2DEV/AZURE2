#include "ScatteringLength.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "JGroup.h"
#include "PPair.h"
#include "ParameterLabel.h"
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace ScatteringLength {

double CoulombX(const PPair *pair) {
  const double z1z2 = (double)pair->GetZ(1) * pair->GetZ(2);
  if (z1z2 <= 0.0) return 0.0;
  // x^2 = 8 Z1 Z2 e^2 mu a / hbar^2 = 8 Z1 Z2 alpha (mu c^2) a / (hbar c)
  return sqrt(8.0 * z1z2 * fstruc * pair->GetRedMass() * uconv * pair->GetChRad() / hbarc);
}

double ZeroEnergyShift(const PPair *pair, int l) {
  const double x = CoulombX(pair);
  if (x <= 0.0) return -(double)l;
  // The zero-energy irregular solution is sqrt(r) K_{2l+1}(x), x = sqrt(8 r / a_B); r d/dr = (x/2) d/dx
  // and K_nu' = -K_{nu-1} - (nu/x) K_nu give S_l(0) = -l - (x/2) K_{2l}/K_{2l+1}.  Scaled Bessel
  // functions: the ratio is unchanged and nothing overflows.
  const double kNum = gsl_sf_bessel_Kn_scaled(2 * l, x);
  const double kDen = gsl_sf_bessel_Kn_scaled(2 * l + 1, x);
  return -(double)l - 0.5 * x * kNum / kDen;
}

Coefficients FormulaCoefficients(const PPair *pair) {
  const double a = pair->GetChRad();
  const double x = CoulombX(pair);
  if (x <= 0.0) return {-a, a};  // a0 = a (1 - M)
  // a0 = -a [ M / (x^2 K1^2) - 2 I1 / (x^2 K1) ] with K1 = e^-x K1s, I1 = e^x I1s.
  const double k1s = gsl_sf_bessel_K1_scaled(x);
  const double i1s = gsl_sf_bessel_I1_scaled(x);
  const double e2x = exp(2.0 * x);
  return {-a * e2x / (x * x * k1s * k1s), 2.0 * a * e2x * i1s / (x * x * k1s)};
}

std::vector<Channel> Channels(CNuc *compound, int pairKey) {
  std::vector<Channel> out;
  if (!compound->IsPairKey(pairKey)) return out;
  const int pairNum = compound->GetPairNumFromKey(pairKey);
  PPair *pair = compound->GetPair(pairNum);
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
    JGroup probe(s, parity);
    probe.SetNucleus(pair->GetNucleus(), pair->GetNucleusZ(), pair->GetNucleusA());
    const int j = compound->IsJGroup(probe);
    if (j > 0) {
      JGroup *jg = compound->GetJGroup(j);
      const int ch = jg->IsChannel(AChannel(0, s, pairNum, 'P'));
      if (ch > 0) {
        c.jGroup = j;
        c.channel = ch;
        c.hasLevels = jg->IsInRMatrix() && jg->NumLevels() > 0;
      }
    }
    out.push_back(c);
  }
  return out;
}

namespace {
/// Fill and invert the level matrices at the pair's threshold.
void ThresholdMatrices(CNuc *compound, const Config &configure, EPoint &point, AMatrixFunc &f) {
  point.CalcEDependentValues(compound, configure, false);
  f.ClearMatrices();
  f.FillMatrices(&point);
  f.InvertMatrices();
}
}  // namespace

Value Evaluate(CNuc *compound, const Config &configure, const Channel &ch) {
  Value v;
  v.ch = ch;
  v.M = complex(0.0, 0.0);
  if (ch.hasLevels) {
    EPoint point(0.0, 0.0, ch.pairKey, ch.pairKey, false, false, false, 0.0, 0, 0);
    AMatrixFunc f(compound, configure);
    ThresholdMatrices(compound, configure, point, f);
    v.M = f.GetUBilinear(ch.jGroup, ch.channel, ch.channel);
  }
  const Coefficients c = FormulaCoefficients(compound->GetPair(ch.pairNum));
  v.a0 = c.c1 * v.M + c.c0;
  return v;
}

bool Adjoint(CNuc *compound, const Config &configure, const Channel &ch, complex aBar, GradAccum &accum,
             const vector_matrix_r *shiftDeriv) {
  if (!ch.hasLevels) return true;  // a0 does not depend on any parameter
  EPoint point(0.0, 0.0, ch.pairKey, ch.pairKey, false, false, false, 0.0, 0, 0);
  AMatrixFunc f(compound, configure);
  ThresholdMatrices(compound, configure, point, f);
  const Coefficients c = FormulaCoefficients(compound->GetPair(ch.pairNum));
  return f.BilinearAdjoint(&point, ch.jGroup, ch.channel, c.c1 * aBar, accum, shiftDeriv);
}

int ReadConstraints(CNuc *compound, const Config &configure, std::vector<Constraint> &constraints) {
  constraints.clear();
  std::ifstream in(configure.configfile.c_str());
  if (!in) return 0;
  std::string line;
  while (!in.eof() && line != "<scatteringLengths>") {
    getline(in, line);
    size_t a = line.find_first_not_of(" \t\r"), b = line.find_last_not_of(" \t\r");
    line = (a == std::string::npos) ? std::string() : line.substr(a, b - a + 1);
  }
  if (line != "<scatteringLengths>") return 0;
  int n = 0;
  auto error = [&](const std::string &what) {
    configure.outStream << "ERROR: <scatteringLengths> line " << n << ": " << what << "." << std::endl;
    return -1;
  };
  while (getline(in, line)) {
    if (line.find("</scatteringLengths>") != std::string::npos) break;
    n++;
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    std::istringstream stm(line);
    Constraint c;
    int active;
    if (!(stm >> active)) continue;  // blank or comment
    if (!(stm >> c.pairKey >> c.s >> c.value >> c.sigma))
      return error("expected 'active pairKey s value sigma [imValue imSigma]'");
    double im, imSigma;
    if (stm >> im) {
      if (!(stm >> imSigma)) return error("an imaginary part needs its uncertainty too");
      c.imValue = im;
      c.imSigma = imSigma;
    }
    std::string extra;
    if (stm >> extra) return error("unexpected '" + extra + "'");
    c.active = (active != 0);
    c.line = n;
    if (c.sigma < 0.0 || c.imSigma < 0.0) return error("uncertainties must not be negative");
    if (!compound->IsPairKey(c.pairKey)) return error("no particle pair with key " + std::to_string(c.pairKey));
    PPair *pair = compound->GetPair(compound->GetPairNumFromKey(c.pairKey));
    if (pair->GetPType() != 0) return error("pair " + std::to_string(c.pairKey) + " is not a particle pair");
    std::vector<Channel> all = Channels(compound, c.pairKey);
    if (c.s < 0.0) {
      if (pair->IsIdentical()) return error("the coherent combination (s = -1) is not defined for an identical pair");
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
        return error(m.str());
      }
    }
    constraints.push_back(c);
  }
  return (int)constraints.size();
}

complex ConstraintValue(CNuc *compound, const Config &configure, const Constraint &c) {
  complex a(0.0, 0.0);
  for (size_t t = 0; t < c.terms.size(); t++) a += c.weights[t] * Evaluate(compound, configure, c.terms[t]).a0;
  return a;
}

int NumResiduals(const std::vector<Constraint> &constraints) {
  int n = 0;
  for (const Constraint &c : constraints) {
    if (!c.active) continue;
    if (c.sigma > 0.0) n++;
    if (c.imSigma > 0.0) n++;
  }
  return n;
}

void Residuals(CNuc *compound, const Config &configure, const std::vector<Constraint> &constraints,
               std::vector<double> &r) {
  for (const Constraint &c : constraints) {
    if (!c.active || (c.sigma <= 0.0 && c.imSigma <= 0.0)) continue;
    complex a = ConstraintValue(compound, configure, c);
    if (c.sigma > 0.0) r.push_back((std::real(a) - c.value) / c.sigma);
    if (c.imSigma > 0.0) r.push_back((std::imag(a) - c.imValue) / c.imSigma);
  }
}

double Chi2(CNuc *compound, const Config &configure, const std::vector<Constraint> &constraints) {
  std::vector<double> r;
  Residuals(compound, configure, constraints, r);
  double chi2 = 0.0;
  for (double x : r) chi2 += x * x;
  return chi2;
}

bool Chi2Gradient(CNuc *compound, const Config &configure, const std::vector<Constraint> &constraints,
                  GradAccum &accum, const vector_matrix_r *shiftDeriv, int residual) {
  int row = 0;
  for (const Constraint &c : constraints) {
    if (!c.active || (c.sigma <= 0.0 && c.imSigma <= 0.0)) continue;
    // d(chi2) = Re(conj(aBar) da):  aBar = 2 r_re/sigma + i 2 r_im/imSigma.  One residual row:
    // aBar = 1/sigma (Re row) or i/imSigma (Im row), since Re(conj(i/s) da) = Im(da)/s.
    complex aBar(0.0, 0.0);
    bool needed = false;
    if (residual < 0) {
      complex a = ConstraintValue(compound, configure, c);
      if (c.sigma > 0.0) aBar += 2.0 * (std::real(a) - c.value) / (c.sigma * c.sigma);
      if (c.imSigma > 0.0) aBar += complex(0.0, 2.0 * (std::imag(a) - c.imValue) / (c.imSigma * c.imSigma));
      needed = true;
    } else {
      if (c.sigma > 0.0) {
        if (row == residual) {
          aBar = 1.0 / c.sigma;
          needed = true;
        }
        row++;
      }
      if (c.imSigma > 0.0) {
        if (row == residual) {
          aBar = complex(0.0, 1.0 / c.imSigma);
          needed = true;
        }
        row++;
      }
    }
    if (!needed) continue;
    for (size_t t = 0; t < c.terms.size(); t++)
      if (!Adjoint(compound, configure, c.terms[t], c.weights[t] * aBar, accum, shiftDeriv)) return false;
    if (residual >= 0) break;
  }
  return true;
}

std::string Report(CNuc *compound, const Config &configure, const std::vector<Constraint> &constraints) {
  std::ostringstream out;
  out << "s-wave scattering lengths at each particle pair's threshold" << std::endl
      << "(a0 = -a [M/(x^2 K1^2) - 2 I1/(x^2 K1)], M = gamma^T A gamma at E = 0; Paneru et al., "
         "PRC 99, 045807 (2019), Eq. 7)" << std::endl;
  out << std::setw(6) << "Pair" << "  " << std::left << std::setw(14) << "Particles" << std::setw(6) << "J^pi"
      << std::right << std::setw(10) << "x" << std::setw(16) << "Re a0 (fm)" << std::setw(16) << "Im a0 (fm)"
      << "  note" << std::endl;
  for (int aa = 1; aa <= compound->NumPairs(); aa++) {
    PPair *pair = compound->GetPair(aa);
    if (pair->GetPType() != 0) continue;
    for (const Channel &ch : Channels(compound, pair->GetPairKey())) {
      Value v = Evaluate(compound, configure, ch);
      out << std::setw(6) << ch.pairKey << "  " << std::left << std::setw(14) << AZURELabel::Pair(compound, aa)
          << std::setw(6) << AZURELabel::JPi(ch.s, ch.parity) << std::right << std::fixed << std::setprecision(4)
          << std::setw(10) << CoulombX(pair) << std::scientific << std::setprecision(6) << std::setw(16)
          << std::real(v.a0) << std::setw(16) << std::imag(v.a0);
      if (!ch.hasLevels) out << "  no level in this channel (hard sphere)";
      out << std::defaultfloat << std::endl;
    }
  }
  bool any = false;
  double total = 0.0;
  for (const Constraint &c : constraints) {
    if (!any) {
      out << "Constraints (<scatteringLengths>):" << std::endl;
      any = true;
    }
    complex a = ConstraintValue(compound, configure, c);
    double chi2 = 0.0;
    if (c.active && c.sigma > 0.0) chi2 += pow((std::real(a) - c.value) / c.sigma, 2);
    if (c.active && c.imSigma > 0.0) chi2 += pow((std::imag(a) - c.imValue) / c.imSigma, 2);
    total += chi2;
    out << "  line " << c.line << ": pair " << c.pairKey << ", "
        << (c.s < 0.0 ? std::string("coherent (spin-weighted)") : "s = " + AZURELabel::Spin(c.s))
        << (c.active ? "" : " (inactive)") << ": model " << std::real(a);
    if (std::imag(a) != 0.0) out << (std::imag(a) < 0 ? " - " : " + ") << fabs(std::imag(a)) << " i";
    out << " fm, measured " << c.value << " +- " << c.sigma;
    if (c.imSigma > 0.0) out << ", Im " << c.imValue << " +- " << c.imSigma;
    out << " fm, chi2 " << chi2 << std::endl;
  }
  if (any) out << "  total scattering-length chi2 " << total << std::endl;
  return out.str();
}

}  // namespace ScatteringLength
