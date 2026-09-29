#include "ThmLineshape.h"

#include "ALevel.h"
#include "AChannel.h"
#include "CNuc.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
#include "ThmExperiment.h"

#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <vector>

namespace {
// E_sF below this (MeV) is clamped: the spectator is at the kinematic limit.
const double kMinEsF = 1.0e-3;
}  // namespace

double ThmLineshape::KsF(double energy) const {
  return std::sqrt(2.0 * MuSF() * std::max(EsF(energy), kMinEsF)) / hbarc;
}

double ThmLineshape::Eta0(double energy) const {
  return Zs * ZF * fstruc * MuSF() / (hbarc * KsF(energy));
}

double ThmLineshape::EtaSB(double energy, int ZB, double mB) const {
  double muSB = ms * mB / (ms + mB) * uconv;
  return Zs * ZB * fstruc * muSB / (hbarc * KsF(energy));
}

double ThmLineshape::EtaSbEstimate(double energy, double eBB, int Zb, double mb, double mB) const {
  if (Zs == 0 || Zb == 0) return 0.0;
  double vs = std::sqrt(2.0 * std::max(EsF(energy), kMinEsF) / MuSF());  // units of c
  double muBB = mb * mB / (mb + mB) * uconv;
  double vb = eBB > 0.0 ? std::sqrt(2.0 * eBB / muBB) * mB / (mb + mB) : 0.0;
  return Zs * Zb * fstruc / std::max(vs, vb);
}

complex ThmLineshapeFactor(double zeta, double x, double Gamma) {
  if (zeta == 0.0) return complex(1.0, 0.0);
  double modulus = std::sqrt(ThmLineshapeFactorSq(zeta, x, Gamma));
  double r = std::hypot(x, 0.5 * std::max(Gamma, 0.0));
  double phase = r > 0.0 ? -zeta * std::log(r) : 0.0;
  return std::polar(modulus, phase);
}

double ThmLineshapeFactorSq(double zeta, double x, double Gamma) {
  if (zeta == 0.0) return 1.0;
  double angle;
  if (Gamma > 0.0)
    angle = std::atan(2.0 * x / Gamma);
  else
    angle = x > 0.0 ? 0.5 * M_PI : x < 0.0 ? -0.5 * M_PI : 0.0;
  return std::exp(2.0 * zeta * angle);
}

double ThmLevelWidth(CNuc *compound, JGroup *jgroup, ALevel *level, const Config &configure) {
  const int numChannels = jgroup->NumChannels();
  const bool useGSL = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const bool rmc = !!(configure.paramMask & Config::USE_RMC_FORMALISM);
  const double energy = level->GetFitE();

  // Per-thread memo: within one evaluation every point asks for the same few
  // levels.  The key is everything the width depends on.
  struct Entry {
    const ALevel *level;
    double e;
    std::vector<double> key;  // gammas, then the channel radii
    double width;
  };
  static const int kMemo = 64;
  thread_local std::vector<Entry> memo;
  thread_local int next = 0;
  std::vector<double> key;
  key.reserve(2 * numChannels);
  for (int ch = 1; ch <= numChannels; ch++) key.push_back(level->GetFitGamma(ch));
  for (int ch = 1; ch <= numChannels; ch++)
    key.push_back(compound->GetPair(jgroup->GetChannel(ch)->GetPairNum())->GetChRad());
  for (const Entry &m : memo)
    if (m.level == level && m.e == energy && m.key == key) return m.width;

  double particle = 0.0, radiative = 0.0, normSum = 0.0;
  for (int ch = 1; ch <= numChannels; ch++) {
    AChannel *channel = jgroup->GetChannel(ch);
    PPair *pair = compound->GetPair(channel->GetPairNum());
    double gamma = level->GetFitGamma(ch);
    if (gamma == 0.0) continue;
    double localEnergy = energy - pair->GetSepE() - pair->GetExE();
    char type = channel->GetRadType();
    if (type == 'P') {
      ChannelFunc channelFunc(pair, useGSL);
      normSum += channelFunc.ShiftDerivative(channel->GetL(), localEnergy) * gamma * gamma;
      if (!channelFunc.IsClosed(localEnergy))
        particle += 2.0 * gamma * gamma * channelFunc.Penetrability(channel->GetL(), localEnergy);
    } else if (type == 'M' || type == 'E') {
      // A ground-state transition parametrizes a moment, not a width.
      if (std::fabs(energy - pair->GetExE()) < 1.0e-3 && jgroup->GetJ() == pair->GetJ(2) &&
          jgroup->GetPi() == pair->GetPi(2))
        continue;
      double pene = rmc ? 1.0 : std::pow(std::fabs(localEnergy) / hbarc, 2.0 * channel->GetL() + 1.0);
      radiative += 2.0 * gamma * gamma * pene;
    }
  }
  double width = (particle + radiative) / (1.0 + normSum);
  if (!(1.0 + normSum > 0.0) || !std::isfinite(width)) width = particle + radiative;

  if ((int)memo.size() < kMemo) {
    memo.push_back(Entry{level, energy, key, width});
  } else {
    memo[next] = Entry{level, energy, key, width};
    next = (next + 1) % kMemo;
  }
  return width;
}

double ThmSpectatorWindow::MeanEs() const {
  double m = 0.0;
  for (size_t k = 0; k < es.size(); k++) m += weight[k] * es[k];
  return m;
}

std::string BuildThmSpectatorWindow(const ThmExperiment &x, double muSx, ThmSpectatorWindow &out) {
  out = ThmSpectatorWindow();
  out.experiment = x.name;
  out.muSx = muSx;
  out.pMin = x.psMin;
  out.pMax = x.psMax;
  std::ostringstream d;
  d.precision(8);
  // Event weight per unit p_s.
  std::function<double(double)> w;
  switch (x.psKind) {
    case ThmExperiment::PS_HULTHEN: {
      const double a2 = x.psA * x.psA, b2 = x.psB * x.psB;
      w = [a2, b2](double p) {
        double q2 = (p / hbarc) * (p / hbarc);
        double phi = 1.0 / (a2 + q2) - 1.0 / (b2 + q2);
        return phi * phi * p * p;
      };
      d << "hulthen a=" << x.psA << " b=" << x.psB << " fm^-1";
      break;
    }
    case ThmExperiment::PS_GAUSS: {
      const double c = 4.0 * std::log(2.0) / (x.psFwhm * x.psFwhm);
      w = [c](double p) { return std::exp(-c * p * p) * p * p; };
      d << "gauss FWHM=" << x.psFwhm << " MeV/c";
      break;
    }
    case ThmExperiment::PS_TABLE: {
      const std::vector<double> &tp = x.psTableP, &tw = x.psTableW;
      w = [&tp, &tw](double p) {
        if (p <= tp.front()) return tw.front();
        if (p >= tp.back()) return tw.back();
        size_t hi = std::upper_bound(tp.begin(), tp.end(), p) - tp.begin(), lo = hi - 1;
        return tw[lo] + (tw[hi] - tw[lo]) * (p - tp[lo]) / (tp[hi] - tp[lo]);
      };
      d << "table " << x.psTable;
      break;
    }
    default:
      return "no ps window";
  }
  d << ", p_s in [" << x.psMin << ", " << x.psMax << "] MeV/c";
  if (!(muSx > 0.0)) return "the s + x reduced mass is not positive";
  if (x.psMax == x.psMin) {
    out.p.push_back(x.psMin);
    out.weight.push_back(1.0);
    d << ", one node";
  } else {
    const int n = x.psNodes;
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
    double total = 0.0;
    for (int i = 0; i < n; i++) {
      double xi, wi;
      gsl_integration_glfixed_point(x.psMin, x.psMax, i, &xi, &wi, t);
      double v = wi * w(xi);
      out.p.push_back(xi);
      out.weight.push_back(v);
      total += v;
    }
    gsl_integration_glfixed_table_free(t);
    if (!(total > 0.0) || !std::isfinite(total))
      return "the weight vanishes at every Gauss-Legendre node of the window (widen it or raise psNodes)";
    for (double &v : out.weight) v /= total;
    d << ", " << n << " Gauss-Legendre nodes";
  }
  for (double p : out.p) out.es.push_back(p * p / (2.0 * muSx));
  out.description = d.str();
  return "";
}
