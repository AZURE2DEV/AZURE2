#include "ThmLineshape.h"

#include "ALevel.h"
#include "AChannel.h"
#include "CNuc.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "JGroup.h"
#include "NuclearPotentialManager.h"
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
  // levels.  The key is everything the width depends on, by value: the memo
  // outlives a CNuc (thread_local), and a later session can get the same
  // level address back (the S_c memo of THMMatrixFunc was keyed by address
  // once and leaked across sessions, 49b6c30).
  struct Entry {
    std::vector<double> key;
    double width;
  };
  static const int kMemo = 64;
  thread_local std::vector<Entry> memo;
  thread_local int next = 0;
  thread_local std::vector<double> key;  // reused: no allocation per call
  key.clear();
  key.push_back(useGSL);
  key.push_back(rmc);
  key.push_back(jgroup->GetJ());
  key.push_back(jgroup->GetPi());
  key.push_back(energy);
  const bool hybridOn = g_config ? g_config->useHybridMethod : false;
  for (int ch = 1; ch <= numChannels; ch++) {
    AChannel *channel = jgroup->GetChannel(ch);
    PPair *pair = compound->GetPair(channel->GetPairNum());
    const bool hybrid = hybridOn && NuclearPotentialManager::instance().isPairEnabled(pair->GetPairKey());
    const double values[] = {level->GetFitGamma(ch), (double)channel->GetL(), (double)channel->GetRadType(),
                             (double)pair->GetZ(1), (double)pair->GetZ(2), pair->GetRedMass(), pair->GetChRad(),
                             pair->GetSepE(), pair->GetExE(), pair->GetJ(2), (double)pair->GetPi(2),
                             hybrid ? (double)NuclearPotentialManager::instance().tagFor(pair->GetPairKey()) : 0.0};
    key.insert(key.end(), values, values + sizeof(values) / sizeof(values[0]));
  }
  for (const Entry &m : memo)
    if (m.key == key) return m.width;

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
    memo.push_back(Entry{key, width});
  } else {
    memo[next] = Entry{key, width};
    next = (next + 1) % kMemo;
  }
  return width;
}

bool ThmSpectatorWindow::NodesAt(double energy, std::vector<Node> &out) const {
  out.clear();
  if (!acc) return false;
  const double esf = acc->EsF(energy);
  if (!(esf > 0.0)) return false;
  const double ksf = std::sqrt(2.0 * acc->sf.mu * esf) / hbarc;
  std::vector<ThmDistortion::AngleNode> an;
  if (!acc->AngleNodes(ksf, an)) return false;
  double total = 0.0;
  for (const ThmDistortion::AngleNode &n : an) {
    if (!(n.w > 0.0)) continue;
    Node k;
    k.theta = n.theta;
    k.p = n.q * hbarc;
    k.weight = n.w * phi2(k.p);
    k.es = k.p * k.p / (2.0 * muSx);
    if (!(k.weight > 0.0)) continue;
    total += k.weight;
    out.push_back(k);
  }
  if (!(total > 0.0) || !std::isfinite(total)) {
    out.clear();
    return false;
  }
  for (Node &k : out) k.weight /= total;
  return true;
}

void ThmSpectatorWindow::Nodes(double energy, std::vector<Node> &out, bool *moved) const {
  if (moved) *moved = false;
  if (NodesAt(energy, out) || dataE.empty()) return;
  // The nearest data energy (every one has accepted directions).
  size_t i = std::lower_bound(dataE.begin(), dataE.end(), energy) - dataE.begin();
  if (i == dataE.size() || (i > 0 && energy - dataE[i - 1] < dataE[i] - energy)) i = i == 0 ? 0 : i - 1;
  if (moved) *moved = true;
  NodesAt(dataE[i], out);
}

double ThmSpectatorWindow::MeanEs(double energy) const {
  std::vector<Node> nodes;
  if (!NodesAt(energy, nodes)) return 0.0;
  double m = 0.0;
  for (const Node &k : nodes) m += k.weight * k.es;
  return m;
}

bool ThmSpectatorWindow::Reach(double energy, double &qLo, double &qHi) const {
  if (!acc) return false;
  const double esf = acc->EsF(energy);
  if (!(esf > 0.0)) return false;
  const double ks = std::sqrt(2.0 * acc->sf.mu * esf) / hbarc, kb = acc->beta * acc->aa.k;
  qLo = std::fabs(ks - kb) * hbarc;
  qHi = (ks + kb) * hbarc;
  return true;
}

double ThmSpectatorWindow::Density(double energy, double q) const {
  double qLo, qHi;
  if (!Reach(energy, qLo, qHi) || q < qLo || q > qHi || q < pMin || q > pMax) return 0.0;
  const ThmDistortion &d = *acc;
  const double ks = std::sqrt(2.0 * d.sf.mu * d.EsF(energy)) / hbarc, kb = d.beta * d.aa.k, qq = q / hbarc;
  const double x = std::max(-1.0, std::min(1.0, (ks * ks + kb * kb - qq * qq) / (2.0 * ks * kb)));
  const double theta = std::acos(d.kin.horseIsBeam ? x : -x);  // c.m. angle to the beam
  double a = 1.0;
  if (!d.angAll) {
    const double vs = hbarc * ks / (d.kin.ms * uconv), g = vs > 0.0 ? d.vcm / vs : HUGE_VAL;
    const double t = (d.angCm ? theta : std::atan2(std::sin(theta), std::cos(theta) + g)) * 180.0 / M_PI;
    if (t < d.angLo || t > d.angHi) return 0.0;
    if (!d.angT.empty()) {
      size_t j = std::upper_bound(d.angT.begin(), d.angT.end(), t) - d.angT.begin();
      j = std::max<size_t>(1, std::min(j, d.angT.size() - 1));
      const double f = (t - d.angT[j - 1]) / (d.angT[j] - d.angT[j - 1]);
      a = d.angW[j - 1] + f * (d.angW[j] - d.angW[j - 1]);
    }
  }
  return a * phi2(q) * q;
}

std::function<double(double)> ThmPsDistribution(const ThmExperiment &x) {
  switch (x.psKind) {
    case ThmExperiment::PS_HULTHEN: {
      const double a2 = x.psA * x.psA, b2 = x.psB * x.psB;
      return [a2, b2](double p) {
        double q2 = (p / hbarc) * (p / hbarc);
        double phi = 1.0 / (a2 + q2) - 1.0 / (b2 + q2);
        return phi * phi;
      };
    }
    case ThmExperiment::PS_GAUSS: {
      const double c = 4.0 * std::log(2.0) / (x.psFwhm * x.psFwhm);
      return [c](double p) { return std::exp(-c * p * p); };
    }
    case ThmExperiment::PS_TABLE: {
      const std::vector<double> tp = x.psTableP, tw = x.psTableW;
      return [tp, tw](double p) {
        if (p <= tp.front()) return tw.front();
        if (p >= tp.back()) return tw.back();
        size_t hi = std::upper_bound(tp.begin(), tp.end(), p) - tp.begin(), lo = hi - 1;
        return tw[lo] + (tw[hi] - tw[lo]) * (p - tp[lo]) / (tp[hi] - tp[lo]);
      };
    }
    default:
      return nullptr;
  }
}

std::string BuildThmSpectatorWindow(const ThmExperiment &x, double muSx, const ThmDistortion::Kinematics &k,
                                    const std::vector<double> &dataE, ThmSpectatorWindow &out) {
  out = ThmSpectatorWindow();
  out.experiment = x.name;
  out.muSx = muSx;
  out.pMin = x.psMin;
  out.pMax = x.psMax;
  out.phi2 = ThmPsDistribution(x);
  if (!out.phi2) return "no ps window";
  if (!(muSx > 0.0)) return "the s + x reduced mass is not positive";
  std::ostringstream d;
  d.precision(8);
  switch (x.psKind) {
    case ThmExperiment::PS_HULTHEN:
      d << "hulthen a=" << x.psA << " b=" << x.psB << " fm^-1";
      break;
    case ThmExperiment::PS_GAUSS:
      d << "gauss FWHM=" << x.psFwhm << " MeV/c";
      break;
    default:
      d << "table " << x.psTable;
      break;
  }
  d << ", p_s in [" << x.psMin << ", " << x.psMax << "] MeV/c";
  // The kinematics and the accepted directions as R(E) and the DW vertex take
  // them (ThmDistortion::Setup), without waves: plane s + F and a + A
  // channels, so that only E_sF > 0 is asked of an energy.
  ThmExperiment xc = x;
  xc.distortion = ThmExperiment::DIST_OPTICAL;
  xc.opticalAA = ThmExperiment::Optical();
  xc.opticalSF = ThmExperiment::Optical();
  xc.opticalAA.kind = xc.opticalSF.kind = 0;
  std::shared_ptr<ThmDistortion> acc = std::make_shared<ThmDistortion>();
  out.dataE = dataE;
  std::sort(out.dataE.begin(), out.dataE.end());
  const double eLo = out.dataE.empty() ? 0.0 : out.dataE.front(), eHi = out.dataE.empty() ? 0.0 : out.dataE.back();
  std::string why = acc->Setup(xc, k, eLo);
  if (!why.empty()) return why;
  // The branches of a lab window up to 0.5 MeV above the data, short of the
  // spectator's threshold (EData::BuildThmGroups for R(E)).
  acc->SetAngleSlots(std::min(eHi + 0.5, acc->eAA - k.bind - 0.5 * acc->EsF(eHi)));
  out.acc = acc;
  for (double e : out.dataE) {
    std::vector<ThmSpectatorWindow::Node> nodes;
    if (out.NodesAt(e, nodes)) continue;
    why = acc->CheckWindow(e);
    return why.empty() ? "at E = " + std::to_string(e) + " MeV the window has no weight (|phi|^2 = 0 there)" : why;
  }
  if (x.angleWindow)
    d << "; spectator directions " << acc->AngleText() << ", " << acc->angNodes << " nodes in cos theta_cm"
      << (acc->angSlots > 1 ? " per branch" : "");
  else if (acc->angNodes == 1)
    d << ", one node";
  else
    d << ", " << acc->angNodes << " Gauss-Legendre nodes in cos theta_cm";
  d << "; weight |phi(p_s)|^2 d cos theta_cm" << (acc->angT.empty() ? "" : " x acceptance");
  out.description = d.str();
  return "";
}
