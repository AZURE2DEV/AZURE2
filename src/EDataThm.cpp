/*
 * The THM experiments of an EData (<thm> experiment[<name>] ...; ThmExperiment.h):
 * building the groups and their physics objects, the profiled scale and
 * background, the coherent-background parameters and the reports and
 * tables (thm_experiments.out, pyazr thm_*).  Members of EData, kept out of
 * EData.cpp so that the classic data handling stays as it is on dev.
 */
#include "CNuc.h"
#include "PPair.h"
#include "ThmFunc.h"
#include "ThmVertexBoundary.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "EData.h"
#include "Minuit2/MnUserParameters.h"
#include "GSLException.h"
#include "NuclearPotentialManager.h"
#include "ThmLineshape.h"
#include "ThmReports.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmAngular.h"
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

void EData::FillThmCoherentFromParams(const vector_r &p, CNuc *theCNuc) {
  if (thmCoherentParams_.empty() || thmCoherentParamOffset_ < 0) return;
  const size_t n = thmCoherentParams_.size();
  if (p.size() < (size_t)thmCoherentParamOffset_ + n) return;
  std::vector<double> v(n);
  for (size_t k = 0; k < n; k++) v[k] = thmCoherentParams_[k].value = p[thmCoherentParamOffset_ + k];
  if (theCNuc) theCNuc->SetThmCoherentValues(v);
}

// ---------------------------------------------------------------------------
// THM experiments (<thm> experiment[<name>] ...): see ThmExperiment.h.

namespace {

const double kAmu = 931.49410242;  // MeV/u (CODATA 2018)

// The lowest and highest c.m. energy of the points of a group's segments.
void ThmGroupEnergyRange(EData &data, const EData::ThmGroup &group, double &eLo, double &eHi) {
  eLo = 1.0e300;
  eHi = -1.0e300;
  for (int s : group.segments)
    for (int p = 1; p <= data.GetSegment(s)->NumPoints(); p++) {
      eLo = std::min(eLo, data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
      eHi = std::max(eHi, data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
    }
}

// What the kinematics of an experiment line give (ThmResolveKinematics).
struct ThmReaction {
  ThmDistortion::Kinematics dk;  ///< for distortion=coulomb|optical and the DW vertex
  PPair *pair = nullptr;         ///< the entrance pair x + A
  int pairKey = 0;
  double mX = 0.0;               ///< u, the pair nucleus that is x
  double bind = 0.0;             ///< B(x+s) from the masses (MeV)
};

/*
 * The three-body reaction of experiment x (beam, target, spectator, Ebeam):
 * the entrance pair x + A of its segments, which of beam and target is the
 * Trojan horse a = x + s, the kinematics, and the startup lines that say so.
 * False (and the ERROR written) if they do not describe the segments.
 */
bool ThmResolveKinematics(EData &data, const Config &configure, CNuc *theCNuc, const ThmExperiment &x,
                          const EData::ThmGroup &group, const std::string &where, ThmReaction &r) {
  ThmDistortion::Kinematics &dk = r.dk;
  int pairKey = data.GetSegment(group.segments[0])->GetEntranceKey();
  for (int s : group.segments)
    if (data.GetSegment(s)->GetEntranceKey() != pairKey) {
      configure.outStream << where << "beam/target/spectator describe one reaction, but its segments "
                             "have different entrance pairs."
                          << std::endl;
      return false;
    }
  PPair *pair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(pairKey));
  const int Z[2] = {pair->GetZ(1), pair->GetZ(2)};
  const double M[2] = {pair->GetM(1), pair->GetM(2)};
  const ThmNuclide &b = x.beam, &t = x.target, &sp = x.spectator;
  ThmReactionKinematics rk;
  const std::string why = ThmResolveReaction(b, t, sp, x.beamEnergy, Z, M, rk);
  if (!why.empty()) {
    configure.outStream << where << why << std::endl;
    return false;
  }
  const ThmNuclide &th = rk.horse, &nA = rk.nucleusA;
  const double mX = rk.mX, bind = rk.bind, exa = rk.exa;
  dk.Za = th.Z;
  dk.ZA = nA.Z;
  dk.Zs = sp.Z;
  dk.Zx = th.Z - sp.Z;
  dk.Aa = th.A;
  dk.AA = nA.A;
  dk.As = sp.A;
  dk.Ax = th.A - sp.A;
  dk.ma = th.mass;
  dk.mA = nA.mass;
  dk.ms = sp.mass;
  dk.mx = mX;
  dk.horseIsBeam = rk.horseIsBeam;
  dk.mBeam = b.mass;
  dk.mTarget = t.mass;
  dk.beamEnergy = x.beamEnergy;
  dk.bind = bind;
  std::ostringstream k;
  k.precision(6);
  k << "  " << b.name << " + " << t.name << " at " << x.beamEnergy << " MeV (lab), Trojan horse "
    << th.name << " = x + " << sp.name << ", B(x+s) = " << bind << " MeV; quasi-free E(x+A) = " << exa
    << " MeV, E_qf = E(x+A) - B = " << exa - bind << " MeV.";
  configure.outStream << k.str() << std::endl;
  // The vertex takes B from the entrance pair's channel lines (field 32),
  // the kinematics from the masses; say so when they disagree.
  if (std::fabs(pair->GetBindingEnergy() - bind) > 1.0e-3) {
    std::ostringstream w;
    w.precision(6);
    w << "WARNING: <thm> experiment[" << x.name << "]: B(x+s) from the masses of " << th.name << " = x + "
      << sp.name << " is " << bind << " MeV, but the entrance pair " << pairKey
      << " carries B = " << pair->GetBindingEnergy()
      << " MeV (field 32 of its channel lines); the THM vertex uses field 32, the kinematics of this "
         "experiment (the quasi-free energy above, E_sF of the line shape and of the distortion "
         "factor) use the masses.";
    configure.outStream << w.str() << std::endl;
  }
  r.pair = pair;
  r.pairKey = pairKey;
  r.mX = mX;
  r.bind = bind;
  return true;
}

// Spectator-momentum window (ThmLineshape.h ThmSpectatorWindow) of a ps= line.
bool ThmBuildWindow(EData &data, const Config &configure, const ThmExperiment &x, EData::ThmGroup &group,
                    const std::string &where, const ThmReaction &r) {
  if (configure.thm.SpectatorEnergy(r.pairKey) != 0.0) {
    configure.outStream << where << "a ps window and spectatorEnergy both set the spectator motion of "
                           "entrance pair "
                        << r.pairKey << "; use one (ps=delta keeps spectatorEnergy)." << std::endl;
    return false;
  }
  // With vertexModel=dw the DW vertex averages over the accepted
  // directions itself (ThmDwVertex); the plane-wave nodes are not used.
  if (x.vertexDW) return true;
  const ThmNuclide &sp = x.spectator;
  double muSx = r.mX * sp.mass / (r.mX + sp.mass) * kAmu;
  std::vector<double> dataE;
  for (int s : group.segments)
    for (int p = 1; p <= data.GetSegment(s)->NumPoints(); p++)
      dataE.push_back(data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
  std::shared_ptr<ThmSpectatorWindow> window = std::make_shared<ThmSpectatorWindow>();
  std::string why = BuildThmSpectatorWindow(x, muSx, r.dk, dataE, *window);
  if (!why.empty()) {
    configure.outStream << where << "ps: " << why << "." << std::endl;
    return false;
  }
  std::vector<ThmSpectatorWindow::Node> lo, hi;
  window->NodesAt(window->dataE.front(), lo);
  window->NodesAt(window->dataE.back(), hi);
  auto range = [](const std::vector<ThmSpectatorWindow::Node> &n) {
    double a = n.front().p, b = a;
    for (const ThmSpectatorWindow::Node &k : n) a = std::min(a, k.p), b = std::max(b, k.p);
    std::ostringstream t;
    t.precision(6);
    t << a << "-" << b;
    return t.str();
  };
  std::ostringstream w;
  w.precision(6);
  w << "  Spectator-momentum window: " << window->description << "; mu_sx = " << muSx
    << " MeV; at the lowest point (E = " << window->dataE.front() << " MeV) " << lo.size()
    << " node(s), |p_s| = " << range(lo) << " MeV/c, <T_s> = " << window->MeanEs(window->dataE.front())
    << " MeV; at the highest (E = " << window->dataE.back() << " MeV) " << hi.size()
    << " node(s), |p_s| = " << range(hi) << " MeV/c, <T_s> = " << window->MeanEs(window->dataE.back())
    << " MeV.";
  configure.outStream << w.str() << std::endl;
  group.window = window;
  for (int s : group.segments) data.GetSegment(s)->SetThmSpectatorWindow(window);
  return true;
}

// Coulomb line shape of the spectator (ThmLineshape.h), lineshape=on.  The
// level energies and widths it uses are the observed ones of Brune.
bool ThmBuildLineshape(EData &data, const Config &configure, CNuc *theCNuc, const ThmExperiment &x,
                       EData::ThmGroup &group, const std::string &where, const ThmReaction &r) {
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    configure.outStream << where << "lineshape=on uses the observed level energies and widths "
                           "as the resonance poles; it needs the Brune parameterization."
                        << std::endl;
    return false;
  }
  const ThmNuclide &b = x.beam, &t = x.target, &sp = x.spectator;
  PPair *pair = r.pair;
  std::shared_ptr<ThmLineshape> shape = std::make_shared<ThmLineshape>();
  shape->experiment = x.name;
  shape->spectator = sp.name;
  shape->Zs = sp.Z;
  shape->ms = sp.mass;
  shape->ZF = pair->GetZ(1) + pair->GetZ(2);
  shape->mF = pair->GetM(1) + pair->GetM(2);
  shape->eAA = x.beamEnergy * t.mass / (b.mass + t.mass);
  shape->bind = r.bind;
  // Every data point must leave the spectator some energy.
  double eMax = -1.0e300;
  for (int s : group.segments)
    for (int p = 1; p <= data.GetSegment(s)->NumPoints(); p++)
      eMax = std::max(eMax, data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
  if (!(shape->EsF(eMax) > 0.0)) {
    configure.outStream << where << "lineshape=on: at E = " << eMax
                        << " MeV the spectator has no energy left (E_sF = E_aA - B - E = "
                        << shape->eAA << " - " << r.bind << " - " << eMax << " MeV <= 0); check Ebeam."
                        << std::endl;
    return false;
  }
  std::ostringstream l;
  l.precision(6);
  l << "  Coulomb line shape on: E_sF = E_aA - B - E with E_aA = " << shape->eAA << " MeV, eta_0 = "
    << shape->Eta0(eMax) << " at the highest point energy (E = " << eMax << " MeV)"
    << (sp.Z == 0 ? "; the spectator is neutral, so N_C = 1." : ".");
  configure.outStream << l.str() << std::endl;
  for (int s : group.segments) {
    int key = data.GetSegment(s)->GetExitKey();
    bool seen = false;
    for (const ThmLineshape::Exit &e : shape->exits) seen = seen || e.pairKey == key;
    if (seen || !theCNuc->IsPairKey(key)) continue;
    PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(key));
    int light = exitPair->GetM(1) <= exitPair->GetM(2) ? 1 : 2;
    ThmLineshape::Exit e;
    e.pairKey = key;
    e.Zb = exitPair->GetZ(light);
    e.ZB = exitPair->GetZ(3 - light);
    e.mb = exitPair->GetM(light);
    e.mB = exitPair->GetM(3 - light);
    e.q = pair->GetSepE() + pair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
    shape->exits.push_back(e);
  }
  group.lineshape = shape;
  for (int s : group.segments) data.GetSegment(s)->SetThmLineshape(shape);
  return true;
}

// Distorted-wave entrance vertex (ThmDwVertex.h): replaces M_l in the HOES
// amplitude of every segment of the experiment; R(E) is not applied.
bool ThmBuildDwVertex(EData &data, const Config &configure, CNuc *theCNuc, const ThmExperiment &x,
                      EData::ThmGroup &group, const std::string &where, const ThmDistortion::Kinematics &dk) {
  int pairKey = data.GetSegment(group.segments[0])->GetEntranceKey();
  int pairNum = theCNuc->GetPairNumFromKey(pairKey);
  PPair *pair = theCNuc->GetPair(pairNum);
  // coulombIntegral=1 is refused with it by Config::ReadThmBlock
  // (CheckThmCoulombConsistency).
  const std::string refused =
      CheckThmDwVertexOptions(configure.thm.coherentL, configure.thm.SpectatorEnergy(pairKey), pairKey);
  if (!refused.empty()) {
    configure.outStream << where << refused << std::endl;
    return false;
  }
  std::vector<int> ls;
  for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
    JGroup *jg = theCNuc->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int ch = 1; ch <= jg->NumChannels(); ch++)
      if (jg->GetChannel(ch)->GetPairNum() == pairNum) ls.push_back(jg->GetChannel(ch)->GetL());
  }
  std::vector<double> energies;
  for (int s : group.segments)
    for (int p = 1; p <= data.GetSegment(s)->NumPoints(); p++)
      energies.push_back(data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
  double eLo, eHi;
  ThmGroupEnergyRange(data, group, eLo, eHi);
  std::shared_ptr<ThmDwVertex> v = std::make_shared<ThmDwVertex>();
  double eAA = dk.beamEnergy * dk.mTarget / (dk.mBeam + dk.mTarget);
  if (!(eAA - dk.bind - eHi > 0.0)) {
    configure.outStream << where << "vertexModel=dw: at E = " << eHi
                        << " MeV the spectator has no energy left (E_sF = E_aA - B - E = " << eAA << " - "
                        << dk.bind << " - " << eHi << " MeV <= 0); check Ebeam." << std::endl;
    return false;
  }
  double gridHi = std::min(eHi + 0.3, eAA - dk.bind - 0.5 * (eAA - dk.bind - eHi));
  std::string why = ls.empty() ? std::string("the entrance pair has no channel in the R matrix")
                               : v->Build(x, dk, pair->GetChRad(), ls, eLo - 0.3, gridHi, energies);
  if (!why.empty()) {
    configure.outStream << where << "vertexModel=dw: " << why << "." << std::endl;
    return false;
  }
  std::ostringstream l;
  l.precision(6);
  l << "  Entrance vertex: " << v->description << ".\n"
    << "  alpha = m_A/m_F = " << v->alpha << ", beta = m_s/m_a = " << v->beta << ", k_aA = " << v->dist.aa.k
    << " fm^-1, eta_aA = " << v->dist.aa.eta << ", kappa = " << v->dist.kappa << " fm^-1, eta_b = "
    << v->dist.etaB << "; channel radius " << v->radius << " fm, l =";
  for (int lv : v->lvals) l << " " << lv;
  l << "; L <= " << v->laMax << " (a + A), " << v->lsMax << " (s + F); u to " << v->uMax << " fm on "
    << v->uNodes << " x " << v->cNodes << " nodes; grid " << v->gridLo << " to "
    << v->gridLo + (v->nE - 1) * v->gridStep << " MeV (" << v->nE << " energies, " << v->nNodes
    << " node(s) each), " << v->buildSeconds << " s.\n"
    << "  The distortion factor R(E) is not applied: the DW vertex carries the energy dependence.";
  for (const std::string &w : v->dist.warnings) l << "\nWARNING: <thm> experiment[" << x.name << "]: " << w;
  for (int s : group.segments)
    if (data.GetSegment(s)->GetThmWeight())
      l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << data.GetSegment(s)->GetSegmentKey()
        << " also has weight[" << data.GetSegment(s)->GetSegmentKey()
        << "]=; it multiplies the model with the DW vertex.";
  configure.outStream << l.str() << std::endl;
  group.dwVertex = v;
  for (int s : group.segments) data.GetSegment(s)->SetThmDwVertex(v);
  return true;
}

// Distortion factor R(E) (ThmDistortion.h): multiplies the model of every
// segment of the experiment before the folding.
bool ThmBuildDistortion(EData &data, const Config &configure, const ThmExperiment &x, EData::ThmGroup &group,
                        const std::string &where, const ThmDistortion::Kinematics &dk) {
  double eLo, eHi;
  ThmGroupEnergyRange(data, group, eLo, eHi);
  std::shared_ptr<ThmDistortion> d = std::make_shared<ThmDistortion>();
  d->experiment = x.name;
  std::ostringstream l;
  l.precision(6);
  if (x.distortion == ThmExperiment::DIST_TABLE) {
    d->kind = ThmDistortion::TABLE;
    d->table = x.distortionWeights;
    const ThmWeightTable &table = *d->table;
    if (!table.Covers(eLo) || !table.Covers(eHi)) {
      configure.outStream << where << "distortion: the points span E_cm = " << eLo << " to " << eHi
                          << " MeV, beyond the table '" << table.name << "' [" << table.e.front() << ", "
                          << table.e.back() << "] MeV." << std::endl;
      return false;
    }
    d->description = "table " + table.name;
    l << "  Distortion factor: w(E) from the table '" << table.name << "' multiplies the model.";
  } else {
    d->kin = dk;
    d->eAA = dk.beamEnergy * dk.mTarget / (dk.mBeam + dk.mTarget);
    // Every data point must be reachable (the grid then extends 0.5
    // MeV beyond the points, short of the spectator's threshold).
    d->angleKind = x.angleKind == 1 ? ThmDistortion::LAB : x.angleKind == 2 ? ThmDistortion::CM : ThmDistortion::QF;
    d->angle = x.angle;
    d->sf.kind = x.distortion == ThmExperiment::DIST_OPTICAL && x.opticalSF.kind == 0
                     ? ThmDistortion::Channel::PLANE
                     : ThmDistortion::Channel::POINT_COULOMB;
    d->sf.mu = dk.ms * (dk.mx + dk.mA) / (dk.ms + dk.mx + dk.mA) * uconv;
    d->vcm = std::sqrt(2.0 * dk.mBeam * uconv * dk.beamEnergy) / ((dk.mBeam + dk.mTarget) * uconv);
    for (int s : group.segments)
      for (int p = 1; p <= data.GetSegment(s)->NumPoints(); p++) {
        std::string why = d->CheckEnergy(data.GetSegment(s)->GetPoint(p)->GetCMEnergy());
        if (!why.empty()) {
          configure.outStream << where << "distortion: " << why << "." << std::endl;
          return false;
        }
      }
    double gridHi = std::min(eHi + 0.5, d->eAA - dk.bind - 0.5 * d->EsF(eHi));
    d->dataLo = eLo;
    d->dataHi = eHi;
    std::string why = d->Build(x, dk, eLo - 0.5, gridHi, 0.5 * (eLo + eHi));
    // spectatorAngles=: every data point needs an accepted direction.
    for (int s = 0; why.empty() && d->angWindow && s < (int)group.segments.size(); s++)
      for (int p = 1; why.empty() && p <= data.GetSegment(group.segments[s])->NumPoints(); p++)
        why = d->CheckWindow(data.GetSegment(group.segments[s])->GetPoint(p)->GetCMEnergy());
    if (!why.empty()) {
      configure.outStream << where << "distortion: " << why << "." << std::endl;
      return false;
    }
    ThmDistortion::Point lo = d->Evaluate(eLo), hi = d->Evaluate(eHi);
    l << "  Distortion factor R(E), zero-range DWBA: " << d->description << ".\n"
      << "  k_aA = " << d->aa.k << " fm^-1, eta_aA = " << d->aa.eta << ", kappa = " << d->kappa
      << " fm^-1, eta_b = " << d->etaB << ", beta = m_s/m_a = " << d->beta << "; E_ref = " << d->eRef
      << " MeV, grid " << d->gridLo << " to " << d->gridLo + (d->lnR.size() - 1) * d->gridStep
      << " MeV, radial step " << d->h << " fm to " << (d->n - 1) * d->h << " fm, l <= " << d->uAA.size() - 1
      << ".\n"
      << "  R = " << d->R(lo) << " at E = " << eLo << " MeV (E_sF = " << lo.esf << ", eta_sF = " << lo.etasf
      << ", theta_cm = " << lo.thetaCm << " deg), " << d->R(hi) << " at E = " << eHi << " MeV (E_sF = "
      << hi.esf << ", eta_sF = " << hi.etasf << ", theta_cm = " << hi.thetaCm << " deg).";
    if (d->angWindow)
      l << "\n  Spectator directions: " << lo.nodes << " accepted node(s) at the lowest point, " << hi.nodes
        << " at the highest; theta_cm above is their acceptance-weighted mean.";
    if (d->pwSignChange && d->ratioPW)
      l << "\nWARNING: <thm> experiment[" << x.name << "]: the plane-wave amplitude M_PW changes sign on "
           "the grid (a node of the momentum distribution at this angle); R = |M/M_PW|^2 is singular "
           "there (distortionRatio=dw avoids it).";
    for (const std::string &w : d->warnings) l << "\nWARNING: <thm> experiment[" << x.name << "]: " << w;
    if (d->tailWorst > 1.0e-8)
      l << "\nWARNING: <thm> experiment[" << x.name << "]: the radial integrals are cut at r = "
        << (d->n - 1) * d->h << " fm with a remainder up to " << d->tailWorst << " of |M|.";
  }
  for (int s : group.segments)
    if (data.GetSegment(s)->GetThmWeight())
      l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << data.GetSegment(s)->GetSegmentKey()
        << " also has weight[" << data.GetSegment(s)->GetSegmentKey()
        << "]=; both multiply its model (the distortion factor and the weight).";
  configure.outStream << l.str() << std::endl;
  group.distortion = d;
  for (int s : group.segments) data.GetSegment(s)->SetThmDistortion(d);
  return true;
}

// Angular window of the exit pair (ThmAngular.h), theta=: the model of every
// segment is dsigma/dOmega averaged over theta_cm in the window.
bool ThmBuildAngleWindow(EData &data, const Config &configure, CNuc *theCNuc, const ThmExperiment &x,
                         EData::ThmGroup &group, const std::string &where) {
  if (configure.thm.coherentL) {
    configure.outStream << where << "theta= computes the interference of the entrance partial waves "
                           "exactly (at fixed angle they interfere); entranceL=coherent is an approximation "
                           "of the angle-integrated observable and cannot be combined with it."
                        << std::endl;
    return false;
  }
  int maxLp = 0;
  for (int j = 1; j <= theCNuc->NumJGroups(); j++)
    for (int ch = 1; ch <= theCNuc->GetJGroup(j)->NumChannels(); ch++)
      maxLp = std::max(maxLp, theCNuc->GetJGroup(j)->GetChannel(ch)->GetL());
  if (2 * maxLp > ThmAngleWindow::kMaxL) {
    configure.outStream << where << "theta= carries Legendre orders up to " << ThmAngleWindow::kMaxL
                        << "; the model has a channel with l = " << maxLp << "." << std::endl;
    return false;
  }
  std::shared_ptr<ThmAngleWindow> w = std::make_shared<ThmAngleWindow>();
  w->experiment = x.name;
  BuildThmAngleWindow(x.thetaMin, x.thetaMax, *w);
  std::ostringstream a;
  a.precision(6);
  a << "  Angular window: theta_cm = " << x.thetaMin << "-" << x.thetaMax
    << " deg (exit particle 1 relative to 2, from p_xA = entrance particle 1 relative to 2); the model "
       "is the HOES dsigma/dOmega averaged over it (4 pi times it is the angle-integrated cross section "
       "for 0-180).";
  configure.outStream << a.str() << std::endl;
  group.angle = w;
  for (int s : group.segments) data.GetSegment(s)->SetThmAngleWindow(w);
  return true;
}

}  // namespace

/*
 * One group per experiment line, in order; per group, in this order: its
 * segments, the profile, the kinematics, the ps window, the line shape, the
 * DW vertex or the distortion factor, the angular window and the coherent
 * background.  The first refusal ends it (the order of the checks is what a
 * user sees, so it is kept).
 */
int EData::BuildThmGroups(const Config &configure, CNuc *theCNuc, int numLines) {
  thmGroups_.clear();
  thmCoherentParams_.clear();
  thmCoherentParamOffset_ = -1;
  for (const ThmExperiment &x : configure.thm.experiments) {
    const std::string where = "ERROR: <thm> experiment[" + x.name + "]: ";
    ThmGroup group;
    group.name = x.name;
    group.terms = x.backgroundTerms;
    for (int key : x.segments) {
      if (key > numLines) {
        configure.outStream << where << "segment " << key << ": <segmentsData> has only " << numLines
                            << " line(s)." << std::endl;
        return -1;
      }
      int index = 0;
      for (int s = 1; s <= NumSegments(); s++)
        if (GetSegment(s)->GetSegmentKey() == key) index = s;
      if (index == 0) {
        configure.outStream << "WARNING: <thm> experiment[" << x.name << "]: segment line " << key
                            << " of <segmentsData> is not in use; it is left out." << std::endl;
        continue;
      }
      ESegment *segment = GetSegment(index);
      if (!segment->IsTHM()) {
        configure.outStream << where << "segment " << key << " is not a THM segment (isDiff < 10)." << std::endl;
        return -1;
      }
      if (!segment->IsVaryNorm()) {
        configure.outStream << where << "segment " << key
                            << " has a fixed norm; the segments of an experiment share one free "
                               "(profiled) norm, so free it."
                            << std::endl;
        return -1;
      }
      group.segments.push_back(index);
    }
    if (group.segments.empty()) {
      configure.outStream << "WARNING: <thm> experiment[" << x.name << "] has no segment in use; ignored."
                          << std::endl;
      continue;
    }
    int points = 0;
    for (int s : group.segments)
      for (int p = 1; p <= GetSegment(s)->NumPoints(); p++)
        if (GetSegment(s)->GetPoint(p)->GetCMCrossSectionError() != 0.0) points++;
    if (points <= 1 + group.terms) {
      configure.outStream << where << points << " point(s) with an error for " << 1 + group.terms
                          << " profiled linear parameter(s) (norm and background " << ThmExperiment::BackgroundName(group.terms)
                          << "); it needs more." << std::endl;
      return -1;
    }
    group.trivial = group.segments.size() == 1 && group.terms == 0;
    group.pairKey = GetSegment(group.segments[0])->GetEntranceKey();
    for (int s : group.segments)
      if (GetSegment(s)->GetEntranceKey() != group.pairKey) group.pairKey = 0;

    std::ostringstream summary;
    summary << "THM experiment '" << x.name << "': segment" << (group.segments.size() > 1 ? "s " : " ");
    for (size_t k = 0; k < group.segments.size(); k++)
      summary << (k ? "," : "") << GetSegment(group.segments[k])->GetSegmentKey();
    summary << (group.segments.size() > 1 ? " share one profiled norm" : " with a profiled norm")
            << ", background " << ThmExperiment::BackgroundName(group.terms) << ".";
    configure.outStream << summary.str() << std::endl;

    ThmReaction reaction;  // its kinematics are those of distortion=coulomb|optical and the DW vertex
    if (x.hasKinematics) {
      if (!ThmResolveKinematics(*this, configure, theCNuc, x, group, where, reaction)) return -1;
      if (x.psKind != ThmExperiment::PS_DELTA && !ThmBuildWindow(*this, configure, x, group, where, reaction))
        return -1;
      if (x.lineshape && !ThmBuildLineshape(*this, configure, theCNuc, x, group, where, reaction)) return -1;
    }
    if (x.vertexDW) {
      if (!ThmBuildDwVertex(*this, configure, theCNuc, x, group, where, reaction.dk)) return -1;
    } else if (x.distortion != ThmExperiment::DIST_NONE) {
      if (!ThmBuildDistortion(*this, configure, x, group, where, reaction.dk)) return -1;
    }
    if (x.hasTheta && !ThmBuildAngleWindow(*this, configure, theCNuc, x, group, where)) return -1;
    if (!x.cbackground.empty()) {
      // Coherent background (cbackground=, ThmCoherentBackground): one complex
      // amplitude per (J^pi, entrance (s,l), exit (s',l')) combination, whose
      // real and imaginary parts are fit parameters.
      if (configure.thm.coherentL) {
        configure.outStream << where << "cbackground= adds an amplitude per entrance bucket (s, l); "
                               "entranceL=coherent merges the l of a channel spin, so it cannot be combined "
                               "with it."
                            << std::endl;
        return -1;
      }
      if (group.pairKey == 0) {
        configure.outStream << where << "cbackground= needs one entrance pair; the segments have different ones."
                            << std::endl;
        return -1;
      }
      const int aa = theCNuc->GetPairNumFromKey(group.pairKey);
      std::shared_ptr<ThmCoherentBackground> cb = std::make_shared<ThmCoherentBackground>();
      cb->experiment = x.name;
      std::ostringstream c;
      c << "  Coherent background (cbackground=): ";
      for (const ThmExperiment::CoherentTerm &t : x.cbackground) {
        const std::string jpi = ThmSpinText(t.J) + (t.parity > 0 ? "+" : "-");
        const std::string term = "cbackground " + jpi + ":" + std::to_string(t.exitKey) + ": ";
        bool exitUsed = false;
        for (int s : group.segments) exitUsed = exitUsed || GetSegment(s)->GetExitKey() == t.exitKey;
        if (!exitUsed || !theCNuc->IsPairKey(t.exitKey)) {
          configure.outStream << where << term << "no segment of the experiment has exit pair " << t.exitKey << "."
                              << std::endl;
          return -1;
        }
        const int exitNum = theCNuc->GetPairNumFromKey(t.exitKey);
        int j = 0;
        for (int jj = 1; jj <= theCNuc->NumJGroups(); jj++) {
          JGroup *jg = theCNuc->GetJGroup(jj);
          if (jg->IsInRMatrix() && std::fabs(jg->GetJ() - t.J) < 1.0e-6 && jg->GetPi() == t.parity) j = jj;
        }
        if (j == 0) {
          configure.outStream << where << term << "the model has no J^pi = " << jpi
                              << " group with a level in the R matrix." << std::endl;
          return -1;
        }
        JGroup *jg = theCNuc->GetJGroup(j);
        int found = 0;
        for (int in = 1; in <= jg->NumChannels(); in++) {
          AChannel *ci = jg->GetChannel(in);
          if (ci->GetPairNum() != aa) continue;
          if (t.hasChannels && (std::fabs(ci->GetS() - t.s) > 1.0e-6 || ci->GetL() != t.l)) continue;
          for (int out = 1; out <= jg->NumChannels(); out++) {
            AChannel *co = jg->GetChannel(out);
            if (co->GetPairNum() != exitNum) continue;
            if (t.hasChannels && (std::fabs(co->GetS() - t.sp) > 1.0e-6 || co->GetL() != t.lp)) continue;
            for (const ThmCoherentBackground::Combo &o : cb->combos)
              if (o.jGroup == j && o.entrance == in && o.exit == out) {
                configure.outStream << where << term << "the combination (s,l) = (" << ThmSpinText(ci->GetS())
                                    << "," << ci->GetL() << ") -> (s',l') = (" << ThmSpinText(co->GetS()) << ","
                                    << co->GetL() << ") is given twice." << std::endl;
                return -1;
              }
            ThmCoherentBackground::Combo combo;
            combo.jGroup = j;
            combo.entrance = in;
            combo.exit = out;
            combo.form = t.form;
            combo.index = (int)thmCoherentParams_.size();
            cb->combos.push_back(combo);
            const std::string stem = ThmCoherentParamStem(x.name, t.J, t.parity, t.exitKey, ci->GetS(), ci->GetL(),
                                                          co->GetS(), co->GetL());
            static const char *parts[4] = {"re0", "im0", "re1", "im1"};
            for (int k = 0; k < 2 * t.form; k++) {
              ThmCoherentParam par;
              par.name = stem + parts[k];
              par.experiment = x.name;
              par.value = t.hasValues ? t.value[k] : 0.0;
              par.fixed = t.hasValues && t.fixed[k];
              thmCoherentParams_.push_back(par);
            }
            c << (found || cb->combos.size() > 1 ? "; " : "") << jpi << " (" << ThmSpinText(ci->GetS()) << ","
              << ci->GetL() << ") -> " << t.exitKey << " (" << ThmSpinText(co->GetS()) << "," << co->GetL()
              << ") " << (t.form == 1 ? "const" : "linear");
            found++;
          }
        }
        if (!found) {
          configure.outStream << where << term
                              << (t.hasChannels ? "the J^pi group has no entrance channel (s,l) = (" + ThmSpinText(t.s) +
                                                      "," + std::to_string(t.l) + ") with an exit channel (s',l') = (" +
                                                      ThmSpinText(t.sp) + "," + std::to_string(t.lp) + ")."
                                                : std::string("the J^pi group does not couple the entrance pair to the "
                                                              "exit pair."))
                              << std::endl;
          return -1;
        }
      }
      c << ".  The amplitude c(E) M_l is added to the resonant HOES amplitude of each before squaring; "
           "Re and Im of c are fit parameters (cbkg_*).";
      configure.outStream << c.str() << std::endl;
      group.coherent = cb;
      for (int s : group.segments) GetSegment(s)->SetThmCoherent(cb);
    }
    thmGroups_.push_back(group);
  }
  // The coherent backgrounds are the last block of the parameter vector
  // (FillMnParams sets the offset again when it lays the vector out).
  if (!thmCoherentParams_.empty()) {
    int offset = 0;
    for (int j = 1; j <= theCNuc->NumJGroups(); j++)
      offset += theCNuc->GetJGroup(j)->NumLevels() * (1 + theCNuc->GetJGroup(j)->NumChannels());
    for (int s = 1; s <= NumSegments(); s++) {
      if (GetSegment(s)->IsVaryNorm() && !GetSegment(s)->IsProfiledNorm()) offset++;
      offset++;  // energy shift
    }
    thmCoherentParamOffset_ = offset;
    std::vector<double> v;
    for (const ThmCoherentParam &par : thmCoherentParams_) v.push_back(par.value);
    theCNuc->SetThmCoherentValues(v);
  }
  return 0;
}

int EData::ThmGroupOf(int i) {
  for (size_t g = 0; g < thmGroups_.size(); g++) {
    if (thmGroups_[g].trivial) continue;
    const std::vector<int> &s = thmGroups_[g].segments;
    if (std::find(s.begin(), s.end(), i) != s.end()) return GetSegment(i)->IsProfiledNorm() ? (int)g : -1;
  }
  return -1;
}

bool EData::IsLastOfThmGroup(int g, int i) {
  const std::vector<int> &s = thmGroups_[g].segments;
  for (size_t k = s.size(); k-- > 0;)
    if (GetSegment(s[k])->IsProfiledNorm()) return s[k] == i;
  return false;
}

namespace {
// The points of the profiled segments of a group, concatenated in order.
void GatherThmGroup(EData *data, const std::vector<int> &segments, std::vector<double> &m,
                    std::vector<double> &d, std::vector<double> &e, std::vector<double> &energy) {
  m.clear();
  d.clear();
  e.clear();
  energy.clear();
  for (int s : segments) {
    ESegment *seg = data->GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (!pt) continue;
      m.push_back(pt->GetFitCrossSection());
      d.push_back(pt->GetCMCrossSection());
      e.push_back(pt->GetCMCrossSectionError());
      energy.push_back(pt->GetCMEnergy());
    }
  }
}
}  // namespace

double EData::ProfileThmGroup(int g) {
  ThmGroup &group = thmGroups_[g];
  std::vector<double> m, d, e, energy;
  GatherThmGroup(this, group.segments, m, d, e, energy);
  group.profile = SolveThmProfile(m, d, e, energy, group.terms);
  const ThmProfile &p = group.profile;
  double total = 0.0;
  size_t k = 0;
  for (int s : group.segments) {
    ESegment *seg = GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    seg->SetNorm(p.Norm());
    double chi = 0.0;
    for (int q = 1; q <= seg->NumPoints(); q++) {
      if (!seg->GetPoint(q)) continue;
      double r = p.Residual(m[k], d[k], e[k], energy[k]);
      chi += r * r;
      k++;
    }
    seg->SetSegmentChiSquared(chi);
    total += chi;
  }
  return total;
}

void EData::ThmGroupResiduals(int g, int i, std::vector<double> &out) {
  out.clear();
  const ThmProfile &p = thmGroups_[g].profile;
  ESegment *seg = GetSegment(i);
  for (int q = 1; q <= seg->NumPoints(); q++) {
    EPoint *pt = seg->GetPoint(q);
    if (!pt) continue;
    out.push_back(p.Residual(pt->GetFitCrossSection(), pt->GetCMCrossSection(), pt->GetCMCrossSectionError(),
                             pt->GetCMEnergy()));
  }
}

double EData::ThmBackgroundAt(int i, double energy) {
  int g = ThmGroupOf(i);
  if (g < 0 || thmGroups_[g].terms == 0) return 0.0;
  return thmGroups_[g].profile.Background(energy);
}

std::vector<ThmExperimentReport> EData::ThmExperimentReports() {
  std::vector<ThmExperimentReport> out;
  for (size_t g = 0; g < thmGroups_.size(); g++) {
    ThmGroup &group = thmGroups_[g];
    ThmExperimentReport r;
    r.name = group.name;
    r.background = ThmExperiment::BackgroundName(group.terms);
    ThmProfile p;
    if (group.trivial || ThmGroupOf(group.segments[0]) < 0) {
      // Profiled per segment (ESegment::ProfileNormChiSquared) or not at all:
      // the same closed form, recomputed here for its uncertainty.
      std::vector<double> m, d, e, energy;
      GatherThmGroup(this, group.segments, m, d, e, energy);
      p = SolveThmProfile(m, d, e, energy, group.terms);
    } else {
      p = group.profile;
    }
    r.points = 0;
    r.chi2 = 0.0;
    for (int s : group.segments) {
      ESegment *seg = GetSegment(s);
      if (!seg->IsProfiledNorm()) continue;
      r.segments.push_back(seg->GetSegmentKey());
      r.chi2 += seg->GetSegmentChiSquared();
    }
    r.points = p.points;
    p.Reported(r.value, r.covariance);
    r.status = r.segments.empty() ? "not profiled: no segment of it has a free norm" : p.status;
    for (const ThmCoherentParam &c : thmCoherentParams_)
      if (c.experiment == group.name) {
        r.coherentNames.push_back(c.name);
        r.coherentValues.push_back(c.value);
        r.coherentFixed.push_back(c.fixed);
      }
    out.push_back(r);
  }
  return out;
}

void EData::WriteThmExperiments(const Config &configure) {
  std::string file = configure.outputdir + "thm_experiments.out";
  std::ofstream out(file.c_str());
  if (!out) {
    configure.outStream << "Could not write " << file << "." << std::endl;
    return;
  }
  out << "# THM experiments (<thm> experiment[<name>] lines): segments sharing one profiled norm and\n"
         "# a background b(E) = b0 + b1 E + b2 E^2 added to the folded model (model units, E the c.m.\n"
         "# energy of the entrance pair in MeV).  norm multiplies the data, as the segment norms of\n"
         "# normalizations.out.  Uncertainties and covariance are those of the closed-form profile\n"
         "# at fixed R-matrix parameters (inverse normal matrix of the linear least squares, not\n"
         "# scaled by chi2/nu).\n";
  out.precision(10);
  out << std::scientific;
  for (const ThmExperimentReport &r : ThmExperimentReports()) {
    out << "\nexperiment: " << r.name << "\nsegments:";
    for (size_t k = 0; k < r.segments.size(); k++) out << (k ? "," : " ") << r.segments[k];
    out << "\nbackground: " << r.background << "\npoints: " << r.points << "\nchi2: " << r.chi2
        << "\nstatus: " << r.status << "\n";
    const char *names[4] = {"norm", "b0", "b1", "b2"};
    int q = 1 + (r.background == "none" ? 0 : r.background == "const" ? 1 : r.background == "linear" ? 2 : 3);
    for (int i = 0; i < q; i++)
      out << std::left << std::setw(6) << names[i] << std::right << std::setw(18) << r.value[i]
          << std::setw(18) << std::sqrt(std::max(0.0, r.covariance[i * 4 + i])) << "\n";
    out << "covariance:\n";
    for (int i = 0; i < q; i++) {
      out << std::left << std::setw(6) << names[i] << std::right;
      for (int j = 0; j < q; j++) out << std::setw(18) << r.covariance[i * 4 + j];
      out << "\n";
    }
    // Coherent background (cbackground=): its fit parameters (their errors
    // are the fit's, in param.par).
    if (!r.coherentNames.empty()) {
      out << "cbackground: " << r.coherentNames.size() << " parameter(s); c(E) = c0 + c1 E times the vertex M_l "
          << "is added to the resonant amplitude\n";
      for (size_t k = 0; k < r.coherentNames.size(); k++)
        out << "cbkg " << std::left << std::setw(40) << r.coherentNames[k] << std::right << std::setw(18)
            << r.coherentValues[k] << (r.coherentFixed[k] ? "  fixed" : "") << "\n";
    }
    // Spectator-momentum window (ps=...): the nodes.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.window) continue;
      const ThmSpectatorWindow &w = *group.window;
      out << "ps: " << w.description << "; mu_sx = " << w.muSx << " MeV\n"
          << "# The model at E is the average of the HOES cross section over the accepted spectator\n"
          << "# directions, weight |phi(p_s)|^2 d cos(theta_cm) (x the acceptance), incoherent; at fixed E\n"
          << "# the direction fixes p_s, and node k adds T_s = p_s^2/2mu_sx to E + B in the vertex.  The\n"
          << "# nodes at the lowest and the highest point; columns: E (MeV), theta_cm (deg), p_s (MeV/c),\n"
          << "# weight, T_s (MeV).\n";
      for (double e : {w.dataE.front(), w.dataE.back()}) {
        std::vector<ThmSpectatorWindow::Node> nodes;
        w.NodesAt(e, nodes);
        for (const ThmSpectatorWindow::Node &k : nodes)
          out << "ps_node" << std::setw(18) << e << std::setw(18) << k.theta << std::setw(18) << k.p
              << std::setw(18) << k.weight << std::setw(18) << k.es << "\n";
        out << std::left << std::setw(16) << "<T_s>" << std::right << std::setw(18) << e << std::setw(18)
            << w.MeanEs(e) << "\n";
      }
      // Memory of the per-point node tables (EPoint::ThmPsTable): points
      // including the folding sub-points, the stored entrance channels summed
      // over the points (x nodes = entries), and their bytes.
      size_t tablePoints = 0, tableBytes = 0;
      for (int s : group.segments)
        for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
          EPoint *point = GetSegment(s)->GetPoint(p);
          tablePoints++;
          tableBytes += point->ThmPsTableBytes();
          for (int q = 1; q <= point->NumSubPoints(); q++) {
            tablePoints++;
            tableBytes += point->GetSubPoint(q)->ThmPsTableBytes();
          }
        }
      out << "# Node tables of the vertex: points (with folding sub-points), bytes.\n"
          << "ps_table" << std::setw(18) << (double)tablePoints << std::setw(18) << (double)tableBytes << "\n";
    }
    // Distortion factor (distortion=...): R at the lowest point, E_ref and the highest point.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.distortion) continue;
      const ThmDistortion &d = *group.distortion;
      out << "distortion: " << d.description << "\n";
      if (d.kind == ThmDistortion::TABLE) continue;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "# Zero-range DWBA transfer amplitude M(E) = <chi(-)_sF phi_sx chi(+)_aA(beta r)> (Mukhamedzhanov &\n"
          << "# Pang PRC 99 (2019) 064618 eqs. 20-24; Mukhamedzhanov arXiv:2609.04498 eqs. 22-30), M_PW its\n"
          << "# plane-wave limit; the model is multiplied by R(E) before folding (PWA-extracted S* / R).\n"
          << "# E_aA = " << d.eAA << " MeV, B = " << d.kin.bind << " MeV, k_aA = " << d.aa.k << " fm^-1, eta_aA = "
          << d.aa.eta << ", kappa = " << d.kappa << " fm^-1, eta_b = " << d.etaB << ", beta = " << d.beta << "\n"
          << (d.angWindow ? "# spectatorAngles: |M|^2 and |M_PW|^2 are averages over the accepted directions\n"
                            "# (weight d cos theta_cm x acceptance), theta_cm their weighted mean.\n"
                          : "")
          << "# Columns: E, E_sF (MeV), eta_sF, theta_cm (deg), |M|^2, |M_PW|^2, R, l_max.\n";
      for (double e : {eLo, d.eRef, eHi}) {
        ThmDistortion::Point p = d.Evaluate(e);
        out << "distortion_point" << std::setw(18) << e << std::setw(18) << p.esf << std::setw(18) << p.etasf
            << std::setw(18) << p.thetaCm << std::setw(18) << ThmDistortion::M2(p) << std::setw(18)
            << ThmDistortion::MPW2(p) << std::setw(18) << (p.ok ? d.R(p) : 0.0) << std::setw(6) << p.lmax << "\n";
      }
    }
    // Distorted-wave entrance vertex (vertexModel=dw): the Gram entries at the
    // lowest, the middle and the highest point.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.dwVertex) continue;
      const ThmDwVertex &v = *group.dwVertex;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "vertex: " << v.description << "\n"
          << "# Surface term of the prior-form DWBA (Mukhamedzhanov PRC 84 (2011) 044616; Mukhamedzhanov,\n"
          << "# Kadyrov & Pang EPJA 56 (2020) 233 eqs. 28-32) replaces M_l; R(E) is not applied.\n"
          << "# G = (4pi/(2l+1)) sum_m (s_m, d_m)^+ (s_m, d_m), (s_m, d_m) = (S_lm(a), a S_lm'(a))/(4pi phi~(q)),\n"
          << "# |M_l|^2 = c^+ G c with c = (B - 1, -1); plane waves: G11 = j_l(pa)^2, G22 = (pa j_l'(pa))^2.\n"
          << "# k_aA = " << v.dist.aa.k << " fm^-1, eta_aA = " << v.dist.aa.eta << ", kappa = " << v.dist.kappa
          << " fm^-1, alpha = " << v.alpha << ", beta = " << v.beta << ", a = " << v.radius << " fm\n"
          << "# Columns: E (MeV), q (MeV/c), p a, l, G11, G22, Re G12, Im G12 (spectatorAngle direction).\n";
      for (double e : {eLo, 0.5 * (eLo + eHi), eHi}) {
        std::vector<double> w, q, g, gd;
        double qd, pd;
        v.Interpolate(e, w, q, g, gd, qd, pd);
        for (size_t li = 0; li < v.lvals.size(); li++)
          out << "dw_vertex_point" << std::setw(18) << e << std::setw(18) << qd * hbarc << std::setw(18)
              << pd * v.radius << std::setw(4) << v.lvals[li] << std::setw(18) << gd[li * 4] << std::setw(18)
              << gd[li * 4 + 1] << std::setw(18) << gd[li * 4 + 2] << std::setw(18) << gd[li * 4 + 3] << "\n";
      }
    }
    // Coulomb line shape (lineshape=on): the ranges over the experiment's points.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.lineshape) continue;
      const ThmLineshape &ls = *group.lineshape;
      double eLo, eHi;
      ThmGroupEnergyRange(*this, group, eLo, eHi);
      out << "lineshape: on (spectator " << ls.spectator << ", Z_s = " << ls.Zs << ", Z_F = " << ls.ZF
          << "; E_aA = " << ls.eAA << " MeV, B = " << ls.bind << " MeV)\n"
          << "# |N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)] per level, zeta = eta_sB - eta_0\n"
          << "# (Mukhamedzhanov et al. EPJA 56 (2020) 233 eqs. 56-62, case 2: m_B >> m_s, m_b, eta_sb\n"
          << "# neglected); zeta < 0 moves the peaks up in E.  eta_sb is the neglected s-b term, averaged\n"
          << "# over the b direction: the approximation assumes |eta_sb| << 1.  Rows: the value at the\n"
          << "# lowest and at the highest point energy E (MeV); [k] is the exit pair key.\n";
      auto row = [&](const char *name, double a, double b) {
        out << std::left << std::setw(16) << name << std::right << std::setw(18) << a << std::setw(18) << b << "\n";
      };
      row("E", eLo, eHi);
      row("E_sF", ls.EsF(eLo), ls.EsF(eHi));
      row("eta_0", ls.Eta0(eLo), ls.Eta0(eHi));
      for (const ThmLineshape::Exit &e : ls.exits) {
        std::ostringstream key;
        key << "zeta[" << e.pairKey << "]";
        row(key.str().c_str(), ls.Zeta(eLo, e.ZB, e.mB), ls.Zeta(eHi, e.ZB, e.mB));
        key.str("");
        key << "eta_sb[" << e.pairKey << "]";
        row(key.str().c_str(), ls.EtaSbEstimate(eLo, eLo + e.q, e.Zb, e.mb, e.mB),
            ls.EtaSbEstimate(eHi, eHi + e.q, e.Zb, e.mb, e.mB));
      }
    }
  }
}

bool EData::ThmLineshapeTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                              const Config &configure, ThmLineshapeReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->lineshape) {
    why = "THM experiment '" + name + "' has no line shape (lineshape=on)";
    return false;
  }
  const ThmLineshape &ls = *group->lineshape;
  out = ThmLineshapeReport();
  out.experiment = name;
  out.spectator = ls.spectator;
  out.Zs = ls.Zs;
  out.ZF = ls.ZF;
  out.eAA = ls.eAA;
  out.bind = ls.bind;
  out.energy = energies;
  for (double e : energies) {
    out.esf.push_back(ls.EsF(e));
    out.eta0.push_back(ls.Eta0(e));
  }
  int entranceKey = GetSegment(group->segments[0])->GetEntranceKey();
  int aa = compound->GetPairNumFromKey(entranceKey);
  PPair *entrance = compound->GetPair(aa);
  double threshold = entrance->GetSepE() + entrance->GetExE();
  for (const ThmLineshape::Exit &x : ls.exits) {
    ThmLineshapeReport::Exit ex;
    ex.pairKey = x.pairKey;
    ex.Zb = x.Zb;
    ex.ZB = x.ZB;
    ex.mb = x.mb;
    ex.mB = x.mB;
    for (double e : energies) {
      ex.zeta.push_back(ls.Zeta(e, x.ZB, x.mB));
      ex.etaSb.push_back(ls.EtaSbEstimate(e, e + x.q, x.Zb, x.mb, x.mB));
    }
    int exitNum = compound->GetPairNumFromKey(x.pairKey);
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      JGroup *jg = compound->GetJGroup(j);
      if (!jg->IsInRMatrix()) continue;
      bool in = false, outCh = false;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        in = in || jg->GetChannel(ch)->GetPairNum() == aa;
        outCh = outCh || jg->GetChannel(ch)->GetPairNum() == exitNum;
      }
      if (!in || !outCh) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmLineshapeReport::Level lv;
        lv.jgroup = j;
        lv.level = la;
        lv.J = jg->GetJ();
        lv.pi = jg->GetPi();
        lv.energy = level->GetFitE() - threshold;
        lv.width = ThmLevelWidth(compound, jg, level, configure);
        for (size_t k = 0; k < energies.size(); k++)
          lv.nc2.push_back(ThmLineshapeFactorSq(ex.zeta[k], lv.energy - energies[k], lv.width));
        ex.levels.push_back(lv);
      }
    }
    out.exits.push_back(ex);
  }
  return true;
}

bool EData::ThmDistortionTable(const std::string &name, const std::vector<double> &energies, ThmDistortionReport &out,
                               std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->distortion) {
    why = "THM experiment '" + name + "' has no distortion (distortion=coulomb|optical|table:<file>)";
    return false;
  }
  const ThmDistortion &d = *group->distortion;
  out = ThmDistortionReport();
  out.experiment = name;
  out.description = d.description;
  out.energy = energies;
  if (d.kind == ThmDistortion::TABLE) {
    out.kind = "table";
    for (double e : energies) out.rModel.push_back(d.Weight(e));
    return true;
  }
  out.kind = d.kind == ThmDistortion::COULOMB ? "coulomb" : "optical";
  out.eRef = d.eRef;
  out.eAA = d.eAA;
  out.bind = d.kin.bind;
  out.kAA = d.aa.k;
  out.etaAA = d.aa.eta;
  out.kappa = d.kappa;
  out.etaB = d.etaB;
  out.beta = d.beta;
  for (double e : energies) {
    ThmDistortion::Point p;
    std::string bad = d.CheckEnergy(e);
    if (bad.empty()) p = d.Evaluate(e);
    if (!bad.empty() || !p.ok) {
      why = "THM experiment '" + name + "': " + (bad.empty() ? p.why : bad);
      return false;
    }
    out.esf.push_back(p.esf);
    out.ksf.push_back(p.ksf);
    out.etasf.push_back(p.etasf);
    out.thetaCm.push_back(p.thetaCm);
    out.x.push_back(p.x);
    out.q.push_back(p.q);
    out.m2.push_back(ThmDistortion::M2(p));
    out.mpw2.push_back(ThmDistortion::MPW2(p));
    out.r.push_back(d.R(p));
    out.rModel.push_back(d.Weight(e));
    out.lmax.push_back(p.lmax);
  }
  return true;
}

bool EData::ThmVertexTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                           const Config &configure, ThmVertexReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (group->pairKey == 0 || !compound->IsPairKey(group->pairKey)) {
    why = "THM experiment '" + name + "': its segments have different entrance pairs";
    return false;
  }
  const bool useGSL = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int aa = compound->GetPairNumFromKey(group->pairKey);
  PPair *pair = compound->GetPair(aa);
  if (pair->GetPType() != 0) {
    why = "THM experiment '" + name + "': the entrance pair is not a particle pair";
    return false;
  }
  out = ThmVertexReport();
  out.experiment = name;
  out.pairKey = group->pairKey;
  out.bind = pair->GetBindingEnergy();
  out.radius = pair->GetChRad();
  const ThmDwVertex *dw = group->dwVertex.get();
  out.energy = energies;
  if (group->window) {
    const ThmSpectatorWindow &w = *group->window;
    out.window = w.description;
    out.muSx = w.muSx;
    out.windowObject = group->window;
    for (double e : energies) {
      std::vector<ThmSpectatorWindow::Node> nodes;
      bool moved = false;
      w.Nodes(e, nodes, &moved);
      out.reached.push_back(moved ? 0 : 1);
      std::vector<double> p, wt, es, th;
      for (const ThmSpectatorWindow::Node &k : nodes) {
        p.push_back(k.p);
        wt.push_back(k.weight);
        es.push_back(k.es);
        th.push_back(k.theta);
      }
      out.p.push_back(p);
      out.weight.push_back(wt);
      out.es.push_back(es);
      out.theta.push_back(th);
    }
  } else if (dw) {
    out.window = dw->angles ? dw->dist.AngleText() : "delta";
    if (dw->angles) out.muSx = dw->dist.muSx;
    for (double e : energies) out.reached.push_back(dw->Reached(e) ? 1 : 0);
  } else {
    out.window = "delta";
    for (size_t i = 0; i < energies.size(); i++) {
      out.p.push_back({0.0});
      out.weight.push_back({1.0});
      out.es.push_back({configure.thm.SpectatorEnergy(group->pairKey)});
      out.reached.push_back(1);
    }
  }
  const double mu = pair->GetRedMass() * uconv;
  // M_l = (B - 1) j_l - rho j_l' + C_l at E with T_s = es added to E + B (EPoint::CalcEDependentValues).
  struct Pieces {
    double jl = 0.0, rhoDjl = 0.0;
    complex coul = complex(0.0, 0.0);
  };
  auto pieces = [&](int l, double e, double es) {
    Pieces q;
    double b = out.bind + es;
    if (e + b > 0.0) {
      ThmBesselParts(l, mu, e, b, out.radius, q.jl, q.rhoDjl);
      if (configure.thm.coulombIntegral && pair->GetZ(1) * pair->GetZ(2) != 0)
        q.coul = ThmCoulombTerm(pair, l, e, ThmRho(mu, e, b, 1.0), useGSL);
    }
    return q;
  };
  // vertexModel=dw: the Gram matrices of the DW vertex on the grid.
  std::vector<std::vector<double>> dwG, dwGd;
  if (dw) {
    out.model = "dw";
    for (double e : energies) {
      std::vector<double> w, q, g, gd;
      double qd, pd;
      dw->Interpolate(e, w, q, g, gd, qd, pd);
      if (dw->angles) {
        // The directions (at the nearest grid energy); M2 uses the averaged G (one node).
        std::vector<double> aw, aq, th;
        dw->AngleNodesAt(e, aw, aq, th);
        out.dwTheta.push_back(th);
        out.dwAngleQ.push_back(aq);
        out.dwAngleWeight.push_back(aw);
      }
      const double ks = std::sqrt(2.0 * dw->dist.sf.mu * std::max(dw->dist.EsF(e), 0.0)) / hbarc;
      const double ka = dw->dist.aa.k, kb = dw->beta * ka;
      std::vector<double> row;
      for (double qk : q) {
        double qq = qk / hbarc, x = ks > 0.0 ? (ks * ks + kb * kb - qq * qq) / (2.0 * ks * kb) : 1.0;
        x = std::max(-1.0, std::min(1.0, x));
        row.push_back(std::sqrt(std::max(0.0, ka * ka + dw->alpha * dw->alpha * ks * ks - 2.0 * ka * dw->alpha * ks * x)) *
                      out.radius);
      }
      out.rho.push_back(row);
      out.dwQ.push_back(q);
      out.dwWeight.push_back(w);
      out.dwQDelta.push_back(qd * hbarc);
      out.dwPDelta.push_back(pd);
      dwG.push_back(g);
      dwGd.push_back(gd);
    }
  } else
    for (double e : energies) {
      std::vector<double> row;
      for (double es : out.es[out.rho.size()])
        row.push_back(e + out.bind + es > 0.0 ? ThmRho(mu, e, out.bind + es, out.radius) : 0.0);
      out.rho.push_back(row);
    }
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      ThmVertexReport::Channel cr;
      cr.jgroup = j;
      cr.channel = ch;
      cr.J = jg->GetJ();
      cr.pi = jg->GetPi();
      cr.l = c->GetL();
      cr.s = c->GetS();
      // The vertex boundary of the model (THMMatrixFunc::CalculateTHMCrossSection).
      const ThmVertexBoundary boundaryOf(configure, pair, jg, ch);
      // The pieces on the grid: [E][node] and the quasi-free [E].
      std::vector<std::vector<Pieces>> grid(energies.size());
      std::vector<Pieces> qf(energies.size());
      if (!dw)
        for (size_t i = 0; i < energies.size(); i++) {
          for (double es : out.es[i]) grid[i].push_back(pieces(cr.l, energies[i], es));
          qf[i] = pieces(cr.l, energies[i], 0.0);
        }
      const int li = dw ? dw->LIndex(cr.l) : -1;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmVertexReport::Level lr;
        lr.level = la;
        lr.boundary = boundaryOf.OnShell() ? std::nan("") : boundaryOf.Level(level);
        for (size_t i = 0; i < energies.size(); i++) {
          const complex B = boundaryOf.OnShell() ? boundaryOf.OnShellAt(energies[i]) : complex(lr.boundary, 0.0);
          if (dw) {
            const int nl = (int)dw->lvals.size();
            double avg = 0.0;
            for (size_t k = 0; k < out.dwWeight[i].size(); k++)
              avg += out.dwWeight[i][k] * (li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwG[i][(k * nl + li) * 4], B));
            lr.m2.push_back(avg);
            lr.m2qf.push_back(li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwGd[i][li * 4], B));
            // Plane waves at the same kinematics: M_l(p), analytic j_l'.
            double rho = out.dwPDelta[i] * out.radius;
            double jl = gsl_sf_bessel_jl(cr.l, rho);
            double djl = rho > 0.0 ? cr.l * jl / rho - gsl_sf_bessel_jl(cr.l + 1, rho) : (cr.l == 1 ? 1.0 / 3.0 : 0.0);
            lr.m2pw.push_back(std::norm((B - 1.0) * jl - rho * djl));
            continue;
          }
          auto m2 = [&](const Pieces &q) { return std::norm((B - 1.0) * q.jl - q.rhoDjl + q.coul); };
          double avg = 0.0;
          for (size_t k = 0; k < out.es[i].size(); k++) avg += out.weight[i][k] * m2(grid[i][k]);
          lr.m2.push_back(avg);
          lr.m2qf.push_back(m2(qf[i]));
        }
        cr.levels.push_back(lr);
      }
      out.channels.push_back(cr);
    }
  }
  return true;
}
