/*
 * The THM experiments of an WorkingCopy (<thm> experiment[<name>] ...; ThmExperiment.h):
 * building the groups and their physics objects, the profiled scale and
 * background, the coherent-background parameters and the per-experiment
 * profile reports (ThmExperimentReports).  thm_experiments.out and the pyazr
 * thm_* tables are written in ThmReports.cpp.  Members of WorkingCopy, kept out
 * of WorkingCopy.cpp so that the classic data handling stays as it is on dev.
 */
#include "CompoundSystem.h"
#include "Constants.h"
#include "PPair.h"
#include "ThmFunc.h"
#include "ThmVertexBoundary.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "WorkingCopy.h"
#include "ModelPass.h"
#include "Minuit2/MnUserParameters.h"
#include "GSLException.h"
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

void WorkingCopy::FillThmCoherentFromParams(const vector_r &p) {
  if (set_->thmCoherentParams_.empty() || set_->thmCoherentParamOffset_ < 0) return;
  const size_t n = set_->thmCoherentParams_.size();
  if (p.size() < (size_t)set_->thmCoherentParamOffset_ + n) return;
  for (size_t k = 0; k < n; k++) thmCoherentValues_[k] = p[set_->thmCoherentParamOffset_ + k];
}

// ---------------------------------------------------------------------------
// THM experiments (<thm> experiment[<name>] ...): see ThmExperiment.h.

// The lowest and highest c.m. energy of the points of a group's segments.
void WorkingCopy::ThmGroupEnergyRange(WorkingCopy &data, const WorkingCopy::ThmGroup &group, double &eLo, double &eHi) {
  eLo = 1.0e300;
  eHi = -1.0e300;
  for (int s : group.segments)
    for (int p = 1; p <= data.Data().GetSegment(s)->NumPoints(); p++) {
      eLo = std::min(eLo, data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
      eHi = std::max(eHi, data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
    }
}

namespace {

const double kAmu = uconv;  // MeV/u

// What the kinematics of an experiment line give (ThmResolveKinematics).
struct ThmReaction {
  ThmDistortion::Kinematics dk;  ///< for distortion=coulomb|optical and the DW vertex
  const PPair *pair = nullptr;   ///< the entrance pair x + A
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
bool ThmResolveKinematics(WorkingCopy &data, const Config &configure, const CompoundSystem *compound, const ThmExperiment &x,
                          const WorkingCopy::ThmGroup &group, const std::string &where, ThmReaction &r) {
  ThmDistortion::Kinematics &dk = r.dk;
  int pairKey = data.Data().GetSegment(group.segments[0])->GetEntranceKey();
  for (int s : group.segments)
    if (data.Data().GetSegment(s)->GetEntranceKey() != pairKey) {
      configure.outStream << where << "beam/target/spectator describe one reaction, but its segments "
                             "have different entrance pairs; give its data sets one entrance pair, or split them "
                             "into experiments of their own."
                          << std::endl;
      return false;
    }
  const PPair *pair = compound->GetPair(compound->GetPairNumFromKey(pairKey));
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
bool ThmBuildWindow(WorkingCopy &data, const Config &configure, const ThmExperiment &x, WorkingCopy::ThmGroup &group,
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
    for (int p = 1; p <= data.Data().GetSegment(s)->NumPoints(); p++)
      dataE.push_back(data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
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
  for (int s : group.segments) data.Draft().GetSegment(s)->SetThmSpectatorWindow(data, window);
  return true;
}

// Coulomb line shape of the spectator (ThmLineshape.h), lineshape=on.  The
// level energies and widths it uses are the observed ones of Brune.
bool ThmBuildLineshape(WorkingCopy &data, const Config &configure, const CompoundSystem *compound, const ThmExperiment &x,
                       WorkingCopy::ThmGroup &group, const std::string &where, const ThmReaction &r) {
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    configure.outStream << where << "lineshape=on uses the observed level energies and widths "
                           "as the resonance poles; it needs the Brune parameterization: run with it (the "
                           "<config> option use_brune), or set lineshape=off."
                        << std::endl;
    return false;
  }
  const ThmNuclide &b = x.beam, &t = x.target, &sp = x.spectator;
  const PPair *pair = r.pair;
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
    for (int p = 1; p <= data.Data().GetSegment(s)->NumPoints(); p++)
      eMax = std::max(eMax, data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
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
    int key = data.Data().GetSegment(s)->GetExitKey();
    bool seen = false;
    for (const ThmLineshape::Exit &e : shape->exits) seen = seen || e.pairKey == key;
    if (seen || !compound->IsPairKey(key)) continue;
    const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(key));
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
  for (int s : group.segments) data.Draft().GetSegment(s)->SetThmLineshape(data, shape);
  return true;
}

// Distorted-wave entrance vertex (ThmDwVertex.h): replaces M_l in the HOES
// amplitude of every segment of the experiment; R(E) is not applied.
bool ThmBuildDwVertex(WorkingCopy &data, const Config &configure, const CompoundSystem *compound, const ThmExperiment &x,
                      WorkingCopy::ThmGroup &group, const std::string &where, const ThmDistortion::Kinematics &dk) {
  int pairKey = data.Data().GetSegment(group.segments[0])->GetEntranceKey();
  int pairNum = compound->GetPairNumFromKey(pairKey);
  const PPair *pair = compound->GetPair(pairNum);
  // coulombIntegral=1 is refused with it by Config::ReadThmBlock
  // (CheckThmCoulombConsistency).
  const std::string refused =
      CheckThmDwVertexOptions(configure.thm.coherentL, configure.thm.SpectatorEnergy(pairKey), pairKey);
  if (!refused.empty()) {
    configure.outStream << where << refused << std::endl;
    return false;
  }
  std::vector<int> ls;
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int ch = 1; ch <= jg->NumChannels(); ch++)
      if (jg->GetChannel(ch)->GetPairNum() == pairNum) ls.push_back(jg->GetChannel(ch)->GetL());
  }
  std::vector<double> energies;
  for (int s : group.segments)
    for (int p = 1; p <= data.Data().GetSegment(s)->NumPoints(); p++)
      energies.push_back(data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
  double eLo, eHi;
  WorkingCopy::ThmGroupEnergyRange(data, group, eLo, eHi);
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
                               : v->Build(x, dk, pair->GetPairRadius(), ls, eLo - 0.3, gridHi, energies);
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
    if (data.Data().GetSegment(s)->GetThmWeight())
      l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << data.Data().GetSegment(s)->GetSegmentKey()
        << " also has weight[" << data.Data().GetSegment(s)->GetSegmentKey()
        << "]=; it multiplies the model with the DW vertex.";
  configure.outStream << l.str() << std::endl;
  group.dwVertex = v;
  for (int s : group.segments) data.Draft().GetSegment(s)->SetThmDwVertex(data, v);
  return true;
}

// Distortion factor R(E) (ThmDistortion.h): multiplies the model of every
// segment of the experiment before the folding.
bool ThmBuildDistortion(WorkingCopy &data, const Config &configure, const ThmExperiment &x, WorkingCopy::ThmGroup &group,
                        const std::string &where, const ThmDistortion::Kinematics &dk) {
  double eLo, eHi;
  WorkingCopy::ThmGroupEnergyRange(data, group, eLo, eHi);
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
                          << table.e.back() << "] MeV; give a table that covers the points, or drop the points "
                             "outside it." << std::endl;
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
      for (int p = 1; p <= data.Data().GetSegment(s)->NumPoints(); p++) {
        std::string why = d->CheckEnergy(data.Data().GetSegment(s)->GetPoint(data, p)->GetCMEnergy());
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
      for (int p = 1; why.empty() && p <= data.Data().GetSegment(group.segments[s])->NumPoints(); p++)
        why = d->CheckWindow(data.Data().GetSegment(group.segments[s])->GetPoint(data, p)->GetCMEnergy());
    if (!why.empty()) {
      configure.outStream << where << "distortion: " << why << "." << std::endl;
      return false;
    }
    ThmDistortion::Point lo = d->Evaluate(eLo), hi = d->Evaluate(eHi);
    l << "  Distortion factor R(E), zero-range DWBA: " << d->description << ".\n"
      << "  k_aA = " << d->aa.k << " fm^-1, eta_aA = " << d->aa.eta << ", kappa = " << d->kappa
      << " fm^-1, eta_b = " << d->etaB << ", beta = m_s/m_a = " << d->beta << "; E_ref = " << d->eRef
      << " MeV, grid " << d->gridLo << " to " << d->gridLo + (d->lnR.size() - 1) * d->gridStep
      << " MeV, radial step " << d->h << " fm to " << (d->n - 1) * d->h << " fm, l <= " << d->entrance.u.size() - 1
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
    if (data.Data().GetSegment(s)->GetThmWeight())
      l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << data.Data().GetSegment(s)->GetSegmentKey()
        << " also has weight[" << data.Data().GetSegment(s)->GetSegmentKey()
        << "]=; both multiply its model (the distortion factor and the weight).";
  configure.outStream << l.str() << std::endl;
  group.distortion = d;
  for (int s : group.segments) data.Draft().GetSegment(s)->SetThmDistortion(data, d);
  return true;
}

// Angular window of the exit pair (ThmAngular.h), theta=: the model of every
// segment is dsigma/dOmega averaged over theta_cm in the window.
bool ThmBuildAngleWindow(WorkingCopy &data, const Config &configure, const CompoundSystem *compound, const ThmExperiment &x,
                         WorkingCopy::ThmGroup &group, const std::string &where) {
  if (configure.thm.coherentL) {
    configure.outStream << where << "theta= computes the interference of the entrance partial waves "
                           "exactly (at fixed angle they interfere); entranceL=coherent is an approximation "
                           "of the angle-integrated observable and cannot be combined with it."
                        << std::endl;
    return false;
  }
  int maxLp = 0;
  for (int j = 1; j <= compound->NumJGroups(); j++)
    for (int ch = 1; ch <= compound->GetJGroup(j)->NumChannels(); ch++)
      maxLp = std::max(maxLp, compound->GetJGroup(j)->GetChannel(ch)->GetL());
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
  for (int s : group.segments) data.Draft().GetSegment(s)->SetThmAngleWindow(data, w);
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
int WorkingCopy::BuildThmGroups(const Config &configure, const CompoundSystem *compound) {
  Draft().thmGroups_.clear();
  Draft().thmCoherentParams_.clear();
  Draft().thmCoherentParamOffset_ = -1;
  for (const ThmExperiment &x : configure.thm.experiments) {
    const std::string where = "ERROR: <thm> experiment[" + x.name + "]: ";
    ThmGroup group;
    group.name = x.name;
    group.terms = x.backgroundTerms;
    for (int key : x.segments) {
      int index = 0;
      for (int s = 1; s <= Draft().NumSegments(); s++)
        if (Draft().GetSegment(s)->GetSegmentKey() == key) index = s;
      if (index == 0) {
        configure.outStream << "WARNING: <thm> experiment[" << x.name << "]: the data set with key " << key
                            << " is not in use; it is left out." << std::endl;
        continue;
      }
      const Segment *segment = Draft().GetSegment(index);
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
      for (int p = 1; p <= Draft().GetSegment(s)->NumPoints(); p++)
        if (Draft().GetSegment(s)->LabError(p - 1) != 0.0) points++;
    if (points <= 1 + group.terms) {
      configure.outStream << where << points << " point(s) with an error for " << 1 + group.terms
                          << " profiled linear parameter(s) (norm and background " << ThmExperiment::BackgroundName(group.terms)
                          << "); it needs more." << std::endl;
      return -1;
    }
    group.trivial = group.segments.size() == 1 && group.terms == 0;
    group.pairKey = Draft().GetSegment(group.segments[0])->GetEntranceKey();
    for (int s : group.segments)
      if (Draft().GetSegment(s)->GetEntranceKey() != group.pairKey) group.pairKey = 0;

    std::ostringstream summary;
    summary << "THM experiment '" << x.name << "': segment" << (group.segments.size() > 1 ? "s " : " ");
    for (size_t k = 0; k < group.segments.size(); k++)
      summary << (k ? "," : "") << Draft().GetSegment(group.segments[k])->GetSegmentKey();
    summary << (group.segments.size() > 1 ? " share one profiled norm" : " with a profiled norm")
            << ", background " << ThmExperiment::BackgroundName(group.terms) << ".";
    configure.outStream << summary.str() << std::endl;

    ThmReaction reaction;  // its kinematics are those of distortion=coulomb|optical and the DW vertex
    if (x.hasKinematics) {
      if (!ThmResolveKinematics(*this, configure, compound, x, group, where, reaction)) return -1;
      if (x.psKind != ThmExperiment::PS_DELTA && !ThmBuildWindow(*this, configure, x, group, where, reaction))
        return -1;
      if (x.lineshape && !ThmBuildLineshape(*this, configure, compound, x, group, where, reaction)) return -1;
    }
    if (x.vertexDW) {
      if (!ThmBuildDwVertex(*this, configure, compound, x, group, where, reaction.dk)) return -1;
    } else if (x.distortion != ThmExperiment::DIST_NONE) {
      if (!ThmBuildDistortion(*this, configure, x, group, where, reaction.dk)) return -1;
    }
    if (x.hasTheta && !ThmBuildAngleWindow(*this, configure, compound, x, group, where)) return -1;
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
      const int aa = compound->GetPairNumFromKey(group.pairKey);
      std::shared_ptr<ThmCoherentBackground> cb = std::make_shared<ThmCoherentBackground>();
      cb->experiment = x.name;
      std::ostringstream c;
      c << "  Coherent background (cbackground=): ";
      for (const ThmExperiment::CoherentTerm &t : x.cbackground) {
        const std::string jpi = ThmSpinText(t.J) + (t.parity > 0 ? "+" : "-");
        const std::string term = "cbackground " + jpi + ":" + std::to_string(t.exitKey) + ": ";
        bool exitUsed = false;
        for (int s : group.segments) exitUsed = exitUsed || Draft().GetSegment(s)->GetExitKey() == t.exitKey;
        if (!exitUsed || !compound->IsPairKey(t.exitKey)) {
          configure.outStream << where << term << "no segment of the experiment has exit pair " << t.exitKey << "."
                              << std::endl;
          return -1;
        }
        const int exitNum = compound->GetPairNumFromKey(t.exitKey);
        int j = 0;
        for (int jj = 1; jj <= compound->NumJGroups(); jj++) {
          const JGroup *jg = compound->GetJGroup(jj);
          if (jg->IsInRMatrix() && std::fabs(jg->GetJ() - t.J) < 1.0e-6 && jg->GetPi() == t.parity) j = jj;
        }
        if (j == 0) {
          configure.outStream << where << term << "the model has no J^pi = " << jpi
                              << " group with a level in the R matrix." << std::endl;
          return -1;
        }
        const JGroup *jg = compound->GetJGroup(j);
        int found = 0;
        for (int in = 1; in <= jg->NumChannels(); in++) {
          const AChannel *ci = jg->GetChannel(in);
          if (ci->GetPairNum() != aa) continue;
          if (t.hasChannels && (std::fabs(ci->GetS() - t.s) > 1.0e-6 || ci->GetL() != t.l)) continue;
          for (int out = 1; out <= jg->NumChannels(); out++) {
            const AChannel *co = jg->GetChannel(out);
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
            combo.index = (int)Draft().thmCoherentParams_.size();
            cb->combos.push_back(combo);
            const std::string stem = ThmCoherentParamStem(x.name, t.J, t.parity, t.exitKey, ci->GetS(), ci->GetL(),
                                                          co->GetS(), co->GetL());
            static const char *const parts[4] = {"re0", "im0", "re1", "im1"};
            for (int k = 0; k < 2 * t.form; k++) {
              ThmCoherentParam par;
              par.name = stem + parts[k];
              par.experiment = x.name;
              par.value = t.hasValues ? t.value[k] : 0.0;
              par.fixed = t.hasValues && t.fixed[k];
              par.jGroup = j;
              par.entrance = in;
              par.exit = out;
              par.term = k;
              Draft().thmCoherentParams_.push_back(par);
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
      for (int s : group.segments) Draft().GetSegment(s)->SetThmCoherent(*this, cb);
    }
    Draft().thmGroups_.push_back(group);
  }
  // The coherent backgrounds are the last block of the parameter vector
  // (FillMnParams sets the offset again when it lays the vector out).
  if (!Draft().thmCoherentParams_.empty()) {
    int offset = compound->NumParameters();
    for (int s = 1; s <= Draft().NumSegments(); s++) {
      if (Draft().GetSegment(s)->IsVaryNorm() && !Draft().GetSegment(s)->IsProfiledNorm()) offset++;
      offset++;  // energy shift
    }
    Draft().thmCoherentParamOffset_ = offset;
    std::vector<double> v;
    for (const ThmCoherentParam &par : Draft().thmCoherentParams_) v.push_back(par.value);
    thmCoherentValues_ = v;
  }
  return 0;
}

namespace {
// The points of the profiled segments of a group, concatenated in order.
void GatherThmGroup(const WorkingCopy *data, const azure2::Results &results, const std::vector<int> &segments, std::vector<double> &m,
                    std::vector<double> &d, std::vector<double> &e, std::vector<double> &energy) {
  m.clear();
  d.clear();
  e.clear();
  energy.clear();
  for (int s : segments) {
    const Segment *seg = data->Data().GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      const Unit *pt = seg->GetPoint(*data, p);
      if (!pt) continue;
      m.push_back(azure2::ModelInDataFrame(*data, *seg, p - 1, pt->GetModelValue(results)));
      d.push_back(seg->LabValue(p - 1));
      e.push_back(seg->LabError(p - 1));
      energy.push_back(pt->GetCMEnergy());
    }
  }
}
}  // namespace

double WorkingCopy::ProfileThmGroup(int g, azure2::Results &results) {
  const ThmGroup &group = set_->thmGroups_[g];
  std::vector<double> m, d, e, energy;
  GatherThmGroup(this, results, group.segments, m, d, e, energy);
  azure2::Results::ObjectiveValues &out = results.objective;
  const ThmProfile &p = out.thmProfiles[g] = SolveThmProfile(m, d, e, energy, group.terms);
  double total = 0.0;
  size_t k = 0;
  for (int s : group.segments) {
    const Segment *seg = Data().GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    out.norm[s] = p.Norm();
    double chi = 0.0;
    for (int q = 1; q <= seg->NumPoints(); q++) {
      if (!seg->GetPoint(*this, q)) continue;
      double r = p.Residual(m[k], d[k], e[k], energy[k]);
      chi += r * r;
      k++;
    }
    out.segmentChi2[s] = chi;
    total += chi;
  }
  return total;
}

void WorkingCopy::ThmGroupResiduals(int g, int i, const azure2::Results &results, std::vector<double> &out) {
  out.clear();
  const ThmProfile &p = results.objective.thmProfiles[g];
  const Segment *seg = Data().GetSegment(i);
  for (int q = 1; q <= seg->NumPoints(); q++) {
    Unit *pt = seg->GetPoint(*this, q);
    if (!pt) continue;
    out.push_back(p.Residual(azure2::ModelInDataFrame(*this, *seg, q - 1, pt->GetModelValue(results)),
                             seg->LabValue(q - 1), seg->LabError(q - 1), pt->GetCMEnergy()));
  }
}

double WorkingCopy::ThmBackgroundAt(int i, double energy, const azure2::Results &results) {
  int g = Data().ThmGroupOf(i);
  if (g < 0 || set_->thmGroups_[g].terms == 0) return 0.0;
  return results.objective.thmProfiles[g].Background(energy);
}

std::vector<ThmExperimentReport> WorkingCopy::ThmExperimentReports(const azure2::Results &results) const {
  std::vector<ThmExperimentReport> out;
  for (size_t g = 0; g < set_->thmGroups_.size(); g++) {
    const ThmGroup &group = set_->thmGroups_[g];
    ThmExperimentReport r;
    r.name = group.name;
    r.background = ThmExperiment::BackgroundName(group.terms);
    ThmProfile p;
    if (group.trivial || Data().ThmGroupOf(group.segments[0]) < 0) {
      // Profiled per segment (Segment::ProfileNormChiSquared) or not at all:
      // the same closed form, recomputed here for its uncertainty.
      std::vector<double> m, d, e, energy;
      GatherThmGroup(this, results, group.segments, m, d, e, energy);
      p = SolveThmProfile(m, d, e, energy, group.terms);
    } else {
      p = results.objective.thmProfiles[g];
    }
    r.points = 0;
    r.chi2 = 0.0;
    for (int s : group.segments) {
      const Segment *seg = Data().GetSegment(s);
      if (!seg->IsProfiledNorm()) continue;
      r.segments.push_back(seg->GetSegmentKey());
      r.chi2 += results.objective.segmentChi2[s];
    }
    r.points = p.points;
    p.Reported(r.value, r.covariance);
    r.status = r.segments.empty() ? "not profiled: no segment of it has a free norm" : p.status;
    for (int k = 0; k < Data().NumThmCoherentParams(); k++)
      if (const ThmCoherentParam &c = Data().GetThmCoherentParam(k); c.experiment == group.name) {
        r.coherentNames.push_back(c.name);
        r.coherentValues.push_back(thmCoherentValues_[k]);
        r.coherentFixed.push_back(c.fixed);
      }
    out.push_back(r);
  }
  return out;
}
