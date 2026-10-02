#include "THMMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Constants.h"
#include "ChannelFunc.h"
#include "CoulFunc.h"
#include "ShftFunc.h"
#include "ThmLineshape.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmAngular.h"
#include "ThmExperiment.h"
#include "NuclearPotentialManager.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

THMMatrixFunc::THMMatrixFunc(CNuc *compound, const Config &configure) :
  AMatrixFunc(compound, configure) {}

/*!
 * Shift function S_c at a level energy (compound-system excitation, MeV), the
 * value CNuc::CalcShiftFunctions stores per level under Brune.  A per-thread
 * memo keeps it cheap: within one evaluation every point asks for the same
 * few (pair, l, energy) triples.
 *
 * The memo is keyed by what S_c depends on -- charges, reduced mass, channel
 * radius, the hybrid-potential state of the pair, l, the resonance energy and
 * the Coulomb routine -- not by the PPair address: the memo outlives a CNuc
 * (thread_local), a later CNuc in the same process (pyazr sessions one after
 * the other, a radius scan) can get the same address back from the allocator,
 * and a pair's radius can change in place.  Keyed by the address, a session at
 * 6.1 fm read S_c of an earlier 4.1 fm session (18O(p,a) narrow levels: chi2
 * 35.99 instead of 208.6 for the same parameters).
 */
static double ShiftAtLevelEnergy(PPair *pair, int l, double levelEnergy, bool useGSL) {
  struct Entry {
    int z1, z2; double redmass, radius; long hybridTag; int l; double e; bool gsl; double s;
  };
  static const int kMemo = 32;
  thread_local Entry memo[kMemo];
  thread_local int filled = 0, next = 0;
  const int z1 = pair->GetZ(1), z2 = pair->GetZ(2);
  const double redmass = pair->GetRedMass(), radius = pair->GetChRad();
  const bool hybrid = (g_config ? g_config->useHybridMethod : false) &&
      NuclearPotentialManager::instance().isPairEnabled(pair->GetPairKey());
  const long hybridTag = hybrid ? NuclearPotentialManager::instance().tagFor(pair->GetPairKey()) : 0;
  double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
  for (int i = 0; i < filled; i++) {
    const Entry &m = memo[i];
    if (m.z1 == z1 && m.z2 == z2 && m.redmass == redmass && m.radius == radius &&
        m.hybridTag == hybridTag && m.l == l && m.e == resonanceEnergy && m.gsl == useGSL)
      return m.s;
  }
  // Continuous through threshold, and at it (ChannelFunc).
  double s = ChannelFunc(pair, useGSL).Shift(l, resonanceEnergy);
  memo[next] = Entry{z1, z2, redmass, radius, hybridTag, l, resonanceEnergy, useGSL, s};
  next = (next + 1) % kMemo;
  if (filled < kMemo) filled++;
  return s;
}

/*!
 * HOES cross section of the modified R-matrix formalism (arbitrary units).
 * Assumes ClearMatrices/FillMatrices/InvertMatrices have already produced the
 * level matrix A (shared interior). Uses the entrance transfer form factors
 * M_l and exit penetrabilities stored on the point by CalcEDependentValues.
 */
void THMMatrixFunc::CalculateTHMCrossSection(EPoint *point) {
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  PPair *entrancePair = compound()->GetPair(aa);
  PPair *exitPair = compound()->GetPair(exitPairNum);

  // Kinematic factors, fixed by how the HOES data were divided out of the
  // triple cross section (Config::ThmOptions::kinematics).  The entrance c.m.
  // energy sets the compound-system energy; the exit channel energy is that
  // minus the exit pair threshold. Below the exit threshold there is no flux.
  double inEnergy = point->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
  double exitEnergy = inEnergy - exitPair->GetSepE() - exitPair->GetExE();
  double muf = exitPair->GetRedMass() * uconv;  // MeV/c^2
  double kf = (exitEnergy > 0.0) ? std::sqrt(2.0 * muf * exitEnergy) / hbarc
                                 : 0.0;  // fm^-1
  if (kf == 0.0) {
    point->SetFitCrossSection(0.0);
    return;
  }
  double fluxFactor = 1.0;
  switch (configure().thm.kinematics) {
    case Config::ThmOptions::LA_COGNATA:  // k_f/mu_f (mrmpy)
      fluxFactor = kf / muf;
      break;
    case Config::ThmOptions::TRIPLE:
      break;
    case Config::ThmOptions::KF_THREE_BODY:
      fluxFactor = 1.0 / (muf * kf);
      break;
    case Config::ThmOptions::LAMBDA32: {
      // On-shell entrance momentum, as Pizzone et al. take it.  Undefined
      // below the entrance threshold; |E| keeps sub-threshold tails finite.
      double mui = entrancePair->GetRedMass() * uconv;
      double ki = std::sqrt(2.0 * mui * std::max(std::fabs(point->GetCMEnergy()), 1.0e-6)) / hbarc;
      fluxFactor = 1.0 / ki;
      break;
    }
  }

  // Coulomb line shape of the spectator (<thm> experiment[...] lineshape=on,
  // ThmLineshape.h): the exit amplitude of level lambda -- the level that
  // decays to b + B -- carries N_C(zeta, E_lambda - E, Gamma_lambda), inside
  // the coherent level sum.  zeta depends on the energy and the exit pair only.
  const ThmLineshape *lineshape = point->GetThmLineshape();
  double zeta = 0.0;
  if (lineshape) {
    int light = exitPair->GetM(1) <= exitPair->GetM(2) ? 1 : 2;
    zeta = lineshape->Zeta(point->GetCMEnergy(), exitPair->GetZ(3 - light), exitPair->GetM(3 - light));
    if (!(lineshape->EsF(point->GetCMEnergy()) > 0.0) && !lineshape->warned.exchange(true))
      configure().outStream << "WARNING: <thm> experiment[" << lineshape->experiment << "] lineshape: at E = "
                            << point->GetCMEnergy() << " MeV the spectator has no energy left (E_sF <= 0); "
                            << "E_sF = 1 keV is used there (reported once)." << std::endl;
  }

  const ThmSpectatorWindow *window = point->GetThmSpectatorWindow();
  int numNodes = window ? point->NumThmPsNodes() : 0;

  // Distorted-wave entrance vertex (vertexModel=dw, ThmDwVertex.h): per node
  // (one: the spectatorAngle direction or the acceptance average) and entrance l, two
  // incoherent components M^(k)(B) = a_k (B - 1) - d_k replace M_l.
  const ThmDwVertex *dw = point->GetThmDwVertex();
  ThmDwVertex::At dwAt;
  if (dw) {
    dw->Evaluate(point->GetCMEnergy(), dwAt);
    if (dwAt.outside && !dw->warned.exchange(true))
      configure().outStream << "WARNING: <thm> experiment[" << dw->experiment
                            << "]: the DW vertex is evaluated at E = " << point->GetCMEnergy()
                            << " MeV, outside its grid; the end value is used there (reported once)." << std::endl;
    window = nullptr;
    numNodes = dwAt.nodes;
  }

  // Angular window of the exit pair (theta=thmin-thmax; ThmAngular.h): at a
  // fixed angle the J^pi groups and the entrance l of one channel spin
  // interfere, so the partial amplitudes of every J group are collected --
  // x = sqrt(K 2P) e^{i(omega - phi)} (amplitude of the exit channel), the exit
  // phase as in the T matrix of AMatrixFunc::CalculateTMatrix, with the
  // entrance vertex M_l in place of e^{i(omega-phi)} sqrt(P) -- and combined in
  // the Legendre sum after the loop, one set per spectator-window node.
  const ThmAngleWindow *angle = point->GetThmAngleWindow();
  std::vector<ThmPartialWave> waves;
  std::vector<std::vector<complex>> partial(angle ? std::max(numNodes, 1) : 0);

  // Coherent background of the experiment (cbackground=,
  // ThmCoherentBackground): c(E) times the bucket's vertex M_l (no width),
  // added to the resonant amplitude of its combinations before squaring.
  const ThmCoherentBackground *coherent = point->GetThmCoherent();
  const std::vector<double> &coherentValues = compound()->ThmCoherentValues();

  double sigma = 0.0;
  for (int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup *jg = compound()->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    int numLevels = jg->NumLevels();
    int numChannels = jg->NumChannels();
    double spinWeight = 2.0 * jg->GetJ() + 1.0;

    // Active-level index map (a_matrices_ / GetAMatrixElement are indexed by the
    // compacted in-R-matrix levels, so a level not in the R matrix is skipped).
    std::vector<int> act(numLevels + 1, 0);
    int na = 0;
    for (int la = 1; la <= numLevels; la++)
      if (jg->GetLevel(la)->IsInRMatrix()) act[la] = ++na;
    if (na == 0) continue;

    // Entrance vertices v[la] = sum_{c_in} gamma_{la,c_in} M_l, one bucket per
    // entrance channel spin s and orbital l, summed incoherently: with the
    // quantization axis along p_xA only m_l = 0 enters, and once the exit
    // direction is integrated and the spin projections summed,
    // sum_{M,m_s} <s m_s l 0|J M><s m_s l' 0|J M> = delta_{ll'} (2J+1)/(2l+1)
    // removes the l cross terms.  `entranceL=coherent` keeps one bucket per s,
    // as mrmpy does.
    //
    // Vertex boundary (Config::ThmOptions::vertex): the per-level shift
    // S_c(E_lambda) under Brune (mrmpy vertex_boundary="per_level"; Tumino et
    // al. 2021 eq. 51), else the channel boundary constant; `constant` uses the
    // channel constant in either formalism; `onshell` the log-derivative of the
    // outgoing wave L_c(E) = S_c(E) + i P_c(E) (Tribble et al. 2014 eq. 2.76),
    // recovered from L_o = L_c - B_c.
    bool perLevel = (configure().paramMask & Config::USE_BRUNE_FORMALISM) &&
                    configure().thm.vertex == Config::ThmOptions::PER_LEVEL;
    bool onShell = configure().thm.vertex == Config::ThmOptions::ON_SHELL;
    bool constantVertex = configure().thm.vertex == Config::ThmOptions::CONSTANT;
    bool coherentL = configure().thm.coherentL;
    bool hasEntrance = false;
    for (int ch = 1; ch <= numChannels; ch++)
      if (jg->GetChannel(ch)->GetPairNum() == aa) hasEntrance = true;
    if (!hasEntrance) continue;  // this J group does not couple the entrance pair

    // N_C of each level (lineshape=on); the pole is the observed energy and
    // total width of the level (Brune).
    std::vector<complex> nc;
    if (lineshape) {
      nc.assign(numLevels + 1, complex(1.0, 0.0));
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        nc[la] = ThmLineshapeFactor(zeta, level->GetFitE() - inEnergy,
                                    ThmLevelWidth(compound(), jg, level, configure()));
      }
    }

    // Spectator-momentum window (ps=..., ThmLineshape.h ThmSpectatorWindow):
    // the cross section -- not the amplitude -- is averaged over the accepted
    // directions at this energy (nodes p_k) with the normalized event weights
    // w_k (|phi|^2 d cos theta_cm), since each direction is a distinct final
    // state (Mukhamedzhanov et al. 2017 eq. 34).  Without a window one
    // pass with the stored single-node vertex (node -1), as before.
    bool coherentJ = false;
    if (coherent)
      for (const ThmCoherentBackground::Combo &cb : coherent->combos) coherentJ = coherentJ || cb.jGroup == j;

    const int passes = numNodes > 0 ? numNodes : 1;
    for (int pass = 0; pass < passes; pass++) {
      const int node = numNodes > 0 ? pass : -1;
      std::map<std::pair<double, int>, std::vector<complex>> vbys;
      // cbackground=: per bucket its entrance channel and the vertex without a width.
      std::map<std::pair<double, int>, std::pair<int, complex>> unitVertex;
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *c = jg->GetChannel(ch);
        if (c->GetPairNum() != aa) continue;
        // Buckets: (s, l), or (s, 0) for entranceL=coherent; with the DW vertex
        // (s, 2 l + k), k = 0, 1 its two components.
        std::vector<complex> *pwVertex = nullptr;
        if (!dw) {
          pwVertex = &vbys[std::make_pair(c->GetS(), coherentL ? 0 : c->GetL())];
          if (pwVertex->empty()) pwVertex->assign(numLevels + 1, complex(0.0, 0.0));
        }
        complex onShellL = point->GetLoElement(j, ch) + c->GetBoundaryCondition();
        // `constant`: S_c at the lowest level of the J group, whatever the order
        // of the levels in the file (the channel boundary constant of the
        // R-matrix is tied to the first level read, which is not physical).
        complex constantB(0.0, 0.0);
        if (constantVertex) {
          double eMin = 0.0;
          bool found = false;
          for (int la = 1; la <= numLevels; la++) {
            ALevel *level = jg->GetLevel(la);
            if (!level->IsInRMatrix()) continue;
            if (!found || level->GetFitE() < eMin) eMin = level->GetFitE();
            found = true;
          }
          constantB = ShiftAtLevelEnergy(compound()->GetPair(aa), c->GetL(), eMin,
                                         !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC));
        }
        if (dw) {
          // Two components per (s, l): buckets 2 l + k.
          const int li = dw->LIndex(c->GetL());
          std::vector<complex> &vertex1 = vbys[std::make_pair(c->GetS(), 2 * c->GetL() + 1)];
          if (vertex1.empty()) vertex1.assign(numLevels + 1, complex(0.0, 0.0));
          std::vector<complex> &vertex0 = vbys[std::make_pair(c->GetS(), 2 * c->GetL())];
          if (vertex0.empty()) vertex0.assign(numLevels + 1, complex(0.0, 0.0));
          if (li < 0) continue;
          const complex *ak = &dwAt.a[((size_t)pass * dwAt.nl + li) * 2];
          const complex *dk = &dwAt.d[((size_t)pass * dwAt.nl + li) * 2];
          if (coherentJ) {
            complex b = onShell ? onShellL : constantVertex ? constantB : complex(c->GetBoundaryCondition(), 0.0);
            unitVertex[std::make_pair(c->GetS(), 2 * c->GetL())] = std::make_pair(ch, ak[0] * (b - 1.0) - dk[0]);
            unitVertex[std::make_pair(c->GetS(), 2 * c->GetL() + 1)] = std::make_pair(ch, ak[1] * (b - 1.0) - dk[1]);
          }
          for (int la = 1; la <= numLevels; la++) {
            ALevel *level = jg->GetLevel(la);
            if (!level->IsInRMatrix()) continue;
            complex boundary = onShell ? onShellL
                               : constantVertex ? constantB
                               : complex(perLevel ? level->GetShiftFunction(ch)
                                                  : c->GetBoundaryCondition(), 0.0);
            double g = level->GetFitGamma(ch);
            vertex0[la] += g * (ak[0] * (boundary - 1.0) - dk[0]);
            vertex1[la] += g * (ak[1] * (boundary - 1.0) - dk[1]);
          }
          continue;
        }
        if (coherentJ) {
          complex b = onShell ? onShellL : constantVertex ? constantB : complex(c->GetBoundaryCondition(), 0.0);
          unitVertex[std::make_pair(c->GetS(), c->GetL())] =
              std::make_pair(ch, node < 0 ? point->GetThmFormFactor(j, ch, b) : point->GetThmFormFactor(j, ch, b, node));
        }
        for (int la = 1; la <= numLevels; la++) {
          ALevel *level = jg->GetLevel(la);
          if (!level->IsInRMatrix()) continue;
          complex boundary = onShell ? onShellL
                             : constantVertex ? constantB
                             : complex(perLevel ? level->GetShiftFunction(ch)
                                                : c->GetBoundaryCondition(), 0.0);
          (*pwVertex)[la] += level->GetFitGamma(ch) * (node < 0 ? point->GetThmFormFactor(j, ch, boundary)
                                                                : point->GetThmFormFactor(j, ch, boundary, node));
        }
      }

      // Exit channels (incoherent), each with sqrt(2 P) folded in as 2 P outside.
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *c = jg->GetChannel(ch);
        if (c->GetPairNum() != exitPairNum) continue;
        double sqrtPen = point->GetSqrtPenetrability(j, ch);
        double pex = sqrtPen * sqrtPen;  // P_l(E_exit); 0 if closed
        if (pex == 0.0) continue;

        double term = 0.0;
        for (std::map<std::pair<double, int>, std::vector<complex>>::iterator it = vbys.begin();
             it != vbys.end(); ++it) {
          std::vector<complex> &vertex = it->second;
          complex amp(0.0, 0.0);
          for (int la = 1; la <= numLevels; la++) {
            if (!jg->GetLevel(la)->IsInRMatrix()) continue;
            double gEx = jg->GetLevel(la)->GetFitGamma(ch);
            if (std::fabs(gEx) < 1.0e-12) continue;
            if (lineshape) {
              complex gExC = gEx * nc[la];
              for (int lap = 1; lap <= numLevels; lap++) {
                if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
                amp += gExC * this->GetAMatrixElement(j, act[la], act[lap]) * vertex[lap];
              }
              continue;
            }
            for (int lap = 1; lap <= numLevels; lap++) {
              if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
              amp += gEx * this->GetAMatrixElement(j, act[la], act[lap]) * vertex[lap];
            }
          }
          if (coherentJ) {
            std::map<std::pair<double, int>, std::pair<int, complex>>::const_iterator u = unitVertex.find(it->first);
            const ThmCoherentBackground::Combo *cb = u != unitVertex.end() ? coherent->Find(j, u->second.first, ch)
                                                                           : nullptr;
            if (cb && cb->index + 2 * cb->form <= (int)coherentValues.size()) {
              const double *v = &coherentValues[cb->index];
              complex cE(v[0], v[1]);
              if (cb->form == 2) cE += complex(v[2], v[3]) * point->GetCMEnergy();
              amp += cE * u->second.second;
            }
          }
          term += std::norm(amp);  // |amp|^2, incoherent over (s, l)
          if (angle) {
            if (pass == 0) {
              ThmPartialWave w;
              w.J = jg->GetJ();
              w.s = it->first.first;
              w.l = it->first.second;
              w.sp = c->GetS();
              w.lp = c->GetL();
              waves.push_back(w);
            }
            partial[pass].push_back(std::sqrt(fluxFactor * 2.0 * pex) * point->GetExpCoulombPhase(j, ch) *
                                    point->GetExpHardSpherePhase(j, ch) * amp);
          }
        }
        if (angle) continue;
        if (dw)
          sigma += dwAt.weight[pass] * (spinWeight * fluxFactor * 2.0 * pex * term);
        else if (node < 0)
          sigma += spinWeight * fluxFactor * 2.0 * pex * term;
        else
          sigma += point->GetThmPsWeight(node) * (spinWeight * fluxFactor * 2.0 * pex * term);
      }
    }
  }

  if (angle) {
    // sum_spins |F|^2 = (1/pi) sum_L b_L P_L(cos theta), averaged over the
    // window (and over the spectator-window nodes with their weights).
    std::vector<double> b;
    for (size_t k = 0; k < partial.size(); k++) {
      ThmLegendreCoefficients(waves, partial[k], b);
      sigma += (numNodes > 0 && !dw ? point->GetThmPsWeight((int)k) : numNodes > 0 ? dwAt.weight[k] : 1.0) *
               angle->Mean(b);
    }
  }

  // Energy-dependent weight of this segment (<thm> weight[k]=), e.g. the
  // Coulomb-distortion factor R(E) of Mukhamedzhanov & Pang, PRC 99 (2019)
  // 064618.  Applied to every point and folding sub-point, before the folding.
  if (const ThmWeightTable *weight = point->GetThmWeight()) {
    bool outside = false;
    sigma *= (*weight)(point->GetCMEnergy(), &outside);
    if (outside && !weight->warned.exchange(true))
      configure().outStream << "WARNING: <thm> weight '" << weight->name << "' is evaluated at E = "
                            << point->GetCMEnergy() << " MeV, outside its table ["
                            << weight->e.front() << ", " << weight->e.back()
                            << "] MeV; the end value is used there (reported once)." << std::endl;
  }
  // Distortion factor R(E) of the point's THM experiment (distortion=...;
  // ThmDistortion.h), the same convention: it multiplies the model before
  // the folding.
  if (const ThmDistortion *distortion = point->GetThmDistortion()) {
    bool outside = false;
    sigma *= distortion->Weight(point->GetCMEnergy(), &outside);
    if (outside && !distortion->warned.exchange(true))
      configure().outStream << "WARNING: <thm> experiment[" << distortion->experiment
                            << "]: the distortion factor is evaluated at E = " << point->GetCMEnergy()
                            << " MeV, outside its grid; the end value is used there (reported once)." << std::endl;
  }

  point->SetFitCrossSection(sigma);
}
