#include "THMMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Constants.h"
#include "CoulFunc.h"
#include "ShftFunc.h"
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
 */
static double ShiftAtLevelEnergy(PPair *pair, int l, double levelEnergy, bool useGSL) {
  struct Entry { PPair *pair; int l; double e; double s; };
  static const int kMemo = 32;
  thread_local Entry memo[kMemo];
  thread_local int filled = 0, next = 0;
  for (int i = 0; i < filled; i++)
    if (memo[i].pair == pair && memo[i].l == l && memo[i].e == levelEnergy) return memo[i].s;
  double resonanceEnergy = levelEnergy - (pair->GetSepE() + pair->GetExE());
  double s;
  if (resonanceEnergy < 0.0) {
    ShftFunc shift(pair);
    s = shift(l, levelEnergy);
  } else {
    CoulFunc coul(pair, useGSL);
    s = coul.PEShift(l, pair->GetChRad(), resonanceEnergy);
  }
  memo[next] = Entry{pair, l, levelEnergy, s};
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
    std::map<std::pair<double, int>, std::vector<complex>> vbys;
    bool hasEntrance = false;
    for (int ch = 1; ch <= numChannels; ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      hasEntrance = true;
      std::vector<complex> &vertex = vbys[std::make_pair(c->GetS(), coherentL ? 0 : c->GetL())];
      if (vertex.empty()) vertex.assign(numLevels + 1, complex(0.0, 0.0));
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
      for (int la = 1; la <= numLevels; la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        complex boundary = onShell ? onShellL
                           : constantVertex ? constantB
                           : complex(perLevel ? level->GetShiftFunction(ch)
                                              : c->GetBoundaryCondition(), 0.0);
        vertex[la] += level->GetFitGamma(ch) * point->GetThmFormFactor(j, ch, boundary);
      }
    }
    if (!hasEntrance) continue;  // this J group does not couple the entrance pair

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
          for (int lap = 1; lap <= numLevels; lap++) {
            if (!jg->GetLevel(lap)->IsInRMatrix()) continue;
            amp += gEx * this->GetAMatrixElement(j, act[la], act[lap]) * vertex[lap];
          }
        }
        term += std::norm(amp);  // |amp|^2, incoherent over (s, l)
      }
      sigma += spinWeight * fluxFactor * 2.0 * pex * term;
    }
  }

  point->SetFitCrossSection(sigma);
}
