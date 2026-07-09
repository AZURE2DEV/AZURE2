#include "THMMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include "Constants.h"
#include <cmath>
#include <map>
#include <vector>

THMMatrixFunc::THMMatrixFunc(CNuc* compound, const Config& configure)
  : AMatrixFunc(compound, configure) {}

/*!
 * HOES cross section of the modified R-matrix formalism (arbitrary units).
 * Assumes ClearMatrices/FillMatrices/InvertMatrices have already produced the
 * level matrix A (shared interior). Uses the entrance transfer form factors
 * M_l and exit penetrabilities stored on the point by CalcEDependentValues.
 */
void THMMatrixFunc::CalculateTHMCrossSection(EPoint* point) {
  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  PPair* entrancePair = compound()->GetPair(aa);
  PPair* exitPair = compound()->GetPair(exitPairNum);

  // Exit-channel flux factor k_f / mu_f (mrmpy). The entrance c.m. energy sets
  // the compound-system energy; the exit channel energy is that minus the exit
  // pair threshold. Below the exit threshold there is no outgoing flux.
  double inEnergy = point->GetCMEnergy() + entrancePair->GetSepE()
                    + entrancePair->GetExE();
  double exitEnergy = inEnergy - exitPair->GetSepE() - exitPair->GetExE();
  double muf = exitPair->GetRedMass() * uconv;                 // MeV/c^2
  double kf = (exitEnergy > 0.0) ? std::sqrt(2.0 * muf * exitEnergy) / hbarc
                                 : 0.0;                        // fm^-1
  double fluxFactor = kf / muf;

  double sigma = 0.0;
  for(int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup* jg = compound()->GetJGroup(j);
    if(!jg->IsInRMatrix()) continue;
    int numLevels = jg->NumLevels();
    int numChannels = jg->NumChannels();
    double spinWeight = 2.0 * jg->GetJ() + 1.0;

    // Active-level index map (a_matrices_ / GetAMatrixElement are indexed by the
    // compacted in-R-matrix levels, so a level not in the R matrix is skipped).
    std::vector<int> act(numLevels + 1, 0);
    int na = 0;
    for(int la = 1; la <= numLevels; la++)
      if(jg->GetLevel(la)->IsInRMatrix()) act[la] = ++na;
    if(na == 0) continue;

    // Entrance vertices: v_s[la] = sum_{c_in in spin s} gamma_{la,c_in} M_l,
    // coherent over entrance partial waves of the same channel spin, kept in
    // separate (incoherent) buckets per channel spin s.
    std::map<double, std::vector<complex> > vbys;
    bool hasEntrance = false;
    for(int ch = 1; ch <= numChannels; ch++) {
      AChannel* c = jg->GetChannel(ch);
      if(c->GetPairNum() != aa) continue;
      hasEntrance = true;
      double ml = point->GetThmFormFactor(j, ch);
      std::vector<complex>& vertex = vbys[c->GetS()];
      if(vertex.empty()) vertex.assign(numLevels + 1, complex(0.0, 0.0));
      for(int la = 1; la <= numLevels; la++) {
        if(!jg->GetLevel(la)->IsInRMatrix()) continue;
        vertex[la] += jg->GetLevel(la)->GetFitGamma(ch) * ml;
      }
    }
    if(!hasEntrance) continue;   // this J group does not couple the entrance pair

    // Exit channels (incoherent), each with sqrt(2 P) folded in as 2 P outside.
    for(int ch = 1; ch <= numChannels; ch++) {
      AChannel* c = jg->GetChannel(ch);
      if(c->GetPairNum() != exitPairNum) continue;
      double sqrtPen = point->GetSqrtPenetrability(j, ch);
      double pex = sqrtPen * sqrtPen;               // P_l(E_exit); 0 if closed
      if(pex == 0.0) continue;

      double term = 0.0;
      for(std::map<double, std::vector<complex> >::iterator it = vbys.begin();
          it != vbys.end(); ++it) {
        std::vector<complex>& vertex = it->second;
        complex amp(0.0, 0.0);
        for(int la = 1; la <= numLevels; la++) {
          if(!jg->GetLevel(la)->IsInRMatrix()) continue;
          double gEx = jg->GetLevel(la)->GetFitGamma(ch);
          if(std::fabs(gEx) < 1.0e-12) continue;
          for(int lap = 1; lap <= numLevels; lap++) {
            if(!jg->GetLevel(lap)->IsInRMatrix()) continue;
            amp += gEx * this->GetAMatrixElement(j, act[la], act[lap])
                       * vertex[lap];
          }
        }
        term += std::norm(amp);                     // |amp|^2, incoherent over s
      }
      sigma += spinWeight * fluxFactor * 2.0 * pex * term;
    }
  }

  point->SetFitCrossSection(sigma);
}
