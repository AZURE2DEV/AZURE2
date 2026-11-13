// The pathway table of a CNuc (PathwayTable.h, docs/upgrade/cnuc_split.md section 11): its inputs taken from the
// loaded compound, its build in CNuc::Initialize, and its handle.
#include "CNucPathways.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "CNucSystem.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"

namespace {

//! 2j as an integer (CNuc.cpp's Twice).
int Twice(double j) { return (int)std::lround(2.0 * j); }

}  // namespace

std::vector<azure2::PathwayFinalState> ECFinalStates(CNuc &c) {
  std::vector<azure2::PathwayFinalState> out;
  for (int j = 1; j <= c.NumJGroups(); j++) {
    JGroup *theFinalJGroup = c.GetJGroup(j);
    for (int la = 1; la <= theFinalJGroup->NumLevels(); la++) {
      ALevel *theFinalLevel = theFinalJGroup->GetLevel(la);
      if (!theFinalLevel->IsECLevel()) continue;
      azure2::PathwayFinalState f;
      f.state.two_J = Twice(theFinalJGroup->GetJ());
      f.state.parity = theFinalJGroup->GetPi();
      f.state.photon_partition = theFinalLevel->GetECPairNum() - 1;
      f.state.multipoles = theFinalLevel->GetECMultMask();
      for (int chp = 1; chp <= theFinalJGroup->NumChannels(); chp++) {
        AChannel *finalChannel = theFinalJGroup->GetChannel(chp);
        f.state.channels.push_back(azure2::physics::FinalStateChannel{
            finalChannel->GetPairNum() - 1, c.GetPair(finalChannel->GetPairNum())->GetPType() == 0, finalChannel->GetL(),
            Twice(finalChannel->GetS())});
      }
      f.group = j;
      f.level = la;
      out.push_back(f);
    }
  }
  return out;
}


const azure2::PathwayTable &CNuc::Pathways() const {
  if (!pathways_) throw std::logic_error("CNuc::Pathways: the compound nucleus has not been initialized");
  return *pathways_;
}

void CNuc::RequestCaptureAy(int entrancePair, int exitPair) {
  const std::pair<int, int> request(entrancePair, exitPair);
  for (const std::pair<int, int> &r : captureAyRequests_)
    if (r == request) return;
  captureAyRequests_.push_back(request);
}

void CNuc::RequestEntrance(int pair) {
  for (int r : entranceRequests_)
    if (r == pair) return;
  entranceRequests_.push_back(pair);
}

void CNuc::CompareEntranceRequests(const char *where) const {
  std::set<int> flagged, requested(entranceRequests_.begin(), entranceRequests_.end());
  for (int aa = 1; aa <= this->NumPairs(); aa++)
    if (this->GetPair(aa)->IsEntrance()) flagged.insert(aa);
  if (flagged == requested) return;
  std::cerr << "ENTRANCE COMPARISON FAILED at " << where << ": flagged by the data set-up {";
  for (int aa : flagged) std::cerr << " " << aa;
  std::cerr << " }, requested from the specs {";
  for (int aa : requested) std::cerr << " " << aa;
  std::cerr << " }" << std::endl;
  std::abort();
}

/*!
 * Builds the pathway table from the System's structure, the external-capture final states and the configuration:
 * R-matrix capture, the largest Legendre order, the entrance pairs (PPair::IsEntrance, set while the data are read)
 * and the capture-A_y requests.  Run by Initialize after the input transformation; nothing writes the table
 * afterwards.
 */
void CNuc::BuildPathways(const Config &configure) {
  azure2::PathwayTableOptions options;
  options.rmatrix_capture = (configure.paramMask & Config::USE_RMC_FORMALISM) != 0;
  options.max_l_order = configure.maxLOrder;
  this->CompareEntranceRequests("CNuc::BuildPathways");
  for (int aa = 1; aa <= this->NumPairs(); aa++)
    if (this->GetPair(aa)->IsEntrance()) options.entrances.push_back(aa);
  options.capture_ay = captureAyRequests_;
  pathways_ = azure2::BuildPathwayTable(BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure),
                                        ECFinalStates(*this), options);
}
