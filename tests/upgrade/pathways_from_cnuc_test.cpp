// physics/ reaction_pathways (docs/upgrade/cnuc_split.md section 10) against the code it replaced: the engine's
// pathway trees (Decay, KGroup, MGroup, ECMGroup on PPair) built by CNuc::SortPathways, which since step 6a fills them
// from reaction_pathways, and built again by a verbatim copy of CNuc::SortPathways at db972c84 (below, `this->`
// written `c.`), on the same loaded CNuc.  Every field of every node, in order, must be equal (doubles with ==).
// Arguments: .azr files or directories (every *.azr below, outside output/); each is loaded, from its own directory, as AZUREAPI does (external capture on), in the Brune and the formal basis.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "Constants.h"
#include "ECAmplitudeCache.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

void OldSortPathways(CNuc &c, const Config &configure) {
  int DecayNum, KGroupNum, MGroupNum;
  for (int aa = 1; aa <= c.NumPairs(); aa++) {
    if (!c.GetPair(aa)->IsEntrance()) continue;
    for (int ir = 1; ir <= c.NumPairs(); ir++) {
      if (c.GetPair(ir)->GetPType() == 20) continue;
      if (c.GetPair(aa)->GetPType() == 20) {
        for (int l = 0; l < 2; l++) {
          for (int j = 1; j <= c.NumJGroups(); j++) {
            if (!c.GetJGroup(j)->IsInRMatrix()) continue;
            for (int ch = 1; ch <= c.GetJGroup(j)->NumChannels(); ch++) {
              if (c.GetJGroup(j)->GetChannel(ch)->GetPairNum() != aa) continue;
              for (int chp = 1; chp <= c.GetJGroup(j)->NumChannels(); chp++) {
                if (c.GetJGroup(j)->GetChannel(chp)->GetPairNum() != ir ||
                    c.GetJGroup(j)->GetChannel(ch)->GetL() != l) continue;
                Decay NewDecay(ir);
                DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                if (!DecayNum) {
                  c.GetPair(aa)->AddDecay(NewDecay);
                  DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                }
                KGroup NewKGroup(l, 0);
                KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                if (!KGroupNum) {
                  c.GetPair(aa)->GetDecay(DecayNum)->AddKGroup(NewKGroup);
                  KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                }
                MGroup NewMGroup(j, ch, chp);
                MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                if (!MGroupNum) {
                  c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->AddMGroup(NewMGroup);
                  MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                }
              }
            }
          }
        }
      } else if (c.GetPair(ir)->GetPType() == 0 && aa == ir) {
        for (double s = fabs(c.GetPair(aa)->GetJ(1) - c.GetPair(aa)->GetJ(2));
             s <= (c.GetPair(aa)->GetJ(1) + c.GetPair(aa)->GetJ(2)); s += 1.) {
          for (double sp = fabs(c.GetPair(ir)->GetJ(1) - c.GetPair(ir)->GetJ(2));
               sp <= (c.GetPair(ir)->GetJ(1) + c.GetPair(ir)->GetJ(2)); sp += 1.) {
            for (int j = 1; j <= c.NumJGroups(); j++) {
              if (!c.GetJGroup(j)->IsInRMatrix()) continue;
              for (int ch = 1; ch <= c.GetJGroup(j)->NumChannels(); ch++) {
                if (c.GetJGroup(j)->GetChannel(ch)->GetPairNum() != aa) continue;
                for (int chp = 1; chp <= c.GetJGroup(j)->NumChannels(); chp++) {
                  if (c.GetJGroup(j)->GetChannel(chp)->GetPairNum() != ir ||
                      c.GetJGroup(j)->GetChannel(ch)->GetS() != s ||
                      c.GetJGroup(j)->GetChannel(chp)->GetS() != sp) continue;
                  Decay NewDecay(ir);
                  DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                  if (!DecayNum) {
                    c.GetPair(aa)->AddDecay(NewDecay);
                    DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                  }
                  KGroup NewKGroup(s, sp);
                  KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                  if (!KGroupNum) {
                    c.GetPair(aa)->GetDecay(DecayNum)->AddKGroup(NewKGroup);
                    KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                  }
                  MGroup NewMGroup(j, ch, chp);
                  MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                  if (!MGroupNum) {
                    c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->AddMGroup(NewMGroup);
                    MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                  }
                  double statspinfactor = (2. * c.GetJGroup(j)->GetJ() + 1.) *
                      c.GetPair(c.GetJGroup(j)->GetChannel(chp)->GetPairNum())->GetI1I2Factor();
                  c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->GetMGroup(MGroupNum)->SetStatSpinFactor(statspinfactor);
                }
              }
            }
          }
        }
      } else if (c.GetPair(ir)->GetPType() == 0 && aa != ir) {
        // Unobserved Primary, Observed Secondary (UPOS): create 3-param KGroups with sp2 for
        // secondary gamma angular distribution calculations.  The secondary-gamma terms
        // multiply T(s') by T(s2')* (Brune and deBoer, PRC 102, 024628, Eq. 22), so a
        // group with sp != sp2 holds the pathways of both exit channel spins; a group with
        // sp == sp2 holds those of sp alone and is the one every other observable uses.
        for (double s = fabs(c.GetPair(aa)->GetJ(1) - c.GetPair(aa)->GetJ(2));
             s <= (c.GetPair(aa)->GetJ(1) + c.GetPair(aa)->GetJ(2)); s += 1.) {
          for (double sp = fabs(c.GetPair(ir)->GetJ(1) - c.GetPair(ir)->GetJ(2));
               sp <= (c.GetPair(ir)->GetJ(1) + c.GetPair(ir)->GetJ(2)); sp += 1.) {
            for (double sp2 = fabs(c.GetPair(ir)->GetJ(1) - c.GetPair(ir)->GetJ(2));
                 sp2 <= (c.GetPair(ir)->GetJ(1) + c.GetPair(ir)->GetJ(2)); sp2 += 1.) {
              for (int j = 1; j <= c.NumJGroups(); j++) {
                if (!c.GetJGroup(j)->IsInRMatrix()) continue;
                for (int ch = 1; ch <= c.GetJGroup(j)->NumChannels(); ch++) {
                  if (c.GetJGroup(j)->GetChannel(ch)->GetPairNum() != aa) continue;
                  for (int chp = 1; chp <= c.GetJGroup(j)->NumChannels(); chp++) {
                    if (c.GetJGroup(j)->GetChannel(chp)->GetPairNum() != ir ||
                        c.GetJGroup(j)->GetChannel(ch)->GetS() != s ||
                        (c.GetJGroup(j)->GetChannel(chp)->GetS() != sp &&
                         c.GetJGroup(j)->GetChannel(chp)->GetS() != sp2)) continue;
                    Decay NewDecay(ir);
                    DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                    if (!DecayNum) {
                      c.GetPair(aa)->AddDecay(NewDecay);
                      DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                    }
                    KGroup NewKGroup(s, sp, sp2);
                    KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup, true);
                    if (!KGroupNum) {
                      c.GetPair(aa)->GetDecay(DecayNum)->AddKGroup(NewKGroup);
                      KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup, true);
                    }
                    MGroup NewMGroup(j, ch, chp);
                    MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                    if (!MGroupNum) {
                      c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->AddMGroup(NewMGroup);
                      MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                    }
                    double statspinfactor = (2. * c.GetJGroup(j)->GetJ() + 1.) *
                        c.GetPair(c.GetJGroup(j)->GetChannel(chp)->GetPairNum())->GetI1I2Factor();
                    c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->GetMGroup(MGroupNum)->SetStatSpinFactor(statspinfactor);
                  }
                }
              }
            }
          }
        }
      } else if (c.GetPair(ir)->GetPType() == 10 && !(configure.paramMask & Config::USE_RMC_FORMALISM)) {
        for (double s = fabs(c.GetPair(aa)->GetJ(1) - c.GetPair(aa)->GetJ(2));
             s <= (c.GetPair(aa)->GetJ(1) + c.GetPair(aa)->GetJ(2)); s += 1.) {
          for (int j = 1; j <= c.NumJGroups(); j++) {
            if (!c.GetJGroup(j)->IsInRMatrix()) continue;
            for (int ch = 1; ch <= c.GetJGroup(j)->NumChannels(); ch++) {
              if (c.GetJGroup(j)->GetChannel(ch)->GetPairNum() != aa) continue;
              for (int chp = 1; chp <= c.GetJGroup(j)->NumChannels(); chp++) {
                if (c.GetJGroup(j)->GetChannel(chp)->GetPairNum() != ir ||
                    c.GetJGroup(j)->GetChannel(ch)->GetS() != s) continue;
                Decay NewDecay(ir);
                DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                if (!DecayNum) {
                  c.GetPair(aa)->AddDecay(NewDecay);
                  DecayNum = c.GetPair(aa)->IsDecay(NewDecay);
                }
                KGroup NewKGroup(s, 0);
                KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                if (!KGroupNum) {
                  c.GetPair(aa)->GetDecay(DecayNum)->AddKGroup(NewKGroup);
                  KGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->IsKGroup(NewKGroup);
                }
                MGroup NewMGroup(j, ch, chp);
                MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                if (!MGroupNum) {
                  c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->AddMGroup(NewMGroup);
                  MGroupNum = c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->IsMGroup(NewMGroup);
                }
                double statspinfactor = (2. * c.GetJGroup(j)->GetJ() + 1.) *
                    c.GetPair(c.GetJGroup(j)->GetChannel(chp)->GetPairNum())->GetI1I2Factor();
                c.GetPair(aa)->GetDecay(DecayNum)->GetKGroup(KGroupNum)->GetMGroup(MGroupNum)->SetStatSpinFactor(statspinfactor);
              }
            }
          }
        }
      }
    }
  }
  for (int aa = 1; aa <= c.NumPairs(); aa++) {  // loop over all pairs
    PPair *entrancePair = c.GetPair(aa);
    if (entrancePair->GetPType() == 20 || !entrancePair->IsEntrance()) continue;
    for (int j = 1; j <= c.NumJGroups(); j++) {
      JGroup *theFinalJGroup = c.GetJGroup(j);
      for (int la = 1; la <= theFinalJGroup->NumLevels(); la++) {
        ALevel *theFinalLevel = theFinalJGroup->GetLevel(la);
        if (!theFinalLevel->IsECLevel()) continue;
        int decayNum = entrancePair->IsDecay(theFinalLevel->GetECPairNum());         // store RESONANCE decay number to final state
        if (!decayNum) continue;                                                     // if this is a resonance decay...
        for (int k = 1; k <= entrancePair->GetDecay(decayNum)->NumKGroups(); k++) {  // loop over all kgroups for decays to final state
          KGroup *theKGroup = entrancePair->GetDecay(decayNum)->GetKGroup(k);
          for (int chp = 1; chp <= theFinalJGroup->NumChannels(); chp++) {  // loop over all final configurations in the capture state
            AChannel *finalChannel = theFinalJGroup->GetChannel(chp);
            if (c.GetPair(finalChannel->GetPairNum())->GetPType() != 0) continue;  // ensure the configuration is a particle pair
            int chDecayNum = entrancePair->IsDecay(finalChannel->GetPairNum());
            if (!chDecayNum) continue;  // if it is actually a resonance decay...
            // The intermediate (channel-capture) decay may carry the 3-parameter UPOS KGroups,
            // where for each (s,sp) there is one KGroup per secondary-channel spin sp2, all
            // holding identical MGroups.  Summing over all of them would double count the
            // channel-capture pathway.  Deduplicate on the physical MGroup channel triple
            // (JNum,ChNum,ChpNum) -- which uniquely determines (s,sp) -- so genuine sp2
            // duplicates are collapsed while distinct spin channels are all retained.
            std::set<std::tuple<int, int, int>> seenChMGroups;
            for (int kp = 1; kp <= entrancePair->GetDecay(chDecayNum)->NumKGroups(); kp++) {
              if (entrancePair->GetDecay(chDecayNum)->GetKGroup(kp)->GetS() != theKGroup->GetS()) continue;
              for (int mp = 1; mp <= entrancePair->GetDecay(chDecayNum)->GetKGroup(kp)->NumMGroups(); mp++) {
                MGroup *chMGroup = entrancePair->GetDecay(chDecayNum)->GetKGroup(kp)->GetMGroup(mp);
                if (!seenChMGroups.insert(std::make_tuple(chMGroup->GetJNum(),
                                                          chMGroup->GetChNum(),
                                                          chMGroup->GetChpNum()))
                         .second) continue;
                AChannel *chChannel = c.GetJGroup(chMGroup->GetJNum())->GetChannel(chMGroup->GetChNum());
                AChannel *chChannelp = c.GetJGroup(chMGroup->GetJNum())->GetChannel(chMGroup->GetChpNum());
                for (int multL = 1; multL <= maxECMult; multL++) {  // loop over all allowed gamma parities
                  char radType;
                  if (c.GetJGroup(chMGroup->GetJNum())->GetPi() * theFinalJGroup->GetPi() == (int)pow(-1, multL))
                    radType = 'E';
                  else
                    radType = 'M';  // calculate radiation type
                  if (!((radType == 'E' && multL == 1) && (theFinalLevel->GetECMultMask() & isE1)) &&
                      !((radType == 'M' && multL == 1) && (theFinalLevel->GetECMultMask() & isM1)) &&
                      !((radType == 'E' && multL == 2) && (theFinalLevel->GetECMultMask() & isE2))) continue;  // allow only m1,e1,e2
                  if (fabs(c.GetJGroup(chMGroup->GetJNum())->GetJ() - multL) > theFinalJGroup->GetJ() ||
                      theFinalJGroup->GetJ() > c.GetJGroup(chMGroup->GetJNum())->GetJ() + multL) continue;
                  if (!(abs(chChannelp->GetL() - multL) <= finalChannel->GetL() &&
                        finalChannel->GetL() <= chChannelp->GetL() + multL &&
                        fabs(chChannelp->GetS() - finalChannel->GetL()) <= theFinalJGroup->GetJ() &&
                        theFinalJGroup->GetJ() <= chChannelp->GetS() + finalChannel->GetL() &&
                        chChannelp->GetS() == finalChannel->GetS()) &&
                      !(fabs(chChannelp->GetS() - multL) <= finalChannel->GetS() &&
                        finalChannel->GetS() <= chChannelp->GetS() + multL &&
                        fabs(chChannelp->GetL() - finalChannel->GetS()) <= theFinalJGroup->GetJ() &&
                        theFinalJGroup->GetJ() <= chChannelp->GetL() + finalChannel->GetS() &&
                        chChannelp->GetL() == finalChannel->GetL() &&
                        radType == 'M')) continue;  // ensure entrance channel for dc can couple to final state
                  if (chChannel == chChannelp) {
                    ECMGroup newECMGroup(radType, multL, chChannel->GetL(),
                                         c.GetJGroup(chMGroup->GetJNum())->GetJ(), chp, j, la);
                    theKGroup->AddECMGroup(newECMGroup);
                  }
                  int internalChannel = 0;
                  for (int intCh = 1; intCh <= c.GetJGroup(chMGroup->GetJNum())->NumChannels(); intCh++) {
                    if (c.GetJGroup(chMGroup->GetJNum())->GetChannel(intCh)->GetRadType() == radType &&
                        c.GetJGroup(chMGroup->GetJNum())->GetChannel(intCh)->GetL() == multL &&
                        c.GetJGroup(chMGroup->GetJNum())->GetChannel(intCh)->GetPairNum() == theFinalLevel->GetECPairNum()) {
                      internalChannel = intCh;
                      break;
                    }
                  }
                  ECMGroup newECMGroup(radType, multL, chChannel->GetL(), c.GetJGroup(chMGroup->GetJNum())->GetJ(),
                                       chp, j, la, chDecayNum, kp, mp, internalChannel);
                  theKGroup->AddECMGroup(newECMGroup);
                }
              }
            }
          }
        }
      }
    }
  }
}

// One line per node, every field at full precision.
std::vector<std::string> Dump(CNuc &c) {
  std::vector<std::string> out;
  char buf[512];
  for (int a = 1; a <= c.NumPairs(); a++) {
    PPair *p = c.GetPair(a);
    for (int d = 1; d <= p->NumDecays(); d++) {
      Decay *dec = p->GetDecay(d);
      std::snprintf(buf, sizeof buf, "pair %d decay %d -> pair %d: %d groups", a, d, dec->GetPairNum(), dec->NumKGroups());
      out.push_back(buf);
      for (int k = 1; k <= dec->NumKGroups(); k++) {
        KGroup *g = dec->GetKGroup(k);
        std::snprintf(buf, sizeof buf, "  k %d: s %.17g s' %.17g s2' %.17g cross %d, %d internal, %d external", k, g->GetS(),
                      g->GetSp(), g->GetSp2(), (int)g->IsCrossSpin(), g->NumMGroups(), g->NumECMGroups());
        out.push_back(buf);
        for (int m = 1; m <= g->NumMGroups(); m++) {
          MGroup *mg = g->GetMGroup(m);
          std::snprintf(buf, sizeof buf, "    m %d: J-group %d ch %d ch' %d stat %.17g", m, mg->GetJNum(), mg->GetChNum(),
                        mg->GetChpNum(), mg->GetStatSpinFactor());
          out.push_back(buf);
        }
        for (int m = 1; m <= g->NumECMGroups(); m++) {
          ECMGroup *e = g->GetECMGroup(m);
          std::snprintf(buf, sizeof buf,
                        "    ec %d: %c%d l %d J %.17g final ch %d J-group %d level %d cc %d/%d/%d photon ch %d stat %.17g", m,
                        e->GetRadType(), e->GetMult(), e->GetL(), e->GetJ(), e->GetFinalChannel(), e->GetJGroupNum(),
                        e->GetLevelNum(), e->GetChanCapDecay(), e->GetChanCapKGroup(), e->GetChanCapMGroup(),
                        e->GetIntChannelNum(), e->GetStatSpinFactor());
          out.push_back(buf);
        }
      }
    }
  }
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  namespace fs = std::filesystem;
  std::vector<std::string> files;
  for (int i = 1; i < argc; i++) {
    if (fs::is_directory(argv[i])) {
      for (const auto &e : fs::recursive_directory_iterator(argv[i]))
        if (e.is_regular_file() && e.path().extension() == ".azr" &&
            e.path().string().find("/output/") == std::string::npos)
          files.push_back(fs::absolute(e.path()).string());
    } else {
      files.push_back(fs::absolute(argv[i]).string());
    }
  }
  std::sort(files.begin(), files.end());
  const fs::path start = fs::current_path();
  int failures = 0;
  long nodes = 0;
  for (const std::string &file : files)
    for (int formal = 0; formal <= 1; formal++) {
      const std::string tag = file + (formal ? " (formal)" : " (Brune)");
      fs::current_path(fs::path(file).parent_path());
      std::ostringstream sink;
      Config cfg(sink);
      cfg.configfile = fs::path(file).filename().string();
      cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
      if (formal) cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
      g_config = &cfg;
      CNuc compound;
      if (cfg.ReadConfigFile() < 0) {
        std::printf("  FAIL  %s: config\n", tag.c_str());
        failures++;
        continue;
      }
      cfg.fileCheckMask = 0;  // writes nothing
      cfg.screenCheckMask = 0;
      cfg.ResolveMatrixRoute();
      InitializeECAmplitudeCache();
      if (compound.Fill(cfg) == -1) {  // a model the engine refuses (tests of the refusal)
        std::printf("  skip  %s: Fill refuses the model\n", tag.c_str());
        continue;
      }
      // Entrances come from the data segments, which are not read here: every particle and beta pair is one, a
      // superset of what any project's data select.
      for (int a = 1; a <= compound.NumPairs(); a++)
        if (compound.GetPair(a)->GetPType() != 10) compound.GetPair(a)->SetEntrance();
      try {
        compound.Initialize(cfg);
      } catch (const std::exception &x) {
        std::printf("  skip  %s: Initialize throws (%s)\n", tag.c_str(), x.what());
        continue;
      }
      const std::vector<std::string> now = Dump(compound);
      for (int a = 1; a <= compound.NumPairs(); a++) compound.GetPair(a)->ClearDecays();
      OldSortPathways(compound, cfg);
      const std::vector<std::string> before = Dump(compound);
      nodes += (long)now.size();
      size_t n = 0;
      while (n < now.size() && n < before.size() && now[n] == before[n]) n++;
      if (n == now.size() && n == before.size()) {
        std::printf("  ok    %s: %zu nodes\n", tag.c_str(), now.size());
      } else {
        std::printf("  FAIL  %s: first difference at node %zu of %zu / %zu\n    now:    %s\n    before: %s\n", tag.c_str(), n,
                    now.size(), before.size(), n < now.size() ? now[n].c_str() : "(end)",
                    n < before.size() ? before[n].c_str() : "(end)");
        failures++;
      }
    }
  fs::current_path(start);
  if (nodes == 0) {
    std::printf("  FAIL  no pathways at all\n");
    failures++;
  }
  std::printf("%s: %zu files, %ld nodes, %d project(s) differ\n", failures ? "FAIL" : "PASS", files.size(), nodes, failures);
  return failures ? 1 : 0;
}
