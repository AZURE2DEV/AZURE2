// physics/ reaction_pathways (docs/upgrade/cnuc_split.md section 10) against the code it replaced: the engine's
// pathway trees (Decay, KGroup, MGroup, ECMGroup on PPair) built by CNuc::SortPathways, which since step 6a fills them
// from reaction_pathways, and built again by a verbatim copy of CNuc::SortPathways at db972c84 (below, `this->`
// written `c.`), on the same loaded CNuc.  Every field of every node, in order, must be equal (doubles with ==).
// Arguments: .azr files or directories (every *.azr below, outside output/); each is loaded, from its own directory, as AZUREAPI does (external capture on), in the Brune and the formal basis.
#include <algorithm>
#include <cstdarg>
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
#include "Decay.h"
#include "ECAmplitudeCache.h"
#include "ECMGroup.h"
#include "KGroup.h"
#include "KLGroup.h"
#include "MGroup.h"
#include "PathwayTable.h"

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

// The pathway table against the trees: one line per row of either, every field, doubles as hexadecimal floats (so
// equal lines mean equal bits).  The capture-A_y terms are the table's alone (built at set-up for the reactions the
// data sets name); physics_pathway_table and upgrade_capture_ay_requests check them.
void Line(std::vector<std::string> &out, const char *format, ...) {
  char buf[1024];
  va_list args;
  va_start(args, format);
  std::vsnprintf(buf, sizeof buf, format, args);
  va_end(args);
  out.push_back(buf);
}

std::vector<std::string> DumpTrees(CNuc &c) {
  std::vector<std::string> out;
  for (int a = 1; a <= c.NumPairs(); a++) {
    PPair *p = c.GetPair(a);
    for (int d = 1; d <= p->NumDecays(); d++) {
      Decay *dec = p->GetDecay(d);
      Line(out, "entrance %d exit %d -> pair %d: %d groups, interferences %d, %zu pathways, %d orders, %zu diagonal",
           a, d, dec->GetPairNum(), dec->NumKGroups(), (int)dec->TablesBuilt(), dec->InterferencePathways().size(),
           dec->NumInterferenceOrders(), dec->DiagonalPathways().size());
      for (int k = 1; k <= dec->NumKGroups(); k++) {
        KGroup *g = dec->GetKGroup(k);
        Line(out, "  k %d: s %a s' %a s2' %a cross %d, %d internal, %d external", k, g->GetS(), g->GetSp(), g->GetSp2(),
             (int)g->IsCrossSpin(), g->NumMGroups(), g->NumECMGroups());
        for (int m = 1; m <= g->NumMGroups(); m++) {
          MGroup *mg = g->GetMGroup(m);
          Line(out, "    m %d: J-group %d ch %d ch' %d stat %a", m, mg->GetJNum(), mg->GetChNum(), mg->GetChpNum(),
               mg->GetStatSpinFactor());
        }
        for (int m = 1; m <= g->NumECMGroups(); m++) {
          ECMGroup *e = g->GetECMGroup(m);
          std::string through = "direct";
          if (e->IsChannelCapture()) {
            MGroup *cc = p->GetDecay(e->GetChanCapDecay())->GetKGroup(e->GetChanCapKGroup())->GetMGroup(e->GetChanCapMGroup());
            through = "through J-group " + std::to_string(cc->GetJNum()) + " ch " + std::to_string(cc->GetChNum()) +
                      " ch' " + std::to_string(cc->GetChpNum());
          }
          Line(out, "    ec %d: %c%d l %d J %a final ch %d J-group %d level %d cc %d %d/%d/%d %s photon ch %d", m,
               e->GetRadType(), e->GetMult(), e->GetL(), e->GetJ(), e->GetFinalChannel(), e->GetJGroupNum(),
               e->GetLevelNum(), (int)e->IsChannelCapture(), e->GetChanCapDecay(), e->GetChanCapKGroup(),
               e->GetChanCapMGroup(), through.c_str(), e->GetIntChannelNum());
        }
      }
      // The (k, L) groups in order, against the flat list built from them.
      const std::vector<FlatInterference> &flat = dec->FlatInterferences();
      std::size_t i = 0;
      for (int kl = 1; kl <= dec->NumKLGroups(); kl++) {
        const KLGroup *g = dec->GetKLGroup(kl);
        for (int n = 1; n <= g->NumInterferences(); n++, i++) {
          const Interference *itf = g->GetInterference(n);
          const FlatInterference f = i < flat.size() ? flat[i] : FlatInterference{};
          Line(out, "  interference k %d L %d m1 %d m2 %d %s z %a z_upos %a | flat k %d L %d m1 %d m2 %d %c%c z %a "
                    "z_upos %a sp equal %d s1 %d s2 %d",
               g->GetK(), g->GetLOrder(), itf->GetM1(), itf->GetM2(), itf->GetInterferenceType().c_str(),
               itf->GetZ1Z2(), itf->GetZ1Z2_UPOS(), f.k, f.l, f.m1, f.m2, f.r1, f.r2, f.z, f.z_upos, (int)f.spEqual,
               f.s1, f.s2);
        }
      }
      if (i != flat.size()) Line(out, "  flat list of %zu for %zu interferences", flat.size(), i);
      for (const InterferencePathway &q : dec->InterferencePathways()) Line(out, "  pathway k %d m %d %c", q.k, q.m, q.rad);
      for (const DiagonalPathway &q : dec->DiagonalPathways())
        Line(out, "  diagonal k %d m %d s %a l %d stat %a", q.k, q.m, q.s, q.l, q.statSpin);
    }
  }
  return out;
}

std::vector<std::string> DumpTable(CNuc &c) {
  std::vector<std::string> out;
  const azure2::PathwayTable &table = c.Pathways();
  for (int a = 1; a <= c.NumPairs(); a++) {
    const azure2::EntrancePathways &P = table.Entrance(a);
    for (int d = 1; d <= (int)P.exits.size(); d++) {
      const azure2::ExitRow &x = P.exits[d - 1];
      Line(out, "entrance %d exit %d -> pair %d: %d groups, interferences %d, %d pathways, %d orders, %d diagonal", a, d,
           x.pair, x.num_groups, (int)x.has_interferences, x.num_interference_pathways, x.num_orders, x.num_diagonal);
      for (int k = 1; k <= x.num_groups; k++) {
        const azure2::SpinGroupRow &g = P.Group(x, k);
        Line(out, "  k %d: s %a s' %a s2' %a cross %d, %d internal, %d external", k, g.s, g.sp, g.sp2, (int)g.CrossSpin(),
             g.num_internal, g.num_external);
        for (int m = 1; m <= g.num_internal; m++) {
          const azure2::InternalRow &r = P.Internal(g, m);
          Line(out, "    m %d: J-group %d ch %d ch' %d stat %a", m, r.group, r.ch, r.chp, r.stat_spin_factor);
        }
        for (int m = 1; m <= g.num_external; m++) {
          const azure2::ExternalRow &e = P.External(g, m);
          std::string through = "direct";
          if (e.channel_capture) {
            const azure2::InternalRow &cc = P.internal[e.cc_row];
            through = "through J-group " + std::to_string(cc.group) + " ch " + std::to_string(cc.ch) + " ch' " +
                      std::to_string(cc.chp);
            // The row and the position triple name the same pathway.
            const azure2::ExitRow &ccExit = P.exits[e.cc_exit - 1];
            if (P.Group(ccExit, e.cc_group).first_internal + e.cc_pathway - 1 != e.cc_row) through += " (cc_row differs)";
          } else if (e.cc_row != -1) {
            through += " (cc_row set)";
          }
          Line(out, "    ec %d: %c%d l %d J %a final ch %d J-group %d level %d cc %d %d/%d/%d %s photon ch %d", m,
               e.rad_type, e.mult, e.l, e.J, e.final_channel, e.final_group, e.final_level, (int)e.channel_capture,
               e.cc_exit, e.cc_group, e.cc_pathway, through.c_str(), e.photon_channel);
        }
      }
      for (int i = 0; i < x.num_interferences; i++) {
        const azure2::FlatInterference &f = P.interferences[x.first_interference + i];
        const std::string type{f.r1, f.r2};
        Line(out, "  interference k %d L %d m1 %d m2 %d %s z %a z_upos %a | flat k %d L %d m1 %d m2 %d %c%c z %a "
                  "z_upos %a sp equal %d s1 %d s2 %d",
             f.k, f.l, f.m1, f.m2, type.c_str(), f.z, f.z_upos, f.k, f.l, f.m1, f.m2, f.r1, f.r2, f.z, f.z_upos,
             (int)f.spEqual, f.s1, f.s2);
      }
      for (int i = 0; i < x.num_interference_pathways; i++) {
        const azure2::InterferencePathway &q = P.interference_pathways[x.first_interference_pathway + i];
        Line(out, "  pathway k %d m %d %c", q.k, q.m, q.rad);
      }
      for (int i = 0; i < x.num_diagonal; i++) {
        const azure2::DiagonalPathway &q = P.diagonal[x.first_diagonal + i];
        Line(out, "  diagonal k %d m %d s %a l %d stat %a", q.k, q.m, q.s, q.l, q.statSpin);
      }
    }
  }
  return out;
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
  long nodes = 0, tableRows = 0;
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
      {
        const std::vector<std::string> table = DumpTable(compound), trees = DumpTrees(compound);
        tableRows += (long)table.size();
        size_t n = 0;
        while (n < table.size() && n < trees.size() && table[n] == trees[n]) n++;
        if (n != table.size() || n != trees.size()) {
          std::printf("  FAIL  %s: table and trees differ at line %zu of %zu / %zu\n    table: %s\n    trees: %s\n",
                      tag.c_str(), n, table.size(), trees.size(), n < table.size() ? table[n].c_str() : "(end)",
                      n < trees.size() ? trees[n].c_str() : "(end)");
          failures++;
        }
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
  if (tableRows == 0) {
    std::printf("  FAIL  no table rows at all\n");
    failures++;
  }
  if (nodes == 0) {
    std::printf("  FAIL  no pathways at all\n");
    failures++;
  }
  std::printf("%s: %zu files, %ld nodes, %ld table lines, %d project(s) differ\n", failures ? "FAIL" : "PASS", files.size(),
              nodes, tableRows, failures);
  return failures ? 1 : 0;
}
