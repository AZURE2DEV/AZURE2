// Properties the pathway table (PathwayTable.h, docs/upgrade/cnuc_split.md section 11.6 step 1) must have, checked
// from the loaded CNuc without the builder's loops:
//  - every internal row joins two channels of one J-group: the entrance channel of the entrance pair with the group's
//    entrance spin s, the exit channel of the exit pair with spin s' (or s2' of a UPOS group); its System indices are
//    the J-group's channel positions;
//  - every external row names an external-capture final level, a channel of it, and for channel capture an internal
//    row of the same entrance whose exit pair is the final channel's pair, also by its (exit, k, m) position;
//  - every interference, interference pathway and diagonal pathway names rows that exist, with the stated kind;
//  - for a photon exit, a_k = 4 Z1Z2 for the pathway pairs of one spin group: the Legendre coefficient of Seyler and
//    Weller's Eq. (20) is the unpolarized coefficient of the differential cross section (CalculateCaptureAnalyzingPower
//    states the relation); every (k, L, m1, m2) with an interference has its capture-A_y term and the converse.
// Arguments: .azr files or directories (every *.azr below, outside output/), loaded as pathways_from_cnuc_test does
// (every particle and beta pair an entrance, capture A_y requested for every photon exit), Brune basis.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "JGroup.h"
#include "PPair.h"
#include "PathwayTable.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

int failures = 0;
long checks = 0, akChecks = 0;

void Check(bool ok, const std::string &tag, const std::string &what) {
  checks++;
  if (ok) return;
  if (failures < 40) std::printf("  FAIL  %s: %s\n", tag.c_str(), what.c_str());
  failures++;
}

std::string Where(int a, int e, int k, int m) {
  return "entrance " + std::to_string(a) + " exit " + std::to_string(e) + " k " + std::to_string(k) + " m " +
         std::to_string(m);
}

void CheckEntrance(CNuc &c, int a, const std::string &tag) {
  const azure2::EntrancePathways &P = c.Pathways().Entrance(a);
  const bool beta = c.GetPair(a)->GetPType() == 20;
  for (int e = 1; e <= (int)P.exits.size(); e++) {
    const azure2::ExitRow &x = P.exits[e - 1];
    Check(x.pair >= 1 && x.pair <= c.NumPairs() && P.ExitOf(x.pair) == e, tag, "exit " + std::to_string(e) + " pair");
    Check(x.first_group + x.num_groups <= (int)P.groups.size(), tag, "exit " + std::to_string(e) + " groups");
    for (int k = 1; k <= x.num_groups; k++) {
      const azure2::SpinGroupRow &g = P.Group(x, k);
      for (int m = 1; m <= g.num_internal; m++) {
        const azure2::InternalRow &r = P.Internal(g, m);
        const std::string at = Where(a, e, k, m);
        const bool group = r.group >= 1 && r.group <= c.NumJGroups();
        Check(group, tag, at + ": J-group");
        if (!group) continue;
        JGroup *J = c.GetJGroup(r.group);
        const bool channels = r.ch >= 1 && r.ch <= J->NumChannels() && r.chp >= 1 && r.chp <= J->NumChannels();
        Check(channels && J->IsInRMatrix(), tag, at + ": channels of an R-matrix group");
        if (!channels) continue;
        AChannel *in = J->GetChannel(r.ch), *out = J->GetChannel(r.chp);
        Check(in->GetPairNum() == a && out->GetPairNum() == x.pair, tag, at + ": channel pairs");
        Check(r.k_in == r.ch - 1 && r.k_out == r.chp - 1, tag, at + ": block channels");
        if (!beta) {
          Check(in->GetS() == g.s, tag, at + ": entrance channel spin");
          // A photon exit's groups are (s, 0): its channel spin is the entrance spin.
          if (c.GetPair(x.pair)->GetPType() == 0)
            Check(out->GetS() == g.sp || (g.has_sp2 && out->GetS() == g.sp2), tag, at + ": exit channel spin");
          else
            Check(g.sp == 0.0, tag, at + ": photon group (s, 0)");
        } else {
          Check(in->GetL() == (int)std::lround(g.s), tag, at + ": beta group l");
        }
      }
      for (int m = 1; m <= g.num_external; m++) {
        const azure2::ExternalRow &r = P.External(g, m);
        const std::string at = Where(a, e, k, m) + " (external)";
        const bool level = r.final_group >= 1 && r.final_group <= c.NumJGroups() && r.final_level >= 1 &&
                           r.final_level <= c.GetJGroup(r.final_group)->NumLevels();
        Check(level && c.GetJGroup(r.final_group)->GetLevel(r.final_level)->IsECLevel(), tag, at + ": final level");
        if (!level) continue;
        JGroup *F = c.GetJGroup(r.final_group);
        Check(c.GetJGroup(r.final_group)->GetLevel(r.final_level)->GetECPairNum() == x.pair, tag, at + ": photon pair");
        const bool finalChannel = r.final_channel >= 1 && r.final_channel <= F->NumChannels();
        Check(finalChannel, tag, at + ": final channel");
        if (!finalChannel) continue;
        if (!r.channel_capture) {
          Check(r.cc_row == -1 && r.cc_exit == 0, tag, at + ": direct capture names no internal row");
          continue;
        }
        const bool row = r.cc_row >= 0 && r.cc_row < (int)P.internal.size();
        Check(row, tag, at + ": channel-capture row");
        if (!row) continue;
        const azure2::InternalRow &cc = P.internal[r.cc_row];
        Check(c.GetJGroup(cc.group)->GetChannel(cc.chp)->GetPairNum() == F->GetChannel(r.final_channel)->GetPairNum(),
              tag, at + ": channel capture through the final channel's pair");
        const bool position = r.cc_exit >= 1 && r.cc_exit <= (int)P.exits.size() && r.cc_group >= 1 &&
                              r.cc_group <= P.exits[r.cc_exit - 1].num_groups && r.cc_pathway >= 1 &&
                              r.cc_pathway <= P.Group(P.exits[r.cc_exit - 1], r.cc_group).num_internal;
        Check(position && &P.Internal(P.Group(P.exits[r.cc_exit - 1], r.cc_group), r.cc_pathway) == &cc, tag,
              at + ": position triple and row");
        if (r.photon_channel) {
          JGroup *B = c.GetJGroup(cc.group);
          Check(r.photon_channel <= B->NumChannels() && B->GetChannel(r.photon_channel)->GetPairNum() == x.pair &&
                    B->GetChannel(r.photon_channel)->GetRadType() == r.rad_type &&
                    B->GetChannel(r.photon_channel)->GetL() == r.mult,
                tag, at + ": photon channel");
        }
      }
    }
    // Interferences and their pathways.
    int orders = 0;
    std::map<std::tuple<int, int, int, char, int, char>, double> z;  // (L, k, m1, r1, m2, r2) -> Z1Z2
    for (int i = 0; i < x.num_interferences; i++) {
      const azure2::FlatInterference &f = P.interferences[x.first_interference + i];
      const std::string at = "exit " + std::to_string(e) + " interference " + std::to_string(i);
      const bool k = f.k >= 1 && f.k <= x.num_groups;
      Check(k, tag, at + ": group");
      if (!k) continue;
      const azure2::SpinGroupRow &g = P.Group(x, f.k);
      auto exists = [&g](int m, char r) { return m >= 1 && m <= (r == 'R' ? g.num_internal : g.num_external); };
      Check((f.r1 == 'R' || f.r1 == 'E') && (f.r2 == 'R' || f.r2 == 'E') && exists(f.m1, f.r1) && exists(f.m2, f.r2), tag,
            at + ": pathways");
      Check(f.s1 >= 0 && f.s1 < x.num_interference_pathways && f.s2 >= 0 && f.s2 < x.num_interference_pathways, tag,
            at + ": pathway positions");
      if (f.s1 >= 0 && f.s1 < x.num_interference_pathways) {
        const azure2::InterferencePathway &p = P.interference_pathways[x.first_interference_pathway + f.s1];
        Check(p.k == f.k && p.m == f.m1 && p.rad == f.r1, tag, at + ": first pathway position");
      }
      if (f.s2 >= 0 && f.s2 < x.num_interference_pathways) {
        const azure2::InterferencePathway &p = P.interference_pathways[x.first_interference_pathway + f.s2];
        Check(p.k == f.k && p.m == f.m2 && p.rad == f.r2, tag, at + ": second pathway position");
      }
      Check(std::fabs(f.z) > 1e-10 || std::fabs(f.z_upos) > 1e-10, tag, at + ": kept coefficient");
      orders = std::max(orders, f.l + 1);
      z[std::make_tuple(f.l, f.k, f.m1, f.r1, f.m2, f.r2)] = f.z;
    }
    Check(orders == x.num_orders, tag, "exit " + std::to_string(e) + ": number of orders");
    Check(x.has_interferences == !beta, tag, "exit " + std::to_string(e) + ": interferences only for a non-beta entrance");
    for (int i = 0; i < x.num_diagonal; i++) {
      const azure2::DiagonalPathway &d = P.diagonal[x.first_diagonal + i];
      const bool k = d.k >= 1 && d.k <= x.num_groups && d.m >= 1 && d.m <= P.Group(x, d.k).num_internal;
      Check(k && P.Internal(P.Group(x, d.k), d.m).ch == P.Internal(P.Group(x, d.k), d.m).chp, tag,
            "exit " + std::to_string(e) + ": diagonal pathway");
    }
    // a_k = 4 Z1Z2 within one spin group of a photon exit.
    if (beta || c.GetPair(x.pair)->GetPType() != 10 || !x.has_capture_ay || x.num_capture_ay == 0) continue;
    std::map<std::tuple<int, int, int, char, int, char>, double> ak;
    for (int t = 0; t < x.num_capture_ay; t++) {
      const azure2::CaptureAyTerm &q = P.capture_ay[x.first_capture_ay + t];
      if (q.kGroup1 != q.kGroup2 || std::fabs(q.ak) < 1e-10) continue;
      ak[std::make_tuple(q.kOrder, q.kGroup1, q.path1, q.isEC1 ? 'E' : 'R', q.path2, q.isEC2 ? 'E' : 'R')] = q.ak;
    }
    for (const auto &p : z) {
      const auto it = ak.find(p.first);
      akChecks++;
      Check(it != ak.end() && std::fabs(it->second - 4.0 * p.second) <= 1e-12 * std::max(1.0, std::fabs(it->second)), tag,
            "exit " + std::to_string(e) + ": a_k = 4 Z1Z2 at L " + std::to_string(std::get<0>(p.first)));
    }
    for (const auto &p : ak)
      Check(z.count(p.first) == 1, tag,
            "exit " + std::to_string(e) + ": a_k without Z1Z2: L " + std::to_string(std::get<0>(p.first)) + " k " +
                std::to_string(std::get<1>(p.first)) + " " + std::get<3>(p.first) + std::to_string(std::get<2>(p.first)) +
                " " + std::get<5>(p.first) + std::to_string(std::get<4>(p.first)) + " a_k " + std::to_string(p.second));
  }
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
  int projects = 0;
  for (const std::string &file : files) {
    fs::current_path(fs::path(file).parent_path());
    std::ostringstream sink;
    Config cfg(sink);
    cfg.configfile = fs::path(file).filename().string();
    cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
    g_config = &cfg;
    CNuc compound;
    if (cfg.ReadConfigFile() < 0) {
      std::printf("  FAIL  %s: config\n", file.c_str());
      failures++;
      continue;
    }
    cfg.fileCheckMask = 0;
    cfg.screenCheckMask = 0;
    cfg.ResolveMatrixRoute();
    InitializeECAmplitudeCache();
    if (compound.Fill(cfg) == -1) continue;
    for (int a = 1; a <= compound.NumPairs(); a++)
      if (compound.GetPair(a)->GetPType() != 10) compound.GetPair(a)->SetEntrance();
    for (int a = 1; a <= compound.NumPairs(); a++)
      for (int p = 1; p <= compound.NumPairs(); p++)
        if (compound.GetPair(a)->GetPType() != 10 && compound.GetPair(p)->GetPType() == 10) compound.RequestCaptureAy(a, p);
    try {
      compound.Initialize(cfg);
    } catch (const std::exception &x) {
      continue;
    }
    projects++;
    for (int a = 1; a <= compound.NumPairs(); a++) CheckEntrance(compound, a, file);
  }
  fs::current_path(start);
  if (akChecks == 0 || projects == 0) {
    std::printf("  FAIL  nothing checked (%d projects, %ld a_k comparisons)\n", projects, akChecks);
    failures++;
  }
  std::printf("%s: %d projects, %ld checks (%ld a_k = 4 Z1Z2), %d failed\n", failures ? "FAIL" : "PASS", projects, checks,
              akChecks, failures);
  return failures ? 1 : 0;
}
