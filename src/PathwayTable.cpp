// The pathway table of a run (docs/upgrade/cnuc_split.md section 11.4): rows from physics/ reaction_pathways in its
// order, and the lists derived at set-up from the coupling coefficients of physics/ ObservableCoupling.
#include "PathwayTable.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "SystemOfCompound.h"
#include "azure2/physics/ObservableCoupling.h"

namespace azure2 {

namespace phys = azure2::physics;

int EntrancePathways::ExitOf(int exitPair) const {
  for (std::size_t e = 0; e < exits.size(); e++)
    if (exits[e].pair == exitPair) return (int)e + 1;
  return 0;
}

namespace {

/// What the coefficient loops need of one pathway of a spin group, internal or external, in the units the trees
/// gave them (2J, l, l' or L, the exit channel spin of an internal pathway, the radiation type).
struct PathwayQuantumNumbers {
  int two_J;
  int l, lp;            ///< entrance l; exit l' (internal) or multipolarity L (external)
  int two_sp_channel;   ///< twice the exit channel spin of an internal pathway; -2 for an external one
  int w;                ///< 0 for a magnetic photon or a particle exit channel, 1 otherwise (electric, beta)
  bool electric;        ///< the exit is an electric photon (capture A_y)
};

PathwayQuantumNumbers OfInternal(const phys::System &system, const InternalRow &r) {
  const phys::JpiBlock &block = system.block(r.block);
  const phys::Channel &in = system.channel(block.channel(r.k_in));
  const phys::Channel &out = system.channel(block.channel(r.k_out));
  PathwayQuantumNumbers q;
  q.two_J = block.two_J();
  q.l = in.l;
  q.lp = out.l;
  q.two_sp_channel = out.two_s;
  q.w = (out.kind == phys::ChannelKind::Magnetic || out.kind == phys::ChannelKind::Particle) ? 0 : 1;
  q.electric = (out.kind == phys::ChannelKind::Electric);
  return q;
}

PathwayQuantumNumbers OfExternal(const ExternalRow &x) {
  PathwayQuantumNumbers q;
  q.two_J = (int)std::lround(2.0 * x.J);
  q.l = x.l;
  q.lp = x.mult;
  q.two_sp_channel = -2;  // as the trees' -1.: never equal to a channel spin
  q.w = (x.rad_type == 'M') ? 0 : 1;
  q.electric = (x.rad_type == 'E');
  return q;
}

int TwiceOf(double j) { return (int)std::lround(2.0 * j); }

/// Pathway m (1-based, internal first, then external) of group g.
PathwayQuantumNumbers OfPathway(const phys::System &system, const EntrancePathways &P, const SpinGroupRow &g, int m) {
  if (m > g.num_internal) return OfExternal(P.External(g, m - g.num_internal));
  return OfInternal(system, P.Internal(g, m));
}

/// The interferences of exit `exit` (loops: group, order, pathway, pathway, the engine's order since before the
/// table), then its flat pathway list and its diagonal pathways.
void AddInterferences(const phys::System &system, int entrancePair, int maxL, EntrancePathways &P, ExitRow &exit) {
  const phys::Partition &exitPartition = system.partition(exit.pair - 1);
  exit.has_interferences = true;
  exit.first_interference = (int)P.interferences.size();
  exit.first_interference_pathway = (int)P.interference_pathways.size();
  exit.num_orders = 0;
  auto position = [&P, &exit](int k, int m, char rad) {
    for (std::size_t i = exit.first_interference_pathway; i < P.interference_pathways.size(); i++) {
      const InterferencePathway &p = P.interference_pathways[i];
      if (p.k == k && p.m == m && p.rad == rad) return (int)i - exit.first_interference_pathway;
    }
    P.interference_pathways.push_back(InterferencePathway{k, m, rad});
    return (int)P.interference_pathways.size() - 1 - exit.first_interference_pathway;
  };
  for (int k = 1; k <= exit.num_groups; k++) {
    const SpinGroupRow &g = P.Group(exit, k);
    const int n = g.num_internal + g.num_external;
    const int two_s = TwiceOf(g.s), two_sp = TwiceOf(g.sp), two_sp2 = TwiceOf(g.sp2);
    for (int lOrder = 0; lOrder <= maxL; lOrder++)
      for (int m1 = 1; m1 <= n; m1++)
        for (int m2 = 1; m2 <= n; m2++) {
          const PathwayQuantumNumbers a = OfPathway(system, P, g, m1), b = OfPathway(system, P, g, m2);
          if ((a.l + b.l + lOrder) % 2 != 0 || (a.lp + b.lp + a.w + b.w + lOrder) % 2 != 0) continue;
          double z = 0.0;
          if (exitPartition.kind == phys::PartitionKind::Particle) {
            z = phys::particle_interference(a.l, a.lp, a.two_J, b.l, b.lp, b.two_J, two_s, two_sp, lOrder);
            if (g.CrossSpin()) z = 0.;
          } else if (exitPartition.kind == phys::PartitionKind::Photon) {
            z = phys::photon_interference(a.l, a.lp, a.two_J, b.l, b.lp, b.two_J, two_s, exitPartition.two_j[1], lOrder);
          }
          // Eq. (22) of Brune and deBoer, PRC 102, 024628 (2020), light particle coupled first (see
          // physics::upos_interference).
          double zUpos = 0.;
          if (exitPartition.kind == phys::PartitionKind::Particle && entrancePair != exit.pair &&
              a.two_sp_channel == two_sp && b.two_sp_channel == two_sp2)
            zUpos = phys::upos_interference(a.l, a.lp, a.two_J, b.l, b.two_J, two_s, two_sp, two_sp2,
                                            exitPartition.two_j[0], exitPartition.two_j[1], lOrder);
          // A capture cascade: the secondary gamma ray of the capture FINAL state after an unobserved primary
          // capture photon.  A photon summed over its direction and both helicities leaves no interference between
          // multipoles (L, electric / magnetic), and the orientation it leaves in the final state is that left by a
          // spin-0 particle of orbital angular momentum L: Eq. (22) with l' = L, I_b = 0 and s' = s2' = I_B = J_f
          // (the conclusions of Brune and deBoer name this extension).  External-capture pathways are capture
          // amplitudes of a given multipole like any other and interfere coherently.
          if (exitPartition.kind == phys::PartitionKind::Photon && a.lp == b.lp && a.w == b.w) {
            const int two_jf = exitPartition.two_j[1];
            zUpos = phys::upos_interference(a.l, a.lp, a.two_J, b.l, b.two_J, two_s, two_jf, two_jf, 0, two_jf, lOrder);
          }
          if (!(std::fabs(z) > 1e-10 || std::fabs(zUpos) > 1e-10)) continue;
          FlatInterference f;
          f.k = k;
          f.r1 = (m1 > g.num_internal) ? 'E' : 'R';
          f.r2 = (m2 > g.num_internal) ? 'E' : 'R';
          f.m1 = (f.r1 == 'E') ? m1 - g.num_internal : m1;
          f.m2 = (f.r2 == 'E') ? m2 - g.num_internal : m2;
          f.l = lOrder;
          f.z = z;
          f.z_upos = zUpos;
          f.spEqual = (g.sp == g.sp2);
          f.s1 = position(k, f.m1, f.r1);
          f.s2 = position(k, f.m2, f.r2);
          exit.num_orders = std::max(exit.num_orders, lOrder + 1);
          P.interferences.push_back(f);
        }
  }
  exit.num_interferences = (int)P.interferences.size() - exit.first_interference;
  exit.num_interference_pathways = (int)P.interference_pathways.size() - exit.first_interference_pathway;

  exit.first_diagonal = (int)P.diagonal.size();
  for (int k = 1; k <= exit.num_groups; k++) {
    const SpinGroupRow &g = P.Group(exit, k);
    for (int m = 1; m <= g.num_internal; m++) {
      const InternalRow &r = P.Internal(g, m);
      if (r.ch != r.chp) continue;
      const phys::Channel &c = system.channel(system.block(r.block).channel(r.k_in));
      P.diagonal.push_back(DiagonalPathway{k, m, c.two_s / 2., c.l, r.stat_spin_factor});
    }
  }
  exit.num_diagonal = (int)P.diagonal.size() - exit.first_diagonal;
}

/// The capture analyzing-power terms of a photon exit (Seyler and Weller, PRC 20 (1979) 453, Eqs. 20, 21): every
/// ordered pair of its pathways, internal and external, across spin groups, for k = 0 ... maxL.  None for a spinless
/// projectile or a beta entrance.
void AddCaptureAy(const phys::System &system, int entrancePair, int maxL, EntrancePathways &P, ExitRow &exit) {
  exit.has_capture_ay = true;
  exit.first_capture_ay = (int)P.capture_ay.size();
  exit.num_capture_ay = 0;
  const phys::Partition &entrance = system.partition(entrancePair - 1);
  const phys::Partition &exitPartition = system.partition(exit.pair - 1);
  // A beta entrance has no projectile to polarize; the readers refuse it before they look at the terms.
  if (entrance.kind != phys::PartitionKind::Particle) return;
  if (exitPartition.kind != phys::PartitionKind::Photon) return;
  const int two_x = entrance.two_j[0], two_a = entrance.two_j[1], two_c = exitPartition.two_j[1];
  if (two_x <= 0) return;  // no vector polarization to have
  struct Pathway {
    int k, m;
    bool ec;
    phys::CapturePathway t;
  };
  std::vector<Pathway> paths;
  for (int k = 1; k <= exit.num_groups; k++) {
    const SpinGroupRow &g = P.Group(exit, k);
    const int two_s = TwiceOf(g.s);
    for (int m = 1; m <= g.num_internal; m++) {
      const PathwayQuantumNumbers q = OfInternal(system, P.Internal(g, m));
      paths.push_back(Pathway{k, m, false, {q.l, q.lp, q.two_J, two_s, q.electric ? 1 : 0}});
    }
    for (int m = 1; m <= g.num_external; m++) {
      const PathwayQuantumNumbers q = OfExternal(P.External(g, m));
      paths.push_back(Pathway{k, m, true, {q.l, q.lp, q.two_J, two_s, q.electric ? 1 : 0}});
    }
  }
  for (const Pathway &t1 : paths)
    for (const Pathway &t2 : paths)
      for (int k = 0; k <= maxL; k++) {
        const phys::CaptureAyCoefficients c = phys::capture_ay_coefficients(t1.t, t2.t, k, two_x, two_a, two_c);
        if (std::fabs(c.a) < 1.e-10 && std::fabs(c.b) < 1.e-10) continue;
        P.capture_ay.push_back(CaptureAyTerm{k, t1.k, t1.m, t1.ec, t2.k, t2.m, t2.ec, c.a, c.b});
      }
  exit.num_capture_ay = (int)P.capture_ay.size() - exit.first_capture_ay;
}

}  // namespace

std::shared_ptr<const PathwayTable> BuildPathwayTable(const SystemOfCompound &structure,
                                                      const std::vector<PathwayFinalState> &finalStates,
                                                      const PathwayTableOptions &options) {
  const phys::System &system = structure.system;
  const SystemOfCompoundMap &map = structure.map;
  for (std::size_t a = 0; a < map.pairOfPartition.size(); a++)
    if (map.pairOfPartition[a] != (int)a + 1)
      throw std::logic_error("BuildPathwayTable: partition " + std::to_string(a) + " is not pair " +
                             std::to_string(a + 1));
  std::vector<phys::FinalState> states;
  for (const PathwayFinalState &f : finalStates) states.push_back(f.state);
  phys::PathwayOptions pathwayOptions;
  pathwayOptions.rmatrix_capture = options.rmatrix_capture;

  auto table = std::make_shared<PathwayTable>();
  table->entrances_.resize(system.num_partitions());
  for (int aa : options.entrances) {
    if (aa < 1 || aa > (int)system.num_partitions())
      throw std::logic_error("BuildPathwayTable: entrance pair " + std::to_string(aa) + " does not exist");
    EntrancePathways &P = table->entrances_[aa - 1];
    if (!P.exits.empty()) continue;  // named twice
    const std::vector<phys::ExitPathways> exits = phys::reaction_pathways(system, aa - 1, states, pathwayOptions);
    // Internal pathways first, exit by exit (the external rows refer to them by position).
    std::vector<std::vector<int> > groupRow(exits.size());  // index into P.groups of exit e's group g
    for (std::size_t e = 0; e < exits.size(); e++) {
      ExitRow x{};
      x.pair = map.pairOfPartition[exits[e].partition];
      x.first_group = (int)P.groups.size();
      x.num_groups = (int)exits[e].groups.size();
      for (const phys::SpinGroup &g : exits[e].groups) {
        SpinGroupRow row{};
        row.s = g.two_s / 2.;
        row.sp = g.two_sp / 2.;
        row.has_sp2 = g.has_sp2;
        row.sp2 = g.has_sp2 ? g.two_sp2 / 2. : 0.;
        row.first_internal = (int)P.internal.size();
        row.num_internal = (int)g.internal.size();
        for (const phys::InternalPathway &p : g.internal)
          P.internal.push_back(InternalRow{map.groupOfBlock[p.block], map.channelOfBlockChannel[p.block][p.k_in],
                                           map.channelOfBlockChannel[p.block][p.k_out],
                                           p.has_stat_spin_factor ? p.stat_spin_factor : 0.0, (int)p.block,
                                           (int)p.k_in, (int)p.k_out});
        groupRow[e].push_back((int)P.groups.size());
        P.groups.push_back(row);
      }
      P.exits.push_back(x);
    }
    // External pathways, group by group.
    for (std::size_t e = 0; e < exits.size(); e++)
      for (std::size_t g = 0; g < exits[e].groups.size(); g++) {
        SpinGroupRow &row = P.groups[groupRow[e][g]];
        row.first_external = (int)P.external.size();
        row.num_external = (int)exits[e].groups[g].external.size();
        for (const phys::ExternalPathway &x : exits[e].groups[g].external) {
          ExternalRow r{};
          r.rad_type = (x.X == phys::Multipole::E) ? 'E' : 'M';
          r.mult = x.L;
          r.l = x.l;
          r.J = x.two_J / 2.;
          r.final_group = finalStates[x.final_state].group;
          r.final_level = finalStates[x.final_state].level;
          r.final_channel = (int)x.final_channel + 1;
          r.channel_capture = x.channel_capture;
          r.cc_row = -1;
          r.through_row = P.groups[groupRow[x.cc_exit][x.cc_group]].first_internal + (int)x.cc_pathway;
          if (x.channel_capture) {
            r.cc_exit = (int)x.cc_exit + 1;
            r.cc_group = (int)x.cc_group + 1;
            r.cc_pathway = (int)x.cc_pathway + 1;
            r.cc_row = P.groups[groupRow[x.cc_exit][x.cc_group]].first_internal + (int)x.cc_pathway;
            const phys::InternalPathway &cc = exits[x.cc_exit].groups[x.cc_group].internal[x.cc_pathway];
            r.photon_channel = x.photon_channel < 0 ? 0 : map.channelOfBlockChannel[cc.block][x.photon_channel];
          }
          P.external.push_back(r);
        }
      }
    // The derived lists; a beta entrance has no angular distribution.
    if (system.partition(aa - 1).kind == phys::PartitionKind::Beta) continue;
    for (ExitRow &exit : P.exits) AddInterferences(system, aa, options.max_l_order, P, exit);
  }
  // The requested capture-A_y terms, entrance by entrance and exit by exit in the table's order, so that the table
  // does not depend on the order of the requests.  A request for a pair no pathway reaches builds nothing.
  for (int aa = 1; aa <= (int)system.num_partitions(); aa++) {
    EntrancePathways &P = table->entrances_[aa - 1];
    for (ExitRow &exit : P.exits) {
      bool requested = false;
      for (const std::pair<int, int> &request : options.capture_ay)
        requested = requested || (request.first == aa && request.second == exit.pair);
      if (requested) AddCaptureAy(system, aa, options.max_l_order, P, exit);
    }
  }
  return table;
}

}  // namespace azure2
