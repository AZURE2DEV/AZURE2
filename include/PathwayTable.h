#ifndef PATHWAYTABLE_H
#define PATHWAYTABLE_H

#include <memory>
#include <utility>
#include <vector>

#include "azure2/physics/Pathways.h"

struct CNucSystem;

/*!
 * The reaction pathways of a run as one immutable table (docs/upgrade/cnuc_split.md section 11.4): per entrance pair,
 * the exits, their channel-spin groups, the internal and external pathways, and the lists derived from them at
 * set-up (the interferences of the differential cross section, the diagonal pathways of the Coulomb-nuclear
 * interference, the capture analyzing-power terms).  Built once by BuildPathwayTable from physics/ reaction_pathways
 * and the coupling coefficients; nothing writes it afterwards, and every working copy of the compound shares it.
 *
 * Rows are in the order of reaction_pathways (exit -> spin group -> internal, then external pathways), the order of
 * the pathway trees the table replaced (until 2026-11-11), so every sum over pathways adds the same terms in the same
 * order as before and the signature of an intEC file is unchanged.
 * The positions (k, m) -- spin group within its exit and pathway within its group, 1-based, internal and external
 * counted separately -- key the T-matrix storage, the external-capture amplitudes and their cache; the table keeps
 * them.  Engine numbers (J-group, channel, level) are 1-based as CNuc numbers them; System indices are 0-based.
 */
namespace azure2 {

/// One interference term of the differential cross section, in the order of the (k, L) groups.
struct FlatInterference {
  int k;          ///< spin group (1-based)
  int m1, m2;     ///< pathways (1-based) in the list r1, r2 name
  int l;          ///< Legendre order
  char r1, r2;    ///< 'R' internal, 'E' external capture pathway
  double z;       ///< Z1Z2 coefficient (0 for a cross-spin UPOS group)
  double z_upos;  ///< the unobserved-primary, observed-secondary coefficient (0 unless a particle exit != entrance)
  bool spEqual;   ///< the group has s' equal to s2' (s2' = 0 for a group without one)
  int s1, s2;     ///< positions of the two pathways in the exit's interference pathways (0-based)
};

/// A pathway taking part in an interference, in order of first appearance.
struct InterferencePathway {
  int k, m;   ///< spin group and pathway (1-based)
  char rad;   ///< 'R' internal, 'E' external capture
};

/// An internal pathway whose entrance and exit channel are the same.
struct DiagonalPathway {
  int k, m;         ///< spin group and pathway (1-based)
  double s;         ///< channel spin
  int l;            ///< orbital angular momentum
  double statSpin;  ///< the pathway's statistical spin factor
};

/*!
 * One (t, t') pathway pair of a capture analyzing power: R. G. Seyler and H. R. Weller, Phys. Rev. C 20 (1979) 453,
 * Eqs. (20) and (21) in the channel-spin representation,
 *   sigma(theta, phi) = N sum_k [a_k P_k(cos theta) + b_k P_k^1(cos theta) p_y],
 *   A_y(theta) = sum_k b_k P_k^1 / sum_k a_k P_k,
 * with a_k multiplying Re T T'* and b_k Re (i T T'*).  The two pathways may lie in different spin groups (b_k does
 * not need s = s'), so each side carries its own group.
 */
struct CaptureAyTerm {
  int kOrder;   ///< Legendre order k
  int kGroup1;  ///< spin group of the first pathway (1-based)
  int path1;    ///< internal (isEC1 false) or external (true) pathway, 1-based
  bool isEC1;
  int kGroup2;
  int path2;
  bool isEC2;
  double ak;    ///< Eq. (20)
  double bk;    ///< Eq. (21)
};

/// An internal pathway: entrance channel ch and exit channel chp of one J-group.
struct InternalRow {
  int group, ch, chp;         ///< J-group and its channel numbers (CNuc's, 1-based)
  double stat_spin_factor;    ///< (2J + 1) / ((2 I1 + 1)(2 I2 + 1)) of the exit pair as the trees hold it; 0.0 for a
                              ///< beta entrance
  int block, k_in, k_out;     ///< physics::InternalPathway (System indices, 0-based)
};

/// An external (direct- or channel-capture) pathway into an external-capture final state.
struct ExternalRow {
  char rad_type;              ///< 'E' or 'M'
  int mult, l;                ///< multipolarity L; entrance orbital angular momentum l
  double J;                   ///< compound J
  int final_group, final_level, final_channel;   ///< the final state's J-group, level and channel (1-based)
  bool channel_capture;
  int cc_exit, cc_group, cc_pathway;  ///< channel capture: the internal pathway it passes through, as positions (exit,
                                      ///< k, m; 1-based); 0 for direct capture
  int cc_row;                 ///< that pathway as an index into EntrancePathways::internal; -1 for direct capture
  int photon_channel;         ///< J-group channel of the block's photon channel to the final state (1-based), 0 if none
};

/// The internal and external pathways of one pair of channel spins (s, s', and s2' for an unobserved primary).
struct SpinGroupRow {
  double s, sp, sp2;          ///< channel spins; sp2 = 0.0 unless has_sp2
  bool has_sp2;
  int first_internal, num_internal, first_external, num_external;
  /// A UPOS group with s' != s2': it holds the pathways of both exit spins and feeds only the secondary-gamma terms.
  bool CrossSpin() const { return has_sp2 && sp != sp2; }
};

/// The pathways from one entrance to one exit pair.
struct ExitRow {
  int pair;                   ///< exit pair number (1-based)
  int first_group, num_groups;
  bool has_interferences;     ///< false for a beta entrance, which has no angular distribution
  int first_interference, num_interferences, num_orders;  ///< num_orders: one more than the largest Legendre order
  int first_interference_pathway, num_interference_pathways;
  int first_diagonal, num_diagonal;
  bool has_capture_ay;        ///< the capture analyzing-power terms were requested for this exit
  int first_capture_ay, num_capture_ay;
};

/// Every pathway of one entrance pair, in flat arrays linked by offsets and counts.
struct EntrancePathways {
  std::vector<ExitRow> exits;
  std::vector<SpinGroupRow> groups;
  std::vector<InternalRow> internal;
  std::vector<ExternalRow> external;
  std::vector<FlatInterference> interferences;
  std::vector<InterferencePathway> interference_pathways;
  std::vector<DiagonalPathway> diagonal;
  std::vector<CaptureAyTerm> capture_ay;
  /// 1-based position of the exit to this pair (1-based), 0 if the pair is no exit of this entrance.
  int ExitOf(int exitPair) const;
  /// Spin group k (1-based) of exit e (1-based).
  const SpinGroupRow &Group(const ExitRow &e, int k) const { return groups[e.first_group + k - 1]; }
  /// Internal pathway m (1-based) of group g.
  const InternalRow &Internal(const SpinGroupRow &g, int m) const { return internal[g.first_internal + m - 1]; }
  /// External pathway m (1-based) of group g.
  const ExternalRow &External(const SpinGroupRow &g, int m) const { return external[g.first_external + m - 1]; }
};

/// An external-capture final state for BuildPathwayTable: what reaction_pathways needs, and where it is in CNuc.
struct PathwayFinalState {
  physics::FinalState state;
  int group = 0, level = 0;   ///< J-group and level (1-based)
};

/// What the table depends on besides the System and the final states (all fixed before CNuc::Initialize).
struct PathwayTableOptions {
  bool rmatrix_capture = false;   ///< the R-matrix-capture formalism: no photon exit pathways
  int max_l_order = 0;            ///< the largest Legendre order of the interferences and capture-A_y terms
  std::vector<int> entrances;     ///< the entrance pairs (1-based)
  /// The (entrance, exit) pairs (1-based) whose capture analyzing-power terms are built.
  std::vector<std::pair<int, int> > capture_ay;
};

class PathwayTable {
 public:
  /// The pathways from pair `pair` (1-based); empty for a pair that is not an entrance.
  const EntrancePathways &Entrance(int pair) const { return entrances_[pair - 1]; }
  int NumPairs() const { return (int)entrances_.size(); }

 private:
  friend std::shared_ptr<const PathwayTable> BuildPathwayTable(const CNucSystem &,
                                                               const std::vector<PathwayFinalState> &,
                                                               const PathwayTableOptions &);
  std::vector<EntrancePathways> entrances_;  ///< one per pair, in CNuc order
};

/*!
 * The table of a run: reaction_pathways for each entrance, the rows in its order with the engine numbers through the
 * System map, then for each exit of a non-beta entrance the interferences of the spin groups for L = 0 ...
 * max_l_order (particle_interference, photon_interference, upos_interference, kept when |Z| > 1e-10), the diagonal
 * pathways, and the requested capture-A_y terms (capture_ay_coefficients, kept when |a| or |b| >= 1e-10).
 * `structure` is BuildSystemFromCNuc(..., CNucSystemContent::Structure).
 */
std::shared_ptr<const PathwayTable> BuildPathwayTable(const CNucSystem &structure,
                                                      const std::vector<PathwayFinalState> &finalStates,
                                                      const PathwayTableOptions &options);

}  // namespace azure2

#endif
