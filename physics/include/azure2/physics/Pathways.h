#ifndef AZURE2_PHYSICS_PATHWAYS_H
#define AZURE2_PHYSICS_PATHWAYS_H

#include <cstddef>
#include <vector>

#include "azure2/physics/ExternalCapture.h"
#include "azure2/physics/System.h"

namespace azure2 {
namespace physics {

/*!
 * The reaction pathways of one entrance partition (docs/upgrade/cnuc_split.md section 10): for every exit partition,
 * the channel-spin groups and in each the internal pathways (level-matrix elements of one block) and the external
 * (direct- and channel-capture) pathways to the external-capture final states.  The lists are in the order the
 * engine's pathway trees (Decay, KGroup, MGroup, ECMGroup) hold them, and indices into them are 0-based.
 */

/// The multipoles a final state is captured into: bits as the engine's isE1, isM1, isE2.
constexpr unsigned char kCaptureE1 = 1 << 0;
constexpr unsigned char kCaptureM1 = 1 << 1;
constexpr unsigned char kCaptureE2 = 1 << 2;

/// One channel of an external-capture final state.  Only particle channels take part.
struct FinalStateChannel {
  int partition = 0;
  bool particle = true;
  int l = 0;
  int two_s = 0;
};

/// An external-capture final state: a level of a J^pi group that need not take part in the R-matrix.
struct FinalState {
  int two_J = 0;
  int parity = 1;
  int photon_partition = 0;           ///< the photon partition the capture populates
  unsigned char multipoles = 0;       ///< kCaptureE1 | kCaptureM1 | kCaptureE2
  std::vector<FinalStateChannel> channels;  ///< all channels of the state's group, in group order
};

/// An internal pathway: entrance channel k_in and exit channel k_out of block b (block channel indices).
struct InternalPathway {
  std::size_t block = 0;
  std::size_t k_in = 0;
  std::size_t k_out = 0;
  /// (2J + 1) / ((2 I_1' + 1)(2 I_2' + 1)) of the exit partition; set for every exit kind but a beta entrance.
  bool has_stat_spin_factor = false;
  double stat_spin_factor = 0.0;
};

/*!
 * An external pathway into final state f and its final channel c_f, through the compound block of an internal
 * pathway with exit partition = the final channel's partition (the channel-capture group).  Direct capture
 * (channel_capture false) is listed where that internal pathway is elastic in one channel (k_in = k_out); every
 * such internal pathway also gives a channel-capture entry.
 */
struct ExternalPathway {
  Multipole X = Multipole::E;
  int L = 1;
  int l = 0;                       ///< entrance orbital angular momentum (the internal pathway's entrance channel)
  int two_J = 0;                   ///< 2J of the compound block
  std::size_t final_state = 0;
  std::size_t final_channel = 0;   ///< index into FinalState::channels
  bool channel_capture = false;
  // The internal pathway it goes through (also for direct capture, whose entrance and exit channel it is):
  std::size_t cc_exit = 0;         ///< index of the exit partition entry of the internal pathway
  std::size_t cc_group = 0;        ///< its channel-spin group
  std::size_t cc_pathway = 0;      ///< its internal pathway
  int photon_channel = -1;         ///< block channel index of the block's photon channel (X, L) to the final state, or -1
};

/// The pathways of one pair of channel spins s (entrance) and s' (exit), and s2' for an unobserved primary.
struct SpinGroup {
  int two_s = 0;
  int two_sp = 0;
  bool has_sp2 = false;
  int two_sp2 = 0;
  std::vector<InternalPathway> internal;
  std::vector<ExternalPathway> external;
};

struct ExitPathways {
  std::size_t partition = 0;
  std::vector<SpinGroup> groups;
};

struct PathwayOptions {
  bool rmatrix_capture = false;  ///< the R-matrix-capture formalism: no photon exit pathways
};

/*!
 * The pathways from entrance partition alpha (an entrance of the project) of the System (R-matrix blocks only) and
 * the external-capture final states.  Exit kinds (as the engine's SortPathways):
 *  - beta entrance: for l = 0, 1, the pathways of entrance channels with that l to any non-beta exit; group (2l, 0);
 *  - elastic (particle, alpha' = alpha): groups (s, s');
 *  - particle exit alpha' != alpha: groups (s, s', s2'), a pathway in every group whose s' or s2' is its exit spin;
 *  - photon exit (unless options.rmatrix_capture): groups (s, 0).
 * External pathways follow for a non-beta entrance: into each final state with a photon decay, for each group of
 * that decay, through the internal pathways (deduplicated on (block, k_in, k_out)) to each particle channel of the
 * final state of the same entrance channel spin, for each multipole the parities, the mask, the J triangle and
 * couples() allow.
 */
std::vector<ExitPathways> reaction_pathways(const System &system, std::size_t alpha,
                                            const std::vector<FinalState> &final_states,
                                            const PathwayOptions &options);

}  // namespace physics
}  // namespace azure2

#endif
