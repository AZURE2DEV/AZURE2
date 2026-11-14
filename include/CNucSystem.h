#ifndef CNUCSYSTEM_H
#define CNUCSYSTEM_H

#include <vector>

#include "azure2/physics/System.h"

class CNuc;
class Config;

/*!
 * The physics/ System of a loaded CNuc (docs/upgrade/physics_system.md), and where each of its entries came from.
 *
 * Indices in the map are CNuc's own, 1-based:
 *   pairOfPartition[alpha]          CNuc pair number (GetPair) of partition alpha
 *   groupOfBlock[b]                 CNuc J-group number (GetJGroup) of block b
 *   levelOfBlockLevel[b][lambda]    level number within that J-group (JGroup::GetLevel) of block level lambda
 *   channelOfBlockChannel[b][k]     channel number within that J-group (JGroup::GetChannel) of block channel k;
 *                                   the System channel itself is system.block(b).channel(k)
 */
struct CNucSystemMap {
  std::vector<int> pairOfPartition;
  std::vector<int> groupOfBlock;
  std::vector<std::vector<int> > levelOfBlockLevel;
  std::vector<std::vector<int> > channelOfBlockChannel;
};

struct CNucSystem {
  azure2::physics::System system;
  CNucSystemMap map;
};

/// What BuildSystemFromCNuc puts in the blocks besides their channels.
enum class CNucSystemContent {
  /// The level parameters and B_c the CNuc holds; the CNuc must have been initialized (below).
  Parameters,
  /// The structure only: every block with its channels, B_c = 0 and no levels; the map still lists the block's
  /// R-matrix levels.  Any state after CNuc::Fill.  This is CNuc::SystemStructure, what CNuc::CalcBoundaryConditions
  /// and CNuc::TransformIn hand to physics/ level_parameters, which takes the level values as separate arguments.
  Structure,
  /// As Structure, but with every J-group and level, also those outside the R-matrix (a group made only to hold an
  /// external-capture final state): the observed quantities (CNuc::ObservedWidths) are reported for all of them.
  /// Not for the R-matrix conversions, which see only R-matrix blocks and levels.
  AllGroups
};

/*!
 * Builds the System from a CNuc after CNuc::Fill and CNuc::Initialize (boundary conditions set and, unless
 * --no-transform, the input parameters transformed in); with CNucSystemContent::Structure, after CNuc::Fill alone.
 * Reads only.  The engine builds the structure and all groups once, at the end of CNuc::Fill (CNuc::SystemStructure,
 * CNuc::SystemAllGroups), and every reader takes those; the full System is used by the tests.
 *
 * What is taken, and from where:
 *  - nuclei: CNuc::GetNucleus (Z, A), in CNuc order;
 *  - partitions: every PPair in CNuc order (Z, M, J, pi, g, e_light, e_heavy, S); kind from the pair type
 *    (0 particle, 10 photon, 20 beta);
 *  - channels: the distinct (pair, kind, l, s) of the J^pi groups that take part in the R-matrix, in order of first
 *    appearance (group by group, channel by channel); a_c = the pair's channel radius;
 *  - blocks: every J^pi group with IsInRMatrix(), in CNuc order; its channels in group order with their B_c; its
 *    levels with IsInRMatrix(), in group order, with E = ALevel::GetE() and gamma = ALevel::GetGamma(ch): the
 *    parameters the project was loaded with, in the basis AZUREParams::BasisForConfig names (after TransformIn these
 *    are the formal/Brune/Park amplitudes in MeV^1/2, not the widths in eV of the file).  The fitted values
 *    (GetFitE, GetFitGamma), which move during a fit, are not read.
 *
 * Throws std::runtime_error when the CNuc is not in that state (CNucSystemContent::Parameters and the parameters
 * still the file's widths), when a spin is not a multiple of 1/2, or for a pair type or radiation type the System has
 * no kind for; std::invalid_argument from the System's own checks.
 */
CNucSystem BuildSystemFromCNuc(const CNuc &compound, const Config &configure,
                               CNucSystemContent content = CNucSystemContent::Parameters);

#endif
