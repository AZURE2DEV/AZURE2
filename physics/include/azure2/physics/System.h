#ifndef AZURE2_PHYSICS_SYSTEM_H
#define AZURE2_PHYSICS_SYSTEM_H

#include <array>
#include <cstddef>
#include <vector>

#include "azure2/physics/JpiBlock.h"

namespace azure2 {
namespace physics {

/*!
 * The reaction system an R-matrix calculation is made for: the compound nuclei, the partitions alpha (the pairs of
 * nuclei, or nucleus and photon, or nucleus and lepton, that the compound nucleus can break up into), the channels
 * c = (alpha, l, s) and the J^pi blocks of levels.  docs/upgrade/ARCHITECTURE.md section 3 and
 * docs/upgrade/physics_system.md give the physics; this header gives the storage.
 *
 * Plain values, read-only after construction: no object holds a pointer to another.  A channel names its partition
 * by index, a block names its channels by index.  All indices are 0-based.
 *
 * Choices not fixed by ARCHITECTURE.md, and why:
 *
 *  - Spins are stored doubled (two_j, two_s, two_J) as integers: half-integers are exact and coupling code needs no
 *    floating-point comparison.
 *  - The channel radius a_c is held by the channel, as ARCHITECTURE.md says, although today's engine has one radius
 *    per partition (every channel of a partition carries the partition's radius); a later l-dependent radius then
 *    changes no layout.
 *  - The boundary value B_c is held by the J^pi block, per block channel, because today's engine sets it per J^pi
 *    group (JpiBlock.h).
 *  - A channel's identity is (partition, kind, l, two_s): for a photon channel the kind (electric or magnetic) is
 *    fixed by the parities, so the same (alpha, L, s) can be an E-channel in one block and an M-channel in another.
 *  - The compound nuclei are listed (Z, A), and each partition and block names its nucleus: a project with several
 *    compound nuclei has a level matrix block-diagonal in the nucleus.
 *  - The level parameters are stored in the basis the project uses (basis()), not converted: conversions belong to
 *    level_parameters (ARCHITECTURE.md section 5), which does not exist yet.
 *
 * Not included yet (docs/upgrade/physics_system.md lists them): the Convention record, external-capture final states
 * and their ANCs, hybrid-potential (alpha, l, j) channels, the THM binding energy of a Trojan-horse partition.
 *
 * The constructor checks every index and the consistency of kinds and throws std::invalid_argument when one fails.
 */

/// What the level energies E_lambda and amplitudes gamma_{lambda c} of the blocks are.
enum class LevelBasis {
  Formal,  ///< formal R-matrix parameters (Lane and Thomas), with the blocks' boundary values B_c
  Brune,   ///< Brune's alternative parameters (PRC 66, 044611 (2002)): observed energies, B_c = S_c(E_lambda)
  Park     ///< Park's observed-amplitude parameters (PRC 104, 064612 (2021)), on Brune's level energies
};

/// A compound nucleus: charge number Z and mass number A.
struct Nucleus {
  int Z = 0;
  int A = 0;
};

/// What a partition is made of.
enum class PartitionKind {
  Particle,  ///< two nuclei (or nucleon and nucleus)
  Photon,    ///< photon + final state of the compound nucleus (radiative capture)
  Beta       ///< lepton + parent nucleus (beta-delayed particle emission)
};

/*!
 * A partition alpha: two particles, index 0 the light one, index 1 the heavy one, each in a definite state.
 * For a photon partition particle 0 is the photon and particle 1 the final nuclear state; for a beta partition
 * particle 0 is the lepton (Z = -1 for an electron) and particle 1 the parent nucleus.
 */
struct Partition {
  PartitionKind kind = PartitionKind::Particle;
  int nucleus = 0;                          ///< index of the compound nucleus in System::nuclei()
  std::array<int, 2> Z{{0, 0}};             ///< charge numbers
  std::array<double, 2> mass_u{{0., 0.}};   ///< atomic masses (u)
  std::array<int, 2> two_j{{0, 0}};         ///< twice the intrinsic spins
  std::array<int, 2> parity{{1, 1}};        ///< intrinsic parities, +1 or -1
  std::array<double, 2> g_factor{{0., 0.}}; ///< g-factors (magnetic moment / (j nuclear magnetons))
  std::array<double, 2> excitation_MeV{{0., 0.}};  ///< excitation energies e of the two particles (MeV)
  /// Separation energy S_alpha (MeV) of the ground-state pair from the compound nucleus; the channel threshold is
  /// S_alpha + e_0 + e_1 in excitation energy.  For a beta partition: the decay Q value Q_beta.
  double separation_MeV = 0.;
};

/// What a channel carries.
enum class ChannelKind {
  Particle,     ///< particle channel: l is the orbital angular momentum
  Electric,     ///< photon channel, electric multipole EL: l is L
  Magnetic,     ///< photon channel, magnetic multipole ML: l is L
  Fermi,        ///< beta channel, Fermi transition (l = 0)
  GamowTeller   ///< beta channel, Gamow-Teller transition (l = 1)
};

/// A channel c = (alpha, l, s): a partition, an orbital angular momentum (or multipolarity) and a channel spin.
struct Channel {
  int partition = 0;                  ///< index alpha into System::partitions()
  ChannelKind kind = ChannelKind::Particle;
  int l = 0;                          ///< orbital angular momentum l, or multipolarity L of a photon channel
  int two_s = 0;                      ///< twice the channel spin s
  double radius_fm = 0.;              ///< channel radius a_c (fm)
};

class System {
 public:
  System(LevelBasis basis, std::vector<Nucleus> nuclei, std::vector<Partition> partitions,
         std::vector<Channel> channels, std::vector<JpiBlock> blocks);

  LevelBasis basis() const { return basis_; }

  std::size_t num_nuclei() const { return nuclei_.size(); }
  std::size_t num_partitions() const { return partitions_.size(); }
  std::size_t num_channels() const { return channels_.size(); }
  std::size_t num_blocks() const { return blocks_.size(); }
  /// Total number of levels, summed over the blocks.
  std::size_t num_levels() const;

  const Nucleus &nucleus(std::size_t i) const { return nuclei_[i]; }
  const Partition &partition(std::size_t alpha) const { return partitions_[alpha]; }
  const Channel &channel(std::size_t c) const { return channels_[c]; }
  const JpiBlock &block(std::size_t b) const { return blocks_[b]; }

  const std::vector<Nucleus> &nuclei() const { return nuclei_; }
  const std::vector<Partition> &partitions() const { return partitions_; }
  const std::vector<Channel> &channels() const { return channels_; }
  const std::vector<JpiBlock> &blocks() const { return blocks_; }

  /// Index of the channel (partition, kind, l, two_s), or -1 if the System has none.
  int find_channel(int partition, ChannelKind kind, int l, int two_s) const;

 private:
  LevelBasis basis_;
  std::vector<Nucleus> nuclei_;
  std::vector<Partition> partitions_;
  std::vector<Channel> channels_;
  std::vector<JpiBlock> blocks_;
};

}  // namespace physics
}  // namespace azure2

#endif
