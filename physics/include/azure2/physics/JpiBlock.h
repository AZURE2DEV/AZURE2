#ifndef AZURE2_PHYSICS_JPIBLOCK_H
#define AZURE2_PHYSICS_JPIBLOCK_H

#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {

/*!
 * One J^pi block of the level matrix: the levels lambda of one total angular momentum J and parity pi of one
 * compound nucleus, and the channels c they decay to.  The level matrix of a project is block-diagonal in (J, pi,
 * nucleus); each block is solved on its own (docs/upgrade/ARCHITECTURE.md section 3, docs/upgrade/physics_system.md).
 *
 * Contents, all read-only after construction:
 *
 *   two_J        2J (integer; J = two_J / 2)
 *   parity       pi = +1 or -1
 *   nucleus      index of the compound nucleus in System::nuclei()
 *   channels[k]  index c into System::channels() of the block's k-th channel, k = 0 .. n_c - 1
 *   B[k]         boundary value B_c of that channel in this block (dimensionless)
 *   E[lambda]    level energy E_lambda (MeV), the excitation energy E_x of the compound nucleus, lambda = 0 .. n_l - 1
 *   gamma        reduced width amplitudes gamma_{lambda c}, row-major: gamma[lambda * n_c + k]
 *
 * Choices not fixed by ARCHITECTURE.md, and why:
 *
 *  - Angular momenta are stored doubled (two_J), as integers, so that half-integer values are exact and the
 *    angular-coupling code to come (6j and Clebsch-Gordan coefficients) needs no floating-point comparison.
 *  - The boundary value B_c lives in the block, not in System::channels(): today's engine (CompoundSystem::
 *    CalcBoundaryConditions) sets B_c = S_c(E_1) at the first level of each J^pi group, so the same channel
 *    (alpha, l, s) has a different B_c in each block.  A System channel holds what does not depend on J^pi.
 *  - gamma is one contiguous row-major array (level-major): one level's amplitudes to all channels of the block are
 *    adjacent, which is the order in which the level matrix and the reduced-width-amplitude vectors are built.
 *  - Indices are 0-based throughout physics/ (the engine's CompoundSystem counts from 1).
 *  - What E and gamma mean (formal, Brune or Park parameters) is stated once, by System::basis().
 *
 * A block may have no levels (n_l = 0): a hard-sphere block, whose channels scatter through their hard-sphere phases
 * only.  It may also have no channels.
 *
 * The constructor checks the shapes and throws std::invalid_argument on a mismatch; it does not check the indices
 * against a System (System's constructor does).
 */
class JpiBlock {
 public:
  JpiBlock(int two_J, int parity, int nucleus, std::vector<int> channels, std::vector<double> boundary,
           std::vector<double> energies_MeV, std::vector<double> gamma);

  int two_J() const { return two_J_; }
  int parity() const { return parity_; }
  int nucleus() const { return nucleus_; }

  std::size_t num_levels() const { return E_.size(); }
  std::size_t num_channels() const { return channels_.size(); }

  /// Index into System::channels() of the block's k-th channel.
  int channel(std::size_t k) const { return channels_[k]; }
  const std::vector<int> &channels() const { return channels_; }

  /// Boundary value B_c of the block's k-th channel (dimensionless).
  double B(std::size_t k) const { return B_[k]; }
  const std::vector<double> &B() const { return B_; }

  /// Level energy E_lambda (MeV, excitation energy of the compound nucleus).
  double E(std::size_t lambda) const { return E_[lambda]; }
  const std::vector<double> &E() const { return E_; }

  /// Reduced width amplitude gamma_{lambda c} of level lambda in the block's k-th channel.
  double gamma(std::size_t lambda, std::size_t k) const { return gamma_[lambda * channels_.size() + k]; }
  /// All amplitudes, row-major: gamma()[lambda * num_channels() + k].
  const std::vector<double> &gamma() const { return gamma_; }

  /// Position k of System channel c in this block, or -1 if the block does not contain it.
  int find_channel(int c) const;

 private:
  int two_J_;
  int parity_;
  int nucleus_;
  std::vector<int> channels_;
  std::vector<double> B_;
  std::vector<double> E_;
  std::vector<double> gamma_;
};

}  // namespace physics
}  // namespace azure2

#endif
