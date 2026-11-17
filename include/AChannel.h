#ifndef ACHANNEL_H
#define ACHANNEL_H

namespace azure2 {
namespace io {
class LevelRow;
}
}  // namespace azure2

/// An AZURE channel object.

/*!
 * An R-Matrix channel for a given \f$J^\pi\f$ group represents a specfic combination of \f$ \alpha,s,l \f$ couplings.
 *
 * A channel is one way of reaching the group's \f$J^\pi\f$: a particle pair,
 * an orbital angular momentum and a channel spin. Every level in the group
 * carries a reduced width for every channel in it -- which is why the channel
 * set belongs to the group and not to the individual level.
 */

class AChannel {
 public:
  /// Build from a channel line of the nuclear input file; \p pairNum is the pair it decays to.
  AChannel(const azure2::io::LevelRow &, int);
  /// Build directly from couplings: l, s, pair number, radiation type.
  AChannel(int, double, int, char);
  /// 1-based index of this channel's particle pair, as CompoundSystem::GetPair takes it.
  int GetPairNum() const;
  /// Orbital angular momentum of the channel.
  int GetL() const;
  /// Channel spin, the coupling of the two particles' intrinsic spins.
  double GetS() const;
  /// Boundary condition \f$B_c\f$, fixed at the energy of the group's lowest level by the compound system's last
  /// part (CompoundSystem::CalcBoundaryConditions); std::logic_error before it set it.
  double GetBoundaryCondition() const;
  /// Radiation type: 'P' particle, 'E'/'M' electric/magnetic multipole, 'F'/'G' beta decay.
  char GetRadType() const;
  void SetBoundaryCondition(double);
  /// Set the channel's radius a_c (fm) and, from it and the pair's reduced mass (u), its Wigner limit.
  /// Deferred until the PPair is known; a zero radius leaves the limit zero.
  void SetRadius(double channelRadius, double reducedMass);
  /// The channel radius a_c (fm) at which every quantity of this channel is taken: penetrability, shift function
  /// and its derivatives, hard-sphere phase, boundary condition, Wigner limit, the bound-state integrals.  Set at
  /// set-up (CompoundSystem's build); std::logic_error before it.
  double GetRadius() const;
  /// Wigner limit on the *reduced width* \f$\gamma^2\f$ in MeV, i.e.
  /// \f$3\hbar^2 / (2 m_{red} a^2)\f$ (Teichmann-Wigner, with the 3/2 factor)
  /// -- not a bound on the partial width.
  /// Zero until SetRadius has been called.
  double GetWignerLimit() const;

 private:
  int l_;
  int pair_;
  double s_;
  char radtype_;
  double boundary_condition_;
  bool hasBoundaryCondition_ = false;  ///< set by SetBoundaryCondition; GetBoundaryCondition refuses before
  double wigner_limit_;
  double radius_ = 0.0;
  bool hasRadius_ = false;  ///< set by SetRadius; GetRadius refuses before
};

#endif
