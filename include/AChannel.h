#ifndef ACHANNEL_H
#define ACHANNEL_H

class NucLine;

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
  AChannel(NucLine, int);
  /// Build directly from couplings: l, s, pair number, radiation type.
  AChannel(int, double, int, char);
  /// 1-based index of this channel's particle pair, as CNuc::GetPair takes it.
  int GetPairNum() const;
  /// Orbital angular momentum of the channel.
  int GetL() const;
  /// Channel spin, the coupling of the two particles' intrinsic spins.
  double GetS() const;
  /// Boundary condition \f$B_c\f$, fixed at the energy of the group's first level.
  double GetBoundaryCondition() const;
  /// Radiation type: 'P' particle, 'E'/'M' electric/magnetic multipole, 'F'/'G' beta decay.
  char GetRadType() const;
  void SetBoundaryCondition(double);
  /// Compute the Wigner limit from the pair's reduced mass (u) and channel radius (fm).
  /// Deferred until the PPair is known; a zero radius leaves the limit zero.
  void SetWignerLimit(double reducedMass, double channelRadius);
  /// Wigner limit on the *reduced width* \f$\gamma^2\f$ in MeV, i.e.
  /// \f$\hbar^2 / (m_{red} a^2)\f$ -- not a bound on the partial width.
  /// Zero until SetWignerLimit has been called.
  double GetWignerLimit() const;

 private:
  int l_;
  int pair_;
  double s_;
  char radtype_;
  double boundary_condition_;
  double wigner_limit_;
};

#endif
