#ifndef PPAIR_H
#define PPAIR_H

#include "Decay.h"

class NucLine;

/// An AZURE Particle Pair

/*!
 * In R-Matrix theory, the configuration space in the external region is decomposed into combinations of particle
 * pairs, traditionally given by the symbol \f$ \alpha \f$.  In AZURE, these particle pair are represented by a PPair object.
 * PPair objects are containers for vectors of Decay objects.
 */

class PPair {
 public:
  /// Build from a line of the nuclear input file.
  PPair(NucLine);
  /// Is this the entrance pair the calculation is driven from?
  bool IsEntrance() const;
  /// Charge number of particle 1 or 2.
  int GetZ(int) const;
  /// Parity of particle 1 or 2, \f$\pm1\f$.
  int GetPi(int) const;
  /// Pair type: 0 particle+particle, 10 particle+gamma.
  int GetPType() const;
  /// Number of exit pairs reachable from here; nonzero only for an entrance pair.
  int NumDecays() const;
  /// 1-based position of the decay in the Decay vector, or 0 if absent.
  int IsDecay(Decay);
  /// 1-based position of the decay to that exit pair, or 0 if absent.
  int IsDecay(int);
  /// The pair's key as written in the input file -- its 1-based position, and
  /// what an AChannel's pair number refers to.
  int GetPairKey() const;
  /// Charge number of the compound nucleus this pair forms, inferred from the
  /// pair columns: Z1+Z2 for a particle pair, Z2 for a gamma pair, Z2-Z1 for a
  /// beta pair (Z1 = -1 for beta-minus).
  int GetNucleusZ() const;
  /// Mass number of the compound nucleus: round(M1+M2) for a particle pair,
  /// round(M2) otherwise.
  int GetNucleusA() const;
  /// 1-based index of the compound nucleus in the CNuc's nucleus list, set by
  /// CNuc::AddPair; 0 until the pair is added.
  int GetNucleus() const;
  void SetNucleus(int);
  /// Mass of particle 1 or 2, in u.
  double GetM(int) const;
  /// g-factor of particle 1 or 2.
  double GetG(int) const;
  /// Intrinsic spin of particle 1 or 2.
  double GetJ(int) const;
  /// Total excitation energy of the pair (MeV), e1 + e2 -- light plus heavy
  /// particle.  The channel threshold is GetSepE() + GetExE().
  double GetExE() const;
  /// Excitation energy of the light particle (MeV), the optional 32nd column
  /// of a <levels> line (e1; 0 when absent).  Always 0 for gamma / beta pairs.
  double GetLightExE() const;
  /// Excitation energy of the heavy particle (MeV), column 17 (e2); for a
  /// gamma pair the energy of the final level.
  double GetHeavyExE() const;
  /// Separation energy of the pair (MeV). Subtract it, and the excitation
  /// energy, from a level's energy to get that channel's centre-of-mass energy.
  double GetSepE() const;
  /// The pair's default channel radius (fm), column 28 of its <levels> lines -- the matching surface between
  /// internal and external regions for every channel of the pair that has no radius of its own.  A channel's
  /// radius is GetChRad(l, s, J) (or AChannel::GetRadius()); this one is for what belongs to the pair, not to a
  /// channel (a beta pair's Fermi-function radius, printouts, the API's pair records).
  double GetPairRadius() const;
  /// Change the default channel radius. Everything derived from it -- penetrabilities,
  /// shift functions, boundary conditions, external-capture integrals -- must
  /// be recomputed afterwards (CNuc::Fill does, on a rebuild).
  void SetChRad(double);
  /// A radius of one channel or a set of channels of this pair: orbital angular momentum l, and optionally the
  /// channel spin s and the total J (a negative value means any).  From the <channelRadii> block (2026-11-10).
  struct ChannelRadius {
    int l;
    double s;
    double J;
    double radius;
  };
  /// Give the channels (l, s, J) of this pair their own radius; s < 0 or J < 0 means any.  A second call for the
  /// same (l, s, J) replaces the first.
  void SetChannelRadius(int l, double s, double J, double radius);
  /// Radius (fm) of the pair's channel (l, s, J): the most specific matching SetChannelRadius entry (l, s and J
  /// > l and J > l and s > l alone), or the pair's default radius.  Every channel quantity -- penetrability,
  /// shift, hard-sphere phase, boundary condition, Wigner limit, ANC conversion -- is taken at this radius.
  double GetChRad(int l, double s, double J) const;
  /// True when channel (l, s, J) matches an "l and s" entry and an "l and J" entry but no entry naming all
  /// three: which of the two applies would be a matter of the precedence rule, so CNuc refuses it.
  bool ChannelRadiusAmbiguous(int l, double s, double J) const;
  const std::vector<ChannelRadius> &GetChannelRadii() const;
  /// Reduced mass of the pair, in u.
  double GetRedMass() const;
  /// \f$1/[(2I_1+1)(2I_2+1)]\f$, the spin-statistics denominator of the entrance channel.
  double GetI1I2Factor() const;
  /// Returns true if the two particles in the pair are identical
  /// (same Z, mass, intrinsic spin, parity, g-factor and excitation energy,
  /// e1 = e2), for a particle pair only.
  /// Detected automatically in the constructor.
  bool IsIdentical() const;
  /// Boson/fermion sign: +1 for identical bosons (2j even), -1 for
  /// identical fermions. Returns +1 when IsIdentical() is false.
  int GetIdenticalSign() const;
  /// True for an identical pair whose members carry spin (j != 0).  For these
  /// the exchange term of the elastic amplitude depends on the channel spin s
  /// through (-1)^s, so no single sign symmetrizes it -- the spin-0 shortcut
  /// (one Coulomb amplitude f_C(theta) + f_C(pi - theta) for every channel)
  /// does not apply.
  bool HasSpinDependentExchange() const;
  void AddDecay(Decay);
  void SetEntrance();
  /// Decay \p i, 1-based.
  Decay *GetDecay(int);

 private:
  bool entrance_;
  bool ec_entrance_;
  bool is_identical_;
  int identical_sign_;
  int pair_z_[2];
  int pair_pi_[2];
  int pair_ptype_;
  int pair_key_;
  int nucleus_z_;
  int nucleus_a_;
  int nucleus_;
  double pair_m_[2];
  double pair_g_[2];
  double pair_j_[2];
  double pair_light_ex_e_;
  double pair_heavy_ex_e_;
  double pair_sep_e_;
  double pair_ch_rad_;
  std::vector<ChannelRadius> channel_radii_;
  double red_mass_;
  double i1i2factor_;
  std::vector<Decay> decays_;
};

#endif
