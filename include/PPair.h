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
  /// Mass of particle 1 or 2, in u.
  double GetM(int) const;
  /// g-factor of particle 1 or 2.
  double GetG(int) const;
  /// Intrinsic spin of particle 1 or 2.
  double GetJ(int) const;
  /// Excitation energy of the residual nucleus (MeV); nonzero for an excited-state pair.
  double GetExE() const;
  /// Separation energy of the pair (MeV). Subtract it, and the excitation
  /// energy, from a level's energy to get that channel's centre-of-mass energy.
  double GetSepE() const;
  /// Channel radius (fm) -- the matching surface between internal and external regions.
  double GetChRad() const;
  /// Change the channel radius. Everything derived from it -- penetrabilities,
  /// shift functions, boundary conditions, external-capture integrals -- must
  /// be recomputed afterwards.
  void SetChRad(double);
  /// THM binding energy B_xs (MeV) of the transferred particle in the
  /// Trojan-Horse nucleus; enters the half-off-shell momentum of the
  /// entrance transfer form factor. 0 for a conventional pair.
  double GetBindingEnergy() const;
  void SetBindingEnergy(double);
  /// Reduced mass of the pair, in u.
  double GetRedMass() const;
  /// \f$1/[(2I_1+1)(2I_2+1)]\f$, the spin-statistics denominator of the entrance channel.
  double GetI1I2Factor() const;
  /// Returns true if the two particles in the pair are identical
  /// (same Z, mass, intrinsic spin, parity, and excitation energy).
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
  /// \f$1+\delta_{12}\f$: 2 for an identical pair, 1 otherwise.  An identical
  /// pair in the ENTRANCE channel multiplies every cross section out of it --
  /// reaction, capture or elastic, angle-integrated or differential -- by this
  /// factor, sigma = (pi/k^2) sum_J g_J (1+delta_12) sum |T|^2 over the channels
  /// symmetry allows -- for spin-0 bosons the familiar
  /// sigma = (pi/k^2) sum_l (2l+1) [1+(-1)^l] T_l (e.g. Hagino, Rowley and Yao,
  /// EPJ Web Conf. 122 (2016) 07002).  It is the convention of the reciprocity
  /// theorem w_12 k_12^2 sigma_12->34 / (1+delta_12) =
  /// w_34 k_34^2 sigma_34->12 / (1+delta_34), w = (2i_1+1)(2i_2+1)
  /// (Blatt and Weisskopf 1952; Fowler,
  /// Caughlan and Zimmerman, ARA&A 5 (1967) 525) with a symmetric U.  An
  /// identical pair in the EXIT channel contributes no factor: the cross section
  /// counts each reaction once, not each outgoing particle.  The elastic
  /// differential cross section, which by convention counts particles, carries
  /// a further 2 of its own (GenMatrixFunc).
  double GetEntranceSymmetryFactor() const;
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
  double pair_m_[2];
  double pair_g_[2];
  double pair_j_[2];
  double pair_ex_e_;
  double pair_sep_e_;
  double pair_ch_rad_;
  double binding_energy_;
  double red_mass_;
  double i1i2factor_;
  std::vector<Decay> decays_;
};

#endif
