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
  PPair(NucLine);
  bool IsEntrance() const;
  int GetZ(int) const;
  int GetPi(int) const;
  int GetPType() const;
  int NumDecays() const;
  int IsDecay(Decay);
  int IsDecay(int);
  int GetPairKey() const;
  double GetM(int) const;
  double GetG(int) const;
  double GetJ(int) const;
  double GetExE() const;
  double GetSepE() const;
  double GetChRad() const;
  void SetChRad(double);
  /// THM binding energy B_xs (MeV) of the transferred particle in the
  /// Trojan-Horse nucleus; enters the half-off-shell momentum of the
  /// entrance transfer form factor. 0 for a conventional pair.
  double GetBindingEnergy() const;
  void SetBindingEnergy(double);
  double GetRedMass() const;
  double GetI1I2Factor() const;
  /// Returns true if the two particles in the pair are identical
  /// (same Z, mass, intrinsic spin, parity, and excitation energy).
  /// Detected automatically in the constructor.
  bool IsIdentical() const;
  /// Boson/fermion sign: +1 for identical bosons (2j even), -1 for
  /// identical fermions. Returns +1 when IsIdentical() is false.
  int GetIdenticalSign() const;
  void AddDecay(Decay);
  void SetEntrance();
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
