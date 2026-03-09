#ifndef PPAIR_H
#define PPAIR_H

#include "Decay.h"

class NucLine;

///An AZURE Particle Pair

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
  double GetRedMass() const;
  double GetI1I2Factor() const;
  double GetJFinal() const;
  int GetRadFinal() const;
  void AddDecay(Decay);
  void SetEntrance();
  Decay *GetDecay(int);
  bool IsUPOS();
  int GetSecondaryDecayL() const;
  double GetIc() const;
private:
  bool entrance_;
  bool ec_entrance_;
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
  double red_mass_;
  double i1i2factor_;
  double jfinal_;
  double radfinal_;
  int pair_is_upos_;
  int pair_secondary_decay_L;
  double pair_Ic;
  std::vector<Decay> decays_;
};

#endif
