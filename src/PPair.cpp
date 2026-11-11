#include "NucLine.h"
#include "PPair.h"
#include <assert.h>
#include <cmath>

/*!
 * A particle pair object is created from and entry in the nuclear input file.
 */

PPair::PPair(NucLine nucLine) {
  pair_z_[0] = nucLine.z1();
  pair_z_[1] = nucLine.z2();
  pair_m_[0] = nucLine.m1();
  pair_m_[1] = nucLine.m2();
  pair_pi_[0] = nucLine.pi1();
  pair_pi_[1] = nucLine.pi2();
  pair_g_[0] = nucLine.g1();
  pair_g_[1] = nucLine.g2();
  pair_j_[0] = nucLine.j1();
  pair_j_[1] = nucLine.j2();
  // e1 (light) and e2 (heavy) are kept apart; GetExE() is their sum, the
  // channel's excitation above the ground-state threshold.
  pair_light_ex_e_ = nucLine.e1();
  pair_heavy_ex_e_ = nucLine.e2();
  pair_sep_e_ = nucLine.sepE();
  pair_ch_rad_ = nucLine.chRad();
  pair_ptype_ = nucLine.pType();
  pair_key_ = nucLine.ir();
  red_mass_ = nucLine.m1() * nucLine.m2() / (nucLine.m1() + nucLine.m2());
  i1i2factor_ = 1.0 / (2. * nucLine.j1() + 1.0) / (2. * nucLine.j2() + 1.0);
  entrance_ = false;
  ec_entrance_ = false;
  nucleus_ = 0;

  // The compound nucleus the pair forms.  A particle pair is light + heavy
  // (Z1+Z2, M1+M2).  A gamma pair (pType 10) is the nucleus itself with a
  // photon (M1 = 0).  A beta pair (pType 20) is the parent (Z2, M2) emitting
  // an electron of charge Z1 = -1 (or a positron, +1), so the daughter the
  // levels belong to has Z2 - Z1 and the same mass number.
  if (pair_ptype_ == 20) {
    nucleus_z_ = pair_z_[1] - pair_z_[0];
    nucleus_a_ = (int)std::floor(pair_m_[1] + 0.5);
  } else if (pair_ptype_ == 10) {
    nucleus_z_ = pair_z_[1];
    nucleus_a_ = (int)std::floor(pair_m_[1] + 0.5);
  } else {
    nucleus_z_ = pair_z_[0] + pair_z_[1];
    nucleus_a_ = (int)std::floor(pair_m_[0] + pair_m_[1] + 0.5);
  }

  // Auto-detect identical-particle pair.  Two particles are identical when
  // every per-particle input matches: Z, mass, intrinsic spin, parity,
  // g-factor and excitation energy (e1 of the light particle, e2 of the
  // heavy one).  Until 2026-11-05 the light particle had no excitation column
  // and the test required e2 = 0, so a mutual-excitation channel such as
  // 12C(2+)+12C(2+) could not be entered as identical; it is now e1 = e2 =
  // 4.439 MeV with j1 = j2 = 2.  The boson/fermion sign is +1 for integer 2j
  // (bosons) and -1 for half-integer 2j (fermions).
  is_identical_ = (pair_z_[0] == pair_z_[1]) &&
      (std::fabs(pair_m_[0] - pair_m_[1]) < 1.0e-6) &&
      (std::fabs(pair_j_[0] - pair_j_[1]) < 1.0e-6) &&
      (pair_pi_[0] == pair_pi_[1]) &&
      (std::fabs(pair_g_[0] - pair_g_[1]) < 1.0e-9) &&
      (std::fabs(pair_light_ex_e_ - pair_heavy_ex_e_) < 1.0e-9) &&
      (pair_ptype_ == 0);
  if (is_identical_) {
    int two_j = static_cast<int>(2.0 * pair_j_[0] + 0.5);
    identical_sign_ = (two_j % 2 == 0) ? +1 : -1;
  } else {
    identical_sign_ = +1;
  }
}

/*!
 * Returns true if the particle pair is an internal entrance pair, otherwise returns false.
 */

bool PPair::IsEntrance() const {
  return entrance_;
}

/*!
 * Returns the atomic number of the specified particle (1 or 2).
 */

int PPair::GetZ(int particle) const {
  return pair_z_[particle - 1];
}

/*!
 * Returns the parity of the specified particle (1 or 2).
 */

int PPair::GetPi(int particle) const {
  return pair_pi_[particle - 1];
}

/*!
 * Returns the integer particle pair type.  Pair types currently used in AZURE are 0: particle,particle and 10: particle,gamma.
 */

int PPair::GetPType() const {
  return pair_ptype_;
}

/*!
 * Returns the number of decay particle pairs for a given pair.  Size of Decay vector will only be nonzero if PPair object is
 * an entrance pair.
 */

int PPair::NumDecays() const {
  return decays_.size();
}

/*!
 * Tests a given decay particle pair to determine if it is in the Decay vector.  If the decay particle pair exists
 * in the vector, the position in the vector is returned.  Otherwise, the function returns 0.
 */

int PPair::IsDecay(Decay decay) {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumDecays()) {
    if (decay.GetPairNum() == this->GetDecay(c + 1)->GetPairNum()) b = true;
    c++;
  }
  if (b)
    return c;
  else
    return 0;
}

/*!
 * Tests a given particle pair number to determine if there exists a corresponding particle pair decay in the Decay vector.
 * If the object exists, the position in the vector is returned.  Otherwise, the function returns 0.
 */

int PPair::IsDecay(int pairNum) {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumDecays()) {
    if (pairNum == this->GetDecay(c + 1)->GetPairNum()) b = true;
    c++;
  }
  if (b)
    return c;
  else
    return 0;
}

/*!
 * Returns the pair key for the particle pair.
 */

int PPair::GetPairKey() const {
  return pair_key_;
}

int PPair::GetNucleusZ() const {
  return nucleus_z_;
}

int PPair::GetNucleusA() const {
  return nucleus_a_;
}

int PPair::GetNucleus() const {
  return nucleus_;
}

void PPair::SetNucleus(int nucleus) {
  nucleus_ = nucleus;
}

/*!
 * Returns the mass number of the specified particle (1 or 2).
 */

double PPair::GetM(int particle) const {
  return pair_m_[particle - 1];
}

/*!
 * Returns the g-factor of the specified particle (1 or 2).
 */

double PPair::GetG(int particle) const {
  return pair_g_[particle - 1];
}

/*!
 * Returns the total spin of the specified particle (1 or 2).
 */

double PPair::GetJ(int particle) const {
  return pair_j_[particle - 1];
}

/*!
 * Returns the total excitation energy of the particle pair, e1 + e2: the
 * amount by which the channel threshold lies above the ground-state
 * separation energy.  Use it for anything about the channel's energy
 * (thresholds, Q-values, kinematics); GetHeavyExE() / GetLightExE() for
 * anything about one particle's own state.
 */

double PPair::GetExE() const {
  return pair_light_ex_e_ + pair_heavy_ex_e_;
}

/*!
 * Returns the excitation energy of the light particle (column 32, e1).
 */

double PPair::GetLightExE() const {
  return pair_light_ex_e_;
}

/*!
 * Returns the excitation energy of the heavy particle (column 17, e2); for a
 * gamma pair, the energy of the final level the photon feeds.
 */

double PPair::GetHeavyExE() const {
  return pair_heavy_ex_e_;
}

/*!
 * Returns the seperation energy of the particle pair.
 */

double PPair::GetSepE() const {
  return pair_sep_e_;
}

/*!
 * Returns the default channel radius of the particle pair (column 28 of its <levels> lines).
 */

double PPair::GetPairRadius() const {
  return pair_ch_rad_;
}

namespace {
bool SameSpin(double a, double b) { return fabs(a - b) < 1.0e-6; }
// 0: does not match; otherwise 1 + 1 (s named) + 2 (J named): l alone 1, l+s 2, l+J 3, l+s+J 4.
int MatchRank(const PPair::ChannelRadius &e, int l, double s, double J) {
  if (e.l != l) return 0;
  if (e.s >= 0.0 && !SameSpin(e.s, s)) return 0;
  if (e.J >= 0.0 && !SameSpin(e.J, J)) return 0;
  return 1 + (e.s >= 0.0 ? 1 : 0) + (e.J >= 0.0 ? 2 : 0);
}
}  // namespace

void PPair::SetChannelRadius(int l, double s, double J, double radius) {
  if (s < 0.0) s = -1.0;
  if (J < 0.0) J = -1.0;
  for (ChannelRadius &e : channel_radii_)
    if (e.l == l && SameSpin(e.s, s) && SameSpin(e.J, J)) {
      e.radius = radius;
      return;
    }
  channel_radii_.push_back(ChannelRadius{l, s, J, radius});
}

double PPair::GetChRad(int l, double s, double J) const {
  int best = 0;
  double radius = pair_ch_rad_;
  for (const ChannelRadius &e : channel_radii_) {
    const int rank = MatchRank(e, l, s, J);
    if (rank > best) {
      best = rank;
      radius = e.radius;
    }
  }
  return radius;
}

bool PPair::ChannelRadiusAmbiguous(int l, double s, double J) const {
  bool sOnly = false, jOnly = false;
  for (const ChannelRadius &e : channel_radii_) {
    const int rank = MatchRank(e, l, s, J);
    if (rank == 4) return false;
    if (rank == 2) sOnly = true;
    if (rank == 3) jOnly = true;
  }
  return sOnly && jOnly;
}

const std::vector<PPair::ChannelRadius> &PPair::GetChannelRadii() const {
  return channel_radii_;
}

/*!
 * Sets the channel radius of the particle pair.
 */

void PPair::SetChRad(double radius) {
  pair_ch_rad_ = radius;
}

/*!
 * Returns the reduced mass of the particle pair.
 */

double PPair::GetRedMass() const {
  return red_mass_;
}

/*!
 * Returns the factor \f$ \frac{1}{(2I_1+1)(2I_2+1)} \f$ of the particle pair.
 */

double PPair::GetI1I2Factor() const {
  return i1i2factor_;
}

bool PPair::IsIdentical() const {
  return is_identical_;
}

int PPair::GetIdenticalSign() const {
  return identical_sign_;
}

bool PPair::HasSpinDependentExchange() const {
  return is_identical_ && pair_j_[0] > 1.0e-6;
}

/*!
 * Adds a decay particle pair to the Decay vector.
 */

void PPair::AddDecay(Decay decay) {
  decays_.push_back(decay);
}

/*!
 * Sets the particle pair to be an internal entrance pair.
 */

void PPair::SetEntrance() {
  entrance_ = true;
}

/*!
 * Returns a pointer to the decay particle pair specified by a position in the Decay vector.
 */

Decay *PPair::GetDecay(int decayNum) {
  Decay *b = &decays_[decayNum - 1];
  return b;
}
