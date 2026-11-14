#include "AChannel.h"
#include "azure2/io/LevelRows.h"
#include "Constants.h"
#include "azure2/physics/Widths.h"
#include <math.h>
#include <assert.h>


/*!
 * This constructor can be used if a channel is to be created from a specific line of the nuclear input file.
 */

AChannel::AChannel(const azure2::io::LevelRow &nucLine, int pairNum) {
  l_ = nucLine.l();
  s_ = nucLine.s();
  pair_ = pairNum;
  wigner_limit_ = 0.0;  // Initialize to zero, will be calculated later when PPair info is available
  if (nucLine.pType() == 0) {
    radtype_ = 'P';
  } else if (nucLine.pType() == 10) {
    if (nucLine.levelPi() * nucLine.pi2() == pow(-1, nucLine.l()))
      radtype_ = 'E';
    else
      radtype_ = 'M';
  } else if (nucLine.pType() == 20) {
    if (nucLine.l() == 0)
      radtype_ = 'F';
    else
      radtype_ = 'G';
  }
}

/*!
 * This constructor can be used if a channel is to be created directly using specified channel couplings.
 */

AChannel::AChannel(int lValue, double sValue, int pairNum, char radType) {
  l_ = lValue;
  s_ = sValue;
  pair_ = pairNum;
  radtype_ = radType;
  wigner_limit_ = 0.0;  // Initialize to zero, will be calculated later when PPair info is available
};

/*!
 * Returns the pair number, or position of the corresponding particle pair in the PPair vector, for the channel.
 */

int AChannel::GetPairNum() const {
  return pair_;
}

/*!
 * Returns the orbital angular momentum for the channel.
 */

int AChannel::GetL() const {
  return l_;
}

/*!
 * Returns the coupled channel spin for the channel.
 */

double AChannel::GetS() const {
  return s_;
}

/*!
 * Returns the boundary condition for the channel.  The energy of the boundary condition is fixed at the first level given in the nuclear input file for a given \f$ J^\pi \f$ group.
 */

double AChannel::GetBoundaryCondition() const {
  return boundary_condition_;
}

/*!
 * Returns the radiation type for the channel. Radiation types are:
 * - P: Particle Radiation
 * - E: EL Electromagnetic Radiation
 * - M: ML Electromagnetic Radiation
 */

char AChannel::GetRadType() const {
  return radtype_;
}

/*!
 * This function is used to set the boundary condition for the channel.
 */
void AChannel::SetBoundaryCondition(double boundaryCondition) {
  boundary_condition_ = boundaryCondition;
}

/*!
 * This function is used to set the Wigner Limit for the channel.
 * The Wigner Limit (single-particle/sum-rule reduced width) is the standard
 * Lane & Thomas (1958) result: gamma_W^2 = 3 * hbar^2 / (2 * reduced_mass * channel_radius^2).
 */
void AChannel::SetWignerLimit(double reducedMass, double channelRadius) {
  // Wigner limit gamma_W^2 = 3 hbar^2 / (2 mu a^2) (MeV), reduced mass in u, radius in fm (physics/Widths.h).
  // No limit (0, skipped by CompoundSystem::FillMnParams) for a channel without a radius, and for a photon or beta
  // channel: a photon pair has mu = 0, so the formula gave infinity whenever its radius was nonzero, and
  // FillMnParams then set Minuit limits (-inf, inf), which its sine transform turns into NaN (audit
  // 2026-10-25, item 32).  The Wigner limit is a statement about particle channels.
  if (channelRadius == 0 || !(reducedMass > 0.0) || radtype_ != 'P') {
    wigner_limit_ = 0;
    return;
  }

  wigner_limit_ = azure2::physics::wigner_limit(reducedMass, channelRadius);  // physics/Widths.h
}

/*!
 * Returns the Wigner Limit for the channel in MeV.
 * The Wigner Limit represents the energy scale below which the R-Matrix formalism
 * becomes less reliable due to the finite channel radius.
 */
double AChannel::GetWignerLimit() const {
  return wigner_limit_;
}
