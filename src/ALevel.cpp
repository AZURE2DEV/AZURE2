#include "ALevel.h"
#include "azure2/io/LevelRows.h"

/*!
 * This constructor is used when a level object is created from an entry in the nuclear file.
 */

ALevel::ALevel(const azure2::io::LevelRow &nucLine) :
  level_e_(nucLine.levelE()),
  isinrmatrix_(true),
  isECLevel_(false),
  ecMultMask_(0),
  ecPairNum_(0) {
  if (nucLine.levelFix() == 1)
    energyfixed_ = true;
  else
    energyfixed_ = false;
}

/*!
 * This constructor is used when a level object is created using a specific energy.
 */

ALevel::ALevel(double energy) :
  energyfixed_(true),
  level_e_(energy),
  isinrmatrix_(false),
  isECLevel_(false),
  ecMultMask_(0),
  ecPairNum_(0) {};

/*!
 * Returns true if the level energy is to be fixed in the fitting process, otherwise returns false.
 */

bool ALevel::EnergyFixed() const {
  return energyfixed_;
}

/*!
 * Returns true if the level is to be included in the A-/R-Matrix calculation, otherwise returns false.
 * A level may specify a bound state for external capture, but may not be an R-Matrix state (i.e. subthreshold state).
 */

bool ALevel::IsInRMatrix() const {
  return isinrmatrix_;
}

/*!
 * Returns true if the reduced width amplitude for corresponding channel number is to be fixed
 * in the fitting process, otherwise returns false.
 */

bool ALevel::ChannelFixed(int channelNum) const {
  return channelfixed_[channelNum - 1];
}

/*!
 * Returns true if the input width for the given channel was specified as a
 * reduced width amplitude (MeV^(1/2)) rather than a physical partial
 * width/ANC, so the parameter transformations must not convert it.
 */

bool ALevel::GammaIsRWA(int channelNum) const {
  return gamma_is_rwa_[channelNum - 1];
}

bool ALevel::FileZeroWidth(int channelNum) const {
  return zero_input_[channelNum - 1];
}

/*!
 * Returns true if the level is a final state for external capture, otherwise returns false.
 */

bool ALevel::IsECLevel() const {
  return isECLevel_;
}

/*!
 * Returns non-zero only if the level is a final state for external capture.
 */

int ALevel::NumNFIntegrals() const {
  return nf_integrals_.size();
}

/*!
 * Returns the position in the pairs vector corresponding the the external capture level.
 */

int ALevel::GetECPairNum() const {
  return ecPairNum_;
}

/*!
 * Returns the multipolarity mask of external capture gammas to the level.
 */

unsigned char ALevel::GetECMultMask() const {
  return ecMultMask_;
}

/*!
 * The level energy of its <levels> line.
 */

double ALevel::FileE() const {
  return level_e_;
}

/*!
 * The width of the channel in the <levels> lines.
 */

double ALevel::FileGamma(int channelNum) const {
  return gammas_[channelNum - 1];
}

/*!
 * Returns the calculated channel integral  in the denominator of the \f$N_f^{1/2} \f$ term for a given channel number.
 */

double ALevel::GetNFIntegral(int channelNum) const {
  return nf_integrals_[channelNum - 1];
}

/*!
 * Returns the conversion factor from reduced width amplitude to ANC for a given channel number.
 */


double ALevel::GetECConversionFactor(int channelNum) const {
  return ec_conv_factors_[channelNum - 1];
}

void ALevel::PadChannels(int n) {
  while ((int)gammas_.size() < n) {
    gammas_.push_back(0.0);
    channelfixed_.push_back(true);
    gamma_is_rwa_.push_back(false);  // a padded channel has a zero width: the flag is moot
    zero_input_.push_back(true);
  }
}

void ALevel::SetGammaFromLine(int channelNum, const azure2::io::LevelRow &nucLine) {
  this->PadChannels(channelNum);
  gammas_[channelNum - 1] = nucLine.gamma();
  channelfixed_[channelNum - 1] = (nucLine.channelFix() == 1);
  // The THM field 33 (gammaIsRWA) goes with its width to the group channel.
  gamma_is_rwa_[channelNum - 1] = (nucLine.gammaIsRWA() == 1);
  zero_input_[channelNum - 1] = (nucLine.gamma() == 0.0);
}

/*!
 * This function creates and fills a position for the channel integral in the
 * denominator of the \f$N_f^{1/2}\f$ term.  The integral is of the form
 * \f$ \int_a^\infty \left[ \frac{W_c(kr)}{W_c{ka_c}} \right]^2 \f$.
 */

void ALevel::AddNFIntegral(double integral) {
  nf_integrals_.push_back(integral);
}

/*!
 * This function adds a conversion factor from reduced width amplitude to ANC.
 */

void ALevel::AddECConversionFactor(double conversionFactor) {
  ec_conv_factors_.push_back(conversionFactor);
}

/*!
 * The energy (MeV, E_x) with which this level enters external capture as a final state: its Brune energy at setup
 * (CompoundSystem::SetUpECFinalBlocks), the input energy except in the formal basis.
 */

double ALevel::GetECEnergy() const {
  return ec_energy_;
}

void ALevel::SetECEnergy(double energy) {
  ec_energy_ = energy;
}

/*!
 * Sets the external capture parameters for the level.
 */

void ALevel::SetECParams(int pairNum, unsigned char multMask) {
  isECLevel_ = true;
  ecPairNum_ = pairNum;
  ecMultMask_ = multMask;
}

void ALevel::AddID(int id) {
  for (int known : ids_)
    if (known == id) return;
  ids_.push_back(id);
}
