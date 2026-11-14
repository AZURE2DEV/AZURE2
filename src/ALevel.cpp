#include "ALevel.h"
#include "azure2/io/LevelRows.h"
#include <atomic>

std::uint64_t ALevel::NewStamp() {
  static std::atomic<std::uint64_t> counter{0};
  return ++counter;
}

/*!
 * This constructor is used when a level object is created from an entry in the nuclear file.
 */

ALevel::ALevel(const azure2::io::LevelRow &nucLine) :
  level_e_(nucLine.levelE()),
  fitlevel_e_(0.0),
  isinrmatrix_(true),
  sqrt_nf_factor_(1.0),
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
  fitlevel_e_(0.0),
  isinrmatrix_(false),
  sqrt_nf_factor_(1.0),
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

bool ALevel::ZeroInputWidth(int channelNum) const {
  return zero_input_[channelNum - 1];
}

void ALevel::SetZeroInputWidth(int channelNum, bool zero) {
  if (zero_input_[channelNum - 1] == zero) return;
  stamp_ = NewStamp();
  zero_input_[channelNum - 1] = zero;
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
 *  Returns the energy of the level.
 */

double ALevel::GetE() const {
  return level_e_;
}

/*!
 * Returns the internal reduced width amplitude for a given channel number.
 */

double ALevel::GetGamma(int channelNum) const {
  return gammas_[channelNum - 1];
}

/*!
 * Returns the fitted internal reduced width amplitude for a given channel number.
 */

double ALevel::GetFitGamma(int channelNum) const {
  return fitgammas_[channelNum - 1];
}

/*!
 * Returns the fitted energy of the level.
 */

double ALevel::GetFitE() const {
  return fitlevel_e_;
}

/*!
 * Returns the calculated channel integral  in the denominator of the \f$N_f^{1/2} \f$ term for a given channel number.
 */

double ALevel::GetNFIntegral(int channelNum) const {
  return nf_integrals_[channelNum - 1];
}

/*!
 * Returns the \f$N_f^{1/2}\f$ term for the level.
 */

double ALevel::GetSqrtNFFactor() const {
  return sqrt_nf_factor_;
}

/*!
 * Returns the conversion factor from reduced width amplitude to ANC for a given channel number.
 */


double ALevel::GetECConversionFactor(int channelNum) const {
  return ec_conv_factors_[channelNum - 1];
}

/*!
 * Returns the Shift function for the specified channel number calculated at the resonance energy.
 */

double ALevel::GetShiftFunction(int channelNum) const {
  return shifts_[channelNum - 1];
}

void ALevel::PadChannels(int n) {
  stamp_ = NewStamp();
  while ((int)gammas_.size() < n) {
    gammas_.push_back(0.0);
    fitgammas_.push_back(0.0);
    ec_gammas_.push_back(0.0);
    channelfixed_.push_back(true);
    gamma_is_rwa_.push_back(false);  // a padded channel has a zero width: the flag is moot
    zero_input_.push_back(true);
    shifts_.push_back(0.0);
    shift_derivs_.push_back(0.0);
    shift_derivs2_.push_back(0.0);
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
 * This function sets the internal reduced width amplitude for a given channel number.
 */

void ALevel::SetGamma(int channelNum, double reducedWidth) {
  gammas_[channelNum - 1] = reducedWidth;
}

/*!
 * This function sets the level energy.
 */

void ALevel::SetE(double energy) {
  level_e_ = energy;
}

/*!
 * This function sets the fitted internal reduced width amplitude for a given channel number.
 */

void ALevel::SetFitGamma(int channelNum, double reducedWidth) {
  stamp_ = NewStamp();
  fitgammas_[channelNum - 1] = reducedWidth;
}

/*!
 * This function sets the fitted level energy.
 */

void ALevel::SetFitE(double energy) {
  stamp_ = NewStamp();
  fitlevel_e_ = energy;
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
 * This function sets the \f$N_f^{1/2}\f$ term for the level.
 */

void ALevel::SetSqrtNFFactor(double term) {
  sqrt_nf_factor_ = term;
}

/*!
 * This function adds a conversion factor from reduced width amplitude to ANC.
 */

void ALevel::AddECConversionFactor(double conversionFactor) {
  ec_conv_factors_.push_back(conversionFactor);
}

/*!
 * The amplitude with which this level enters external capture as a final state: its Brune (Park) amplitude, the
 * fit amplitude except in the formal basis (CNuc::FillCompoundFromParams; cnuc_split.md section 8).
 */

double ALevel::GetECGamma(int channelNum) const {
  return ec_gammas_[channelNum - 1];
}

void ALevel::SetECGamma(int channelNum, double reducedWidth) {
  ec_gammas_[channelNum - 1] = reducedWidth;
}

/*!
 * The energy (MeV, E_x) with which this level enters external capture as a final state: its Brune energy at setup
 * (CNuc::SetUpECFinalBlocks), the input energy except in the formal basis.
 */

double ALevel::GetECEnergy() const {
  return ec_energy_;
}

void ALevel::SetECEnergy(double energy) {
  ec_energy_ = energy;
}

/*!
 * Sets the value of the shift function calculated at the resonance energy.
 */

void ALevel::SetShiftFunction(int channelNum, double shiftFunction) {
  stamp_ = NewStamp();
  shifts_[channelNum - 1] = shiftFunction;
}

/*!
 * Returns the Park (2021) diagonal overlap \f$J_{\lambda\lambda}\f$ of the level.
 */

double ALevel::GetParkNorm() const {
  return park_norm_;
}

/*!
 * Sets the Park (2021) diagonal overlap \f$J_{\lambda\lambda}\f$ of the level.
 */

void ALevel::SetParkNorm(double parkNorm) {
  stamp_ = NewStamp();
  park_norm_ = parkNorm;
}

/*!
 * Returns dS_c/dE at the level energy (Park formalism only).
 */

double ALevel::GetShiftDerivative(int channelNum) const {
  return shift_derivs_[channelNum - 1];
}

/*!
 * Returns d^2S_c/dE^2 at the level energy (Park formalism only).
 */

double ALevel::GetShiftSecondDerivative(int channelNum) const {
  return shift_derivs2_[channelNum - 1];
}

void ALevel::SetShiftDerivative(int channelNum, double first) {
  shift_derivs_[channelNum - 1] = first;
}

/*!
 * Stores d^2S_c/dE^2; changes the stamp (ChannelSpace::CurrentLevelDerivatives
 * reads it).
 */

void ALevel::SetShiftSecondDerivative(int channelNum, double second) {
  stamp_ = NewStamp();
  shift_derivs2_[channelNum - 1] = second;
}

/*!
 * True when the stored shift functions were computed at exactly this level
 * energy.  They depend on nothing else that changes during a fit (the channel
 * radii and potentials are fixed once the compound nucleus is built).
 */

bool ALevel::ShiftFunctionsValidAt(double energy) const {
  return shift_cache_energy_ == energy;
}

void ALevel::SetShiftFunctionsEnergy(double energy) {
  shift_cache_energy_ = energy;
}

bool ALevel::ShiftSecondDerivativesValidAt(double energy) const {
  return shift_d2_energy_ == energy;
}

void ALevel::SetShiftSecondDerivativesEnergy(double energy) {
  shift_d2_energy_ = energy;
}

void ALevel::InvalidateShiftFunctions() {
  shift_cache_energy_ = std::numeric_limits<double>::quiet_NaN();
  shift_d2_energy_ = std::numeric_limits<double>::quiet_NaN();
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
