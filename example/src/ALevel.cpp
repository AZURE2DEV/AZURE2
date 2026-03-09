#include "ALevel.h"
#include "NucLine.h"

/*!
 * This constructor is used when a level object is created from an entry in the nuclear file.
 */

ALevel::ALevel(NucLine nucLine) :
  level_e_(nucLine.levelE()),fitlevel_e_(0.0), isinrmatrix_(true), sqrt_nf_factor_(1.0), isECLevel_(false),
  ecMultMask_(0), ecPairNum_(0) {
  if(nucLine.levelFix()==1) energyfixed_=true;
  else energyfixed_=false;
}

/*!
 * This constructor is used when a level object is created using a specific energy.
 */

ALevel::ALevel(double energy) :
  energyfixed_(true),level_e_(energy),fitlevel_e_(0.0), isinrmatrix_(false), sqrt_nf_factor_(1.0), isECLevel_(false),
  ecMultMask_(0), ecPairNum_(0) {};

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
  return channelfixed_[channelNum-1];
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
 * Returns the number of iterations required to transform the level from formal to physical parameters.
 */

int ALevel::GetTransformIterations() const {
  return transform_iter_;
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
  return gammas_[channelNum-1];
}

/*!
 * Returns the fitted internal reduced width amplitude for a given channel number.
 */

double ALevel::GetFitGamma(int channelNum) const {
  return fitgammas_[channelNum-1];
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
  return nf_integrals_[channelNum-1];
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
  return ec_conv_factors_[channelNum-1];
}

/*!
 * Returns the physical internal reduced width amplitude for a given channel number.
 */

double ALevel::GetTransformGamma(int channelNum) const {
  return transform_gammas_[channelNum-1];
}

/*!
 * Returns the physical level energy.
 */

double ALevel::GetTransformE() const {
  return transform_e_;
}

/*!
 * Returns the Breit-Wigner partial width for a given channel number.
 */

double ALevel::GetBigGamma(int channelNum) const {
  return big_gammas_[channelNum-1];
}

/*!
 * Returns the Shift function for the specified channel number calculated at the resonance energy.
 */

double ALevel::GetShiftFunction(int channelNum) const {
  return shifts_[channelNum-1];
}

/*!
 * Returns the external portion of the reduced width amplitude for a given channel number.
 */

complex ALevel::GetExternalGamma(int channelNum) const {
  return external_gammas_[channelNum-1];
}

/*!
 * This function adds a position in the width vectors corresponding to a new channel.  
 * The initial reduced width amplitude is set from an entry in the nuclear input file.
 */

void ALevel::AddGamma(NucLine nucLine) {
  double b=nucLine.gamma();
  gammas_.push_back(b);
  fitgammas_.push_back(0.0);
  transform_gammas_.push_back(0.0);
  big_gammas_.push_back(0.0);
  external_gammas_.push_back(complex(0.0,0.0));
  if(nucLine.channelFix()==1) channelfixed_.push_back(true);
  else channelfixed_.push_back(false);
  shifts_.push_back(0.0);
}

/*!
 * This function adds a position in the width vectors corresponding to a new channel.  
 * The initial reduced width amplitude is set directly.
 */

void ALevel::AddGamma(double reducedWidth) {
  gammas_.push_back(reducedWidth);
  fitgammas_.push_back(0.0);
  transform_gammas_.push_back(0.0);
  big_gammas_.push_back(0.0);
  external_gammas_.push_back(complex(0.0,0.0));
  channelfixed_.push_back(false);
}

/*!
 * This function sets the internal reduced width amplitude for a given channel number.
 */

void ALevel::SetGamma(int channelNum, double reducedWidth) {
  gammas_[channelNum-1]=reducedWidth;
}

/*!
 * This function sets the level energy.
 */

void ALevel::SetE(double energy) {
  level_e_=energy;
}

/*!
 * This function sets the fitted internal reduced width amplitude for a given channel number.
 */

void ALevel::SetFitGamma(int channelNum,double reducedWidth) {
  fitgammas_[channelNum-1]=reducedWidth;
}

/*!
 * This function sets the fitted level energy.
 */

void ALevel::SetFitE(double energy) {
  fitlevel_e_=energy;
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
  sqrt_nf_factor_=term;
}

/*!
 * This function adds a conversion factor from reduced width amplitude to ANC.
 */

void ALevel::AddECConversionFactor(double conversionFactor) {
  ec_conv_factors_.push_back(conversionFactor);
}

/*!
 * This function sets the physical reduced width amplitude for a given channel number.
 */

void ALevel::SetTransformGamma(int channelNum,double reducedWidth) {
  transform_gammas_[channelNum-1]=reducedWidth;
}

/*!
 * This function sets the physical level energy.
 */

void ALevel::SetTransformE(double energy) {
  transform_e_=energy;
}

/*!
 * This function sets the Breit-Wigner partial width for a given channel number.
 */

void ALevel::SetBigGamma(int channelNum, double partialWidth) {
  big_gammas_[channelNum-1]=partialWidth;
}

/*!
 * This function sets the number of iterations that were required for the transformation from formal to physical parameters.
 */

void ALevel::SetTransformIterations(int iterations) {
  transform_iter_=iterations;
}

/*!
 * This function sets the external reduced width amplitude for a given channel number.
 */

void ALevel::SetExternalGamma(int channelNum, complex reducedWidth) {
  external_gammas_[channelNum-1]=reducedWidth;
}

/*!
 * Sets the value of the shift function calculated at the resonance energy.
 */

void ALevel::SetShiftFunction(int channelNum, double shiftFunction) {
  shifts_[channelNum-1]=shiftFunction;
}

/*!
 * Sets the external capture parameters for the level.
 */ 

void ALevel::SetECParams(int pairNum, unsigned char multMask) {
  isECLevel_=true;
  ecPairNum_=pairNum;
  ecMultMask_=multMask;
}
