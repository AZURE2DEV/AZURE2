#ifndef AZURECALC_H
#define AZURECALC_H

#include "Minuit2/FCNBase.h"
#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;

///A function class to perform the calculation of the chi-squared value

/*!
 * The AZURECalc function class calculates the cross section based on a 
 * parameter set for all available data, and returns a chi-squared value.
 * This function class is what Minuit calls repeatedly during the fitting
 * process to perform the minimization.
 */

class AZURECalc : public ROOT::Minuit2::FCNBase {
 public:
  /*!
   * The AZURECalc object is created with reference to an EData and CNuc object.
   *. The runtime configurations are also passed through a Config structure.
   */
  AZURECalc(EData* data,CNuc* compound, const Config& configure) : configure_(configure) {
    data_=data;
    compound_=compound;
  };
  
  ~AZURECalc() {};
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  virtual double Up() const {return theErrorDef;};
  /*!
   * Overloaded operator to make the class instance callable as a function. 
   * A Minuit parameter array is passed as the dependent variable.  The function
   * returns the total chi-squared value.
   */
  virtual double operator()(const vector_r&) const;
  
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const {return configure_;};
  /*!
   * Returns a pointer to the EData object.
   */
  EData *data() const {return data_;};
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
 
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  void SetErrorDef(double def) {theErrorDef=def;};
 private:
  const Config &configure_;
  EData *data_;
  CNuc *compound_;
  double theErrorDef;
};

#endif
