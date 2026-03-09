#ifndef REACTIONRATE_H
#define REACTIONRATE_H

#include "Constants.h"

class CNuc;

extern double gsl_reactionrate_integration(double,CNuc*,const Config&,int,int);

///A container structure for a reaction rate.

/*!
 * The RateData container structure hold a temperature and the corresponding reaction rate.
 */ 

class RateData {
 public:
 /// Constructor creates RateData object from a given temperature and rate value.
 RateData(double t, double r) :
  temperature(t),rate(r) {};
 /// This function defines the "less than" operator for use in sorting.
  bool operator<(const RateData& right) const {
    return temperature < right.temperature;
  };
  /// Temperature at which the rate was calculated.
  double temperature;
  /// Reaction rate at corresponding temperature.
  double rate;
};

///A function class to calculate the reaction rate

/*!
 * The ReactionRate function class is used to calculate the reaction rate based on a set of
 * R-Matrix parameters over a range of stellar temperatures.
 */

class ReactionRate {
 public:
  /*!
   * The ReactionRate object is created with reference to a CNuc object, a vector of Minuit parameters,
   * a Config structure, and a set of entrance and exit pair keys.
   */
  ReactionRate(CNuc*, const vector_r&, const Config &, int, int);
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const {return configure_;};
  /*!
   * Returns the entrance pair key.
   */
  int entranceKey() const {return entrance_key_;};
  /*!
   * Returns the exit pair key.
   */
  int exitKey() const {return exit_key_;};
  void CalculateRates();
  void CalculateFileRates();
  /*!
   * Writes the rates to an output file.
   */
  void WriteRates();
 private:
  int entrance_key_;
  int exit_key_;
  CNuc *compound_;
  const Config &configure_;
  std::vector<RateData> rates_;
};



#endif
