#ifndef REACTIONRATE_H
#define REACTIONRATE_H

#include "Constants.h"
#include "LevelInputs.h"
#include "LevelState.h"
#include <vector>

class CNuc;

extern std::vector<double> gsl_reactionrate_breakpoints(const CNuc *, const LevelState &, const LevelInputs &,
                                                        const Config &, int);
extern double gsl_reactionrate_integration(double, const CNuc *, const LevelState &, const Config &, int, int,
                                           const std::vector<double> &);

/// A container structure for a reaction rate.

/*!
 * The RateData container structure hold a temperature and the corresponding reaction rate.
 */

class RateData {
 public:
  /// Constructor creates RateData object from a given temperature and rate value.
  RateData(double t, double r) :
    temperature(t),
    rate(r) {};
  /// This function defines the "less than" operator for use in sorting.
  bool operator<(const RateData &right) const {
    return temperature < right.temperature;
  };
  /// Temperature at which the rate was calculated.
  double temperature;
  /// Reaction rate at corresponding temperature.
  double rate;
};

/// A function class to calculate the reaction rate

/*!
 * The ReactionRate function class is used to calculate the reaction rate based on a set of
 * R-Matrix parameters over a range of stellar temperatures.
 */

class ReactionRate {
 public:
  /*!
   * The ReactionRate object is created with reference to a CNuc object (shared and constant; the caller outlives
   * this object), the caller's levels (`state`, copied into the rate's own and filled from `params`; `inputs`, read),
   * a vector of Minuit parameters, a Config structure, and a set of entrance and exit pair keys.  It warns of large
   * radiative widths at `params` and, under Brune, computes its state's shift functions.
   */
  ReactionRate(const CNuc *compound, const LevelState &state, const LevelInputs &inputs, const vector_r &params,
               const Config &, int, int);
  /*!
   * Returns a pointer to the CNuc object.
   */
  const CNuc *compound() const { return compound_; };
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const { return configure_; };
  /*!
   * Returns the entrance pair key.
   */
  int entranceKey() const { return entrance_key_; };
  /*!
   * Returns the exit pair key; -1 is the total radiative capture, the sum over every capture
   * (gamma) pair of the entrance pair's compound nucleus (CNuc::CapturePairKeys).
   */
  int exitKey() const { return exit_key_; };
  /// Reaction rate on the temperature grid, writing reactionrates.out.
  void CalculateRates();
  /// As above, on temperatures read from a file.
  void CalculateFileRates();
  /*!
   * Writes the rates to an output file.
   */
  void WriteRates();
  /// The rates computed (temperature in GK, N_A<sigma v> in cm^3 mol^-1 s^-1), in the order WriteRates writes them.
  const std::vector<RateData> &Rates() const { return rates_; }
  /// The rate's levels: the caller's state filled at the rate's parameters (each thread evaluates a copy of it).
  const LevelState &State() const { return state_; }

 private:
  int entrance_key_;
  int exit_key_;
  const CNuc *compound_;
  const LevelInputs &inputs_;  ///< the caller's
  const Config &configure_;
  LevelState state_;  ///< the rate's own
  std::vector<RateData> rates_;
};


#endif
