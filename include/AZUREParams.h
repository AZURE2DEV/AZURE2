#ifndef AZUREPARAMS_H
#define AZUREPARAMS_H

#include "Minuit2/MnUserParameters.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "Constants.h"

class Config;

/// A container class to hold Minuit parameters in AZURE

/*!
 * The AZUREParams class holds the Minuit parameters determined in the fit.
 * The class also has member functions corresponding to reading and writing of
 * the parameters and their errors.
 */

class AZUREParams {
 public:
  ROOT::Minuit2::MnUserParameters &GetMinuitParams();
  /// Read starting parameters from the configured external parameter file.
  void ReadUserParameters(const Config &);
  /// Read starting parameters from a named file.
  void ReadUserParameters(const std::string &);
  /// Write the current parameters to param.sav.
  void WriteUserParameters(const Config &, bool);
  /// Write Minos asymmetric errors to param.errors.
  void WriteParameterErrors(const std::vector<std::pair<double, double>> &, const Config &);

 private:
  ROOT::Minuit2::MnUserParameters params_;
};

#endif
