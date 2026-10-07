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

class CNuc;

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

  /// Which amplitudes a parameter file holds: the comment-style line
  /// `#parametrization` that WriteUserParameters puts first.  Files written
  /// before the line existed carry none and read back as kBasisUnknown.
  enum Basis { kBasisUnknown = -1, kBasisStandard = 0, kBasisBrune = 1, kBasisPark = 2 };
  /// The basis the last ReadUserParameters found in its file.
  int BasisTag() const { return basisTag_; }
  /// After ReadUserParameters: convert the amplitudes to this run's basis if
  /// the file was written in the other alternative one (Brune <-> Park), and
  /// say so.  False when no conversion exists (standard <-> alternative, or a
  /// Park level with J <= 0), after printing why.
  bool ReconcileBasis(CNuc *, const Config &);
  /// The basis a run writes, from its runtime options.
  static int BasisForConfig(const Config &);
  static const char *BasisName(int);

 private:
  ROOT::Minuit2::MnUserParameters params_;
  int basisTag_ = kBasisUnknown;
};

#endif
