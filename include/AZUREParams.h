#ifndef AZUREPARAMS_H
#define AZUREPARAMS_H

#include "Minuit2/MnUserParameters.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "Constants.h"
#include "InferenceTypes.h"

class Config;

/// A container class to hold Minuit parameters in AZURE

/*!
 * The AZUREParams class holds the Minuit parameters determined in the fit and,
 * beside each, its ParameterRecord: what the parameter is and what it belongs
 * to.  The class also has member functions corresponding to reading and writing
 * of the parameters and their errors.
 */

class CompoundSystem;
class LevelState;

class AZUREParams {
 public:
  ROOT::Minuit2::MnUserParameters &GetMinuitParams();
  const ROOT::Minuit2::MnUserParameters &GetMinuitParams() const { return params_; }
  /// Append a parameter (value and Minuit step in the unit of its kind, ParamKindTable()) with its record.  The
  /// creation sites (CompoundSystem::FillMnParams, WorkingCopy::FillMnParams) are the only callers: the record is written where the
  /// parameter is made, the one place that knows what it is.
  void Add(const std::string &name, double value, double step, ParameterRecord record);
  /// What parameter i (0-based, full vector) is; throws std::logic_error if the vector holds a parameter added
  /// without its record.
  const ParameterRecord &Record(int i) const;
  const std::vector<ParameterRecord> &Records() const { return records_; }
  ParamKind Kind(int i) const { return Record(i).kind; }
  /// Index of the parameter whose name is exactly `name`, or -1: how a name in a file (param.sav, priors,
  /// <parameterSettings>) is resolved to its parameter.  Whole names only, never a prefix.
  int IndexOfName(const std::string &name) const;
  /// Read starting parameters from the configured external parameter file.  False (after the error) when no name of
  /// the file is a parameter of the run; a partial match is a warning (experiment/ParameterFileNames.h).
  bool ReadUserParameters(const Config &);
  /// Read starting parameters from a named file; false as above.
  bool ReadUserParameters(const std::string &);
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
  bool ReconcileBasis(const CompoundSystem *, const Config &);
  /// After ReadUserParameters and ReconcileBasis: the --ignore-externals mask of the photon channels the file
  /// named follows the file (docs/upgrade/external_capture.md section 3).  A photon amplitude the file gives as
  /// exactly 0 is "no gamma decay of the level to that final state": masked (no external part) and fixed, whatever
  /// <levels> says; one it gives as non-zero is not masked.  Brune's and Park's bases only (in the formal basis
  /// the file's amplitudes are not the Brune levels', on which the mask is defined).  Says what it changed.
  void ApplyFileWidthMask(const CompoundSystem *, LevelState &, const Config &);
  /// The basis a run writes, from its runtime options.
  static int BasisForConfig(const Config &);
  static const char *BasisName(int);

 private:
  /// The body of both ReadUserParameters overloads; "could not read" goes to out.  With
  /// segmentEntries false (a run without data: extrapolation, reaction rates) the file's rows for
  /// data-segment parameters (ParamKindRow::segmentField: normalizations, energy shifts, thicknesses)
  /// are not applied.  `otherRunNames`: the names the project gives to parameters of another run (a format-3
  /// project's data sets in an extrapolation run and the reverse), neither applied nor unknown.
  bool ApplyUserParameterFile(const std::string &filename, const std::string &pathToOpen, std::ostream &out, bool segmentEntries = true,
                              const std::vector<std::string> &otherRunNames = {});
  ROOT::Minuit2::MnUserParameters params_;
  std::vector<ParameterRecord> records_;  ///< one per entry of params_, in its order
  std::vector<int> kindCount_;            ///< parameters of each kind so far (the next ordinal)
  int basisTag_ = kBasisUnknown;
  std::vector<int> fileSet_;  ///< the parameters the last parameter file set, by index
};

#endif
