#ifndef PARAMETERLIMITSMANAGER_H
#define PARAMETERLIMITSMANAGER_H

#include <string>
#include <map>
#include <vector>
#include <Minuit2/MnUserParameters.h>
#include "AZUREParams.h"

// Forward declarations
class Config;
class CNuc;
class LevelInputs;
class LevelState;
class WorkingCopy;

// Structure to hold parameter settings
struct ParameterSetting {
  std::string name;  // Parameter name
  double lowerLimit;
  double upperLimit;
  double error;
  double nominalValue;         // Nominal/central value for nuisance parameters
  double nominalValueReduced;  // Reduced value for width parameters
  double errorReduced;         // Reduced error for width parameters
  double fitError;             // Error from fitting results
  bool useAsNuisance;
  std::string category;
  int minuitIndex;  // Minuit2 parameter index
  /// Stored under a parameter's own name (the one-token form): matched by that name or not at all.
  bool parameterNameForm = false;
  int line = 0;  ///< 1-based line of the row in the project file
};

class ParameterLimitsManager {
 public:
  /// `state` and `inputs` are the caller's levels (the run's, the session's), which outlive this object: a conversion
  /// of a physical limit works on copies of them.
  ParameterLimitsManager(const Config *config, const CNuc *compound, const LevelState &state, const LevelInputs &inputs,
                         WorkingCopy *data, AZUREParams *params);

  // Read parameter settings from AZURE2 file
  bool ReadParameterSettings();

  // Apply all parameter settings to the parameters (limits; nuisance priors for NuisancePrior)
  void ApplyAllParameterSettings(AZUREParams &params);

  // Check if parameter is marked as nuisance by Minuit2 index
  bool IsNuisanceParameterByIndex(int nonFixedIndex) const;

  /*!
   * Gaussian prior on the parameter at position actualIndex of the FULL Minuit
   * parameter vector (fixed parameters included -- the vector the chi-squared
   * function, its gradient and CNuc::FillCompoundFromParams all index).
   * Returns false when that parameter carries no nuisance prior; otherwise sets
   * nominal and error in fit units (reduced amplitudes for widths).  Segment
   * normalizations and energy shifts always return false: their priors come
   * from the segment definition and are added by the segment code.
   *
   * Callers that hold a full-vector index should use this rather than the
   * ...ByIndex() lookups, whose non-fixed numbering follows the fixed flags of
   * the parameters this manager was applied to (which, with a param.par file,
   * need not match a freshly filled parameter list).
   */
  /// As NuisancePrior, in the parameter's own units also for a width (the physical width, eV): an MCMC that samples
  /// physical parameters.
  bool NuisancePriorPhysical(int actualIndex, double &nominal, double &error) const;
  bool NuisancePrior(int actualIndex, double &nominal, double &error) const;

 private:
  const Config *config_;
  const CNuc *compound_;
  const LevelState &state_;
  const LevelInputs &inputs_;
  WorkingCopy *data_;
  AZUREParams *params_;
  std::map<std::string, ParameterSetting> parameterSettings_;
  /// Settings entry that governs each enumerated non-fixed parameter, or NULL.
  /// Built by BuildIndexMap() and consulted by every ...ByIndex() lookup, so
  /// that a stale minuit_index in the file cannot attach one parameter's
  /// settings to a different parameter.
  std::vector<ParameterSetting *> indexToSetting_;
  /// Nuisance-prior settings entry per full-vector parameter index (NULL for
  /// none, and for segment norms/shifts).  Built by BuildIndexMap().
  std::vector<ParameterSetting *> priorByActualIndex_;

  /// Whether the nuisance prior of each full-vector parameter is on a level width (its settings are partial widths,
  /// eV, converted to reduced-width amplitudes).  Built by BuildIndexMap().
  std::vector<bool> priorIsWidth_;

  // Internal functions
  /// The parameter at full-vector index paramIndex set to a physical value (a level width in eV), as the fit
  /// parameter (the reduced-width amplitude) that value corresponds to.  Called for widths only; an index outside
  /// the level parameters returns the value unchanged.
  double ConvertPhysicalLimitToReduced(double physicalLimit, int paramIndex) const;
  void ApplyParameterSettingByIndex(int nonFixedIndex, int actualIndex, AZUREParams &params);
  /*!
   * The name under which the GUI writes a parameter into the <parameterSettings>
   * block: a level energy as "Level 2 Energy (MeV)", a width as "Level 2
   * Channel 3 Width (eV)" (from the parameter's record), every other parameter
   * under its own name.
   */
  static std::string SettingNameFor(const ParameterRecord &record, const std::string &name);
  /*!
   * Resolves each enumerated non-fixed parameter to its settings entry, by name
   * where possible and by the stored minuit_index only for entries no name
   * matched.  Must be called before any ...ByIndex() lookup.
   */
  void BuildIndexMap(const AZUREParams &params);
  /// W_PARAMETER_SETTING_UNKNOWN for each row stored under a name that is no parameter of the project
  /// (experiment::UnknownParameterSettings): such a row is ignored.
  void WarnUnknownSettings(const AZUREParams &params) const;
  ParameterSetting *SettingForIndex(int nonFixedIndex) const;
};

#endif