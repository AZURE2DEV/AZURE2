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
class EData;

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
};

class ParameterLimitsManager {
 public:
  ParameterLimitsManager(const Config *config, CNuc *compound, EData *data, AZUREParams *params);

  // Read parameter settings from AZURE2 file
  bool ReadParameterSettings();

  // Apply all parameter settings to Minuit parameters
  void ApplyAllParameterSettings(ROOT::Minuit2::MnUserParameters &p);

  // Check if parameter is marked as nuisance
  bool IsNuisanceParameter(const std::string &paramName) const;

  // Check if parameter is marked as nuisance by Minuit2 index
  bool IsNuisanceParameterByIndex(int nonFixedIndex) const;

  // Get parameter error for nuisance calculation
  double GetParameterError(const std::string &paramName) const;

  // Get converted nominal value (physical to reduced for width parameters)
  double GetConvertedNominalValue(const std::string &paramName) const;

  // Get converted error (physical to reduced for width parameters)
  double GetConvertedError(const std::string &paramName) const;

  // Get converted nominal value by Minuit2 index (physical to reduced for width parameters)
  double GetConvertedNominalValueByIndex(int nonFixedIndex) const;

  // Get converted error by Minuit2 index (physical to reduced for width parameters)
  double GetConvertedErrorByIndex(int nonFixedIndex) const;

 private:
  const Config *config_;
  CNuc *compound_;
  EData *data_;
  AZUREParams *params_;
  std::map<std::string, ParameterSetting> parameterSettings_;
  /// Settings entry that governs each enumerated non-fixed parameter, or NULL.
  /// Built by BuildIndexMap() and consulted by every ...ByIndex() lookup, so
  /// that a stale minuit_index in the file cannot attach one parameter's
  /// settings to a different parameter.
  std::vector<ParameterSetting *> indexToSetting_;

  // Internal functions
  int FindParameterIndex(const std::string &paramName) const;
  double ConvertPhysicalLimitToReduced(double physicalLimit, const std::string &paramName) const;
  void ApplyParameterSetting(const std::string &paramName, ROOT::Minuit2::MnUserParameters &p);
  void ApplyParameterSettingByIndex(int nonFixedIndex, int actualIndex, ROOT::Minuit2::MnUserParameters &p);
  /*!
   * Translates a Minuit parameter name into the name the GUI writes into the
   * <parameterSettings> block ("energy_2" -> "Level 2 Energy (MeV)",
   * "width_2_3" -> "Level 2 Channel 3 Width (eV)"; segment norms and energy
   * shifts already share their names).
   */
  static std::string SettingNameForMinuitName(const std::string &minuitName);
  /*!
   * Resolves each enumerated non-fixed parameter to its settings entry, by name
   * where possible and by the stored minuit_index only for entries no name
   * matched.  Must be called before any ...ByIndex() lookup.
   */
  void BuildIndexMap(const ROOT::Minuit2::MnUserParameters &p);
  ParameterSetting *SettingForIndex(int nonFixedIndex) const;
};

#endif