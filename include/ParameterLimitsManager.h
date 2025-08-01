#ifndef PARAMETERLIMITSMANAGER_H
#define PARAMETERLIMITSMANAGER_H

#include <string>
#include <map>
#include <Minuit2/MnUserParameters.h>
#include "AZUREParams.h"

// Forward declarations
class Config;
class CNuc;
class EData;

// Structure to hold parameter settings
struct ParameterSetting {
  double lowerLimit;
  double upperLimit;
  double error;
  double nominalValue;  // Nominal/central value for nuisance parameters
  double nominalValueReduced; // Reduced value for width parameters
  double errorReduced; // Reduced error for width parameters
  bool useAsNuisance;
  std::string category;
};

class ParameterLimitsManager {
public:
  ParameterLimitsManager(const Config* config, CNuc* compound, EData* data, AZUREParams* params);
  
  // Read parameter settings from AZURE2 file
  bool ReadParameterSettings();
  
  // Apply all parameter settings to Minuit parameters
  void ApplyAllParameterSettings(ROOT::Minuit2::MnUserParameters& p);
  
  // Check if parameter is marked as nuisance
  bool IsNuisanceParameter(const std::string& paramName) const;
  
  // Get parameter error for nuisance calculation
  double GetParameterError(const std::string& paramName) const;
  
  // Get converted nominal value (physical to reduced for width parameters)
  double GetConvertedNominalValue(const std::string& paramName) const;
  
  // Get converted error (physical to reduced for width parameters)
  double GetConvertedError(const std::string& paramName) const;
  
private:
  const Config* config_;
  CNuc* compound_;
  EData* data_;
  AZUREParams* params_;
  std::map<std::string, ParameterSetting> parameterSettings_;
  
  // Internal functions
  int FindParameterIndex(const std::string& paramName) const;
  double ConvertPhysicalLimitToReduced(double physicalLimit, const std::string& paramName) const;
  void ApplyParameterSetting(const std::string& paramName, ROOT::Minuit2::MnUserParameters& p);
};

#endif