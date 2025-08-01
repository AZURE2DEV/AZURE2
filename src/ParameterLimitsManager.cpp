#include "ParameterLimitsManager.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

ParameterLimitsManager::ParameterLimitsManager(const Config* config, CNuc* compound, EData* data, AZUREParams* params)
    : config_(config), compound_(compound), data_(data), params_(params) {
}

bool ParameterLimitsManager::ReadParameterSettings() {
  parameterSettings_.clear();
  
  std::ifstream file(config_->configfile.c_str());
  if(!file.is_open()) {
    return false;
  }
  
  std::string line;
  bool inParameterSettings = false;
  
  while(std::getline(file, line)) {
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t\r\n");
    if(start == std::string::npos) continue;
    line = line.substr(start);
    size_t end = line.find_last_not_of(" \t\r\n");
    if(end != std::string::npos) line = line.substr(0, end + 1);
    
    if(line == "<parameterSettings>") {
      inParameterSettings = true;
      continue;
    } else if(line == "</parameterSettings>" || line.substr(0, 1) == "<") {
      inParameterSettings = false;
      continue;
    }
    
    if(inParameterSettings && !line.empty() && line.substr(0, 1) != "#") {
      std::istringstream iss(line);
      std::string name, category;
      double value, lowerLimit, upperLimit, error;
      int nuisance;
      
      if(iss >> name >> value >> lowerLimit >> upperLimit >> error >> nuisance >> category) {
        ParameterSetting setting;
        setting.lowerLimit = lowerLimit;
        setting.upperLimit = upperLimit;
        setting.error = error;
        setting.nominalValue = value;  // Store the nominal value from config file
        setting.useAsNuisance = (nuisance == 1);
        setting.category = category;
        
        parameterSettings_[name] = setting;
      }
    }
  }
  
  file.close();
  return true;
}

void ParameterLimitsManager::ApplyAllParameterSettings(ROOT::Minuit2::MnUserParameters& p) {  
  for(int i = 0; i < p.Params().size(); i++) {
    std::string paramName = p.Parameter(i).GetName();
    ApplyParameterSetting(paramName, p);
  }
}

void ParameterLimitsManager::ApplyParameterSetting(const std::string& paramName, ROOT::Minuit2::MnUserParameters& p) {
  auto it = parameterSettings_.find(paramName);

  if(it == parameterSettings_.end()) {
    return; // No settings found for this parameter
  }
  
  if(it != parameterSettings_.end()) {
    ParameterSetting& setting = it->second;
    
    // Apply limits if they are not both zero (0,0 means free parameter)
    if(setting.lowerLimit != 0.0 || setting.upperLimit != 0.0) {
      double lower = setting.lowerLimit;
      double upper = setting.upperLimit;

      // For width parameters (physical values need conversion to reduced)
      // Width parameters have format: j=%d_la=%d_ch=%d_rwa
      if(paramName.find("width") != std::string::npos && setting.category == "level") {
        // Convert physical width limits to reduced width limits
        lower = ConvertPhysicalLimitToReduced(setting.lowerLimit, paramName);
        upper = ConvertPhysicalLimitToReduced(setting.upperLimit, paramName);
      }

      // If nuisance parameter, convert nominal value and error
      if(setting.useAsNuisance) {
        setting.nominalValueReduced = ConvertPhysicalLimitToReduced(setting.nominalValue, paramName);
        double errorSmall = std::abs(setting.nominalValue - ConvertPhysicalLimitToReduced(setting.nominalValue - setting.error, paramName));
        double errorLarge = std::abs(ConvertPhysicalLimitToReduced(setting.nominalValue + setting.error, paramName) - setting.nominalValueReduced);
        setting.errorReduced = std::max(errorSmall, errorLarge);
      }

      // Ensure lower < upper, swap if necessary
      if(lower > upper) {
        std::swap(lower, upper);
      }
      
      // Set limits in Minuit
      if(p.Parameter(p.Index(paramName)).HasLimits()) {
        p.RemoveLimits(paramName);
      }
      p.SetLimits(paramName, lower, upper);
    }
  }
}

int ParameterLimitsManager::FindParameterIndex(const std::string& paramName) const {
  AZUREParams tempParams;
  compound_->FillMnParams(tempParams.GetMinuitParams());
  if(data_) data_->FillMnParams(tempParams.GetMinuitParams());
  
  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); ++i) {
    if (tempParams.GetMinuitParams().Parameter(i).GetName() == paramName) {
      return i;
    }
  }
  return -1;
}

double ParameterLimitsManager::ConvertPhysicalLimitToReduced(double physicalLimit, const std::string& paramName) const {
  int paramIndex = FindParameterIndex(paramName);
  if (paramIndex == -1) {
    return physicalLimit;
  }

  CNuc* clone = compound_->Clone();

  // Fill clone with current parameters
  clone->FillCompoundFromParams(params_->GetMinuitParams().Params());
  clone->CalcShiftFunctions(*config_);
  clone->TransformOut(*config_);
  
  try {
    vector_r currentParams = clone->GetTransformParams(*config_);
    
    currentParams[paramIndex] = physicalLimit;
    clone->FillCompoundFromParamsPhysical(currentParams);
    
    if (!clone->TransformIn(*config_)) {
      delete clone;
      return physicalLimit;
    }

    AZUREParams tempParams;
    clone->FillMnParams(tempParams.GetMinuitParams());

    vector_r reducedParams = tempParams.GetMinuitParams().Params();
    delete clone;
    
    if (paramIndex < reducedParams.size()) {
      double convertedLimit = reducedParams[paramIndex];
      return convertedLimit;
    } else {
      return physicalLimit;
    }
  } catch(...) {
    delete clone;
    return physicalLimit;
  }
}

bool ParameterLimitsManager::IsNuisanceParameter(const std::string& paramName) const {
  auto it = parameterSettings_.find(paramName);
  if(it != parameterSettings_.end()) {
    return it->second.useAsNuisance;
  }
  return false;
}

double ParameterLimitsManager::GetParameterError(const std::string& paramName) const {
  auto it = parameterSettings_.find(paramName);
  if(it != parameterSettings_.end()) {
    return it->second.error;
  }
  return 0.0;
}

double ParameterLimitsManager::GetConvertedNominalValue(const std::string& paramName) const {
  auto it = parameterSettings_.find(paramName);
  if(it != parameterSettings_.end()) {
    const ParameterSetting& setting = it->second;
    
    // For width parameters (physical values need conversion to reduced)
    // Width parameters have format: j=%d_la=%d_ch=%d_rwa
    if(paramName.find("width") != std::string::npos && setting.category == "level") {
      // Convert physical nominal value to reduced
      return setting.nominalValueReduced;
    } else {
      // For non-width parameters, use nominal value as-is
      return setting.nominalValue;
    }
  }
  return 0.0;
}

double ParameterLimitsManager::GetConvertedError(const std::string& paramName) const {
  auto it = parameterSettings_.find(paramName);
  if(it != parameterSettings_.end()) {
    const ParameterSetting& setting = it->second;
    
    // For width parameters (physical values need conversion to reduced)
    // Width parameters have format: j=%d_la=%d_ch=%d_rwa
    if(paramName.find("width") != std::string::npos && setting.category == "level") {
      return setting.errorReduced; // Return reduced error for width parameters
    } else {
      // For non-width parameters, use error as-is
      return setting.error;
    }
  }
  return 0.0;
}