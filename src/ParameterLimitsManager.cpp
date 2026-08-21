#include "ParameterLimitsManager.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <regex>

inline double custom_stod(const std::string &str, std::size_t *pos = nullptr) {
  std::size_t i = 0;

  // Skip leading whitespace
  while (i < str.size() && std::isspace(static_cast<unsigned char>(str[i]))) {
    ++i;
  }

  if (i == str.size()) {
    throw std::invalid_argument("custom_stod: empty string or only whitespace");
  }

  // Handle sign
  int sign = 1;
  if (str[i] == '+') {
    ++i;
  } else if (str[i] == '-') {
    sign = -1;
    ++i;
  }

  // Parse integer part
  double value = 0.0;
  bool has_digits = false;
  while (i < str.size() && std::isdigit(static_cast<unsigned char>(str[i]))) {
    has_digits = true;
    value = value * 10.0 + (str[i] - '0');
    ++i;
  }

  // Parse fractional part
  if (i < str.size() && str[i] == '.') {
    ++i;
    double factor = 0.1;
    while (i < str.size() && std::isdigit(static_cast<unsigned char>(str[i]))) {
      has_digits = true;
      value += (str[i] - '0') * factor;
      factor *= 0.1;
      ++i;
    }
  }

  if (!has_digits) {
    throw std::invalid_argument("custom_stod: no digits found in input");
  }

  // Parse scientific notation (e/E)
  if (i < str.size() && (str[i] == 'e' || str[i] == 'E')) {
    ++i;
    int exp_sign = 1;
    if (i < str.size() && (str[i] == '+' || str[i] == '-')) {
      if (str[i] == '-') exp_sign = -1;
      ++i;
    }
    int exponent = 0;
    bool exp_digits = false;
    while (i < str.size() && std::isdigit(static_cast<unsigned char>(str[i]))) {
      exp_digits = true;
      exponent = exponent * 10 + (str[i] - '0');
      ++i;
    }
    if (!exp_digits) {
      throw std::invalid_argument("custom_stod: malformed exponent");
    }

    double pow10 = 1.0;
    int abs_exp = exponent;
    while (abs_exp--) {
      pow10 *= 10.0;
    }
    if (exp_sign > 0)
      value *= pow10;
    else
      value /= pow10;
  }

  if (pos) {
    *pos = i;
  }

  return sign * value;
}

ParameterLimitsManager::ParameterLimitsManager(const Config *config, CNuc *compound, EData *data, AZUREParams *params) :
  config_(config),
  compound_(compound),
  data_(data),
  params_(params) {
}

bool ParameterLimitsManager::ReadParameterSettings() {
  parameterSettings_.clear();

  std::ifstream file(config_->configfile.c_str());
  if (!file.is_open()) {
    return false;
  }

  std::string line;
  bool inParameterSettings = false;

  while (std::getline(file, line)) {
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) continue;
    line = line.substr(start);
    size_t end = line.find_last_not_of(" \t\r\n");
    if (end != std::string::npos) line = line.substr(0, end + 1);

    if (line == "<parameterSettings>") {
      inParameterSettings = true;
      continue;
    } else if (line == "</parameterSettings>" || line.substr(0, 1) == "<") {
      inParameterSettings = false;
      continue;
    }

    if (inParameterSettings && !line.empty() && line.substr(0, 1) != "#") {
      std::istringstream iss(line);
      std::vector<std::string> parts;
      std::string part;

      // Split the line by spaces
      while (iss >> part) {
        parts.push_back(part);
      }

      if (parts.size() == 9) {
        // Format: "segment_1_norm 1 0 0 1 0 0 norm 102" (9 parts)
        ParameterSetting setting;
        setting.name = parts[0];
        setting.nominalValue = custom_stod(parts[1]);
        setting.lowerLimit = custom_stod(parts[2]);
        setting.upperLimit = custom_stod(parts[3]);
        setting.error = custom_stod(parts[4]);
        setting.fitError = custom_stod(parts[5]);  // fitError field
        setting.useAsNuisance = (std::stoi(parts[6]) == 1);
        setting.category = parts[7];
        setting.minuitIndex = std::stoi(parts[8]);  // Store Minuit2 index

        parameterSettings_[setting.name] = setting;
      } else if (parts.size() == 12) {
        // Format: "Level 2 Energy (MeV) 11.2914 0 0 0.01 0 0 level 0" (12 parts)
        ParameterSetting setting;
        setting.name = parts[0] + " " + parts[1] + " " + parts[2] + " " + parts[3];  // "Level 2 Energy (MeV)"
        setting.nominalValue = custom_stod(parts[4]);
        setting.lowerLimit = custom_stod(parts[5]);
        setting.upperLimit = custom_stod(parts[6]);
        setting.error = custom_stod(parts[7]);
        setting.fitError = custom_stod(parts[8]);
        setting.useAsNuisance = (std::stoi(parts[9]) == 1);
        setting.category = parts[10];
        setting.minuitIndex = std::stoi(parts[11]);  // Store Minuit2 index

        parameterSettings_[setting.name] = setting;
      } else if (parts.size() == 14) {
        // Format: "Level 2 Channel 6 Width (eV) -0.000119735 0 0 -1.19735e-05 0 0 level 1" (14 parts)
        ParameterSetting setting;
        setting.name = parts[0] + " " + parts[1] + " " + parts[2] + " " + parts[3] + " " + parts[4] + " " + parts[5];  // "Level 2 Channel 6 Width (eV)"
        setting.nominalValue = custom_stod(parts[6]);
        setting.lowerLimit = custom_stod(parts[7]);
        setting.upperLimit = custom_stod(parts[8]);
        setting.error = custom_stod(parts[9]);
        setting.fitError = custom_stod(parts[10]);
        setting.useAsNuisance = (std::stoi(parts[11]) == 1);
        setting.category = parts[12];
        setting.minuitIndex = std::stoi(parts[13]);  // Store Minuit2 index

        parameterSettings_[setting.name] = setting;
      }
    }
  }

  file.close();
  return true;
}

std::string ParameterLimitsManager::SettingNameForMinuitName(const std::string &minuitName) {
  if (minuitName.compare(0, 7, "energy_") == 0)
    return "Level " + minuitName.substr(7) + " Energy (MeV)";
  if (minuitName.compare(0, 6, "width_") == 0) {
    size_t sep = minuitName.find('_', 6);
    if (sep != std::string::npos)
      return "Level " + minuitName.substr(6, sep - 6) +
          " Channel " + minuitName.substr(sep + 1) + " Width (eV)";
  }
  // Normalizations and energy shifts are stored under their Minuit names already.
  return minuitName;
}

/*!
 * Associates every enumerated non-fixed parameter with the settings entry that
 * describes it.
 *
 * The names are authoritative: the minuit_index recorded in the file is written
 * by the GUI from its own view of which parameters are free, and goes stale as
 * soon as a normalization or energy shift is freed or fixed in the Segments tab
 * without the Fitting tab rebuilding its list.  Trusting it then hands one
 * parameter's limits and nuisance prior to a different parameter.  The index is
 * still honoured, but only for entries whose name matched nothing -- which keeps
 * older files working without letting a stale index override a name match.
 */

void ParameterLimitsManager::BuildIndexMap(const ROOT::Minuit2::MnUserParameters &p) {
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < (int)p.Params().size(); i++) {
    if (!p.Parameter(i).IsFixed() || p.Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  indexToSetting_.assign(nonFixedToActualIndex.size(), NULL);
  std::map<const ParameterSetting *, bool> claimed;

  for (size_t nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    std::string settingName =
        SettingNameForMinuitName(p.Parameter(nonFixedToActualIndex[nonFixedIndex]).GetName());
    std::map<std::string, ParameterSetting>::iterator it = parameterSettings_.find(settingName);
    if (it != parameterSettings_.end()) {
      indexToSetting_[nonFixedIndex] = &it->second;
      claimed[&it->second] = true;
    }
  }

  for (size_t nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    if (indexToSetting_[nonFixedIndex]) continue;
    for (std::map<std::string, ParameterSetting>::iterator it = parameterSettings_.begin();
         it != parameterSettings_.end(); ++it) {
      if (it->second.minuitIndex == (int)nonFixedIndex && !claimed.count(&it->second)) {
        indexToSetting_[nonFixedIndex] = &it->second;
        claimed[&it->second] = true;
        break;
      }
    }
  }
}

ParameterSetting *ParameterLimitsManager::SettingForIndex(int nonFixedIndex) const {
  if (nonFixedIndex < 0 || nonFixedIndex >= (int)indexToSetting_.size()) return NULL;
  return indexToSetting_[nonFixedIndex];
}

void ParameterLimitsManager::ApplyAllParameterSettings(ROOT::Minuit2::MnUserParameters &p) {
  // Build mapping from non-fixed parameter index to actual parameter index
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < p.Params().size(); i++) {
    if (!p.Parameter(i).IsFixed() || p.Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  BuildIndexMap(p);

  // Apply settings using non-fixed parameter indices
  for (int nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    ApplyParameterSettingByIndex(nonFixedIndex, actualIndex, p);
  }
}

void ParameterLimitsManager::ApplyParameterSetting(const std::string &paramName, ROOT::Minuit2::MnUserParameters &p) {
  // Try to find parameter settings by the GUI name format
  // Since we now store parameters with their full GUI names like "Level 2 Energy (MeV)"
  // we need to find the setting that corresponds to the Minuit parameter name like "energy_2"

  ParameterSetting *setting = nullptr;
  std::string matchedName = "";

  // Search through all parameter settings to find a match
  for (auto &entry : parameterSettings_) {
    const std::string &guiName = entry.first;

    // Try to match Minuit parameter names to GUI parameter names

    // Energy parameters: paramName="energy_2" should match guiName="Level 2 Energy (MeV)"
    if (paramName.find("energy_") == 0) {
      std::string levelNum = paramName.substr(7);  // Extract number after "energy_"
      std::string expectedName = "Level " + levelNum + " Energy (MeV)";
      if (guiName == expectedName) {
        setting = &entry.second;
        matchedName = guiName;
        break;
      }
    }

    // Width parameters: paramName="width_2_3" should match guiName="Level 2 Channel 3 Width (eV)"
    else if (paramName.find("width_") == 0) {
      size_t firstUnderscore = paramName.find("_", 6);  // Find underscore after "width_"
      if (firstUnderscore != std::string::npos) {
        std::string levelNum = paramName.substr(6, firstUnderscore - 6);
        std::string channelNum = paramName.substr(firstUnderscore + 1);
        std::string expectedName = "Level " + levelNum + " Channel " + channelNum + " Width (eV)";
        if (guiName == expectedName) {
          setting = &entry.second;
          matchedName = guiName;
          break;
        }
      }
    }

    // Normalization parameters: paramName="segment_1_norm" should match guiName="segment_1_norm"
    else if (paramName.find("_norm") != std::string::npos || paramName.find("_energy_shift") != std::string::npos) {
      if (guiName == paramName) {
        setting = &entry.second;
        matchedName = guiName;
        break;
      }
    }
  }

  if (!setting) {
    return;  // No settings found for this parameter
  }

  // Apply limits if they are not both zero (0,0 means free parameter)
  if (setting->lowerLimit != 0.0 || setting->upperLimit != 0.0) {
    double lower = setting->lowerLimit;
    double upper = setting->upperLimit;

    // For width parameters (physical values need conversion to reduced)
    if ((matchedName.find("width") != std::string::npos || matchedName.find("Width") != std::string::npos) && setting->category == "level") {
      // Convert physical width limits to reduced width limits
      lower = ConvertPhysicalLimitToReduced(setting->lowerLimit, paramName);
      upper = ConvertPhysicalLimitToReduced(setting->upperLimit, paramName);
    }

    // If nuisance parameter, convert nominal value and error
    if (setting->useAsNuisance && (matchedName.find("width") != std::string::npos || matchedName.find("Width") != std::string::npos)) {
      setting->nominalValueReduced = ConvertPhysicalLimitToReduced(setting->nominalValue, paramName);
      double errorSmall = std::abs(setting->nominalValue - ConvertPhysicalLimitToReduced(setting->nominalValue - setting->error, paramName));
      double errorLarge = std::abs(ConvertPhysicalLimitToReduced(setting->nominalValue + setting->error, paramName) - setting->nominalValueReduced);
      setting->errorReduced = std::max(errorSmall, errorLarge);
    }

    // Ensure lower < upper, swap if necessary
    if (lower > upper) {
      std::swap(lower, upper);
    }

    // Set limits in Minuit
    if (p.Parameter(p.Index(paramName)).HasLimits()) {
      p.RemoveLimits(paramName);
    }
    p.SetLimits(paramName, lower, upper);
  }
}

void ParameterLimitsManager::ApplyParameterSettingByIndex(int nonFixedIndex, int actualIndex, ROOT::Minuit2::MnUserParameters &p) {
  ParameterSetting *setting = SettingForIndex(nonFixedIndex);

  if (!setting) {
    return;  // No settings found for this parameter index
  }
  const std::string &matchedName = setting->name;

  // Apply limits if they are not both zero (0,0 means free parameter)
  if (setting->lowerLimit != 0.0 || setting->upperLimit != 0.0) {
    double lower = setting->lowerLimit;
    double upper = setting->upperLimit;

    // For width parameters (physical values need conversion to reduced)
    if ((matchedName.find("width") != std::string::npos || matchedName.find("Width") != std::string::npos) && setting->category == "level") {
      // Convert physical width limits to reduced width limits
      std::string paramName = p.Parameter(actualIndex).GetName();
      lower = ConvertPhysicalLimitToReduced(setting->lowerLimit, paramName);
      upper = ConvertPhysicalLimitToReduced(setting->upperLimit, paramName);
    }

    // Ensure lower < upper, swap if necessary
    if (lower > upper) {
      std::swap(lower, upper);
    }

    // Set limits in Minuit using actual parameter index
    std::string paramName = p.Parameter(actualIndex).GetName();
    if (p.Parameter(actualIndex).HasLimits()) {
      p.RemoveLimits(paramName);
    }
    p.SetLimits(paramName, lower, upper);
  }

  // If nuisance parameter, convert nominal value and error
  if ((matchedName.find("segment") == std::string::npos) && setting->useAsNuisance) {
    std::string paramName = p.Parameter(actualIndex).GetName();
    setting->nominalValueReduced = ConvertPhysicalLimitToReduced(setting->nominalValue, paramName);
    double errorSmall = std::abs(setting->nominalValueReduced - ConvertPhysicalLimitToReduced(setting->nominalValue - setting->error, paramName));
    double errorLarge = std::abs(ConvertPhysicalLimitToReduced(setting->nominalValue + setting->error, paramName) - setting->nominalValueReduced);
    setting->errorReduced = std::max(errorSmall, errorLarge);
  }
}

int ParameterLimitsManager::FindParameterIndex(const std::string &paramName) const {
  AZUREParams tempParams;
  compound_->FillMnParams(tempParams.GetMinuitParams(), config_);
  if (data_) data_->FillMnParams(tempParams.GetMinuitParams());

  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); ++i) {
    if (tempParams.GetMinuitParams().Parameter(i).GetName() == paramName) {
      return i;
    }
  }
  return -1;
}

double ParameterLimitsManager::ConvertPhysicalLimitToReduced(double physicalLimit, const std::string &paramName) const {
  int paramIndex = FindParameterIndex(paramName);
  if (paramIndex == -1) {
    return physicalLimit;
  }

  CNuc *clone = compound_->Clone();

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
    clone->FillMnParams(tempParams.GetMinuitParams(), config_);

    vector_r reducedParams = tempParams.GetMinuitParams().Params();
    delete clone;

    if (paramIndex < reducedParams.size()) {
      double convertedLimit = reducedParams[paramIndex];
      return convertedLimit;
    } else {
      return physicalLimit;
    }
  } catch (...) {
    delete clone;
    return physicalLimit;
  }
}

bool ParameterLimitsManager::IsNuisanceParameter(const std::string &paramName) const {
  // First try to find by Minuit non-fixed index if we can determine it
  int actualIndex = FindParameterIndex(paramName);
  if (actualIndex >= 0) {
    // Convert actual index to non-fixed index
    AZUREParams tempParams;
    compound_->FillMnParams(tempParams.GetMinuitParams(), config_);
    if (data_) data_->FillMnParams(tempParams.GetMinuitParams());

    int nonFixedIndex = 0;
    for (int i = 0; i <= actualIndex && i < tempParams.GetMinuitParams().Params().size(); i++) {
      if (!tempParams.GetMinuitParams().Parameter(i).IsFixed()) {
        if (i == actualIndex) {
          // Found the non-fixed index corresponding to actualIndex
          for (const auto &entry : parameterSettings_) {
            if (entry.second.minuitIndex == nonFixedIndex) {
              return entry.second.useAsNuisance;
            }
          }
          break;
        }
        nonFixedIndex++;
      }
    }
  }

  // Fallback to name-based search for compatibility
  for (const auto &entry : parameterSettings_) {
    const std::string &guiName = entry.first;

    // Energy parameters: paramName="energy_2" should match guiName="Level 2 Energy (MeV)"
    if (paramName.find("energy_") == 0) {
      std::string levelNum = paramName.substr(7);
      std::string expectedName = "Level " + levelNum + " Energy (MeV)";
      if (guiName == expectedName) {
        return entry.second.useAsNuisance;
      }
    }

    // Width parameters: paramName="width_2_3" should match guiName="Level 2 Channel 3 Width (eV)"
    else if (paramName.find("width_") == 0) {
      size_t firstUnderscore = paramName.find("_", 6);
      if (firstUnderscore != std::string::npos) {
        std::string levelNum = paramName.substr(6, firstUnderscore - 6);
        std::string channelNum = paramName.substr(firstUnderscore + 1);
        std::string expectedName = "Level " + levelNum + " Channel " + channelNum + " Width (eV)";
        if (guiName == expectedName) {
          return entry.second.useAsNuisance;
        }
      }
    }

    // Direct match for norm/shift parameters
    else if (guiName == paramName) {
      return entry.second.useAsNuisance;
    }
  }

  return false;
}

bool ParameterLimitsManager::IsNuisanceParameterByIndex(int nonFixedIndex) const {
  const ParameterSetting *setting = SettingForIndex(nonFixedIndex);
  return setting ? setting->useAsNuisance : false;
}

double ParameterLimitsManager::GetParameterError(const std::string &paramName) const {
  // First try to find by Minuit non-fixed index if we can determine it
  int actualIndex = FindParameterIndex(paramName);
  if (actualIndex >= 0) {
    // Convert actual index to non-fixed index
    AZUREParams tempParams;
    compound_->FillMnParams(tempParams.GetMinuitParams(), config_);
    if (data_) data_->FillMnParams(tempParams.GetMinuitParams());

    int nonFixedIndex = 0;
    for (int i = 0; i <= actualIndex && i < tempParams.GetMinuitParams().Params().size(); i++) {
      if (!tempParams.GetMinuitParams().Parameter(i).IsFixed()) {
        if (i == actualIndex) {
          // Found the non-fixed index corresponding to actualIndex
          for (const auto &entry : parameterSettings_) {
            if (entry.second.minuitIndex == nonFixedIndex) {
              return entry.second.error;
            }
          }
          break;
        }
        nonFixedIndex++;
      }
    }
  }

  // Fallback to name-based search for compatibility
  for (const auto &entry : parameterSettings_) {
    const std::string &guiName = entry.first;

    // Energy parameters
    if (paramName.find("energy_") == 0) {
      std::string levelNum = paramName.substr(7);
      std::string expectedName = "Level " + levelNum + " Energy (MeV)";
      if (guiName == expectedName) {
        return entry.second.error;
      }
    }
    // Width parameters
    else if (paramName.find("width_") == 0) {
      size_t firstUnderscore = paramName.find("_", 6);
      if (firstUnderscore != std::string::npos) {
        std::string levelNum = paramName.substr(6, firstUnderscore - 6);
        std::string channelNum = paramName.substr(firstUnderscore + 1);
        std::string expectedName = "Level " + levelNum + " Channel " + channelNum + " Width (eV)";
        if (guiName == expectedName) {
          return entry.second.error;
        }
      }
    }
    // Direct match for norm/shift parameters
    else if (guiName == paramName) {
      return entry.second.error;
    }
  }

  return 0.0;
}

double ParameterLimitsManager::GetConvertedNominalValue(const std::string &paramName) const {
  auto it = parameterSettings_.find(paramName);
  if (it != parameterSettings_.end()) {
    const ParameterSetting &setting = it->second;

    // For width parameters (physical values need conversion to reduced)
    // Width parameters have format: j=%d_la=%d_ch=%d_rwa
    if (paramName.find("width") != std::string::npos && setting.category == "level") {
      // Convert physical nominal value to reduced
      return setting.nominalValueReduced;
    } else {
      // For non-width parameters, use nominal value as-is
      return setting.nominalValue;
    }
  }
  return 0.0;
}

double ParameterLimitsManager::GetConvertedError(const std::string &paramName) const {
  auto it = parameterSettings_.find(paramName);
  if (it != parameterSettings_.end()) {
    const ParameterSetting &setting = it->second;

    // For width parameters (physical values need conversion to reduced)
    // Width parameters have format: j=%d_la=%d_ch=%d_rwa
    if (paramName.find("width") != std::string::npos && setting.category == "level") {
      return setting.errorReduced;  // Return reduced error for width parameters
    } else {
      // For non-width parameters, use error as-is
      return setting.error;
    }
  }
  return 0.0;
}

double ParameterLimitsManager::GetConvertedNominalValueByIndex(int nonFixedIndex) const {
  // Find parameter settings by matching Minuit2 non-fixed index
  const ParameterSetting *setting = SettingForIndex(nonFixedIndex);
  if (!setting) return 0.0;

  // For width parameters (physical values need conversion to reduced)
  if ((setting->name.find("width") != std::string::npos || setting->name.find("Width") != std::string::npos) && setting->category == "level")
    return setting->nominalValueReduced;  // Return reduced nominal value for width parameters
  // For non-width parameters, use nominal value as-is
  return setting->nominalValue;
}

double ParameterLimitsManager::GetConvertedErrorByIndex(int nonFixedIndex) const {
  // Find parameter settings by matching Minuit2 non-fixed index
  const ParameterSetting *setting = SettingForIndex(nonFixedIndex);
  if (!setting) return 0.0;

  // For width parameters (physical values need conversion to reduced)
  if ((setting->name.find("width") != std::string::npos || setting->name.find("Width") != std::string::npos) && setting->category == "level")
    return setting->errorReduced;  // Return reduced error for width parameters
  // For non-width parameters, use error as-is
  return setting->error;
}