#include "ParameterLimitsManager.h"
#include "Config.h"
#include "CNuc.h"
#include "WorkingCopy.h"
#include "AZUREParams.h"
#include "azure2/experiment/Build.h"
#include "azure2/experiment/ParameterFileNames.h"
#include "azure2/experiment/ProjectReader.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cctype>
#include <limits>

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

  // nan / inf (any case, optional sign already consumed): a fitError written from a singular covariance
  // (audit 2026-10-25, item 30) used to throw here and end the CLI run; they are legal doubles.
  {
    std::string rest = str.substr(i);
    for (std::size_t q = 0; q < rest.size(); q++) rest[q] = (char)std::tolower((unsigned char)rest[q]);
    if (rest.compare(0, 3, "nan") == 0 || rest.compare(0, 3, "inf") == 0) {
      const bool isNan = rest.compare(0, 3, "nan") == 0;
      i += 3;
      if (!isNan && rest.compare(0, 8, "infinity") == 0) i += 5;
      if (pos) *pos = i;
      return isNan ? std::numeric_limits<double>::quiet_NaN() : sign * std::numeric_limits<double>::infinity();
    }
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

ParameterLimitsManager::ParameterLimitsManager(const Config *config, const CNuc *compound, const LevelState &state,
                                               const LevelInputs &inputs, WorkingCopy *data, AZUREParams *params) :
  config_(config),
  compound_(compound),
  state_(state),
  inputs_(inputs),
  data_(data),
  params_(params) {
}

bool ParameterLimitsManager::ReadParameterSettings() {
  parameterSettings_.clear();

  {
    std::ifstream test(config_->configfile.c_str());
    if (!test.is_open()) return false;
  }
  int lineNumber = 0;  // the section's opening tag, then each line
  std::unique_ptr<std::istream> section = config_->SectionStream("parameterSettings", &lineNumber);
  if (!section) return true;  // no settings
  std::istream &file = *section;
  lineNumber--;

  std::string line;
  bool inParameterSettings = false;

  while (std::getline(file, line)) {
    lineNumber++;
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
        setting.parameterNameForm = true;
        setting.line = lineNumber;

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

  return true;
}

std::string ParameterLimitsManager::SettingNameFor(const ParameterRecord &record, const std::string &name) {
  if (record.kind == ParamKind::LevelEnergy) return "Level " + std::to_string(record.levelNumber) + " Energy (MeV)";
  if (record.kind == ParamKind::Gamma)
    return "Level " + std::to_string(record.levelNumber) + " Channel " + std::to_string(record.channel) + " Width (eV)";
  // Every other parameter is stored under its own name.
  return name;
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
 *
 * An entry whose name matched nothing is not always from an older file: it can
 * name a parameter that no longer exists (a segment made inactive, a level
 * renumbered or removed), and its stale index then points at an unrelated
 * parameter.  12C+a (2026-10-16): 33 such bindings, among them a level energy
 * given segment_113_norm's nuisance flag, which switched off the analytic G2.
 * So the index fallback only pairs compatible entries: an entry stored under a
 * parameter's own name (the one-token form the GUI writes for normalizations,
 * energy shifts and every parameter other than a level's) never falls back, and
 * a "Level N ..." entry only binds to a parameter of level N of the same kind
 * (energy or width, from the parameter's record; channel numbers are not
 * compared, older files numbered them globally).  Entries in any other format
 * fall back as before.
 */
static bool SettingFitsParameter(const ParameterSetting &setting, const ParameterRecord &param) {
  if (setting.parameterNameForm) return false;  // a parameter's own name: matched by name or not at all
  const std::string &name = setting.name;
  if (name.compare(0, 6, "Level ") != 0) return true;  // unknown format: keep the old behaviour
  const size_t sp = name.find(' ', 6);
  if (sp == std::string::npos) return true;
  const std::string level = name.substr(6, sp - 6);
  const bool isEnergy = name.find(" Energy") != std::string::npos;
  const bool isWidth = name.find(" Width") != std::string::npos;
  if (isEnergy) return param.kind == ParamKind::LevelEnergy && std::to_string(param.levelNumber) == level;
  if (isWidth) return param.kind == ParamKind::Gamma && std::to_string(param.levelNumber) == level;
  return true;
}


void ParameterLimitsManager::BuildIndexMap(const AZUREParams &params) {
  const ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // The GUI's numbering: the free parameters and every data-segment parameter, fixed or not.
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < (int)p.Params().size(); i++) {
    if (!p.Parameter(i).IsFixed() || params.Record(i).Row().segmentField) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  indexToSetting_.assign(nonFixedToActualIndex.size(), NULL);
  std::map<const ParameterSetting *, bool> claimed;

  for (size_t nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    const int i = nonFixedToActualIndex[nonFixedIndex];
    std::map<std::string, ParameterSetting>::iterator it =
        parameterSettings_.find(SettingNameFor(params.Record(i), p.GetName(i)));
    if (it != parameterSettings_.end()) {
      indexToSetting_[nonFixedIndex] = &it->second;
      claimed[&it->second] = true;
    }
  }

  for (size_t nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    if (indexToSetting_[nonFixedIndex]) continue;
    for (std::map<std::string, ParameterSetting>::iterator it = parameterSettings_.begin();
         it != parameterSettings_.end(); ++it) {
      if (it->second.minuitIndex == (int)nonFixedIndex && !claimed.count(&it->second) &&
          SettingFitsParameter(it->second, params.Record(nonFixedToActualIndex[nonFixedIndex]))) {
        indexToSetting_[nonFixedIndex] = &it->second;
        claimed[&it->second] = true;
        break;
      }
    }
  }

  // The same association keyed by the full-vector index, for the nuisance
  // prior.  The parameters with a prior of their own (ParamKindRow::ownPrior: norms, shifts, abundances,
  // electron-screening potentials, fitted target densities, Doppler temperatures and sample thicknesses) are left
  // out, as the chi-squared always has: they carry their own priors (segment errors, <abundances>, the screening,
  // densityfit, doppler and transmission / yield blocks).
  priorByActualIndex_.assign(p.Params().size(), NULL);
  priorIsWidth_.assign(p.Params().size(), false);
  for (size_t nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size(); nonFixedIndex++) {
    ParameterSetting *setting = indexToSetting_[nonFixedIndex];
    if (!setting || !setting->useAsNuisance) continue;
    const int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    if (params.Record(actualIndex).Row().ownPrior) continue;
    priorByActualIndex_[actualIndex] = setting;
    priorIsWidth_[actualIndex] = params.Kind(actualIndex) == ParamKind::Gamma && setting->category == "level";
  }
}

void ParameterLimitsManager::WarnUnknownSettings(const AZUREParams &params) const {
  std::vector<std::string> known =
      azure2::experiment::ParameterNamesOf(config_->project->built);
  for (unsigned int i = 0; i < params.GetMinuitParams().Params().size(); i++) known.push_back(params.GetMinuitParams().GetName(i));
  std::vector<azure2::experiment::ParameterSettingRow> rows;
  for (const auto &entry : parameterSettings_)
    if (entry.second.parameterNameForm) rows.push_back({entry.second.name, entry.second.line});
  std::sort(rows.begin(), rows.end(), [](const azure2::experiment::ParameterSettingRow &a,
                                         const azure2::experiment::ParameterSettingRow &b) { return a.line < b.line; });
  for (const azure2::experiment::Diagnostic &d : azure2::experiment::UnknownParameterSettings(config_->configfile, rows, known))
    config_->outStream << "WARNING: " << d.ToText() << std::endl;
}

ParameterSetting *ParameterLimitsManager::SettingForIndex(int nonFixedIndex) const {
  if (nonFixedIndex < 0 || nonFixedIndex >= (int)indexToSetting_.size()) return NULL;
  return indexToSetting_[nonFixedIndex];
}

void ParameterLimitsManager::ApplyAllParameterSettings(AZUREParams &params) {
  const ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // Build mapping from non-fixed parameter index to actual parameter index (as BuildIndexMap)
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < (int)p.Params().size(); i++) {
    if (!p.Parameter(i).IsFixed() || params.Record(i).Row().segmentField) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  BuildIndexMap(params);
  WarnUnknownSettings(params);

  // Apply settings using non-fixed parameter indices
  for (int nonFixedIndex = 0; nonFixedIndex < (int)nonFixedToActualIndex.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    ApplyParameterSettingByIndex(nonFixedIndex, actualIndex, params);
  }
}

void ParameterLimitsManager::ApplyParameterSettingByIndex(int nonFixedIndex, int actualIndex, AZUREParams &params) {
  ParameterSetting *setting = SettingForIndex(nonFixedIndex);

  if (!setting) {
    return;  // No settings found for this parameter index
  }
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // A level width's settings are partial widths (eV), converted to the reduced-width amplitude the fit varies.
  const bool isWidth = params.Kind(actualIndex) == ParamKind::Gamma && setting->category == "level";

  // Apply limits if they are not both zero (0,0 means free parameter)
  if (setting->lowerLimit != 0.0 || setting->upperLimit != 0.0) {
    double lower = setting->lowerLimit;
    double upper = setting->upperLimit;

    // For width parameters (physical values need conversion to reduced)
    if (isWidth) {
      lower = ConvertPhysicalLimitToReduced(setting->lowerLimit, actualIndex);
      upper = ConvertPhysicalLimitToReduced(setting->upperLimit, actualIndex);
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

  // A nuisance prior on a width is given as a partial width (eV) and converted to the reduced-width amplitude; on any
  // other parameter it is used as given (NuisancePrior reads the reduced values of a width only).
  if (setting->useAsNuisance && isWidth) {
    setting->nominalValueReduced = ConvertPhysicalLimitToReduced(setting->nominalValue, actualIndex);
    double errorSmall = std::abs(setting->nominalValueReduced - ConvertPhysicalLimitToReduced(setting->nominalValue - setting->error, actualIndex));
    double errorLarge = std::abs(ConvertPhysicalLimitToReduced(setting->nominalValue + setting->error, actualIndex) - setting->nominalValueReduced);
    setting->errorReduced = std::max(errorSmall, errorLarge);
  } else if (setting->useAsNuisance) {
    setting->nominalValueReduced = setting->nominalValue;
    setting->errorReduced = setting->error;
  }
}

double ParameterLimitsManager::ConvertPhysicalLimitToReduced(double physicalLimit, int paramIndex) const {
  // On copies of the caller's levels, filled with the current parameters.
  LevelState state(state_);
  LevelInputs inputs(inputs_);
  compound_->FillCompoundFromParams(params_->GetMinuitParams().Params(), state);
  compound_->CalcShiftFunctions(*config_, state);
  const ObservedLevelTable observed = compound_->TransformOut(*config_, state, inputs);

  try {
    vector_r currentParams = compound_->GetTransformParams(state, observed);
    // Only the level parameters (and abundances) are in the physical vector; the callers convert widths only.
    if (paramIndex < 0 || paramIndex >= (int)currentParams.size()) return physicalLimit;

    currentParams[paramIndex] = physicalLimit;
    compound_->FillCompoundFromParamsPhysical(currentParams, state, inputs);

    if (!compound_->TransformIn(*config_, state, inputs)) return physicalLimit;

    AZUREParams tempParams;
    compound_->FillMnParams(tempParams, state, inputs, config_);

    vector_r reducedParams = tempParams.GetMinuitParams().Params();

    if (paramIndex < (int)reducedParams.size()) {
      double convertedLimit = reducedParams[paramIndex];
      return convertedLimit;
    } else {
      return physicalLimit;
    }
  } catch (...) {
    return physicalLimit;
  }
}

bool ParameterLimitsManager::IsNuisanceParameterByIndex(int nonFixedIndex) const {
  const ParameterSetting *setting = SettingForIndex(nonFixedIndex);
  return setting ? setting->useAsNuisance : false;
}

bool ParameterLimitsManager::NuisancePriorPhysical(int actualIndex, double &nominal, double &error) const {
  if (actualIndex < 0 || actualIndex >= (int)priorByActualIndex_.size()) return false;
  const ParameterSetting *setting = priorByActualIndex_[actualIndex];
  if (!setting) return false;
  nominal = setting->nominalValue;
  error = setting->error;
  return error > 0.0;
}

bool ParameterLimitsManager::NuisancePrior(int actualIndex, double &nominal, double &error) const {
  if (actualIndex < 0 || actualIndex >= (int)priorByActualIndex_.size()) return false;
  const ParameterSetting *setting = priorByActualIndex_[actualIndex];
  if (!setting) return false;
  const bool isWidth = priorIsWidth_[actualIndex];
  nominal = isWidth ? setting->nominalValueReduced : setting->nominalValue;
  error = isWidth ? setting->errorReduced : setting->error;
  return error > 0.0;
}
