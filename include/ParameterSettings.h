#ifndef PARAMETERSETTINGS_H
#define PARAMETERSETTINGS_H

#include <string>
#include <Minuit2/MnUserParameters.h>

// Forward declarations
class Config;
class CNuc;
class EData;

// Function to apply parameter settings to Minuit parameter
void ApplyParameterSettings(const std::string& paramName, ROOT::Minuit2::MnUserParameters& p, 
                          const Config* config = nullptr, CNuc* compound = nullptr, EData* data = nullptr);

#endif