#include <algorithm>
#include "AZUREAPI.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include "ParameterLabel.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "CoulFunc.h"
#include "ECIntegral.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>

AZUREAPI::~AZUREAPI() {
  // A session used to be a process, so the OS reclaimed these; in-process the
  // compound nucleus and its data are ~10 MB per model and must be returned.
  if (compound_ != nullptr) delete compound_;
  if (data_ != nullptr) delete data_;
}

int AZUREAPI::Initialize() {
  // Initialize caches for performance
  InitializeCoulFuncCache();
  InitializeECAmplitudeCache();

  configure().paramMask |= Config::USE_EXTERNAL_CAPTURE;

  std::string file;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA)
    file = configure().outputdir + "intEC.dat";
  else
    file = configure().outputdir + "intEC.extrap";

  std::ifstream in(file.c_str());
  if (!in)
    configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  else
    configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  configure().integralsfile = file;

  // Initialize EC Integral caching system

  std::string cacheFile;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
    cacheFile = configure().outputdir + "intEC_cache.dat";
  } else {
    cacheFile = configure().outputdir + "intEC_cache.extrap";
  }

  // Initialize() is called again on every mode switch, so the previous model
  // has to go or it is leaked.  (This used to crash: the constructor left
  // data_/compound_ uninitialized, so the first call deleted a garbage
  // pointer.  They are null-initialized now.)  Same order as SetRadius.
  if (compound_ != nullptr) {
    delete compound_;
    compound_ = nullptr;
  }
  if (data_ != nullptr) {
    delete data_;
    data_ = nullptr;
  }

  data_ = new EData();
  compound_ = new CNuc();

  // configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if (compound()->Fill(configure()) == -1) {
    // configure().outStream << "Could not fill compound nucleus from file." << std::endl;
    return -1;
  } else if (compound()->NumPairs() == 0 || compound()->NumJGroups() == 0) {
    // configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl;
    return -1;
  }
  if ((configure().screenCheckMask | configure().fileCheckMask) &
      Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if (!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    // Fill the data object from the segments and data file
    //   Compound object is passed to the function for pair key verification and
    //   center of mass conversions, s-factor conversions, etc.
    // configure().outStream << "Filling Data Structures..." << std::endl;
    if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if (data()->Fill(configure(), compound()) == -1) {
        // configure().outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (data()->NumSegments() == 0) {
        // configure().outStream << "There is no data provided." << std::endl;
        return -1;
      }
    } else {
      if (data()->MakePoints(configure(), compound()) == -1) {
        // configure().outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (data()->NumSegments() == 0) {
        // configure().outStream << "Extrapolation segments produce no data." << std::endl;
        return -1;
      }
    }
    if ((configure().fileCheckMask | configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());
  } else {
    if (!compound()->IsPairKey(configure().rateParams.entrancePair) || !compound()->IsPairKey(configure().rateParams.exitPair)) {
      // configure().outStream << "Reaction rate pairs do not exist in compound nucleus." << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  // Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << std::endl
                          << "Calculation was aborted." << std::endl;
    return -1;
  }

  UpdateParameters();

  if (data()->Initialize(compound(), configure()) == -1) return -1;

  return 0;
}

bool AZUREAPI::UpdateParameters() {
  all_.clear();
  all_rwa_.clear();
  names_.clear();
  fixed_.clear();
  values_.clear();
  values_rwa_.clear();

  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams(), &configure());
  data()->FillMnParams(params.GetMinuitParams());

  compound()->FillCompoundFromParams(params.GetMinuitParams().Params());

  compound()->CalcShiftFunctions(configure());
  compound()->TransformOut(configure());

  for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
    // if( !params.GetMinuitParams().Parameter(i).IsFixed( ) ){
    names_.push_back(params.GetMinuitParams().Parameter(i).GetName());
    //}
    all_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    fixed_.push_back(params.GetMinuitParams().Parameter(i).IsFixed());
    if (!fixed_.back()) {
      values_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    }
  }

  all_ = compound()->GetTransformParams(configure());

  // Norms and shifts are missing from all_ and need to be added using the all_rwa_ values
  for (int i = all_.size(); i < all_rwa_.size(); ++i) {
    all_.push_back(all_rwa_[i]);
  }

  for (int i = 0; i < all_.size(); ++i) {
    // std::cout << "all_[" << i << "] = " << all_[i] << std::endl;
    // std::cout << "fixed_[" << i << "] = " << fixed_[i] << std::endl;
    // std::cout << "names_[" << i << "] = " << names_[i] << std::endl;
    if (!fixed_[i]) {
      values_.push_back(all_[i]);
    }
  }

  return true;
}

int AZUREAPI::UpdateSegments(vector_r &p) {
  calculatedConv_.clear();
  calculatedEnergies_.clear();
  calculatedAngles_.clear();
  calculatedAngularDists_.clear();
  calculatedSegments_.clear();
  calculatedSegmentsE1_.clear();
  calculatedSegmentsE2_.clear();
  calculatedExcitationEnergies_.clear();

  int k = 0;
  vector_r params_ = all_;
  for (int i = 0; i < all_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = p[k];
      ++k;
    }
  }

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(params_);
  bool isValid = localCompound->TransformIn(configure());

  if (!isValid) return 0;

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams().Params());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment> &segments = localData->GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    std::vector<EPoint> &data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for (int k = 0; k < data.size(); ++k) {
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back(theoreticalValue);
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    } else {
      // Regular segment calculation (existing logic)
      for (int k = 0; k < data.size(); ++k) {
        if (!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound, configure());
          } catch (GSLException &e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back(data[k].GetFitCrossSection());
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    }

    calculatedConv_.push_back(conv);
    calculatedSegments_.push_back(cross);
    calculatedSegmentsE1_.push_back(crossE1);
    calculatedSegmentsE2_.push_back(crossE2);
    calculatedEnergies_.push_back(energies);
    calculatedAngles_.push_back(angles);
    calculatedAngularDists_.push_back(angularDists);
    calculatedExcitationEnergies_.push_back(excitationEnergies);
  }

  return calculatedSegments_.size();
}

int AZUREAPI::UpdateSegmentsRWA(vector_r &p) {
  calculatedConv_.clear();
  calculatedEnergies_.clear();
  calculatedAngles_.clear();
  calculatedAngularDists_.clear();
  calculatedSegments_.clear();
  calculatedSegmentsE1_.clear();
  calculatedSegmentsE2_.clear();
  calculatedExcitationEnergies_.clear();

  int k = 0;
  vector_r params_ = all_rwa_;
  for (int i = 0; i < all_rwa_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = p[k];
      ++k;
    }
  }

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment> &segments = localData->GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    std::vector<EPoint> &data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for (int k = 0; k < data.size(); ++k) {
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back(theoreticalValue);
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    } else {
      // Regular segment calculation (existing logic)
      for (int k = 0; k < data.size(); ++k) {
        if (!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound, configure());
          } catch (GSLException &e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back(data[k].GetFitCrossSection());
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    }

    calculatedConv_.push_back(conv);
    calculatedSegments_.push_back(cross);
    calculatedSegmentsE1_.push_back(crossE1);
    calculatedSegmentsE2_.push_back(crossE2);
    calculatedEnergies_.push_back(energies);
    calculatedAngles_.push_back(angles);
    calculatedAngularDists_.push_back(angularDists);
    calculatedExcitationEnergies_.push_back(excitationEnergies);
  }

  return calculatedSegments_.size();
}

int AZUREAPI::UpdateSegmentsAllRWA(vector_r &p) {
  calculatedConv_.clear();
  calculatedEnergies_.clear();
  calculatedAngles_.clear();
  calculatedAngularDists_.clear();
  calculatedSegments_.clear();
  calculatedSegmentsE1_.clear();
  calculatedSegmentsE2_.clear();
  calculatedExcitationEnergies_.clear();

  int k = 0;
  for (int i = 0; i < all_rwa_.size(); ++i) {
    all_rwa_[i] = p[i];
    if (!fixed_[i]) {
      values_[k] = p[i];
      ++k;
    }
  }
  vector_r params_ = all_rwa_;

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment> &segments = localData->GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    std::vector<EPoint> &data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for (int k = 0; k < data.size(); ++k) {
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back(theoreticalValue);
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    } else {
      // Regular segment calculation (existing logic)
      for (int k = 0; k < data.size(); ++k) {
        if (!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound, configure());
          } catch (GSLException &e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back(data[k].GetFitCrossSection());
        crossE1.push_back(data[k].GetFitE1CrossSection());
        crossE2.push_back(data[k].GetFitE2CrossSection());
        angles.push_back(data[k].GetCMAngle());
        energies.push_back(data[k].GetCMEnergy());
        conv.push_back(data[k].GetSFactorConversion());
        excitationEnergies.push_back(data[k].GetExcitationEnergy());
        angularDists.push_back((double)data[k].GetNumAngularDists());
        for (int a = 0; a < data[k].GetNumAngularDists(); ++a)
          angularDists.push_back(data[k].GetAngularDist(a));
      }
    }

    calculatedConv_.push_back(conv);
    calculatedSegments_.push_back(cross);
    calculatedSegmentsE1_.push_back(crossE1);
    calculatedSegmentsE2_.push_back(crossE2);
    calculatedEnergies_.push_back(energies);
    calculatedAngles_.push_back(angles);
    calculatedAngularDists_.push_back(angularDists);
    calculatedExcitationEnergies_.push_back(excitationEnergies);
  }

  return calculatedSegments_.size();
}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformRWAParameters(const vector_r &p) const {
  // p holds the non-fixed parameters only, so the loop has to run over the
  // full parameter array -- running it to p.size() leaves every parameter
  // beyond that index at its .azr value, which silently returns stale widths
  // for the levels at the end of the file.  A shorter p (e.g. only the
  // leading R-matrix block, with the norm/shift tail sliced off) is still
  // accepted: the free parameters it does not cover keep their all_rwa_ value.
  vector_r params = all_rwa_;
  int k = 0;
  for (int i = 0; i < params.size(); ++i) {
    if (!fixed_[i]) {
      if (k < p.size()) params[i] = p[k];
      ++k;
    }
  }

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  localCompound->FillCompoundFromParams(params);

  localCompound->TransformOut(configure());

  vector_r transformedParams = compound()->GetTransformParams(configure());

  // Get only non fixed parameters
  vector_r transformed;
  k = 0;
  for (int i = 0; i < transformedParams.size(); ++i) {
    if (!fixed_[i]) {
      transformed.push_back(transformedParams[i]);
      ++k;
    }
  }

  return transformed;
}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformAllRWAParameters(const vector_r &p, bool includeFixed) const {
  // p holds every parameter, fixed ones included; a short p updates a prefix.
  vector_r params = all_rwa_;
  for (int i = 0; i < p.size() && i < params.size(); ++i) {
    params[i] = p[i];
  }

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  localCompound->FillCompoundFromParams(params);

  localCompound->TransformOut(configure());

  vector_r transformedParams = compound()->GetTransformParams(configure());

  if (includeFixed) return transformedParams;
  // Get only non fixed parameters
  vector_r transformed;
  for (int i = 0; i < transformedParams.size(); ++i) {
    if (!fixed_[i]) {
      transformed.push_back(transformedParams[i]);
    }
  }

  return transformed;
}

bool AZUREAPI::CalculateExternalCapture() {
  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  data()->CalculateECAmplitudes(compound(), configure());
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  return true;
}

int AZUREAPI::UpdateData() {
  dataEnergies_.clear();
  dataAngles_.clear();
  dataSegments_.clear();
  dataSegmentsErrors_.clear();
  dataConv_.clear();
  dataExcitationEnergies_.clear();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment> &segments = localData->GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    std::vector<EPoint> &data = segments[i].GetPoints();

    std::vector<double> energies, angles, cross, crossErr, conv, excitationEnergies;
    for (int k = 0; k < data.size(); ++k) {
      energies.push_back(data[k].GetCMEnergy());
      angles.push_back(data[k].GetCMAngle());
      cross.push_back(data[k].GetCMCrossSection());
      crossErr.push_back(data[k].GetCMCrossSectionError());
      conv.push_back(data[k].GetSFactorConversion());
      excitationEnergies.push_back(data[k].GetExcitationEnergy());
    }

    dataEnergies_.push_back(energies);
    dataSegments_.push_back(cross);
    dataAngles_.push_back(angles);
    dataSegmentsErrors_.push_back(crossErr);
    dataConv_.push_back(conv);
    dataExcitationEnergies_.push_back(excitationEnergies);
  }

  return nSegments;
}

void AZUREAPI::UpdateNorms() {
  norms_.clear();
  normsErrors_.clear();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment> &segments = localData->GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    double norm = segments[i].GetNominalNorm();
    double normErr = segments[i].GetNormError();

    norms_.push_back(norm);
    normsErrors_.push_back(normErr);
  }
}

// Set AZURE2 to calculate data points
void AZUREAPI::SetData() {
  configure().paramMask |= Config::CALCULATE_WITH_DATA;
}

// Set AZURE2 to calculate extrapolations
void AZUREAPI::SetExtrap() {
  configure().paramMask &= ~Config::CALCULATE_WITH_DATA;
}

// Set the channel radius of one particle pair (1-based) and rebuild.
//
// The radius enters the penetrabilities, shift functions, boundary conditions,
// Wigner limits, the ANC <-> reduced-width conversion and the lower limit of
// every external-capture integral, so nothing can be reused: the compound
// nucleus and data objects are rebuilt from the .azr with the override, the
// EC-integral cache is bypassed (USE_PREVIOUS_INTEGRALS cleared) so the
// integrals are recomputed on the new radius rather than read back from
// output/intEC*, and the parameter bookkeeping -- values, names, fixed flags
// and the transformed physical parameters -- is refilled to match.
bool AZUREAPI::SetRadius(int idx, double r) {
  return RebuildImpl(&idx, &r);
}

bool AZUREAPI::WriteOutputFiles() {
  if (data_ == nullptr) return false;
  data()->WriteOutputFiles(configure(), false, nullptr);
  return true;
}

std::vector<ThmExperimentReport> AZUREAPI::GetThmExperiments() const {
  if (data_ == nullptr) return std::vector<ThmExperimentReport>();
  return data()->ThmExperimentReports();
}

bool AZUREAPI::GetThmLineshape(const std::string &name, const vector_r &energies, ThmLineshapeReport &out,
                               std::string &why) {
  if (data_ == nullptr || compound_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmLineshapeTable(name, energies, compound(), configure(), out, why);
}

bool AZUREAPI::GetThmVertex(const std::string &name, const vector_r &energies, ThmVertexReport &out,
                            std::string &why) {
  if (data_ == nullptr || compound_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmVertexTable(name, energies, compound(), configure(), out, why);
}

bool AZUREAPI::GetThmDistortion(const std::string &name, const vector_r &energies, ThmDistortionReport &out,
                                std::string &why) {
  if (data_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmDistortionTable(name, energies, out, why);
}

bool AZUREAPI::Rebuild() {
  return RebuildImpl(nullptr, nullptr);
}

// Shared body of SetRadius and Rebuild.  Everything the two do is identical
// apart from whether CNuc::Fill is handed a radius override.
bool AZUREAPI::RebuildImpl(const int *idx, const double *r) {
  if (compound_ != nullptr) delete compound_;
  if (data_ != nullptr) delete data_;

  compound_ = new CNuc;
  data_ = new EData;

  if (idx && r) {
    std::pair<int, double> pair = std::make_pair(*idx, *r);
    if (compound()->Fill(configure(), pair) == -1) return false;
  } else {
    if (compound()->Fill(configure()) == -1) return false;
  }
  if (compound()->NumPairs() == 0 || compound()->NumJGroups() == 0) return false;

  // Mirror the branch taken at startup: data mode reads the data segments,
  // extrapolation mode builds its points from <segmentsTest>.
  if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
    if (data()->Fill(configure(), compound()) == -1) return false;
  } else {
    if (data()->MakePoints(configure(), compound()) == -1) return false;
  }
  if (data()->NumSegments() == 0) return false;

  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    configure().outStream << e.what() << std::endl;
    return false;
  }
  if (data()->Initialize(compound(), configure()) == -1) {
    configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    return false;
  }
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  // Without this the client keeps the pre-change parameter vector, names,
  // fixed flags and Wigner limits -- silently mismatched with the new model.
  UpdateParameters();

  return true;
}

// Map a packed non-fixed RWA parameter vector into the full all_rwa_ layout.
static vector_r MapPackedToFull(const vector_r &packed, const vector_r &all_rwa,
                                const std::vector<bool> &fixed) {
  vector_r full = all_rwa;
  int k = 0;
  for (int i = 0; i < (int)all_rwa.size(); ++i)
    if (!fixed[i] && k < (int)packed.size()) full[i] = packed[k++];
  return full;
}

void AZUREAPI::FillFromFullRWA(const vector_r &full) const {
  CNuc *lc = compound();
  EData *ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());
}

// The data term of AZURECalc::operator(), on the filled canonical objects.
double AZUREAPI::EvaluateFilledChi2(vector_r *res) const {
  CNuc *lc = compound();
  EData *ld = data();
  if (res) res->clear();
  double chiSquared = 0.0;
  // Rows of the segments of a THM experiment, filled once it is profiled.
  std::vector<std::vector<std::pair<int, size_t>>> groupRows(ld->NumThmGroups());
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *segment = ld->GetSegment(i);
    if (!segment) continue;
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
      EPoint *point = segment->GetPoint(pid + 1);
      if (point) point->SetFitCrossSection(th);
    }

    // A THM experiment: one profile (shared norm, background) over its
    // segments, at the last of them, as the CLI.
    const int g = ld->ThmGroupOf(i);
    if (g >= 0) {
      if (res) {
        groupRows[g].push_back(std::make_pair(i, res->size()));
        for (int pid = 0; pid < segment->NumPoints(); pid++)
          if (segment->GetPoint(pid + 1)) res->push_back(0.0);
      }
      if (ld->IsLastOfThmGroup(g, i)) {
        chiSquared += ld->ProfileThmGroup(g);
        if (res) {
          std::vector<double> r;
          for (const std::pair<int, size_t> &row : groupRows[g]) {
            ld->ThmGroupResiduals(g, row.first, r);
            std::copy(r.begin(), r.end(), res->begin() + row.second);
          }
        }
      }
      continue;
    }

    // A free THM norm is the arbitrary HOES scale: profiled to its optimum
    // (which also stores n* on the segment), exactly as the CLI does.
    double segChi = 0.0;
    const bool profiled = segment->IsProfiledNorm();
    if (profiled) segChi = segment->ProfileNormChiSquared();

    const double norm = segment->GetNorm();
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      EPoint *point = segment->GetPoint(pid + 1);
      if (!point) continue;
      double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * norm;
      double error = point->GetCMCrossSectionError() * norm;
      double r = (error != 0.0) ? residual / error : 0.0;
      if (!profiled) segChi += r * r;
      if (res) res->push_back(r);
    }

    segment->SetSegmentChiSquared(segChi);
    chiSquared += segChi;
  }
  return chiSquared;
}

double AZUREAPI::CalculateChi2RWA(const vector_r &rwaParams) const {
  // One request at a time: operate on the canonical compound/data in place
  // (re-filled here from these parameters) instead of cloning.
  FillFromFullRWA(MapPackedToFull(rwaParams, all_rwa_, fixed_));
  double chiSquared = EvaluateFilledChi2(nullptr);
  // Park formalism: the J > 0 wall, part of what a fit must minimize.
  if (configure().paramMask & Config::USE_PARK_FORMALISM) chiSquared += compound()->ParkNormPenalty();
  return chiSquared;
}

vector_r AZUREAPI::CalculateResidualsRWA(const vector_r &params) const {
  FillFromFullRWA(MapPackedToFull(params, all_rwa_, fixed_));
  vector_r res;
  EvaluateFilledChi2(&res);
  return res;
}

vector_r AZUREAPI::GetCurrentNorms() const {
  vector_r out;
  int prevKey = -1;
  std::vector<ESegment> &segments = data()->GetSegments();
  for (size_t i = 0; i < segments.size(); ++i) {
    int newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    out.push_back(segments[i].GetNorm());
  }
  return out;
}

vector_r AZUREAPI::ParkNorms(const vector_r &rwaParams) const {
  vector_r params_ = all_rwa_;
  for (int i = 0, k = 0; i < (int)all_rwa_.size(); ++i)
    if (!fixed_[i] && k < (int)rwaParams.size()) params_[i] = rwaParams[k++];
  CNuc *localCompound = compound();
  localCompound->FillCompoundFromParams(params_);
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  vector_r out;
  for (int j = 1; j <= localCompound->NumJGroups(); j++) {
    JGroup *jg = localCompound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++)
      if (jg->GetLevel(la)->IsInRMatrix()) out.push_back(jg->GetLevel(la)->GetParkNorm());
  }
  return out;
}

double AZUREAPI::CalculateChi2Physical(const vector_r &physicalParams) const {
  int k = 0;
  vector_r params_ = all_;
  for (int i = 0; i < all_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = physicalParams[k];
      ++k;
    }
  }

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled below from these parameters) instead of cloning.
  CNuc *localCompound = compound();
  EData *localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(params_);
  bool isValid = localCompound->TransformIn(configure());

  if (!isValid) return 0;

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());
  // Norms and shifts sit at the same offsets in the physical vector as in the
  // Minuit one (all_ is the transformed R-matrix block plus all_rwa_'s tail).
  localData->FillNormsFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams().Params());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  return EvaluateFilledChi2(nullptr);
}

std::vector<AZUREAPI::THMRows> AZUREAPI::ComputeTHMRows(const vector_r &full,
                                                        const ParamIndexMap &pmap) const {
  return ::ComputeTHMRows(compound(), data(), configure(), full, pmap);
}


bool AZUREAPI::Chi2GradEGammaNorm(const vector_r &full, vector_r &gradFull,
                                  double &chi2Out) const {
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  // Only one API request runs at a time, so operate on the canonical
  // compound/data in place (like UpdateSegments) rather than cloning; every
  // entry point re-fills from its own parameters before use.
  CNuc *lc = compound();
  EData *ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (brune) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  GradAccum accum;
  accum.Init(lc);

  // chi2 = sum (fit - data*n)^2/(cmErr*n)^2.  d(chi2)/d(model) = 2 r / err^2;
  // the same model gives the data-term norm gradient, accumulated per segment.
  std::vector<double> normData(ld->NumSegments() + 1, 0.0);
  double chi2 = 0.0;
  FitBarFn fb = [&](ESegment *seg, int i, int pid, double model) -> double {
    EPoint *pt = seg->GetPoint(pid + 1);
    if (!pt) return 0.0;
    double norm = seg->GetNorm();
    double dataval = pt->GetCMCrossSection();
    double cmErr = pt->GetCMCrossSectionError();
    double r = model - dataval * norm;
    double err = cmErr * norm;
    if (err == 0.0) return 0.0;
    // chi2 falls out of the same residual the gradient uses (so value and
    // gradient stay consistent). Runs in the parallel point loop, so guard it.
#pragma omp atomic
    chi2 += (r * r) / (err * err);
    if (seg->IsVaryNorm() && norm != 0.0 && i >= 1 && i < (int)normData.size()) {
      double e2 = cmErr * cmErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[i], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm) - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[i] += dNorm;
    }
    return 2.0 * r / (err * err);
  };

  bool ok = AccumulateEGammaGradient(lc, ld, configure(), pmap, sdp, fb, accum);
  if (ok) {
    if (configure().paramMask & Config::USE_PARK_FORMALISM) {
      AddParkPenaltyGradient(lc, configure(), accum);
      chi2 += lc->ParkNormPenalty();
    }
    accum.Scatter(pmap, gradFull);
    for (int s = 1; s <= ld->NumSegments(); s++) {
      ESegment *seg = ld->GetSegment(s);
      if (!seg || !seg->IsVaryNorm()) continue;
      int idx = pmap.NormIndex(s);
      if (idx >= 0 && idx < (int)gradFull.size()) gradFull[idx] = normData[s];
    }
    chi2Out = chi2;
  }
  return ok;
}

vector_r AZUREAPI::CalculateChi2GradRWA(const vector_r &params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // The analytic gradient pass runs the full forward model at every point, so
  // chi2 comes out as a byproduct -- no separate forward pass. Only fall back to
  // a standalone chi2 evaluation if the analytic path bails.
  vector_r gradFull(all_rwa_.size(), 0.0);
  double chi2 = 0.0;
  bool eg = Chi2GradEGammaNorm(full, gradFull, chi2);
  if (!eg) chi2 = CalculateChi2RWA(params);

  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed_);

  // THM (HOES) segments are skipped by the adjoint, so neither their chi2 nor
  // their gradient is in what it returned.  Add both: the chi2 with the norm
  // profiled as the CLI does, the gradient as 2 J^T r over the THM rows (whose
  // J carries the dependence of the profiled scale; at the profile optimum
  // that term drops out of the gradient anyway, sum_i r_i m_i/e_i = 0).
  if (eg) {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    const int nCols = pmap.NumPacked();
    for (const THMRows &tr : thm) {
      ESegment *seg = data()->GetSegment(tr.segment);
      double c = 0.0;
      if (data()->ThmGroupOf(tr.segment) >= 0) {
        // A THM experiment, profiled by ComputeTHMRows: its residuals.
        for (double r : tr.r) c += r * r;
      } else if (seg->IsProfiledNorm()) {
        c = seg->ProfileNormChiSquared();
      } else {
        for (double r : tr.r) c += r * r;
      }
      seg->SetSegmentChiSquared(c);
      chi2 += c;
      for (int col = 0; col < nCols; col++) {
        double g = 0.0;
        for (size_t i = 0; i < tr.r.size(); i++) g += tr.r[i] * tr.J[i * nCols + col];
        gradFull[pmap.PackedToFull(col)] += 2.0 * g;
      }
    }
  }

  // Finite differences for energy shifts (and the whole block if the analytic
  // path bailed), using the scalar chi-squared.
  for (int f = 0; f < (int)full.size() && f < pmap.NumFull(); f++) {
    if (fixed_[f]) continue;
    ParamKind kind = pmap.Desc(f).kind;
    // THM coherent backgrounds (cbackground=) are columns of the THM rows above.
    if (eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::Norm ||
               kind == ParamKind::ThmCoherent))
      continue;
    int packed = pmap.FullToPacked(f);
    if (packed < 0 || packed >= (int)params.size()) continue;
    double x0 = params[packed];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = params;
    pp[packed] = x0 + h;
    vector_r pm = params;
    pm[packed] = x0 - h;
    gradFull[f] = (CalculateChi2RWA(pp) - CalculateChi2RWA(pm)) / (2.0 * h);
  }

  vector_r out;
  out.reserve(1 + params.size());
  out.push_back(chi2);
  for (int i = 0; i < (int)all_rwa_.size(); i++)
    if (!fixed_[i]) out.push_back(gradFull[i]);
  return out;
}

vector_r AZUREAPI::CalculateResidualJacobianRWA(const vector_r &params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled here from these parameters) instead of cloning.
  CNuc *lc = compound();
  EData *ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  vector_r residuals, jacobian;
  int nCols = 0;
  // THM rows are left zero by the adjoint and filled here: residuals with the
  // profiled norm, Jacobian by central differences of the HOES model plus the
  // analytic dependence of the profiled scale (see ComputeTHMRows).
  bool ok = ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, residuals, jacobian, nCols,
                                    /*skipTHM=*/true);

  if (!ok) return vector_r{-1.0};

  {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    for (const THMRows &tr : thm) {
      for (size_t i = 0; i < tr.r.size(); i++) {
        const size_t row = (size_t)tr.firstRow + i;
        if (row >= residuals.size()) continue;
        residuals[row] = tr.r[i];
        for (int c = 0; c < nCols; c++) jacobian[row * (size_t)nCols + c] = tr.J[i * nCols + c];
      }
    }
  }

  // ---- Energy-shift columns, by central differences -----------------------
  //
  // The adjoint covers level energies, reduced widths and normalizations. An
  // energy shift is a different derivative: it translates the energy axis of a
  // whole segment, so what is wanted is d(model)/dE rather than d(model)/d(a
  // parameter). Every energy-dependent quantity moves with it --- the level
  // matrix, the penetrabilities and shift functions, the Coulomb and
  // hard-sphere phases, the external-capture amplitudes, and, where a segment
  // carries target integration, the sub-point grid itself (which brings
  // boundary terms with it). AZURE2 applies a shift by rebuilding all of that
  // in UpdatePointEnergiesWithShift, so there is no cheap analytic route
  // through the existing forward code.
  //
  // Returning these columns as zero, which is what this function used to do,
  // is the worst of the options: a least-squares driver handed a zero column
  // simply never moves that parameter, converges, and reports success. Finite
  // differences cost two residual evaluations per free shift and are correct.
  //
  // This mirrors what CalculateChi2GradRWA already does for the scalar
  // gradient. The differenced residuals are the forward ones of
  // EvaluateFilledChi2, which emits one row per point in the order
  // ComputeResidualJacobian assigns them (THM rows with their profiled norm).
  {
    const size_t nRes = residuals.size();
    vector_r rPlus, rMinus;
    for (int f = 0; f < pmap.NumFull(); f++) {
      if (fixed_[f]) continue;
      if (pmap.Desc(f).kind != ParamKind::EnergyShift &&
          pmap.Desc(f).kind != ParamKind::EnergyShiftSqrt) continue;
      const int packed = pmap.FullToPacked(f);
      if (packed < 0 || packed >= (int)params.size()) continue;

      const double x0 = params[packed];
      // A few eV: far inside any quoted beam-energy calibration uncertainty,
      // far outside the noise of the forward model.
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);

      // The forward residuals: the same rows as ComputeResidualJacobian's,
      // THM rows included with their profiled norm, and no adjoint to pay for.
      auto residualsAt = [&](double value, vector_r &out) -> bool {
        vector_r pk = params;
        pk[packed] = value;
        FillFromFullRWA(MapPackedToFull(pk, all_rwa_, fixed_));
        EvaluateFilledChi2(&out);
        return true;
      };

      const bool okp = residualsAt(x0 + h, rPlus);
      const bool okm = residualsAt(x0 - h, rMinus);
      if (okp && okm && rPlus.size() == nRes && rMinus.size() == nRes) {
        for (size_t r = 0; r < nRes; r++)
          jacobian[r * (size_t)nCols + (size_t)packed] =
              (rPlus[r] - rMinus[r]) / (2.0 * h);
      }
    }
    // Leave the model at the parameters that were asked for.
    FillFromFullRWA(full);
    EvaluateFilledChi2(nullptr);
  }

  vector_r out;
  out.reserve(2 + residuals.size() + jacobian.size());
  out.push_back((double)residuals.size());
  out.push_back((double)nCols);
  out.insert(out.end(), residuals.begin(), residuals.end());
  out.insert(out.end(), jacobian.begin(), jacobian.end());
  return out;
}

vector_r AZUREAPI::CalculateModelGradientsRWA(const vector_r &params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // Same preparation as CalculateResidualJacobianRWA: one request at a time, so
  // the canonical compound/data are re-filled in place from these parameters.
  CNuc *lc = compound();
  EData *ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);

  // A band is sensitive only to the R-matrix (and THM coherent-background)
  // parameters, and covariance.dat spans exactly those columns
  // (BandPackedColumns) -- so reduce each full packed row to them here
  // rather than shipping zero columns for every normalization.
  const std::vector<int> rc = BandPackedColumns(pmap);
  const int nCols = (int)rc.size();

  // THM (HOES) rows come from central differences of the HOES model (the
  // adjoint differentiates the T-matrix observable); see ComputeTHMRows.
  std::map<EPoint *, vector_r> grad;
  if (!ComputeModelGradients(lc, ld, configure(), pmap, sdp, grad, /*skipTHM=*/true))
    return vector_r{-1.0};
  {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    const int nPacked = pmap.NumPacked();
    for (const THMRows &tr : thm) {
      ESegment *seg = ld->GetSegment(tr.segment);
      size_t i = 0;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        EPoint *pt = seg->GetPoint(pid + 1);
        if (!pt) continue;
        std::map<EPoint *, vector_r>::iterator it = grad.find(pt);
        if (it != grad.end() && i < tr.m.size())
          for (int c = 0; c < nPacked && c < (int)it->second.size(); c++)
            it->second[c] = tr.Jm[i * nPacked + c];
        i++;
      }
    }
  }

  // Walk the segments exactly as UpdateSegments does, so row k of segment s
  // lines up with GET_CALCULATED_SEGMENT's point k: segments sharing a segment
  // key collapse to one calculated segment, and only the first of them is used.
  std::vector<int> counts;
  vector_r rows;
  int prevKey = -1;
  std::vector<ESegment> &segments = ld->GetSegments();
  for (int i = 0; i < (int)segments.size(); ++i) {
    const int newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;

    std::vector<EPoint> &points = segments[i].GetPoints();
    int n = 0;
    for (int k = 0; k < (int)points.size(); ++k) {
      std::map<EPoint *, vector_r>::const_iterator it = grad.find(&points[k]);
      if (it == grad.end()) continue;
      const vector_r &f = it->second;
      for (int c = 0; c < nCols; ++c)
        rows.push_back(rc[c] < (int)f.size() ? f[rc[c]] : 0.0);
      ++n;
    }
    counts.push_back(n);
  }

  vector_r out;
  out.reserve(2 + counts.size() + rows.size());
  out.push_back((double)counts.size());
  out.push_back((double)nCols);
  for (int i = 0; i < (int)counts.size(); ++i) out.push_back((double)counts[i]);
  out.insert(out.end(), rows.begin(), rows.end());
  return out;
}

// Function to get indeces of non-fixed normalization parameters
vector_r AZUREAPI::GetNormalizationIndices() {
  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for (int i = 0; i < totalParams; ++i) {
    if (!fixed_[i]) {
      // A segment norm (by form: a cbkg name can contain "norm")
      if (AZURELabel::IsNormName(names_[i])) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}

// Function to get indeces of non-fixed energy shift parameters
vector_r AZUREAPI::GetEnergyShiftIndices() {
  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for (int i = 0; i < totalParams; ++i) {
    if (!fixed_[i]) {
      // A segment energy shift (by form: a cbkg name can contain "shift")
      if (AZURELabel::IsEnergyShiftName(names_[i])) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}

// Structured metadata for every parameter.  The parameters are walked in the
// exact order CNuc::FillMnParams (energies + widths) then EData::FillMnParams
// (norms, energy shifts, sqrt(E) shift coefficients, then the THM coherent
// backgrounds) emit them, so the records line up one-to-one with
// names_ / all_ / fixed_.  See AZUREAPI.h for the field layout.
vector_r AZUREAPI::GetParameterInfo() const {
  vector_r info;

  // Running global parameter index; reads the fixed flag and physical value
  // that UpdateParameters() already cached.
  int gi = 0;
  auto push = [&](double type, double jgroup, double J, double parity,
                  double level, double levelE, double channel, double L,
                  double S, double pair, double radtype, double segKey,
                  double wignerLimit, double inputIsRWA) {
    info.push_back(type);
    info.push_back(jgroup);
    info.push_back(J);
    info.push_back(parity);
    info.push_back(level);
    info.push_back(levelE);
    info.push_back(channel);
    info.push_back(L);
    info.push_back(S);
    info.push_back(pair);
    info.push_back(radtype);
    info.push_back((gi < (int)fixed_.size() && fixed_[gi]) ? 1.0 : 0.0);
    info.push_back(gi < (int)all_.size() ? all_[gi] : 0.0);
    info.push_back(segKey);
    info.push_back(wignerLimit);
    info.push_back(inputIsRWA);
    ++gi;
  };

  // R-matrix parameters: one energy then one width per channel, per level.
  CNuc *nuc = compound();
  for (int j = 1; j <= nuc->NumJGroups(); ++j) {
    JGroup *jg = nuc->GetJGroup(j);
    double J = jg->GetJ();
    double parity = jg->GetPi();
    for (int la = 1; la <= jg->NumLevels(); ++la) {
      ALevel *level = jg->GetLevel(la);
      double levelE = level->GetE();
      // energy parameter
      push(0, j, J, parity, la, levelE, -1, -1, -1, -1, -1, -1, -1, -1);
      // width parameters (one per channel)
      for (int ch = 1; ch <= jg->NumChannels(); ++ch) {
        AChannel *chan = jg->GetChannel(ch);
        push(1, j, J, parity, la, levelE, ch, chan->GetL(), chan->GetS(),
             chan->GetPairNum(), (double)chan->GetRadType(), -1,
             chan->GetWignerLimit(),
             level->GammaIsRWA(ch) ? 1.0 : 0.0);
      }
    }
  }

  // Normalization parameters: one per segment with IsVaryNorm(), except a
  // profiled (free THM) norm, which EData::FillMnParams leaves out too.
  std::vector<ESegment> &segments = data()->GetSegments();
  for (size_t s = 0; s < segments.size(); ++s) {
    if (segments[s].IsVaryNorm() && !segments[s].IsProfiledNorm())
      push(2, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1,
           segments[s].GetSegmentKey(), -1, -1);
  }

  // Energy-shift parameters: one per segment (always emitted).
  for (size_t s = 0; s < segments.size(); ++s) {
    push(3, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1,
         segments[s].GetSegmentKey(), -1, -1);
  }

  // sqrt(E) energy-shift coefficients: one per segment (always emitted).
  for (size_t s = 0; s < segments.size(); ++s) {
    push(4, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1,
         segments[s].GetSegmentKey(), -1, -1);
  }

  // THM coherent backgrounds (cbackground=), last: Re c0, Im c0 [, Re c1,
  // Im c1] per combination, with its J group, entrance channel (channel, L,
  // S) and exit pair number (pair).
  for (int g = 0; g < data()->NumThmGroups(); ++g) {
    const ThmCoherentBackground *cb = data()->GetThmGroup(g).coherent.get();
    if (!cb) continue;
    for (const ThmCoherentBackground::Combo &c : cb->combos) {
      JGroup *jg = nuc->GetJGroup(c.jGroup);
      AChannel *in = jg->GetChannel(c.entrance);
      for (int k = 0; k < 2 * c.form; ++k)
        push(5, c.jGroup, jg->GetJ(), jg->GetPi(), -1, 0, c.entrance, in->GetL(), in->GetS(),
             jg->GetChannel(c.exit)->GetPairNum(), -1, -1, -1, -1);
    }
  }

  return info;
}

vector_r AZUREAPI::GetPairsInfo() const {
  vector_r info;
  CNuc *nuc = compound();

  // One record per pair, in 1-based pair-number order so field "pair" of
  // GetParameterInfo() indexes directly into this list.
  for (int p = 1; p <= nuc->NumPairs(); ++p) {
    PPair *pair = nuc->GetPair(p);
    info.push_back(p);
    info.push_back(pair->GetPairKey());
    info.push_back(pair->GetPType());
    info.push_back(pair->IsEntrance() ? 1.0 : 0.0);
    info.push_back(pair->GetJ(1));
    info.push_back(pair->GetPi(1));
    info.push_back(pair->GetZ(1));
    info.push_back(pair->GetM(1));
    info.push_back(pair->GetJ(2));
    info.push_back(pair->GetPi(2));
    info.push_back(pair->GetZ(2));
    info.push_back(pair->GetM(2));
    info.push_back(pair->GetSepE());
    info.push_back(pair->GetExE());
    info.push_back(pair->GetChRad());
    info.push_back(pair->GetI1I2Factor());
    info.push_back(pair->GetBindingEnergy());
  }

  return info;
}
// ---------------------------------------------------------------------------
//  Diagnostics: the external region, and the caches that make it affordable
// ---------------------------------------------------------------------------

/*!
 * Coulomb wave functions, penetrability, shift function and hard-sphere phase
 * on a requested energy grid.
 *
 * These are ordinary R-matrix quantities that AZURE2 has always computed
 * internally; what is new is being able to ask for them.  The values follow the
 * run's own configuration, so the same call returns the accurate Coulomb
 * routine's answer, GSL's, or the Numerov solution through a nuclear potential,
 * according to how the calculation was set up -- which is how one sees what the
 * hybrid model actually does to the external region.
 */
vector_r AZUREAPI::GetCoulombFunctions(const vector_r &request) const {
  vector_r out;
  if (request.size() < 4) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int lValue = static_cast<int>(std::lround(request[1]));
  double radius = request[2];
  int nE = static_cast<int>(std::lround(request[3]));
  if (nE < 0 || request.size() < static_cast<size_t>(4 + nE)) return out;

  if (!compound_ || !compound_->IsPairKey(pairKey)) return out;
  PPair *pair = compound_->GetPair(compound_->GetPairNumFromKey(pairKey));
  if (!pair) return out;
  if (radius <= 0.0) radius = pair->GetChRad();  // the channel radius by default

  CoulFunc coul(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC));

  out.push_back(static_cast<double>(nE));
  for (int i = 0; i < nE; ++i) {
    double energy = request[4 + i];
    double F = 0., dF = 0., G = 0., dG = 0., P = 0., S = 0., delta = 0.;
    if (energy > 0.0) {
      try {
        CoulWaves w = coul(lValue, radius, energy);
        F = w.F;
        dF = w.dF;
        G = w.G;
        dG = w.dG;
        P = coul.Penetrability(lValue, radius, energy);
        S = coul.PEShift(lValue, radius, energy);
        delta = -std::atan2(w.F, w.G);
      } catch (...) {
        // A failed evaluation returns zeros for that energy rather than
        // aborting the whole grid; the caller can see which points are missing.
        F = dF = G = dG = P = S = delta = 0.;
      }
    }
    out.push_back(F);
    out.push_back(dF);
    out.push_back(G);
    out.push_back(dG);
    out.push_back(P);
    out.push_back(S);
    out.push_back(delta);
  }
  return out;
}

/*!
 * External-capture radial integrals on a requested energy grid.
 *
 * Walks exactly the pathway structure EPoint::CalculateECAmplitudes walks --
 * every EC level whose final pair matches, every KGroup of the entrance pair's
 * decay, every ECMGroup within it -- and evaluates the integral at each
 * requested energy instead of at a data point's energy.  The quantum numbers of
 * each pathway come back with it, so no bookkeeping is needed on the caller's
 * side.
 */
vector_r AZUREAPI::GetECIntegrals(const vector_r &request) const {
  vector_r out;
  if (request.size() < 2) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int nE = static_cast<int>(std::lround(request[1]));
  if (nE <= 0 || request.size() < static_cast<size_t>(2 + nE)) return out;
  std::vector<double> energies(request.begin() + 2, request.begin() + 2 + nE);

  if (!compound_ || !compound_->IsPairKey(pairKey)) return out;
  int aa = compound_->GetPairNumFromKey(pairKey);
  PPair *entrancePair = compound_->GetPair(aa);
  if (!entrancePair || !entrancePair->IsEntrance()) return out;

  // One block per pathway: 6 descriptors then 2*nE numbers.
  std::vector<vector_r> blocks;

  for (int j = 1; j <= compound_->NumJGroups(); j++) {
    for (int la = 1; la <= compound_->GetJGroup(j)->NumLevels(); la++) {
      ALevel *ecLevel = compound_->GetJGroup(j)->GetLevel(la);
      if (!ecLevel->IsECLevel()) continue;
      int ir = ecLevel->GetECPairNum();
      if (ir < 1 || ir > compound_->NumPairs()) continue;

      for (int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
        KGroup *theKGroup = entrancePair->GetDecay(ir)->GetKGroup(k);
        for (int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
          ECMGroup *g = theKGroup->GetECMGroup(ecm);

          AChannel *finalChannel =
              compound_->GetJGroup(j)->GetChannel(g->GetFinalChannel());
          PPair *finalPair = compound_->GetPair(finalChannel->GetPairNum());

          int liValue;
          double siValue;
          if (g->IsChannelCapture()) {
            MGroup *cc = entrancePair->GetDecay(g->GetChanCapDecay())
                             ->GetKGroup(g->GetChanCapKGroup())
                             ->GetMGroup(g->GetChanCapMGroup());
            liValue = compound_->GetJGroup(cc->GetJNum())
                          ->GetChannel(cc->GetChpNum())
                          ->GetL();
            siValue = compound_->GetJGroup(cc->GetJNum())
                          ->GetChannel(cc->GetChpNum())
                          ->GetS();
          } else {
            liValue = g->GetL();
            siValue = theKGroup->GetS();
          }

          vector_r block;
          block.push_back(static_cast<double>(liValue));
          block.push_back(static_cast<double>(finalChannel->GetL()));
          block.push_back(2.0 * siValue);
          block.push_back(2.0 * finalChannel->GetS());
          block.push_back(static_cast<double>(g->GetMult()));
          block.push_back(g->GetRadType() == 'E' ? 1.0 : 0.0);

          ECIntegral theECIntegral(finalPair, configure());
          double levelEnergy = ecLevel->GetE();
          for (int i = 0; i < nE; ++i) {
            double inEnergy = energies[i] + entrancePair->GetSepE() + entrancePair->GetExE();
            complex v(0., 0.);
            if (energies[i] > 0.0) {
              try {
                v = theECIntegral(liValue, finalChannel->GetL(),
                                  siValue, finalChannel->GetS(),
                                  g->GetJ(), compound_->GetJGroup(j)->GetJ(),
                                  g->GetMult(), g->GetRadType(),
                                  inEnergy, levelEnergy,
                                  g->IsChannelCapture());
              } catch (...) {
                v = complex(0., 0.);
              }
            }
            block.push_back(real(v));
            block.push_back(imag(v));
          }
          blocks.push_back(std::move(block));
        }
      }
    }
  }

  out.push_back(static_cast<double>(blocks.size()));
  out.push_back(static_cast<double>(nE));
  for (const vector_r &b : blocks) out.insert(out.end(), b.begin(), b.end());
  return out;
}

/*!
 * Coulomb-function cache counters, aggregated over threads.
 */
vector_r AZUREAPI::GetCacheStats() const {
  vector_r out(6, 0.0);
  if (!g_coulFuncCache) return out;
  CoulFuncCache::Stats s = g_coulFuncCache->GetStats();
  out[0] = static_cast<double>(s.queries);
  out[1] = static_cast<double>(s.hits);
  out[2] = static_cast<double>(s.entries);
  out[3] = static_cast<double>(s.keys);
  out[4] = static_cast<double>(s.disabledKeys);
  out[5] = static_cast<double>(s.threads);
  return out;
}
