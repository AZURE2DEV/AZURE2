// AZUREAPI, what a session reports of a calculation: the arrays of the calculated points for a parameter vector in
// each basis (ReportSegments, UpdateSegments*), the data and the norms as set up (UpdateData, UpdateNorms,
// GetCurrentNorms), the output files, the tables of a THM experiment, and the model's own numbers at a vector
// (scattering lengths, Park's norms).  One unit because each fills the session's level state at a vector, runs the
// forward calculation and hands arrays to the binding; none of them is the objective.
#include <algorithm>
#include "Objective.h"
#include "ParameterLimitsManager.h"
#include "AZUREAPI.h"
#include "CompoundSystemSetUp.h"
#include "AdaptiveExtrapolation.h"
#include "ChannelFunc.h"
#include "ModelPass.h"
#include "ScatteringLength.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include <stdexcept>
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "PPair.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "PriorCovarianceFile.h"
#include "AZURECalc.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "Profile.h"
#include "InputFiles.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>
#include "AZUREAPIInternal.h"

// The segments the UpdateSegments* loops report: the first of each run of equal keys.
static std::vector<const Segment *> ReportedSegments(const std::vector<Segment> &segments) {
  std::vector<const Segment *> out;
  int prevKey = -1;
  for (const Segment &seg : segments) {
    if (seg.GetSegmentKey() == prevKey) continue;
    prevKey = seg.GetSegmentKey();
    out.push_back(&seg);
  }
  return out;
}


// The arrays UpdateSegments* report, for the segments of ReportedSegments: per point the model value
// (Results::model, the value the residual compares with the datum: the E1 / E2 part of an E1- or E2-only segment, the
// combination of a SUM, RATIO, sample or transform), the E1 and E2 parts and the angular-distribution coefficients
// of its base point, its c.m. energy and angle, S-factor conversion and excitation energy.  A point whose
// calculation fails is NaN (Unit::MarkFailed), and the arrays keep one entry per point.  Without data: the arrays
// are emptied.
int AZUREAPI::ReportSegments(const CompoundSystem *localCompound, WorkingCopy *localData) {
  calculatedConv_.clear();
  calculatedEnergies_.clear();
  calculatedAngles_.clear();
  calculatedAngularDists_.clear();
  calculatedSegments_.clear();
  calculatedSegmentsE1_.clear();
  calculatedSegmentsE2_.clear();
  calculatedExcitationEnergies_.clear();
  if (!localData) return 0;
  AZUREProfile::ScopedReport profileReport(configure().profile, "api segments");
  const std::vector<const Segment *> segments = ReportedSegments(localData->Data().GetSegments());
  const bool pre = PrecalculateAll(*localData, segments, localCompound, state(), configure(), results());
  for (const Segment *segment : segments) {
    const UnitRange<Unit> data = segment->GetPoints(*localData);
    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;
    // Angular-distribution coefficients, self-describing per point: a count followed by that many Legendre
    // coefficients (a zero for a point outside an angular-distribution segment).
    std::vector<double> angularDists;
    for (int k = 0; k < (int)data.size(); ++k) {
      cross.push_back(azure2::EvaluateModelPoint(*localData, *segment, k, localCompound, state(), configure(), results(),
                                                 pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate));
      crossE1.push_back(data[k].GetFitE1CrossSection(results()));
      crossE2.push_back(data[k].GetFitE2CrossSection(results()));
      angles.push_back(data[k].GetCMAngle());
      energies.push_back(data[k].GetCMEnergy());
      conv.push_back(data[k].GetSFactorConversion());
      excitationEnergies.push_back(data[k].GetExcitationEnergy());
      const std::vector<double> coefficients = segment->LegendreCoefficients(data[k], results());
      angularDists.push_back((double)coefficients.size());
      angularDists.insert(angularDists.end(), coefficients.begin(), coefficients.end());
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

int AZUREAPI::UpdateSegments(vector_r &p) {
  int k = 0;
  vector_r params_ = all_;
  for (int i = 0; i < all_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = p[k];
      ++k;
    }
  }

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  if (localCompound->GetIsospin().Active()) {
    // The transformation from observed widths mixes the levels of a J^pi group, so the observed values of derived
    // isospin members would leak into the free ones (IsospinTies.h): the fit-basis calls take such a project.
    configure().outStream << "ERROR: a project with isospin levels takes parameters in the fit's basis "
                          << "(UpdateSegmentsRWA / CalculateChi2RWA), not observed ones." << std::endl;
    ReportSegments(nullptr, nullptr);
    return -1;
  }
  localCompound->FillCompoundFromParamsPhysical(params_, state(), inputs());
  localData->FillAbundancesFromVector(params_);
  bool isValid = localCompound->TransformIn(configure(), state(), inputs());

  if (!isValid) {
    ReportSegments(nullptr, nullptr);  // nothing reported for an invalid vector (the arrays of the last call go)
    return 0;
  }

  localCompound->FillMnParams(params, state(), inputs());
  localData->FillMnParams(params);
  localCompound->FillCompoundFromParams(params.GetMinuitParams().Params(), state());
  // Norms are the same in both bases (model curves do not use them, but the
  // segments then hold the vector's values, as after CalculateChi2Physical).
  localData->FillNormsFromParams(params_);
  // After the compound: the sub-point grids are laid out from its fit values.
  localData->FillEnergyShiftsFromParams(params_, results(), localCompound, &state(), &inputs(), &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());

  return ReportSegments(localCompound, localData);
}

int AZUREAPI::UpdateSegmentsRWA(vector_r &p) {
  int k = 0;
  vector_r params_ = all_rwa_;
  for (int i = 0; i < all_rwa_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = p[k];
      ++k;
    }
  }

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_, state());
  localData->FillEnergyShiftsFromParams(params_, results(), localCompound, &state(), &inputs(), &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());

  return ReportSegments(localCompound, localData);
}

int AZUREAPI::UpdateSegmentsAllRWA(vector_r &p) {
  int k = 0;
  for (int i = 0; i < all_rwa_.size(); ++i) {
    all_rwa_[i] = p[i];
    if (!fixed_[i]) {
      values_[k] = p[i];
      ++k;
    }
  }
  vector_r params_ = all_rwa_;

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_, state());
  localData->FillEnergyShiftsFromParams(params_, results(), localCompound, &state(), &inputs(), &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());

  return ReportSegments(localCompound, localData);
}


int AZUREAPI::UpdateData() {
  dataEnergies_.clear();
  dataAngles_.clear();
  dataSegments_.clear();
  dataSegmentsErrors_.clear();
  dataConv_.clear();
  dataExcitationEnergies_.clear();

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  const std::vector<Segment> &segments = localData->Data().GetSegments();
  for (int i = 0; i < segments.size(); ++i) {
    newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    ++nSegments;

    const UnitRange<Unit> data = segments[i].GetPoints(*localData);

    std::vector<double> energies, angles, cross, crossErr, conv, excitationEnergies;
    for (int k = 0; k < data.size(); ++k) {
      energies.push_back(data[k].GetCMEnergy());
      angles.push_back(data[k].GetCMAngle());
      cross.push_back(segments[i].CMValue(localData->State(segments[i]), k));
      crossErr.push_back(segments[i].CMError(localData->State(segments[i]), k));
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

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey = -1;
  int prevKey = -1;
  int nSegments = 0;

  const std::vector<Segment> &segments = localData->Data().GetSegments();
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

bool AZUREAPI::WriteOutputFiles() {
  if (data_ == nullptr) return false;
  // chiSquared.out reports the Park, abundance and scattering-length totals from the WorkingCopy object; the CLI's
  // chi2 path stores them there, the API's CalculateChi2RWA does not (pyazr adds those terms in Python), so
  // through pyazr the file carried 0 or stale values (audit 2026-10-25, item 30).  Evaluate them at the current
  // parameters before writing.
  const CompoundSystem *localCompound = compound();
  if (localCompound) {
    const vector_r none;
    ObjectiveContext context{localCompound, state(), inputs(), data(), results(), none, configure()};
    Objective::Engine().ModelChi2(context);
  }
  data()->WriteOutputFiles(configure(), results(), false, nullptr);
  return true;
}

std::vector<ThmExperimentReport> AZUREAPI::GetThmExperiments() const {
  if (data_ == nullptr) return std::vector<ThmExperimentReport>();
  return data()->ThmExperimentReports(results());
}

bool AZUREAPI::GetThmLineshape(const std::string &name, const vector_r &energies, ThmLineshapeReport &out,
                               std::string &why) {
  if (data_ == nullptr || compound_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmLineshapeTable(name, energies, compound(), state(), configure(), out, why);
}

bool AZUREAPI::GetThmVertex(const std::string &name, const vector_r &energies, ThmVertexReport &out,
                            std::string &why) {
  if (data_ == nullptr || compound_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmVertexTable(name, energies, compound(), state(), configure(), out, why);
}

bool AZUREAPI::GetThmDistortion(const std::string &name, const vector_r &energies, ThmDistortionReport &out,
                                std::string &why) {
  if (data_ == nullptr) {
    why = "no data loaded";
    return false;
  }
  return data()->ThmDistortionTable(name, energies, out, why);
}

vector_r AZUREAPI::CalculateScatteringLengthsRWA(const vector_r &rwaParams) {
  vector_r params_ = all_rwa_;
  int k = 0;
  for (int i = 0; i < (int)all_rwa_.size(); ++i)
    if (!fixed_[i]) params_[i] = rwaParams[k++];
  const CompoundSystem *localCompound = compound();
  localCompound->FillCompoundFromParams(params_, state());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());
  vector_r out;
  for (int aa = 1; aa <= localCompound->NumPairs(); aa++) {
    const PPair *pair = localCompound->GetPair(aa);
    if (pair->GetPType() != 0) continue;
    for (const ScatteringLength::Channel &ch : ScatteringLength::Channels(localCompound, pair->GetPairKey())) {
      ScatteringLength::Value v = ScatteringLength::Evaluate(localCompound, state(), configure(), ch);
      out.insert(out.end(), {(double)ch.pairKey, ch.s, (double)ch.parity, std::real(v.a0), std::imag(v.a0),
                             ch.hasLevels ? 1.0 : 0.0});
    }
  }
  return out;
}


vector_r AZUREAPI::GetCurrentNorms() const {
  vector_r out;
  int prevKey = -1;
  const std::vector<Segment> &segments = data()->Data().GetSegments();
  for (size_t i = 0; i < segments.size(); ++i) {
    int newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;
    out.push_back(data()->NormUsed((int)i + 1, results()));
  }
  return out;
}

vector_r AZUREAPI::ParkNorms(const vector_r &rwaParams) {
  vector_r params_ = all_rwa_;
  for (int i = 0, k = 0; i < (int)all_rwa_.size(); ++i)
    if (!fixed_[i] && k < (int)rwaParams.size()) params_[i] = rwaParams[k++];
  const CompoundSystem *localCompound = compound();
  localCompound->FillCompoundFromParams(params_, state());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());
  vector_r out;
  for (int j = 1; j <= localCompound->NumJGroups(); j++) {
    const JGroup *jg = localCompound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++)
      if (jg->GetLevel(la)->IsInRMatrix()) out.push_back(state().ParkNorm(j, la));
  }
  return out;
}

vector_r AZUREAPI::ObservedLevelRowsAllRWA(const vector_r &p) {
  if (p.size() != all_rwa_.size())
    throw std::invalid_argument("observed levels: expected every parameter (" + std::to_string(all_rwa_.size()) +
                                " values, fixed ones included), got " + std::to_string(p.size()));
  const vector_r &params = p;
  const CompoundSystem *localCompound = compound();
  localCompound->FillCompoundFromParams(params, state());
  const ObservedLevelTable table = localCompound->TransformOut(configure(), state(), inputs());
  vector_r out;
  for (int j = 1; j <= localCompound->NumJGroups(); j++)
    for (int la = 1; la <= localCompound->GetJGroup(j)->NumLevels(); la++) {
      const ObservedLevel &level = table.Level(j, la);
      for (int ch = 1; ch <= localCompound->GetJGroup(j)->NumChannels(); ch++)
        out.insert(out.end(), {(double)j, (double)la, (double)ch, level.E, (double)level.kind[ch - 1], level.W[ch - 1],
                               level.value[ch - 1], level.P[ch - 1], level.mu_u[ch - 1], level.a_fm[ch - 1]});
    }
  return out;
}
