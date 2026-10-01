#include "AZUREOutput.h"
#include "CNuc.h"
#include "PPair.h"
#include "ChannelFunc.h"
#include "ThmFunc.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "EData.h"
#include "ECAmplitudeCache.h"
#include "ExtrapLine.h"
#include "SegLine.h"
#include "AdaptiveIntegrationGrid.h"
#include "Minuit2/MnUserParameters.h"
#include "GSLException.h"
#include "NuclearPotentialManager.h"
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <unordered_map>

/*!
 * The EData object has a private attribute containing the number of iterations needed to find the best fit parameters.
 * At creation, this attribute is set to 0.
 */

namespace {

/// The <thm> weight table named for segment line `key`, or null.
std::shared_ptr<const ThmWeightTable> ThmWeightFor(
    const std::map<int, std::shared_ptr<const ThmWeightTable>> &tables, int key) {
  auto it = tables.find(key);
  return it == tables.end() ? nullptr : it->second;
}

/*!
 * Checks the <thm> weight[k]= (or weightTest[k]=) keys against the segments
 * just read: k must be a line of the block, a segment in use must be THM, and
 * every one of its points must lie inside the table (folding sub-points may
 * reach beyond it; they get the end value).  A key naming an inactive or
 * unusable line is ignored with a warning.  Returns -1 after an ERROR line.
 */
int CheckThmWeights(const Config &configure,
                    const std::map<int, std::shared_ptr<const ThmWeightTable>> &tables,
                    const char *keyName, const char *blockName, int numLines, EData *data) {
  for (const auto &entry : tables) {
    const int key = entry.first;
    const ThmWeightTable &table = *entry.second;
    if (key > numLines) {
      configure.outStream << "ERROR: <thm> " << keyName << "[" << key << "]: " << blockName << " has only "
                          << numLines << " line(s)." << std::endl;
      return -1;
    }
    ESegment *segment = nullptr;
    for (int s = 1; s <= data->NumSegments(); s++)
      if (data->GetSegment(s)->GetSegmentKey() == key) segment = data->GetSegment(s);
    if (!segment) {
      configure.outStream << "WARNING: <thm> " << keyName << "[" << key << "]: segment line " << key
                          << " of " << blockName << " is not in use; the weight is ignored." << std::endl;
      continue;
    }
    if (!segment->IsTHM()) {
      configure.outStream << "ERROR: <thm> " << keyName << "[" << key << "]: segment " << key << " of "
                          << blockName << " is not a THM segment (isDiff < 10)." << std::endl;
      return -1;
    }
    for (int p = 1; p <= segment->NumPoints(); p++) {
      double energy = segment->GetPoint(p)->GetCMEnergy();
      if (!table.Covers(energy)) {
        configure.outStream << "ERROR: <thm> " << keyName << "[" << key << "]: segment " << key
                            << " has a point at E_cm = " << energy << " MeV, outside the table '"
                            << table.name << "' [" << table.e.front() << ", " << table.e.back() << "] MeV."
                            << std::endl;
        return -1;
      }
    }
  }
  return 0;
}

}  // namespace

EData::EData() {
  iterations_ = 0;
  normParamOffset_ = 0;
  energyShiftParamOffset_ = 0;
  isFit_ = true;
  isErrorAnalysis_ = false;
  ecReadPos_ = std::streampos(0);
  ecUsePrevious_ = false;
}

/*!
 * Returns the number of segment objects in the ESegment vector.
 */

int EData::NumSegments() const {
  return segments_.size();
}

/*!
 * This function fills the data object with segments from the segment data files.
 * After a segment is created, the ESegment::Fill method is called for that segment.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int EData::Fill(const Config &configure, CNuc *theCNuc) {
  std::ifstream in(configure.configfile.c_str());
  if (!in) return -1;
  std::string line = "";
  while (line != "<segmentsData>" && !in.eof()) getline(in, line);
  if (line != "<segmentsData>") return -1;
  line = "";
  int numTotalSegments = 0;
  while (!in.eof() && line != "</segmentsData>") {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (!in.eof() && line != "</segmentsData>") {
      std::istringstream stm;
      stm.str(line);
      SegLine segment(stm);
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      numTotalSegments++;
      if (segment.isActive() == 1) {
        ESegment NewSegment(segment);

        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()))->SetEntrance();
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            for (int i = 1; i <= theCNuc->NumPairs(); i++) {
              if (theCNuc->GetPair(i)->GetPType() == 10) {
                isValidTotal = true;
                break;
              }
            }
          }
          if (isValidTotal || theCNuc->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightBySegment, numTotalSegments));
            this->AddSegment(NewSegment);

            // Fill the segment with data first
            if (this->GetSegment(this->NumSegments())->Fill(theCNuc, this, configure) == -1) {
              configure.outStream << "WARNING: Could Not Fill Segment #" << this->NumSegments()
                                  << " from file." << std::endl;
              this->DeleteLastSegment();
            } else if (this->GetSegment(this->NumSegments())->NumPoints() == 0) {
              configure.outStream << "WARNING: Segment #" << numTotalSegments
                                  << " is empty and will not be used." << std::endl;
              this->DeleteLastSegment();
            } else {
              // Handle advanced segments - parse components and add them to the segment
              // This must be done AFTER the segment is filled with data points
              if (NewSegment.IsAdvanced()) {
                // Set the operation type from the segment data
                OperationType operation = (NewSegment.GetLegacyOperationType() == 0) ? SUM : RATIO;
                ESegment *filledSegment = this->GetSegment(this->NumSegments());
                filledSegment->SetOperationType(operation);

                // Parse components and create full segment copies
                std::string componentsStr = NewSegment.GetComponentsList();
                if (!componentsStr.empty()) {
                  // Components are stored as "Entrance: X, Exit: Y;Entrance: A, Exit: B;..." or with optional "Angle: Z" and "Scaling: S"
                  std::istringstream stream(componentsStr);
                  std::string component;
                  while (std::getline(stream, component, ';')) {
                    if (!component.empty()) {
                      // Parse "Entrance: X, Exit: Y" format (with optional ", Angle: Z" and ", Scaling: S")
                      size_t entrancePos = component.find("Entrance: ");
                      size_t exitPos = component.find("Exit: ");
                      size_t anglePos = component.find("Angle: ");
                      size_t scalingPos = component.find("Scaling: ");

                      double componentScaling = 1.0;
                      if (scalingPos != std::string::npos) {
                        try {
                          componentScaling = std::stod(component.substr(scalingPos + 9));
                        } catch (...) {
                          componentScaling = 1.0;
                        }
                      }

                      if (entrancePos != std::string::npos && exitPos != std::string::npos) {
                        entrancePos += 10;  // Length of "Entrance: "
                        size_t commaPos = component.find(", Exit: ");
                        if (commaPos != std::string::npos) {
                          int entranceKey = std::stoi(component.substr(entrancePos, commaPos - entrancePos));
                          exitPos += 6;  // Length of "Exit: "

                          // Check if there's an angle specification
                          ESegment *componentSegment = nullptr;
                          if (anglePos != std::string::npos) {
                            // Parse the exit key up to ", Angle:"
                            size_t angleCommaPos = component.find(", Angle: ");
                            int exitKey = std::stoi(component.substr(exitPos, angleCommaPos - exitPos));

                            // Parse the angle value (std::stod stops at the trailing ", Scaling: ...")
                            anglePos += 7;  // Length of "Angle: "
                            double fixedAngle = std::stod(component.substr(anglePos));

                            // Create component segment with fixed angle
                            componentSegment = this->CreateComponentSegment(*filledSegment, entranceKey, exitKey, fixedAngle);
                          } else {
                            // No angle specified, use standard method
                            int exitKey = std::stoi(component.substr(exitPos));
                            componentSegment = this->CreateComponentSegment(*filledSegment, entranceKey, exitKey);
                          }

                          if (componentSegment) {
                            componentSegment->SetComponentScaling(componentScaling);
                            filledSegment->AddComponentSegment(componentSegment);
                          }
                        }
                      }
                    }
                  }
                }
              }

              int thisSegmentNum = this->NumSegments();
              if (this->GetSegment(thisSegmentNum)->IsTotalCapture()) {
                ESegment *baseSegment = this->GetSegment(thisSegmentNum);
                int numCapturePairs = 0;

                // Set operation type to SUM for total capture
                baseSegment->SetOperationType(SUM);

                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(theCNuc->GetPair(i)->GetPairKey());
                    } else {
                      // For additional capture pairs, create component segments instead of separate segments
                      ESegment *componentSegment = this->CreateComponentSegment(*baseSegment,
                                                                                baseSegment->GetEntranceKey(),
                                                                                theCNuc->GetPair(i)->GetPairKey());
                      if (componentSegment) {
                        baseSegment->AddComponentSegment(componentSegment);
                      }
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
            }
          } else {
            if (NewSegment.GetExitKey() == -1) {
              configure.outStream << "WARNING: Total capture specified but no capture pair exists."
                                  << std::endl;
            } else {
              configure.outStream << "WARNING: Pair key " << NewSegment.GetExitKey()
                                  << " not in compound nucleus." << std::endl;
            }
          }
        } else
          configure.outStream << "WARNING: Pair key " << NewSegment.GetEntranceKey()
                              << " not in compound nucleus." << std::endl;
      }
    }
  }

  if (line != "</segmentsData>") return -1;

  in.close();

  if (CheckThmWeights(configure, configure.thm.weightBySegment, "weight", "<segmentsData>",
                      numTotalSegments, this) != 0)
    return -1;
  if (BuildThmGroups(configure, theCNuc, numTotalSegments) != 0) return -1;

  if (this->NumSegments() > 0) {
    if (this->ReadTargetEffectsFile(configure, theCNuc) == -1) return -1;
    this->MapData();
  }

  return 0;
}

/*!
 * If the AZURE calculation is not data driven, this function is called in place
 * of the EData::Fill function to create points at specified energies and angles.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int EData::MakePoints(const Config &configure, CNuc *theCNuc) {
  std::ifstream in(configure.configfile.c_str());
  if (!in) return -1;
  std::string line = "";
  while (line != "<segmentsTest>" && !in.eof()) getline(in, line);
  if (line != "<segmentsTest>") return -1;
  line = "";
  int numTotalSegments = 0;
  while (!in.eof() && line != "</segmentsTest>") {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (!in.eof() && line != "</segmentsTest>") {
      std::istringstream stm;
      stm.str(line);
      ExtrapLine segment(stm);
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      numTotalSegments++;
      if (segment.isActive() == 1) {
        ESegment NewSegment(segment);
        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            for (int i = 1; i <= theCNuc->NumPairs(); i++) {
              if (theCNuc->GetPair(i)->GetPType() == 10) {
                isValidTotal = true;
                break;
              }
            }
          }
          // See ESegment::Fill: isDiff 8 has no photon implementation.  A test
          // segment costs no chi2, but it would write a column of zeros to
          // AZUREOut and read as a prediction, so drop it with a reason too.
          const bool polProductCapture =
              NewSegment.IsPolarizationProduct() && theCNuc->IsPairKey(NewSegment.GetExitKey()) &&
              theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetExitKey()))->GetPType() == 10;
          if (polProductCapture) {
            configure.outStream
                << "WARNING: Test segment #" << numTotalSegments
                << " is Polarization x Cross Section (isDiff 8) with a capture exit channel,"
                << " which is not implemented; it will not be used." << std::endl;
          }
          if (!polProductCapture && (isValidTotal || theCNuc->IsPairKey(NewSegment.GetExitKey()))) {
            NewSegment.SetSegmentKey(numTotalSegments);
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightByTestSegment, numTotalSegments));
            this->AddSegment(NewSegment);
            ESegment *theSegment = this->GetSegment(this->NumSegments());

            theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()))->SetEntrance();
            PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()));
            PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetExitKey()));
            double aStep = theSegment->GetAStep();
            double eStep = theSegment->GetEStep();
            for (double angle = theSegment->GetMinAngle();
                 angle <= theSegment->GetMaxAngle(); angle += aStep) {
              for (double energy = theSegment->GetMinEnergy();
                   energy <= theSegment->GetMaxEnergy(); energy += eStep) {
                EPoint NewPoint(angle, energy, theSegment);
                theSegment->AddPoint(NewPoint);
                EPoint *thePoint = theSegment->GetPoint(theSegment->NumPoints());
                thePoint->SetParentData(this);
                if (entrancePair->GetPType() == 20)
                  thePoint->ConvertDecayEnergy(exitPair);
                else if (!theSegment->IsCMDifferential())
                  thePoint->ConvertLabEnergy(entrancePair);
                else if (theSegment->IsCMDifferential())
                  thePoint->ConvertLabEnergy(entrancePair);
                if (exitPair->GetPType() == 0 && theSegment->IsDifferential() &&
                    !theSegment->IsPhase() && !theSegment->IsAngularDist() && !theSegment->IsCMDifferential()) {
                  if (theSegment->GetEntranceKey() == theSegment->GetExitKey()) {
                    thePoint->ConvertLabAngle(entrancePair);
                  } else {
                    thePoint->ConvertLabAngle(entrancePair, exitPair, configure);
                  }
                  thePoint->ConvertCrossSection(entrancePair, exitPair);
                }
                if (exitPair->GetPType() == 10 && theSegment->IsDifferential() &&
                    !theSegment->IsPhase() && !theSegment->IsAngularDist() && !theSegment->IsCMDifferential()) {
                  thePoint->ConvertLabAngleGammas(entrancePair);
                  thePoint->ConvertCrossSectionGammas(entrancePair);
                }
                if (eStep == 0.0) break;
              }
              if (aStep == 0.0) break;
            }

            // Handle advanced segments (sum/ratio of components) - AFTER points are added
            if (NewSegment.IsAdvanced()) {
              int operation = NewSegment.GetOperationType();
              theSegment->SetOperationType((OperationType)operation);

              // Parse components and create full segment copies
              std::string componentsStr = NewSegment.GetComponentsList();
              if (!componentsStr.empty()) {
                // Components are stored as "Entrance: X, Exit: Y;Entrance: A, Exit: B;..." or with optional "Angle: Z" and "Scaling: S"
                std::istringstream stream(componentsStr);
                std::string component;
                int componentCount = 0;
                while (std::getline(stream, component, ';')) {
                  if (!component.empty()) {
                    // Parse "Entrance: X, Exit: Y" format (with optional ", Angle: Z" and ", Scaling: S")
                    size_t entrancePos = component.find("Entrance: ");
                    size_t exitPos = component.find("Exit: ");
                    size_t anglePos = component.find("Angle: ");
                    size_t scalingPos = component.find("Scaling: ");

                    double componentScaling = 1.0;
                    if (scalingPos != std::string::npos) {
                      try {
                        componentScaling = std::stod(component.substr(scalingPos + 9));
                      } catch (...) {
                        componentScaling = 1.0;
                      }
                    }

                    if (entrancePos != std::string::npos && exitPos != std::string::npos) {
                      entrancePos += 10;  // Length of "Entrance: "
                      size_t commaPos = component.find(", Exit: ");
                      if (commaPos != std::string::npos) {
                        int entranceKey = std::stoi(component.substr(entrancePos, commaPos - entrancePos));
                        exitPos += 6;  // Length of "Exit: "

                        // Check if there's an angle specification
                        ESegment *componentSegment = nullptr;
                        if (anglePos != std::string::npos) {
                          // Parse the exit key up to ", Angle:"
                          size_t angleCommaPos = component.find(", Angle: ");
                          int exitKey = std::stoi(component.substr(exitPos, angleCommaPos - exitPos));

                          // Parse the angle value (std::stod stops at the trailing ", Scaling: ...")
                          anglePos += 7;  // Length of "Angle: "
                          double fixedAngle = std::stod(component.substr(anglePos));

                          // Create component segment with fixed angle
                          componentSegment = this->CreateComponentSegment(*theSegment, entranceKey, exitKey, fixedAngle);
                        } else {
                          // No angle specified, use standard method
                          int exitKey = std::stoi(component.substr(exitPos));
                          componentSegment = this->CreateComponentSegment(*theSegment, entranceKey, exitKey);
                        }

                        if (componentSegment) {
                          componentSegment->SetComponentScaling(componentScaling);
                          theSegment->AddComponentSegment(componentSegment);
                          componentCount++;
                        }
                      }
                    }
                  }
                }
              } else {
              }
            }
            if (theSegment->NumPoints() == 0) {
              configure.outStream << "WARNING: Extrapolation segment #" << numTotalSegments
                                  << " is empty and will not be used." << std::endl;
              this->DeleteLastSegment();
            } else {
              int thisSegmentNum = this->NumSegments();
              if (this->GetSegment(thisSegmentNum)->IsTotalCapture()) {
                ESegment *baseSegment = this->GetSegment(thisSegmentNum);
                int numCapturePairs = 0;

                // Set operation type to SUM for total capture
                baseSegment->SetOperationType(SUM);

                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(theCNuc->GetPair(i)->GetPairKey());
                    } else {
                      // For additional capture pairs, create component segments instead of separate segments
                      ESegment *componentSegment = this->CreateComponentSegment(*baseSegment,
                                                                                baseSegment->GetEntranceKey(),
                                                                                theCNuc->GetPair(i)->GetPairKey());
                      if (componentSegment) {
                        baseSegment->AddComponentSegment(componentSegment);
                      }
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
            }
          } else {
            if (NewSegment.GetExitKey() == -1) {
              configure.outStream << "WARNING: Total capture specified but no capture pair exists."
                                  << std::endl;
            } else {
              configure.outStream << "WARNING: Pair key " << NewSegment.GetExitKey()
                                  << " not in compound nucleus." << std::endl;
            }
          }
        } else
          configure.outStream << "WARNING: Pair key " << NewSegment.GetEntranceKey()
                              << " not in compound nucleus." << std::endl;
      }
    }
  }

  if (line != "</segmentsTest>") return -1;

  in.close();

  if (CheckThmWeights(configure, configure.thm.weightByTestSegment, "weightTest", "<segmentsTest>",
                      numTotalSegments, this) != 0)
    return -1;

  if (this->NumSegments() > 0) {
    if (this->ReadTargetEffectsFile(configure, theCNuc) == -1) return -1;
    this->MapData();
  }

  return 0;
}

/*!
 * Returns the number of fit iterations needed to minimize the parameters to the data.
 */

int EData::Iterations() const {
  return iterations_;
}

/*!
 * Returns the number of TargetEffect objects contained in the present object.
 */

int EData::NumTargetEffects() const {
  return targetEffects_.size();
}

/*!
 * Returns the offset of the normalization paramters in the Minuit parameter vector.
 */
int EData::GetNormParamOffset() const {
  return normParamOffset_;
}

/*!
 * Returns the offset of the energy shift parameters in the Minuit parameter vector.
 */
int EData::GetEnergyShiftParamOffset() const {
  return energyShiftParamOffset_;
}

namespace {

/*
Builds the sub-points of a target-effect point on [endEnergy, startEnergy] (CM)
and records on the point how the grid was built, so EPoint::RefreshSubPointGrid
can rebuild it the same way from the current fit parameters.
*/
void FillSubPoints(EPoint *point, ESegment *segment, TargetEffect *targetEffect, double startEnergy,
                   double endEnergy, double cmConversion, int entranceKey, bool refreshable,
                   CNuc *compound, const Config &configure) {
  std::vector<double> energyGrid;
  int numPoints = targetEffect->NumSubPoints();
  SubGridSpec spec;
  if (configure.useAdaptiveGrid) {
    AdaptiveIntegrationGrid::GridConfig gridConfig;
    gridConfig.maxPoints = numPoints;
    gridConfig.entranceKey = entranceKey;
    gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
    gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
    gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
    gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
    AdaptiveIntegrationGrid gridGenerator(gridConfig);
    std::vector<AdaptiveIntegrationGrid::ResonanceInfo> anchors = gridGenerator.Anchors(compound);
    energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, anchors);
    spec.refreshable = refreshable;
    spec.startEnergy = startEnergy;
    spec.endEnergy = endEnergy;
    spec.baseEnergyStep = gridConfig.baseEnergyStep;
    spec.resonanceWidthMultiplier = gridConfig.resonanceWidthMultiplier;
    spec.pointsPerWidth = gridConfig.pointsPerWidth;
    spec.entranceKey = entranceKey;
    spec.cmConversion = cmConversion;
    for (const AdaptiveIntegrationGrid::ResonanceInfo &r : gridGenerator.AnchorsInReach(startEnergy, endEnergy, anchors)) {
      spec.anchors.push_back(r.id);
      spec.anchors.push_back(r.energy);
      spec.anchors.push_back(r.particleWidth);
    }
  } else {
    double step = (startEnergy - endEnergy) / numPoints;
    for (int i = 0; i <= numPoints; i++)
      energyGrid.push_back(startEnergy - i * step);
  }

  for (size_t i = 0; i < energyGrid.size(); i++) {
    double subEnergy = energyGrid[i];
    EPoint subPoint(point->GetCMAngle(), subEnergy, segment);
    if (targetEffect->IsTargetIntegration()) {
      double stoppingPower = cmConversion * targetEffect->GetStoppingPowerEq()->Evaluate(configure, subEnergy / cmConversion);
      subPoint.SetStoppingPower(stoppingPower);
    }
    point->AddSubPoint(subPoint);
  }
  point->SetSubGridSpec(spec);
}

}  // namespace

/*!
 * Reads the target effects input file and creates the TargetEffect objects
 * to be applied to the data.
 */

int EData::ReadTargetEffectsFile(const Config &configure, CNuc *compound) {
  std::ifstream in(configure.configfile.c_str());
  if (!in) return -1;
  std::string line = "";
  while (line != "<targetInt>" && !in.eof()) getline(in, line);
  if (line != "<targetInt>") return -1;
  line = "";
  while (line != "</targetInt>" && !in.eof()) {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (line != "</targetInt>" && !in.eof()) {
      std::istringstream stm;
      stm.str(line);
      TargetEffect targetEffect(stm, configure);
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      if (targetEffect.IsActive()) {
        // Every listed segment gets its own copy of the effect.  The effect's
        // energies (sigma, beam profile) are converted lab -> c.m. below with
        // the segment's own factor; one object shared by "1,2" was converted
        // once per listed segment, and both segments were then folded with
        // the twice-converted sigma (17 keV c.m. became 16.1 keV for 18O+p).
        // A copy per segment also lets segments of different entrance pairs
        // share one line.
        std::vector<int> segmentsList = targetEffect.GetSegmentsList();
        for (int i = 1; i <= segmentsList.size(); i++) {
          if (this->IsSegmentKey(segmentsList[i - 1])) {
            ESegment *segment = this->GetSegmentFromKey(segmentsList[i - 1]);
            if (segment) {
              this->AddTargetEffect(targetEffect);
              segment->SetTargetEffectNum(this->NumTargetEffects());
              // If segment has components, set the target effect to these as
              // well -- again a copy each, since each is converted below.
              if (segment->HasComponents()) {
                for (auto component : segment->GetComponentSegments()) {
                  this->AddTargetEffect(targetEffect);
                  component->SetTargetEffectNum(this->NumTargetEffects());
                }
              }
            }
          }
        }
      }
    }
  }
  if (line != "</targetInt>") return -1;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(segment->GetEntranceKey()));
    PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(segment->GetExitKey()));
    double cmConversion;
    if (entrancePair->GetPType() == 20)
      cmConversion = (exitPair->GetM(1) + exitPair->GetM(2)) / exitPair->GetM(2);
    else
      cmConversion = entrancePair->GetM(2) / (entrancePair->GetM(1) + entrancePair->GetM(2));

    if (segment->IsTargetEffect()) {
      TargetEffect *targetEffect = this->GetTargetEffect(segment->GetTargetEffectNum());
      double sigma = targetEffect->GetSigma();
      targetEffect->SetSigma(cmConversion * sigma);
      if (targetEffect->IsBeamProfile()) targetEffect->ConvertBeamProfileToCM(cmConversion);

      for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
        // An effect restricted to energy ranges leaves points outside them
        // completely untouched: they never receive the effect number, so the
        // whole downstream machinery treats them as ordinary points.
        double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy());
        if (blendWeight <= 0.) continue;
        point->SetTargetEffectNum(segment->GetTargetEffectNum());
        point->SetTargetBlendWeight(blendWeight);

        if (targetEffect->IsSubPointEffect()) {
          double forwardDepth = 0.0;
          double backwardDepth = 0.0;

          if (targetEffect->IsTargetIntegration()) {
            double totalM = entrancePair->GetM(1) + entrancePair->GetM(2);
            double targetThickness = cmConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure);
            point->SetTargetThickness(targetThickness);
            if (targetEffect->IsConvolution() || targetEffect->IsConvCoefficients()) {
              if (targetEffect->IsConvCoefficients()) {
                backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure) * 5.0;
                forwardDepth = targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure) * 5.0;
              } else {
                backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->GetSigma() * 5.0;
                forwardDepth = targetEffect->convolutionRange * targetEffect->GetSigma() * 5.0;
              }
              // Add straggling range extension at the deepest layer
              if (targetEffect->IsStraggling()) {
                double targetThickness_keV = targetThickness * 1000.0;  // Convert MeV to keV
                double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                double stragglingSigma = stragglingSigma_keV / 1000.0;  // Convert back to MeV
                double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                backwardDepth += stragglingRange;
              }
            } else {
              backwardDepth = targetThickness;
              forwardDepth = 0.0;
              // Add straggling range extension for target integration only
              if (targetEffect->IsStraggling()) {
                double targetThickness_keV = targetThickness * 1000.0;
                double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                double stragglingSigma = stragglingSigma_keV / 1000.0;
                double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                backwardDepth += stragglingRange;
              }
            }
          }

          else if (targetEffect->IsConvolution() || targetEffect->IsConvCoefficients()) {
            double convRange = segment->IsTHM() ? targetEffect->thmConvolutionRange
                                                : targetEffect->convolutionRange;
            if (targetEffect->IsConvCoefficients()) {
              backwardDepth = convRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure);
              forwardDepth = convRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure);
            } else {
              backwardDepth = convRange * targetEffect->GetSigma();
              forwardDepth = convRange * targetEffect->GetSigma();
            }
          }

          // Generate adaptive integration grid
          double startEnergy = point->GetCMEnergy() + forwardDepth;
          double endEnergy;

          /*
          THM (HOES) segments carry data below the entrance threshold, so their
          sub-point grids may extend to negative energies as long as the
          compound-system energy Ecm+SepE+ExE stays positive.
          */
          double minSubEnergy = segment->IsTHM() ? TargetEffect::minIntegrationEnergy - (entrancePair->GetSepE() + entrancePair->GetExE()) : TargetEffect::minIntegrationEnergy;

          // Check if target thickness is zero (for convolution-only cases)
          if (targetEffect->IsTargetIntegration()) {
            double targetThickness = point->GetTargetThickness();
            if (targetThickness < 1.0e-10) {
              // Target thickness is effectively zero - integrate down to the floor energy
              endEnergy = minSubEnergy;
              // Set density to 1e24 to prevent division issues
              targetEffect->SetDensity(1.0e24);
            } else {
              // Normal case - use backward depth
              endEnergy = point->GetCMEnergy() - backwardDepth;
              // Safety check: if backwardDepth > energy, set endEnergy to minimum
              if (endEnergy < minSubEnergy) {
                endEnergy = minSubEnergy;
              }
            }
          } else {
            // Convolution or ConvCoefficients - use backward depth
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if (endEnergy < minSubEnergy) {
              endEnergy = minSubEnergy;
            }
          }

          if (targetEffect->IsBeamProfile()) {
            // The kernel is an absolute beam profile: sample where the beam
            // is, narrowed to the point's energy window (plus the resolution
            // tails) when it has one, not around the point's own energy.
            double low, high;
            targetEffect->BeamProfileSupport(low, high);
            double s = targetEffect->GetBeamTpcSigma();
            // The resolution tails: 4 s, or for a THM (HOES) point the 5 s of
            // its Gaussian fold -- a narrow peak 4-5 s from a point dominates
            // the tail there (-2 % at 81 keV from a 289 eV level, 18O(p,a) THM).
            double tails = segment->IsTHM() ? TargetEffect::thmConvolutionRange : 4.0;
            if (point->HasBinWindow()) {
              low = std::max(low, point->GetBinLowCM() - tails * s);
              high = std::min(high, point->GetBinHighCM() + tails * s);
            }
            // The floor of the Gaussian fold: a THM point may sit below the entrance
            // threshold (this was +1 keV, and the fold of a point at -45 keV came
            // out 66 times too large).
            if (low < minSubEnergy) low = minSubEnergy;
            if (high <= low) high = low + 0.001;
            startEnergy = high;
            endEnergy = low;
            point->SetPhotoKinematics(entrancePair->GetSepE() - exitPair->GetExE(),
                                      (entrancePair->GetM(1) + entrancePair->GetM(2)) * uconv);
          }
          // The grid of a main segment may be rebuilt later when a fit moves a
          // narrow level (EPoint::RefreshSubPointGrid); a beam-profile window
          // is absolute and is left as built.
          FillSubPoints(&*point, &*segment, targetEffect, startEnergy, endEnergy, cmConversion,
                        segment->GetEntranceKey(), !targetEffect->IsBeamProfile(), compound, configure);
        }
      }
    }
    if (segment->HasComponents()) {
      for (auto component : segment->GetComponentSegments()) {
        if (component->IsTargetEffect()) {
          TargetEffect *targetEffect = this->GetTargetEffect(component->GetTargetEffectNum());
          double sigma = targetEffect->GetSigma();
          targetEffect->SetSigma(cmConversion * sigma);
          if (targetEffect->IsBeamProfile()) targetEffect->ConvertBeamProfileToCM(cmConversion);

          for (EPointIterator point = component->GetPoints().begin(); point < component->GetPoints().end(); point++) {
            double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy());
            if (blendWeight <= 0.) continue;
            point->SetTargetEffectNum(component->GetTargetEffectNum());
            point->SetTargetBlendWeight(blendWeight);

            if (targetEffect->IsSubPointEffect()) {
              double forwardDepth = 0.0;
              double backwardDepth = 0.0;

              if (targetEffect->IsTargetIntegration()) {
                double totalM = entrancePair->GetM(1) + entrancePair->GetM(2);
                double targetThickness = cmConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure);
                point->SetTargetThickness(targetThickness);
                if (targetEffect->IsConvolution() || targetEffect->IsConvCoefficients()) {
                  if (targetEffect->IsConvCoefficients()) {
                    backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure) * 5.0;
                    forwardDepth = targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure) * 5.0;
                  } else {
                    backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->GetSigma() * 5.0;
                    forwardDepth = targetEffect->convolutionRange * targetEffect->GetSigma() * 5.0;
                  }
                  // Add straggling range extension at the deepest layer (component segments)
                  if (targetEffect->IsStraggling()) {
                    double targetThickness_keV = targetThickness * 1000.0;  // Convert MeV to keV
                    double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                    double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                    double stragglingSigma = stragglingSigma_keV / 1000.0;  // Convert back to MeV
                    double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                    backwardDepth += stragglingRange;
                  }
                } else {
                  backwardDepth = targetThickness;
                  forwardDepth = 0.0;
                  // Add straggling range extension for target integration only (component segments)
                  if (targetEffect->IsStraggling()) {
                    double targetThickness_keV = targetThickness * 1000.0;
                    double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                    double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                    double stragglingSigma = stragglingSigma_keV / 1000.0;
                    double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                    backwardDepth += stragglingRange;
                  }
                }
              }

              else if (targetEffect->IsConvolution() || targetEffect->IsConvCoefficients()) {
                double convRange = component->IsTHM() ? targetEffect->thmConvolutionRange
                                                      : targetEffect->convolutionRange;
                if (targetEffect->IsConvCoefficients()) {
                  backwardDepth = convRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure);
                  forwardDepth = convRange * targetEffect->CalculateSigma(point->GetLabEnergy(), configure);
                } else {
                  backwardDepth = convRange * targetEffect->GetSigma();
                  forwardDepth = convRange * targetEffect->GetSigma();
                }
              }

              // Generate adaptive integration grid for component segment
              double startEnergy = point->GetCMEnergy() + forwardDepth;
              double endEnergy;

              /*
              THM (HOES) segments carry data below the entrance threshold, so their
              sub-point grids may extend to negative energies as long as the
              compound-system energy Ecm+SepE+ExE stays positive. Note this uses the
              component's own entrance pair, since an advanced (sum/ratio) component
              can have a different entrance channel than its master segment.
              */
              PPair *componentEntrancePair = compound->GetPair(compound->GetPairNumFromKey(component->GetEntranceKey()));
              double minSubEnergy = component->IsTHM() ? TargetEffect::minIntegrationEnergy - (componentEntrancePair->GetSepE() + componentEntrancePair->GetExE()) : TargetEffect::minIntegrationEnergy;

              // Check if target thickness is zero (for convolution-only cases)
              if (targetEffect->IsTargetIntegration()) {
                double targetThickness = point->GetTargetThickness();
                if (targetThickness < 1.0e-10) {
                  // Target thickness is effectively zero - integrate down to the floor energy
                  endEnergy = minSubEnergy;
                  // Set density to 1e24 to prevent division issues
                  targetEffect->SetDensity(1.0e24);
                } else {
                  // Normal case - use backward depth
                  endEnergy = point->GetCMEnergy() - backwardDepth;
                  // Safety check: if backwardDepth > energy, set endEnergy to minimum
                  if (endEnergy < minSubEnergy) {
                    endEnergy = minSubEnergy;
                  }
                }
              } else {
                // Convolution or ConvCoefficients - use backward depth
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if (endEnergy < minSubEnergy) {
                  endEnergy = minSubEnergy;
                }
              }

              if (targetEffect->IsBeamProfile()) {
                double low, high;
                targetEffect->BeamProfileSupport(low, high);
                double s = targetEffect->GetBeamTpcSigma();
                double tails = component->IsTHM() ? TargetEffect::thmConvolutionRange : 4.0;
                if (point->HasBinWindow()) {
                  low = std::max(low, point->GetBinLowCM() - tails * s);
                  high = std::min(high, point->GetBinHighCM() + tails * s);
                }
                if (low < minSubEnergy) low = minSubEnergy;  // as above
                if (high <= low) high = low + 0.001;
                startEnergy = high;
                endEnergy = low;
                point->SetPhotoKinematics(entrancePair->GetSepE() - exitPair->GetExE(),
                                          (entrancePair->GetM(1) + entrancePair->GetM(2)) * uconv);
              }
              // A component's points are converted again after this
              // (InitializeComponentSegments), so its grid is never rebuilt.
              FillSubPoints(&*point, &*component, targetEffect, startEnergy, endEnergy, cmConversion,
                            segment->GetEntranceKey(), false, compound, configure);
            }
          }
        }
      }
    }
  }
  return 0;
}

/*!
 * Returns true if the data is to be fit, otherwise returns false.  Used in the AZURECalc function class
 * to determine if a clone of the CNuc and EData objects should be made for thread safety.
 */

bool EData::IsFit() const {
  return isFit_;
}

/*!
 * Returns true if the call to function is for error analysis via Minos, otherwise returns false.  Used in the AZURECalc function
 * class to suppress transformation and file output during error analysis.
 */

bool EData::IsErrorAnalysis() const {
  return isErrorAnalysis_;
}

/*!
 * Sets the boolean indicating if the data is to be fit by AZURECalc function class. Used in AZUREMain function
 * class before calls to Minuit and AZURECalc.
 */

/*!
 * Returns true if the specified segment key exists corresponds to a segment in the ESegment vector,
 * otherwise returns false.
 */

bool EData::IsSegmentKey(int segmentKey) {
  bool isKey = false;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->GetSegmentKey() == segmentKey) {
      isKey = true;
      break;
    }
  }
  return isKey;
}

/*!
 * Sets an internal variable specifying if the data is to be fit by Minuit.  Needed to determine cloning behavior in AZURECalc for
 * thread safety.
 */

void EData::SetFit(bool fit) {
  isFit_ = fit;
}

/*!
 * Sets the boolean indicating if the call to the function is for error analysis via Minos.
 */

void EData::SetErrorAnalysis(bool errorAnalysis) {
  isErrorAnalysis_ = errorAnalysis;
}

/*!
 * This function updates the number of fit iterations per iteration during the fitting process.
 */

void EData::Iterate() {
  iterations_++;
}

/*!
 * This function sets the number of iterations to zero.
 */

void EData::ResetIterations() {
  iterations_ = 0;
}

/*!
 * This function is identical in role to the EPoint::Initialize function, except that it initializes
 * and entire EData object instead of a single EPoint object.
 */

int EData::Initialize(CNuc *compound, const Config &configure) {
  // Identical-particle phase-shift segment validation: warn once per segment
  // if the requested L is forbidden by Bose/Fermi symmetry on an identical
  // pair. The phase branch in GenMatrixFunc would silently report delta = 0
  // for these (no matching channel), which is easy to misread as "no nuclear
  // interaction" rather than "this partial wave does not exist".
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (!seg->IsPhase()) continue;
    int aa = compound->GetPairNumFromKey(seg->GetEntranceKey());
    if (aa == 0) continue;
    PPair *pp = compound->GetPair(aa);
    if (!pp->IsIdentical()) continue;
    int segL = seg->GetL();
    bool allowedFound = false;
    for (int j = 1; j <= compound->NumJGroups() && !allowedFound; j++) {
      JGroup *jg = compound->GetJGroup(j);
      if (jg->GetJ() != seg->GetJ()) continue;
      for (int ch = 1; ch <= jg->NumChannels() && !allowedFound; ch++) {
        AChannel *channel = jg->GetChannel(ch);
        if (channel->GetRadType() != 'P') continue;
        if (channel->GetPairNum() != aa) continue;
        if (channel->GetL() == segL) allowedFound = true;
      }
    }
    if (!allowedFound) {
      configure.outStream << "**WARNING: Phase-shift segment #" << s
                          << " requests J=" << seg->GetJ()
                          << ", L=" << segL
                          << " for identical pair (Z=" << pp->GetZ(1)
                          << ", A=" << pp->GetM(1)
                          << "); no matching partial wave exists "
                          << "(likely forbidden by Bose/Fermi symmetry, "
                          << "which requires L+S even). The fit will report delta = 0 at every energy."
                          << std::endl;
    }
  }

  // Calculate channel lo-matrix and channel penetrability for each channel at each local energy
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Lo-Matrix, Phases, and Penetrabilities..." << std::endl;
  if (this->CalcEDependentValues(compound, configure) == -1) return -1;
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_ENERGY_DEP)
    this->PrintEDependentValues(configure, compound);
  // Calculate legendre polynomials for each data point
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Legendre Polynomials..." << std::endl;
  this->CalcLegendreP(configure.maxLOrder, compound);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_LEGENDRE)
    this->PrintLegendreP(configure);

  // Calculate Coulomb Amplitudes
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Coulomb Amplitudes..." << std::endl;
  this->CalcCoulombAmplitude(compound);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_COUL_AMPLITUDES) {
    this->PrintCoulombAmplitude(configure, compound);
  }

  // Calculate new ec amplitudes
  if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
    if (!(configure.paramMask & Config::USE_API))
      configure.outStream << "Calculating External Capture Amplitudes..." << std::endl;
    if (this->CalculateECAmplitudes(compound, configure) == -1) return -1;
  }

  // Initialize component segments - ensure their entrance/exit pairs are properly set up
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Initializing Component Segments..." << std::endl;
  if (this->InitializeComponentSegments(compound, configure) == -1) return -1;

  return 0;
}

/*!
 * Adds a segment to the ESegment vector.
 */

void EData::AddSegment(ESegment segment) {
  segments_.push_back(segment);
}

/*!
 * Prints the data point after the object is filled or points are created.
 */

void EData::PrintData(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_DATA) {
    std::string outfile = configure.checkdir + "data.chk";
    fbuffer.open(outfile.c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_DATA)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_DATA) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_DATA)) {
    out << std::endl
        << "************************************" << std::endl
        << "*            Segments              *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
        << std::setw(17) << "Segment Key #"
        << std::setw(17) << "Entrance Key #"
        << std::setw(13) << "Exit Key #"
        << std::setw(12) << "Min Energy"
        << std::setw(12) << "Max Energy"
        << std::setw(11) << "Min Angle"
        << std::setw(11) << "Max Angle"
        << std::setw(25) << "Data File"
        << std::endl;
    for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
      out << std::setw(11) << segment - GetSegments().begin() + 1
          << std::setw(17) << segment->GetSegmentKey()
          << std::setw(17) << segment->GetEntranceKey()
          << std::setw(13) << segment->GetExitKey()
          << std::setw(12) << segment->GetMinEnergy()
          << std::setw(12) << segment->GetMaxEnergy()
          << std::setw(11) << segment->GetMinAngle()
          << std::setw(11) << segment->GetMaxAngle()
          << std::setw(25) << segment->GetDataFile()
          << std::endl;
    }
    out << std::endl
        << "************************************" << std::endl
        << "*               Data               *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
        << std::setw(14) << "Data Point #"
        << std::setw(15) << "Lab Energy"
        << std::setw(15) << "CM Energy"
        << std::setw(15) << "Angle"
        << std::setw(20) << "Cross Section"
        << std::setw(22) << "Cross Section Error"
        << std::setw(12) << "Map Point"
        << std::setw(18) << "# of Subpoints"
        << std::setw(18) << "Low Sub Energy"
        << std::setw(18) << "High Sub Energy"
        << std::endl;
    for (EDataIterator data = begin(); data != end(); data++) {
      out << std::setw(11) << data.segment() - GetSegments().begin() + 1
          << std::setw(14) << data.point() - (data.segment()->GetPoints()).begin() + 1
          << std::setw(15) << data.point()->GetLabEnergy()
          << std::setw(15) << data.point()->GetCMEnergy()
          << std::setw(15) << data.point()->GetCMAngle()
          << std::setw(20) << data.point()->GetCMCrossSection()
          << std::setw(22) << data.point()->GetCMCrossSectionError();
      if (data.point()->IsMapped()) {
        EnergyMap map = data.point()->GetMap();
        char tempMap[25];
        snprintf(tempMap, sizeof(tempMap), "(%d,%d)", map.segment, map.point);
        out << std::setw(12) << tempMap << std::endl;
      } else
        out << std::setw(12) << "Not Mapped"
            << std::setw(18) << data.point()->NumSubPoints();
      if (data.point()->IsTargetEffect() &&
          (data.point()->GetParentData()->GetTargetEffect(data.point()->GetTargetEffectNum())->IsConvolution() ||
           data.point()->GetParentData()->GetTargetEffect(data.point()->GetTargetEffectNum())->IsTargetIntegration() ||
           data.point()->GetParentData()->GetTargetEffect(data.point()->GetTargetEffectNum())->IsBeamProfile())) {
        out << std::setw(18) << data.point()->GetSubPoint(data.point()->NumSubPoints())->GetCMEnergy()
            << std::setw(18) << data.point()->GetSubPoint(1)->GetCMEnergy();
      }
      out << std::endl;
      if (data.point() == data.segment()->GetPoints().end() - 1) out << std::endl;
    }
  } else
    configure.outStream << "Could not write data check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Calls EPoint::CalcLegendreP for each point in the entire EData object.
 */

void EData::CalcLegendreP(int maxL, CNuc *theCNuc) {
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    TargetEffect *effect = (segment->IsTargetEffect() &&
                            this->GetTargetEffect(segment->GetTargetEffectNum())->IsQCoefficients())
        ? this->GetTargetEffect(segment->GetTargetEffectNum())
        : NULL;
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) {
      EPoint *point = segment->GetPoint(i);
      point->CalcLegendreP(maxL, theCNuc, effect);
    }
  }
}

/*!
 * Prints the Legendre polynomials for each point in the EData object.
 */

void EData::PrintLegendreP(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_LEGENDRE) {
    std::string outfile = configure.checkdir + "legendre.chk";
    fbuffer.open(outfile.c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_LEGENDRE)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_LEGENDRE) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_LEGENDRE)) {
    out << std::endl
        << "************************************" << std::endl
        << "*       Legendre Polynomials       *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "Segment #"
        << std::setw(10) << "Point #"
        << std::setw(15) << "CM Energy"
        << std::setw(15) << "Angle"
        << std::setw(5) << "L"
        << std::setw(15) << "Leg. Poly." << std::endl;
    for (EDataIterator data = begin(); data != end(); data++) {
      for (int lOrder = 0; lOrder <= data.point()->GetMaxLOrder(); lOrder++) {
        out << std::setw(10) << data.segment() - GetSegments().begin() + 1
            << std::setw(10) << data.point() - (data.segment()->GetPoints()).begin() + 1
            << std::setw(15) << data.point()->GetCMEnergy()
            << std::setw(15) << data.point()->GetCMAngle()
            << std::setw(5) << lOrder
            << std::setw(15) << data.point()->GetLegendreP(lOrder) << std::endl;
      }
    }
  } else
    configure.outStream << "Could not write legendre polynomials check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Calls EPoint::CalcEDependentValues for each point in the entire EData object.
 */

int EData::CalcEDependentValues(CNuc *theCNuc, const Config &configure) {
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    bool localStop = false;
#pragma omp parallel for shared(localStop, configure)
    for (int i = 1; i <= segment->NumPoints(); i++) {
      if (configure.stopFlag || localStop) continue;
      EPoint *point = segment->GetPoint(i);
      if (!(point->IsMapped())) {
        try {
          point->CalcEDependentValues(theCNuc, configure);
        } catch (GSLException e) {
#pragma omp critical
          {
            std::cout << point->GetLabEnergy() << "\t" << point->GetLabAngle() << std::endl;
            configure.outStream << e.what() << std::endl;
            localStop = true;
          }
        }
      }
    }
    if (configure.stopFlag || localStop) return -1;
  }
  return 0;
}

/*!
 * Prints the values calculated by EPoint::CalcEDependentValues for each point in the entire EData object.
 */

void EData::PrintEDependentValues(const Config &configure, CNuc *theCNuc) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_ENERGY_DEP) {
    std::string outfile = configure.checkdir + "lomatrixandpene.chk";
    fbuffer.open(outfile.c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_ENERGY_DEP)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_ENERGY_DEP) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_ENERGY_DEP)) {
    out << std::endl
        << "************************************" << std::endl
        << "*  Lo Matrix and Penetrabilities   *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "Seg_#"
        << std::setw(10) << "Point_#"
        << std::setw(5) << "j"
        << std::setw(5) << "ch"
        << std::setw(5) << "l"
        << std::setw(15) << "E_chan"
        << std::setw(15) << "sqrt_pene"
        << std::setw(25) << "Lo"
        << std::setw(25) << "expHSP" << std::endl;
    for (EDataIterator data = begin(); data != end(); data++) {
      double inEnergy = data.point()->GetCMEnergy() + theCNuc->GetPair(theCNuc->GetPairNumFromKey(data.segment()->GetEntranceKey()))->GetSepE();
      for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
        if (theCNuc->GetJGroup(j)->IsInRMatrix()) {
          JGroup *theJGroup = theCNuc->GetJGroup(j);
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
            AChannel *theChannel = theJGroup->GetChannel(ch);
            PPair *thePair = theCNuc->GetPair(theChannel->GetPairNum());
            int lValue = theChannel->GetL();
            double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
            complex loElement = data.point()->GetLoElement(j, ch);
            complex expHSP = data.point()->GetExpHardSpherePhase(j, ch);
            out << std::setw(10) << data.segment() - GetSegments().begin() + 1
                << std::setw(10) << data.point() - (data.segment()->GetPoints()).begin() + 1
                << std::setw(5) << j
                << std::setw(5) << ch
                << std::setw(5) << lValue
                << std::setw(15) << localEnergy
                << std::setw(15) << data.point()->GetSqrtPenetrability(j, ch)
                << std::setw(25) << loElement.real() << " " << loElement.imag()
                << std::setw(25) << expHSP.real() << " " << expHSP.imag() << std::endl;
          }
        }
      }
    }
  } else
    configure.outStream << "Could not write lo-matrix and penetrabilities check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Calls EPoint::CalcCoulombAmplitude for each point in the entire EData object.
 */

void EData::CalcCoulombAmplitude(CNuc *theCNuc) {
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) {
      EPoint *point = segment->GetPoint(i);
      point->CalcCoulombAmplitude(theCNuc);
    }
  }
}

/*!
 * Prints the values calculated by EPoint::CalcCoulombAmplitude for each point in the entire EData object.
 */

void EData::PrintCoulombAmplitude(const Config &configure, CNuc *theCNuc) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES) {
    std::string outfile = configure.checkdir + "coulombamplitudes.chk";
    fbuffer.open(outfile.c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES)) {
    out << std::endl
        << "************************************" << std::endl
        << "*        Coulomb Amplitudes        *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "segment #"
        << std::setw(10) << "point #"
        << std::setw(10) << "aa"
        << std::setw(15) << "cmenergy"
        << std::setw(15) << "angle"
        << std::setw(25) << "coulomb amplitude"
        << std::endl;
    for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
      if (segment->GetEntranceKey() == segment->GetExitKey()) {
        for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
          out << std::setw(10) << segment - GetSegments().begin() + 1
              << std::setw(10) << point - segment->GetPoints().begin() + 1
              << std::setw(10) << theCNuc->GetPairNumFromKey(segment->GetEntranceKey())
              << std::setw(15) << point->GetCMEnergy()
              << std::setw(15) << point->GetCMAngle()
              << std::setw(25) << point->GetCoulombAmplitude()
              << std::endl;
        }
      }
    }
  } else
    configure.outStream << "Could not write coulomb amplitudes check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Writes the output files for the calculation.  The output files are all in center of mass frame, and contain
 * columns for energy, angle, calculated cross section, calculated s-factor, experimental cross section and error
 * and experimental s-factor and error.
 */

// One ".band" line: energy, excitation, angle, xs, d(xs), S-factor, d(S-factor).
static void WriteBandLine(std::ostream &o, double energy, double excitation,
                          double angle, double xs, double dxs, double conv) {
  o << std::setw(18) << std::scientific << energy
    << std::setw(18) << std::scientific << excitation
    << std::setw(18) << std::scientific << angle
    << std::setw(18) << std::scientific << xs
    << std::setw(18) << std::scientific << dxs
    << std::setw(18) << std::scientific << xs * conv
    << std::setw(18) << std::scientific << dxs * conv
    << std::endl;
}

void EData::WriteOutputFiles(const Config &configure, bool isFit, const BandData *band) {
  AZUREOutput output(configure.outputdir);
  std::ofstream chiOut;
  if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA)) {
    std::string chiOutFile = configure.outputdir + "chiSquared.out";
    chiOut.open(chiOutFile.c_str());
  }
  chiOut << "Segment#," << " Chi-Squared, " << " N, " << " Norm, " << " Norm-Chi-Squared " << std::endl;

  if (!(configure.paramMask & Config::CALCULATE_WITH_DATA)) output.SetExtrap();
  bool isVaryNorm = false;
  double totalChiSquared = 0.;
  double totalNormChiSquared = 0.;
  double totalN = 0.;
  // use kinflag.dat to control output option
  std::ifstream kinoption;
  kinoption.open("kinflag.dat");
  int kinflag = 0;  // 0 cm output, 1 lab output
  if (kinoption) {
    kinoption >> kinflag;
    kinoption.close();
  }
  if (kinflag != 0) configure.outStream << "Using alternate output format..." << std::endl;

  // THM experiments: refresh their profile from the models the points hold
  // (a snapshot written during a fit has not been profiled), so that the norm
  // and the background written are those of these models.
  if (configure.paramMask & Config::CALCULATE_WITH_DATA)
    for (int g = 0; g < NumThmGroups(); g++)
      if (!thmGroups_[g].trivial) ProfileThmGroup(g);

  // When a covariance is available, write a sibling ".band" file per output file,
  // with the same block/point structure so the GUI can pair them.
  bool writeBand = band && (configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) && !band->grad.empty();
  std::map<std::string, std::ofstream> bandFiles;
  auto bandStreamFor = [&](int aa, int ir) -> std::ofstream * {
    std::string name = configure.outputdir + "AZUREOut_aa=" + std::to_string(aa);
    if (ir == -1)
      name += "_TOTAL_CAPTURE";
    else
      name += "_R=" + std::to_string(ir);
    name += output.IsExtrap() ? ".extrap.band" : ".out.band";
    std::ofstream &f = bandFiles[name];
    if (!f.is_open()) {
      f.open(name.c_str());
      f.precision(10);
    }
    return f.is_open() ? &f : nullptr;
  };
  // Look up a point's parameter-sensitivity row, or nullptr if absent.
  auto gradFor = [&](EPoint *p) -> const std::vector<double> * {
    if (!band) return nullptr;
    std::map<EPoint *, std::vector<double>>::const_iterator it = band->grad.find(p);
    return (it == band->grad.end()) ? nullptr : &it->second;
  };

  ESegmentIterator firstSumIterator = GetSegments().end();
  for (ESegmentIterator segment = GetSegments().begin();
       segment < GetSegments().end(); segment++) {
    if (segment->IsTotalCapture()) {
      firstSumIterator = segment;
      // With component segments, no need to skip - total capture is handled within the segment
    }
    if (segment->IsVaryNorm()) isVaryNorm = true;
    int aa = segment->GetEntranceKey();
    int ir = segment->GetExitKey();
    std::filebuf *buf;
    if (firstSumIterator != GetSegments().end())
      buf = output(aa, -1);
    else {
      if (segment->IsAngularDist() &&
          !(configure.paramMask & Config::CALCULATE_WITH_DATA))
        buf = output(aa, ir, true);
      else
        buf = output(aa, ir);
    }
    std::ostream out(buf);
    // A THM experiment with a background: the fitted curve is the model plus b(E).
    const int segmentIndex = (int)(segment - GetSegments().begin()) + 1;
    const bool thmBackground = !output.IsExtrap() && ThmGroupOf(segmentIndex) >= 0 &&
                               thmGroups_[ThmGroupOf(segmentIndex)].terms > 0;
    ESegmentIterator thisSegment = segment;
    if (firstSumIterator != GetSegments().end()) thisSegment = firstSumIterator;

    // Band stream for this segment (skip angular-coefficient extrap files, which
    // have no scalar cross section; other segments get a block for alignment).
    bool bandThisSegment = writeBand && !(segment->IsAngularDist() && output.IsExtrap());
    std::ofstream *bandOut = bandThisSegment ? bandStreamFor(aa, (firstSumIterator != GetSegments().end()) ? -1 : ir) : nullptr;
    // Point gradient, summed across total-capture components like the value is.
    auto pointBandGrad = [&](EPointIterator point) -> std::vector<double> {
      std::vector<double> g;
      const std::vector<double> *g0 = gradFor(&*point);
      if (g0) g = *g0;
      if (firstSumIterator != GetSegments().end()) {
        int pointIndex = point - segment->GetPoints().begin() + 1;
        for (ESegmentIterator it = firstSumIterator; it < segment; it++) {
          const std::vector<double> *gi = gradFor(it->GetPoint(pointIndex));
          if (!gi) continue;
          if (g.empty())
            g = *gi;
          else
            for (size_t k = 0; k < g.size() && k < gi->size(); k++) g[k] += (*gi)[k];
        }
      }
      return g;
    };

    if (kinflag == 0) {
      for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
        out.precision(10);
        if (segment->IsAngularDist()) {
          out << std::setw(18) << std::scientific << point->GetCMEnergy();
          for (int i = 0; i < point->GetNumAngularDists(); i++) out << std::setw(18) << point->GetAngularDist(i);
          out << std::endl;
          if (bandOut) WriteBandLine(*bandOut, point->GetCMEnergy(), point->GetExcitationEnergy(),
                                     point->GetCMAngle(), 0., 0., 0.);
        } else {
          double fitCrossSection = point->GetFitCrossSection();
          if (thmBackground) fitCrossSection += ThmBackgroundAt(segmentIndex, point->GetCMEnergy());
          if (firstSumIterator != GetSegments().end()) {
            int pointIndex = point - segment->GetPoints().begin() + 1;
            for (ESegmentIterator it = firstSumIterator; it < segment; it++)
              fitCrossSection += it->GetPoint(pointIndex)->GetFitCrossSection();
          }
          out << std::setw(18) << std::scientific << point->GetCMEnergy()
              << std::setw(18) << std::scientific << point->GetExcitationEnergy()
              << std::setw(18) << std::scientific << point->GetCMAngle()
              << std::setw(18) << std::scientific << fitCrossSection
              << std::setw(18) << std::scientific << fitCrossSection * point->GetSFactorConversion();
          if (!output.IsExtrap()) {
            double dataNorm = thisSegment->GetNorm();
            out << std::setw(18) << std::scientific << point->GetCMCrossSection() * dataNorm
                << std::setw(18) << std::scientific << point->GetCMCrossSectionError() * dataNorm
                << std::setw(18) << std::scientific << point->GetCMCrossSection() * dataNorm * point->GetSFactorConversion()
                << std::setw(18) << std::scientific << point->GetCMCrossSectionError() * dataNorm * point->GetSFactorConversion()
                << std::endl;
          } else
            out << std::endl;
          if (bandOut) {
            std::vector<double> g = pointBandGrad(point);
            double dxs = g.empty() ? 0. : band->dXS(g);
            WriteBandLine(*bandOut, point->GetCMEnergy(), point->GetExcitationEnergy(),
                          point->GetCMAngle(), fitCrossSection, dxs, point->GetSFactorConversion());
          }
        }
      }
    }
    if (kinflag == 1) {
      for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
        out.precision(10);
        if (segment->IsAngularDist()) {
          out << std::setw(18) << std::scientific << point->GetCMEnergy();
          for (int i = 0; i < point->GetNumAngularDists(); i++) out << std::setw(18) << point->GetAngularDist(i);
          out << std::endl;
          if (bandOut) WriteBandLine(*bandOut, point->GetLabEnergy(), point->GetExcitationEnergy(),
                                     point->GetLabAngle(), 0., 0., 0.);
        } else {
          double fitCrossSection = (point->GetFitCrossSection() +
                                    (thmBackground ? ThmBackgroundAt(segmentIndex, point->GetCMEnergy()) : 0.0)) /
                                   point->GetCrossSectionKinFactor();
          if (firstSumIterator != GetSegments().end()) {
            int pointIndex = point - segment->GetPoints().begin() + 1;
            for (ESegmentIterator it = firstSumIterator; it < segment; it++)
              fitCrossSection += it->GetPoint(pointIndex)->GetFitCrossSection();
          }
          out << std::setw(18) << std::scientific << point->GetLabEnergy()
              << std::setw(18) << std::scientific << point->GetExcitationEnergy()
              << std::setw(18) << std::scientific << point->GetLabAngle()
              << std::setw(18) << std::scientific << fitCrossSection
              << std::setw(18) << std::scientific << fitCrossSection * point->GetSFactorConversion();
          if (!output.IsExtrap()) {
            double dataNorm = thisSegment->GetNorm();
            out << std::setw(18) << std::scientific << point->GetLabCrossSection() * dataNorm
                << std::setw(18) << std::scientific << point->GetLabCrossSectionError() * dataNorm
                << std::setw(18) << std::scientific << point->GetLabCrossSection() * dataNorm * point->GetSFactorConversion()
                << std::setw(18) << std::scientific << point->GetLabCrossSectionError() * dataNorm * point->GetSFactorConversion()
                << std::endl;
          } else
            out << std::endl;
          if (bandOut) {
            std::vector<double> g = pointBandGrad(point);
            double kin = point->GetCrossSectionKinFactor();
            double dxs = g.empty() ? 0. : band->dXS(g) / kin;
            WriteBandLine(*bandOut, point->GetLabEnergy(), point->GetExcitationEnergy(),
                          point->GetLabAngle(), fitCrossSection, dxs, point->GetSFactorConversion());
          }
        }
      }
    }
    if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA)) {
      // The normalization penalty as the fit sees it (AZURECalc::operator()):
      // the deviation from the nominal normalization in units of its uncertainty,
      // which the segment stores as a percentage.  This used to add a literal 1
      // per segment, so the reported total was just the number of segments.
      double normChiSquared = 0.;
      double normError = thisSegment->GetNominalNorm() / 100. * thisSegment->GetNormError();
      if (normError != 0.) normChiSquared =
                               pow((thisSegment->GetNorm() - thisSegment->GetNominalNorm()) / normError, 2.0);
      totalChiSquared += (thisSegment->GetSegmentChiSquared());
      totalNormChiSquared += normChiSquared;
      totalN += thisSegment->NumPoints();
      chiOut << thisSegment->GetSegmentKey()
             << ","
             << thisSegment->GetSegmentChiSquared()
             << ","
             << thisSegment->NumPoints()
             << ","
             << thisSegment->GetNorm()
             << ","
             << normChiSquared
             << std::endl;
    }
    out << std::endl
        << std::endl;
    out.flush();
    if (bandOut) {
      *bandOut << std::endl
               << std::endl;
      bandOut->flush();
    }
    firstSumIterator = GetSegments().end();
  }


  if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA)) {
    chiOut << "Total-Chi-Squared: "
           << totalChiSquared
           << " Total-Norm-Chi-Squared: "
           << totalNormChiSquared
           << " Total-N: "
           << totalN
           << std::endl
           << std::endl;
    chiOut.flush();
    chiOut.close();
  }
  if (!isFit && (configure.paramMask & Config::CALCULATE_WITH_DATA) && NumThmGroups() > 0)
    WriteThmExperiments(configure);
  if (isVaryNorm) {
    std::string outputfile = configure.outputdir + "normalizations.out";
    std::ofstream out(outputfile.c_str());
    if (out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
        if (segment->IsVaryNorm()) out << segment->GetSegmentKey() << ","
                                       << segment->GetDataFile() << ","
                                       << segment->GetMinAngle() << ","
                                       << segment->GetMaxAngle() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetMinEnergy() << ","
            << segment->GetMaxEnergy() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetNorm() << ","
            << segment->GetEnergyShift() << std::endl;
      }
      out.flush();
      out.close();
    } else
      configure.outStream << "Could not write normalization file." << std::endl;
  }

  // Check if any energy shifts are varied
  bool isVaryEnergyShift = false;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryEnergyShift()) {
      isVaryEnergyShift = true;
      break;
    }
  }

  // Write energy shifts file
  if (isVaryEnergyShift) {
    std::string outputfile = configure.outputdir + "shifts.out";
    std::ofstream out(outputfile.c_str());
    if (out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
        if (segment->IsVaryEnergyShift()) out << segment->GetSegmentKey() << ","
                                              << segment->GetDataFile() << ","
                                              << segment->GetMinAngle() << ","
                                              << segment->GetMaxAngle() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetMinEnergy() << ","
            << segment->GetMaxEnergy() << ",";
        if (segment->IsCMDifferential() == true) {
          out << "CM,";
        } else
          out << "Lab,";
        out << segment->GetNorm() << ","
            << segment->GetEnergyShift() << std::endl;
      }
      out.flush();
      out.close();
    } else
      configure.outStream << "Could not write energy shifts file." << std::endl;
  }
}

namespace {

/*!
 * 64-bit FNV-1a over a stream of values, each written as text.  Doubles are
 * printed at fixed precision (%.12e) rather than hashed bit for bit, so the
 * signature survives a value that went through a text file and back, while
 * still moving for any change that matters to the integrals.
 */
class ECSignatureHash {
 public:
  void Add(const char *text) {
    for (const char *c = text; *c; c++) {
      hash_ ^= (unsigned char)*c;
      hash_ *= 1099511628211ULL;
    }
    // A separator, so that ("1","23") and ("12","3") hash differently.
    hash_ ^= (unsigned char)';';
    hash_ *= 1099511628211ULL;
  }
  void Add(double value) {
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "%.12e", value);
    Add(buffer);
  }
  void Add(int value) {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%d", value);
    Add(buffer);
  }
  std::string Hex() const {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%016" PRIx64, hash_);
    return buffer;
  }

 private:
  uint64_t hash_ = 14695981039346656037ULL;
};

/*!
 * Walks the external-capture amplitudes of one set of segments in exactly the
 * order CalculateECAmplitudes and InitializeComponentSegments write them, and
 * returns how many there are.  With a hash, it also feeds in what each block of
 * amplitudes depends on, and collects the pairs involved so their properties
 * can be added afterwards.
 */
template <class Segments>
long long WalkECAmplitudes(Segments &segments, CNuc *theCNuc, ECSignatureHash *hash,
                           std::set<int> *pairs) {
  long long count = 0;
  for (auto segment = segments.begin(); segment != segments.end(); segment++) {
    int aa = theCNuc->GetPairNumFromKey(segment->GetEntranceKey());
    if (theCNuc->GetPair(aa)->GetPType() == 20) continue;
    if (!theCNuc->GetPair(aa)->IsEntrance()) continue;
    PPair *entrancePair = theCNuc->GetPair(aa);
    for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
      for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
        if (!theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) continue;
        ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
        int ir = theCNuc->GetPairNumFromKey(segment->GetExitKey());
        if (ecLevel->GetECPairNum() != ir) continue;
        Decay *decay = entrancePair->GetDecay(ir);
        if (hash) {
          // Which block this is, and the final state it captures to.  The
          // level energy enters every integral (EPoint::CalculateECAmplitudes);
          // the ANC does not -- it multiplies the amplitude later -- so it is
          // deliberately left out and can be varied without invalidating.
          hash->Add("segment");
          hash->Add(segment->GetSegmentKey());
          hash->Add(segment->GetEntranceKey());
          hash->Add(segment->GetExitKey());
          hash->Add(theCNuc->GetJGroup(j)->GetJ());
          hash->Add(ecLevel->GetE());
          pairs->insert(aa);
          for (int k = 1; k <= decay->NumKGroups(); k++) {
            KGroup *theKGroup = decay->GetKGroup(k);
            hash->Add(theKGroup->GetS());
            for (int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
              ECMGroup *theECMGroup = theKGroup->GetECMGroup(ecm);
              AChannel *finalChannel = theCNuc->GetJGroup(j)->GetChannel(theECMGroup->GetFinalChannel());
              hash->Add(theECMGroup->GetL());
              hash->Add(theECMGroup->GetMult());
              hash->Add((int)theECMGroup->GetRadType());
              hash->Add(theECMGroup->GetJ());
              hash->Add(finalChannel->GetL());
              hash->Add(finalChannel->GetS());
              hash->Add(theCNuc->GetPair(finalChannel->GetPairNum())->GetPairKey());
              pairs->insert(finalChannel->GetPairNum());
              hash->Add(theECMGroup->IsChannelCapture() ? 1 : 0);
              if (theECMGroup->IsChannelCapture()) {
                MGroup *chanCap = entrancePair->GetDecay(theECMGroup->GetChanCapDecay())->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
                AChannel *initial = theCNuc->GetJGroup(chanCap->GetJNum())->GetChannel(chanCap->GetChpNum());
                hash->Add(initial->GetL());
                hash->Add(initial->GetS());
              }
            }
          }
        }
        for (auto point = segment->GetPoints().begin(); point != segment->GetPoints().end(); point++) {
          if (point->IsMapped()) continue;
          for (int k = 1; k <= decay->NumKGroups(); k++)
            for (int ecm = 1; ecm <= decay->GetKGroup(k)->NumECMGroups(); ecm++)
              count += 1 + (long long)point->NumSubPoints();
          if (hash) {
            // The energies the amplitudes are evaluated at.  The angle does
            // not enter an EC amplitude, so it is not part of the signature.
            // The lab energy is included as well as the c.m. one because a
            // component segment's c.m. energy is converted again, with its own
            // entrance pair, after this point (InitializeComponentSegments).
            hash->Add("point");
            hash->Add(point->GetLabEnergy());
            hash->Add(point->GetCMEnergy());
            hash->Add(point->NumSubPoints());
            for (auto subPoint = point->GetSubPoints().begin(); subPoint != point->GetSubPoints().end(); subPoint++)
              hash->Add(subPoint->GetCMEnergy());
          }
        }
      }
    }
  }
  return count;
}

}  // namespace

/*!
 * Counts the external-capture amplitudes this model expects to read from, or
 * write to, an intEC file: those of the data segments, which
 * CalculateECAmplitudes writes, followed by those of the component segments,
 * which InitializeComponentSegments appends.  The walk mirrors both loops
 * exactly; if either changes, WalkECAmplitudes must change with it.
 *
 * A count only catches a file built for a different number of points.  The
 * signature (ECSignature) is what catches the same number of points at
 * different energies.
 */

long long EData::CountECAmplitudes(CNuc *theCNuc, const Config &configure) {
  return WalkECAmplitudes(GetSegments(), theCNuc, nullptr, nullptr) +
      WalkECAmplitudes(componentSegments_, theCNuc, nullptr, nullptr);
}

/*!
 * The signature of the external-capture integrals this calculation needs.
 *
 * An intEC file records amplitudes at the sub-point energies of one grid and
 * carries no record of which grid.  Energy straggling, a target thickness or
 * the adaptive-grid settings move those energies without changing how many
 * there are; a channel radius, the Coulomb-function routine or the hybrid
 * potential change the integrand at the same energies.  Either way the count
 * check passes and the stale amplitudes are reused, silently.  The signature
 * covers all of it: every energy an amplitude is evaluated at, the structure
 * of each capture block, the level energy of each final state, the entrance
 * and final pairs (masses, charges, spins, separation and excitation
 * energies, channel radius, hybrid potential), and the options that select
 * how the Coulomb functions and the integrals are computed.
 *
 * Energy shifts: the amplitudes are computed once, here, at the energies the
 * points have at initialization, and are not recomputed when a fit moves a
 * segment's energy shift (ESegment::UpdatePointEnergiesWithShift refreshes the
 * penetrabilities and phases, not the EC amplitudes).  The file is therefore
 * the initial state only, and the signature is taken over the energies at the
 * time the file is written or read, which is the same stage of the same
 * calculation in both cases.
 */

std::string EData::ECSignature(CNuc *theCNuc, const Config &configure) {
  ECSignatureHash hash;
  // Bump the version if what goes into the signature changes, so files written
  // under the old definition are recomputed rather than wrongly matched.
  hash.Add("AZURE2 EC signature v1");
  hash.Add((configure.paramMask & Config::USE_GSL_COULOMB_FUNC) ? 1 : 0);
  hash.Add((configure.paramMask & Config::CALCULATE_REACTION_RATE) ? 1 : 0);
  hash.Add(configure.useHybridMethod ? 1 : 0);

  std::set<int> pairs;
  hash.Add("data segments");
  WalkECAmplitudes(GetSegments(), theCNuc, &hash, &pairs);
  hash.Add("component segments");
  WalkECAmplitudes(componentSegments_, theCNuc, &hash, &pairs);

  for (std::set<int>::const_iterator it = pairs.begin(); it != pairs.end(); ++it) {
    PPair *pair = theCNuc->GetPair(*it);
    hash.Add("pair");
    hash.Add(pair->GetPairKey());
    hash.Add(pair->GetPType());
    for (int i = 1; i <= 2; i++) {
      hash.Add(pair->GetZ(i));
      hash.Add(pair->GetM(i));
      hash.Add(pair->GetJ(i));
      hash.Add(pair->GetPi(i));
      hash.Add(pair->GetG(i));
    }
    hash.Add(pair->GetSepE());
    hash.Add(pair->GetExE());
    hash.Add(pair->GetChRad());
    // The hybrid potential bends the Coulomb functions of this pair, and so
    // both the entrance-channel phases and the final-state integrand.
    // (CoulFunc applies it when the global switch and the pair's own
    // setting are both on.)
    const NuclearPotentialManager &potentials = NuclearPotentialManager::instance();
    bool hybrid = configure.useHybridMethod && potentials.isPairEnabled(pair->GetPairKey());
    hash.Add(hybrid ? 1 : 0);
    if (hybrid) {
      NuclearPotentialSetting setting = potentials.getSetting(pair->GetPairKey());
      hash.Add(setting.type.c_str());
      hash.Add(setting.V0);
      hash.Add(setting.R);
      hash.Add(setting.a);
      hash.Add(setting.r0);
    }
  }
  return hash.Hex();
}

std::string EData::ECSignaturePath(const std::string &integralsFile) {
  return integralsFile + ".sig";
}

/*!
 * If external capture amplitudes are to be calculated, EPoint::CalculateECAmplitudes is
 * called for each point with a corresponding external capture component in the EData object.
 * Otherwise, the amplitudes are read from the specified file.
 */

int EData::CalculateECAmplitudes(CNuc *theCNuc, const Config &configure) {
  std::ifstream in;
  std::ofstream out;
  std::string outputfile;
  if (configure.paramMask & Config::CALCULATE_WITH_DATA)
    outputfile = configure.outputdir + "intEC.dat";
  else
    outputfile = configure.outputdir + "intEC.extrap";

  // An intEC file records amplitudes at the sub-point energies of the grid it
  // was built for.  Reusing one built for a different grid, or for different
  // inputs at the same energies, returns amplitudes that belong to another
  // calculation -- silently, and with a result that looks like physics.  Two
  // checks guard the reuse: the number of amplitudes in the file, and the
  // signature recorded beside it (ECSignature).  If either disagrees, say so
  // and recompute rather than proceed.
  ecSignature_ = ECSignature(theCNuc, configure);
  ecOutputFile_ = outputfile;
  bool usePrevious = (configure.paramMask & Config::USE_PREVIOUS_INTEGRALS) != 0;
  if (usePrevious) {
    long long expected = CountECAmplitudes(theCNuc, configure);
    std::ifstream check(configure.integralsfile.c_str());
    if (!check) {
      configure.outStream << "Could not open external capture file '"
                          << configure.integralsfile
                          << "'; the integrals will be recalculated." << std::endl;
      usePrevious = false;
    } else {
      long long found = 0;
      complex dummy(0.0, 0.0);
      while (check >> dummy) ++found;
      check.close();
      if (found != expected) {
        configure.outStream << "WARNING: '" << configure.integralsfile
                            << "' holds " << found << " external capture amplitudes but this"
                            << " calculation needs " << expected << "." << std::endl
                            << "         The file belongs to a different set of data segments or"
                            << " integration points." << std::endl
                            << "         Recalculating the integrals." << std::endl;
        usePrevious = false;
      }
    }
  }
  if (usePrevious) {
    std::string sigFile = ECSignaturePath(configure.integralsfile);
    std::ifstream sigIn(sigFile.c_str());
    std::string recorded, line;
    while (sigIn && std::getline(sigIn, line)) {
      std::istringstream fields(line);
      std::string tag, version, value;
      if (fields >> tag >> version >> value && tag == "AZURE2-EC-SIGNATURE") recorded = version + " " + value;
    }
    if (recorded.empty()) {
      // A file from before signatures existed (or one copied without its
      // sidecar): there is no record of what it was computed for, and the
      // cases that matter are exactly the ones the count cannot see.
      configure.outStream << "WARNING: '" << configure.integralsfile
                          << "' has no signature file ('" << sigFile << "')," << std::endl
                          << "         so there is no record of the energies and inputs its"
                          << " integrals were computed for." << std::endl
                          << "         Recalculating the integrals." << std::endl;
      usePrevious = false;
    } else if (recorded != "1 " + ecSignature_) {
      configure.outStream << "WARNING: '" << configure.integralsfile
                          << "' was computed for a different integration grid or different inputs" << std::endl
                          << "         (signature " << recorded.substr(recorded.find(' ') + 1)
                          << ", this calculation needs " << ecSignature_ << ")." << std::endl
                          << "         The number of amplitudes matches, but their energies or integrands do not:"
                          << std::endl
                          << "         straggling, a target thickness, the integration-grid settings, a channel"
                          << std::endl
                          << "         radius, a level energy, the Coulomb-function routine or the hybrid potential"
                          << std::endl
                          << "         may have changed.  Recalculating the integrals." << std::endl;
      usePrevious = false;
    }
  }
  ecUsePrevious_ = usePrevious;

  if (usePrevious)
    in.open(configure.integralsfile.c_str());
  else {
    // The signature is written only once the file is complete, after the
    // component segments are appended (InitializeComponentSegments).  Remove
    // the old one first, so that a run interrupted part way leaves a file
    // with no signature -- which is recomputed -- rather than one whose old
    // signature still vouches for it.
    std::remove(ECSignaturePath(outputfile).c_str());
    out.open(outputfile.c_str());
    if (!out) configure.outStream << "Could not write to EC Amplitude File." << std::endl;
  }
  int sumSegmentI = 0;
  int numSumSegments = 0;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsTotalCapture()) {
      numSumSegments = segment->IsTotalCapture();
      sumSegmentI = 1;  // Start at 1 since total capture is now handled within the segment
    } else {
      sumSegmentI = 0;
      numSumSegments = 0;
    }
    char segmentKeyOut[256];
    if (segment->IsTotalCapture())
      snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d (Total: %d)", segment->GetSegmentKey(), numSumSegments);
    else
      snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d", segment->GetSegmentKey());
    int aa = theCNuc->GetPairNumFromKey(segment->GetEntranceKey());
    if (theCNuc->GetPair(aa)->GetPType() == 20) continue;
    if (theCNuc->GetPair(aa)->IsEntrance()) {
      PPair *entrancePair = theCNuc->GetPair(aa);
      for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
        for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
          if (theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
            ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
            int ir = theCNuc->GetPairNumFromKey(segment->GetExitKey());
            if (ecLevel->GetECPairNum() == ir) {
              if (!usePrevious) {
                if (!(configure.paramMask & Config::USE_API)) {
                  configure.outStream << "\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [                         ] 0%";
                  configure.outStream.flush();
                }
                int numPoints = segment->NumPoints();
                int pointIndex = 0;
                time_t startTime = time(NULL);
                bool localStop = false;
#pragma omp parallel for shared(configure, localStop)
                for (int i = 1; i <= numPoints; i++) {
                  if (configure.stopFlag || localStop) continue;
                  EPoint *point = segment->GetPoint(i);
                  if (!(point->IsMapped())) {
                    try {
                      point->CalculateECAmplitudes(theCNuc, configure);
                    } catch (GSLException e) {
#pragma omp critical
                      {
                        configure.outStream << e.what() << std::endl;
                        localStop = true;
                      }
                    }
                  }
                  ++pointIndex;
                  if (difftime(time(NULL), startTime) > 0.25) {
                    startTime = time(NULL);
                    std::string progress = " [";
                    double percent = 0.;
                    for (int j = 1; j <= 25; j++) {
                      if (pointIndex >= percent * numPoints && percent < 1.) {
                        percent += 0.04;
                        progress += '*';
                      } else
                        progress += ' ';
                    }
                    progress += "] ";
                    if (!(configure.paramMask & Config::USE_API))
                      configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                          << std::setw(0) << progress << percent * 100 << '%';
                    configure.outStream.flush();
                  }
                }
                if (configure.stopFlag || localStop) {
                  if (out.is_open()) out.close();
                  if (in.is_open()) in.close();
                  return -1;
                }
                if (!(configure.paramMask & Config::USE_API))
                  configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [*************************] 100%" << std::endl;
              }
              for (EPointIterator point = segment->GetPoints().begin();
                   point < segment->GetPoints().end(); point++) {
                if (!(point->IsMapped())) {
                  for (int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
                    for (int ecm = 1; ecm <= entrancePair->GetDecay(ir)->GetKGroup(k)->NumECMGroups(); ecm++) {
                      if (!usePrevious) {
                        if (out.is_open()) out << point->GetECAmplitude(k, ecm) << std::endl;
                        for (EPointIterator subPoint = point->GetSubPoints().begin();
                             subPoint < point->GetSubPoints().end(); subPoint++)
                          if (out.is_open()) out << subPoint->GetECAmplitude(k, ecm) << std::endl;
                      } else {
                        complex ecAmplitude(0.0, 0.0);
                        in >> ecAmplitude;
                        point->AddECAmplitude(k, ecm, ecAmplitude);
                        for (EPointIterator subPoint = point->GetSubPoints().begin();
                             subPoint < point->GetSubPoints().end(); subPoint++) {
                          ecAmplitude = complex(0.0, 0.0);
                          in >> ecAmplitude;
                          subPoint->AddECAmplitude(k, ecm, ecAmplitude);
                        }
                      }
                      for (EPointMapIterator mappedPoint = point->GetMappedPoints().begin();
                           mappedPoint < point->GetMappedPoints().begin(); mappedPoint++) {
                        (*mappedPoint)->AddECAmplitude(k, ecm, point->GetECAmplitude(k, ecm));
                        for (int i = 1; i <= point->NumSubPoints(); i++) {
                          (*mappedPoint)->GetSubPoint(i)->AddECAmplitude(k, ecm, point->GetSubPoint(i)->GetECAmplitude(k, ecm));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (sumSegmentI == numSumSegments) {
      sumSegmentI = 0;
      numSumSegments = 0;
    }
  }
  if (out.is_open()) {
    out.flush();
    out.close();
  }
  // Remember where the regular-segment integrals end so the component-segment
  // integrals (appended to the same file) can be resumed from here on read.
  if (in.is_open()) {
    ecReadPos_ = in.tellg();
    in.close();
  }
  return 0;
}

/*!
 * Initialize component segments by ensuring their entrance/exit pairs are properly set up
 * in the compound nucleus for calculations.
 */

int EData::InitializeComponentSegments(CNuc *theCNuc, const Config &configure) {
  // Initialize all component segments with the same process as regular segments
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Initializing Component Segments (" << componentSegments_.size() << " components)..." << std::endl;

  // First, ensure entrance/exit pairs are properly set up in the compound nucleus
  for (auto &componentSegment : componentSegments_) {
    int entranceKey = componentSegment.GetEntranceKey();
    int exitKey = componentSegment.GetExitKey();

    // Ensure entrance pair is set as entrance if it exists
    if (theCNuc->IsPairKey(entranceKey)) {
      int pairNum = theCNuc->GetPairNumFromKey(entranceKey);
      theCNuc->GetPair(pairNum)->SetEntrance();
    }

    // Check if exit pair is valid (either -1 for total capture or existing pair)
    bool isValidExit = false;
    if (exitKey == -1) {
      // Check for total capture validity
      for (int i = 1; i <= theCNuc->NumPairs(); i++) {
        if (theCNuc->GetPair(i)->GetPType() == 10) {
          isValidExit = true;
          break;
        }
      }
    } else {
      isValidExit = theCNuc->IsPairKey(exitKey);
    }

    if (!isValidExit) {
      configure.outStream << "Warning: Invalid exit key " << exitKey
                          << " for component segment" << std::endl;
    }
  }

  // Initialize component segments using the same process as regular segments
  // Apply pair-dependent conversions for component segments
  for (auto &componentSegment : componentSegments_) {
    int entranceKey = componentSegment.GetEntranceKey();
    int exitKey = componentSegment.GetExitKey();

    PPair *entrancePair = nullptr;
    PPair *exitPair = nullptr;

    if (theCNuc->IsPairKey(entranceKey)) {
      entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(entranceKey));
    }
    if (theCNuc->IsPairKey(exitKey)) {
      exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(exitKey));
    }

    // Apply the same conversions that are done for test segment points
    for (int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint *point = componentSegment.GetPoint(i);

      if (entrancePair && exitPair) {
        if (entrancePair->GetPType() == 20) {
          point->ConvertDecayEnergy(exitPair);
        } else if (!componentSegment.IsCMDifferential()) {
          point->ConvertLabEnergy(entrancePair);
        } else if (componentSegment.IsCMDifferential()) {
          point->ConvertLabEnergy(entrancePair);
        }

        if (exitPair->GetPType() == 0 && componentSegment.IsDifferential() &&
            !componentSegment.IsPhase() && !componentSegment.IsAngularDist() && !componentSegment.IsCMDifferential()) {
          if (entranceKey == exitKey) {
            point->ConvertLabAngle(entrancePair);
          } else {
            point->ConvertLabAngle(entrancePair, exitPair, configure);
          }
          point->ConvertCrossSection(entrancePair, exitPair);
        }

        if (exitPair->GetPType() == 10 && componentSegment.IsDifferential() &&
            !componentSegment.IsPhase() && !componentSegment.IsAngularDist() && !componentSegment.IsCMDifferential()) {
          point->ConvertLabAngleGammas(entrancePair);
          point->ConvertCrossSectionGammas(entrancePair);
        }
      }
    }
  }

  // CalcEDependentValues for component segments
  for (auto &componentSegment : componentSegments_) {
    bool localStop = false;
    for (int i = 1; i <= componentSegment.NumPoints(); i++) {
      if (configure.stopFlag || localStop) continue;
      EPoint *point = componentSegment.GetPoint(i);
      if (!(point->IsMapped())) {
        try {
          point->CalcEDependentValues(theCNuc, configure);
        } catch (GSLException e) {
          configure.outStream << "Component segment calculation error: " << e.what() << std::endl;
          localStop = true;
        }
      }
    }
    if (configure.stopFlag || localStop) return -1;
  }

  // CalcLegendreP for component segments
  for (auto &componentSegment : componentSegments_) {
    // Get TargetEffect with Q-coefficients if applicable
    TargetEffect *effect = NULL;
    if (componentSegment.IsTargetEffect()) {
      TargetEffect *te = this->GetTargetEffect(componentSegment.GetTargetEffectNum());
      if (te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    for (int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint *point = componentSegment.GetPoint(i);
      point->CalcLegendreP(configure.maxLOrder, theCNuc, effect);
    }
  }

  // CalcCoulombAmplitude for component segments
  for (auto &componentSegment : componentSegments_) {
    for (int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint *point = componentSegment.GetPoint(i);
      point->CalcCoulombAmplitude(theCNuc);
    }
  }

  // Calculate EC amplitudes for component segments if external capture is enabled.
  // The integrals are saved to (or read back from) the same intEC.dat/intEC.extrap
  // file as the regular segments.  When previously calculated integrals are reused,
  // the (expensive) component-segment EC calculation is skipped entirely, which is
  // what makes subsequent calculations much faster.
  if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
    // Follow the decision CalculateECAmplitudes made for the data segments.
    // Reading the flag again here would resume reading a file that was just
    // rejected, and rewritten without these amplitudes, from its start.
    bool usePrevious = ecUsePrevious_;
    std::ifstream in;
    std::ofstream out;
    if (usePrevious) {
      in.open(configure.integralsfile.c_str());
      // Resume reading right after the regular-segment integrals recorded earlier.
      if (in.is_open() && ecReadPos_ > std::streampos(0)) in.seekg(ecReadPos_);
    } else {
      std::string outputfile;
      if (configure.paramMask & Config::CALCULATE_WITH_DATA)
        outputfile = configure.outputdir + "intEC.dat";
      else
        outputfile = configure.outputdir + "intEC.extrap";
      // Append the component-segment integrals after the regular-segment ones.
      out.open(outputfile.c_str(), std::ios::app);
      if (!out) configure.outStream << "Could not append to EC Amplitude File." << std::endl;
      if (!(configure.paramMask & Config::USE_API))
        configure.outStream << "Calculating EC Amplitudes for Component Segments..." << std::endl;
    }
    for (auto &componentSegment : componentSegments_) {
      int aa = theCNuc->GetPairNumFromKey(componentSegment.GetEntranceKey());
      if (theCNuc->GetPair(aa)->GetPType() != 20 && theCNuc->GetPair(aa)->IsEntrance()) {
        PPair *entrancePair = theCNuc->GetPair(aa);
        for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
          for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
            if (theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
              ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
              int ir = theCNuc->GetPairNumFromKey(componentSegment.GetExitKey());
              if (ecLevel->GetECPairNum() == ir) {
                if (!usePrevious) {
                  char segmentKeyOut[256];
                  snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d", componentSegment.GetSegmentKey());
                  if (!(configure.paramMask & Config::USE_API)) {
                    configure.outStream << "\tSegment #" << std::setw(12) << segmentKeyOut
                                        << std::setw(0) << " [                         ] 0%";
                    configure.outStream.flush();
                  }
                  int numPoints = componentSegment.NumPoints();
                  int pointIndex = 0;
                  time_t startTime = time(NULL);
                  bool localStop = false;
#pragma omp parallel for shared(configure, localStop)
                  for (int i = 1; i <= numPoints; i++) {
                    if (configure.stopFlag || localStop) continue;
                    EPoint *point = componentSegment.GetPoint(i);
                    if (!(point->IsMapped())) {
                      try {
                        point->CalculateECAmplitudes(theCNuc, configure);
                      } catch (GSLException e) {
#pragma omp critical
                        {
                          configure.outStream << e.what() << std::endl;
                          localStop = true;
                        }
                      } catch (...) {
#pragma omp critical
                        {
                          configure.outStream << "Warning: EC amplitude calculation failed for component segment point" << std::endl;
                        }
                      }
                    }
                    ++pointIndex;
                    if (difftime(time(NULL), startTime) > 0.25) {
                      startTime = time(NULL);
                      std::string progress = " [";
                      double percent = 0.;
                      for (int k = 1; k <= 25; k++) {
                        if (pointIndex >= percent * numPoints && percent < 1.) {
                          percent += 0.04;
                          progress += '*';
                        } else
                          progress += ' ';
                      }
                      progress += "] ";
                      if (!(configure.paramMask & Config::USE_API))
                        configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                            << std::setw(0) << progress << percent * 100 << '%';
                      configure.outStream.flush();
                    }
                  }
                  if (configure.stopFlag || localStop) {
                    if (out.is_open()) out.close();
                    if (in.is_open()) in.close();
                    return -1;
                  }
                  if (!(configure.paramMask & Config::USE_API))
                    configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                        << std::setw(0) << " [*************************] 100%" << std::endl;
                }
                // Write the freshly calculated amplitudes to the file, or read the
                // previously saved ones back in.  The iteration order matches the
                // regular-segment loop in CalculateECAmplitudes exactly.
                for (EPointIterator point = componentSegment.GetPoints().begin();
                     point < componentSegment.GetPoints().end(); point++) {
                  if (!(point->IsMapped())) {
                    for (int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
                      for (int ecm = 1; ecm <= entrancePair->GetDecay(ir)->GetKGroup(k)->NumECMGroups(); ecm++) {
                        if (!usePrevious) {
                          if (out.is_open()) out << point->GetECAmplitude(k, ecm) << std::endl;
                          for (EPointIterator subPoint = point->GetSubPoints().begin();
                               subPoint < point->GetSubPoints().end(); subPoint++)
                            if (out.is_open()) out << subPoint->GetECAmplitude(k, ecm) << std::endl;
                        } else {
                          complex ecAmplitude(0.0, 0.0);
                          in >> ecAmplitude;
                          point->AddECAmplitude(k, ecm, ecAmplitude);
                          for (EPointIterator subPoint = point->GetSubPoints().begin();
                               subPoint < point->GetSubPoints().end(); subPoint++) {
                            ecAmplitude = complex(0.0, 0.0);
                            in >> ecAmplitude;
                            subPoint->AddECAmplitude(k, ecm, ecAmplitude);
                          }
                        }
                        for (EPointMapIterator mappedPoint = point->GetMappedPoints().begin();
                             mappedPoint < point->GetMappedPoints().begin(); mappedPoint++) {
                          (*mappedPoint)->AddECAmplitude(k, ecm, point->GetECAmplitude(k, ecm));
                          for (int i = 1; i <= point->NumSubPoints(); i++) {
                            (*mappedPoint)->GetSubPoint(i)->AddECAmplitude(k, ecm, point->GetSubPoint(i)->GetECAmplitude(k, ecm));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (out.is_open()) {
      out.flush();
      out.close();
      // The file is complete: record what it was computed for.  (Only when
      // CalculateECAmplitudes started it, which is what set the signature.)
      std::ofstream sig;
      if (!ecOutputFile_.empty() && !ecSignature_.empty()) sig.open(ECSignaturePath(ecOutputFile_).c_str());
      if (sig.is_open()) {
        sig << "# External-capture integral signature for " << ecOutputFile_ << std::endl
            << "# (FNV-1a 64 of the energies, pairs and options the integrals depend on;"
            << " see EData::ECSignature)" << std::endl
            << "AZURE2-EC-SIGNATURE 1 " << ecSignature_ << std::endl;
      }
    }
    if (in.is_open()) in.close();
  }

  return 0;
}

/*!
 * Creates a component segment as a full copy of the base segment with different entrance/exit keys
 */
ESegment *EData::CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey) {
  // Create a copy of the base segment with different entrance/exit keys
  // Store component segments separately to avoid iterator invalidation

  // Add a copy of the base segment to the component segments vector
  this->componentSegments_.push_back(baseSegment);
  ESegment *componentSegment = &componentSegments_.back();

  // Set the new entrance/exit keys for the component segment
  componentSegment->SetEntranceKey(entranceKey);
  componentSegment->SetExitKey(exitKey);

  // Set a unique segment key to avoid cache conflicts
  // Use negative keys for component segments to distinguish from regular segments
  static int componentSegmentKeyCounter = -1000;
  componentSegment->SetSegmentKey(componentSegmentKeyCounter--);

  // Clear any existing components to avoid circular references
  componentSegment->ClearComponents();

  // The component segment now has the same data points as the base segment
  // but with different entrance/exit keys, so it will calculate different
  // theoretical cross sections when initialized

  return componentSegment;
}

/*!
 * Creates a component segment with a fixed angle for ratio denominators
 */
ESegment *EData::CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey, double fixedAngle) {
  // Create a copy of the base segment with different entrance/exit keys
  this->componentSegments_.push_back(baseSegment);
  ESegment *componentSegment = &componentSegments_.back();

  // Set the new entrance/exit keys for the component segment
  componentSegment->SetEntranceKey(entranceKey);
  componentSegment->SetExitKey(exitKey);

  // Set the fixed angle by constraining min and max to the same value
  componentSegment->SetMinAngle(fixedAngle);
  componentSegment->SetMaxAngle(fixedAngle);

  // CRITICAL: Override all point angles to the fixed angle
  // This ensures calculations use the fixed angle, not the original data angles
  std::vector<EPoint> &points = componentSegment->GetPoints();
  for (auto &point : points) {
    point.SetCMAngle(fixedAngle);
    point.SetLabAngle(fixedAngle);  // Also set lab angle to be consistent
  }

  // Set a unique segment key to avoid cache conflicts
  static int componentSegmentKeyCounter = -1000;
  componentSegment->SetSegmentKey(componentSegmentKeyCounter--);

  // Clear any existing components to avoid circular references
  componentSegment->ClearComponents();

  return componentSegment;
}

/*!
 * This function determined what points should be mapped to another to reduce
 * redundant calculations at like energies.
 */

void EData::MapData() {
  // Disable all point mapping - just clear existing mappings
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
      point->ClearMapping();
      point->ClearLocalMappedPoints();
    }
  }
}

/*!
 * Adds a TargetEffect object to the vector contained within the present object.
 */

void EData::AddTargetEffect(TargetEffect targetEffect) {
  targetEffects_.push_back(targetEffect);
}

/*!
 * Sets the normalization parameter offset in the parameter vector.
 */

void EData::SetNormParamOffset(int offset) {
  normParamOffset_ = offset;
}

/*!
 * Sets the energy shift parameter offset in the parameter vector.
 */

void EData::SetEnergyShiftParamOffset(int offset) {
  energyShiftParamOffset_ = offset;
}

/*!
 * Fills the Minuit parameter array from initial values in the EData object.
 */

void EData::FillMnParams(ROOT::Minuit2::MnUserParameters &p) {
  SetNormParamOffset(p.Params().size());
  char varname[50];
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      snprintf(varname, sizeof(varname), "segment_%d_norm", segment->GetSegmentKey());
      p.Add(varname, segment->GetNorm(), segment->GetNorm() * 0.05);
      p.SetLowerLimit(varname, 0.0);

      // Parameter settings will be applied by ParameterLimitsManager during fit
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }

  // Add energy shift parameters FOR ALL SEGMENTS (like normalization)
  SetEnergyShiftParamOffset(p.Params().size());
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    // Always add energy shift parameter, regardless of IsVaryEnergyShift()
    snprintf(varname, sizeof(varname), "segment_%d_energy_shift", segment->GetSegmentKey());
    double stepSize = (segment->GetEnergyShiftError() > 0.0) ? segment->GetEnergyShiftError() * 0.01 : 0.0005;
    p.Add(varname, segment->GetEnergyShift(), stepSize);

    if (!segment->IsVaryEnergyShift()) {
      p.Fix(varname);  // Fix parameter if not varying
    }
  }
}


/*!
 * Deletes the last segment from the segment vector.
 */

void EData::DeleteLastSegment() {
  segments_.pop_back();
}

/*!
 * Fills the Normalizations from the Minuit parameter array.
 */

void EData::FillNormsFromParams(const vector_r &p) {
  int i = GetNormParamOffset();
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      segment->SetNorm(p[i]);
      i++;
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
}

/*!
 * Fills the Energy Shifts from the Minuit parameter array.
 */

void EData::FillEnergyShiftsFromParams(const vector_r &p, EData *data, CNuc *theCNuc, const Config *configure) {
  int i = GetEnergyShiftParamOffset();
  int k = 0;
  bool anyEnergyShifted = false;

  if (data) {
    for (ESegmentIterator segment = data->GetSegments().begin(); segment < data->GetSegments().end(); segment++) {
      k++;

      // Always apply energy shift since parameters are always present for all segments
      if (segment->IsVaryEnergyShift() || p[i] != 0.0) {
        // Check if energy is the same, if so, continue
        if (segment->GetLastEnergyShift() == p[i]) {
          i++;
          continue;
        }

        segment->SetEnergyShift(p[i]);
        segment->SetLastEnergyShift(p[i]);
        segment->UpdatePointEnergiesWithShift(theCNuc, configure);
        anyEnergyShifted = true;

        // Apply the same energy shift to ALL component segments of this master segment
        if (segment->HasComponents()) {
          for (ESegment *componentSegment : segment->GetComponentSegments()) {
            if (componentSegment) {
              componentSegment->SetEnergyShift(segment->GetEnergyShift());
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
            }
          }
        }
      }

      // Apply the same energy shift to all component segments in total capture segments
      if (segment->IsTotalCapture() && segment->HasComponents()) {
        if (segment->IsVaryEnergyShift() || p[i] != 0.0) {
          // Apply energy shift to all component segments
          const std::vector<ESegment *> &componentSegments = segment->GetComponentSegments();
          for (ESegment *componentSegment : componentSegments) {
            if (componentSegment) {
              // Check if energy is the same, if so, skip the recompute.
              // NOTE: do NOT advance i here. All components of a total-capture
              // segment share this segment's single energy-shift parameter p[i];
              // i is advanced exactly once per segment at the end of the loop.
              // Advancing it inside the component loop desyncs the
              // parameter-to-segment mapping for every subsequent segment.
              if (componentSegment->GetLastEnergyShift() == p[i]) {
                continue;
              }

              componentSegment->SetEnergyShift(p[i]);
              componentSegment->SetLastEnergyShift(p[i]);
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
              anyEnergyShifted = true;
            }
          }
        }
      }

      // Always increment i since energy shift parameters are now present for ALL segments
      i++;
    }

    // Critical fix: Rebuild energy-based mapping system after energy shifts
    if (anyEnergyShifted) {
      data->MapData();
    }

    return;
  }
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    // Always apply energy shift since parameters are always present for all segments
    if (segment->IsVaryEnergyShift() || p[i] != 0.0) {
      segment->SetEnergyShift(p[i]);
      segment->UpdatePointEnergiesWithShift();
      anyEnergyShifted = true;

      // CRITICAL FIX: Apply the same energy shift to ALL component segments of this master segment
      if (segment->HasComponents()) {
        for (ESegment *componentSegment : segment->GetComponentSegments()) {
          if (componentSegment) {
            componentSegment->SetEnergyShift(segment->GetEnergyShift());
            componentSegment->UpdatePointEnergiesWithShift();
          }
        }
      }
    }

    // CRITICAL FIX: Apply the same energy shift to all component segments in total capture segments
    if (segment->IsTotalCapture() && segment->HasComponents()) {
      if (segment->IsVaryEnergyShift() || p[i] != 0.0) {
        // Apply energy shift to all component segments
        const std::vector<ESegment *> &componentSegments = segment->GetComponentSegments();
        for (ESegment *componentSegment : componentSegments) {
          if (componentSegment) {
            componentSegment->SetEnergyShift(p[i]);
            componentSegment->UpdatePointEnergiesWithShift();
            anyEnergyShifted = true;
          }
        }
      }
    }

    // Always increment i since energy shift parameters are now present for ALL segments
    i++;
  }

  // Critical fix: Rebuild energy-based mapping system after energy shifts
  if (anyEnergyShifted) {
    this->MapData();
  }
}

/*!
 * Returns a pointer to a segment specified by a position in the ESegment vector.
 */

ESegment *EData::GetSegment(int segmentNum) {
  ESegment *b = &segments_[segmentNum - 1];
  return b;
}

/*!
 * Returns a pointer to a segment based on the segment key, as opposed to a position in the ESegment vector.
 */

ESegment *EData::GetSegmentFromKey(int segmentKey) {
  int segmentNumber = 1;
  while (segmentNumber <= this->NumSegments()) {
    if (segmentKey == this->GetSegment(segmentNumber)->GetSegmentKey())
      break;
    else
      segmentNumber++;
  }
  if (segmentNumber <= this->NumSegments())
    return this->GetSegment(segmentNumber);
  else
    return NULL;
}

/*!
 * Creates a new copy of the EData object in memory and returns a pointer to the new object.
 * Used in AZURECalc function class for thread safety.
 */

EData *EData::Clone() const {
  EData *dataCopy = new EData();

  // Explicitly copy all member variables
  dataCopy->iterations_ = this->iterations_;
  dataCopy->normParamOffset_ = this->normParamOffset_;
  dataCopy->energyShiftParamOffset_ = this->energyShiftParamOffset_;
  dataCopy->isFit_ = this->isFit_;
  dataCopy->isErrorAnalysis_ = this->isErrorAnalysis_;
  dataCopy->targetEffects_ = this->targetEffects_;
  dataCopy->segments_ = this->segments_;
  dataCopy->componentSegments_ = this->componentSegments_;
  dataCopy->thmGroups_ = this->thmGroups_;

  // Build a mapping from component segment keys to cloned component segments
  std::unordered_map<int, ESegment *> clonedComponentMap;
  for (size_t i = 0; i < dataCopy->componentSegments_.size(); i++) {
    int segmentKey = dataCopy->componentSegments_[i].GetSegmentKey();
    clonedComponentMap[segmentKey] = &(dataCopy->componentSegments_[i]);
  }

  // Update component segment pointers in cloned regular segments
  for (size_t i = 0; i < dataCopy->segments_.size(); i++) {
    auto &segment = dataCopy->segments_[i];
    const auto &originalSegment = this->segments_[i];

    if (originalSegment.HasComponents()) {
      // Clear the old component pointers (they point to original EData)
      segment.ClearComponents();

      // Preserve the operation type
      segment.SetOperationType(originalSegment.GetOperationType());

      // Add the remapped component segments by finding them in the cloned componentSegments_
      for (const auto &originalComponent : originalSegment.GetComponentSegments()) {
        // Find this component in the original componentSegments_ to get its position
        int componentIndex = -1;
        for (size_t j = 0; j < this->componentSegments_.size(); j++) {
          if (&this->componentSegments_[j] == originalComponent) {
            componentIndex = j;
            break;
          }
        }

        // Add the corresponding cloned component
        if (componentIndex >= 0 && componentIndex < (int)dataCopy->componentSegments_.size()) {
          segment.AddComponentSegment(&dataCopy->componentSegments_[componentIndex]);
        }
      }
    }
  }

  for (EDataIterator data = dataCopy->begin(); data != dataCopy->end(); data++) {
    data.point()->SetParentData(dataCopy);
    data.point()->ClearLocalMappedPoints();
    data.point()->ClearECAmplitudes();
  }
  // Also fix parentData for component segment points
  for (size_t i = 0; i < dataCopy->componentSegments_.size(); i++) {
    for (int j = 1; j <= dataCopy->componentSegments_[i].NumPoints(); j++) {
      dataCopy->componentSegments_[i].GetPoint(j)->SetParentData(dataCopy);
      dataCopy->componentSegments_[i].GetPoint(j)->ClearLocalMappedPoints();
      dataCopy->componentSegments_[i].GetPoint(j)->ClearECAmplitudes();
    }
  }
  for (EDataIterator data = dataCopy->begin(); data != dataCopy->end(); data++) {
    if (data.point()->IsMapped()) {
      EnergyMap pointMap = data.point()->GetMap();
      dataCopy->GetSegment(pointMap.segment)->GetPoint(pointMap.point)->AddLocalMappedPoint(&*data.point());
    }
  }

  return dataCopy;
}

/*!
 * Returns a pointer to the specified TargetEffect object.
 */

TargetEffect *EData::GetTargetEffect(int effectNumber) {
  TargetEffect *temp;
  if (effectNumber <= targetEffects_.size())
    temp = &targetEffects_[effectNumber - 1];
  else
    return temp = NULL;
  return temp;
}

/*!
 * Returns an EDataIterator referring to the first data point in the set.
 */

EDataIterator EData::begin() {
  return EDataIterator(&segments_);
}

/*!
 * Returns an EDataIterator referring to one object past the last data point in the set.
 */

EDataIterator EData::end() {
  EDataIterator it(&segments_);
  return it.SetEnd();
}

/*!
 * Returns a reference to the vector of ESegment objects.
 */

std::vector<ESegment> &EData::GetSegments() {
  return segments_;
}

// ---------------------------------------------------------------------------
// THM experiments (<thm> experiment[<name>] ...): see ThmExperiment.h.

int EData::BuildThmGroups(const Config &configure, CNuc *theCNuc, int numLines) {
  thmGroups_.clear();
  const double amu = 931.49410242;  // MeV/u (CODATA 2018)
  for (const ThmExperiment &x : configure.thm.experiments) {
    const std::string where = "ERROR: <thm> experiment[" + x.name + "]: ";
    ThmGroup group;
    group.name = x.name;
    group.terms = x.backgroundTerms;
    for (int key : x.segments) {
      if (key > numLines) {
        configure.outStream << where << "segment " << key << ": <segmentsData> has only " << numLines
                            << " line(s)." << std::endl;
        return -1;
      }
      int index = 0;
      for (int s = 1; s <= NumSegments(); s++)
        if (GetSegment(s)->GetSegmentKey() == key) index = s;
      if (index == 0) {
        configure.outStream << "WARNING: <thm> experiment[" << x.name << "]: segment line " << key
                            << " of <segmentsData> is not in use; it is left out." << std::endl;
        continue;
      }
      ESegment *segment = GetSegment(index);
      if (!segment->IsTHM()) {
        configure.outStream << where << "segment " << key << " is not a THM segment (isDiff < 10)." << std::endl;
        return -1;
      }
      if (!segment->IsVaryNorm()) {
        configure.outStream << where << "segment " << key
                            << " has a fixed norm; the segments of an experiment share one free "
                               "(profiled) norm, so free it."
                            << std::endl;
        return -1;
      }
      group.segments.push_back(index);
    }
    if (group.segments.empty()) {
      configure.outStream << "WARNING: <thm> experiment[" << x.name << "] has no segment in use; ignored."
                          << std::endl;
      continue;
    }
    int points = 0;
    for (int s : group.segments)
      for (int p = 1; p <= GetSegment(s)->NumPoints(); p++)
        if (GetSegment(s)->GetPoint(p)->GetCMCrossSectionError() != 0.0) points++;
    if (points <= 1 + group.terms) {
      configure.outStream << where << points << " point(s) with an error for " << 1 + group.terms
                          << " profiled linear parameter(s) (norm and background " << ThmExperiment::BackgroundName(group.terms)
                          << "); it needs more." << std::endl;
      return -1;
    }
    group.trivial = group.segments.size() == 1 && group.terms == 0;
    group.pairKey = GetSegment(group.segments[0])->GetEntranceKey();
    for (int s : group.segments)
      if (GetSegment(s)->GetEntranceKey() != group.pairKey) group.pairKey = 0;

    std::ostringstream summary;
    summary << "THM experiment '" << x.name << "': segment" << (group.segments.size() > 1 ? "s " : " ");
    for (size_t k = 0; k < group.segments.size(); k++)
      summary << (k ? "," : "") << GetSegment(group.segments[k])->GetSegmentKey();
    summary << (group.segments.size() > 1 ? " share one profiled norm" : " with a profiled norm")
            << ", background " << ThmExperiment::BackgroundName(group.terms) << ".";
    configure.outStream << summary.str() << std::endl;

    {
      ThmDistortion::Kinematics dk;  // for distortion=coulomb|optical
      if (x.hasKinematics) {
        // The entrance pair x + A of the THM segments, and which of beam and
        // target is the Trojan horse a = x + s.
        int pairKey = GetSegment(group.segments[0])->GetEntranceKey();
        for (int s : group.segments)
          if (GetSegment(s)->GetEntranceKey() != pairKey) {
            configure.outStream << where << "beam/target/spectator describe one reaction, but its segments "
                                   "have different entrance pairs."
                                << std::endl;
            return -1;
          }
        PPair *pair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(pairKey));
        int Z[2] = {pair->GetZ(1), pair->GetZ(2)};
        int A[2] = {(int)std::lround(pair->GetM(1)), (int)std::lround(pair->GetM(2))};
        const ThmNuclide &b = x.beam, &t = x.target, &sp = x.spectator;
        // horse: 0 beam, 1 target; other: index of the pair nucleus that is the other one.
        int horse = -1, other = -1;
        for (int h = 0; h < 2 && horse < 0; h++) {
          const ThmNuclide &th = h == 0 ? b : t, &tg = h == 0 ? t : b;
          for (int k = 0; k < 2; k++)
            if (tg.Z == Z[k] && tg.A == A[k] && th.Z - sp.Z == Z[1 - k] && th.A - sp.A == A[1 - k]) {
              horse = h;
              other = k;
              break;
            }
        }
        if (horse < 0) {
          configure.outStream << where << "beam " << b.name << " + target " << t.name << " with spectator "
                              << sp.name << " does not give the entrance pair of its segments (Z,A) = ("
                              << Z[0] << "," << A[0] << ") + (" << Z[1] << "," << A[1]
                              << "): one of beam/target must be a nucleus of the pair and the other the "
                                 "second nucleus plus the spectator."
                              << std::endl;
          return -1;
        }
        const ThmNuclide &th = horse == 0 ? b : t, &nA = horse == 0 ? t : b;
        const ThmNuclide *tabX = ThmNuclide::Find(th.Z - sp.Z, th.A - sp.A);
        double mX = tabX ? tabX->mass : pair->GetM(2 - other);  // the pair nucleus that is x
        double mA = nA.mass;
        double bind = (mX + sp.mass - th.mass) * amu;
        dk.Za = th.Z;
        dk.ZA = nA.Z;
        dk.Zs = sp.Z;
        dk.Zx = th.Z - sp.Z;
        dk.Aa = th.A;
        dk.AA = nA.A;
        dk.As = sp.A;
        dk.Ax = th.A - sp.A;
        dk.ma = th.mass;
        dk.mA = nA.mass;
        dk.ms = sp.mass;
        dk.mx = mX;
        dk.horseIsBeam = horse == 0;
        dk.mBeam = b.mass;
        dk.mTarget = t.mass;
        dk.beamEnergy = x.beamEnergy;
        dk.bind = bind;
        // Quasi-free x + A energy: the spectator keeps the Trojan horse's
        // velocity (horse = beam) or stays at rest (horse = target).
        double exa = horse == 0 ? x.beamEnergy * mX / th.mass * mA / (mX + mA) : x.beamEnergy * mX / (mA + mX);
        std::ostringstream k;
        k.precision(6);
        k << "  " << b.name << " + " << t.name << " at " << x.beamEnergy << " MeV (lab), Trojan horse "
          << th.name << " = x + " << sp.name << ", B(x+s) = " << bind << " MeV; quasi-free E(x+A) = " << exa
          << " MeV, E_qf = E(x+A) - B = " << exa - bind << " MeV.";
        configure.outStream << k.str() << std::endl;
        // The vertex takes B from the entrance pair's channel lines (field 32),
        // the kinematics from the masses; say so when they disagree.
        if (std::fabs(pair->GetBindingEnergy() - bind) > 1.0e-3) {
          std::ostringstream w;
          w.precision(6);
          w << "WARNING: <thm> experiment[" << x.name << "]: B(x+s) from the masses of " << th.name << " = x + "
            << sp.name << " is " << bind << " MeV, but the entrance pair " << pairKey
            << " carries B = " << pair->GetBindingEnergy()
            << " MeV (field 32 of its channel lines); the THM vertex uses field 32, the kinematics of this "
               "experiment (the quasi-free energy above, E_sF of the line shape and of the distortion "
               "factor) use the masses.";
          configure.outStream << w.str() << std::endl;
        }

        if (x.psKind != ThmExperiment::PS_DELTA) {
          // Spectator-momentum window (ThmLineshape.h ThmSpectatorWindow).
          if (configure.thm.SpectatorEnergy(pairKey) != 0.0) {
            configure.outStream << where << "a ps window and spectatorEnergy both set the spectator motion of "
                                   "entrance pair "
                                << pairKey << "; use one (ps=delta keeps spectatorEnergy)." << std::endl;
            return -1;
          }
          double muSx = mX * sp.mass / (mX + sp.mass) * amu;
          std::shared_ptr<ThmSpectatorWindow> window = std::make_shared<ThmSpectatorWindow>();
          std::string why = BuildThmSpectatorWindow(x, muSx, *window);
          if (!why.empty()) {
            configure.outStream << where << "ps: " << why << "." << std::endl;
            return -1;
          }
          std::ostringstream w;
          w.precision(6);
          w << "  Spectator-momentum window: " << window->description << "; mu_sx = " << muSx
            << " MeV, T_s = p_s^2/2mu_sx from " << window->es.front() << " to " << window->es.back()
            << " MeV, <T_s> = " << window->MeanEs() << " MeV.";
          configure.outStream << w.str() << std::endl;
          group.window = window;
          // With vertexModel=dw the DW vertex puts the nodes on the reachable
          // part of the window itself (ThmDwVertex); the plane-wave nodes are
          // not used.
          if (!x.vertexDW)
            for (int s : group.segments) GetSegment(s)->SetThmSpectatorWindow(window);
        }

        if (x.lineshape) {
          // Coulomb line shape of the spectator (ThmLineshape.h).  The level
          // energies and widths it uses are the observed ones of Brune.
          if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
            configure.outStream << where << "lineshape=on uses the observed level energies and widths "
                                   "as the resonance poles; it needs the Brune parameterization."
                                << std::endl;
            return -1;
          }
          std::shared_ptr<ThmLineshape> shape = std::make_shared<ThmLineshape>();
          shape->experiment = x.name;
          shape->spectator = sp.name;
          shape->Zs = sp.Z;
          shape->ms = sp.mass;
          shape->ZF = pair->GetZ(1) + pair->GetZ(2);
          shape->mF = pair->GetM(1) + pair->GetM(2);
          shape->eAA = x.beamEnergy * t.mass / (b.mass + t.mass);
          shape->bind = bind;
          // Every data point must leave the spectator some energy.
          double eMax = -1.0e300;
          for (int s : group.segments)
            for (int p = 1; p <= GetSegment(s)->NumPoints(); p++)
              eMax = std::max(eMax, GetSegment(s)->GetPoint(p)->GetCMEnergy());
          if (!(shape->EsF(eMax) > 0.0)) {
            configure.outStream << where << "lineshape=on: at E = " << eMax
                                << " MeV the spectator has no energy left (E_sF = E_aA - B - E = "
                                << shape->eAA << " - " << bind << " - " << eMax << " MeV <= 0); check Ebeam."
                                << std::endl;
            return -1;
          }
          std::ostringstream l;
          l.precision(6);
          l << "  Coulomb line shape on: E_sF = E_aA - B - E with E_aA = " << shape->eAA << " MeV, eta_0 = "
            << shape->Eta0(eMax) << " at the highest point energy (E = " << eMax << " MeV)"
            << (sp.Z == 0 ? "; the spectator is neutral, so N_C = 1." : ".");
          configure.outStream << l.str() << std::endl;
          for (int s : group.segments) {
            int key = GetSegment(s)->GetExitKey();
            bool seen = false;
            for (const ThmLineshape::Exit &e : shape->exits) seen = seen || e.pairKey == key;
            if (seen || !theCNuc->IsPairKey(key)) continue;
            PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(key));
            int light = exitPair->GetM(1) <= exitPair->GetM(2) ? 1 : 2;
            ThmLineshape::Exit e;
            e.pairKey = key;
            e.Zb = exitPair->GetZ(light);
            e.ZB = exitPair->GetZ(3 - light);
            e.mb = exitPair->GetM(light);
            e.mB = exitPair->GetM(3 - light);
            e.q = pair->GetSepE() + pair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
            shape->exits.push_back(e);
          }
          group.lineshape = shape;
          for (int s : group.segments) GetSegment(s)->SetThmLineshape(shape);
        }
      }
      if (x.vertexDW) {
        // Distorted-wave entrance vertex (ThmDwVertex.h): replaces M_l in the
        // HOES amplitude of every segment of the experiment; R(E) is not applied.
        int pairKey = GetSegment(group.segments[0])->GetEntranceKey();
        int pairNum = theCNuc->GetPairNumFromKey(pairKey);
        PPair *pair = theCNuc->GetPair(pairNum);
        // coulombIntegral=1 is refused with it by Config::ReadThmBlock
        // (CheckThmCoulombConsistency).
        if (configure.thm.coherentL) {
          configure.outStream << where << "vertexModel=dw sums the entrance partial waves (and their projections) "
                                 "incoherently, as the angle-integrated observable requires; entranceL=coherent "
                                 "cannot be combined with it."
                              << std::endl;
          return -1;
        }
        if (configure.thm.SpectatorEnergy(pairKey) != 0.0) {
          configure.outStream << where << "with vertexModel=dw the spectator kinematics come from Ebeam and the "
                                 "spectator direction; spectatorEnergy for entrance pair "
                              << pairKey << " must be 0." << std::endl;
          return -1;
        }
        std::vector<int> ls;
        for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
          JGroup *jg = theCNuc->GetJGroup(j);
          if (!jg->IsInRMatrix()) continue;
          for (int ch = 1; ch <= jg->NumChannels(); ch++)
            if (jg->GetChannel(ch)->GetPairNum() == pairNum) ls.push_back(jg->GetChannel(ch)->GetL());
        }
        std::vector<double> energies;
        double eLo = 1.0e300, eHi = -1.0e300;
        for (int s : group.segments)
          for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
            double e = GetSegment(s)->GetPoint(p)->GetCMEnergy();
            energies.push_back(e);
            eLo = std::min(eLo, e);
            eHi = std::max(eHi, e);
          }
        std::shared_ptr<ThmDwVertex> v = std::make_shared<ThmDwVertex>();
        double eAA = dk.beamEnergy * dk.mTarget / (dk.mBeam + dk.mTarget);
        if (!(eAA - dk.bind - eHi > 0.0)) {
          configure.outStream << where << "vertexModel=dw: at E = " << eHi
                              << " MeV the spectator has no energy left (E_sF = E_aA - B - E = " << eAA << " - "
                              << dk.bind << " - " << eHi << " MeV <= 0); check Ebeam." << std::endl;
          return -1;
        }
        double gridHi = std::min(eHi + 0.3, eAA - dk.bind - 0.5 * (eAA - dk.bind - eHi));
        std::string why = ls.empty() ? std::string("the entrance pair has no channel in the R matrix")
                                     : v->Build(x, dk, pair->GetChRad(), ls, eLo - 0.3, gridHi, energies);
        if (!why.empty()) {
          configure.outStream << where << "vertexModel=dw: " << why << "." << std::endl;
          return -1;
        }
        std::ostringstream l;
        l.precision(6);
        l << "  Entrance vertex: " << v->description << ".\n"
          << "  alpha = m_A/m_F = " << v->alpha << ", beta = m_s/m_a = " << v->beta << ", k_aA = " << v->dist.aa.k
          << " fm^-1, eta_aA = " << v->dist.aa.eta << ", kappa = " << v->dist.kappa << " fm^-1, eta_b = "
          << v->dist.etaB << "; channel radius " << v->radius << " fm, l =";
        for (int lv : v->lvals) l << " " << lv;
        l << "; L <= " << v->laMax << " (a + A), " << v->lsMax << " (s + F); u to " << v->uMax << " fm on "
          << v->uNodes << " x " << v->cNodes << " nodes; grid " << v->gridLo << " to "
          << v->gridLo + (v->nE - 1) * v->gridStep << " MeV (" << v->nE << " energies, " << v->nNodes
          << " node(s) each), " << v->buildSeconds << " s.\n"
          << "  The distortion factor R(E) is not applied: the DW vertex carries the energy dependence.";
        for (const std::string &w : v->dist.warnings)
          l << "\nWARNING: <thm> experiment[" << x.name << "]: " << w;
        for (int s : group.segments)
          if (GetSegment(s)->GetThmWeight())
            l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << GetSegment(s)->GetSegmentKey()
              << " also has weight[" << GetSegment(s)->GetSegmentKey()
              << "]=; it multiplies the model with the DW vertex.";
        configure.outStream << l.str() << std::endl;
        group.dwVertex = v;
        for (int s : group.segments) GetSegment(s)->SetThmDwVertex(v);
      } else if (x.distortion != ThmExperiment::DIST_NONE) {
        // Distortion factor R(E) (ThmDistortion.h): multiplies the model of
        // every segment of the experiment before the folding.
        double eLo = 1.0e300, eHi = -1.0e300;
        for (int s : group.segments)
          for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
            eLo = std::min(eLo, GetSegment(s)->GetPoint(p)->GetCMEnergy());
            eHi = std::max(eHi, GetSegment(s)->GetPoint(p)->GetCMEnergy());
          }
        std::shared_ptr<ThmDistortion> d = std::make_shared<ThmDistortion>();
        d->experiment = x.name;
        std::ostringstream l;
        l.precision(6);
        if (x.distortion == ThmExperiment::DIST_TABLE) {
          d->kind = ThmDistortion::TABLE;
          d->table = x.distortionWeights;
          const ThmWeightTable &table = *d->table;
          if (!table.Covers(eLo) || !table.Covers(eHi)) {
            configure.outStream << where << "distortion: the points span E_cm = " << eLo << " to " << eHi
                                << " MeV, beyond the table '" << table.name << "' [" << table.e.front() << ", "
                                << table.e.back() << "] MeV." << std::endl;
            return -1;
          }
          d->description = "table " + table.name;
          l << "  Distortion factor: w(E) from the table '" << table.name << "' multiplies the model.";
        } else {
          d->kin = dk;
          d->eAA = dk.beamEnergy * dk.mTarget / (dk.mBeam + dk.mTarget);
          // Every data point must be reachable (the grid then extends 0.5
          // MeV beyond the points, short of the spectator's threshold).
          d->angleKind = x.angleKind == 1 ? ThmDistortion::LAB : x.angleKind == 2 ? ThmDistortion::CM : ThmDistortion::QF;
          d->angle = x.angle;
          d->sf.kind = x.distortion == ThmExperiment::DIST_OPTICAL && x.opticalSF.kind == 0
                           ? ThmDistortion::Channel::PLANE
                           : ThmDistortion::Channel::POINT_COULOMB;
          d->sf.mu = dk.ms * (dk.mx + dk.mA) / (dk.ms + dk.mx + dk.mA) * uconv;
          d->vcm = std::sqrt(2.0 * dk.mBeam * uconv * dk.beamEnergy) / ((dk.mBeam + dk.mTarget) * uconv);
          for (int s : group.segments)
            for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
              std::string why = d->CheckEnergy(GetSegment(s)->GetPoint(p)->GetCMEnergy());
              if (!why.empty()) {
                configure.outStream << where << "distortion: " << why << "." << std::endl;
                return -1;
              }
            }
          double gridHi = std::min(eHi + 0.5, d->eAA - dk.bind - 0.5 * d->EsF(eHi));
          d->dataLo = eLo;
          d->dataHi = eHi;
          std::string why = d->Build(x, dk, eLo - 0.5, gridHi, 0.5 * (eLo + eHi));
          if (!why.empty()) {
            configure.outStream << where << "distortion: " << why << "." << std::endl;
            return -1;
          }
          ThmDistortion::Point lo = d->Evaluate(eLo), hi = d->Evaluate(eHi);
          l << "  Distortion factor R(E), zero-range DWBA: " << d->description << ".\n"
            << "  k_aA = " << d->aa.k << " fm^-1, eta_aA = " << d->aa.eta << ", kappa = " << d->kappa
            << " fm^-1, eta_b = " << d->etaB << ", beta = m_s/m_a = " << d->beta << "; E_ref = " << d->eRef
            << " MeV, grid " << d->gridLo << " to " << d->gridLo + (d->lnR.size() - 1) * d->gridStep
            << " MeV, radial step " << d->h << " fm to " << (d->n - 1) * d->h << " fm, l <= " << d->uAA.size() - 1
            << ".\n"
            << "  R = " << d->R(lo) << " at E = " << eLo << " MeV (E_sF = " << lo.esf << ", eta_sF = " << lo.etasf
            << ", theta_cm = " << lo.thetaCm << " deg), " << d->R(hi) << " at E = " << eHi << " MeV (E_sF = "
            << hi.esf << ", eta_sF = " << hi.etasf << ", theta_cm = " << hi.thetaCm << " deg).";
          if (d->pwSignChange && d->ratioPW)
            l << "\nWARNING: <thm> experiment[" << x.name << "]: the plane-wave amplitude M_PW changes sign on "
                 "the grid (a node of the momentum distribution at this angle); R = |M/M_PW|^2 is singular "
                 "there (distortionRatio=dw avoids it).";
          for (const std::string &w : d->warnings) l << "\nWARNING: <thm> experiment[" << x.name << "]: " << w;
          if (d->tailWorst > 1.0e-8)
            l << "\nWARNING: <thm> experiment[" << x.name << "]: the radial integrals are cut at r = "
              << (d->n - 1) * d->h << " fm with a remainder up to " << d->tailWorst << " of |M|.";
        }
        for (int s : group.segments)
          if (GetSegment(s)->GetThmWeight())
            l << "\nWARNING: <thm> experiment[" << x.name << "]: segment " << GetSegment(s)->GetSegmentKey()
              << " also has weight[" << GetSegment(s)->GetSegmentKey()
              << "]=; both multiply its model (the distortion factor and the weight).";
        configure.outStream << l.str() << std::endl;
        group.distortion = d;
        for (int s : group.segments) GetSegment(s)->SetThmDistortion(d);
      }
    }
    if (x.hasTheta) {
      // Angular window of the exit pair (ThmAngular.h): the model of every
      // segment is dsigma/dOmega averaged over theta_cm in the window.
      if (configure.thm.coherentL) {
        configure.outStream << where << "theta= computes the interference of the entrance partial waves "
                               "exactly (at fixed angle they interfere); entranceL=coherent is an approximation "
                               "of the angle-integrated observable and cannot be combined with it."
                            << std::endl;
        return -1;
      }
      int maxLp = 0;
      for (int j = 1; j <= theCNuc->NumJGroups(); j++)
        for (int ch = 1; ch <= theCNuc->GetJGroup(j)->NumChannels(); ch++)
          maxLp = std::max(maxLp, theCNuc->GetJGroup(j)->GetChannel(ch)->GetL());
      if (2 * maxLp > ThmAngleWindow::kMaxL) {
        configure.outStream << where << "theta= carries Legendre orders up to " << ThmAngleWindow::kMaxL
                            << "; the model has a channel with l = " << maxLp << "." << std::endl;
        return -1;
      }
      std::shared_ptr<ThmAngleWindow> w = std::make_shared<ThmAngleWindow>();
      w->experiment = x.name;
      BuildThmAngleWindow(x.thetaMin, x.thetaMax, *w);
      std::ostringstream a;
      a.precision(6);
      a << "  Angular window: theta_cm = " << x.thetaMin << "-" << x.thetaMax
        << " deg (exit particle 1 relative to 2, from p_xA = entrance particle 1 relative to 2); the model "
           "is the HOES dsigma/dOmega averaged over it (4 pi times it is the angle-integrated cross section "
           "for 0-180).";
      configure.outStream << a.str() << std::endl;
      group.angle = w;
      for (int s : group.segments) GetSegment(s)->SetThmAngleWindow(w);
    }
    thmGroups_.push_back(group);
  }
  return 0;
}

int EData::ThmGroupOf(int i) {
  for (size_t g = 0; g < thmGroups_.size(); g++) {
    if (thmGroups_[g].trivial) continue;
    const std::vector<int> &s = thmGroups_[g].segments;
    if (std::find(s.begin(), s.end(), i) != s.end()) return GetSegment(i)->IsProfiledNorm() ? (int)g : -1;
  }
  return -1;
}

bool EData::IsLastOfThmGroup(int g, int i) {
  const std::vector<int> &s = thmGroups_[g].segments;
  for (size_t k = s.size(); k-- > 0;)
    if (GetSegment(s[k])->IsProfiledNorm()) return s[k] == i;
  return false;
}

namespace {
// The points of the profiled segments of a group, concatenated in order.
void GatherThmGroup(EData *data, const std::vector<int> &segments, std::vector<double> &m,
                    std::vector<double> &d, std::vector<double> &e, std::vector<double> &energy) {
  m.clear();
  d.clear();
  e.clear();
  energy.clear();
  for (int s : segments) {
    ESegment *seg = data->GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (!pt) continue;
      m.push_back(pt->GetFitCrossSection());
      d.push_back(pt->GetCMCrossSection());
      e.push_back(pt->GetCMCrossSectionError());
      energy.push_back(pt->GetCMEnergy());
    }
  }
}
}  // namespace

double EData::ProfileThmGroup(int g) {
  ThmGroup &group = thmGroups_[g];
  std::vector<double> m, d, e, energy;
  GatherThmGroup(this, group.segments, m, d, e, energy);
  group.profile = SolveThmProfile(m, d, e, energy, group.terms);
  const ThmProfile &p = group.profile;
  double total = 0.0;
  size_t k = 0;
  for (int s : group.segments) {
    ESegment *seg = GetSegment(s);
    if (!seg->IsProfiledNorm()) continue;
    seg->SetNorm(p.Norm());
    double chi = 0.0;
    for (int q = 1; q <= seg->NumPoints(); q++) {
      if (!seg->GetPoint(q)) continue;
      double r = p.Residual(m[k], d[k], e[k], energy[k]);
      chi += r * r;
      k++;
    }
    seg->SetSegmentChiSquared(chi);
    total += chi;
  }
  return total;
}

void EData::ThmGroupResiduals(int g, int i, std::vector<double> &out) {
  out.clear();
  const ThmProfile &p = thmGroups_[g].profile;
  ESegment *seg = GetSegment(i);
  for (int q = 1; q <= seg->NumPoints(); q++) {
    EPoint *pt = seg->GetPoint(q);
    if (!pt) continue;
    out.push_back(p.Residual(pt->GetFitCrossSection(), pt->GetCMCrossSection(), pt->GetCMCrossSectionError(),
                             pt->GetCMEnergy()));
  }
}

double EData::ThmBackgroundAt(int i, double energy) {
  int g = ThmGroupOf(i);
  if (g < 0 || thmGroups_[g].terms == 0) return 0.0;
  return thmGroups_[g].profile.Background(energy);
}

std::vector<ThmExperimentReport> EData::ThmExperimentReports() {
  std::vector<ThmExperimentReport> out;
  for (size_t g = 0; g < thmGroups_.size(); g++) {
    ThmGroup &group = thmGroups_[g];
    ThmExperimentReport r;
    r.name = group.name;
    r.background = ThmExperiment::BackgroundName(group.terms);
    ThmProfile p;
    if (group.trivial || ThmGroupOf(group.segments[0]) < 0) {
      // Profiled per segment (ESegment::ProfileNormChiSquared) or not at all:
      // the same closed form, recomputed here for its uncertainty.
      std::vector<double> m, d, e, energy;
      GatherThmGroup(this, group.segments, m, d, e, energy);
      p = SolveThmProfile(m, d, e, energy, group.terms);
    } else {
      p = group.profile;
    }
    r.points = 0;
    r.chi2 = 0.0;
    for (int s : group.segments) {
      ESegment *seg = GetSegment(s);
      if (!seg->IsProfiledNorm()) continue;
      r.segments.push_back(seg->GetSegmentKey());
      r.chi2 += seg->GetSegmentChiSquared();
    }
    r.points = p.points;
    p.Reported(r.value, r.covariance);
    r.status = r.segments.empty() ? "not profiled: no segment of it has a free norm" : p.status;
    out.push_back(r);
  }
  return out;
}

void EData::WriteThmExperiments(const Config &configure) {
  std::string file = configure.outputdir + "thm_experiments.out";
  std::ofstream out(file.c_str());
  if (!out) {
    configure.outStream << "Could not write " << file << "." << std::endl;
    return;
  }
  out << "# THM experiments (<thm> experiment[<name>] lines): segments sharing one profiled norm and\n"
         "# a background b(E) = b0 + b1 E + b2 E^2 added to the folded model (model units, E the c.m.\n"
         "# energy of the entrance pair in MeV).  norm multiplies the data, as the segment norms of\n"
         "# normalizations.out.  Uncertainties and covariance are those of the closed-form profile\n"
         "# at fixed R-matrix parameters (inverse normal matrix of the linear least squares, not\n"
         "# scaled by chi2/nu).\n";
  out.precision(10);
  out << std::scientific;
  for (const ThmExperimentReport &r : ThmExperimentReports()) {
    out << "\nexperiment: " << r.name << "\nsegments:";
    for (size_t k = 0; k < r.segments.size(); k++) out << (k ? "," : " ") << r.segments[k];
    out << "\nbackground: " << r.background << "\npoints: " << r.points << "\nchi2: " << r.chi2
        << "\nstatus: " << r.status << "\n";
    const char *names[4] = {"norm", "b0", "b1", "b2"};
    int q = 1 + (r.background == "none" ? 0 : r.background == "const" ? 1 : r.background == "linear" ? 2 : 3);
    for (int i = 0; i < q; i++)
      out << std::left << std::setw(6) << names[i] << std::right << std::setw(18) << r.value[i]
          << std::setw(18) << std::sqrt(std::max(0.0, r.covariance[i * 4 + i])) << "\n";
    out << "covariance:\n";
    for (int i = 0; i < q; i++) {
      out << std::left << std::setw(6) << names[i] << std::right;
      for (int j = 0; j < q; j++) out << std::setw(18) << r.covariance[i * 4 + j];
      out << "\n";
    }
    // Spectator-momentum window (ps=...): the nodes.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.window) continue;
      const ThmSpectatorWindow &w = *group.window;
      out << "ps: " << w.description << "; mu_sx = " << w.muSx << " MeV\n"
          << "# The model at E is the average of the HOES cross section over the spectator momentum\n"
          << "# p_s, weight |phi(p_s)|^2 p_s^2 (a table: its w(p_s)), incoherent; node k adds\n"
          << "# T_s = p_s^2/2mu_sx to E + B in the vertex.  Columns: p_s (MeV/c), weight, T_s (MeV).\n";
      for (size_t k = 0; k < w.p.size(); k++)
        out << "ps_node" << std::setw(18) << w.p[k] << std::setw(18) << w.weight[k] << std::setw(18) << w.es[k]
            << "\n";
      out << std::left << std::setw(16) << "<T_s>" << std::right << std::setw(18) << w.MeanEs() << "\n";
      // Memory of the per-point node tables (EPoint::ThmPsTable): points
      // including the folding sub-points, the stored entrance channels summed
      // over the points (x nodes = entries), and their bytes.
      size_t tablePoints = 0, tableBytes = 0;
      for (int s : group.segments)
        for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
          EPoint *point = GetSegment(s)->GetPoint(p);
          tablePoints++;
          tableBytes += point->ThmPsTableBytes();
          for (int q = 1; q <= point->NumSubPoints(); q++) {
            tablePoints++;
            tableBytes += point->GetSubPoint(q)->ThmPsTableBytes();
          }
        }
      out << "# Node tables of the vertex: points (with folding sub-points), bytes.\n"
          << "ps_table" << std::setw(18) << (double)tablePoints << std::setw(18) << (double)tableBytes << "\n";
    }
    // Distortion factor (distortion=...): R at the lowest point, E_ref and the highest point.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.distortion) continue;
      const ThmDistortion &d = *group.distortion;
      out << "distortion: " << d.description << "\n";
      if (d.kind == ThmDistortion::TABLE) continue;
      double eLo = 1.0e300, eHi = -1.0e300;
      for (int s : group.segments)
        for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
          eLo = std::min(eLo, GetSegment(s)->GetPoint(p)->GetCMEnergy());
          eHi = std::max(eHi, GetSegment(s)->GetPoint(p)->GetCMEnergy());
        }
      out << "# Zero-range DWBA transfer amplitude M(E) = <chi(-)_sF phi_sx chi(+)_aA(beta r)> (Mukhamedzhanov &\n"
          << "# Pang PRC 99 (2019) 064618 eqs. 20-24; Mukhamedzhanov arXiv:2609.04498 eqs. 22-30), M_PW its\n"
          << "# plane-wave limit; the model is multiplied by R(E) before folding (PWA-extracted S* / R).\n"
          << "# E_aA = " << d.eAA << " MeV, B = " << d.kin.bind << " MeV, k_aA = " << d.aa.k << " fm^-1, eta_aA = "
          << d.aa.eta << ", kappa = " << d.kappa << " fm^-1, eta_b = " << d.etaB << ", beta = " << d.beta << "\n"
          << "# Columns: E, E_sF (MeV), eta_sF, theta_cm (deg), |M|^2, |M_PW|^2, R, l_max.\n";
      for (double e : {eLo, d.eRef, eHi}) {
        ThmDistortion::Point p = d.Evaluate(e);
        out << "distortion_point" << std::setw(18) << e << std::setw(18) << p.esf << std::setw(18) << p.etasf
            << std::setw(18) << p.thetaCm << std::setw(18) << std::norm(p.m) << std::setw(18) << p.mpw * p.mpw
            << std::setw(18) << (p.ok ? d.R(p) : 0.0) << std::setw(6) << p.lmax << "\n";
      }
    }
    // Distorted-wave entrance vertex (vertexModel=dw): the Gram entries at the
    // lowest, the middle and the highest point.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.dwVertex) continue;
      const ThmDwVertex &v = *group.dwVertex;
      double eLo = 1.0e300, eHi = -1.0e300;
      for (int s : group.segments)
        for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
          eLo = std::min(eLo, GetSegment(s)->GetPoint(p)->GetCMEnergy());
          eHi = std::max(eHi, GetSegment(s)->GetPoint(p)->GetCMEnergy());
        }
      out << "vertex: " << v.description << "\n"
          << "# Surface term of the prior-form DWBA (Mukhamedzhanov PRC 84 (2011) 044616; Mukhamedzhanov,\n"
          << "# Kadyrov & Pang EPJA 56 (2020) 233 eqs. 28-32) replaces M_l; R(E) is not applied.\n"
          << "# G = (4pi/(2l+1)) sum_m (s_m, d_m)^+ (s_m, d_m), (s_m, d_m) = (S_lm(a), a S_lm'(a))/(4pi phi~(q)),\n"
          << "# |M_l|^2 = c^+ G c with c = (B - 1, -1); plane waves: G11 = j_l(pa)^2, G22 = (pa j_l'(pa))^2.\n"
          << "# k_aA = " << v.dist.aa.k << " fm^-1, eta_aA = " << v.dist.aa.eta << ", kappa = " << v.dist.kappa
          << " fm^-1, alpha = " << v.alpha << ", beta = " << v.beta << ", a = " << v.radius << " fm\n"
          << "# Columns: E (MeV), q (MeV/c), p a, l, G11, G22, Re G12, Im G12 (spectatorAngle direction).\n";
      for (double e : {eLo, 0.5 * (eLo + eHi), eHi}) {
        std::vector<double> w, q, g, gd;
        double qd, pd;
        v.Interpolate(e, w, q, g, gd, qd, pd);
        for (size_t li = 0; li < v.lvals.size(); li++)
          out << "dw_vertex_point" << std::setw(18) << e << std::setw(18) << qd * hbarc << std::setw(18)
              << pd * v.radius << std::setw(4) << v.lvals[li] << std::setw(18) << gd[li * 4] << std::setw(18)
              << gd[li * 4 + 1] << std::setw(18) << gd[li * 4 + 2] << std::setw(18) << gd[li * 4 + 3] << "\n";
      }
    }
    // Coulomb line shape (lineshape=on): the ranges over the experiment's points.
    for (const ThmGroup &group : thmGroups_) {
      if (group.name != r.name || !group.lineshape) continue;
      const ThmLineshape &ls = *group.lineshape;
      double eLo = 1.0e300, eHi = -1.0e300;
      for (int s : group.segments)
        for (int p = 1; p <= GetSegment(s)->NumPoints(); p++) {
          eLo = std::min(eLo, GetSegment(s)->GetPoint(p)->GetCMEnergy());
          eHi = std::max(eHi, GetSegment(s)->GetPoint(p)->GetCMEnergy());
        }
      out << "lineshape: on (spectator " << ls.spectator << ", Z_s = " << ls.Zs << ", Z_F = " << ls.ZF
          << "; E_aA = " << ls.eAA << " MeV, B = " << ls.bind << " MeV)\n"
          << "# |N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)] per level, zeta = eta_sB - eta_0\n"
          << "# (Mukhamedzhanov et al. EPJA 56 (2020) 233 eqs. 56-62, case 2: m_B >> m_s, m_b, eta_sb\n"
          << "# neglected); zeta < 0 moves the peaks up in E.  eta_sb is the neglected s-b term, averaged\n"
          << "# over the b direction: the approximation assumes |eta_sb| << 1.  Rows: the value at the\n"
          << "# lowest and at the highest point energy E (MeV); [k] is the exit pair key.\n";
      auto row = [&](const char *name, double a, double b) {
        out << std::left << std::setw(16) << name << std::right << std::setw(18) << a << std::setw(18) << b << "\n";
      };
      row("E", eLo, eHi);
      row("E_sF", ls.EsF(eLo), ls.EsF(eHi));
      row("eta_0", ls.Eta0(eLo), ls.Eta0(eHi));
      for (const ThmLineshape::Exit &e : ls.exits) {
        std::ostringstream key;
        key << "zeta[" << e.pairKey << "]";
        row(key.str().c_str(), ls.Zeta(eLo, e.ZB, e.mB), ls.Zeta(eHi, e.ZB, e.mB));
        key.str("");
        key << "eta_sb[" << e.pairKey << "]";
        row(key.str().c_str(), ls.EtaSbEstimate(eLo, eLo + e.q, e.Zb, e.mb, e.mB),
            ls.EtaSbEstimate(eHi, eHi + e.q, e.Zb, e.mb, e.mB));
      }
    }
  }
}

bool EData::ThmLineshapeTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                              const Config &configure, ThmLineshapeReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->lineshape) {
    why = "THM experiment '" + name + "' has no line shape (lineshape=on)";
    return false;
  }
  const ThmLineshape &ls = *group->lineshape;
  out = ThmLineshapeReport();
  out.experiment = name;
  out.spectator = ls.spectator;
  out.Zs = ls.Zs;
  out.ZF = ls.ZF;
  out.eAA = ls.eAA;
  out.bind = ls.bind;
  out.energy = energies;
  for (double e : energies) {
    out.esf.push_back(ls.EsF(e));
    out.eta0.push_back(ls.Eta0(e));
  }
  int entranceKey = GetSegment(group->segments[0])->GetEntranceKey();
  int aa = compound->GetPairNumFromKey(entranceKey);
  PPair *entrance = compound->GetPair(aa);
  double threshold = entrance->GetSepE() + entrance->GetExE();
  for (const ThmLineshape::Exit &x : ls.exits) {
    ThmLineshapeReport::Exit ex;
    ex.pairKey = x.pairKey;
    ex.Zb = x.Zb;
    ex.ZB = x.ZB;
    ex.mb = x.mb;
    ex.mB = x.mB;
    for (double e : energies) {
      ex.zeta.push_back(ls.Zeta(e, x.ZB, x.mB));
      ex.etaSb.push_back(ls.EtaSbEstimate(e, e + x.q, x.Zb, x.mb, x.mB));
    }
    int exitNum = compound->GetPairNumFromKey(x.pairKey);
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      JGroup *jg = compound->GetJGroup(j);
      if (!jg->IsInRMatrix()) continue;
      bool in = false, outCh = false;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        in = in || jg->GetChannel(ch)->GetPairNum() == aa;
        outCh = outCh || jg->GetChannel(ch)->GetPairNum() == exitNum;
      }
      if (!in || !outCh) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmLineshapeReport::Level lv;
        lv.jgroup = j;
        lv.level = la;
        lv.J = jg->GetJ();
        lv.pi = jg->GetPi();
        lv.energy = level->GetFitE() - threshold;
        lv.width = ThmLevelWidth(compound, jg, level, configure);
        for (size_t k = 0; k < energies.size(); k++)
          lv.nc2.push_back(ThmLineshapeFactorSq(ex.zeta[k], lv.energy - energies[k], lv.width));
        ex.levels.push_back(lv);
      }
    }
    out.exits.push_back(ex);
  }
  return true;
}

bool EData::ThmDistortionTable(const std::string &name, const std::vector<double> &energies, ThmDistortionReport &out,
                               std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (!group->distortion) {
    why = "THM experiment '" + name + "' has no distortion (distortion=coulomb|optical|table:<file>)";
    return false;
  }
  const ThmDistortion &d = *group->distortion;
  out = ThmDistortionReport();
  out.experiment = name;
  out.description = d.description;
  out.energy = energies;
  if (d.kind == ThmDistortion::TABLE) {
    out.kind = "table";
    for (double e : energies) out.rModel.push_back(d.Weight(e));
    return true;
  }
  out.kind = d.kind == ThmDistortion::COULOMB ? "coulomb" : "optical";
  out.eRef = d.eRef;
  out.eAA = d.eAA;
  out.bind = d.kin.bind;
  out.kAA = d.aa.k;
  out.etaAA = d.aa.eta;
  out.kappa = d.kappa;
  out.etaB = d.etaB;
  out.beta = d.beta;
  for (double e : energies) {
    ThmDistortion::Point p;
    std::string bad = d.CheckEnergy(e);
    if (bad.empty()) p = d.Evaluate(e);
    if (!bad.empty() || !p.ok) {
      why = "THM experiment '" + name + "': " + (bad.empty() ? p.why : bad);
      return false;
    }
    out.esf.push_back(p.esf);
    out.ksf.push_back(p.ksf);
    out.etasf.push_back(p.etasf);
    out.thetaCm.push_back(p.thetaCm);
    out.x.push_back(p.x);
    out.q.push_back(p.q);
    out.m2.push_back(std::norm(p.m));
    out.mpw2.push_back(p.mpw * p.mpw);
    out.r.push_back(d.R(p));
    out.rModel.push_back(d.Weight(e));
    out.lmax.push_back(p.lmax);
  }
  return true;
}

bool EData::ThmVertexTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                           const Config &configure, ThmVertexReport &out, std::string &why) {
  const ThmGroup *group = nullptr;
  for (const ThmGroup &g : thmGroups_)
    if (g.name == name) group = &g;
  if (!group) {
    why = "no THM experiment '" + name + "' in use";
    return false;
  }
  if (group->pairKey == 0 || !compound->IsPairKey(group->pairKey)) {
    why = "THM experiment '" + name + "': its segments have different entrance pairs";
    return false;
  }
  const bool useGSL = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int aa = compound->GetPairNumFromKey(group->pairKey);
  PPair *pair = compound->GetPair(aa);
  if (pair->GetPType() != 0) {
    why = "THM experiment '" + name + "': the entrance pair is not a particle pair";
    return false;
  }
  out = ThmVertexReport();
  out.experiment = name;
  out.pairKey = group->pairKey;
  out.bind = pair->GetBindingEnergy();
  out.radius = pair->GetChRad();
  if (group->window) {
    const ThmSpectatorWindow &w = *group->window;
    out.window = w.description;
    out.muSx = w.muSx;
    out.p = w.p;
    out.weight = w.weight;
    out.es = w.es;
  } else {
    out.window = "delta";
    out.p.push_back(0.0);
    out.weight.push_back(1.0);
    out.es.push_back(configure.thm.SpectatorEnergy(group->pairKey));
  }
  out.energy = energies;
  const ThmDwVertex *dw = group->dwVertex.get();
  const double mu = pair->GetRedMass() * uconv;
  // M_l = (B - 1) j_l - rho j_l' + C_l at E with T_s = es added to E + B (EPoint::CalcEDependentValues).
  struct Pieces {
    double jl = 0.0, rhoDjl = 0.0;
    complex coul = complex(0.0, 0.0);
  };
  auto pieces = [&](int l, double e, double es) {
    Pieces q;
    double b = out.bind + es;
    if (e + b > 0.0) {
      ThmBesselParts(l, mu, e, b, out.radius, q.jl, q.rhoDjl);
      if (configure.thm.coulombIntegral && pair->GetZ(1) * pair->GetZ(2) != 0)
        q.coul = ThmCoulombTerm(pair, l, e, ThmRho(mu, e, b, 1.0), useGSL);
    }
    return q;
  };
  // vertexModel=dw: the Gram matrices of the DW vertex on the grid.
  std::vector<std::vector<double>> dwG, dwGd;
  if (dw) {
    out.model = "dw";
    for (double e : energies) {
      std::vector<double> w, q, g, gd;
      double qd, pd;
      dw->Interpolate(e, w, q, g, gd, qd, pd);
      const double ks = std::sqrt(2.0 * dw->dist.sf.mu * std::max(dw->dist.EsF(e), 0.0)) / hbarc;
      const double ka = dw->dist.aa.k, kb = dw->beta * ka;
      std::vector<double> row;
      for (double qk : q) {
        double qq = qk / hbarc, x = ks > 0.0 ? (ks * ks + kb * kb - qq * qq) / (2.0 * ks * kb) : 1.0;
        x = std::max(-1.0, std::min(1.0, x));
        row.push_back(std::sqrt(std::max(0.0, ka * ka + dw->alpha * dw->alpha * ks * ks - 2.0 * ka * dw->alpha * ks * x)) *
                      out.radius);
      }
      out.rho.push_back(row);
      out.dwQ.push_back(q);
      out.dwWeight.push_back(w);
      out.dwQDelta.push_back(qd * hbarc);
      out.dwPDelta.push_back(pd);
      dwG.push_back(g);
      dwGd.push_back(gd);
    }
  } else
    for (double e : energies) {
      std::vector<double> row;
      for (double es : out.es) row.push_back(e + out.bind + es > 0.0 ? ThmRho(mu, e, out.bind + es, out.radius) : 0.0);
      out.rho.push_back(row);
    }
  // Vertex boundary, as THMMatrixFunc::CalculateTHMCrossSection chooses it.
  const bool perLevel = (configure.paramMask & Config::USE_BRUNE_FORMALISM) &&
                        configure.thm.vertex == Config::ThmOptions::PER_LEVEL;
  const bool onShell = configure.thm.vertex == Config::ThmOptions::ON_SHELL;
  const bool constantVertex = configure.thm.vertex == Config::ThmOptions::CONSTANT;
  const double threshold = pair->GetSepE() + pair->GetExE();
  ChannelFunc channelFunc(pair, useGSL);
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa) continue;
      ThmVertexReport::Channel cr;
      cr.jgroup = j;
      cr.channel = ch;
      cr.J = jg->GetJ();
      cr.pi = jg->GetPi();
      cr.l = c->GetL();
      cr.s = c->GetS();
      double eMin = 0.0;
      bool found = false;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        if (!found || level->GetFitE() < eMin) eMin = level->GetFitE();
        found = true;
      }
      // The pieces on the grid: [E][node] and the quasi-free [E].
      std::vector<std::vector<Pieces>> grid(energies.size());
      std::vector<Pieces> qf(energies.size());
      if (!dw)
        for (size_t i = 0; i < energies.size(); i++) {
          for (double es : out.es) grid[i].push_back(pieces(cr.l, energies[i], es));
          qf[i] = pieces(cr.l, energies[i], 0.0);
        }
      const int li = dw ? dw->LIndex(cr.l) : -1;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        ThmVertexReport::Level lr;
        lr.level = la;
        double fixedB = constantVertex ? channelFunc.Shift(cr.l, eMin - threshold)
                        : perLevel     ? channelFunc.Shift(cr.l, level->GetFitE() - threshold)
                                       : c->GetBoundaryCondition();
        lr.boundary = onShell ? std::nan("") : fixedB;
        for (size_t i = 0; i < energies.size(); i++) {
          complex B(fixedB, 0.0);
          if (onShell) B = complex(channelFunc.Shift(cr.l, energies[i]), channelFunc.Penetrability(cr.l, energies[i]));
          if (dw) {
            const int nl = (int)dw->lvals.size();
            double avg = 0.0;
            for (size_t k = 0; k < out.dwWeight[i].size(); k++)
              avg += out.dwWeight[i][k] * (li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwG[i][(k * nl + li) * 4], B));
            lr.m2.push_back(avg);
            lr.m2qf.push_back(li < 0 ? 0.0 : ThmDwVertex::Vertex2Factored(&dwGd[i][li * 4], B));
            // Plane waves at the same kinematics: M_l(p), analytic j_l'.
            double rho = out.dwPDelta[i] * out.radius;
            double jl = gsl_sf_bessel_jl(cr.l, rho);
            double djl = rho > 0.0 ? cr.l * jl / rho - gsl_sf_bessel_jl(cr.l + 1, rho) : (cr.l == 1 ? 1.0 / 3.0 : 0.0);
            lr.m2pw.push_back(std::norm((B - 1.0) * jl - rho * djl));
            continue;
          }
          auto m2 = [&](const Pieces &q) { return std::norm((B - 1.0) * q.jl - q.rhoDjl + q.coul); };
          double avg = 0.0;
          for (size_t k = 0; k < out.es.size(); k++) avg += out.weight[k] * m2(grid[i][k]);
          lr.m2.push_back(avg);
          lr.m2qf.push_back(m2(qf[i]));
        }
        cr.levels.push_back(lr);
      }
      out.channels.push_back(cr);
    }
  }
  return true;
}
