#include "AZUREOutput.h"
#include "CNuc.h"
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
  energyShiftSqrtParamOffset_ = 0;
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
      if (segment.isDiff() == 8 || segment.isDiff() == 18) {
        // Polarization x Cross Section was removed: the outgoing polarization of
        // A(a,b)B is the analyzing power of the inverse reaction B(b,a)A.  18 is
        // the same observable on a THM segment (isDiff offset 10).
        configure.outStream
            << "ERROR: Data segment #" << numTotalSegments << " uses isDiff " << segment.isDiff() << " (Polarization x"
            << " Cross Section), which is no longer supported." << std::endl
            << "       Divide the data by the differential cross section and enter the"
            << " polarization as an Analyzing Power (isDiff 7) segment on the inverse"
            << " channel (entrance = the measured exit pair, exit = the measured entrance pair)."
            << std::endl;
        return -1;
      }
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

  if (ReadPriorCentres(configure) != 0) return -1;

  if (CheckThmWeights(configure, configure.thm.weightBySegment, "weight", "<segmentsData>",
                      numTotalSegments, this) != 0)
    return -1;
  // THM (HOES) segments are built and validated with Brune's level matrix;
  // Park's parametrization (--use-park) has not been checked for them.
  if (configure.paramMask & Config::USE_PARK_FORMALISM)
    for (int s = 1; s <= this->NumSegments(); s++)
      if (this->GetSegment(s)->IsTHM()) {
        configure.outStream << "ERROR: Data segment " << this->GetSegment(s)->GetSegmentKey()
                            << " is a THM segment (isDiff >= 10); THM is not available with --use-park."
                            << std::endl;
        return -1;
      }
  if (BuildThmGroups(configure, theCNuc, numTotalSegments) != 0) return -1;

  if (this->NumSegments() > 0) {
    if (this->ReadTargetEffectsFile(configure, theCNuc) == -1) return -1;
    this->MapData();
  }

  return 0;
}

/*!
 * Reads the explicit prior centres of <parameterSettings> (see EData.h) and
 * moves the nominal normalization / energy shift of the segments they name.
 * Only rows of exactly three tokens whose second token is "prior_centre" are
 * read; every other row is left to ParameterLimitsManager, as before.  A row
 * for a segment key with no active segment is reported and ignored (an
 * inactive line keeps its key); a malformed row, or a normalization centre
 * that is not positive, fails the run.
 */

int EData::ReadPriorCentres(const Config &configure) {
  std::ifstream in(configure.configfile.c_str());
  if (!in) return 0;
  std::string line;
  bool inside = false;
  while (getline(in, line)) {
    const size_t a = line.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) continue;
    const size_t b = line.find_last_not_of(" \t\r\n");
    const std::string trimmed = line.substr(a, b - a + 1);
    if (trimmed == "<parameterSettings>") {
      inside = true;
      continue;
    }
    if (trimmed[0] == '<') {
      inside = false;
      continue;
    }
    if (!inside || trimmed[0] == '#') continue;
    std::istringstream stm(trimmed);
    std::vector<std::string> tok;
    std::string t;
    while (stm >> t) tok.push_back(t);
    if (tok.size() != 3 || tok[1] != "prior_centre") continue;

    const std::string &name = tok[0];
    bool isNorm = false;
    std::string number;
    const std::string normTail = "_norm", shiftTail = "_energy_shift";
    if (name.compare(0, 8, "segment_") == 0 && name.size() > 8 + normTail.size() &&
        name.compare(name.size() - normTail.size(), normTail.size(), normTail) == 0) {
      isNorm = true;
      number = name.substr(8, name.size() - 8 - normTail.size());
    } else if (name.compare(0, 8, "segment_") == 0 && name.size() > 8 + shiftTail.size() &&
               name.compare(name.size() - shiftTail.size(), shiftTail.size(), shiftTail) == 0) {
      number = name.substr(8, name.size() - 8 - shiftTail.size());
    }
    const bool digits = !number.empty() && number.find_first_not_of("0123456789") == std::string::npos;
    char *end = NULL;
    const double centre = strtod(tok[2].c_str(), &end);
    if (!digits || end == tok[2].c_str() || *end != '\0' || !std::isfinite(centre)) {
      configure.outStream << "ERROR: <parameterSettings> row \"" << trimmed
                          << "\" is not \"segment_N_norm prior_centre value\" or "
                             "\"segment_N_energy_shift prior_centre value\"." << std::endl;
      return -1;
    }
    if (isNorm && centre <= 0.) {
      configure.outStream << "ERROR: <parameterSettings> row \"" << trimmed
                          << "\": a normalization prior centre must be positive." << std::endl;
      return -1;
    }
    const int key = atoi(number.c_str());
    bool found = false;
    for (int s = 1; s <= this->NumSegments(); s++) {
      ESegment *segment = this->GetSegment(s);
      if (segment->GetSegmentKey() != key) continue;
      found = true;
      if (isNorm)
        segment->SetNominalNorm(centre);
      else
        segment->SetNominalEnergyShift(centre);
    }
    if (!found)
      configure.outStream << "WARNING: <parameterSettings> row \"" << trimmed
                          << "\" names no active data segment; ignored." << std::endl;
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
      if (segment.isDiff() == 8 || segment.isDiff() == 18) {
        configure.outStream
            << "ERROR: Test segment #" << numTotalSegments << " uses isDiff " << segment.isDiff() << " (Polarization x"
            << " Cross Section), which is no longer supported; use an Analyzing Power (isDiff 7)"
            << " segment on the inverse channel instead." << std::endl;
        return -1;
      }
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
          if (isValidTotal || theCNuc->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            if (NewSegment.IsTHM() && (configure.paramMask & Config::USE_PARK_FORMALISM)) {
              configure.outStream << "ERROR: Test segment " << numTotalSegments
                                  << " is a THM segment (isDiff >= 10); THM is not available with --use-park."
                                  << std::endl;
              return -1;
            }
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
 * Returns the index at which the sqrt(E) energy-shift coefficients start in
 * the Minuit parameter vector (one per segment, after the energy shifts).
 */

int EData::GetEnergyShiftSqrtParamOffset() const {
  return energyShiftSqrtParamOffset_;
}

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
            if (segment && segment->IsTHM() && targetEffect.IsTargetIntegration()) {
              // A THM energy is the x + A energy reconstructed from the
              // three-body kinematics, not a beam energy slowing down in the
              // target; its sub-points also go below the entrance threshold,
              // where a stopping power is evaluated at E < 0 (NaN).
              configure.outStream << "ERROR: <targetInt>: segment " << segmentsList[i - 1]
                                  << " is a THM segment: target integration does not apply to it "
                                     "(use the Gaussian convolution for the THM energy resolution)."
                                  << std::endl;
              return -1;
            }
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
      targetEffect->ConvertSigmaToCM(cmConversion);
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
                backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure) * 5.0;
                forwardDepth = targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure) * 5.0;
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
              backwardDepth = convRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure);
              forwardDepth = convRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure);
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
          targetEffect->ConvertSigmaToCM(cmConversion);
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
                    backwardDepth = targetThickness + targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure) * 5.0;
                    forwardDepth = targetEffect->convolutionRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure) * 5.0;
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
                  backwardDepth = convRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure);
                  forwardDepth = convRange * targetEffect->CalculateSigma(point->GetCMEnergy(), configure);
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
           << totalN;
    // Park formalism: the J > 0 wall, when a level sits beyond it.
    if (configure.paramMask & Config::USE_PARK_FORMALISM)
      chiOut << " Total-Park-Chi-Squared: " << parkPenalty_;
    chiOut << std::endl
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

  // Check if any energy shifts (constant or sqrt(E) term) are varied
  bool isVaryEnergyShift = false;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt()) {
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
      // One complete row per segment with a free shift term; sqrt_shift is the
      // coefficient b of E' = E + shift + b*sqrt(E/MeV), in MeV^1/2.
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift, sqrt_shift" << std::endl;
      for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
        if (!(segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt())) continue;
        out << segment->GetSegmentKey() << ","
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
            << segment->GetEnergyShift() << ","
            << segment->GetEnergyShiftSqrt() << std::endl;
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

void EData::SetEnergyShiftSqrtParamOffset(int offset) {
  energyShiftSqrtParamOffset_ = offset;
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

  // sqrt(E) energy-shift coefficients, again one per segment so the layout
  // (and every parameter index downstream) depends only on the segment count.
  // A segment without a "sqrtshift" block carries a fixed 0 here, so a model
  // written before the term existed fits exactly as before.
  SetEnergyShiftSqrtParamOffset(p.Params().size());
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    snprintf(varname, sizeof(varname), "segment_%d_energy_shift_sqrt", segment->GetSegmentKey());
    double stepSize = (segment->GetEnergyShiftSqrtError() > 0.0) ? segment->GetEnergyShiftSqrtError() * 0.01 : 0.0005;
    p.Add(varname, segment->GetEnergyShiftSqrt(), stepSize);
    if (!segment->IsVaryEnergyShiftSqrt()) p.Fix(varname);
  }

  // THM coherent backgrounds (cbackground=): Re/Im of each amplitude, last.
  thmCoherentParamOffset_ = thmCoherentParams_.empty() ? -1 : (int)p.Params().size();
  for (const ThmCoherentParam &c : thmCoherentParams_) {
    p.Add(c.name, c.value, c.value != 0.0 ? 0.1 * std::fabs(c.value) : 0.1);
    if (c.fixed) p.Fix(c.name);
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
  // The THM coherent backgrounds, the block after the shifts (none without
  // cbackground=): every caller fills the shifts right after the norms, with
  // the compound the model will run on.
  FillThmCoherentFromParams(p, theCNuc);
  int i = GetEnergyShiftParamOffset();
  int j = GetEnergyShiftSqrtParamOffset();
  int k = 0;
  bool anyEnergyShifted = false;

  if (data) {
    for (ESegmentIterator segment = data->GetSegments().begin(); segment < data->GetSegments().end(); segment++) {
      k++;

      // Both terms are present for every segment: the constant shift at p[i]
      // and the sqrt(E) coefficient at p[j].  Either moving the points means
      // the segment must be re-shifted (the mapping E -> E + a + b*sqrt(E) is
      // applied from the original energies, so one update covers both).
      const double shiftConst = p[i];
      const double shiftSqrt = (j >= 0 && j < (int)p.size()) ? p[j] : 0.0;
      const bool shiftActive = segment->IsVaryEnergyShift() || shiftConst != 0.0 ||
                               segment->IsVaryEnergyShiftSqrt() || shiftSqrt != 0.0;
      if (shiftActive) {
        // Check if energy is the same, if so, continue
        if (segment->GetLastEnergyShift() == shiftConst && segment->GetLastEnergyShiftSqrt() == shiftSqrt) {
          i++;
          j++;
          continue;
        }

        segment->SetEnergyShift(shiftConst);
        segment->SetLastEnergyShift(shiftConst);
        segment->SetEnergyShiftSqrt(shiftSqrt);
        segment->SetLastEnergyShiftSqrt(shiftSqrt);
        segment->UpdatePointEnergiesWithShift(theCNuc, configure);
        anyEnergyShifted = true;

        // Apply the same energy shift to ALL component segments of this master segment
        if (segment->HasComponents()) {
          for (ESegment *componentSegment : segment->GetComponentSegments()) {
            if (componentSegment) {
              componentSegment->SetEnergyShift(segment->GetEnergyShift());
              componentSegment->SetEnergyShiftSqrt(segment->GetEnergyShiftSqrt());
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
            }
          }
        }
      }

      // Apply the same energy shift to all component segments in total capture segments
      if (segment->IsTotalCapture() && segment->HasComponents()) {
        if (shiftActive) {
          // Apply energy shift to all component segments
          const std::vector<ESegment *> &componentSegments = segment->GetComponentSegments();
          for (ESegment *componentSegment : componentSegments) {
            if (componentSegment) {
              // Check if energy is the same, if so, skip the recompute.
              // NOTE: do NOT advance i or j here. All components of a total-capture
              // segment share this segment's single pair of shift parameters;
              // i and j are advanced exactly once per segment at the end of the
              // loop.  Advancing them inside the component loop desyncs the
              // parameter-to-segment mapping for every subsequent segment.
              if (componentSegment->GetLastEnergyShift() == shiftConst &&
                  componentSegment->GetLastEnergyShiftSqrt() == shiftSqrt) {
                continue;
              }

              componentSegment->SetEnergyShift(shiftConst);
              componentSegment->SetLastEnergyShift(shiftConst);
              componentSegment->SetEnergyShiftSqrt(shiftSqrt);
              componentSegment->SetLastEnergyShiftSqrt(shiftSqrt);
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
              anyEnergyShifted = true;
            }
          }
        }
      }

      // Always increment i and j since both shift parameters are present for ALL segments
      i++;
      j++;
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
  dataCopy->energyShiftSqrtParamOffset_ = this->energyShiftSqrtParamOffset_;
  dataCopy->isFit_ = this->isFit_;
  dataCopy->isErrorAnalysis_ = this->isErrorAnalysis_;
  dataCopy->targetEffects_ = this->targetEffects_;
  dataCopy->segments_ = this->segments_;
  dataCopy->componentSegments_ = this->componentSegments_;
  dataCopy->thmGroups_ = this->thmGroups_;
  dataCopy->thmCoherentParams_ = this->thmCoherentParams_;
  dataCopy->thmCoherentParamOffset_ = this->thmCoherentParamOffset_;

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
