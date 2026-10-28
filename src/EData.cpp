#include "AZUREOutput.h"
#include "CNuc.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "EData.h"
#include "ECAmplitudeCache.h"
#include "ECTable.h"
#include "ExtrapLine.h"
#include "SegLine.h"
#include "AdaptiveIntegrationGrid.h"
#include "Minuit2/MnUserParameters.h"
#include "GSLException.h"
#include <gsl/gsl_linalg.h>
#include <map>
#include <tuple>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cmath>
#include <limits>
#include "NuclearPotentialManager.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <unordered_map>

/*!
 * A sample segment ("sample <1|name>", alias "natural") adds abundance-weighted terms of
 * several compound nuclei, so it is accepted where a SUM of cross sections is
 * meaningful -- angle-integrated, differential and capture cross sections,
 * plain or SUM composites -- and refused for a RATIO, a phase shift, an
 * analyzing power and a test-mode angular-distribution segment, which are not
 * additive.  It also needs the <abundances> block: without it every weight is
 * 1 and the token would be silently meaningless, which is exactly what an
 * older build would do with it too.
 */

/*!
 * An analyzing-power segment (isDiff 7) needs a spin-1/2 projectile: the
 * particle route (PolarizationFunc) is written for one and returned 0 for any
 * other, and the capture route (CNuc::CalcCaptureAnalyzingPower) ran for any
 * projectile spin but returned 3/(x+1) <S_y>, which is the Madison A_y only
 * at x = 1/2 (for a deuteron 1.5 times Ohlsen's Cartesian A_y).  Both are
 * refused here instead (2026-10-26, audit item 9); tensor observables and the
 * spin-1 vector analyzing power are not implemented.
 */
/*!
 * An angular-coefficient segment (isDiff 9) is the Legendre expansion of a reaction's c.m.
 * angular distribution.  Refused: an elastic pair (the expansion of GenMatrixFunc leaves the
 * Coulomb amplitude out), total capture, composites, samples, "total", UPOS, and a free
 * normalization on normalised coefficients (a_k / a_0 has no scale; use "absolute 1" for A_k).
 */
static bool CoefficientsAllowed(CNuc *theCNuc, const ESegment &segment, int number, const Config &configure) {
  if (!segment.IsCoefficients()) return true;
  std::string why;
  if (segment.GetExitKey() == -1)
    why = "total capture has no angular distribution";
  else if (segment.GetExitKey() == segment.GetEntranceKey())
    why = "an elastic distribution includes the Coulomb amplitude, which the Legendre expansion leaves out";
  else if (segment.HasComponents() || segment.IsAdvanced())
    why = "a composite (SUM / RATIO) is not supported";
  else if (segment.IsNatural())
    why = "a sample segment is not supported";
  else if (segment.IsTotal())
    why = "\"total 1\" has no angular distribution";
  else if (segment.IsUPOS())
    why = "a secondary-gamma (UPOS) distribution is not supported";
  else if (!segment.IsCoefficientAbsolute() && segment.IsVaryNorm())
    why = "normalised coefficients a_k / a_0 have no scale: fix the normalization at 1, or give the absolute A_k with \"absolute 1\"";
  if (why.empty()) return true;
  configure.outStream << "ERROR: segment #" << number << " (angular coefficients): " << why << "." << std::endl;
  return false;
}

static bool AnalyzingPowerAllowed(CNuc *theCNuc, const ESegment &segment, int number, bool isTest,
                                  const Config &configure) {
  if (!segment.IsAnalyzingPower()) return true;
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetEntranceKey()));
  if (fabs(entrancePair->GetJ(1) - 0.5) < 1.e-6) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << number
                      << ": an analyzing power needs a spin-1/2 projectile; entrance pair " << segment.GetEntranceKey()
                      << " has a spin-" << entrancePair->GetJ(1) << " projectile (tensor observables and the spin-1 vector"
                      << " analyzing power are not implemented)." << std::endl;
  return false;
}

/*!
 * A "total 1" segment is the total cross section of its entrance pair by the optical theorem
 * (GenMatrixFunc::CalculateBareCrossSection): angle-integrated, its exit pair its entrance pair, and a neutral pair
 * (the Coulomb total diverges) of distinguishable particles.  The R-matrix-for-capture formalism defines T
 * differently, so it is refused there.
 */
/*!
 * A UPOS (secondary-gamma) segment given as a lab differential is converted with the kinematics of particle b, the
 * light ejectile of the exit pair, which is right for a gamma measured in coincidence with it and not otherwise
 * (segments.rst, UPOS).  Say so instead of converting silently (audit 2026-10-25, item 31).
 */
static void UPOSLabNote(const ESegment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !segment.IsDifferential() || segment.IsCMDifferential() || segment.IsPhase()) return;
  configure.outStream << "WARNING: " << (isTest ? "test segment #" : "segment #") << key
                      << " (UPOS, lab differential): the lab angle is converted with the kinematics of particle b, "
                      << "the light ejectile of the exit pair (not the gamma's); give c.m. angles if that is not meant."
                      << std::endl;
}

/*!
 * Total capture is the exit key -1 together with the total-capture data type (3 on a data line, 4 on a test
 * line): the reader expands such a segment into one term per capture pair.  Either half alone -- the key -1
 * with another type, or the type with a real exit pair -- left a segment whose exit pair does not exist:
 * mode 3 with "1 1 -1 ... 0" read pairs_[-1] and crashed in the external-capture lookup (13N), and a model
 * without external capture gave meaningless numbers (audit 2026-10-25, item 32).
 */
static bool TotalCaptureCodeOk(const ESegment &segment, int key, bool isTest, const Config &configure) {
  const bool keyTotal = segment.GetExitKey() == -1;
  const bool codeTotal = segment.IsTotalCapture() != 0;
  if (keyTotal == codeTotal) return true;
  const char *what = isTest ? "test segment #" : "segment #";
  if (keyTotal) {
    configure.outStream << "ERROR: " << what << key << ": the exit key -1 (total capture) needs the total-capture data type ("
                        << (isTest ? "4 on a test line" : "3 on a data line") << ")." << std::endl;
    return false;
  }
  // The type with a real exit key has always meant total capture: the expansion below replaces the exit key with
  // the first capture pair and adds the others as components.  Say so instead of doing it silently.
  configure.outStream << "WARNING: " << what << key << ": the total-capture data type with exit pair " << segment.GetExitKey()
                      << ": the exit key is ignored; the segment is the total capture over every capture pair of the "
                      << "entrance nucleus (write the exit key -1 to say so)." << std::endl;
  return true;
}

static bool TotalAllowed(CNuc *theCNuc, ESegment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsTotal()) return true;
  const char *what = isTest ? "test segment" : "segment";
  auto refuse = [&](const std::string &why) {
    configure.outStream << "ERROR: " << what << " #" << key << " (total 1): " << why << std::endl;
    return false;
  };
  if (segment.GetExitKey() != segment.GetEntranceKey())
    return refuse("the total cross section of a pair needs the segment's exit pair to be its entrance pair.");
  if (segment.IsDifferential() || segment.IsPhase() || segment.IsAnalyzingPower() || segment.IsAngularDist() ||
      segment.IsUPOS())
    return refuse("a total cross section is angle-integrated (data type 0).");
  if (!theCNuc->IsPairKey(segment.GetEntranceKey())) return true;  // reported where the pair is checked
  PPair *pair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetEntranceKey()));
  if (pair->GetPType() != 0 || pair->GetZ(1) * pair->GetZ(2) != 0)
    return refuse("the optical theorem needs a neutral particle pair (the total cross section of charged particles "
                  "diverges).");
  if (pair->IsIdentical()) return refuse("not for a pair of identical particles.");
  if (configure.paramMask & Config::USE_RMC_FORMALISM)
    return refuse("not with the R-matrix-for-capture formalism (--use-rmc).");
  return true;
}

static bool NaturalAllowed(CNuc *theCNuc, ESegment &segment, int key, bool isTest, const Config &configure) {
  const char *what = isTest ? "test segment" : "segment";
  const std::string &token = segment.GetNaturalToken();
  if (token == "0") return true;   // not a natural segment
  if (!theCNuc->HasAbundanceBlock()) {
    configure.outStream << "ERROR: " << what << " #" << key << " is a sample (abundance-weighted) segment, but the project has no "
                        << "<abundances> block (one line per compound nucleus: [<sample>:] <Z> <A> <abundance> <vary> <sigma>)." << std::endl;
    return false;
  }
  // "1" = the unnamed sample (or the only one); anything else names a sample of the block.
  int sample = (token == "1") ? theCNuc->DefaultSample() : theCNuc->FindSample(token);
  if (sample == 0) {
    configure.outStream << "ERROR: " << what << " #" << key << ": sample \"" << token
                        << "\" is not in the <abundances> block";
    if (token == "1") configure.outStream << " (it has several named samples and no unnamed one: name one on the segment line)";
    configure.outStream << "." << std::endl;
    return false;
  }
  segment.SetNaturalSample(sample);
  if (segment.IsAdvanced() && segment.GetLegacyOperationType() != 0) {
    configure.outStream << "ERROR: " << what << " #" << key << ": a sample segment cannot be a RATIO." << std::endl;
    return false;
  }
  if (segment.IsPhase() || segment.IsAnalyzingPower() || (isTest && segment.IsAngularDist())) {
    configure.outStream << "ERROR: " << what << " #" << key << ": a sample segment must be a cross section "
                        << "(angle-integrated, differential or capture), not a phase shift, analyzing power or "
                        << "angular-distribution coefficients." << std::endl;
    return false;
  }
  return true;
}

/*!
 * The decomposition of a sample segment, one file per segment since
 * AZUREOut files are named by pair keys only: the energy, the total the fit
 * compares with the data, and the abundance-weighted term of every compound
 * nucleus (ESegment::CalculateTheoreticalCrossSection).  Same frame as the
 * AZUREOut file (kinflag.dat), one line per point.
 */
static void WriteNaturalDecomposition(const Config &configure, bool isExtrap, int kinflag,
                                      const std::vector<std::string> &nucleusLabels, ESegment *segment) {
  std::string name = configure.outputdir + "AZUREOut_sample_seg=" + std::to_string(segment->GetSegmentKey()) +
                     (isExtrap ? ".extrap" : ".out");
  std::ofstream out(name.c_str());
  if (!out) return;
  out << "# Segment " << segment->GetSegmentKey() << ": abundance-weighted sum over compound nuclei (sample "
      << segment->GetNaturalToken() << ")" << std::endl;
  out << "# E_" << (kinflag ? "lab" : "cm") << " [MeV]   total";
  for (size_t k = 0; k < nucleusLabels.size(); k++) out << "   " << nucleusLabels[k];
  out << std::endl;
  out.precision(10);
  for (int i = 1; i <= segment->NumPoints(); i++) {
    EPoint *point = segment->GetPoint(i);
    double kin = (kinflag != 0) ? point->GetCrossSectionKinFactor() : 1.0;
    if (kin == 0.0) kin = 1.0;
    out << std::setw(18) << std::scientific << (kinflag != 0 ? point->GetLabEnergy() : point->GetCMEnergy())
        << std::setw(18) << std::scientific << point->GetFitCrossSection() / kin;
    const std::vector<double> &terms = point->GetNaturalTerms();
    for (size_t k = 0; k < nucleusLabels.size(); k++)
      out << std::setw(18) << std::scientific << (k < terms.size() ? terms[k] / kin : 0.0);
    out << std::endl;
  }
}


/*!
 * The EData object has a private attribute containing the number of iterations needed to find the best fit parameters.
 * At creation, this attribute is set to 0.
 */

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
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) {
        if (!segment.parseError().empty())
          configure.outStream << "ERROR: Data segment #" << numTotalSegments + 1 << ": " << segment.parseError() << "."
                              << std::endl;
        return -1;
      }
      numTotalSegments++;
      if (segment.isDiff() == 8) {
        // Polarization x Cross Section was removed: the outgoing polarization of
        // A(a,b)B is the analyzing power of the inverse reaction B(b,a)A.
        configure.outStream
            << "ERROR: Data segment #" << numTotalSegments << " uses isDiff 8 (Polarization x"
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
          PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          entrancePair->SetEntrance();
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!CoefficientsAllowed(theCNuc, NewSegment, numTotalSegments, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalCaptureCodeOk(NewSegment, numTotalSegments, false, configure)) return -1;
          UPOSLabNote(NewSegment, numTotalSegments, false, configure);
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            // total capture: every capture pair of the ENTRANCE pair's compound nucleus
            for (int i = 1; i <= theCNuc->NumPairs(); i++) {
              if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == entrancePair->GetNucleus()) {
                isValidTotal = true;
                break;
              }
            }
          }
          if (!isValidTotal && theCNuc->IsPairKey(NewSegment.GetExitKey()) &&
              theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus() != entrancePair->GetNucleus()) {
            configure.outStream << "ERROR: segment #" << numTotalSegments << ": entrance pair " << NewSegment.GetEntranceKey()
                                << " and exit pair " << NewSegment.GetExitKey() << " belong to different compound nuclei ("
                                << theCNuc->NucleusLabel(entrancePair->GetNucleus()) << " and "
                                << theCNuc->NucleusLabel(theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus())
                                << "); a reaction stays within one nucleus." << std::endl;
            return -1;
          }
          if (isValidTotal || theCNuc->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            this->AddSegment(NewSegment);

            // Fill the segment with data first
            if (this->GetSegment(this->NumSegments())->Fill(theCNuc, this, configure) == -1) {
              if (covarianceError_ || kinematicsError_) return -1;  // never skipped silently
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
                            if (!this->CheckComponentPairs(theCNuc, *filledSegment, entranceKey, exitKey, configure)) return -1;
                            componentSegment = this->CreateComponentSegment(*filledSegment, entranceKey, exitKey, fixedAngle);
                          } else {
                            // No angle specified, use standard method
                            int exitKey = std::stoi(component.substr(exitPos));
                            if (!this->CheckComponentPairs(theCNuc, *filledSegment, entranceKey, exitKey, configure)) return -1;
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

                const int entranceNucleus = theCNuc->GetPair(theCNuc->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == entranceNucleus) {
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

  if (this->ResolveShiftGroups(configure) == -1) return -1;
  this->ResolveNormGroups(configure);
  if (this->ResolveCovarianceGroups(configure) == -1) return -1;

  if (this->NumSegments() > 0) {
    if (this->ReadTargetEffectsFile(configure, theCNuc) == -1) return -1;
    this->MapData(theCNuc);
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
      if (segment.hasCovariance()) {
        configure.outStream << "ERROR: Test segment #" << numTotalSegments
                            << " carries a covariance block; a data covariance belongs on data segments only."
                            << std::endl;
        return -1;
      }
      if (segment.isDiff() == 9) {
        configure.outStream << "ERROR: Test segment #" << numTotalSegments
                            << " uses isDiff 9 (angular coefficients), a data-segment type; a test segment"
                            << " computes the coefficients with isDiff 3 and a maximum order." << std::endl;
        return -1;
      }
      if (segment.isDiff() == 8) {
        configure.outStream
            << "ERROR: Test segment #" << numTotalSegments << " uses isDiff 8 (Polarization x"
            << " Cross Section), which is no longer supported; use an Analyzing Power (isDiff 7)"
            << " segment on the inverse channel instead." << std::endl;
        return -1;
      }
      if (segment.isActive() == 1) {
        ESegment NewSegment(segment);
        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          PPair *segEntrance = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalCaptureCodeOk(NewSegment, numTotalSegments, true, configure)) return -1;
          UPOSLabNote(NewSegment, numTotalSegments, true, configure);
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            // total capture: every capture pair of the ENTRANCE pair's compound nucleus
            for (int i = 1; i <= theCNuc->NumPairs(); i++) {
              if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == segEntrance->GetNucleus()) {
                isValidTotal = true;
                break;
              }
            }
          }
          if (!isValidTotal && theCNuc->IsPairKey(NewSegment.GetExitKey()) &&
              theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus() != segEntrance->GetNucleus()) {
            configure.outStream << "ERROR: test segment #" << numTotalSegments << ": entrance pair " << NewSegment.GetEntranceKey()
                                << " and exit pair " << NewSegment.GetExitKey() << " belong to different compound nuclei ("
                                << theCNuc->NucleusLabel(segEntrance->GetNucleus()) << " and "
                                << theCNuc->NucleusLabel(theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus())
                                << "); a reaction stays within one nucleus." << std::endl;
            return -1;
          }
          if (isValidTotal || theCNuc->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            this->AddSegment(NewSegment);
            ESegment *theSegment = this->GetSegment(this->NumSegments());

            theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()))->SetEntrance();
            PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()));
            PPair *exitPair = theCNuc->GetPairFromKey(theSegment->GetExitKey());  // nullptr for total capture (-1)
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
                if (entrancePair->GetPType() == 20 && exitPair)
                  thePoint->ConvertDecayEnergy(exitPair);
                else if (!theSegment->IsCMDifferential())
                  thePoint->ConvertLabEnergy(entrancePair);
                else if (theSegment->IsCMDifferential())
                  thePoint->ConvertLabEnergy(entrancePair);
                // A RATIO is dimensionless: no cross-section conversion (as ESegment::Fill for data segments;
                // the stored factor would otherwise divide the ratio in the lab output, kinflag.dat = 1).
                const bool isRatio = theSegment->IsAdvanced() && theSegment->GetLegacyOperationType() == 1;
                if (exitPair && exitPair->GetPType() == 0 && theSegment->IsDifferential() &&
                    !theSegment->IsPhase() && !theSegment->IsAngularDist() && !theSegment->IsCMDifferential()) {
                  const bool inRange = (theSegment->GetEntranceKey() == theSegment->GetExitKey())
                                           ? thePoint->ConvertLabAngle(entrancePair)
                                           : thePoint->ConvertLabAngle(entrancePair, exitPair, configure);
                  if (!inRange) {
                    configure.outStream << "ERROR: " << "test segment #" << theSegment->GetSegmentKey() << ": lab angle " << thePoint->GetLabAngle() << " deg at E_lab = " << thePoint->GetLabEnergy()
    << " MeV: " << thePoint->GetLabAngleWhy() << "." << std::endl;
                    return -1;
                  }
                  if (!isRatio) thePoint->ConvertCrossSection(entrancePair, exitPair);
                }
                if (exitPair && exitPair->GetPType() == 10 && theSegment->IsDifferential() &&
                    !theSegment->IsPhase() && !theSegment->IsAngularDist() && !theSegment->IsCMDifferential()) {
                  thePoint->ConvertLabAngleGammas(entrancePair);
                  if (!isRatio) thePoint->ConvertCrossSectionGammas(entrancePair);
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
                          if (!this->CheckComponentPairs(theCNuc, *theSegment, entranceKey, exitKey, configure)) return -1;
                          componentSegment = this->CreateComponentSegment(*theSegment, entranceKey, exitKey, fixedAngle);
                        } else {
                          // No angle specified, use standard method
                          int exitKey = std::stoi(component.substr(exitPos));
                          if (!this->CheckComponentPairs(theCNuc, *theSegment, entranceKey, exitKey, configure)) return -1;
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

                const int entranceNucleus = theCNuc->GetPair(theCNuc->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == entranceNucleus) {
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

  if (this->NumSegments() > 0) {
    if (this->ReadTargetEffectsFile(configure, theCNuc) == -1) return -1;
    this->MapData(theCNuc);
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

/*!
 * Returns the index at which the sqrt(E) energy-shift coefficients start in
 * the Minuit parameter vector (one per segment, after the energy shifts).
 */

int EData::GetEnergyShiftSqrtParamOffset() const {
  return energyShiftSqrtParamOffset_;
}

int EData::GetEnergyShiftFnParamOffset() const {
  return energyShiftFnParamOffset_;
}

/*!
 * Reads the target effects input file and creates the TargetEffect objects
 * to be applied to the data.
 */

namespace {
EPoint::SubGridSpec MakeSubGridSpec(const AdaptiveIntegrationGrid::GridConfig &gc, double start, double end,
                                    const TargetEffect *effect, double conversion) {
  EPoint::SubGridSpec spec;
  spec.valid = true;
  spec.startEnergy = start;
  spec.endEnergy = end;
  spec.maxPoints = gc.maxPoints;
  spec.entranceKey = gc.entranceKey;
  spec.baseEnergyStep = gc.baseEnergyStep;
  spec.resonanceWidthMultiplier = gc.resonanceWidthMultiplier;
  spec.pointsPerWidth = gc.pointsPerWidth;
  spec.targetIntegration = effect->IsTargetIntegration();
  spec.stoppingConversion = conversion;
  return spec;
}
}  // namespace

EData::GridNucleus::GridNucleus() {}
EData::GridNucleus::GridNucleus(const GridNucleus &) {}
EData::GridNucleus &EData::GridNucleus::operator=(const GridNucleus &) { return *this; }
EData::GridNucleus::~GridNucleus() {}

namespace {
// Does this segment need external-capture amplitudes at its sub-points?  They are
// integrals computed once at initialization (CalculateECAmplitudes), so its grid
// cannot be re-placed.
bool SegmentNeedsEC(ESegment &segment, CNuc *theCNuc, const Config &configure) {
  if (!(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) return false;
  const int ir = theCNuc->GetPairNumFromKey(segment.GetExitKey());
  for (int j = 1; j <= theCNuc->NumJGroups(); j++)
    for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
      ALevel *level = theCNuc->GetJGroup(j)->GetLevel(la);
      if (level->IsECLevel() && level->GetECPairNum() == ir) return true;
    }
  return false;
}
}  // namespace

/*!
 * Re-places the adaptive sub-point grids (DeBoer 2026-10-19, Coulomb table stage 3).
 *
 * ReadTargetEffectsFile lays each grid out once, from the level energies and
 * widths in the input file.  A fit that moves a narrow level then integrates it
 * on a grid built for another position, and the objective depends on which file
 * the session started from (natural neon refit1: 25,773.96 from the start file,
 * 25,776.08 from the fitted one, same parameters).  Here the grid is laid out
 * again for the current parameters: the fit values are converted to observed
 * energies and widths, as a saved file would hold them (CNuc::TransformOut on a
 * private copy of the compound nucleus), and fed to the same generator with the
 * same window, so the grid of a parameter set is the one a session started from
 * its saved file would place.  Only points whose grid changed get new
 * sub-points; their Coulomb functions are table lookups.
 *
 * Energy shifts: the sub-points are laid out at the unshifted energies and moved
 * with the data by the shift step; the resonances are offset by the point's c.m.
 * shift so the dense part of the grid lands on each resonance after the move.
 *
 * A grid differs from the current one if a node moved by more than 1e-12
 * relative: the round trip from the file's observed values through the formal
 * parameters and back changes them in the last digits only.
 */
int EData::ReplaceSubPointGrids(CNuc *theCNuc, const Config &configure, const std::vector<ESegment *> &shiftActive,
                                std::vector<ESegment *> &reshift) {
  if (!configure.replaceSubPointGrids || !configure.useAdaptiveGrid || !theCNuc) return 0;
  const auto t0 = std::chrono::steady_clock::now();

  struct Item { ESegment *seg; int i; double offset; bool shifted; };
  std::vector<Item> work;
  static std::atomic<bool> ecNoted{false};
  auto collect = [&](ESegment *seg) {
    bool any = false;
    for (int i = 1; i <= seg->NumPoints() && !any; i++) any = seg->GetPoint(i)->GetSubGridSpec().valid;
    if (!any) return;
    if (SegmentNeedsEC(*seg, theCNuc, configure)) {
      if (!ecNoted.exchange(true))
        configure.outStream << "Note: the sub-point grids of segments with external capture are not re-placed "
                               "(their capture integrals are computed once)." << std::endl;
      return;
    }
    const bool shifted = std::find(shiftActive.begin(), shiftActive.end(), seg) != shiftActive.end();
    PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(seg->GetEntranceKey()));
    PPair *exitPair = theCNuc->GetPairFromKey(seg->GetExitKey());  // nullptr for total capture (-1)
    // The c.m./lab factor the shift step uses for sub-points (ESegment::UpdatePointEnergyWithShift).
    const double cmFactor = (entrancePair->GetPType() == 20 && exitPair)
        ? (exitPair->GetM(1) + exitPair->GetM(2)) / exitPair->GetM(2)
        : entrancePair->GetM(2) / (entrancePair->GetM(1) + entrancePair->GetM(2));
    for (int i = 1; i <= seg->NumPoints(); i++) {
      EPoint *point = seg->GetPoint(i);
      if (!point->GetSubGridSpec().valid) continue;
      // the parent's shift, which is what moves every sub-point (ESegment::UpdatePointEnergyWithShift)
      const double offset = shifted ? cmFactor * seg->TotalEnergyShift(point->GetOriginalEnergy()) : 0.0;
      work.push_back({seg, i, offset, shifted});
    }
  };
  for (ESegment &seg : segments_) collect(&seg);
  for (ESegment &seg : componentSegments_) collect(&seg);
  if (work.empty()) return 0;

  // The observed parameters of the current fit values, on the private copy.
  if (!gridNuc_.nuc) gridNuc_.nuc.reset(theCNuc->Clone());
  CNuc *g = gridNuc_.nuc.get();
  for (int j = 1; j <= theCNuc->NumJGroups(); j++)
    for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
      ALevel *from = theCNuc->GetJGroup(j)->GetLevel(la);
      ALevel *to = g->GetJGroup(j)->GetLevel(la);
      to->SetFitE(from->GetFitE());
      for (int ch = 1; ch <= theCNuc->GetJGroup(j)->NumChannels(); ch++) to->SetFitGamma(ch, from->GetFitGamma(ch));
    }
  g->TransformOut(configure);
  g->FillCompoundFromParamsPhysical(g->GetTransformParams(configure));
  const auto t1 = std::chrono::steady_clock::now();

  std::atomic<int> changed{0};
  std::atomic<long> newSubPoints{0};
  // For AZURE_GRID_TIMING: per work item, -1 = a different number of nodes, else the
  // largest node move (MeV) of a grid that changed.
  std::vector<double> moved(work.size(), 0.0);
  auto replace = [&](const Item &it, size_t q) {
    EPoint *point = it.seg->GetPoint(it.i);
    const EPoint::SubGridSpec &spec = point->GetSubGridSpec();
    AdaptiveIntegrationGrid::GridConfig gc;
    gc.maxPoints = spec.maxPoints;
    gc.entranceKey = spec.entranceKey;
    gc.inputWidthsArePhysical = true;
    gc.baseEnergyStep = spec.baseEnergyStep;
    gc.resonanceWidthMultiplier = spec.resonanceWidthMultiplier;
    gc.pointsPerWidth = spec.pointsPerWidth;
    gc.resonanceOffset = it.offset;
    std::vector<double> grid = AdaptiveIntegrationGrid(gc).GenerateGrid(spec.startEnergy, spec.endEnergy, g);
    std::vector<EPoint> &subs = point->GetSubPoints();
    bool same = (subs.size() == grid.size());
    for (size_t k = 0; same && k < grid.size(); k++)
      same = std::fabs(subs[k].GetOriginalEnergy() - grid[k]) <= 1.0e-12 * std::max(1.0, std::fabs(grid[k]));
    if (same) return;
    if (subs.size() != grid.size())
      moved[q] = -1.0;
    else
      for (size_t k = 0; k < grid.size(); k++) moved[q] = std::max(moved[q], std::fabs(subs[k].GetOriginalEnergy() - grid[k]));
    TargetEffect *effect = this->GetTargetEffect(it.seg->GetTargetEffectNum());
    TargetEffect *qEffect = effect->IsQCoefficients() ? effect : nullptr;
    subs.clear();
    for (double e : grid) {
      EPoint sub(point->GetCMAngle(), e, it.seg);
      sub.SetCoefficientOrder(point->GetCoefficientOrder());
      sub.SetCoefficientAbsolute(point->GetCoefficientAbsolute());
      if (spec.targetIntegration)
        sub.SetStoppingPower(spec.stoppingConversion *
                             effect->GetStoppingPowerEq()->Evaluate(configure, e / spec.stoppingConversion));
      sub.CalcLegendreP(configure.maxLOrder, theCNuc, qEffect);
      // A shifted segment's sub-points are moved and energised by the shift step.
      if (!it.shifted) {
        sub.CalcEDependentValues(theCNuc, configure);
        sub.CalcCoulombAmplitude(theCNuc);
      }
      point->AddSubPoint(sub);
    }
    changed++;
    newSubPoints += (long)grid.size();
  };
#ifdef _OPENMP
  const bool serial = omp_in_parallel();
#else
  const bool serial = true;
#endif
  if (serial) {
    for (size_t q = 0; q < work.size(); q++) replace(work[q], q);
  } else {
#pragma omp parallel for schedule(dynamic, 4)
    for (int q = 0; q < (int)work.size(); q++) replace(work[q], (size_t)q);
  }
  if (changed > 0)
    for (const Item &it : work)
      if (it.shifted && std::find(reshift.begin(), reshift.end(), it.seg) == reshift.end()) reshift.push_back(it.seg);

  if (std::getenv("AZURE_GRID_TIMING")) {
    static std::atomic<int> calls{0};
    const int c = ++calls;
    if (c <= 20 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      int resized = 0;
      double maxMove = 0.0;
      for (double m : moved) {
        if (m < 0.0) resized++;
        else maxMove = std::max(maxMove, m);
      }
      std::cerr << "[grid] call " << c << ": " << work.size() << " grids, " << changed.load() << " re-placed ("
                << newSubPoints.load() << " sub-points; " << resized << " with a new node count, largest node move "
                << maxMove * 1.0e6 << " eV otherwise); observed parameters " << ms(t0, t1) << " ms, grids "
                << ms(t1, std::chrono::steady_clock::now()) << " ms" << std::endl;
    }
  }
  return changed.load();
}

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
      if (targetEffect.HasScreeningError()) {
        std::cout << "ERROR: <targetInt> line \"" << line << "\": a screening block is 'screening Ue_eV vary sigma_eV'"
                  << std::endl;
        return -1;
      }
      if (stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      if (targetEffect.IsActive()) {
        this->AddTargetEffect(targetEffect);
        TargetEffect *thisTargetEffect = this->GetTargetEffect(this->NumTargetEffects());
        std::vector<int> segmentsList = thisTargetEffect->GetSegmentsList();
        for (int i = 1; i <= segmentsList.size(); i++) {
          if (this->IsSegmentKey(segmentsList[i - 1])) {
            ESegment *segment = this->GetSegmentFromKey(segmentsList[i - 1]);
            if (segment) {
              segment->SetTargetEffectNum(this->NumTargetEffects());
              // If segment has components, set the target effect to these as well
              if (segment->HasComponents()) {
                for (auto component : segment->GetComponentSegments()) {
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
    PPair *exitPair = compound->GetPairFromKey(segment->GetExitKey());  // nullptr for total capture (-1)
    double cmConversion;
    if (entrancePair->GetPType() == 20 && exitPair)
      cmConversion = (exitPair->GetM(1) + exitPair->GetM(2)) / exitPair->GetM(2);
    else
      cmConversion = entrancePair->GetM(2) / (entrancePair->GetM(1) + entrancePair->GetM(2));

    if (segment->IsTargetEffect()) {
      TargetEffect *targetEffect = this->GetTargetEffect(segment->GetTargetEffectNum());
      // the line's energies stay in the lab frame; every point converts them with its own pair's factor
      const int bohrZ1 = (entrancePair->GetPType() == 20) ? 0 : entrancePair->GetZ(1);
      const int bohrZ2 = entrancePair->GetZ(2);
      if (targetEffect->IsUdr()) {
        double eMin = 1.0e300, eMax = 0.0;
        for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
          eMin = std::min(eMin, point->GetLabEnergy());
          eMax = std::max(eMax, point->GetLabEnergy());
        }
        if (eMax > 0.0 && !targetEffect->LoadUdrTable(eMin, eMax, configure.outStream)) return -1;
      }

      for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
        // An effect restricted to energy ranges leaves points outside them
        // completely untouched: they never receive the effect number, so the
        // whole downstream machinery treats them as ordinary points.
        double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy());
        if (blendWeight <= 0.) continue;
        point->SetTargetEffectNum(segment->GetTargetEffectNum());
        point->SetTargetBlendWeight(blendWeight);
        point->SetTargetKinematics(cmConversion, bohrZ1, bohrZ2);

        if (targetEffect->IsSubPointEffect()) {
          double forwardDepth = 0.0;
          double backwardDepth = 0.0;

          double targetThickness = 0.0;
          if (targetEffect->IsTargetIntegration()) {
            if (targetEffect->IsThickTarget()) {
              // an infinitely thick target: the point integrates from its energy down to 1 keV
              targetThickness = std::max(point->GetCMEnergy() - 0.001, 0.0);
            } else {
              targetThickness = cmConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure);
              if (!(targetThickness > 0.0)) {
                configure.outStream << "ERROR: target integration of segment #" << segment->GetSegmentKey()
                                    << " at E_lab = " << point->GetLabEnergy() << " MeV: the stopping cross section times the"
                                    << " density is " << targetThickness / cmConversion << " (not positive).  A density of 0"
                                    << " means an infinitely thick target; otherwise the stopping power must be positive"
                                    << " at every point." << std::endl;
                return -1;
              }
            }
            point->SetTargetThickness(targetThickness);
          }
          // How far the grid reaches above and below the point for this effect (beam spread, target, straggling).
          targetEffect->GridExtent(point->GetCMEnergy(), targetThickness, configure, cmConversion, bohrZ1, bohrZ2, forwardDepth, backwardDepth);

          // Generate adaptive integration grid
          double startEnergy = point->GetCMEnergy() + forwardDepth;
          double endEnergy;

          if (targetEffect->IsTargetIntegration()) {
            // the target (or, for a thick target, everything down to 1 keV), plus the spread below it
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if (endEnergy < 0.001) {
              endEnergy = 0.001;
            }
          } else {
            // Convolution or ConvCoefficients - use backward depth
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if (endEnergy < 0.001) {
              endEnergy = 0.001;
            }
          }

          if (targetEffect->IsBeamProfile()) {
            // The kernel is an absolute beam profile: sample where the beam
            // is, narrowed to the point's energy window (plus the resolution
            // tails) when it has one, not around the point's own energy.
            double low, high;
            targetEffect->BeamProfileSupport(low, high);  // lab
            low *= cmConversion; high *= cmConversion;
            double s = cmConversion * targetEffect->GetBeamTpcSigma();
            if (point->HasBinWindow()) {
              low = std::max(low, point->GetBinLowCM() - 4.0 * s);
              high = std::min(high, point->GetBinHighCM() + 4.0 * s);
            }
            if (low < 0.001) low = 0.001;
            if (high <= low) high = low + 0.001;
            startEnergy = high;
            endEnergy = low;
            point->SetPhotoKinematics(entrancePair->GetSepE() - exitPair->GetExE(),
                                      (entrancePair->GetM(1) + entrancePair->GetM(2)) * uconv);
          }
          if (targetEffect->IsUdr()) {
            // A delay tau maps the point's nominal flight time t0 to a true
            // flight time t0 - tau: the kernel's support in delay becomes a
            // window in true energy, extending above the point for positive
            // delays (moderation) and below it for negative ones.
            double mass = entrancePair->GetM(1) * uconv;
            std::shared_ptr<const UdrKernel> kernel = targetEffect->BuildUdrKernel(point->GetLabEnergy(), mass);
            point->SetUdrKernel(kernel, mass);
            double t0 = UdrKinematics::TimeOfFlight(point->GetLabEnergy(), mass, targetEffect->GetUdrFlightPath());
            double eHighLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMax, mass, targetEffect->GetUdrFlightPath());
            double eLowLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMin, mass, targetEffect->GetUdrFlightPath());
            double cmPerLab = point->GetCMEnergy() / point->GetLabEnergy();
            if (eHighLab <= point->GetLabEnergy()) eHighLab = point->GetLabEnergy() * 1.0001;
            if (eLowLab <= 0.0 || eLowLab >= point->GetLabEnergy()) eLowLab = point->GetLabEnergy() * 0.9999;
            startEnergy = eHighLab * cmPerLab;
            endEnergy = std::max(eLowLab * cmPerLab, 0.001);
          }
          std::vector<double> energyGrid;
          int numPoints = targetEffect->NumSubPoints();
          if (configure.useAdaptiveGrid) {
            AdaptiveIntegrationGrid::GridConfig gridConfig;
            gridConfig.maxPoints = numPoints;
            gridConfig.entranceKey = segment->GetEntranceKey();
            gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
            gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
            gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
            gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
            AdaptiveIntegrationGrid gridGenerator(gridConfig);
            energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
            point->SetSubGridSpec(MakeSubGridSpec(gridConfig, startEnergy, endEnergy, targetEffect, cmConversion));
          } else {
            double step = (startEnergy - endEnergy) / numPoints;
            for (int i = 0; i <= numPoints; i++)
              energyGrid.push_back(startEnergy - i * step);
          }

          for (size_t i = 0; i < energyGrid.size(); i++) {
            double subEnergy = energyGrid[i];
            EPoint subPoint(point->GetCMAngle(), subEnergy, &*segment);
            subPoint.SetCoefficientOrder(point->GetCoefficientOrder());
            subPoint.SetCoefficientAbsolute(point->GetCoefficientAbsolute());
            if (targetEffect->IsTargetIntegration()) {
              double stoppingPower = cmConversion * targetEffect->GetStoppingPowerEq()->Evaluate(configure, subEnergy / cmConversion);
              subPoint.SetStoppingPower(stoppingPower);
            }
            point->AddSubPoint(subPoint);
          }
        }
      }
    }
    if (segment->HasComponents()) {
      for (auto component : segment->GetComponentSegments()) {
        // The component's grid, kinematics and lab<->c.m. factor belong to the
        // pair the component calculates -- for a natural target another
        // compound nucleus with its own masses -- not to the base segment's.
        // The points were copied from the base segment with its c.m.
        // energies; convert them for this pair before the grid is laid out
        // (EData::Initialize repeats this, idempotently, from the lab energy).
        PPair *compEntrance = compound->GetPair(compound->GetPairNumFromKey(component->GetEntranceKey()));
        PPair *compExit = compound->GetPairFromKey(component->GetExitKey());
        double compConversion;
        if (compEntrance->GetPType() == 20 && compExit)
          compConversion = (compExit->GetM(1) + compExit->GetM(2)) / compExit->GetM(2);
        else
          compConversion = compEntrance->GetM(2) / (compEntrance->GetM(1) + compEntrance->GetM(2));
        for (EPointIterator point = component->GetPoints().begin(); point < component->GetPoints().end(); point++) {
          if (compEntrance->GetPType() == 20)
            point->ConvertDecayEnergy(compExit);
          else
            point->ConvertLabEnergy(compEntrance);
        }
        if (component->IsTargetEffect()) {
          TargetEffect *targetEffect = this->GetTargetEffect(component->GetTargetEffectNum());
          const int bohrZ1 = (compEntrance->GetPType() == 20) ? 0 : compEntrance->GetZ(1);
          const int bohrZ2 = compEntrance->GetZ(2);
          if (targetEffect->IsUdr()) {
            double eMin = 1.0e300, eMax = 0.0;
            for (EPointIterator point = component->GetPoints().begin(); point < component->GetPoints().end(); point++) {
              eMin = std::min(eMin, point->GetLabEnergy());
              eMax = std::max(eMax, point->GetLabEnergy());
            }
            if (eMax > 0.0 && !targetEffect->LoadUdrTable(eMin, eMax, configure.outStream)) return -1;
          }

          for (EPointIterator point = component->GetPoints().begin(); point < component->GetPoints().end(); point++) {
            double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy());
            if (blendWeight <= 0.) continue;
            point->SetTargetEffectNum(component->GetTargetEffectNum());
            point->SetTargetBlendWeight(blendWeight);
            point->SetTargetKinematics(compConversion, bohrZ1, bohrZ2);

            if (targetEffect->IsSubPointEffect()) {
              double forwardDepth = 0.0;
              double backwardDepth = 0.0;

              double targetThickness = 0.0;
              if (targetEffect->IsTargetIntegration()) {
                if (targetEffect->IsThickTarget()) {
                  targetThickness = std::max(point->GetCMEnergy() - 0.001, 0.0);
                } else {
                  targetThickness = compConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure);
                  if (!(targetThickness > 0.0)) {
                    configure.outStream << "ERROR: target integration of a component of segment #" << segment->GetSegmentKey()
                                        << " at E_lab = " << point->GetLabEnergy() << " MeV: the stopping cross section times"
                                        << " the density is " << targetThickness / compConversion << " (not positive)." << std::endl;
                    return -1;
                  }
                }
                point->SetTargetThickness(targetThickness);
              }
              // How far the grid reaches above and below the point for this effect (beam spread, target, straggling).
              targetEffect->GridExtent(point->GetCMEnergy(), targetThickness, configure, compConversion, bohrZ1, bohrZ2, forwardDepth, backwardDepth);

              // Generate adaptive integration grid for component segment
              double startEnergy = point->GetCMEnergy() + forwardDepth;
              double endEnergy;

              if (targetEffect->IsTargetIntegration()) {
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if (endEnergy < 0.001) {
                  endEnergy = 0.001;
                }
              } else {
                // Convolution or ConvCoefficients - use backward depth
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if (endEnergy < 0.001) {
                  endEnergy = 0.001;
                }
              }

              if (targetEffect->IsBeamProfile()) {
                double low, high;
                targetEffect->BeamProfileSupport(low, high);  // lab
                low *= compConversion; high *= compConversion;
                double s = compConversion * targetEffect->GetBeamTpcSigma();
                if (point->HasBinWindow()) {
                  low = std::max(low, point->GetBinLowCM() - 4.0 * s);
                  high = std::min(high, point->GetBinHighCM() + 4.0 * s);
                }
                if (low < 0.001) low = 0.001;
                if (high <= low) high = low + 0.001;
                startEnergy = high;
                endEnergy = low;
                point->SetPhotoKinematics(compEntrance->GetSepE() - compExit->GetExE(),
                                          (compEntrance->GetM(1) + compEntrance->GetM(2)) * uconv);
              }
              if (targetEffect->IsUdr()) {
                double mass = compEntrance->GetM(1) * uconv;
                std::shared_ptr<const UdrKernel> kernel = targetEffect->BuildUdrKernel(point->GetLabEnergy(), mass);
                point->SetUdrKernel(kernel, mass);
                double t0 = UdrKinematics::TimeOfFlight(point->GetLabEnergy(), mass, targetEffect->GetUdrFlightPath());
                double eHighLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMax, mass, targetEffect->GetUdrFlightPath());
                double eLowLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMin, mass, targetEffect->GetUdrFlightPath());
                double cmPerLab = point->GetCMEnergy() / point->GetLabEnergy();
                if (eHighLab <= point->GetLabEnergy()) eHighLab = point->GetLabEnergy() * 1.0001;
                if (eLowLab <= 0.0 || eLowLab >= point->GetLabEnergy()) eLowLab = point->GetLabEnergy() * 0.9999;
                startEnergy = eHighLab * cmPerLab;
                endEnergy = std::max(eLowLab * cmPerLab, 0.001);
              }
              std::vector<double> energyGrid;
              int numPoints = targetEffect->NumSubPoints();
              if (configure.useAdaptiveGrid) {
                AdaptiveIntegrationGrid::GridConfig gridConfig;
                gridConfig.maxPoints = numPoints;
                // The component's own entrance pair: its resonances, on its energy axis, are what the grid must resolve.
                gridConfig.entranceKey = component->GetEntranceKey();
                gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
                gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
                gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
                gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
                AdaptiveIntegrationGrid gridGenerator(gridConfig);
                energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
                point->SetSubGridSpec(MakeSubGridSpec(gridConfig, startEnergy, endEnergy, targetEffect, compConversion));
              } else {
                double step = (startEnergy - endEnergy) / numPoints;
                for (int i = 0; i <= numPoints; i++)
                  energyGrid.push_back(startEnergy - i * step);
              }

              for (size_t i = 0; i < energyGrid.size(); i++) {
                double subEnergy = energyGrid[i];
                EPoint subPoint(point->GetCMAngle(), subEnergy, &*component);
                subPoint.SetCoefficientOrder(point->GetCoefficientOrder());
                subPoint.SetCoefficientAbsolute(point->GetCoefficientAbsolute());
                if (targetEffect->IsTargetIntegration()) {
                  double stoppingPower = compConversion * targetEffect->GetStoppingPowerEq()->Evaluate(configure, subEnergy / compConversion);
                  subPoint.SetStoppingPower(stoppingPower);
                }
                point->AddSubPoint(subPoint);
              }
            }
          }
        }
      }
    }
  }
  screeningEffects_.clear();
  for (int n = 1; n <= this->NumTargetEffects(); n++)
    if (this->GetTargetEffect(n)->HasScreening()) screeningEffects_.push_back(n);
  this->ApplyScreening();
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
  // Labels of the compound nuclei, for the sample decomposition files.
  nucleusLabels_.clear();
  for (int k = 1; k <= compound->NumNuclei(); k++) nucleusLabels_.push_back(compound->NucleusLabel(k));
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

  // A phase-shift segment is one partial wave: (J, l) and the channel spin s.  Up to 2026-10-27 the
  // diagonal U-matrix elements of every channel spin with the segment's (J, l) were summed and half
  // the argument reported, which is not a phase shift of anything when two spins contribute (p + 3He
  // 1P1 and 3P1, p + 7Be J = 1 l = 1 with s = 1 and 2; audit 2026-10-25, tier 3 entry 4).  A segment
  // whose (J, l) has several channel spins in the level scheme must say which with "spin <s>", and a
  // given s must exist.
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (!seg->IsPhase()) continue;
    int aa = compound->GetPairNumFromKey(seg->GetEntranceKey());
    if (aa == 0) continue;
    std::vector<double> spins;
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      JGroup *jg = compound->GetJGroup(j);
      if (jg->GetJ() != seg->GetJ()) continue;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        AChannel *channel = jg->GetChannel(ch);
        if (channel->GetRadType() != 'P' || channel->GetPairNum() != aa || channel->GetL() != seg->GetL()) continue;
        bool seen = false;
        for (double v : spins) if (std::fabs(v - channel->GetS()) < 1.0e-6) seen = true;
        if (!seen) spins.push_back(channel->GetS());
      }
    }
    const double want = seg->GetPhaseS();
    if (want >= 0.0) {
      bool found = false;
      for (double v : spins) if (std::fabs(v - want) < 1.0e-6) found = true;
      if (!found) {
        configure.outStream << "ERROR: phase-shift segment #" << seg->GetSegmentKey() << " (J = " << seg->GetJ() << ", l = " << seg->GetL()
                            << ") asks for channel spin s = " << want << ", but the level scheme has ";
        if (spins.empty()) configure.outStream << "no channel of that (J, l) for its entrance pair." << std::endl;
        else {
          configure.outStream << "only s =";
          for (double v : spins) configure.outStream << " " << v;
          configure.outStream << " for it." << std::endl;
        }
        return -1;
      }
    } else if (spins.size() > 1) {
      configure.outStream << "ERROR: phase-shift segment #" << seg->GetSegmentKey() << " (J = " << seg->GetJ() << ", l = " << seg->GetL()
                          << ") has " << spins.size() << " channel spins in the level scheme (s =";
      for (double v : spins) configure.outStream << " " << v;
      configure.outStream << "): a phase shift belongs to one partial wave, so add \"spin <s>\" to the segment line"
                          << " (older builds summed the S-matrix elements of all of them, which is not a phase shift)." << std::endl;
      return -1;
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
        // Same columns as an unmapped point: a mapped point has its own
        // (identical) sub-point grid.
        out << std::setw(12) << tempMap
            << std::setw(18) << data.point()->NumSubPoints();
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
      PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(data.segment()->GetEntranceKey()));
      double inEnergy = data.point()->GetCMEnergy() + entrancePair->GetSepE();
      for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
        if (theCNuc->GetJGroup(j)->IsInRMatrix() && theCNuc->GetJGroup(j)->GetNucleus() == entrancePair->GetNucleus()) {
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
          double fitCrossSection = point->GetFitCrossSection() / point->GetCrossSectionKinFactor();
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
    if (segment->IsNatural()) WriteNaturalDecomposition(configure, output.IsExtrap(), kinflag, nucleusLabels_, &*segment);
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
    // The energy-shift penalties the fit minimised (ESegment::AddShiftPenalty), when
    // any shift is a fit parameter.  Without this the file reported less than the
    // minimiser's own total by exactly this amount.
    bool anyShiftVaried = false;
    double totalShiftChiSquared = 0.;
    for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
      if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt()) anyShiftVaried = true;
      if (segment->OwnsShiftFn())
        for (int q = 0; q < segment->NumShiftFnParams(); q++)
          if (segment->IsVaryShiftFn(q)) anyShiftVaried = true;
      segment->AddShiftPenalty(totalShiftChiSquared);
    }
    if (anyShiftVaried) chiOut << " Total-Shift-Chi-Squared: " << totalShiftChiSquared;
    // Isotopic abundances: their priors and the sum-to-one term.
    if (abundancePenalty_ != 0.0)
      chiOut << " Total-Abundance-Chi-Squared: " << abundancePenalty_;
    // Electron-screening priors, when a screening potential has one.
    {
      bool anyPrior = false;
      for (int k = 1; k <= this->NumScreeningParams(); k++)
        if (this->GetTargetEffect(this->GetScreeningEffect(k))->GetScreeningSigma() > 0.0) anyPrior = true;
      if (anyPrior) chiOut << " Total-Screening-Chi-Squared: " << this->ScreeningPenalty();
    }
    // Measured scattering lengths (<scatteringLengths>; table in scatteringLengths.out).
    if (scatteringLengthChi2_ != 0.0)
      chiOut << " Total-Scattering-Length-Chi-Squared: " << scatteringLengthChi2_;
    chiOut << std::endl
           << std::endl;
    chiOut.flush();
    chiOut.close();
    if (this->HasCovariance()) this->WriteCovarianceSummary(configure);
  }
  if (isVaryNorm) {
    std::string outputfile = configure.outputdir + "normalizations.out";
    std::ofstream out(outputfile.c_str());
    if (out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      // One complete row per segment whose normalization is varied (a norm-group member keeps its own
      // flag).  Before 2026-10-27 only the first four fields were conditional, so every other segment
      // left a row of six fields without key, file or angles.
      for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
        if (!segment->IsVaryNorm()) continue;
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
    if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt() || segment->IsShiftGroupFollower() ||
        segment->HasShiftFn()) {
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
        if (!(segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt() || segment->IsShiftGroupFollower() ||
              segment->HasShiftFn()))
          continue;
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
      // Shift functions: one line per function, after the per-segment rows.  The
      // lines start with '#' so that readers of the rows above skip them.
      bool header = false;
      double totalShiftChi = 0.;
      for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++)
        segment->AddShiftPenalty(totalShiftChi);
      for (int pos = 0; pos < (int)segments_.size(); pos++) {
        ESegment &owner = segments_[pos];
        if (!owner.OwnsShiftFn() || !owner.HasShiftFn()) continue;
        if (!header) {
          out << "# shift functions: dE(x) = shift + f(x), x = original lab energy of a point (MeV), dE in MeV" << std::endl
              << "# segment_key_#, member_keys, f(x), shift, a0 a1 ..., x_low, x_high, dE(x_low), dE(x_high), penalty"
              << std::endl;
          header = true;
        }
        std::ostringstream keys;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (int m = 0; m < (int)segments_.size(); m++) {
          if (segments_[m].GetShiftFnSourcePos() != pos) continue;
          keys << (m == pos ? "" : " ") << segments_[m].GetSegmentKey();
          for (int pid = 1; pid <= segments_[m].NumPoints(); pid++) {
            EPoint *pt = segments_[m].GetPoint(pid);
            const double x = pt->GetOriginalEnergy() > 0 ? pt->GetOriginalEnergy() : pt->GetLabEnergy();
            lo = std::min(lo, x);
            hi = std::max(hi, x);
          }
        }
        // The function's own penalty: the owner's shift penalties minus those of its
        // constant and sqrt(E) terms.
        double pen = 0.0, penOwn = 0.0;
        owner.AddShiftPenalty(pen);
        if (owner.IsVaryEnergyShift() && owner.GetEnergyShiftError() != 0.0)
          penOwn = pow(owner.GetEnergyShift() / owner.GetEnergyShiftError(), 2.0);  // from 0, as AddShiftPenalty
        out << "# " << owner.GetSegmentKey() << "," << keys.str() << "," << owner.GetShiftFnExpression() << ","
            << owner.GetEnergyShift() << ",";
        for (int q = 0; q < owner.NumShiftFnParams(); q++) out << (q ? " " : "") << owner.GetShiftFnValues()[q];
        out << "," << lo << "," << hi << "," << owner.TotalEnergyShift(lo) << "," << owner.TotalEnergyShift(hi) << ","
            << pen - penOwn << std::endl;
      }
      if (header) out << "# Total-Shift-Chi-Squared: " << totalShiftChi << std::endl;
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
        // Every point carries its own amplitudes, mapped or not (see
        // CalculateECAmplitudes), so the intEC file does not depend on mapping.
        for (auto point = segment->GetPoints().begin(); point != segment->GetPoints().end(); point++) {
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
  // v2 (2026-10-23): the long-wavelength switch is part of the signature, and
  // the M1 amplitudes changed (sqrt(16/3) prefactor, no E1 charge function), so
  // every v1 file is recomputed once.
  // v3 (2026-10-26): the M1 spin-term sign for a cluster whose partner has
  // half-integer spin (ECIntegral.cpp), so every v2 file is recomputed once.
  // v4 (2026-10-26): the M1 orbital term is dropped on channel-spin-flip pathways (ECIntegral.cpp), so every
  // v3 file with such pathways is recomputed once.
  hash.Add("AZURE2 EC signature v4");
  hash.Add((configure.paramMask & Config::USE_GSL_COULOMB_FUNC) ? 1 : 0);
  hash.Add((configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX) ? 1 : 0);
  // The convection current changes exact-form electric terms only; long-wavelength signatures stay as in v2.
  if (configure.ecConvectionCurrent && !(configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX))
    hash.Add("convection current");
  hash.Add((configure.paramMask & Config::CALCULATE_REACTION_RATE) ? 1 : 0);
  hash.Add(configure.useHybridMethod ? 1 : 0);
  // Relativistic kinematics (relativistic=1 in <config>) change eta, rho and the flux factor of every amplitude.
  // Added only when on, so every non-relativistic v4 file keeps its signature and is not recomputed; a relativistic
  // run no longer reads back a file written without the switch (2026-10-27, checklist t3-25).
  if (configure.relativistic) hash.Add("relativistic");

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
  ECTable::ResetCounters();

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
                  // Mapped points compute their own amplitudes too: the
                  // amplitudes are a one-off setup cost, and the intEC file
                  // then has the same layout with and without point mapping.
                  {
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
                {
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
  if (!usePrevious && configure.useECTable && !(configure.paramMask & Config::USE_API) &&
      ECTable::Lookups() + ECTable::DirectCalls() > 0)
    configure.outStream << "\t" << ECTable::Report() << std::endl;
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
          const bool inRange = (entranceKey == exitKey) ? point->ConvertLabAngle(entrancePair)
                                                        : point->ConvertLabAngle(entrancePair, exitPair, configure);
          if (!inRange) {
            configure.outStream << "ERROR: " << "component (" << entranceKey << " -> " << exitKey << ") of segment #" << componentSegment.GetSegmentKey() << ": lab angle " << point->GetLabAngle() << " deg at E_lab = " << point->GetLabEnergy()
    << " MeV: " << point->GetLabAngleWhy() << "." << std::endl;
            return -1;
          }
          point->ConvertCrossSection(entrancePair, exitPair);
        }

        if (exitPair->GetPType() == 10 && componentSegment.IsDifferential() &&
            !componentSegment.IsPhase() && !componentSegment.IsAngularDist() && !componentSegment.IsCMDifferential()) {
          point->ConvertLabAngleGammas(entrancePair);
          point->ConvertCrossSectionGammas(entrancePair);
        }
        // The sub-points of a target effect were laid out (ReadTargetEffectsFile) from the copied base
        // point's c.m. angle, before this conversion with the component's own pair; put them at the
        // converted angle (CalcLegendreP below recurses into them).  Until 2026-10-27 a component with a
        // different exit pair, or a fixed-angle RATIO denominator (raw lab angle), was integrated at the
        // wrong angle (audit 2026-10-25, item 32).
        for (int k = 1; k <= point->NumSubPoints(); k++) {
          point->GetSubPoint(k)->SetLabAngle(point->GetLabAngle());
          point->GetSubPoint(k)->SetCMAngle(point->GetCMAngle());
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
      if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
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
                    if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
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
                  if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
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
/*!
 * Checks a component's pair keys before the component is created: both must
 * exist and belong to one compound nucleus (a reaction does not change the
 * compound nucleus), and the component's entrance pair is flagged as an
 * entrance HERE, while CNuc::Initialize -- which builds the reaction pathways
 * of every entrance pair -- still lies ahead.  Flagging it only in
 * EData::Initialize, as was done, left a pair used by components alone (the
 * second isotope of a natural target) without pathways: GetDecay(0).
 * Returns false, with the error printed, when the component cannot be made.
 */
bool EData::CheckComponentPairs(CNuc *theCNuc, const ESegment &baseSegment, int entranceKey, int exitKey,
                                const Config &configure) {
  if (!theCNuc->IsPairKey(entranceKey) || !theCNuc->IsPairKey(exitKey)) {
    configure.outStream << "ERROR: segment #" << baseSegment.GetSegmentKey() << ": component " << entranceKey
                        << " -> " << exitKey << " names a pair key that is not in the compound nucleus." << std::endl;
    componentError_ = true;
    return false;
  }
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(entranceKey));
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(exitKey));
  if (entrancePair->GetNucleus() != exitPair->GetNucleus()) {
    configure.outStream << "ERROR: segment #" << baseSegment.GetSegmentKey() << ": component " << entranceKey
                        << " -> " << exitKey << " connects two compound nuclei ("
                        << theCNuc->NucleusLabel(entrancePair->GetNucleus()) << " and "
                        << theCNuc->NucleusLabel(exitPair->GetNucleus()) << ")." << std::endl;
    componentError_ = true;
    return false;
  }
  entrancePair->SetEntrance();
  return true;
}

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
  componentSegment->SetNatural(false);   // the weighting is the base segment's job
  componentSegment->ApplyTotalToPoints();

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
  componentSegment->SetNatural(false);   // the weighting is the base segment's job
  componentSegment->ApplyTotalToPoints();

  return componentSegment;
}

namespace {

// True when the calculation of `b` can be taken from that of `a`: both are at
// the same energy, with the same entrance and exit pair, so they share the
// T-matrix of the point and -- checked here, sub-point by sub-point -- of the
// integration grid.  Only what enters through GenMatrixFunc::
// CalculateCrossSection on the point itself (angle, Coulomb amplitude,
// observable type) may differ.  Whether `a` and `b` may sit in different
// segments is the caller's question (ShiftsAgree).
bool SharesCalculation(EPoint &a, EPoint &b, bool sameExit = true) {
  if (a.GetEntranceKey() != b.GetEntranceKey() || (sameExit && a.GetExitKey() != b.GetExitKey())) return false;
  if (a.GetOriginalEnergy() != b.GetOriginalEnergy() || a.GetLabEnergy() != b.GetLabEnergy() ||
      a.GetCMEnergy() != b.GetCMEnergy())
    return false;
  if (a.IsTargetEffect() != b.IsTargetEffect() || a.GetTargetEffectNum() != b.GetTargetEffectNum() ||
      a.GetTargetBlendWeight() != b.GetTargetBlendWeight())
    return false;
  if (a.IsAnalyzingPower() != b.IsAnalyzingPower()) return false;
  if (a.IsAngularDist() != b.IsAngularDist()) return false;
  if (a.IsCoefficient() != b.IsCoefficient()) return false;
  if (a.IsCoefficient() && a.GetCoefficientAbsolute() != b.GetCoefficientAbsolute()) return false;
  if (a.NumSubPoints() != b.NumSubPoints()) return false;
  std::vector<EPoint> &sa = a.GetSubPoints();
  std::vector<EPoint> &sb = b.GetSubPoints();
  for (size_t i = 0; i < sa.size(); i++) {
    if (sa[i].GetOriginalEnergy() != sb[i].GetOriginalEnergy() ||
        sa[i].GetCMEnergy() != sb[i].GetCMEnergy() ||
        sa[i].GetLabEnergy() != sb[i].GetLabEnergy() ||
        sa[i].GetStoppingPower() != sb[i].GetStoppingPower())
      return false;
  }
  return true;
}

// The energy shift of a segment as a function of the fit parameters: which
// free parameters move it, and its current fixed form.  Two segments whose
// shifts agree here are at the same energy for EVERY value of the free shift
// parameters, not only the current one: a finite-difference step on a shift
// (AZURECalc's segment-local shift derivatives) moves both together, so a
// mapping between them stays valid without being rebuilt.
struct ShiftForm {
  int constOwner = -1;  // segment position of the free constant shift, -1 if fixed
  bool sqrtFree = false;
  int fnOwner = -1;     // segment position of the shift function with a free parameter, -1 if none
  double constShift = 0.0, sqrtShift = 0.0;
  bool hasFn = false;
  std::string fnExpression;
  std::vector<double> fnValues;
};

ShiftForm ShiftFormOf(std::vector<ESegment> &segments, int pos) {
  ESegment &seg = segments[pos];
  ShiftForm f;
  const int lead = seg.GetShiftLeaderPos() >= 0 ? seg.GetShiftLeaderPos() : pos;
  if (segments[lead].IsVaryEnergyShift()) f.constOwner = lead;
  f.sqrtFree = seg.IsVaryEnergyShiftSqrt();
  f.constShift = seg.GetEnergyShift();
  f.sqrtShift = seg.GetEnergyShiftSqrt();
  f.hasFn = seg.HasShiftFn();
  if (f.hasFn) {
    f.fnExpression = seg.GetShiftFnExpression();
    f.fnValues = seg.GetShiftFnValues();
    const int src = seg.GetShiftFnSourcePos();
    if (src >= 0)
      for (int j = 0; j < segments[src].NumShiftFnParams(); j++)
        if (segments[src].IsVaryShiftFn(j)) f.fnOwner = src;
  }
  return f;
}

// A free sqrt(E) coefficient is the segment's own parameter, so a segment with
// one agrees with no other segment.
bool ShiftsAgree(const ShiftForm &a, const ShiftForm &b) {
  if (a.sqrtFree || b.sqrtFree) return false;
  return a.constOwner == b.constOwner && a.fnOwner == b.fnOwner && a.constShift == b.constShift &&
         a.sqrtShift == b.sqrtShift && a.hasFn == b.hasFn && a.fnExpression == b.fnExpression &&
         a.fnValues == b.fnValues;
}

}  // namespace

/*!
 * Maps points onto other points at the same energy, so the energy-dependent
 * part of the calculation (the energy-dependent values, and the level matrix
 * and T-matrix of the point and of each target-effect sub-point) is done once
 * per energy; each mapped point then only forms its own observable (angle,
 * Coulomb amplitude, A_y) from it and integrates its own target effect
 * (EPoint::Calculate).
 *
 * Points are mapped within a segment and across segments.  Two points share a
 * calculation (SharesCalculation) when they have the same entrance and exit
 * pair, the same original, lab and c.m. energies, the same target effect and
 * sub-point grid, and the same kind of observable (cross section or A_y,
 * angular-distribution coefficients or not).  Across segments their energy
 * shifts must also agree as functions of the fit parameters (ShiftsAgree):
 * the same free parameters and the same fixed values.  Then the two segments
 * are at the same energy whatever a fit or a finite-difference step does to
 * the shifts.  Segments whose shifts merely coincide at the current values are
 * mapped too, but only while they coincide: every change of a shift or of a
 * sub-point grid rebuilds the mapping (EData::FillEnergyShiftsFromParams).
 * Normalizations play no part: they scale the data, not the model.
 *
 * The host of a group is its first point in data order (segment, then point
 * index): every loop that calculates the whole data set walks it in that order,
 * so a host is calculated before the points that read from it
 * (ESegment::CalculateTheoreticalCrossSection; ESegment::
 * UpdatePointEnergiesWithShift for the energy-dependent values, within a
 * segment).  The parallel loops (AZURECalc::PrecalculatePoints) hand out hosts
 * only, each of which calculates its mapped points.  The analytic gradient
 * (AZUREGrad.cpp) works on the same groups: one forward solve per group.
 *
 * Segments with components (SUM/RATIO, total capture) are not mapped; a
 * shifted segment without a particle exit (external capture) only within itself.
 * AZURE_NO_POINT_MAPPING=1 in the environment switches mapping off entirely;
 * AZURE_NO_CROSS_SEGMENT_MAPPING=1 keeps it within each segment only.
 */

void EData::MapData(CNuc *theCNuc) {
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++) {
      point->ClearMapping();
      point->ClearLocalMappedPoints();
    }
  }
  for (ESegment &component : componentSegments_) {
    for (EPointIterator point = component.GetPoints().begin(); point < component.GetPoints().end(); point++) {
      point->ClearMapping();
      point->ClearLocalMappedPoints();
    }
  }
  const char *off = std::getenv("AZURE_NO_POINT_MAPPING");
  if (off && *off && std::string(off) != "0") return;
  const char *crossOff = std::getenv("AZURE_NO_CROSS_SEGMENT_MAPPING");
  const bool crossSegment = !(crossOff && *crossOff && std::string(crossOff) != "0");

  std::vector<ESegment> &segments = GetSegments();
  const int numSegments = (int)segments.size();
  std::vector<ShiftForm> shiftForms(numSegments);
  for (int s = 0; s < numSegments; s++) shiftForms[s] = ShiftFormOf(segments, s);
  // External-capture amplitudes at a shifted energy may be interpolated from
  // the energies of the segment's own points (ECAmplitudeCache), so a shifted
  // capture segment is mapped within itself only.  Without the compound
  // nucleus at hand, every shifted segment is.
  std::vector<char> ownOnly(numSegments, 0);
  for (int s = 0; s < numSegments; s++) {
    const ShiftForm &f = shiftForms[s];
    const bool shifted = f.constOwner >= 0 || f.sqrtFree || f.fnOwner >= 0 || f.constShift != 0.0 ||
                         f.sqrtShift != 0.0 || f.hasFn;
    if (!shifted) continue;
    if (!theCNuc) {
      ownOnly[s] = 1;
      continue;
    }
    const int exitPair = theCNuc->GetPairNumFromKey(segments[s].GetExitKey());
    if (exitPair < 1 || theCNuc->GetPair(exitPair)->GetPType() != 0) ownOnly[s] = 1;
  }

  // Every point of every mappable segment, in data order.
  struct Ref { int seg, point; };
  std::vector<Ref> refs;
  for (int s = 0; s < numSegments; s++) {
    if (segments[s].HasComponents() || segments[s].IsTotalCapture()) continue;
    for (int p = 0; p < segments[s].NumPoints(); p++) refs.push_back({s, p});
  }
  auto pointOf = [&segments](const Ref &r) -> EPoint & { return segments[r.seg].GetPoints()[r.point]; };

  // Group by pair and energy; a stable sort keeps data order inside a group,
  // so the first point of a group is its first point in data order.
  std::stable_sort(refs.begin(), refs.end(), [&](const Ref &x, const Ref &y) {
    EPoint &a = pointOf(x), &b = pointOf(y);
    if (a.GetEntranceKey() != b.GetEntranceKey()) return a.GetEntranceKey() < b.GetEntranceKey();
    if (a.GetExitKey() != b.GetExitKey()) return a.GetExitKey() < b.GetExitKey();
    if (a.GetOriginalEnergy() != b.GetOriginalEnergy()) return a.GetOriginalEnergy() < b.GetOriginalEnergy();
    return a.GetCMEnergy() < b.GetCMEnergy();
  });

  std::vector<Ref> hosts;
  const int n = (int)refs.size();
  int start = 0;
  while (start < n) {
    EPoint &first = pointOf(refs[start]);
    int end = start + 1;
    while (end < n) {
      EPoint &e = pointOf(refs[end]);
      if (e.GetEntranceKey() != first.GetEntranceKey() || e.GetExitKey() != first.GetExitKey() ||
          e.GetOriginalEnergy() != first.GetOriginalEnergy() || e.GetCMEnergy() != first.GetCMEnergy())
        break;
      end++;
    }
    hosts.clear();
    for (int k = start; k < end; k++) {
      EPoint &p = pointOf(refs[k]);
      bool mapped = false;
      for (const Ref &h : hosts) {
        if (h.seg != refs[k].seg && !(crossSegment && !ownOnly[h.seg] && !ownOnly[refs[k].seg] &&
                                      ShiftsAgree(shiftForms[h.seg], shiftForms[refs[k].seg])))
          continue;
        EPoint &hp = pointOf(h);
        if (SharesCalculation(hp, p)) {
          p.SetMap(h.seg + 1, h.point + 1);
          hp.AddLocalMappedPoint(&p);
          mapped = true;
          break;
        }
      }
      if (!mapped) hosts.push_back(refs[k]);
    }
    start = end;
  }

  MapComponentSegments(theCNuc);
}

/*!
 * Point mapping inside a segment with components (SUM, RATIO, total capture): each component's point i is mapped onto
 * the base segment's point i.  A component has its own exit pair but the base's entrance pair, energies, sub-point
 * grid and energy shift, so its level matrix is the host's: the host's calculation builds only the component's
 * T-matrix and observable (EPoint::Calculate).  Before 2026-10-24 every term was a calculation of its own, and a
 * neutron total cross section summed over nine exit pairs cost nine level-matrix solutions per point (13C+a: 87 % of
 * an evaluation).
 *
 * Such a family is not mapped with other segments.  A component re-energises its own points on an energy shift
 * (OwnsEDependentValues).  A beta-delayed entrance is left out: its energy-dependent values depend on the exit pair.
 * AZURE_NO_COMPONENT_MAPPING=1 switches this off (and AZURE_NO_POINT_MAPPING=1 all mapping).
 */
void EData::MapComponentSegments(CNuc *theCNuc) {
  if (!theCNuc) return;
  const char *off = std::getenv("AZURE_NO_COMPONENT_MAPPING");
  if (off && *off && std::string(off) != "0") return;
  std::vector<ESegment> &segments = GetSegments();
  for (int s = 0; s < (int)segments.size(); s++) {
    ESegment &seg = segments[s];
    if (!seg.HasComponents() && !seg.IsTotalCapture()) continue;
    const int entrancePair = theCNuc->GetPairNumFromKey(seg.GetEntranceKey());
    if (entrancePair < 1 || theCNuc->GetPair(entrancePair)->GetPType() == 20) continue;
    std::vector<EPoint> &points = seg.GetPoints();
    const int n = (int)points.size();
    // Each component's point i onto the base's point i (components in turn, so a host's mapped points come grouped
    // by exit and EPoint::Calculate builds each exit's T-matrix once).  The base's points stay unmapped among
    // themselves, as before: the external-capture amplitudes and their file (intEC.dat) skip mapped points, so
    // mapping them would change that file's layout.
    for (ESegment *component : seg.GetComponentSegments()) {
      if (!component || component->NumPoints() != n) continue;
      for (int p = 0; p < n; p++) {
        EPoint &host = points[p];
        EPoint *cp = component->GetPoint(p + 1);
        if (!cp || !SharesCalculation(host, *cp, false)) continue;
        cp->SetMap(s + 1, p + 1);
        cp->SetOwnsEDependentValues(true);
        host.AddLocalMappedPoint(cp);
      }
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

void EData::SetEnergyShiftFnParamOffset(int offset) {
  energyShiftFnParamOffset_ = offset;
}

/*!
 * Links the segments of each shift group ("shiftgroup <label>" on the segment
 * line) to the group's first segment in file order: that segment keeps its own
 * constant-shift parameter (value, uncertainty, vary flag, penalty), and every
 * other member takes that parameter's value instead of its own, which is fixed
 * in the fit.  A note per group lists the members.
 */

int EData::ResolveShiftGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    ESegment &seg = segments_[pos];
    const std::string label = seg.GetShiftGroup();
    if (label.empty()) continue;
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    ESegment &leader = segments_[it->second];
    seg.SetShiftGroupFollower(it->second, leader.GetEnergyShift());
    leader.SetHasShiftGroupFollowers(true);
  }
  for (const auto &g : leaders) {
    const ESegment &leader = segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (const ESegment &seg : segments_)
      if (seg.GetShiftLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.outStream << "Note: shift group '" << g.first << "' has only segment #" << leader.GetSegmentKey()
                          << "; nothing is shared." << std::endl;
      continue;
    }
    configure.outStream << "Note: shift group '" << g.first << "': segments #" << leader.GetSegmentKey() << " +"
                        << members.str() << " share one constant energy shift (segment #"
                        << leader.GetSegmentKey() << "'s parameter, " << (leader.IsVaryEnergyShift() ? "free" : "fixed")
                        << ")." << std::endl;
  }
  return ResolveShiftFunctions(configure);
}

namespace {

// Indices j of the parameters "aj" an Equation expression refers to.  The Equation
// class reads "a" followed by digits as a parameter (asin, acos, atan are not).
std::set<int> ShiftFnParamRefs(const std::string &expression) {
  std::set<int> refs;
  for (size_t k = 0; k + 1 < expression.size(); k++) {
    if (expression[k] != 'a' || !std::isdigit((unsigned char)expression[k + 1])) continue;
    size_t e = k + 1;
    while (e < expression.size() && std::isdigit((unsigned char)expression[e])) e++;
    refs.insert(std::atoi(expression.substr(k + 1, e - k - 1).c_str()));
    k = e - 1;
  }
  return refs;
}

}  // namespace

/*!
 * Shift functions ("shiftfn" block, design note 10-13-26): the segment declaring one
 * owns its parameters; it and every follower of its shift group apply
 * dE(x) = a + f(x; a0..), x the point's original lab energy.  Rejected: a follower
 * declaring a function (the group's calibration belongs to its first segment), a
 * function together with a sqrt(E) term on any member (two descriptions of one
 * calibration), and an expression referring to an undeclared parameter or leaving a
 * declared one unused.  A shift-function group ("shiftfngroup <label>") also applies
 * the function of its first segment to every member, but each member keeps its own
 * constant shift; such a group may not overlap a shift group.  For each function the
 * free shift terms (each free constant a of its members, and every varied aj) are
 * checked for degeneracy over the members' data energies; a term with no effect there
 * is an error.  One log line per function.
 */

int EData::ResolveShiftFunctions(const Config &configure) {
  // Shift-function groups: the first segment with a label (file order) declares the
  // function; the others apply it and keep their own constant shifts.
  std::vector<int> fnGroupLeader(segments_.size(), -1);
  {
    std::vector<std::pair<std::string, int>> fnLeaders;  // label, leader position
    for (int pos = 0; pos < (int)segments_.size(); pos++) {
      const ESegment &seg = segments_[pos];
      const std::string label = seg.GetShiftFnGroup();
      if (label.empty()) continue;
      if (!seg.GetShiftGroup().empty()) {
        configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " is in shift group '"
                            << seg.GetShiftGroup() << "' and in shift-function group '" << label
                            << "': use one of them (shiftgroup shares the constant shift as well)." << std::endl;
        return -1;
      }
      auto it = std::find_if(fnLeaders.begin(), fnLeaders.end(),
                             [&label](const std::pair<std::string, int> &e) { return e.first == label; });
      if (it == fnLeaders.end()) {
        fnLeaders.push_back({label, pos});
        if (!seg.OwnsShiftFn()) {
          configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " is the first segment of "
                              << "shift-function group '" << label << "' but declares no shiftfn: put the "
                              << "group's shiftfn block on it." << std::endl;
          return -1;
        }
        continue;
      }
      if (seg.OwnsShiftFn()) {
        configure.outStream << "ERROR: Data segment #" << seg.GetSegmentKey() << " declares a shiftfn, but it is not "
                            << "the first segment of shift-function group '" << label << "' (segment #"
                            << segments_[it->second].GetSegmentKey() << " is): one function per group." << std::endl;
        return -1;
      }
      fnGroupLeader[pos] = it->second;
    }
  }
  int fnStart = 0;
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    ESegment &owner = segments_[pos];
    if (!owner.OwnsShiftFn()) continue;
    const int key = owner.GetSegmentKey();
    const std::string expression = owner.GetDeclaredShiftFn();
    const int n = (int)owner.GetNominalShiftFnValues().size();
    if (owner.IsShiftGroupFollower()) {
      configure.outStream << "ERROR: Data segment #" << key << " declares a shiftfn, but it is not the first segment of "
                          << "shift group '" << owner.GetShiftGroup() << "' (segment #"
                          << segments_[owner.GetShiftLeaderPos()].GetSegmentKey()
                          << " is): put the block on that segment; it applies to the whole group." << std::endl;
      return -1;
    }
    std::vector<int> members{pos};
    for (int f = 0; f < (int)segments_.size(); f++)
      if (segments_[f].GetShiftLeaderPos() == pos || fnGroupLeader[f] == pos) members.push_back(f);
    for (int m : members)
      if (segments_[m].HasSqrtShiftBlock()) {
        configure.outStream << "ERROR: Data segment #" << segments_[m].GetSegmentKey() << " has a sqrtshift block and "
                            << "shift function \"" << expression << "\" (segment #" << key
                            << "): two descriptions of one energy calibration.  Write the sqrt(E) term into the function, "
                            << "e.g. a" << n << "*sqrt(x), and drop the sqrtshift block." << std::endl;
        return -1;
      }
    const std::set<int> refs = ShiftFnParamRefs(expression);
    for (int j : refs)
      if (j >= n) {
        configure.outStream << "ERROR: Data segment #" << key << ": shift function \"" << expression << "\" uses a" << j
                            << " but declares " << n << " parameter" << (n == 1 ? "" : "s") << " (a0..a" << n - 1 << ")."
                            << std::endl;
        return -1;
      }
    for (int j = 0; j < n; j++)
      if (!refs.count(j)) {
        configure.outStream << "ERROR: Data segment #" << key << ": shift function \"" << expression << "\" declares a"
                            << j << " but never uses it." << std::endl;
        return -1;
      }
    for (int m : members) segments_[m].InstallShiftFn(owner, pos, configure);
    owner.SetShiftFnParamStart(fnStart);
    fnStart += n;

    // Degeneracy check: rows = the members' distinct original lab energies (per member,
    // since members of a shift-function group have their own constants); columns =
    // d(dE)/d(term) for every free constant (one per constant owner) and every varied aj.
    std::vector<std::pair<int, double>> rows;  // constant owner position, energy
    for (int m : members) {
      const int constOwner = segments_[m].GetShiftLeaderPos() >= 0 ? segments_[m].GetShiftLeaderPos() : m;
      std::vector<double> e;
      for (int pid = 1; pid <= segments_[m].NumPoints(); pid++) {
        EPoint *pt = segments_[m].GetPoint(pid);
        e.push_back(pt->GetOriginalEnergy() > 0 ? pt->GetOriginalEnergy() : pt->GetLabEnergy());
      }
      std::sort(e.begin(), e.end());
      e.erase(std::unique(e.begin(), e.end()), e.end());
      for (double x : e) rows.push_back({constOwner, x});
    }
    std::vector<double> energies;
    for (const auto &r : rows) energies.push_back(r.second);
    std::sort(energies.begin(), energies.end());
    energies.erase(std::unique(energies.begin(), energies.end()), energies.end());
    std::vector<std::string> names;
    std::vector<std::vector<double>> cols;
    std::vector<int> constOwners;
    for (const auto &r : rows)
      if (std::find(constOwners.begin(), constOwners.end(), r.first) == constOwners.end()) constOwners.push_back(r.first);
    for (int c : constOwners) {
      if (!segments_[c].IsVaryEnergyShift()) continue;
      names.push_back(constOwners.size() == 1 ? std::string("a") : "a#" + std::to_string(segments_[c].GetSegmentKey()));
      std::vector<double> col(rows.size());
      for (size_t r = 0; r < rows.size(); r++) col[r] = (rows[r].first == c) ? 1.0 : 0.0;
      cols.push_back(col);
    }
    ESegment probe = owner;  // a scratch copy: the owner's own values stay untouched
    for (int j = 0; j < n; j++) {
      if (!owner.IsVaryShiftFn(j)) continue;
      const double v0 = owner.GetShiftFnValues()[j], h = 1.0e-6 * (std::fabs(v0) + 1.0);
      std::vector<double> c(rows.size());
      for (size_t r = 0; r < rows.size(); r++) {
        probe.SetShiftFnValue(j, v0 + h);
        const double fp = probe.ShiftFnTerm(rows[r].second);
        probe.SetShiftFnValue(j, v0 - h);
        c[r] = (fp - probe.ShiftFnTerm(rows[r].second)) / (2.0 * h);
      }
      probe.SetShiftFnValue(j, v0);
      double norm = 0.0;
      for (double v : c) norm += v * v;
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        configure.outStream << "ERROR: Data segment #" << key << ": in shift function \"" << expression << "\", a" << j
                            << " has no effect at this group's energies (" << energies.front() << "-" << energies.back()
                            << " MeV); check the expression." << std::endl;
        return -1;
      }
      names.push_back("a" + std::to_string(j));
      cols.push_back(c);
    }
    std::ostringstream note;
    note << "Note: shift function on segment #" << key;
    if (members.size() > 1) {
      const bool fnGroup = !owner.GetShiftFnGroup().empty();
      note << " (" << (fnGroup ? "shift-function group '" + owner.GetShiftFnGroup() : "group '" + owner.GetShiftGroup())
           << "', segments #" << key;
      for (size_t m = 1; m < members.size(); m++) note << " " << segments_[members[m]].GetSegmentKey();
      note << (fnGroup ? "; each keeps its own constant shift)" : ")");
    }
    note << ": dE = a + " << expression << ", x = lab energy (MeV);";
    for (int j = 0; j < n; j++) {
      const double err = owner.GetShiftFnErrors()[j];
      note << " a" << j << " " << owner.GetNominalShiftFnValues()[j]
           << (owner.IsVaryShiftFn(j) ? (err != 0.0 ? " free (penalty width " + std::to_string(err) + ")" : " free, unpenalised")
                                      : " fixed")
           << (j + 1 < n ? "," : ".");
    }
    const int m = (int)rows.size(), k = (int)cols.size();
    if (k >= 2) {
      // Condition number of the column-normalised matrix (scale-free).
      double cond = std::numeric_limits<double>::infinity();
      if (m >= k) {
        gsl_matrix *A = gsl_matrix_alloc(m, k);
        gsl_matrix *V = gsl_matrix_alloc(k, k);
        gsl_vector *S = gsl_vector_alloc(k), *work = gsl_vector_alloc(k);
        for (int c = 0; c < k; c++) {
          double norm = 0.0;
          for (int r = 0; r < m; r++) norm += cols[c][r] * cols[c][r];
          norm = std::sqrt(norm);
          for (int r = 0; r < m; r++) gsl_matrix_set(A, r, c, cols[c][r] / norm);
        }
        gsl_set_error_handler_off();
        if (gsl_linalg_SV_decomp(A, V, S, work) == 0 && gsl_vector_get(S, k - 1) > 0.0)
          cond = gsl_vector_get(S, 0) / gsl_vector_get(S, k - 1);
        gsl_matrix_free(A);
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_vector_free(work);
      }
      note << " Condition number of the free terms (";
      for (int c = 0; c < k; c++) note << names[c] << (c + 1 < k ? ", " : "");
      note << ") over " << m << " member energ" << (m == 1 ? "y" : "ies") << " " << energies.front() << "-" << energies.back()
           << " MeV: " << cond << ".";
      configure.outStream << note.str() << std::endl;
      if (!(cond < 1.0e8))
        configure.outStream << "WARNING: shift function on segment #" << key << ": its free terms cannot be told apart "
                            << "over this group's energies (condition number " << cond << "); fix one of them or "
                            << "constrain it with a penalty." << std::endl;
    } else {
      configure.outStream << note.str() << std::endl;
    }
  }
  return 0;
}

/*!
 * Writes each shift-group follower's parameter as the value it actually uses,
 * its group leader's, so that param.sav / param.fit read back unchanged.
 */

void EData::SyncShiftGroupParams(ROOT::Minuit2::MnUserParameters &p) const {
  const int off = energyShiftParamOffset_;
  if (off < 0) return;
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    const int lead = segments_[pos].GetShiftLeaderPos();
    if (lead < 0 || off + pos >= (int)p.Params().size() || off + lead >= (int)p.Params().size()) continue;
    p.SetValue(off + pos, p.Value(off + lead));
  }
}

/*!
 * Links the segments of each normalization group ("normgroup <label>" on the
 * segment line) to the group's first segment in file order: that segment keeps
 * its own norm (value, uncertainty, vary flag, penalty), and every other member
 * takes its value instead.  A follower whose line varies its norm keeps that
 * parameter, fixed, so the parameter layout does not depend on the groups.  A
 * note per group lists the members.
 */

void EData::ResolveNormGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    ESegment &seg = segments_[pos];
    const std::string label = seg.GetNormGroup();
    if (label.empty()) continue;
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    ESegment &leader = segments_[it->second];
    seg.SetNormGroupFollower(it->second, leader.GetNorm(), leader.GetNominalNorm());
    leader.SetHasNormGroupFollowers(true);
  }
  for (const auto &g : leaders) {
    const ESegment &leader = segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (const ESegment &seg : segments_)
      if (seg.GetNormLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.outStream << "Note: normalization group '" << g.first << "' has only segment #"
                          << leader.GetSegmentKey() << "; nothing is shared." << std::endl;
      continue;
    }
    configure.outStream << "Note: normalization group '" << g.first << "': segments #" << leader.GetSegmentKey()
                        << " +" << members.str() << " share one normalization (segment #" << leader.GetSegmentKey()
                        << "'s, " << (leader.IsVaryNorm() ? "free" : "fixed") << ")." << std::endl;
  }
}

int EData::NormOwner(int s) const {
  if (s < 1 || s > (int)segments_.size()) return s;
  const int lead = segments_[s - 1].GetNormLeaderPos();
  return (lead >= 0) ? lead + 1 : s;
}

std::shared_ptr<DataCovariance::File> EData::CovarianceFileFor(const std::string &path, const Config &configure) {
  const std::string key = DataCovariance::CanonicalPath(path);
  auto it = covFiles_.find(key);
  if (it != covFiles_.end()) return it->second;
  std::string err;
  std::shared_ptr<DataCovariance::File> f = DataCovariance::Read(path, err);
  if (!f) {
    configure.outStream << "ERROR: Covariance file: " << err << "." << std::endl;
    covarianceError_ = true;
    return nullptr;
  }
  configure.outStream << "Covariance file " << path << ": " << f->N << " rows over " << f->files.size()
                      << " data file(s)" << (f->rel.empty() ? std::string() : ", " + std::to_string(f->rel.size()) +
                                                                                   " relvector term(s)")
                      << "." << std::endl;
  covFiles_[key] = f;
  return f;
}

/*!
 * One group per covariance file: its segments' points in the file's row order, the
 * absolute part of their submatrix, the relvector fractions on their rows, and the
 * Cholesky factor of M = D^-1 C D^-1.  Without relvector terms the factor is final here;
 * with them it is checked here with the data in place of the model and frozen at the
 * first evaluation (DataCovariance::Group::EnsureFactored).
 */
int EData::ResolveCovarianceGroups(const Config &configure) {
  covGroups_.clear();
  std::map<std::string, int> groupOf;
  struct Item { int row, seg, pid; };
  std::vector<std::vector<Item>> items;
  std::vector<std::shared_ptr<DataCovariance::File>> files;
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (seg->GetCovarianceFile().empty()) continue;
    const std::string key = DataCovariance::CanonicalPath(seg->GetCovarianceFile());
    auto fit = covFiles_.find(key);
    if (fit == covFiles_.end()) return -1;  // Fill reported it
    auto git = groupOf.find(key);
    int g;
    if (git == groupOf.end()) {
      g = (int)items.size();
      groupOf[key] = g;
      items.emplace_back();
      files.push_back(fit->second);
      auto grp = std::make_shared<DataCovariance::Group>();
      grp->path = seg->GetCovarianceFile();
      covGroups_.push_back(grp);
    } else {
      g = git->second;
    }
    seg->SetCovarianceGroup(g);
    covGroups_[g]->numSegments++;
    const std::string canonicalData = DataCovariance::CanonicalPath(seg->GetDataFile());
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      const int row = fit->second->RowOf(canonicalData, seg->GetPoint(pid + 1)->GetFileRecord());
      if (row < 0) return -1;  // cannot happen: Fill checked the file and the record count
      items[g].push_back({row, s, pid});
    }
  }
  if (covGroups_.empty()) {
    covFiles_.clear();
    return 0;
  }
  const std::string modelFile = configure.outputdir + "covarianceModel.out";
  bool anyRelative = false;
  for (size_t g = 0; g < covGroups_.size(); g++) {
    DataCovariance::Group &grp = *covGroups_[g];
    const DataCovariance::File &f = *files[g];
    std::vector<Item> &it = items[g];
    std::stable_sort(it.begin(), it.end(), [](const Item &a, const Item &b) { return a.row < b.row; });
    for (size_t k = 1; k < it.size(); k++)
      if (it[k].row == it[k - 1].row) {
        configure.outStream << "ERROR: Row " << it[k].row + 1 << " of covariance file " << grp.path
                            << " is fitted twice, in segments #" << this->GetSegment(it[k - 1].seg)->GetSegmentKey()
                            << " and #" << this->GetSegment(it[k].seg)->GetSegmentKey()
                            << " (overlapping energy or angle windows on one data file)." << std::endl;
        return -1;
      }
    const int n = (int)it.size();
    grp.n = n;
    grp.N = f.N;
    grp.members.resize(n);
    grp.rows.resize(n);
    grp.sigma.resize(n);
    grp.data.resize(n);
    grp.label.resize(n);
    for (int a = 0; a < n; a++) {
      grp.members[a] = {it[a].seg, it[a].pid};
      grp.rows[a] = it[a].row;
      ESegment *seg = this->GetSegment(it[a].seg);
      EPoint *pt = seg->GetPoint(it[a].pid + 1);
      grp.sigma[a] = pt->GetLabCrossSectionError();
      grp.data[a] = pt->GetLabCrossSection();
      std::ostringstream lab;
      lab << seg->GetSegmentKey() << " " << pt->GetLabEnergy() << " " << pt->GetLabAngle();
      grp.label[a] = lab.str();
    }
    grp.modelFile = modelFile;
    // B (diagonal and blocks) factored by connected components; it must be positive
    // definite by itself on the fitted rows.
    const int bad = grp.SetUp(f);
    if (bad >= 0) {
      configure.outStream << "ERROR: Covariance file " << grp.path << ": its diagonal and block terms are not positive"
                          << " definite on the fitted rows (at file row " << grp.rows[bad] + 1
                          << "). Every fitted row needs a variance of its own from a diagonal or block term;"
                          << " vector and relvector terms alone are not enough." << std::endl;
      return -1;
    }
    if (grp.HasRelative())
      anyRelative = true;  // the low-rank part waits for the model of the first evaluation
    else
      grp.Factor(grp.data);
    configure.outStream << "Covariance " << grp.path << ": " << n << " of " << f.N << " rows fitted, in "
                        << grp.numSegments << " segment(s); " << grp.NumComponents()
                        << " diagonal/block component(s) (largest " << grp.LargestComponent() << " rows), "
                        << grp.absolute.size() << " vector term(s)"
                        << (grp.HasRelative() ? ", " + std::to_string(grp.rel.size()) +
                                                    " relvector term(s) on the model frozen at the first evaluation"
                                              : std::string())
                        << "." << std::endl;
    if (f.numSigmaMismatch > 0)
      configure.outStream << "WARNING: " << f.numSigmaMismatch << " point(s) of " << grp.path
                          << " have an error column more than 1 % from sqrt(C_ii) (worst "
                          << 100.0 * f.worstSigmaMismatch << " %; first: " << f.firstSigmaMismatch
                          << "); sqrt(C_ii) is used." << std::endl;
  }
  if (anyRelative) std::ofstream(modelFile.c_str(), std::ios::trunc);
  covFiles_.clear();
  return 0;
}

double EData::CovarianceChi2Correction(std::vector<double> *segChis) {
  double correction = 0.0;
  std::vector<double> r, z;
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    r.assign(g.n, 0.0);
    double diag = 0.0;
    for (int a = 0; a < g.n; a++) {
      ESegment *seg = this->GetSegment(g.members[a].first);
      EPoint *pt = seg->GetPoint(g.members[a].second + 1);
      const double norm = seg->GetNorm();
      const double err = pt->GetCMCrossSectionError() * norm;
      r[a] = err != 0.0 ? (pt->GetFitCrossSection() - pt->GetCMCrossSection() * norm) / err : 0.0;
      diag += r[a] * r[a];
    }
    g.EnsureFactored(r);
    z = r;
    g.Whiten(z.data());
    double q = 0.0;
    for (int a = 0; a < g.n; a++) q += z[a] * z[a];
    g.WhitenTranspose(z.data());
    correction += q - diag;
    for (int a = 0; a < g.n; a++) {
      ESegment *seg = this->GetSegment(g.members[a].first);
      seg->SetCovarianceSeed(g.members[a].second, z[a]);
    }
    // Per-segment chi2: the exact partition sum_i r_i (M^-1 r)_i of the group's chi2.
    std::map<int, double> part;
    for (int a = 0; a < g.n; a++) part[g.members[a].first] += r[a] * z[a];
    for (const auto &ps : part) {
      this->GetSegment(ps.first)->SetSegmentChiSquared(ps.second);
      if (segChis && ps.first < (int)segChis->size()) (*segChis)[ps.first] = ps.second;
    }
  }
  return correction;
}

void EData::WhitenResidualRows(std::vector<double> &residuals, std::vector<double> *jacobian, int nCols) {
  if (covGroups_.empty()) return;
  std::vector<int> rowStart(this->NumSegments() + 2, 0);
  int row = 0;
  for (int s = 1; s <= this->NumSegments(); s++) {
    rowStart[s] = row;
    row += this->GetSegment(s)->NumPoints();
  }
  if (row > (int)residuals.size()) return;  // not the data-row layout
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    std::vector<int> rows(g.n);
    for (int a = 0; a < g.n; a++) rows[a] = rowStart[g.members[a].first] + g.members[a].second;
    std::vector<double> y(g.n);
    for (int a = 0; a < g.n; a++) y[a] = residuals[rows[a]];
    g.EnsureFactored(y);  // the first evaluation of a run with relvector terms freezes the model here
    g.Whiten(y.data());
    for (int a = 0; a < g.n; a++) residuals[rows[a]] = y[a];
    if (!jacobian || nCols <= 0) continue;
    std::vector<double> &J = *jacobian;
#pragma omp parallel for schedule(dynamic, 4)
    for (int c = 0; c < nCols; c++) {
      std::vector<double> col(g.n);
      bool any = false;
      for (int a = 0; a < g.n; a++) any |= (col[a] = J[(size_t)rows[a] * nCols + c]) != 0.0;
      if (!any) continue;
      g.Whiten(col.data());
      for (int a = 0; a < g.n; a++) J[(size_t)rows[a] * nCols + c] = col[a];
    }
  }
}

/*!
 * covariance.out: per covariance group, the fitted rows, the chi2 r^T M^-1 r of the
 * current fit values and, for comparison, the sum of r_i^2 (the chi2 the same points
 * would have without their correlations).
 */
void EData::WriteCovarianceSummary(const Config &configure) {
  std::ofstream out((configure.outputdir + "covariance.out").c_str());
  if (!out) return;
  out << "# file, fitted rows, rows, segments, chi2 (r^T M^-1 r), sum of r_i^2 (correlations ignored), relvector terms"
      << std::endl;
  out.precision(12);
  out << std::scientific;
  std::vector<double> r;
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    r.assign(g.n, 0.0);
    double diag = 0.0;
    for (int a = 0; a < g.n; a++) {
      ESegment *seg = this->GetSegment(g.members[a].first);
      EPoint *pt = seg->GetPoint(g.members[a].second + 1);
      const double norm = seg->GetNorm();
      const double err = pt->GetCMCrossSectionError() * norm;
      r[a] = err != 0.0 ? (pt->GetFitCrossSection() - pt->GetCMCrossSection() * norm) / err : 0.0;
      diag += r[a] * r[a];
    }
    g.EnsureFactored(r);
    g.Whiten(r.data());
    double q = 0.0;
    for (int a = 0; a < g.n; a++) q += r[a] * r[a];
    out << g.path << ", " << g.n << ", " << g.N << ", " << g.numSegments << ", " << q << ", " << diag << ", "
        << g.rel.size() << std::endl;
  }
}

void EData::RefreezeCovariance() {
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) gp->Refreeze();
}

void EData::PrepareCovarianceSeeds(CNuc *theCNuc, const Config &configure) {
  if (covGroups_.empty()) return;
  std::vector<std::pair<ESegment *, int>> work;
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (seg->GetCovarianceGroup() < 0) continue;
    for (int pid = 0; pid < seg->NumPoints(); pid++) work.push_back({seg, pid});
  }
  // Each point on its own, without point mapping (as the gradient's own point loops): a
  // mapped point would otherwise read a host value that may not be calculated yet.
#pragma omp parallel for schedule(dynamic) if (!omp_in_parallel())
  for (int k = 0; k < (int)work.size(); k++) {
    ScopedPointMappingOff noPointMapping;
    ESegment *seg = work[k].first;
    const double th = seg->CalculateTheoreticalCrossSection(work[k].second, theCNuc, configure, this);
    seg->GetPoint(work[k].second + 1)->SetFitCrossSection(th);
  }
  this->CovarianceChi2Correction();
}


/*!
 * Writes each norm-group follower's (fixed) norm parameter as the value it
 * actually uses, its group leader's, so that param.sav / param.fit read back
 * unchanged.
 */

void EData::SyncNormGroupParams(ROOT::Minuit2::MnUserParameters &p) const {
  const int off = GetNormParamOffset();
  if (off < 0) return;
  std::vector<int> index(segments_.size(), -1);  // norm parameter of each segment, -1 = none
  int i = off;
  for (int pos = 0; pos < (int)segments_.size(); pos++)
    if (segments_[pos].IsVaryNorm()) index[pos] = i++;
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    const int lead = segments_[pos].GetNormLeaderPos();
    if (lead < 0 || index[pos] < 0 || index[pos] >= (int)p.Params().size()) continue;
    const double value = (index[lead] >= 0 && index[lead] < (int)p.Params().size()) ? p.Value(index[lead])
                                                                                     : segments_[lead].GetNorm();
    p.SetValue(index[pos], value);
  }
}

/*!
 * Fills the Minuit parameter array from initial values in the EData object.
 */

void EData::FillMnParams(ROOT::Minuit2::MnUserParameters &p) {
  SetNormParamOffset(p.Params().size());
  char varname[50];
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryNorm()) {
      snprintf(varname, sizeof(varname), "segment_%d_norm", segment->GetSegmentKey());
      p.Add(varname, segment->GetNorm(), segment->GetNorm() * 0.05);
      p.SetLowerLimit(varname, 0.0);
      if (segment->IsNormGroupFollower()) p.Fix(varname);  // the group's norm is its leader's parameter

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

  // Shift-function parameters: only segments that declare a "shiftfn" block have
  // any, so a model without one gets no new parameter at all.
  SetEnergyShiftFnParamOffset(p.Params().size());
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (!segment->OwnsShiftFn() || !segment->HasShiftFn()) continue;
    for (int j = 0; j < segment->NumShiftFnParams(); j++) {
      snprintf(varname, sizeof(varname), "segment_%d_energy_shift_fn%d", segment->GetSegmentKey(), j);
      const double err = segment->GetShiftFnErrors()[j];
      p.Add(varname, segment->GetShiftFnValues()[j], (err > 0.0) ? err * 0.01 : 0.0005);
      if (!segment->IsVaryShiftFn(j)) p.Fix(varname);
    }
  }

  // Electron-screening potentials (eV): one per target effect with a "screening" block, none otherwise.
  screeningParamOffset_ = p.Params().size();
  for (int k = 1; k <= this->NumScreeningParams(); k++) {
    const TargetEffect *te = this->GetTargetEffect(this->GetScreeningEffect(k));
    snprintf(varname, sizeof(varname), "screening_%d", this->GetScreeningEffect(k));
    const double ue = te->GetScreeningUe();
    p.Add(varname, ue, std::max(10.0, 0.05 * std::fabs(ue)));
    p.SetLowerLimit(varname, 0.0);
    if (!te->IsScreeningVaried()) p.Fix(varname);
  }
}

void EData::ApplyScreening() {
  if (screeningEffects_.empty()) return;
  std::vector<int> paramOf(this->NumTargetEffects() + 1, 0);
  for (int k = 1; k <= this->NumScreeningParams(); k++) paramOf[this->GetScreeningEffect(k)] = k;
  auto apply = [&](ESegment &seg) {
    const int te = seg.GetTargetEffectNum();
    if (te <= 0 || te >= (int)paramOf.size() || paramOf[te] == 0) return;
    const double ueMeV = this->GetTargetEffect(te)->GetScreeningUe() * 1.0e-6;
    for (int i = 1; i <= seg.NumPoints(); i++) {
      EPoint *point = seg.GetPoint(i);
      point->SetScreening(ueMeV, paramOf[te]);
      for (int j = 1; j <= point->NumSubPoints(); j++) point->GetSubPoint(j)->SetScreening(ueMeV, paramOf[te]);
    }
  };
  for (ESegment &seg : segments_) apply(seg);
  for (ESegment &seg : componentSegments_) apply(seg);
}

void EData::FillScreeningFromParams(const vector_r &p) {
  if (screeningEffects_.empty() || screeningParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumScreeningParams(); k++) {
    const int idx = screeningParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the screening block
    this->GetTargetEffect(this->GetScreeningEffect(k))->SetScreeningUe(p[idx]);
  }
  this->ApplyScreening();
}

double EData::ScreeningPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumScreeningParams(); k++)
    chi2 += const_cast<EData *>(this)->GetTargetEffect(this->GetScreeningEffect(k))->ScreeningPenalty();
  return chi2;
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
    if (segment->IsVaryNorm()) {
      segment->SetNorm(p[i]);
      i++;
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
  // Norm-group followers take their leader's norm (a follower's own parameter is fixed).
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++)
    if (segment->IsNormGroupFollower()) segment->SetNorm(GetSegments()[segment->GetNormLeaderPos()].GetNorm());
}

/*!
 * Fills the Energy Shifts from the Minuit parameter array.
 */

void EData::FillEnergyShiftsFromParams(const vector_r &p, EData *data, CNuc *theCNuc, const Config *configure) {
  // Screening first: sub-points rebuilt below (ReplaceSubPointGrids) inherit their parent's U_e.
  this->FillScreeningFromParams(p);
  int i = GetEnergyShiftParamOffset();
  int j = GetEnergyShiftSqrtParamOffset();
  int k = 0;
  bool anyEnergyShifted = false;

  if (data) {
    // Re-energising a point (ESegment::UpdatePointEnergyWithShift: Coulomb
    // functions of the point and its sub-points) touches only that point, so
    // the shift values are set here in the original order and the points of
    // the segments that moved are re-energised afterwards in parallel, each
    // segment once, in its final state.  A line search moves every shift at once; serially that cost
    // ~195 s per evaluation on the 11B+a mixed grid against 0.5 s for the
    // points.  Inside an enclosing parallel region (Minuit's numerical
    // gradient) the loop stays serial, as before.
    std::vector<ESegment *> toUpdate;
    auto queueUpdate = [&toUpdate](ESegment *s) {
      if (std::find(toUpdate.begin(), toUpdate.end(), s) == toUpdate.end()) toUpdate.push_back(s);
    };
    // Segments (and their components) whose points carry an energy shift, moved or not.
    std::vector<ESegment *> shiftActiveSegments;
    for (ESegmentIterator segment = data->GetSegments().begin(); segment < data->GetSegments().end(); segment++) {
      k++;

      // Both terms are present for every segment: the constant shift at p[i]
      // and the sqrt(E) coefficient at p[j].  Either moving the points means
      // the segment must be re-shifted (the mapping E -> E + a + b*sqrt(E) is
      // applied from the original energies, so one update covers both).
      // A shift-group follower reads its group leader's constant shift.
      const int lead = segment->GetShiftLeaderPos();
      const double shiftConst = (lead >= 0) ? p[GetEnergyShiftParamOffset() + lead] : p[i];
      const double shiftSqrt = (j >= 0 && j < (int)p.size()) ? p[j] : 0.0;
      // A shift function's parameters sit in the shift-function block, at the
      // owning segment's (this one or its group leader) start.
      std::vector<double> fnValues;
      if (segment->HasShiftFn()) {
        const int start = data->GetSegment(segment->GetShiftFnSourcePos() + 1)->GetShiftFnParamStart();
        for (int q = 0; q < segment->NumShiftFnParams(); q++) {
          const int idx = GetEnergyShiftFnParamOffset() + start + q;
          fnValues.push_back((start >= 0 && idx < (int)p.size()) ? p[idx] : segment->GetShiftFnValues()[q]);
        }
      }
      const bool shiftActive = segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 ||
                               segment->IsVaryEnergyShiftSqrt() || shiftSqrt != 0.0 || segment->HasShiftFn();
      if (shiftActive) {
        shiftActiveSegments.push_back(&*segment);
        if (segment->HasComponents())
          for (ESegment *componentSegment : segment->GetComponentSegments())
            if (componentSegment) shiftActiveSegments.push_back(componentSegment);
        // Check if energy is the same, if so, continue
        if (segment->GetLastEnergyShift() == shiftConst && segment->GetLastEnergyShiftSqrt() == shiftSqrt &&
            (!segment->HasShiftFn() || segment->GetLastShiftFnValues() == fnValues)) {
          i++;
          j++;
          continue;
        }

        segment->SetEnergyShift(shiftConst);
        segment->SetLastEnergyShift(shiftConst);
        segment->SetEnergyShiftSqrt(shiftSqrt);
        segment->SetLastEnergyShiftSqrt(shiftSqrt);
        for (int q = 0; q < (int)fnValues.size(); q++) segment->SetShiftFnValue(q, fnValues[q]);
        if (segment->HasShiftFn()) segment->SetLastShiftFnValues(fnValues);
        queueUpdate(&*segment);
        anyEnergyShifted = true;

        // Apply the same energy shift to ALL component segments of this master segment
        if (segment->HasComponents()) {
          for (ESegment *componentSegment : segment->GetComponentSegments()) {
            if (componentSegment) {
              componentSegment->SetEnergyShift(segment->GetEnergyShift());
              componentSegment->SetEnergyShiftSqrt(segment->GetEnergyShiftSqrt());
              componentSegment->CopyShiftFnFrom(*segment);
              queueUpdate(componentSegment);
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
                  componentSegment->GetLastEnergyShiftSqrt() == shiftSqrt &&
                  (!segment->HasShiftFn() || componentSegment->GetLastShiftFnValues() == fnValues)) {
                continue;
              }

              componentSegment->SetEnergyShift(shiftConst);
              componentSegment->SetLastEnergyShift(shiftConst);
              componentSegment->SetEnergyShiftSqrt(shiftSqrt);
              componentSegment->SetLastEnergyShiftSqrt(shiftSqrt);
              if (segment->HasShiftFn()) {
                componentSegment->CopyShiftFnFrom(*segment);
                componentSegment->SetLastShiftFnValues(fnValues);
              }
              queueUpdate(componentSegment);
              anyEnergyShifted = true;
            }
          }
        }
      }

      // Always increment i and j since both shift parameters are present for ALL segments
      i++;
      j++;
    }

    // The sub-point grids, for the levels and shifts just set.  Sub-points of
    // shifted segments are moved and energised by the step below.
    if (theCNuc && configure) {
      std::vector<ESegment *> reshift;
      if (data->ReplaceSubPointGrids(theCNuc, *configure, shiftActiveSegments, reshift) > 0) anyEnergyShifted = true;
      for (ESegment *s : reshift) queueUpdate(s);
    }

    if (!toUpdate.empty()) {
      // One work item per point (with its sub-points), largest first, so a
      // segment with thousands of sub-points per point does not serialize.
      // A point mapped onto another (EData::MapData) copies its values from
      // that host, so it is re-energised with it, on the same thread, after
      // it (ESegment::UpdatePointGroupEnergyWithShift); it gets no item of its
      // own.  A mapped point costs little (a copy, its Coulomb amplitudes).
      struct Item { double cost; ESegment *seg; int i; };
      std::vector<Item> work;
      for (ESegment *s : toUpdate)
        for (int i = 0; i < s->NumPoints(); i++) {
          if (s->ShiftHost(i)) continue;
          EPoint *pt = s->GetPoint(i + 1);
          const double nsp = pt ? pt->NumSubPoints() : 0;
          const double nMapped = pt ? pt->NumLocalMappedPoints() : 0;
          work.push_back({(1.0 + nsp) * (1.0 + 0.1 * nMapped), s, i});
        }
      std::stable_sort(work.begin(), work.end(), [](const Item &x, const Item &y) { return x.cost > y.cost; });
#ifdef _OPENMP
      const bool serial = omp_in_parallel();
#else
      const bool serial = true;
#endif
      // No parallel construct at all when already inside a parallel region: an
      // inactive nested team renumbers its thread 0, and per-thread state
      // indexed by thread number would then be shared between the outer threads.
      if (serial) {
        for (size_t q = 0; q < work.size(); q++) work[q].seg->UpdatePointGroupEnergyWithShift(work[q].i, theCNuc, configure);
      } else {
#pragma omp parallel for schedule(dynamic, 1)
        for (int q = 0; q < (int)work.size(); q++) work[q].seg->UpdatePointGroupEnergyWithShift(work[q].i, theCNuc, configure);
      }
    }

    // Critical fix: Rebuild energy-based mapping system after energy shifts
    if (anyEnergyShifted) {
      data->MapData(theCNuc);
    }

    return;
  }
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    const int lead = segment->GetShiftLeaderPos();
    const double shiftConst = (lead >= 0) ? p[GetEnergyShiftParamOffset() + lead] : p[i];
    if (segment->HasShiftFn()) {
      const int start = GetSegment(segment->GetShiftFnSourcePos() + 1)->GetShiftFnParamStart();
      for (int q = 0; q < segment->NumShiftFnParams(); q++) {
        const int idx = GetEnergyShiftFnParamOffset() + start + q;
        if (start >= 0 && idx < (int)p.size()) segment->SetShiftFnValue(q, p[idx]);
      }
    }
    // Always apply energy shift since parameters are always present for all segments
    if (segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 || segment->HasShiftFn()) {
      segment->SetEnergyShift(shiftConst);
      segment->UpdatePointEnergiesWithShift();
      anyEnergyShifted = true;

      // CRITICAL FIX: Apply the same energy shift to ALL component segments of this master segment
      if (segment->HasComponents()) {
        for (ESegment *componentSegment : segment->GetComponentSegments()) {
          if (componentSegment) {
            componentSegment->SetEnergyShift(segment->GetEnergyShift());
            componentSegment->CopyShiftFnFrom(*segment);
            componentSegment->UpdatePointEnergiesWithShift();
          }
        }
      }
    }

    // CRITICAL FIX: Apply the same energy shift to all component segments in total capture segments
    if (segment->IsTotalCapture() && segment->HasComponents()) {
      if (segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 || segment->HasShiftFn()) {
        // Apply energy shift to all component segments
        const std::vector<ESegment *> &componentSegments = segment->GetComponentSegments();
        for (ESegment *componentSegment : componentSegments) {
          if (componentSegment) {
            componentSegment->SetEnergyShift(shiftConst);
            componentSegment->CopyShiftFnFrom(*segment);
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
    this->MapData(theCNuc);
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
  dataCopy->energyShiftFnParamOffset_ = this->energyShiftFnParamOffset_;
  dataCopy->screeningParamOffset_ = this->screeningParamOffset_;
  dataCopy->screeningEffects_ = this->screeningEffects_;
  dataCopy->isFit_ = this->isFit_;
  dataCopy->isErrorAnalysis_ = this->isErrorAnalysis_;
  dataCopy->targetEffects_ = this->targetEffects_;
  dataCopy->segments_ = this->segments_;
  dataCopy->componentSegments_ = this->componentSegments_;
  dataCopy->covGroups_ = this->covGroups_;  // shared: one factor (and one frozen model) for all copies

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
  // Component points mapped onto their base segment's points (MapComponentSegments), in the order MapData added
  // them: components in turn.
  for (size_t i = 0; i < dataCopy->segments_.size(); i++)
    for (ESegment *component : dataCopy->segments_[i].GetComponentSegments())
      for (int j = 1; component && j <= component->NumPoints(); j++) {
        EPoint *cp = component->GetPoint(j);
        if (!cp->IsMapped()) continue;
        EnergyMap pointMap = cp->GetMap();
        dataCopy->GetSegment(pointMap.segment)->GetPoint(pointMap.point)->AddLocalMappedPoint(cp);
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
