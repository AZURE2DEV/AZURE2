#include "InputFiles.h"
#include "AZUREOutput.h"
#include "Profile.h"
#include "CNuc.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "EData.h"
#include "ModelPass.h"
#include "ECAmplitudeCache.h"
#include "ECTable.h"
#include "AdaptiveIntegrationGrid.h"
#include "Minuit2/MnUserParameters.h"
#include "AZUREParams.h"
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
#include <cmath>
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
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

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

/*!
 * The UPOS coefficient (CNuc, Brune and deBoer 2020 Eq. 22) is written for a secondary gamma ray from the HEAVY
 * partner of the exit pair: its final-state spin is J2.  With the light particle excited (column 32 of the
 * <levels> line, 2026-11-05) the gamma ray would come from the light particle, which that formula does not
 * describe -- refuse instead of computing the heavy partner's correlation silently.
 */
static bool UPOSAllowed(CNuc *theCNuc, const ESegment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !theCNuc->IsPairKey(segment.GetExitKey())) return true;
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetExitKey()));
  if (exitPair->GetLightExE() == 0.) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << key << " (UPOS): exit pair "
                      << segment.GetExitKey() << " has an excited light particle (column 32, " << exitPair->GetLightExE()
                      << " MeV); the secondary-gamma correlation is implemented for the decay of the heavy partner only."
                      << std::endl;
  return false;
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

// A beta-delayed spectrum (entrance pair of type 20) is the energy distribution of one decay pair: its energies are
// converted with that pair's masses (ESegment::ConvertPointEnergyToCM), so a segment without one (exit -1) is refused.
static bool BetaDelayedAllowed(CNuc *theCNuc, ESegment &segment, int key, bool isTest, const Config &configure) {
  if (!theCNuc->IsPairKey(segment.GetEntranceKey()) || segment.GetExitKey() != -1) return true;
  if (theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetEntranceKey()))->GetPType() != 20) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment" : "segment") << " #" << key << ": a beta-delayed "
                      << "spectrum (entrance pair " << segment.GetEntranceKey() << " of type 20) needs its decay pair as"
                      << " the exit pair, not -1 (total capture)." << std::endl;
  return false;
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

int EData::Fill(const Config &configure, CNuc *theCNuc, const azure2::experiment::DataSpecs &specs) {
  const int numTotalSegments_ = (int)specs.datasets.items.size();
  for (const auto &item : specs.datasets.items) {
    {
      const azure2::experiment::DatasetSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        ESegment NewSegment(segment);

        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          entrancePair->SetEntrance();
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!CoefficientsAllowed(theCNuc, NewSegment, numTotalSegments, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!BetaDelayedAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!UPOSAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
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
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightBySegment, numTotalSegments));
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

                // One full segment copy per component (spec: entrance, exit, optional fixed angle, scaling)
                for (const azure2::experiment::ComponentSpec &c : NewSegment.GetComponents()) {
                  if (!this->CheckComponentPairs(theCNuc, *filledSegment, c.entranceKey, c.exitKey, configure)) return -1;
                  ESegment *componentSegment =
                      c.hasAngle ? this->CreateComponentSegment(*filledSegment, c.entranceKey, c.exitKey, c.angle)
                                 : this->CreateComponentSegment(*filledSegment, c.entranceKey, c.exitKey);
                  if (componentSegment) {
                    componentSegment->SetComponentScaling(c.scaling);
                    filledSegment->AddComponentSegment(componentSegment);
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
              if (!this->FinishTransformedSegment(theCNuc, this->GetSegment(thisSegmentNum), numTotalSegments, false, configure))
                return -1;
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

  if (specs.datasets.failed) {
    if (!specs.datasets.error.empty()) configure.outStream << specs.datasets.error << std::endl;
    return -1;
  }

  if (this->ResolveShiftGroups(configure) == -1) return -1;
  this->ResolveNormGroups(configure);
  if (this->ResolveThicknessGroups(configure) == -1) return -1;
  if (this->ResolveCovarianceGroups(configure) == -1) return -1;

  if (ApplyPriorCentres(configure, specs) != 0) return -1;

  if (CheckThmWeights(configure, configure.thm.weightBySegment, "weight", "<segmentsData>",
                      numTotalSegments_, this) != 0)
    return -1;
  if (BuildThmGroups(configure, theCNuc, numTotalSegments_) != 0) return -1;

  if (this->NumSegments() > 0) {
    if (this->SetUpTargetEffects(configure, theCNuc, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(theCNuc, configure)) return -1;
    this->MapData(theCNuc);
  }

  return 0;
}

/*!
 * Moves the nominal normalization / energy shift of the segments whose specs give a prior centre (after the groups
 * are resolved, so that a group's follower keeps its own centre).  A centre for a data set that has no segment (not
 * active, or skipped at load) is reported and ignored.  A refusal of the reader (specs.priorCentreError) is reported
 * after the centres before it are applied, and fails the run.  Segment keys are unique in segments_, so at most one
 * segment matches a spec.
 */

int EData::ApplyPriorCentres(const Config &configure, const azure2::experiment::DataSpecs &specs) {
  for (const auto &item : specs.datasets.items) {
    const azure2::experiment::DatasetSpec &d = item.spec;
    if (!d.hasNormCentre && !d.hasShiftCentre) continue;
    ESegment *segment = nullptr;
    for (int s = 1; s <= this->NumSegments(); s++)
      if (this->GetSegment(s)->GetSegmentKey() == d.key) segment = this->GetSegment(s);
    if (d.hasNormCentre) {
      if (segment)
        segment->SetNominalNorm(d.normCentre);
      else
        configure.outStream << "WARNING: <parameterSettings> row \"" << d.normCentreSource
                            << "\" names no active data segment; ignored." << std::endl;
    }
    if (d.hasShiftCentre) {
      if (segment)
        segment->SetNominalEnergyShift(d.shiftCentre);
      else
        configure.outStream << "WARNING: <parameterSettings> row \"" << d.shiftCentreSource
                            << "\" names no active data segment; ignored." << std::endl;
    }
  }
  if (!specs.priorCentreError.empty()) {
    configure.outStream << specs.priorCentreError << std::endl;
    return -1;
  }
  return 0;
}

/*!
 * If the AZURE calculation is not data driven, this function is called in place
 * of the EData::Fill function to create points at specified energies and angles.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int EData::MakePoints(const Config &configure, CNuc *theCNuc, const azure2::experiment::ExtrapolationSpecs &specs) {
  const int numTotalSegments_ = (int)specs.extrapolations.items.size();
  for (const auto &item : specs.extrapolations.items) {
    {
      const azure2::experiment::ExtrapolationSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        ESegment NewSegment(segment);
        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          PPair *segEntrance = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!BetaDelayedAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!UPOSAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
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
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightByTestSegment, numTotalSegments));
            this->AddSegment(NewSegment);
            ESegment *theSegment = this->GetSegment(this->NumSegments());

            theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()))->SetEntrance();
            PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()));
            PPair *exitPair = theCNuc->GetPairFromKey(theSegment->GetExitKey());  // nullptr for total capture (-1)
            for (double angle : segment.angles) {
              for (double energy : segment.energies) {
                EPoint NewPoint(angle, energy, theSegment);
                theSegment->AddPoint(NewPoint);
                EPoint *thePoint = theSegment->GetPoint(theSegment->NumPoints());
                thePoint->SetParentData(this);
                // A ratio is dimensionless: no cross-section conversion (as ESegment::Fill for data segments;
                // the stored factor would otherwise divide the ratio in the lab output, kinflag.dat = 1).
                if (!theSegment->ConvertPointToCM(thePoint, entrancePair, exitPair, configure,
                                                  !theSegment->IsDimensionlessRatio(), ESegment::BeyondRange::Stop)) {
                  configure.outStream << "ERROR: " << "test segment #" << theSegment->GetSegmentKey() << ": lab angle " << thePoint->GetLabAngle() << " deg at E_lab = " << thePoint->GetLabEnergy()
                                      << " MeV: " << thePoint->GetLabAngleWhy() << "." << std::endl;
                  return -1;
                }
              }
            }

            // Handle advanced segments (sum/ratio of components) - AFTER points are added
            if (NewSegment.IsAdvanced()) {
              int operation = NewSegment.GetOperationType();
              theSegment->SetOperationType((OperationType)operation);

              // One full segment copy per component (spec: entrance, exit, optional fixed angle, scaling)
              for (const azure2::experiment::ComponentSpec &c : NewSegment.GetComponents()) {
                if (!this->CheckComponentPairs(theCNuc, *theSegment, c.entranceKey, c.exitKey, configure)) return -1;
                ESegment *componentSegment =
                    c.hasAngle ? this->CreateComponentSegment(*theSegment, c.entranceKey, c.exitKey, c.angle)
                               : this->CreateComponentSegment(*theSegment, c.entranceKey, c.exitKey);
                if (componentSegment) {
                  componentSegment->SetComponentScaling(c.scaling);
                  theSegment->AddComponentSegment(componentSegment);
                }
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
              if (!this->FinishTransformedSegment(theCNuc, this->GetSegment(thisSegmentNum), numTotalSegments, true, configure))
                return -1;
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

  if (specs.extrapolations.failed) {
    if (!specs.extrapolations.error.empty()) configure.outStream << specs.extrapolations.error << std::endl;
    return -1;
  }

  if (CheckThmWeights(configure, configure.thm.weightByTestSegment, "weightTest", "<segmentsTest>",
                      numTotalSegments_, this) != 0)
    return -1;

  if (this->NumSegments() > 0) {
    if (this->SetUpTargetEffects(configure, theCNuc, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(theCNuc, configure)) return -1;
    this->MapData(theCNuc);
  }

  return 0;
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
 * This function counts one more evaluation of the objective and returns the count before it. The read and the
 * increment are one atomic operation, so evaluations made concurrently from the threads of Minuit2's numerical
 * gradient each receive a distinct number. Relaxed order suffices: the counter orders no other memory.
 */

int EData::Iterate() {
  return iterations_.fetch_add(1, std::memory_order_relaxed);
}

/*!
 * This function is identical in role to the EPoint::Initialize function, except that it initializes
 * and entire EData object instead of a single EPoint object.
 */

int EData::Initialize(CNuc *compound, const Config &configure) {
  // Labels of the compound nuclei, for the sample decomposition files.
  nucleusLabels_.clear();
  for (int k = 1; k <= compound->NumNuclei(); k++) nucleusLabels_.push_back(compound->NucleusLabel(k));
  numNuclei_ = compound->NumNuclei();
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
  // Every point of this stage has its rows now; nothing reads the channel table until the energies change, and the
  // table (about 64 bytes a row) is not kept beside the rows it filled (docs/upgrade/channel_table.md section 2).
  channelTable_.reset();
  // Where the result record of every point and sub-point lies (Results.h); laid out again when the points change
  // (BeginEvaluation).
  this->LayoutResults(nullptr);

  return 0;
}

/*!
 * Adds a segment to the ESegment vector.
 */

void EData::AddSegment(ESegment segment) {
  segments_.push_back(segment);
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
 * Calls EPoint::CalcEDependentValues for each point in the entire EData object.
 */

/*!
 * The channel table of the given points and their sub-points, built before the threads that fill them start
 * (docs/upgrade/channel_table.md).  Rows the current table holds are copied from it.
 */
std::shared_ptr<const ChannelTable> EData::BuildChannelTable(CNuc *theCNuc, const Config &configure,
                                                             const std::vector<EPoint *> &points) {
  // Each point's keys are made unique before they are added: its slots repeat (pair, l, E_c) once per channel spin
  // and J^pi block, and the full list was the peak memory of the build (natural_neon: 11 M keys).
  std::vector<ChannelTable::Key> keys, own;
  std::size_t requested = 0;
  for (EPoint *point : points) {
    own.clear();
    point->ChannelKeys(theCNuc, own);
    requested += own.size();
    std::sort(own.begin(), own.end());
    own.erase(std::unique(own.begin(), own.end()), own.end());
    keys.insert(keys.end(), own.begin(), own.end());
  }
  std::shared_ptr<const ChannelTable> table = ChannelTable::Build(theCNuc, configure, std::move(keys), channelTable_.get());
  if (std::getenv("AZURE_PROFILE"))
    configure.outStream << "Channel table: " << requested << " rows asked for, " << table->Size() << " distinct, "
                        << table->Computed() << " computed" << std::endl;
  return table;
}

int EData::CalcEDependentValues(CNuc *theCNuc, const Config &configure) {
  {
    std::vector<EPoint *> points;
    for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++)
      for (int i = 1; i <= segment->NumPoints(); i++)
        if (!segment->GetPoint(i)->IsMapped()) points.push_back(segment->GetPoint(i));
    try {
      channelTable_ = BuildChannelTable(theCNuc, configure, points);
    } catch (std::exception &e) {
      configure.outStream << "Channel functions: " << e.what() << std::endl;
      return -1;
    }
  }
  ChannelRowSource rows;
  rows.tables.push_back(channelTable_.get());
  const ChannelRowSource *table = &rows;
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    bool localStop = false;
#pragma omp parallel for shared(localStop, configure)
    for (int i = 1; i <= segment->NumPoints(); i++) {
      if (configure.stopFlag || localStop) continue;
      EPoint *point = segment->GetPoint(i);
      if (!(point->IsMapped())) {
        try {
          point->CalcEDependentValues(theCNuc, configure, true, true, table);
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
        // The base segment's sequence, with the component's own pairs; a component's value is always converted.
        if (!componentSegment.ConvertPointToCM(point, entrancePair, exitPair, configure, true, ESegment::BeyondRange::Stop)) {
          configure.outStream << "ERROR: " << "component (" << entranceKey << " -> " << exitKey << ") of segment #" << componentSegment.GetSegmentKey() << ": lab angle " << point->GetLabAngle() << " deg at E_lab = " << point->GetLabEnergy()
                              << " MeV: " << point->GetLabAngleWhy() << "." << std::endl;
          return -1;
        }
        // The sub-points of a target effect were laid out (SetUpTargetEffects) from the copied base
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

  // CalcEDependentValues for component segments, in parallel as for the data segments
  // (EData::CalcEDependentValues).  Each point writes only itself, its sub-points and its local
  // mapped points that do not own their values; those are never computed here (the condition
  // below) and each has exactly one host, so no two threads write the same point.  Serial, this
  // loop was 9 of the 17.6 s of a natural_neon calculation on 24 threads.
  {
    std::vector<EPoint *> points;
    for (auto &componentSegment : componentSegments_)
      for (int i = 1; i <= componentSegment.NumPoints(); i++)
        if (!(componentSegment.GetPoint(i)->IsMapped()) || componentSegment.GetPoint(i)->OwnsEDependentValues())
          points.push_back(componentSegment.GetPoint(i));
    if (!points.empty()) {
      try {
        channelTable_ = BuildChannelTable(theCNuc, configure, points);
      } catch (std::exception &e) {
        configure.outStream << "Channel functions: " << e.what() << std::endl;
        return -1;
      }
    }
  }
  ChannelRowSource componentRows;
  if (channelTable_) componentRows.tables.push_back(channelTable_.get());
  const ChannelRowSource *componentTable = &componentRows;
  for (auto &componentSegment : componentSegments_) {
    bool localStop = false;
#pragma omp parallel for shared(localStop, configure)
    for (int i = 1; i <= componentSegment.NumPoints(); i++) {
      if (configure.stopFlag || localStop) continue;
      EPoint *point = componentSegment.GetPoint(i);
      if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
        try {
          point->CalcEDependentValues(theCNuc, configure, true, true, componentTable);
        } catch (GSLException e) {
#pragma omp critical
          {
            configure.outStream << "Component segment calculation error: " << e.what() << std::endl;
            localStop = true;
          }
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
      out.precision(17);  // as the data-segment amplitudes (CalculateECAmplitudes)
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
                  std::atomic<int> pointIndex(0);
                  time_t startTime = time(NULL);
                  bool localStop = false;
#pragma omp parallel for shared(configure, localStop, pointIndex, startTime)
                  for (int i = 1; i <= numPoints; i++) {
                    if (configure.stopFlag || localStop) continue;
                    EPoint *point = componentSegment.GetPoint(i);
                    if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
                      try {
                        point->CalculateECAmplitudes(theCNuc, configure, channelTable_.get());
                      } catch (GSLException e) {
#pragma omp critical(ec_progress)
                        {
                          configure.outStream << e.what() << std::endl;
                          localStop = true;
                        }
                      } catch (...) {
#pragma omp critical(ec_progress)
                        {
                          configure.outStream << "Warning: EC amplitude calculation failed for component segment point" << std::endl;
                        }
                      }
                    }
                    ++pointIndex;
                    ECProgressRedraw(configure, segmentKeyOut, pointIndex, numPoints, startTime);
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
    return false;
  }
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(entranceKey));
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(exitKey));
  if (entrancePair->GetNucleus() != exitPair->GetNucleus()) {
    configure.outStream << "ERROR: segment #" << baseSegment.GetSegmentKey() << ": component " << entranceKey
                        << " -> " << exitKey << " connects two compound nuclei ("
                        << theCNuc->NucleusLabel(entrancePair->GetNucleus()) << " and "
                        << theCNuc->NucleusLabel(exitPair->GetNucleus()) << ")." << std::endl;
    return false;
  }
  entrancePair->SetEntrance();
  return true;
}

/*!
 * A transmission or yield segment (ESegment::IsTransformed), once its components exist: refuse what the transform
 * cannot mean, make every component a plain cross section, and give a yield its hidden total-cross-section term
 * (entrance pair to itself, "total 1", by the optical theorem), the s_t of Y0 = (s / s_t)(1 - exp(-n s_t)).
 */
bool EData::FinishTransformedSegment(CNuc *theCNuc, ESegment *seg, int key, bool isTest, const Config &configure) {
  if (!seg || !seg->IsTransformed()) return true;
  const std::string what = std::string(isTest ? "test segment #" : "segment #") + std::to_string(key) + " (" +
                           (seg->IsTransmission() ? "transmission" : "yield") + "): ";
  auto refuse = [&](const std::string &why) {
    configure.outStream << "ERROR: " << what << why << std::endl;
    return false;
  };
  if (seg->IsTransmission() && seg->IsYield()) return refuse("a line is a transmission or a yield, not both.");
  if (seg->IsDifferential() || seg->IsPhase() || seg->IsAnalyzingPower() || seg->IsAngularDist() || seg->IsUPOS() ||
      seg->IsCoefficients() || seg->GetCrossSectionComponent() != 0)
    return refuse("the transform acts on an angle-integrated cross section (data type 0, or 3 for total capture).");
  if (seg->IsAdvanced() && seg->GetLegacyOperationType() != 0) return refuse("not on a RATIO.");
  PPair *entrance = theCNuc->GetPair(theCNuc->GetPairNumFromKey(seg->GetEntranceKey()));
  if (entrance->GetPType() != 0 || entrance->GetZ(1) * entrance->GetZ(2) != 0)
    return refuse("a sample transmission or self-shielded yield needs a neutral entrance pair (a neutron beam).");
  for (ESegment *comp : seg->GetComponentSegments())
    if (comp) comp->ClearTransform();
  if (seg->IsTransmission()) return true;

  if (seg->IsNatural())
    return refuse("not on a sample (abundance-weighted) segment: the total of a mixture is not implemented for Y0.");
  if (seg->GetExitKey() == seg->GetEntranceKey() || seg->IsTotal())
    return refuse("Y0 is a capture or reaction yield: the exit pair must differ from the entrance pair (no \"total 1\").");
  for (ESegment *comp : seg->GetComponentSegments()) {
    if (!comp) continue;
    if (comp->GetEntranceKey() != seg->GetEntranceKey() || comp->GetExitKey() == comp->GetEntranceKey())
      return refuse("every term of a yield must be a reaction of the segment's entrance pair (pair " +
                    std::to_string(seg->GetEntranceKey()) + ") to another exit pair.");
  }
  if (entrance->IsIdentical()) return refuse("not for a pair of identical particles.");
  if (configure.paramMask & Config::USE_RMC_FORMALISM)
    return refuse("the total cross section of Y0 needs the optical theorem, not available with --use-rmc.");
  ESegment *total = this->CreateComponentSegment(*seg, seg->GetEntranceKey(), seg->GetEntranceKey());
  if (!total) return refuse("could not create the total-cross-section term.");
  total->ClearTransform();
  total->SetComponentScaling(1.0);
  total->SetYieldTotalTerm();
  total->ApplyTotalToPoints();
  if (!seg->HasComponents()) seg->SetOperationType(SUM);
  seg->AddComponentSegment(total);
  return true;
}

/*!
 * After the target effects are read: a transmission or yield segment whose points carry a sub-point effect sums its
 * terms at every sub-point before the transform, so every term must sit on the base point's sub-point grid.  Terms of
 * one entrance pair do (their grids are laid out alike); the terms of a sample with several compound nuclei do not,
 * and an effect blended in at range edges or skipped by an automatic tolerance is not one kernel on transformed
 * values.  Refused here rather than computed wrongly.
 */
bool EData::CheckTransformedSegments(const Config &configure) {
  for (ESegment &seg : this->GetSegments()) {
    if (!seg.IsTransformed()) continue;
    const std::string what = "segment #" + std::to_string(seg.GetSegmentKey()) + " (" +
                             (seg.IsTransmission() ? "transmission" : "yield") + "): ";
    for (int p = 1; p <= seg.NumPoints(); p++) {
      EPoint *base = seg.GetPoint(p);
      if (!ESegment::HasSubPointIntegration(base)) {
        for (ESegment *comp : seg.GetComponentSegments())
          if (comp && p <= comp->NumPoints() && ESegment::HasSubPointIntegration(comp->GetPoint(p))) {
            configure.outStream << "ERROR: " << what << "a term is integrated over sub-points at E_lab "
                                << base->GetLabEnergy() << " MeV where the segment's own point is not." << std::endl;
            return false;
          }
        continue;
      }
      TargetEffect *te = this->GetTargetEffect(base->GetTargetEffectNum());
      // A thin-target integration blends with sigma(E_bar) inside its kernel (EPoint::TargetFallbackImpl), which is
      // still one linear kernel on the transformed values; other kernels blend with a plain point outside it.
      const bool thinTarget = te && te->IsTargetIntegration() && !te->IsThickTarget();
      if ((base->GetTargetBlendWeight() < 1.0 && !thinTarget) || (te && te->GetAutoTolerance() > 0.0)) {
        configure.outStream << "ERROR: " << what << "its target effect blends in at range edges or has an automatic "
                            << "tolerance; the transform needs the full kernel at every point (drop the blend width "
                            << "and the tolerance of that <targetInt> line)." << std::endl;
        return false;
      }
      for (ESegment *comp : seg.GetComponentSegments()) {
        EPoint *cp = (comp && p <= comp->NumPoints()) ? comp->GetPoint(p) : nullptr;
        bool same = cp && cp->NumSubPoints() == base->NumSubPoints();
        for (int j = 1; same && j <= base->NumSubPoints(); j++) {
          const double a = base->GetSubPoint(j)->GetCMEnergy(), b = cp->GetSubPoint(j)->GetCMEnergy();
          if (std::fabs(a - b) > 1.0e-9 * std::max(1.0, std::fabs(a))) same = false;
        }
        if (!same) {
          configure.outStream << "ERROR: " << what << "its terms do not share the base point's sub-point grid at E_lab "
                              << base->GetLabEnergy() << " MeV (a sample with several compound nuclei under a "
                              << "sub-point target effect is not supported yet: the transform needs every term at the "
                              << "same energies)." << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

/*!
 * A component of a SUM, RATIO, sample or total-capture segment: a copy of the base segment, its points included,
 * with its own pairs, kept in componentSegments_ (a deque: the base segments' pointers stay valid).  A fixed angle
 * (deg; a RATIO denominator at one angle) puts every point there, lab and c.m. slots alike.
 */
ESegment *EData::CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey,
                                        std::optional<double> fixedAngle) {
  this->componentSegments_.push_back(baseSegment);
  ESegment *componentSegment = &componentSegments_.back();
  componentSegment->SetEntranceKey(entranceKey);
  componentSegment->SetExitKey(exitKey);
  if (fixedAngle) {
    componentSegment->SetMinAngle(*fixedAngle);
    componentSegment->SetMaxAngle(*fixedAngle);
    for (EPoint &point : componentSegment->GetPoints()) {
      point.SetCMAngle(*fixedAngle);
      point.SetLabAngle(*fixedAngle);
    }
  }
  // A key of its own (negative: not a segment of the file), from the data set's counter, for the caches keyed by it
  componentSegment->SetSegmentKey(nextComponentKey_--);
  componentSegment->ClearComponents();  // no circular references
  componentSegment->SetNatural(false);  // the weighting is the base segment's job
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
 * (azure2::EvaluateModelPoint; ESegment::
 * UpdatePointEnergiesWithShift for the energy-dependent values, within a
 * segment).  The parallel loops (AZURECalc::PrecalculatePoints) hand out hosts
 * only, each of which calculates its mapped points.  The analytic gradient
 * (AZUREGrad.cpp) works on the same groups: one forward solve per group.
 *
 * Segments with components (SUM/RATIO, total capture) are not mapped; a
 * shifted segment without a particle exit (external capture) and a THM segment
 * only within themselves.
 * AZURE_NO_POINT_MAPPING=1 in the environment switches mapping off entirely;
 * AZURE_NO_CROSS_SEGMENT_MAPPING=1 keeps it within each segment only.
 */

void EData::MapData(CNuc *theCNuc) {
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    for (EPointIterator point = segment->GetPoints().begin(); point < segment->GetPoints().end(); point++)
      point->ClearMapping();
  }
  for (ESegment &component : componentSegments_)
    for (EPointIterator point = component.GetPoints().begin(); point < component.GetPoints().end(); point++)
      point->ClearMapping();
  // (host, mapped point) in the order each host calculates its mapped points
  std::vector<std::pair<azure2::PointRef, azure2::PointRef>> mapping;
  const char *off = std::getenv("AZURE_NO_POINT_MAPPING");
  if (off && *off && std::string(off) != "0") return plan_.Build(*this, mapping);
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
    // A THM (HOES) segment is mapped within itself only, and no other point is
    // mapped onto one of its points.  Its model is THMMatrixFunc's, from the
    // segment's own options (observable, vertex, kinematics, windows, line
    // shape, distortion, experiment, coherent background, 5-sigma fold), which
    // SharesCalculation does not compare; its gradient rows (ComputeTHMRows,
    // the THM-only differences of AZURECalc::Gradient) calculate THM segments
    // alone, so a host in another segment would not be up to date there.
    if (segments[s].IsTHM()) {
      ownOnly[s] = 1;
      continue;
    }
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
        if (SharesCalculation(pointOf(h), p)) {
          p.SetMap(h.seg + 1, h.point + 1);
          mapping.push_back({{azure2::PointRef::kSegments, h.seg, h.point},
                             {azure2::PointRef::kSegments, refs[k].seg, refs[k].point}});
          mapped = true;
          break;
        }
      }
      if (!mapped) hosts.push_back(refs[k]);
    }
    start = end;
  }

  MapComponentSegments(theCNuc, mapping);
  plan_.Build(*this, mapping);
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
void EData::MapComponentSegments(CNuc *theCNuc, std::vector<std::pair<azure2::PointRef, azure2::PointRef>> &mapping) {
  if (!theCNuc) return;
  const char *off = std::getenv("AZURE_NO_COMPONENT_MAPPING");
  if (off && *off && std::string(off) != "0") return;
  std::unordered_map<const ESegment *, int> componentPos;
  for (int c = 0; c < (int)componentSegments_.size(); c++) componentPos[&componentSegments_[c]] = c;
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
      if (!component || component->NumPoints() != n || !componentPos.count(component)) continue;
      for (int p = 0; p < n; p++) {
        EPoint *cp = component->GetPoint(p + 1);
        if (!cp || !SharesCalculation(points[p], *cp, false)) continue;
        cp->SetMap(s + 1, p + 1);
        cp->SetOwnsEDependentValues(true);
        mapping.push_back({{azure2::PointRef::kSegments, s, p}, {azure2::PointRef::kComponents, componentPos[component], p}});
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
      configure.Note("Note: shift group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: shift group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) + " +" +
                   members.str() + " share one constant energy shift (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s parameter, " +
                   (leader.IsVaryEnergyShift() ? "free" : "fixed") + ").");
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
        int svd;
        {
          GslQuiet quiet;  // scoped: a process-wide handler swap would outlive this check
          svd = gsl_linalg_SV_decomp(A, V, S, work);
        }
        if (svd == 0 && gsl_vector_get(S, k - 1) > 0.0)
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
      configure.Note(note.str());
      if (!(cond < 1.0e8))
        configure.outStream << "WARNING: shift function on segment #" << key << ": its free terms cannot be told apart "
                            << "over this group's energies (condition number " << cond << "); fix one of them or "
                            << "constrain it with a penalty." << std::endl;
    } else {
      configure.Note(note.str());
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
      configure.Note("Note: normalization group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: normalization group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) +
                   " +" + members.str() + " share one normalization (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s, " + (leader.IsVaryNorm() ? "free" : "fixed") + ").");
  }
}

int EData::NormOwner(int s) const {
  if (s < 1 || s > (int)segments_.size()) return s;
  const int lead = segments_[s - 1].GetNormLeaderPos();
  return (lead >= 0) ? lead + 1 : s;
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
  // Thickness-group followers likewise: their (fixed) segment_<key>_thickness is written as the leader's value.
  if (thicknessParamOffset_ >= 0) {
    for (int k = 1; k <= this->NumThicknessParams(); k++) {
      const ESegment &seg = segments_[thicknessSegments_[k - 1] - 1];
      const int idx = thicknessParamOffset_ + k - 1;
      if (!seg.IsThicknessGroupFollower() || idx >= (int)p.Params().size()) continue;
      const int lk = this->ThicknessParamOfSegment(thicknessSegments_[k - 1]);
      const int lidx = (lk > 0) ? thicknessParamOffset_ + lk - 1 : -1;
      p.SetValue(idx, (lidx >= 0 && lidx < (int)p.Params().size()) ? p.Value(lidx)
                                                                     : segments_[seg.GetThicknessLeaderPos()].GetSampleThickness());
    }
  }
}

/*!
 * Fills the Minuit parameter array from initial values in the EData object.
 */

namespace {

// The name of a target effect's parameter: from its target in format 3 (<target>.areal_density), else the old
// <kind>_<effect number>.
std::string EffectParameterName(const TargetEffect &te, const char *oldKind, const char *field, int effect) {
  return te.GetTargetName().empty() ? std::string(oldKind) + "_" + std::to_string(effect) : te.GetTargetName() + field;
}

}  // namespace

void EData::FillMnParams(AZUREParams &params) {
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // 1-based position of a segment in the list (the records' and ParamIndexMap's numbering).
  auto SegmentPosition = [this](ESegmentIterator it) { return (int)(it - GetSegments().begin()) + 1; };
  SetNormParamOffset(p.Params().size());
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      const std::string &varname = segment->GetParameterNames().norm;
      params.Add(varname, segment->GetNorm(), segment->GetNorm() * 0.05,
                 ParameterRecord::OfSegment(ParamKind::Norm, SegmentPosition(segment), segment->GetSegmentKey()));
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
    const std::string &varname = segment->GetParameterNames().energyShift;
    double stepSize = (segment->GetEnergyShiftError() > 0.0) ? segment->GetEnergyShiftError() * 0.01 : 0.0005;
    params.Add(varname, segment->GetEnergyShift(), stepSize,
               ParameterRecord::OfSegment(ParamKind::EnergyShift, SegmentPosition(segment), segment->GetSegmentKey()));

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
    const std::string &varname = segment->GetParameterNames().energyShiftSqrt;
    double stepSize = (segment->GetEnergyShiftSqrtError() > 0.0) ? segment->GetEnergyShiftSqrtError() * 0.01 : 0.0005;
    params.Add(varname, segment->GetEnergyShiftSqrt(), stepSize,
               ParameterRecord::OfSegment(ParamKind::EnergyShiftSqrt, SegmentPosition(segment), segment->GetSegmentKey()));
    if (!segment->IsVaryEnergyShiftSqrt()) p.Fix(varname);
  }

  // Shift-function parameters: only segments that declare a "shiftfn" block have
  // any, so a model without one gets no new parameter at all.
  SetEnergyShiftFnParamOffset(p.Params().size());
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (!segment->OwnsShiftFn() || !segment->HasShiftFn()) continue;
    for (int j = 0; j < segment->NumShiftFnParams(); j++) {
      const std::string &varname = segment->GetParameterNames().shiftFunction.at(j);
      const double err = segment->GetShiftFnErrors()[j];
      params.Add(varname, segment->GetShiftFnValues()[j], (err > 0.0) ? err * 0.01 : 0.0005,
                 ParameterRecord::ShiftFnTerm(SegmentPosition(segment), segment->GetSegmentKey(), j));
      if (!segment->IsVaryShiftFn(j)) p.Fix(varname);
    }
  }

  // Electron-screening potentials (eV): one per target effect with a "screening" block, none otherwise.
  screeningParamOffset_ = p.Params().size();
  for (int k = 1; k <= this->NumScreeningParams(); k++) {
    const TargetEffect *te = this->GetScreeningTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "screening", ".screening", this->GetScreeningEffect(k));
    const double ue = te->GetScreeningUe();
    params.Add(varname, ue, std::max(10.0, 0.05 * std::fabs(ue)),
               ParameterRecord::OfEffect(ParamKind::Screening, this->GetScreeningEffect(k)));
    p.SetLowerLimit(varname, 0.0);
    if (!te->IsScreeningVaried()) p.Fix(varname);
  }

  // Fitted target densities (atoms/cm^2): one per densityfit line or density group, none otherwise.  Bounded by
  // the density the sub-point grids span (TargetEffect::GetLayoutDensity) and a hundredth of the line's.
  densityParamOffset_ = p.Params().size();
  for (int k = 1; k <= this->NumDensityParams(); k++) {
    const TargetEffect *te = this->GetDensityTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "density", ".areal_density", this->GetDensityEffect(k));
    const double n = te->GetDensity();
    params.Add(varname, n, (te->GetDensitySigma() > 0.0) ? 0.1 * te->GetDensitySigma() : 0.01 * n,
               ParameterRecord::OfEffect(ParamKind::Density, this->GetDensityEffect(k)));
    p.SetLimits(varname, 0.01 * te->GetDensityNominal(), te->GetLayoutDensity());
    if (!te->IsDensityVaried()) p.Fix(varname);
  }

  // Fitted Doppler temperatures (K): one per line with the fit form 'doppler T vary sigma [Tmax]', none otherwise.
  // Bounded by the temperature the grids span (TargetEffect::GetDopplerLayoutTemperature) and T0/4.
  dopplerParamOffset_ = p.Params().size();
  for (int k = 1; k <= this->NumDopplerParams(); k++) {
    const TargetEffect *te = this->GetDopplerTargetEffect(k);
    const std::string varname = EffectParameterName(*te, "doppler", ".doppler", this->GetDopplerEffect(k));
    const double t = te->GetDopplerTemperature();
    params.Add(varname, t, (te->GetDopplerSigma() > 0.0) ? 0.1 * te->GetDopplerSigma() : 0.01 * t,
               ParameterRecord::OfEffect(ParamKind::Doppler, this->GetDopplerEffect(k)));
    p.SetLimits(varname, te->GetDopplerMinimum(), te->GetDopplerLayoutTemperature());
    if (!te->IsDopplerVaried()) p.Fix(varname);
  }

  // Fitted sample thicknesses (atoms/barn): one per data segment whose transmission / yield block has the fit form
  // 'n vary sigma', none otherwise.  A thickness-group follower's is fixed (its value is its leader's).
  thicknessParamOffset_ = p.Params().size();
  for (int k = 1; k <= this->NumThicknessParams(); k++) {
    ESegment *seg = this->GetSegment(this->GetThicknessSegment(k));
    const std::string &varname = seg->GetParameterNames().thickness;
    const double n = seg->GetSampleThickness();
    params.Add(varname, n, (seg->GetThicknessSigma() > 0.0) ? 0.1 * seg->GetThicknessSigma() : 0.01 * n,
               ParameterRecord::OfSegment(ParamKind::Thickness, this->GetThicknessSegment(k), seg->GetSegmentKey()));
    p.SetLowerLimit(varname, 0.0);
    if (!seg->IsThicknessVaried()) p.Fix(varname);
  }

  // THM coherent backgrounds (cbackground=): Re/Im of each amplitude, last.
  thmCoherentParamOffset_ = thmCoherentParams_.empty() ? -1 : (int)p.Params().size();
  for (const ThmCoherentParam &c : thmCoherentParams_) {
    params.Add(c.name, c.value, c.value != 0.0 ? 0.1 * std::fabs(c.value) : 0.1,
               ParameterRecord::ThmCoherent(c.experiment, c.jGroup, c.entrance, c.exit, c.term));
    if (c.fixed) p.Fix(c.name);
  }
}

int EData::ResolveThicknessGroups(const Config &configure) {
  std::vector<std::pair<std::string, int>> leaders;  // label, leader position (0-based)
  for (int pos = 0; pos < (int)segments_.size(); pos++) {
    ESegment &seg = segments_[pos];
    const std::string label = seg.GetThicknessGroup();
    if (label.empty()) continue;
    if (!seg.IsTransformed()) {
      configure.outStream << "ERROR: segment #" << seg.GetSegmentKey() << ": thicknessgroup '" << label
                          << "' needs the segment's own transmission or yield block (the group shares the sample "
                          << "thickness of transmission / yield segments)." << std::endl;
      return -1;
    }
    auto it = std::find_if(leaders.begin(), leaders.end(),
                           [&label](const std::pair<std::string, int> &e) { return e.first == label; });
    if (it == leaders.end()) {
      leaders.push_back({label, pos});
      continue;
    }
    seg.SetThicknessGroupFollower(it->second, segments_[it->second].GetSampleThickness());
  }
  for (const auto &g : leaders) {
    const ESegment &leader = segments_[g.second];
    std::ostringstream members;
    int n = 0;
    for (const ESegment &seg : segments_)
      if (seg.GetThicknessLeaderPos() == g.second) {
        members << " " << seg.GetSegmentKey();
        n++;
      }
    if (n == 0) {
      configure.Note("Note: thickness group '" + g.first + "' has only segment #" +
                     std::to_string(leader.GetSegmentKey()) + "; nothing is shared.");
      continue;
    }
    configure.Note("Note: thickness group '" + g.first + "': segments #" + std::to_string(leader.GetSegmentKey()) +
                   " +" + members.str() + " share one sample thickness (segment #" +
                   std::to_string(leader.GetSegmentKey()) + "'s, " +
                   (leader.HasThicknessFit() ? "parameter " + leader.GetParameterNames().thickness
                                             : std::string("fixed")) +
                   ").");
  }
  thicknessSegments_.clear();
  for (int pos = 0; pos < (int)segments_.size(); pos++)
    if (segments_[pos].HasThicknessFit()) thicknessSegments_.push_back(pos + 1);
  return 0;
}

int EData::ThicknessParamOfSegment(int pos) const {
  if (pos < 1 || pos > (int)segments_.size()) return 0;
  const int lead = segments_[pos - 1].GetThicknessLeaderPos();
  const int owner = (lead >= 0) ? lead + 1 : pos;
  for (int k = 1; k <= this->NumThicknessParams(); k++)
    if (thicknessSegments_[k - 1] == owner) return k;
  return 0;
}

void EData::FillThicknessFromParams(const vector_r &p) {
  if (thicknessSegments_.empty() || thicknessParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumThicknessParams(); k++) {
    const int idx = thicknessParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the thickness block
    ESegment &seg = segments_[thicknessSegments_[k - 1] - 1];
    if (seg.IsThicknessGroupFollower()) continue;
    seg.SetSampleThickness(std::max(0.0, p[idx]));
  }
  for (ESegment &seg : segments_)
    if (seg.IsThicknessGroupFollower())
      seg.SetSampleThickness(segments_[seg.GetThicknessLeaderPos()].GetSampleThickness());
}

double EData::ThicknessPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumThicknessParams(); k++) chi2 += segments_[thicknessSegments_[k - 1] - 1].ThicknessPenalty();
  return chi2;
}

int EData::DopplerParamOfEffect(int n) const {
  for (int k = 1; k <= this->NumDopplerParams(); k++)
    for (int c : dopplerCopies_[k - 1])
      if (c == n) return k;
  return 0;
}

void EData::FillDopplerFromParams(const vector_r &p, const Config *configure) {
  if (dopplerEffects_.empty() || dopplerParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumDopplerParams(); k++) {
    const int idx = dopplerParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the Doppler block
    TargetEffect *te = this->GetDopplerTargetEffect(k);
    double v = p[idx];
    // The grids and the integration windows span the layout temperature; below T0/4 the Gauss-Legendre pieces
    // (laid out for T0) would be too wide for the kernel.
    const double vMax = te->GetDopplerLayoutTemperature(), vMin = te->GetDopplerMinimum();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        (configure ? configure->outStream : std::cout)
            << "WARNING: doppler_" << this->GetDopplerEffect(k) << " = " << v << " K is outside [" << vMin << ", "
            << vMax << "] K, the range its sub-point grids and integration pieces span; the calculation uses the"
            << " nearest end.  Give the line a wider range with 'doppler T vary sigma Tmax'. (Printed once.)"
            << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n : dopplerCopies_[k - 1]) this->GetTargetEffect(n)->SetDopplerTemperature(v);
  }
}

double EData::DopplerPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumDopplerParams(); k++)
    chi2 += this->GetDopplerTargetEffect(k)->DopplerPenalty();
  return chi2;
}

/*!
 * Groups the densityfit / densitygroup lines into parameters.  SetUpTargetEffects gives every segment a line
 * lists its own copy of the effect (and every component of such a segment one more), so the effects are visited
 * by line (TargetEffect::GetLineNumber): the first copy of a line decides, the other copies follow it.  A
 * line that lists no existing segment has no copy and so no parameter (dev makes one that moves nothing).
 */
int EData::SetUpDensityParams(const Config &configure) {
  densityEffects_.clear();
  densityLeadCopies_.clear();
  densityParamOfEffect_.assign(this->NumTargetEffects(), 0);
  std::map<int, int> firstCopyOfLine;       // line number -> its first copy (target effect number)
  std::map<std::string, int> groupLeader;   // label -> first copy of the group's first line
  for (int n = 1; n <= this->NumTargetEffects(); n++) {
    TargetEffect *te = this->GetTargetEffect(n);
    const int line = te->GetLineNumber();
    auto seen = firstCopyOfLine.find(line);
    if (seen != firstCopyOfLine.end()) {
      // Another copy of a line already grouped: the same parameter, and the group's values if it joined one.
      const TargetEffect *lineFirst = this->GetTargetEffect(seen->second);
      if (lineFirst->HasDensityFit()) te->JoinDensityGroup(*lineFirst);
      densityParamOfEffect_[n - 1] = densityParamOfEffect_[seen->second - 1];
      continue;
    }
    firstCopyOfLine[line] = n;
    const std::string &group = te->GetDensityGroup();
    if (!group.empty()) {
      auto it = groupLeader.find(group);
      if (it != groupLeader.end()) {
        const TargetEffect *lead = this->GetTargetEffect(it->second);
        const int leadLine = lead->GetLineNumber();
        if (std::fabs(te->GetDensity() - lead->GetDensity()) > 1.0e-9 * lead->GetDensity()) {
          std::ostringstream note;
          note << "Note: <targetInt> effect " << line << " (density group " << group << ") has density "
               << te->GetDensity() << "; the group's density is its first line's, " << lead->GetDensity()
               << " (effect " << leadLine << ").";
          configure.Note(note.str());
        }
        if (te->HasDensityFit() && (te->IsDensityVaried() != lead->IsDensityVaried() ||
                                    std::fabs(te->GetDensitySigma() - lead->GetDensitySigma()) > 1.0e-9 * lead->GetDensity()))
          configure.Note("Note: <targetInt> effect " + std::to_string(line) + " (density group " + group +
                         "): its densityfit block is ignored, the group's is its first line's (effect " +
                         std::to_string(leadLine) + ").");
        te->JoinDensityGroup(*lead);
        densityParamOfEffect_[n - 1] = densityParamOfEffect_[it->second - 1];
        continue;
      }
      if (!te->HasDensityFit()) {
        configure.outStream << "ERROR: density group " << group << ": its first <targetInt> line (effect " << line
                            << ") has no 'densityfit vary sigma' block; the group's parameter and prior are read from it."
                            << std::endl;
        return -1;
      }
      groupLeader[group] = n;
    }
    if (!te->HasDensityFit()) continue;
    densityEffects_.push_back(line);
    densityLeadCopies_.push_back(n);
    densityParamOfEffect_[n - 1] = (int)densityEffects_.size();
  }
  return 0;
}

void EData::FillDensitiesFromParams(const vector_r &p, const Config &configure) {
  if (densityEffects_.empty() || densityParamOffset_ < 0) return;
  std::vector<char> changed(this->NumTargetEffects(), 0);
  bool any = false;
  for (int k = 1; k <= this->NumDensityParams(); k++) {
    const int idx = densityParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the density block
    const TargetEffect *lead = this->GetDensityTargetEffect(k);
    double v = p[idx];
    // The sub-point grids span the layout density: beyond it the depth integral would run off the grid.
    const double vMax = lead->GetLayoutDensity(), vMin = 0.01 * lead->GetDensityNominal();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        configure.outStream << "WARNING: density_" << this->GetDensityEffect(k) << " = " << v << " is outside ["
                            << vMin << ", " << vMax << "], the range its sub-point grids span; the calculation uses the"
                            << " nearest end.  Give the line a wider range with 'densityfit vary sigma max'."
                            << " (Printed once.)" << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n = 1; n <= this->NumTargetEffects(); n++) {
      if (densityParamOfEffect_[n - 1] != k) continue;
      TargetEffect *te = this->GetTargetEffect(n);
      if (te->GetDensity() == v) continue;
      te->SetDensity(v);
      changed[n - 1] = 1;
      any = true;
    }
  }
  if (!any) return;
  auto update = [&](ESegment &seg) {
    if (!seg.IsTargetEffect()) return;
    const int n = seg.GetTargetEffectNum();
    if (n < 1 || n > (int)changed.size() || !changed[n - 1]) return;
    TargetEffect *te = this->GetTargetEffect(n);
    for (EPointIterator point = seg.GetPoints().begin(); point < seg.GetPoints().end(); point++) {
      if (point->GetTargetEffectNum() != n || !(point->GetTargetThickness() > 0.0)) continue;
      point->SetTargetThickness(point->GetCMFactor() * te->TargetThickness(point->GetTargetLabEnergy(), configure));
      // E_bar of the sigma(E_bar) fallback moves with the density too (EPoint::SetTargetFallbackOnly)
      if (!te->IsThickTarget())
        point->SetTargetFallbackDepth(point->GetCMFactor() *
                                      te->MeanEnergyLoss(point->GetTargetLabEnergy(), configure, te->GetDensity()));
    }
  };
  for (ESegment &seg : segments_) update(seg);
  for (ESegment &seg : componentSegments_) update(seg);
}

double EData::DensityPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumDensityParams(); k++)
    chi2 += this->GetDensityTargetEffect(k)->DensityPenalty();
  return chi2;
}

void EData::ApplyScreening() {
  if (screeningEffects_.empty()) return;
  std::vector<int> paramOf(this->NumTargetEffects() + 1, 0);
  for (int k = 1; k <= this->NumScreeningParams(); k++)
    for (int n : screeningCopies_[k - 1]) paramOf[n] = k;
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
    for (int n : screeningCopies_[k - 1]) this->GetTargetEffect(n)->SetScreeningUe(p[idx]);
  }
  this->ApplyScreening();
}

double EData::ScreeningPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumScreeningParams(); k++)
    chi2 += this->GetScreeningTargetEffect(k)->ScreeningPenalty();
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
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
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

void EData::FillEnergyShiftsFromParams(const vector_r &p, azure2::Results &results, CNuc *theCNuc,
                                       const Config *configure) {
  this->FillEnergyShiftsFromParamsImpl(p, theCNuc, configure);
  this->BeginEvaluation(results);
}

template <class F>
void EData::ForEachResultUnit(F f) {
  const int nNatural = numNuclei_;
  auto segment = [&](ESegment &seg) {
    const int nat = seg.IsNatural() ? nNatural : 0;
    for (EPoint &point : seg.GetPoints()) {
      f(point, point.ResultAngularCapacity(), nat);
      for (EPoint &sub : point.GetSubPoints()) f(sub, sub.ResultAngularCapacity(), 0);
    }
  };
  for (ESegment &seg : segments_) segment(seg);
  for (ESegment &seg : componentSegments_) segment(seg);
}

void EData::LayoutResults(azure2::Results *results) {
  size_t total = 0;
  ForEachResultUnit([&](EPoint &, int nAng, int nNat) { total += azure2::Results::RecordSize(nAng, nNat); });
  std::vector<double> next(results ? total : 0, 0.0);
  size_t offset = 0;
  ForEachResultUnit([&](EPoint &point, int nAng, int nNat) {
    if (results) {
      const size_t from = point.ResultSlot();
      const int fromAng = point.ResultAngularSlots(), fromNat = point.ResultNaturalSlots();
      if (from != EPoint::kNoRecord && from + azure2::Results::RecordSize(fromAng, fromNat) <= results->slot.size())
        azure2::CopyRecord(results->slot.data() + from, fromAng, fromNat, next.data() + offset, nAng, nNat);
    }
    point.SetResultRecord(offset, nAng, nNat);
    offset += azure2::Results::RecordSize(nAng, nNat);
  });
  resultsSize_ = total;
  if (results) results->slot.swap(next);
  modelSize_ = 0;
  for (ESegment &seg : segments_)
    for (EPoint &point : seg.GetPoints()) {
      point.SetModelRecord(modelSize_, point.ResultAngularCapacity());
      modelSize_ += azure2::Results::ModelRecordSize(point.ResultAngularCapacity());
    }
}

void EData::BeginEvaluation(azure2::Results &results) {
  size_t offset = 0;
  bool laidOut = true;
  ForEachResultUnit([&](EPoint &point, int nAng, int nNat) {
    if (!laidOut) return;
    laidOut = point.ResultSlot() == offset && point.ResultAngularSlots() == nAng && point.ResultNaturalSlots() == nNat;
    offset += azure2::Results::RecordSize(nAng, nNat);
  });
  if (!laidOut || offset != resultsSize_) this->LayoutResults(&results);
  this->SizeResults(results);
  if (azure2::Results::PoisonRequested()) results.Poison();
}

void EData::FillEnergyShiftsFromParamsImpl(const vector_r &p, CNuc *theCNuc, const Config *configure) {
  // The THM coherent backgrounds, the block after the shifts (none without
  // cbackground=): every caller fills the shifts right after the norms, with
  // the compound the model will run on.
  FillThmCoherentFromParams(p, theCNuc);
  // Screening first: sub-points rebuilt below (ReplaceSubPointGrids) inherit their parent's U_e.
  this->FillScreeningFromParams(p);
  // Fitted target densities: the thickness of each point of their lines (the grids stay put).
  if (configure) this->FillDensitiesFromParams(p, *configure);
  // Fitted Doppler temperatures: only the kernel values move (the grids and pieces stay put).
  this->FillDopplerFromParams(p, configure);
  // Fitted sample thicknesses: the n of the transmission / yield transforms.
  this->FillThicknessFromParams(p);
  int i = GetEnergyShiftParamOffset();
  int j = GetEnergyShiftSqrtParamOffset();
  bool anyEnergyShifted = false;

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
  for (ESegmentIterator segment = this->GetSegments().begin(); segment < this->GetSegments().end(); segment++) {
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
      const int start = this->GetSegment(segment->GetShiftFnSourcePos() + 1)->GetShiftFnParamStart();
      for (int q = 0; q < segment->NumShiftFnParams(); q++) {
        const int idx = GetEnergyShiftFnParamOffset() + start + q;
        fnValues.push_back((start >= 0 && idx < (int)p.size()) ? p[idx] : segment->GetShiftFnValues()[q]);
      }
    }
    // A segment whose points carry a shift of an earlier evaluation (last shift not NaN) stays active, so a fixed
    // shift set back to 0 moves them back (upgrade_history: such a copy gives a fresh copy's bits).
    const bool shiftActive = segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 || segment->HasShiftFn() ||
                             segment->IsVaryEnergyShiftSqrt() || shiftSqrt != 0.0 || !std::isnan(segment->GetLastEnergyShift());
    if (shiftActive) {
      shiftActiveSegments.push_back(&*segment);
      if (segment->HasComponents())
        for (ESegment *componentSegment : segment->GetComponentSegments())
          if (componentSegment) shiftActiveSegments.push_back(componentSegment);
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

    i++;
    j++;
  }

  // The sub-point grids, for the levels and shifts just set.  Sub-points of
  // shifted segments are moved and energised by the step below.
  if (theCNuc && configure) {
    std::vector<ESegment *> reshift;
    if (this->ReplaceSubPointGrids(theCNuc, *configure, shiftActiveSegments, reshift) > 0) anyEnergyShifted = true;
    for (ESegment *s : reshift) queueUpdate(s);
  }

  if (!toUpdate.empty()) {
    if (theCNuc && configure)
      for (ESegment *s : toUpdate) s->BeginShiftRefill(theCNuc, *configure);
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
        const double nMapped = pt ? pt->PlanMapped().size() : 0;
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
    if (theCNuc && configure)
      for (ESegment *s : toUpdate) s->CommitShiftTable(theCNuc, *configure);
  }

  // The energy-based point mapping (MapData) follows the shifted energies.
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
  dataCopy->iterations_.store(this->iterations_.load(std::memory_order_relaxed), std::memory_order_relaxed);
  dataCopy->normParamOffset_ = this->normParamOffset_;
  dataCopy->energyShiftParamOffset_ = this->energyShiftParamOffset_;
  dataCopy->energyShiftSqrtParamOffset_ = this->energyShiftSqrtParamOffset_;
  dataCopy->energyShiftFnParamOffset_ = this->energyShiftFnParamOffset_;
  dataCopy->screeningParamOffset_ = this->screeningParamOffset_;
  dataCopy->screeningEffects_ = this->screeningEffects_;
  dataCopy->screeningCopies_ = this->screeningCopies_;
  dataCopy->densityParamOffset_ = this->densityParamOffset_;
  dataCopy->densityEffects_ = this->densityEffects_;
  dataCopy->densityLeadCopies_ = this->densityLeadCopies_;
  dataCopy->densityParamOfEffect_ = this->densityParamOfEffect_;
  dataCopy->dopplerParamOffset_ = this->dopplerParamOffset_;
  dataCopy->dopplerEffects_ = this->dopplerEffects_;
  dataCopy->dopplerCopies_ = this->dopplerCopies_;
  dataCopy->thicknessParamOffset_ = this->thicknessParamOffset_;
  dataCopy->thicknessSegments_ = this->thicknessSegments_;
  dataCopy->isFit_ = this->isFit_;
  dataCopy->isErrorAnalysis_ = this->isErrorAnalysis_;
  dataCopy->targetEffects_ = this->targetEffects_;
  dataCopy->segments_ = this->segments_;
  dataCopy->componentSegments_ = this->componentSegments_;
  dataCopy->plan_ = this->plan_;
  dataCopy->thmGroups_ = this->thmGroups_;
  dataCopy->thmCoherentParams_ = this->thmCoherentParams_;
  dataCopy->thmCoherentParamOffset_ = this->thmCoherentParamOffset_;
  dataCopy->covGroups_ = this->covGroups_;  // shared: one factor (and one frozen model) for all copies
  dataCopy->nodeCache_ = this->nodeCache_;  // shared: a re-placed node is computed and held once
  dataCopy->priorCov_ = this->priorCov_;    // shared, read-only
  dataCopy->numNuclei_ = this->numNuclei_;

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

  // The points of the copy belong to it (a fitted density or Doppler temperature reaches them through it).  The plan
  // and the points' task indices are positions, valid in the copy as they are: no point holds a pointer to another.
  for (ESegment &segment : dataCopy->segments_)
    for (EPoint &point : segment.GetPoints()) point.SetParentData(dataCopy), point.ClearECAmplitudes();
  for (ESegment &segment : dataCopy->componentSegments_)
    for (EPoint &point : segment.GetPoints()) point.SetParentData(dataCopy), point.ClearECAmplitudes();
  // The copy's layout: its points carry this data set's offsets, which are the copy's own when the copy has the
  // same points; the values are in the buffers of the callers that evaluate the two.
  dataCopy->LayoutResults(nullptr);
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

/*!
 * One parallel pass over the points of an evaluation, shared by the CLI
 * (AZURECalc::PrecalculatePoints) and the API's calculations (AZUREAPI).
 */

int EData::CalculatePointsInParallel(const std::vector<ESegment *> &segments, CNuc *theCNuc,
                                     const Config &configure, azure2::Results &results) {
#ifdef _OPENMP
  if (omp_in_parallel()) return -1;
#else
  return -1;
#endif
  if (ScopedPointMappingOff::Active()) return -1;
  if (std::getenv("AZURE_SERIAL_CHI2")) return -1;
  std::vector<EPoint *> work;
  std::set<ESegment *> seen;
  auto add = [&](ESegment *seg) {
    if (!seg || !seen.insert(seg).second) return;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (pt && !pt->IsMapped()) work.push_back(pt);
    }
  };
  for (ESegment *seg : segments) {
    if (!seg) continue;
    add(seg);
    for (ESegment *comp : seg->GetComponentSegments()) add(comp);
  }
  // A host also calculates every point mapped onto it (EData::MapData), so
  // its work grows with that number: one Fowler-type angular distribution at
  // a single energy is one item.  Largest items first, one at a time, so such
  // an item does not end up last on one thread.  Each item is calculated the
  // same way whatever the order.
  std::vector<double> cost(work.size());
  for (size_t k = 0; k < work.size(); k++)
    cost[k] = (1.0 + work[k]->NumSubPoints()) * (1.0 + work[k]->PlanMapped().size());
  std::vector<size_t> order(work.size());
  for (size_t k = 0; k < order.size(); k++) order[k] = k;
  std::stable_sort(order.begin(), order.end(), [&cost](size_t a, size_t b) { return cost[a] > cost[b]; });
  std::vector<EPoint *> sorted(work.size());
  for (size_t k = 0; k < order.size(); k++) sorted[k] = work[order[k]];
  work.swap(sorted);
  int failed = 0;
  AZUREProfile::Timer profileTimer(AZUREProfile::kPoints);
  // One item per task of the plan: its host calculates itself and the points the plan maps onto it.
  auto runTask = [&](EPoint *host) {
    try {
      host->Calculate(theCNuc, configure, results);
    } catch (const std::exception &e) {
      // NaN on the point and the points mapped onto it, reported once (EPoint::MarkFailed); 0 on the point alone,
      // silently, before 2026-11-01 -- and the CLI's objective takes these values (PrecalculatePoints).
      host->MarkFailed(results, "point", configure, e.what());
      return 1;
    } catch (...) {
      host->MarkFailed(results, "point", configure, "unknown exception.");
      return 1;
    }
    return 0;
  };
  if (azure2::CalculationPlan::RunTimeCheckRequested() && !work.empty()) {
    azure2::CalculationPlan::RunChecked(*work[0]->GetParentData(), work, results, [&](EPoint *h) { failed += runTask(h); });
    return failed;
  }
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : failed)
  for (int k = 0; k < (int)work.size(); k++) failed += runTask(work[k]);
  return failed;
}
