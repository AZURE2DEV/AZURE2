#include "InputFiles.h"
#include "AZUREOutput.h"
#include "Profile.h"
#include "CNuc.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "WorkingCopy.h"
#include "ModelPass.h"
#include "ECAmplitudeCache.h"
#include "ECTable.h"
#include "AdaptiveIntegrationGrid.h"
#include "Minuit2/MnUserParameters.h"
#include "AZUREParams.h"
#include "GSLException.h"
#include "ThmVertexBoundary.h"
#include <gsl/gsl_linalg.h>
#include <map>
#include <tuple>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cmath>
#include <limits>
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
static bool CoefficientsAllowed(CNuc *theCNuc, const Segment &segment, int number, const Config &configure) {
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

static bool AnalyzingPowerAllowed(CNuc *theCNuc, const Segment &segment, int number, bool isTest,
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
static void UPOSLabNote(const Segment &segment, int key, bool isTest, const Config &configure) {
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
static bool TotalCaptureCodeOk(const Segment &segment, int key, bool isTest, const Config &configure) {
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
static bool UPOSAllowed(CNuc *theCNuc, const Segment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !theCNuc->IsPairKey(segment.GetExitKey())) return true;
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetExitKey()));
  if (exitPair->GetLightExE() == 0.) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << key << " (UPOS): exit pair "
                      << segment.GetExitKey() << " has an excited light particle (column 32, " << exitPair->GetLightExE()
                      << " MeV); the secondary-gamma correlation is implemented for the decay of the heavy partner only."
                      << std::endl;
  return false;
}

static bool TotalAllowed(CNuc *theCNuc, Segment &segment, int key, bool isTest, const Config &configure) {
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
// converted with that pair's masses (Segment::ConvertPointEnergyToCM), so a segment without one (exit -1) is refused.
static bool BetaDelayedAllowed(CNuc *theCNuc, Segment &segment, int key, bool isTest, const Config &configure) {
  if (!theCNuc->IsPairKey(segment.GetEntranceKey()) || segment.GetExitKey() != -1) return true;
  if (theCNuc->GetPair(theCNuc->GetPairNumFromKey(segment.GetEntranceKey()))->GetPType() != 20) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment" : "segment") << " #" << key << ": a beta-delayed "
                      << "spectrum (entrance pair " << segment.GetEntranceKey() << " of type 20) needs its decay pair as"
                      << " the exit pair, not -1 (total capture)." << std::endl;
  return false;
}

static bool NaturalAllowed(CNuc *theCNuc, Segment &segment, int key, bool isTest, const Config &configure) {
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
 * The WorkingCopy object has a private attribute containing the number of iterations needed to find the best fit parameters.
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
 * Checks the THM weight tables of the records (a data set's or extrapolation's thm_weight, ThmOptions) against the
 * points just read: every point of the data set must lie inside its table (folding sub-points may reach beyond it;
 * they get the end value).  The reading has refused a weight on a record that is not THM; the weight of a record not
 * in use is ignored with a warning.  Returns -1 after an ERROR line.
 */
int CheckThmWeights(const Config &configure,
                    const std::map<int, std::shared_ptr<const ThmWeightTable>> &tables,
                    const char *keyName, WorkingCopy *data) {
  for (const auto &entry : tables) {
    const int key = entry.first;
    const ThmWeightTable &table = *entry.second;
    const Segment *segment = nullptr;
    for (int s = 1; s <= data->Data().NumSegments(); s++)
      if (data->Data().GetSegment(s)->GetSegmentKey() == key) segment = data->Data().GetSegment(s);
    if (!segment) {
      configure.outStream << "WARNING: <thm> " << keyName << "[" << key << "]: the data set or extrapolation with key "
                          << key << " is not in use; the weight is ignored." << std::endl;
      continue;
    }
    for (int p = 1; p <= segment->NumPoints(); p++) {
      double energy = segment->GetPoint(*data, p)->GetCMEnergy();
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

WorkingCopy::WorkingCopy() : draft_(std::make_shared<DataSet>()) {
  set_ = draft_;
  isFit_ = true;
  isErrorAnalysis_ = false;
  ecReadPos_ = std::streampos(0);
  ecUsePrevious_ = false;
}

/*!
 * This function fills the data object with segments from the segment data files.
 * After a segment is created, the Segment::Fill method is called for that segment.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int WorkingCopy::Fill(const Config &configure, CNuc *theCNuc, const azure2::experiment::DataSpecs &specs) {
  const int numTotalSegments_ = (int)specs.datasets.items.size();
  for (const auto &item : specs.datasets.items) {
    {
      const azure2::experiment::DatasetSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        Segment NewSegment(segment);

        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          entrancePair->SetEntrance();
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!CoefficientsAllowed(theCNuc, NewSegment, numTotalSegments, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!BetaDelayedAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!UPOSAllowed(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!ThmVertexAllowsPair(theCNuc, NewSegment, numTotalSegments, false, configure)) return -1;
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
            if (Draft().GetSegment(Draft().NumSegments())->Fill(theCNuc, this, configure) == -1) {
              if (covarianceError_ || kinematicsError_) return -1;  // never skipped silently
              configure.outStream << "WARNING: Could Not Fill Segment #" << Draft().NumSegments()
                                  << " from file." << std::endl;
              this->DeleteLastSegment();
            } else if (Draft().GetSegment(Draft().NumSegments())->NumPoints() == 0) {
              configure.outStream << "WARNING: Segment #" << numTotalSegments
                                  << " is empty and will not be used." << std::endl;
              this->DeleteLastSegment();
            } else {
              // Handle advanced segments - parse components and add them to the segment
              // This must be done AFTER the segment is filled with data points
              if (NewSegment.IsAdvanced()) {
                // Set the operation type from the segment data
                OperationType operation = (NewSegment.GetLegacyOperationType() == 0) ? SUM : RATIO;
                Segment *filledSegment = Draft().GetSegment(Draft().NumSegments());
                filledSegment->SetOperationType(operation);

                // One leaf per component (spec: entrance, exit, optional fixed angle, scaling)
                for (const azure2::experiment::ComponentSpec &c : NewSegment.GetComponents()) {
                  if (!this->CheckComponentPairs(theCNuc, *filledSegment, c.entranceKey, c.exitKey, configure)) return -1;
                  this->AddLeaf(*filledSegment, c.entranceKey, c.exitKey,
                                c.hasAngle ? std::optional<double>(c.angle) : std::nullopt, c.scaling);
                }
              }

              int thisSegmentNum = Draft().NumSegments();
              if (Draft().GetSegment(thisSegmentNum)->IsTotalCapture()) {
                Segment *baseSegment = Draft().GetSegment(thisSegmentNum);
                int numCapturePairs = 0;

                // Set operation type to SUM for total capture
                baseSegment->SetOperationType(SUM);

                const int entranceNucleus = theCNuc->GetPair(theCNuc->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == entranceNucleus) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(*this, theCNuc->GetPair(i)->GetPairKey());
                    } else {
                      // Each further capture pair is a leaf of the segment
                      this->AddLeaf(*baseSegment, baseSegment->GetEntranceKey(), theCNuc->GetPair(i)->GetPairKey(),
                                    std::nullopt, 1.0);
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
              if (!this->FinishTransformedSegment(theCNuc, Draft().GetSegment(thisSegmentNum), numTotalSegments, false, configure))
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

  if (CheckThmWeights(configure, configure.thm.weightBySegment, "weight", this) != 0)
    return -1;
  if (BuildThmGroups(configure, theCNuc) != 0) return -1;

  if (Draft().NumSegments() > 0) {
    if (this->SetUpTargetEffects(configure, theCNuc, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(theCNuc, configure)) return -1;
    this->MapData(theCNuc);
  }

  return 0;
}

/*!
 * If the AZURE calculation is not data driven, this function is called in place
 * of the WorkingCopy::Fill function to create points at specified energies and angles.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int WorkingCopy::MakePoints(const Config &configure, CNuc *theCNuc, const azure2::experiment::ExtrapolationSpecs &specs) {
  const int numTotalSegments_ = (int)specs.extrapolations.items.size();
  for (const auto &item : specs.extrapolations.items) {
    {
      const azure2::experiment::ExtrapolationSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        Segment NewSegment(segment);
        if (theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
          PPair *segEntrance = theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          if (!AnalyzingPowerAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!NaturalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!BetaDelayedAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!UPOSAllowed(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!ThmVertexAllowsPair(theCNuc, NewSegment, numTotalSegments, true, configure)) return -1;
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
            Segment *theSegment = Draft().GetSegment(Draft().NumSegments());

            theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()))->SetEntrance();
            PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()));
            PPair *exitPair = theCNuc->GetPairFromKey(theSegment->GetExitKey());  // nullptr for total capture (-1)
            theSegment->SetGridRows(*this, segment.energies, segment.angles);
            for (int i = 0; i < theSegment->NumPoints(); i++) {
              // A ratio is dimensionless: no cross-section conversion (as Segment::Fill for data segments;
              // the stored factor would otherwise divide the ratio in the lab output, kinflag.dat = 1).
              azure2::experiment::AngleLimit lim;
              if (!theSegment->ConvertRowToCM(*this, 0, i, entrancePair, exitPair, !theSegment->IsDimensionlessRatio(),
                                              Segment::BeyondRange::Stop, lim)) {
                const Unit *thePoint = theSegment->GetPoint(*this, i + 1);
                configure.outStream << "ERROR: " << "test segment #" << theSegment->GetSegmentKey() << ": lab angle " << thePoint->GetLabAngle() << " deg at E_lab = " << thePoint->GetLabEnergy()
                                    << " MeV: " << lim.why << "." << std::endl;
                return -1;
              }
            }
            theSegment->KeepLoadFrame(*this, 0);

            // Handle advanced segments (sum/ratio of components) - AFTER points are added
            if (NewSegment.IsAdvanced()) {
              int operation = NewSegment.GetOperationType();
              theSegment->SetOperationType((OperationType)operation);

              // One leaf per component (spec: entrance, exit, optional fixed angle, scaling)
              for (const azure2::experiment::ComponentSpec &c : NewSegment.GetComponents()) {
                if (!this->CheckComponentPairs(theCNuc, *theSegment, c.entranceKey, c.exitKey, configure)) return -1;
                this->AddLeaf(*theSegment, c.entranceKey, c.exitKey,
                              c.hasAngle ? std::optional<double>(c.angle) : std::nullopt, c.scaling);
              }
            }
            if (theSegment->NumPoints() == 0) {
              configure.outStream << "WARNING: Extrapolation segment #" << numTotalSegments
                                  << " is empty and will not be used." << std::endl;
              this->DeleteLastSegment();
            } else {
              int thisSegmentNum = Draft().NumSegments();
              if (Draft().GetSegment(thisSegmentNum)->IsTotalCapture()) {
                Segment *baseSegment = Draft().GetSegment(thisSegmentNum);
                int numCapturePairs = 0;

                // Set operation type to SUM for total capture
                baseSegment->SetOperationType(SUM);

                const int entranceNucleus = theCNuc->GetPair(theCNuc->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= theCNuc->NumPairs(); i++) {
                  if (theCNuc->GetPair(i)->GetPType() == 10 && theCNuc->GetPair(i)->GetNucleus() == entranceNucleus) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(*this, theCNuc->GetPair(i)->GetPairKey());
                    } else {
                      // Each further capture pair is a leaf of the segment
                      this->AddLeaf(*baseSegment, baseSegment->GetEntranceKey(), theCNuc->GetPair(i)->GetPairKey(),
                                    std::nullopt, 1.0);
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
              if (!this->FinishTransformedSegment(theCNuc, Draft().GetSegment(thisSegmentNum), numTotalSegments, true, configure))
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

  if (CheckThmWeights(configure, configure.thm.weightByTestSegment, "weightTest", this) != 0)
    return -1;

  if (Draft().NumSegments() > 0) {
    if (this->SetUpTargetEffects(configure, theCNuc, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(theCNuc, configure)) return -1;
    this->MapData(theCNuc);
  }

  return 0;
}

/*!
 * Returns true if the data is to be fit, otherwise returns false.  Used in the AZURECalc function class
 * to determine if a clone of the CNuc and WorkingCopy objects should be made for thread safety.
 */

bool WorkingCopy::IsFit() const {
  return isFit_;
}

/*!
 * Returns true if the call to function is for error analysis via Minos, otherwise returns false.  Used in the AZURECalc function
 * class to suppress transformation and file output during error analysis.
 */

bool WorkingCopy::IsErrorAnalysis() const {
  return isErrorAnalysis_;
}

/*!
 * Sets the boolean indicating if the data is to be fit by AZURECalc function class. Used in AZUREMain function
 * class before calls to Minuit and AZURECalc.
 */

/*!
 * Sets an internal variable specifying if the data is to be fit by Minuit.  Needed to determine cloning behavior in AZURECalc for
 * thread safety.
 */

void WorkingCopy::SetFit(bool fit) {
  isFit_ = fit;
}

/*!
 * Sets the boolean indicating if the call to the function is for error analysis via Minos.
 */

void WorkingCopy::SetErrorAnalysis(bool errorAnalysis) {
  isErrorAnalysis_ = errorAnalysis;
}

/*!
 * This function counts one more evaluation of the objective and returns the count before it. The read and the
 * increment are one atomic operation, so evaluations made concurrently from the threads of Minuit2's numerical
 * gradient each receive a distinct number. Relaxed order suffices: the counter orders no other memory.
 */

int WorkingCopy::Iterate() {
  return iterations_.n.fetch_add(1, std::memory_order_relaxed);
}

/*!
 * This function is identical in role to the Unit::Initialize function, except that it initializes
 * and entire WorkingCopy object instead of a single Unit object.
 */

int WorkingCopy::Initialize(CNuc *compound, const Config &configure) {
  // The last step of the set-up: the data set is written once more (the nuclei, the leaves' frames at load) and then
  // frozen.  A second Initialize (a session set up again) works on a new data set.
  if (!draft_) {
    draft_ = std::make_shared<DataSet>(*set_);
    set_ = draft_;
  }
  // Labels of the compound nuclei, for the sample decomposition files.
  Draft().nucleusLabels_.clear();
  for (int k = 1; k <= compound->NumNuclei(); k++) Draft().nucleusLabels_.push_back(compound->NucleusLabel(k));
  Draft().numNuclei_ = compound->NumNuclei();
  // Identical-particle phase-shift segment validation: warn once per segment
  // if the requested L is forbidden by Bose/Fermi symmetry on an identical
  // pair. The phase branch in GenMatrixFunc would silently report delta = 0
  // for these (no matching channel), which is easy to misread as "no nuclear
  // interaction" rather than "this partial wave does not exist".
  for (int s = 1; s <= Draft().NumSegments(); s++) {
    const Segment *seg = Draft().GetSegment(s);
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
  for (int s = 1; s <= Draft().NumSegments(); s++) {
    const Segment *seg = Draft().GetSegment(s);
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

  // The leaves: their pairs, frames and tables
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Initializing Component Segments..." << std::endl;
  if (this->InitializeLeaves(compound, configure) == -1) return -1;
  // Every point of this stage has its rows now; nothing reads the channel table until the energies change, and the
  // table (about 64 bytes a row) is not kept beside the rows it filled (docs/upgrade/channel_table.md section 2).
  channelTable_.reset();
  // Where the result record of every point and sub-point lies (Results.h); laid out again when the points change
  // (BeginEvaluation).
  this->LayoutResults(nullptr);
  this->Freeze();

  return 0;
}

/*!
 * Adds a segment to the Segment vector.
 */

void WorkingCopy::AddSegment(Segment segment) {
  // The segment's state in this copy starts at the segment's start values, with its own term and no points yet.
  SegmentState state;
  state.norm = segment.StartNorm();
  state.energyShift = segment.StartEnergyShift();
  state.energyShiftSqrt = segment.StartEnergyShiftSqrt();
  state.thickness = segment.StartThickness();
  state.terms.resize(1);
  segment.SetPos((int)Draft().segments_.size());
  Draft().segments_.push_back(std::move(segment));
  segmentStates_.push_back(std::move(state));
}

// A copy shares the data set; one copied while the set-up still writes it gets the data set as it stands, frozen.
WorkingCopy::WorkingCopy(const WorkingCopy &other)
    : set_(other.draft_ ? std::make_shared<const DataSet>(*other.draft_) : other.set_),
      segmentStates_(other.segmentStates_),
      effectValues_(other.effectValues_),
      thmCoherentValues_(other.thmCoherentValues_),
      units_(other.units_),
      pendingSubPoints_(other.pendingSubPoints_),
      plan_(other.plan_),
      covFiles_(other.covFiles_),
      covarianceError_(other.covarianceError_),
      kinematicsError_(other.kinematicsError_),
      iterations_(other.iterations_),
      isFit_(other.isFit_),
      isErrorAnalysis_(other.isErrorAnalysis_),
      ecReadPos_(other.ecReadPos_),
      ecUsePrevious_(other.ecUsePrevious_),
      ecSignature_(other.ecSignature_),
      ecOutputFile_(other.ecOutputFile_),
      resultsSize_(other.resultsSize_),
      modelSize_(other.modelSize_),
      gridNuc_(other.gridNuc_),
      nodeCache_(other.nodeCache_),
      channelTable_(other.channelTable_) {}

DataSet &WorkingCopy::Draft() {
  if (!draft_) throw std::logic_error("WorkingCopy::Draft: the data set is frozen after the set-up (ReviseDataSet changes it)");
  return *draft_;
}

void WorkingCopy::SetPriorCovariance(std::shared_ptr<const PriorCovariance> prior) {
  ReviseDataSet([&](DataSet &d) { d.priorCov_ = std::move(prior); });
}

/*!
 * Calls Unit::CalcLegendreP for each point in the entire WorkingCopy object.
 */

void WorkingCopy::CalcLegendreP(int maxL, CNuc *theCNuc) {
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    const TargetEffect *effect = (segment->IsTargetEffect() &&
                            Data().GetTargetEffect(segment->GetTargetEffectNum())->IsQCoefficients())
        ? Data().GetTargetEffect(segment->GetTargetEffectNum())
        : NULL;
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) {
      Unit *point = segment->GetPoint(*this, i);
      point->CalcLegendreP(this, maxL, theCNuc, effect);
    }
  }
}

/*!
 * Calls Unit::CalcEDependentValues for each point in the entire WorkingCopy object.
 */

/*!
 * The channel table of the given points and their sub-points, built before the threads that fill them start
 * (docs/upgrade/channel_table.md).  Rows the current table holds are copied from it.
 */
std::shared_ptr<const ChannelTable> WorkingCopy::BuildChannelTable(CNuc *theCNuc, const Config &configure,
                                                             const std::vector<Unit *> &points) {
  // Each point's keys are made unique before they are added: its slots repeat (pair, l, E_c) once per channel spin
  // and J^pi block, and the full list was the peak memory of the build (natural_neon: 11 M keys).
  std::vector<ChannelTable::Key> keys, own;
  std::size_t requested = 0;
  for (Unit *point : points) {
    own.clear();
    point->ChannelKeys(this, theCNuc, own);
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

int WorkingCopy::CalcEDependentValues(CNuc *theCNuc, const Config &configure) {
  {
    std::vector<Unit *> points;
    for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++)
      for (int i = 1; i <= segment->NumPoints(); i++)
        if (!segment->GetPoint(*this, i)->IsMapped()) points.push_back(segment->GetPoint(*this, i));
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
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    bool localStop = false;
#pragma omp parallel for shared(localStop, configure)
    for (int i = 1; i <= segment->NumPoints(); i++) {
      if (configure.stopFlag || localStop) continue;
      Unit *point = segment->GetPoint(*this, i);
      if (!(point->IsMapped())) {
        try {
          point->CalcEDependentValues(this, theCNuc, configure, true, this, table);
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
 * Calls Unit::CalcCoulombAmplitude for each point in the entire WorkingCopy object.
 */

void WorkingCopy::CalcCoulombAmplitude(CNuc *theCNuc) {
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) {
      Unit *point = segment->GetPoint(*this, i);
      point->CalcCoulombAmplitude(this, theCNuc);
    }
  }
}

/*!
 * Sets up the leaves of every segment (Segment::AddLeaf: the components of SUM, RATIO, sample and total-capture
 * segments, the total term of a yield): their pairs in the compound nucleus, their frames with their own pairs, their
 * energy-dependent values, Legendre polynomials, Coulomb amplitudes and external-capture amplitudes (appended to the
 * intEC file after those of the segments).  The run log still calls them component segments.
 */

int WorkingCopy::InitializeLeaves(CNuc *theCNuc, const Config &configure) {
  // The leaves go through the steps of the segments' own points
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Initializing Component Segments (" << Draft().leafRefs_.size() << " components)..." << std::endl;

  // First, ensure entrance/exit pairs are properly set up in the compound nucleus
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    int entranceKey = Draft().segments_[leaf.first].TermEntranceKey(leaf.second);
    int exitKey = Draft().segments_[leaf.first].TermExitKey(leaf.second);

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

  // The frame of each leaf with its own pairs
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    int entranceKey = Draft().segments_[leaf.first].TermEntranceKey(leaf.second);
    int exitKey = Draft().segments_[leaf.first].TermExitKey(leaf.second);

    PPair *entrancePair = nullptr;
    PPair *exitPair = nullptr;

    if (theCNuc->IsPairKey(entranceKey)) {
      entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(entranceKey));
    }
    if (theCNuc->IsPairKey(exitKey)) {
      exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(exitKey));
    }

    // Apply the same conversions that are done for test segment points
    for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++) {
      Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);

      if (entrancePair && exitPair) {
        // The base segment's sequence, with the component's own pairs; a component's value is always converted.
        azure2::experiment::AngleLimit lim;
        if (!Draft().segments_[leaf.first].ConvertRowToCM(*this, leaf.second, i - 1, entrancePair, exitPair, true,
                                                          Segment::BeyondRange::Stop, lim)) {
          configure.outStream << "ERROR: " << "component (" << entranceKey << " -> " << exitKey << ") of segment #" << Draft().segments_[leaf.first].TermKey(leaf.second) << ": lab angle " << point->GetLabAngle() << " deg at E_lab = " << point->GetLabEnergy()
                              << " MeV: " << lim.why << "." << std::endl;
          return -1;
        }
        // The sub-points of a target effect were laid out (SetUpTargetEffects) from the copied base
        // point's c.m. angle, before this conversion with the component's own pair; put them at the
        // converted angle (CalcLegendreP below recurses into them).  Until 2026-10-27 a component with a
        // different exit pair, or a fixed-angle RATIO denominator (raw lab angle), was integrated at the
        // wrong angle (audit 2026-10-25, item 32).
        for (int k = 1; k <= point->NumSubPoints(); k++) {
          point->GetSubPoint(*this, k)->SetLabAngle(point->GetLabAngle());
          point->GetSubPoint(*this, k)->SetCMAngle(point->GetCMAngle());
        }
      }
    }
    if (entrancePair && exitPair) Draft().segments_[leaf.first].KeepLoadFrame(*this, leaf.second);
  }

  // CalcEDependentValues for the leaves, in parallel as for the data segments
  // (WorkingCopy::CalcEDependentValues).  Each point writes only itself, its sub-points and its local
  // mapped points that do not own their values; those are never computed here (the condition
  // below) and each has exactly one host, so no two threads write the same point.  Serial, this
  // loop was 9 of the 17.6 s of a natural_neon calculation on 24 threads.
  {
    std::vector<Unit *> points;
    for (const std::pair<int, int> &leaf : Draft().leafRefs_)
      for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++)
        if (!(Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i)->IsMapped()) || Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i)->OwnsEDependentValues())
          points.push_back(Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i));
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
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    bool localStop = false;
#pragma omp parallel for shared(localStop, configure)
    for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++) {
      if (configure.stopFlag || localStop) continue;
      Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);
      if (!(point->IsMapped()) || point->OwnsEDependentValues()) {
        try {
          point->CalcEDependentValues(this, theCNuc, configure, true, this, componentTable);
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

  // CalcLegendreP for the leaves
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    // Get TargetEffect with Q-coefficients if applicable
    TargetEffect *effect = NULL;
    if (Draft().segments_[leaf.first].IsTargetEffect()) {
      TargetEffect *te = Draft().GetTargetEffect(Draft().segments_[leaf.first].TermTargetEffectNum(leaf.second));
      if (te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++) {
      Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);
      point->CalcLegendreP(this, configure.maxLOrder, theCNuc, effect);
    }
  }

  // CalcCoulombAmplitude for the leaves
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++) {
      Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);
      point->CalcCoulombAmplitude(this, theCNuc);
    }
  }

  // Calculate EC amplitudes for the leaves if external capture is enabled.
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
    if (this->LeafECAmplitudes(theCNuc, configure, usePrevious, in, out) == -1) return -1;
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
 * Checks a component's pair keys before the component is created: both must
 * exist and belong to one compound nucleus (a reaction does not change the
 * compound nucleus), and the component's entrance pair is flagged as an
 * entrance HERE, while CNuc::Initialize -- which builds the reaction pathways
 * of every entrance pair -- still lies ahead.  Flagging it only in
 * WorkingCopy::Initialize, as was done, left a pair used by components alone (the
 * second isotope of a natural target) without pathways: GetDecay(0).
 * Returns false, with the error printed, when the component cannot be made.
 */
bool WorkingCopy::CheckComponentPairs(CNuc *theCNuc, const Segment &baseSegment, int entranceKey, int exitKey,
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
 * A transmission or yield segment (Segment::IsTransformed), once its components exist: refuse what the transform
 * cannot mean, make every component a plain cross section, and give a yield its hidden total-cross-section term
 * (entrance pair to itself, "total 1", by the optical theorem), the s_t of Y0 = (s / s_t)(1 - exp(-n s_t)).
 */
bool WorkingCopy::FinishTransformedSegment(CNuc *theCNuc, Segment *seg, int key, bool isTest, const Config &configure) {
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
  if (seg->IsTransmission()) return true;

  if (seg->IsNatural())
    return refuse("not on a sample (abundance-weighted) segment: the total of a mixture is not implemented for Y0.");
  if (seg->GetExitKey() == seg->GetEntranceKey() || seg->IsTotal())
    return refuse("Y0 is a capture or reaction yield: the exit pair must differ from the entrance pair (no \"total 1\").");
  for (int t = 1; t < seg->NumTerms(); t++) {
    if (seg->TermEntranceKey(t) != seg->GetEntranceKey() || seg->TermExitKey(t) == seg->TermEntranceKey(t))
      return refuse("every term of a yield must be a reaction of the segment's entrance pair (pair " +
                    std::to_string(seg->GetEntranceKey()) + ") to another exit pair.");
  }
  if (entrance->IsIdentical()) return refuse("not for a pair of identical particles.");
  if (configure.paramMask & Config::USE_RMC_FORMALISM)
    return refuse("the total cross section of Y0 needs the optical theorem, not available with --use-rmc.");
  if (!seg->HasComponents()) seg->SetOperationType(SUM);
  seg->MakeYieldTotalLeaf(*this, this->AddLeaf(*seg, seg->GetEntranceKey(), seg->GetEntranceKey(), std::nullopt, 1.0));
  return true;
}

/*!
 * After the target effects are read: a transmission or yield segment whose points carry a sub-point effect sums its
 * terms at every sub-point before the transform, so every term must sit on the base point's sub-point grid.  Terms of
 * one entrance pair do (their grids are laid out alike); the terms of a sample with several compound nuclei do not,
 * and an effect blended in at range edges or skipped by an automatic tolerance is not one kernel on transformed
 * values.  Refused here rather than computed wrongly.
 */
bool WorkingCopy::CheckTransformedSegments(const Config &configure) {
  for (const Segment &seg : Draft().GetSegments()) {
    if (!seg.IsTransformed()) continue;
    const std::string what = "segment #" + std::to_string(seg.GetSegmentKey()) + " (" +
                             (seg.IsTransmission() ? "transmission" : "yield") + "): ";
    for (int p = 1; p <= seg.NumPoints(); p++) {
      Unit *base = seg.GetPoint(*this, p);
      if (!Segment::HasSubPointIntegration(*this, base)) {
        for (int t = 1; t < seg.NumTerms(); t++)
          if (p <= (int)seg.TermPoints(*this, t).size() && Segment::HasSubPointIntegration(*this, seg.GetTermPoint(*this, t, p))) {
            configure.outStream << "ERROR: " << what << "a term is integrated over sub-points at E_lab "
                                << base->GetLabEnergy() << " MeV where the segment's own point is not." << std::endl;
            return false;
          }
        continue;
      }
      TargetEffect *te = Draft().GetTargetEffect(base->GetTargetEffectNum());
      // A thin-target integration blends with sigma(E_bar) inside its kernel (Unit::TargetFallbackImpl), which is
      // still one linear kernel on the transformed values; other kernels blend with a plain point outside it.
      const bool thinTarget = te && te->IsTargetIntegration() && !te->IsThickTarget();
      if ((base->GetTargetBlendWeight() < 1.0 && !thinTarget) || (te && te->GetAutoTolerance() > 0.0)) {
        configure.outStream << "ERROR: " << what << "its target effect blends in at range edges or has an automatic "
                            << "tolerance; the transform needs the full kernel at every point (drop the blend width "
                            << "and the auto_tolerance of its fold_windows)." << std::endl;
        return false;
      }
      for (int t = 1; t < seg.NumTerms(); t++) {
        Unit *cp = p <= (int)seg.TermPoints(*this, t).size() ? seg.GetTermPoint(*this, t, p) : nullptr;
        bool same = cp && cp->NumSubPoints() == base->NumSubPoints();
        for (int j = 1; same && j <= base->NumSubPoints(); j++) {
          const double a = base->GetSubPoint(*this, j)->GetCMEnergy(), b = cp->GetSubPoint(*this, j)->GetCMEnergy();
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
 * A leaf of a SUM, RATIO, sample, total-capture or yield segment (Segment::AddLeaf): its pairs, a fixed angle (deg;
 * a RATIO denominator at one angle) or none, its scale, and a key of its own from the data set's counter (negative:
 * not a segment of the file; the caches and the intEC file are keyed by it).  Returns its term.
 */
int WorkingCopy::AddLeaf(Segment &base, int entranceKey, int exitKey, std::optional<double> fixedAngle, double scaling) {
  const int t = base.AddLeaf(*this, entranceKey, exitKey, fixedAngle, scaling, Draft().nextComponentKey_--);
  Draft().leafRefs_.push_back({(int)(&base - Draft().segments_.data()), t});
  return t;
}

namespace {

// True when the calculation of `b` can be taken from that of `a`: both are at
// the same energy, with the same entrance and exit pair, so they share the
// T-matrix of the point and -- checked here, sub-point by sub-point -- of the
// integration grid.  Only what enters through GenMatrixFunc::
// CalculateCrossSection on the point itself (angle, Coulomb amplitude,
// observable type) may differ.  Whether `a` and `b` may sit in different
// segments is the caller's question (ShiftsAgree).
bool SharesCalculation(WorkingCopy &data, Unit &a, Unit &b, bool sameExit = true) {
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
  const UnitRange<Unit> sa = a.GetSubPoints(data);
  const UnitRange<Unit> sb = b.GetSubPoints(data);
  for (int i = 0; i < sa.size(); i++) {
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

ShiftForm ShiftFormOf(const WorkingCopy &data, int pos) {
  const std::vector<Segment> &segments = data.Data().GetSegments();
  const Segment &seg = segments[pos];
  const SegmentState &state = data.State(seg);
  ShiftForm f;
  const int lead = seg.GetShiftLeaderPos() >= 0 ? seg.GetShiftLeaderPos() : pos;
  if (segments[lead].IsVaryEnergyShift()) f.constOwner = lead;
  f.sqrtFree = seg.IsVaryEnergyShiftSqrt();
  f.constShift = state.energyShift;
  f.sqrtShift = state.energyShiftSqrt;
  f.hasFn = seg.HasShiftFn();
  if (f.hasFn) {
    f.fnExpression = seg.GetShiftFnExpression();
    f.fnValues = state.shiftFnValues;
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
 * (Unit::Calculate).
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
 * sub-point grid rebuilds the mapping (WorkingCopy::FillEnergyShiftsFromParams).
 * Normalizations play no part: they scale the data, not the model.
 *
 * The host of a group is its first point in data order (segment, then point
 * index): every loop that calculates the whole data set walks it in that order,
 * so a host is calculated before the points that read from it
 * (azure2::EvaluateModelPoint; Segment::
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

void WorkingCopy::MapData(CNuc *theCNuc) {
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++)
      point->ClearMapping();
  }
  for (const std::pair<int, int> &leaf : Data().leafRefs_)
    for (Unit &point : Data().segments_[leaf.first].TermPoints(*this, leaf.second)) point.ClearMapping();
  // (host, mapped point) in the order each host calculates its mapped points
  std::vector<std::pair<azure2::PointRef, azure2::PointRef>> mapping;
  const char *off = std::getenv("AZURE_NO_POINT_MAPPING");
  if (off && *off && std::string(off) != "0") return plan_.Build(*this, mapping);
  const char *crossOff = std::getenv("AZURE_NO_CROSS_SEGMENT_MAPPING");
  const bool crossSegment = !(crossOff && *crossOff && std::string(crossOff) != "0");

  const std::vector<Segment> &segments = Data().GetSegments();
  const int numSegments = (int)segments.size();
  std::vector<ShiftForm> shiftForms(numSegments);
  for (int s = 0; s < numSegments; s++) shiftForms[s] = ShiftFormOf(*this, s);
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
  auto pointOf = [this, &segments](const Ref &r) -> Unit & { return segments[r.seg].GetPoints(*this)[r.point]; };

  // Group by pair and energy; a stable sort keeps data order inside a group,
  // so the first point of a group is its first point in data order.
  std::stable_sort(refs.begin(), refs.end(), [&](const Ref &x, const Ref &y) {
    Unit &a = pointOf(x), &b = pointOf(y);
    if (a.GetEntranceKey() != b.GetEntranceKey()) return a.GetEntranceKey() < b.GetEntranceKey();
    if (a.GetExitKey() != b.GetExitKey()) return a.GetExitKey() < b.GetExitKey();
    if (a.GetOriginalEnergy() != b.GetOriginalEnergy()) return a.GetOriginalEnergy() < b.GetOriginalEnergy();
    return a.GetCMEnergy() < b.GetCMEnergy();
  });

  std::vector<Ref> hosts;
  const int n = (int)refs.size();
  int start = 0;
  while (start < n) {
    Unit &first = pointOf(refs[start]);
    int end = start + 1;
    while (end < n) {
      Unit &e = pointOf(refs[end]);
      if (e.GetEntranceKey() != first.GetEntranceKey() || e.GetExitKey() != first.GetExitKey() ||
          e.GetOriginalEnergy() != first.GetOriginalEnergy() || e.GetCMEnergy() != first.GetCMEnergy())
        break;
      end++;
    }
    hosts.clear();
    for (int k = start; k < end; k++) {
      Unit &p = pointOf(refs[k]);
      bool mapped = false;
      for (const Ref &h : hosts) {
        if (h.seg != refs[k].seg && !(crossSegment && !ownOnly[h.seg] && !ownOnly[refs[k].seg] &&
                                      ShiftsAgree(shiftForms[h.seg], shiftForms[refs[k].seg])))
          continue;
        if (SharesCalculation(*this, pointOf(h), p)) {
          p.SetMap(h.seg + 1, h.point + 1);
          mapping.push_back({segments[h.seg].TermFirstUnit(*this, 0) + h.point,
                             segments[refs[k].seg].TermFirstUnit(*this, 0) + refs[k].point});
          mapped = true;
          break;
        }
      }
      if (!mapped) hosts.push_back(refs[k]);
    }
    start = end;
  }

  MapLeaves(theCNuc, mapping);
  plan_.Build(*this, mapping);
}

/*!
 * Point mapping inside a segment with components (SUM, RATIO, total capture): each component's point i is mapped onto
 * the base segment's point i.  A component has its own exit pair but the base's entrance pair, energies, sub-point
 * grid and energy shift, so its level matrix is the host's: the host's calculation builds only the component's
 * T-matrix and observable (Unit::Calculate).  Before 2026-10-24 every term was a calculation of its own, and a
 * neutron total cross section summed over nine exit pairs cost nine level-matrix solutions per point (13C+a: 87 % of
 * an evaluation).
 *
 * Such a family is not mapped with other segments.  A component re-energises its own points on an energy shift
 * (OwnsEDependentValues).  A beta-delayed entrance is left out: its energy-dependent values depend on the exit pair.
 * AZURE_NO_COMPONENT_MAPPING=1 switches this off (and AZURE_NO_POINT_MAPPING=1 all mapping).
 */
void WorkingCopy::MapLeaves(CNuc *theCNuc, std::vector<std::pair<azure2::PointRef, azure2::PointRef>> &mapping) {
  if (!theCNuc) return;
  const char *off = std::getenv("AZURE_NO_COMPONENT_MAPPING");
  if (off && *off && std::string(off) != "0") return;
  std::map<std::pair<int, int>, int> leafPos;  // (segment, term) -> position in the plan's leaf list
  for (int c = 0; c < (int)Data().leafRefs_.size(); c++) leafPos[Data().leafRefs_[c]] = c;
  const std::vector<Segment> &segments = Data().GetSegments();
  for (int s = 0; s < (int)segments.size(); s++) {
    const Segment &seg = segments[s];
    if (!seg.HasComponents() && !seg.IsTotalCapture()) continue;
    const int entrancePair = theCNuc->GetPairNumFromKey(seg.GetEntranceKey());
    if (entrancePair < 1 || theCNuc->GetPair(entrancePair)->GetPType() == 20) continue;
    const UnitRange<Unit> points = seg.GetPoints(*this);
    const int n = points.size();
    // Each component's point i onto the base's point i (components in turn, so a host's mapped points come grouped
    // by exit and Unit::Calculate builds each exit's T-matrix once).  The base's points stay unmapped among
    // themselves, as before: the external-capture amplitudes and their file (intEC.dat) skip mapped points, so
    // mapping them would change that file's layout.
    for (int t = 1; t < seg.NumTerms(); t++) {
      const UnitRange<Unit> leafPoints = seg.TermPoints(*this, t);
      if (leafPoints.size() != n || !leafPos.count({s, t})) continue;
      for (int p = 0; p < n; p++) {
        Unit *cp = &leafPoints[p];
        if (!SharesCalculation(*this, points[p], *cp, false)) continue;
        cp->SetMap(s + 1, p + 1);
        cp->SetOwnsEDependentValues(true);
        mapping.push_back({seg.TermFirstUnit(*this, 0) + p, seg.TermFirstUnit(*this, t) + p});
      }
    }
  }
}

/*!
 * Adds a TargetEffect object to the vector contained within the present object.
 */

void WorkingCopy::AddTargetEffect(TargetEffect targetEffect) {
  TargetEffectValues values;
  values.screeningUe = targetEffect.GetScreeningUe();
  values.density = targetEffect.GetDensity();
  values.dopplerTemperature = targetEffect.GetDopplerTemperature();
  values.stoppingScale = targetEffect.GetStoppingScaleNominal();
  Draft().targetEffects_.push_back(std::move(targetEffect));
  effectValues_.push_back(values);
}


/*!
 * Deletes the last segment from the segment vector.
 */

void WorkingCopy::DeleteLastSegment() {
  Draft().segments_.pop_back();
  segmentStates_.pop_back();
}

template <class F>
void WorkingCopy::ForEachResultUnit(F f) {
  const int nNatural = set_->numNuclei_;
  auto segment = [&](const Segment &seg, int t) {
    const int nat = (t == 0 && seg.IsNatural()) ? nNatural : 0;
    for (Unit &point : seg.TermPoints(*this, t)) {
      f(point, point.ResultAngularCapacity(), nat);
      for (Unit &sub : point.GetSubPoints(*this)) f(sub, sub.ResultAngularCapacity(), 0);
    }
  };
  for (const Segment &seg : set_->segments_) segment(seg, 0);
  for (const std::pair<int, int> &leaf : set_->leafRefs_) segment(set_->segments_[leaf.first], leaf.second);
}

void WorkingCopy::LayoutResults(azure2::Results *results) {
  size_t total = 0;
  ForEachResultUnit([&](Unit &, int nAng, int nNat) { total += azure2::Results::RecordSize(nAng, nNat); });
  std::vector<double> next(results ? total : 0, 0.0);
  size_t offset = 0;
  ForEachResultUnit([&](Unit &point, int nAng, int nNat) {
    if (results) {
      const size_t from = point.ResultSlot();
      const int fromAng = point.ResultAngularSlots(), fromNat = point.ResultNaturalSlots();
      if (from != Unit::kNoRecord && from + azure2::Results::RecordSize(fromAng, fromNat) <= results->slot.size())
        azure2::CopyRecord(results->slot.data() + from, fromAng, fromNat, next.data() + offset, nAng, nNat);
    }
    point.SetResultRecord(offset, nAng, nNat);
    offset += azure2::Results::RecordSize(nAng, nNat);
  });
  resultsSize_ = total;
  if (results) results->slot.swap(next);
  modelSize_ = 0;
  for (const Segment &seg : set_->segments_)
    for (Unit &point : seg.GetPoints(*this)) {
      point.SetModelRecord(modelSize_, point.ResultAngularCapacity());
      modelSize_ += azure2::Results::ModelRecordSize(point.ResultAngularCapacity());
    }
}

void WorkingCopy::BeginEvaluation(azure2::Results &results) {
  size_t offset = 0;
  bool laidOut = true;
  ForEachResultUnit([&](Unit &point, int nAng, int nNat) {
    if (!laidOut) return;
    laidOut = point.ResultSlot() == offset && point.ResultAngularSlots() == nAng && point.ResultNaturalSlots() == nNat;
    offset += azure2::Results::RecordSize(nAng, nNat);
  });
  if (!laidOut || offset != resultsSize_) this->LayoutResults(&results);
  this->SizeResults(results);
  if (azure2::Results::PoisonRequested()) results.Poison();
}

/*!
 * Creates a new copy of the WorkingCopy object in memory and returns a pointer to the new object.
 * Used in AZURECalc function class for thread safety.
 */

WorkingCopy *WorkingCopy::Clone() const {
  // A plain copy: the units, the plan and the points' task indices are values and positions, valid in the copy as
  // they are (no unit points to another or to its data set); the measurements, the term settings, the covariance
  // groups, the node tables and the prior covariance are shared and immutable (or, the two caches, shared on
  // purpose).  A copy has no grid nucleus of its own (GridNucleus).
  WorkingCopy *dataCopy = new WorkingCopy(*this);
  // The copy's points (not their sub-points) start without the external-capture amplitudes of the set-up: a working
  // copy evaluates them at its points' energies (Unit::GetECAmplitudeWithShift).
  for (azure2::PointRef u : dataCopy->PlanUnits()) dataCopy->units_[u].ClearECAmplitudes();
  return dataCopy;
}
