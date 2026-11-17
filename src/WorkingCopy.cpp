#include "InputFiles.h"
#include "AZUREOutput.h"
#include "Profile.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "WorkingCopy.h"
#include "azure2/support/Environment.h"
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
#include <cstdlib>
#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <sstream>
#include "azure2/support/FloatingPointPrecision.h"
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
 * other, and the capture route (CompoundSystem::CalcCaptureAnalyzingPower) ran for any
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
static bool CoefficientsAllowed(const CompoundSystem *compound, const Segment &segment, int number, const Config &configure) {
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

static bool AnalyzingPowerAllowed(const CompoundSystem *compound, const Segment &segment, int number, bool isTest,
                                  const Config &configure) {
  if (!segment.IsAnalyzingPower()) return true;
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(segment.GetEntranceKey()));
  if (const azure2::physics::SpinObservable *observable = segment.GetSpinObservable()) {
    // A general spin observable: particle pairs on both sides, and every particle its operator names with the spin
    // the operator needs (physics/SpinObservable.h, Check).
    const char *what = isTest ? "test segment #" : "segment #";
    const PPair *exitPair = compound->IsPairKey(segment.GetExitKey())
                                ? compound->GetPair(compound->GetPairNumFromKey(segment.GetExitKey()))
                                : nullptr;
    std::string why;
    if (!exitPair)
      why = "needs one entrance and one exit pair";
    else if (entrancePair->GetPType() != 0 || exitPair->GetPType() != 0)
      why = "needs particle pairs on both sides (for capture only A_y, quantity analyzing_power, is computed)";
    else
      observable->Check(entrancePair->GetJ(1), entrancePair->GetJ(2), exitPair->GetJ(1), exitPair->GetJ(2), why);
    if (why.empty()) return true;
    configure.outStream << "ERROR: " << what << number << ": spin observable \"" << observable->token << "\": " << why
                        << "." << std::endl;
    return false;
  }
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
 * The UPOS coefficient (CompoundSystem, Brune and deBoer 2020 Eq. 22) is written for a secondary gamma ray from the HEAVY
 * partner of the exit pair: its final-state spin is J2.  With the light particle excited (column 32 of the
 * <levels> line, 2026-11-05) the gamma ray would come from the light particle, which that formula does not
 * describe -- refuse instead of computing the heavy partner's correlation silently.
 */
static bool UPOSAllowed(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !compound->IsPairKey(segment.GetExitKey())) return true;
  const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(segment.GetExitKey()));
  if (exitPair->GetLightExE() == 0.) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << key << " (UPOS): exit pair "
                      << segment.GetExitKey() << " has an excited light particle (column 32, " << exitPair->GetLightExE()
                      << " MeV); the secondary-gamma correlation is implemented for the decay of the heavy partner only."
                      << std::endl;
  return false;
}

static bool TotalAllowed(const CompoundSystem *compound, Segment &segment, int key, bool isTest, const Config &configure) {
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
  if (!compound->IsPairKey(segment.GetEntranceKey())) return true;  // reported where the pair is checked
  const PPair *pair = compound->GetPair(compound->GetPairNumFromKey(segment.GetEntranceKey()));
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
static bool BetaDelayedAllowed(const CompoundSystem *compound, Segment &segment, int key, bool isTest, const Config &configure) {
  if (!compound->IsPairKey(segment.GetEntranceKey()) || segment.GetExitKey() != -1) return true;
  if (compound->GetPair(compound->GetPairNumFromKey(segment.GetEntranceKey()))->GetPType() != 20) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment" : "segment") << " #" << key << ": a beta-delayed "
                      << "spectrum (entrance pair " << segment.GetEntranceKey() << " of type 20) needs its decay pair as"
                      << " the exit pair, not -1 (total capture)." << std::endl;
  return false;
}

static bool NaturalAllowed(const AbundanceSamples &samples, Segment &segment, int key, bool isTest, const Config &configure) {
  const char *what = isTest ? "test segment" : "segment";
  const std::string &token = segment.GetNaturalToken();
  if (token == "0") return true;   // not a natural segment
  if (!samples.Present()) {
    configure.outStream << "ERROR: " << what << " #" << key << " is a sample (abundance-weighted) segment, but the project has no "
                        << "<abundances> block (one line per compound nucleus: [<sample>:] <Z> <A> <abundance> <vary> <sigma>)." << std::endl;
    return false;
  }
  // "1" = the unnamed sample (or the only one); anything else names a sample of the block.
  int sample = (token == "1") ? samples.Default() : samples.Find(token);
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

int WorkingCopy::Fill(const Config &configure, const CompoundSystem *compound, const LevelInputs &inputs, const azure2::experiment::DataSpecs &specs) {
  const int numTotalSegments_ = (int)specs.datasets.items.size();
  for (const auto &item : specs.datasets.items) {
    {
      const azure2::experiment::DatasetSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        Segment NewSegment(segment);

        if (compound->IsPairKey(NewSegment.GetEntranceKey())) {
          const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          if (!AnalyzingPowerAllowed(compound, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!CoefficientsAllowed(compound, NewSegment, numTotalSegments, configure)) return -1;
          if (!NaturalAllowed(Data().Abundances(), NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalAllowed(compound, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!BetaDelayedAllowed(compound, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!UPOSAllowed(compound, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!ThmVertexAllowsPair(compound, NewSegment, numTotalSegments, false, configure)) return -1;
          if (!TotalCaptureCodeOk(NewSegment, numTotalSegments, false, configure)) return -1;
          UPOSLabNote(NewSegment, numTotalSegments, false, configure);
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            // total capture: every capture pair of the ENTRANCE pair's compound nucleus
            for (int i = 1; i <= compound->NumPairs(); i++) {
              if (compound->GetPair(i)->GetPType() == 10 && compound->GetPair(i)->GetNucleus() == entrancePair->GetNucleus()) {
                isValidTotal = true;
                break;
              }
            }
          }
          if (!isValidTotal && compound->IsPairKey(NewSegment.GetExitKey()) &&
              compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus() != entrancePair->GetNucleus()) {
            configure.outStream << "ERROR: segment #" << numTotalSegments << ": entrance pair " << NewSegment.GetEntranceKey()
                                << " and exit pair " << NewSegment.GetExitKey() << " belong to different compound nuclei ("
                                << compound->NucleusLabel(entrancePair->GetNucleus()) << " and "
                                << compound->NucleusLabel(compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus())
                                << "); a reaction stays within one nucleus." << std::endl;
            return -1;
          }
          if (isValidTotal || compound->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightBySegment, numTotalSegments));
            this->AddSegment(NewSegment);

            // Fill the segment with data first
            if (Draft().GetSegment(Draft().NumSegments())->Fill(compound, this, configure) == -1) {
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
                  if (!this->CheckComponentPairs(compound, *filledSegment, c.entranceKey, c.exitKey, configure)) return -1;
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

                const int entranceNucleus = compound->GetPair(compound->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= compound->NumPairs(); i++) {
                  if (compound->GetPair(i)->GetPType() == 10 && compound->GetPair(i)->GetNucleus() == entranceNucleus) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(*this, compound->GetPair(i)->GetPairKey());
                    } else {
                      // Each further capture pair is a leaf of the segment
                      this->AddLeaf(*baseSegment, baseSegment->GetEntranceKey(), compound->GetPair(i)->GetPairKey(),
                                    std::nullopt, 1.0);
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
              if (!this->FinishTransformedSegment(compound, Draft().GetSegment(thisSegmentNum), numTotalSegments, false, configure))
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
  if (BuildThmGroups(configure, compound) != 0) return -1;

  if (Draft().NumSegments() > 0) {
    if (this->SetUpTargetEffects(configure, compound, inputs, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(compound, configure)) return -1;
    this->BuildPlan();
  }

  return 0;
}

/*!
 * If the AZURE calculation is not data driven, this function is called in place
 * of the WorkingCopy::Fill function to create points at specified energies and angles.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int WorkingCopy::MakePoints(const Config &configure, const CompoundSystem *compound, const LevelInputs &inputs, const azure2::experiment::ExtrapolationSpecs &specs) {
  const int numTotalSegments_ = (int)specs.extrapolations.items.size();
  for (const auto &item : specs.extrapolations.items) {
    {
      const azure2::experiment::ExtrapolationSpec &segment = item.spec;
      const int numTotalSegments = segment.key;
      if (!item.warning.empty()) configure.outStream << "WARNING: " << item.warning << std::endl;
      if (!item.note.empty()) configure.Note(item.note);
      if (segment.active) {
        Segment NewSegment(segment);
        if (compound->IsPairKey(NewSegment.GetEntranceKey())) {
          const PPair *segEntrance = compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetEntranceKey()));
          if (!AnalyzingPowerAllowed(compound, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!NaturalAllowed(Data().Abundances(), NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalAllowed(compound, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!BetaDelayedAllowed(compound, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!UPOSAllowed(compound, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!ThmVertexAllowsPair(compound, NewSegment, numTotalSegments, true, configure)) return -1;
          if (!TotalCaptureCodeOk(NewSegment, numTotalSegments, true, configure)) return -1;
          UPOSLabNote(NewSegment, numTotalSegments, true, configure);
          bool isValidTotal = false;
          if (NewSegment.GetExitKey() == -1) {
            // total capture: every capture pair of the ENTRANCE pair's compound nucleus
            for (int i = 1; i <= compound->NumPairs(); i++) {
              if (compound->GetPair(i)->GetPType() == 10 && compound->GetPair(i)->GetNucleus() == segEntrance->GetNucleus()) {
                isValidTotal = true;
                break;
              }
            }
          }
          if (!isValidTotal && compound->IsPairKey(NewSegment.GetExitKey()) &&
              compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus() != segEntrance->GetNucleus()) {
            configure.outStream << "ERROR: test segment #" << numTotalSegments << ": entrance pair " << NewSegment.GetEntranceKey()
                                << " and exit pair " << NewSegment.GetExitKey() << " belong to different compound nuclei ("
                                << compound->NucleusLabel(segEntrance->GetNucleus()) << " and "
                                << compound->NucleusLabel(compound->GetPair(compound->GetPairNumFromKey(NewSegment.GetExitKey()))->GetNucleus())
                                << "); a reaction stays within one nucleus." << std::endl;
            return -1;
          }
          if (isValidTotal || compound->IsPairKey(NewSegment.GetExitKey())) {
            NewSegment.SetSegmentKey(numTotalSegments);
            if (NewSegment.IsTHM()) NewSegment.SetThmWeight(ThmWeightFor(configure.thm.weightByTestSegment, numTotalSegments));
            this->AddSegment(NewSegment);
            Segment *theSegment = Draft().GetSegment(Draft().NumSegments());

            const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(theSegment->GetEntranceKey()));
            const PPair *exitPair = compound->GetPairFromKey(theSegment->GetExitKey());  // nullptr for total capture (-1)
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
                if (!this->CheckComponentPairs(compound, *theSegment, c.entranceKey, c.exitKey, configure)) return -1;
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

                const int entranceNucleus = compound->GetPair(compound->GetPairNumFromKey(baseSegment->GetEntranceKey()))->GetNucleus();
                for (int i = 1; i <= compound->NumPairs(); i++) {
                  if (compound->GetPair(i)->GetPType() == 10 && compound->GetPair(i)->GetNucleus() == entranceNucleus) {
                    numCapturePairs++;
                    if (numCapturePairs == 1) {
                      // For the first capture pair, set the base segment's exit key
                      baseSegment->SetExitKey(*this, compound->GetPair(i)->GetPairKey());
                    } else {
                      // Each further capture pair is a leaf of the segment
                      this->AddLeaf(*baseSegment, baseSegment->GetEntranceKey(), compound->GetPair(i)->GetPairKey(),
                                    std::nullopt, 1.0);
                    }
                  }
                }
                baseSegment->SetIsTotalCapture(numCapturePairs);
              }
              if (!this->FinishTransformedSegment(compound, Draft().GetSegment(thisSegmentNum), numTotalSegments, true, configure))
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
    if (this->SetUpTargetEffects(configure, compound, inputs, specs.targetEffects) == -1) return -1;
    if (!this->CheckTransformedSegments(configure) || !this->BuildPredictions(compound, configure)) return -1;
    this->BuildPlan();
  }

  return 0;
}

/*!
 * Returns true if the data is to be fit, otherwise returns false.  Used in the AZURECalc function class
 * to determine if a clone of the CompoundSystem and WorkingCopy objects should be made for thread safety.
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

int WorkingCopy::Initialize(const CompoundSystem *compound, const Config &configure) {
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
    const PPair *pp = compound->GetPair(aa);
    if (!pp->IsIdentical()) continue;
    int segL = seg->GetL();
    bool allowedFound = false;
    for (int j = 1; j <= compound->NumJGroups() && !allowedFound; j++) {
      const JGroup *jg = compound->GetJGroup(j);
      if (jg->GetJ() != seg->GetJ()) continue;
      for (int ch = 1; ch <= jg->NumChannels() && !allowedFound; ch++) {
        const AChannel *channel = jg->GetChannel(ch);
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
      const JGroup *jg = compound->GetJGroup(j);
      if (jg->GetJ() != seg->GetJ()) continue;
      for (int ch = 1; ch <= jg->NumChannels(); ch++) {
        const AChannel *channel = jg->GetChannel(ch);
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
  // The readers of angle-free slots with one angle share their rows (RowSet.h).
  this->BuildRowSets(compound, configure);
  // Every point of this stage has its rows now; nothing reads the channel table until the energies change, and the
  // table (about 64 bytes a row) is not kept beside the rows it filled (docs/upgrade/channel_table.md section 2).
  channelTable_.reset();
  // Where the record of every point and sub-point lies (Results.h), sized for the compound's reactions.
  capacity_ = RecordCapacityOf(*compound);
  this->LayoutResults(nullptr);
  // The stage's node table, from the units as the set-up placed them (NodeTable.h).
  nodes_ = azure2::BuildNodeTable(*this, compound, configure);
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
      abundances_(other.abundances_),
      units_(other.units_),
      slots_(other.slots_),
      slotBlock_(other.slotBlock_),
      slotWriter_(other.slotWriter_),
      slotReaders_(other.slotReaders_),
      pendingBlocks_(other.pendingBlocks_),
      pendingOwnBlock_(other.pendingOwnBlock_),
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
      capacity_(other.capacity_),
      rowSets_(other.rowSets_),
      contractedSize_(other.contractedSize_),
      gridLevels_(other.gridLevels_),
      meshAnchors_(other.meshAnchors_),
      meshSupports_(other.meshSupports_),
      meshCentres_(other.meshCentres_),
      channelTable_(other.channelTable_), nodes_(other.nodes_) {}

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

void WorkingCopy::CalcLegendreP(int maxL, const CompoundSystem *compound) {
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    const TargetEffect *effect = (segment->IsTargetEffect() &&
                            Data().GetTargetEffect(segment->GetTargetEffectNum())->IsQCoefficients())
        ? Data().GetTargetEffect(segment->GetTargetEffectNum())
        : NULL;
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) { const azure2::support::ExtendedPrecisionScope fpPrecision;
      Unit *point = segment->GetPoint(*this, i);
      point->CalcLegendreP(this, maxL, compound, effect);
    }
  }
}

/*!
 * Calls Unit::CalcEDependentValues for each point in the entire WorkingCopy object.
 */

/*!
 * Calls Unit::CalcCoulombAmplitude for each point in the entire WorkingCopy object.
 */

void WorkingCopy::CalcCoulombAmplitude(const CompoundSystem *compound) {
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
#pragma omp parallel for
    for (int i = 1; i <= segment->NumPoints(); i++) { const azure2::support::ExtendedPrecisionScope fpPrecision;
      Unit *point = segment->GetPoint(*this, i);
      point->CalcCoulombAmplitude(this, compound);
    }
  }
}

/*!
 * Sets up the leaves of every segment (Segment::AddLeaf: the components of SUM, RATIO, sample and total-capture
 * segments, the total term of a yield): their pairs in the compound nucleus, their frames with their own pairs, their
 * energy-dependent values, Legendre polynomials, Coulomb amplitudes and external-capture amplitudes (appended to the
 * intEC file after those of the segments).  The run log still calls them component segments.
 */

int WorkingCopy::InitializeLeaves(const CompoundSystem *compound, const Config &configure) {
  // The leaves go through the steps of the segments' own points
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Initializing Component Segments (" << Draft().leafRefs_.size() << " components)..." << std::endl;

  // First, check the leaves' exit pairs in the compound nucleus (their entrance pairs were requested at set-up)
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    int exitKey = Draft().segments_[leaf.first].TermExitKey(leaf.second);

    // Check if exit pair is valid (either -1 for total capture or existing pair)
    bool isValidExit = false;
    if (exitKey == -1) {
      // Check for total capture validity
      for (int i = 1; i <= compound->NumPairs(); i++) {
        if (compound->GetPair(i)->GetPType() == 10) {
          isValidExit = true;
          break;
        }
      }
    } else {
      isValidExit = compound->IsPairKey(exitKey);
    }

    if (!isValidExit) {
      configure.outStream << "Warning: Invalid exit key " << exitKey
                          << " for component segment" << std::endl;
    }
  }

  // The frame of each leaf with its own pairs
  std::vector<int> angled;  // the points whose slots take their converted angles
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    int entranceKey = Draft().segments_[leaf.first].TermEntranceKey(leaf.second);
    int exitKey = Draft().segments_[leaf.first].TermExitKey(leaf.second);

    const PPair *entrancePair = nullptr;
    const PPair *exitPair = nullptr;

    if (compound->IsPairKey(entranceKey)) {
      entrancePair = compound->GetPair(compound->GetPairNumFromKey(entranceKey));
    }
    if (compound->IsPairKey(exitKey)) {
      exitPair = compound->GetPair(compound->GetPairNumFromKey(exitKey));
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
        // A slot another point of the group writes takes that point's angles; SeparateSlotsByAngle below gives
        // a point whose converted angles differ slots of its own.
        // Slots holding Legendre coefficients have no angle (Site::ValueFor): nothing to put.
        if (point->AngleFree()) continue;
        if (point->NumSubPoints()) angled.push_back(static_cast<int>(point - &UnitAt(0)));
        for (int k = 1; k <= point->NumSubPoints(); k++) {
          if (!WritesSlot(*point, k)) continue;
          point->GetSubPoint(*this, k)->SetLabAngle(point->GetLabAngle());
          point->GetSubPoint(*this, k)->SetCMAngle(point->GetCMAngle());
        }
      }
    }
    if (entrancePair && exitPair) Draft().segments_[leaf.first].KeepLoadFrame(*this, leaf.second);
  }
  this->SeparateSlotsByAngle(angled);

  // CalcEDependentValues for the leaves, in parallel as for the data segments (WorkingCopy::CalcEDependentValues).
  // Each point writes only itself and its sub-points.  Serial, this loop was 9 of the 17.6 s of a natural_neon
  // calculation on 24 threads.
  {
    std::vector<std::vector<Unit *>> batches;
    for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
      batches.emplace_back();
      for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++)
        batches.back().push_back(Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i));
    }
    if (FillChannelRows(batches, compound, configure, [&configure](const Unit &, const GSLException &e) {
          configure.outStream << "Component segment calculation error: " << e.what() << std::endl;
        }) == -1)
      return -1;
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
      point->CalcLegendreP(this, configure.maxLOrder, compound, effect);
    }
  }

  // CalcCoulombAmplitude for the leaves
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    for (int i = 1; i <= (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size(); i++) {
      Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);
      point->CalcCoulombAmplitude(this, compound);
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
      in.open(configure.OpenPath(configure.integralsfile).c_str());
      // Resume reading right after the regular-segment integrals recorded earlier.
      if (in.is_open() && ecReadPos_ > std::streampos(0)) in.seekg(ecReadPos_);
    } else {
      std::string outputfile;
      if (configure.paramMask & Config::CALCULATE_WITH_DATA)
        outputfile = configure.outputdir + "intEC.dat";
      else
        outputfile = configure.outputdir + "intEC.extrap";
      // Append the component-segment integrals after the regular-segment ones.
      out.open(configure.OpenPath(outputfile).c_str(), std::ios::app);
      if (!out) configure.outStream << "Could not append to EC Amplitude File." << std::endl;
      out.precision(17);  // as the data-segment amplitudes (CalculateECAmplitudes)
      if (!(configure.paramMask & Config::USE_API))
        configure.outStream << "Calculating EC Amplitudes for Component Segments..." << std::endl;
    }
    if (this->LeafECAmplitudes(compound, configure, usePrevious, in, out) == -1) return -1;
    if (out.is_open()) {
      out.flush();
      out.close();
      // The file is complete: record what it was computed for.  (Only when
      // CalculateECAmplitudes started it, which is what set the signature.)
      std::ofstream sig;
      if (!ecOutputFile_.empty() && !ecSignature_.empty()) sig.open(configure.OpenPath(ECSignaturePath(ecOutputFile_)).c_str());
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
 * compound nucleus).  The component's entrance pair is an entrance of the run
 * (RequestDataEntrances, RequestExtrapolationEntrances: every term's entrance), so
 * CompoundSystem::Initialize builds its pathways.
 * Returns false, with the error printed, when the component cannot be made.
 */
bool WorkingCopy::CheckComponentPairs(const CompoundSystem *compound, const Segment &baseSegment, int entranceKey, int exitKey,
                                const Config &configure) {
  if (!compound->IsPairKey(entranceKey) || !compound->IsPairKey(exitKey)) {
    configure.outStream << "ERROR: segment #" << baseSegment.GetSegmentKey() << ": component " << entranceKey
                        << " -> " << exitKey << " names a pair key that is not in the compound nucleus." << std::endl;
    return false;
  }
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(entranceKey));
  const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(exitKey));
  if (entrancePair->GetNucleus() != exitPair->GetNucleus()) {
    configure.outStream << "ERROR: segment #" << baseSegment.GetSegmentKey() << ": component " << entranceKey
                        << " -> " << exitKey << " connects two compound nuclei ("
                        << compound->NucleusLabel(entrancePair->GetNucleus()) << " and "
                        << compound->NucleusLabel(exitPair->GetNucleus()) << ")." << std::endl;
    return false;
  }
  return true;
}

/*!
 * A transmission or yield segment (Segment::IsTransformed), once its components exist: refuse what the transform
 * cannot mean, make every component a plain cross section, and give a yield its hidden total-cross-section term
 * (entrance pair to itself, "total 1", by the optical theorem), the s_t of Y0 = (s / s_t)(1 - exp(-n s_t)).
 */
bool WorkingCopy::FinishTransformedSegment(const CompoundSystem *compound, Segment *seg, int key, bool isTest, const Config &configure) {
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
  const PPair *entrance = compound->GetPair(compound->GetPairNumFromKey(seg->GetEntranceKey()));
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
 * one entrance pair do (their grids are laid out alike); the terms of a sample with several compound nuclei do not.
 * Refused here rather than computed wrongly.
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

/*!
 * The calculation plan (CalculationPlan.h): one task per point, built and checked.
 */
void WorkingCopy::BuildPlan() { plan_.Build(*this); }

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
      f(point, point.ResultAngularCapacity(capacity_), nat);
      for (int i = 1; i <= point.NumSubPoints(); i++)
        if (WritesSlot(point, i))
          f(*point.GetSubPoint(*this, i), point.GetSubPoint(*this, i)->ResultAngularCapacity(capacity_), 0);
    }
  };
  for (const Segment &seg : set_->segments_) segment(seg, 0);
  for (const std::pair<int, int> &leaf : set_->leafRefs_) segment(set_->segments_[leaf.first], leaf.second);
}

void WorkingCopy::LayoutResults(azure2::Results *results) {
  size_t total = 0;
  ForEachResultUnit([&](Site &, int nAng, int nNat) { total += azure2::Results::RecordSize(nAng, nNat); });
  std::vector<double> next(results ? total : 0, 0.0);
  size_t offset = 0;
  ForEachResultUnit([&](Site &point, int nAng, int nNat) {
    if (results) {
      const size_t from = point.ResultSlot();
      const int fromAng = point.ResultAngularSlots(), fromNat = point.ResultNaturalSlots();
      if (from != Site::kNoRecord && from + azure2::Results::RecordSize(fromAng, fromNat) <= results->slot.size())
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
      point.SetModelRecord(modelSize_, point.ModelAngularCapacity());
      modelSize_ += azure2::Results::ModelRecordSize(point.ModelAngularCapacity());
    }
}

void WorkingCopy::BeginEvaluation(azure2::Results &results) {
  size_t offset = 0;
  bool laidOut = true;
  ForEachResultUnit([&](Site &point, int nAng, int nNat) {
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
  // purpose).  A copy has no grid levels of its own (GridLevels).
  WorkingCopy *dataCopy = new WorkingCopy(*this);
  // The copy's points (not their sub-points) start without the external-capture amplitudes of the set-up: a working
  // copy evaluates them at its points' energies (Unit::GetECAmplitudeWithShift).
  for (azure2::PointRef u : dataCopy->PlanUnits()) dataCopy->units_[u].ClearECAmplitudes();
  return dataCopy;
}
