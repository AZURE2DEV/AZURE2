// WorkingCopy, reading the data: the segments of the data sets (Fill) and of the extrapolations (MakePoints), the
// gates that refuse a segment a quantity it cannot carry, and the check of the <thm> weight tables against the
// points read.  The set-up of the points read is WorkingCopy.cpp.
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
