// The set-up checks of a data set with a secondary gamma ray (UposChecks.h), called by WorkingCopy::Fill and MakePoints.
#include "UposChecks.h"

#include <iostream>
#include <string>

#include "CompoundSystem.h"
#include "Config.h"
#include "PPair.h"
#include "Segment.h"

/*!
 * A UPOS (secondary-gamma) segment given as a lab differential is converted with the kinematics of particle b, the
 * light ejectile of the exit pair, which is right for a gamma measured in coincidence with it and not otherwise
 * (segments.rst, UPOS).  Say so instead of converting silently (audit 2026-10-25, item 31).
 */
void UPOSLabNote(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !segment.IsDifferential() || segment.IsCMDifferential() || segment.IsPhase()) return;
  // A capture cascade is converted as a capture gamma ray (emitted by the recoiling compound system; the recoil from
  // the primary photon is neglected), which is what the secondary gamma ray is.
  if (compound->IsPairKey(segment.GetExitKey()) &&
      compound->GetPair(compound->GetPairNumFromKey(segment.GetExitKey()))->GetPType() == 10)
    return;
  configure.outStream << "WARNING: " << (isTest ? "test segment #" : "segment #") << key
                      << " (UPOS, lab differential): the lab angle is converted with the kinematics of particle b, "
                      << "the light ejectile of the exit pair (not the gamma's); give c.m. angles if that is not meant."
                      << std::endl;
}

/*!
 * The UPOS coefficient (CompoundSystem, Brune and deBoer 2020 Eq. 22) is written for a secondary gamma ray from the HEAVY
 * partner of the exit pair: its final-state spin is J2.  With the light particle excited (column 32 of the
 * <levels> line, 2026-11-05) the gamma ray would come from the light particle, which that formula does not
 * describe -- refuse instead of computing the heavy partner's correlation silently.
 */
bool UPOSAllowed(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure) {
  if (!segment.IsUPOS() || !compound->IsPairKey(segment.GetExitKey())) return true;
  const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(segment.GetExitKey()));
  if (exitPair->GetPType() == 10) {
    // A capture cascade (the secondary gamma ray of the capture final state, primary photon unobserved).
    std::string why;
    if (configure.paramMask & Config::USE_RMC_FORMALISM)
      why = "not with the R-matrix-for-capture formalism (--use-rmc)";
    else if (segment.GetLineWindow() > 0.0)
      why = "not together with a line window (the line-shape integration has no secondary-gamma form yet)";
    else if (segment.IsAnalyzingPower())
      why = "an analyzing power of the secondary gamma ray is not implemented";
    if (why.empty()) return true;
    configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << key
                        << " (UPOS on a capture pair, the secondary gamma ray of a capture cascade): " << why << "."
                        << std::endl;
    return false;
  }
  if (exitPair->GetLightExE() == 0.) return true;
  configure.outStream << "ERROR: " << (isTest ? "test segment #" : "segment #") << key << " (UPOS): exit pair "
                      << segment.GetExitKey() << " has an excited light particle (column 32, " << exitPair->GetLightExE()
                      << " MeV); the secondary-gamma correlation is implemented for the decay of the heavy partner only."
                      << std::endl;
  return false;
}
