#include "AdaptiveIntegrationGrid.h"
#include "JGroup.h"
#include "ALevel.h"
#include "PPair.h"
#include "CoulFunc.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {
bool DebugGridEnabled() {
  static const bool enabled = (std::getenv("AZR_DEBUG_GRID") != nullptr);
  return enabled;
}
}  // namespace

/*!
 * \brief Constructor
 */
AdaptiveIntegrationGrid::AdaptiveIntegrationGrid(const GridConfig &config) :
  config_(config) {
}

/*!
 * \brief Generate adaptive grid for target integration
 *
 * Algorithm:
 * 1. Add fine symmetric grids around each resonance FIRST
 *    - Track resonance regions (±3Γ around each resonance)
 *    - Add symmetric points with fine step = Γ/pointsPerWidth
 *    - Handle partial coverage when resonances near boundaries
 * 2. Add coarse grid ONLY in smooth regions (outside resonance regions)
 *    - Step size = baseEnergyStep (default 50 keV)
 *    - Skip any points falling within resonance regions
 * 3. Sort and remove duplicates
 *
 * Note: No overlap between fine and coarse grids - resonances get exclusive coverage.
 */
std::vector<double> AdaptiveIntegrationGrid::GenerateGrid(double startEnergy, double endEnergy, CNuc *compound) {
  std::vector<double> grid;

  // Ensure proper ordering (startEnergy > endEnergy for backward integration)
  if (startEnergy < endEnergy) {
    std::swap(startEnergy, endEnergy);
  }

  // Handle trivial case
  double totalRange = startEnergy - endEnergy;
  if (totalRange <= 0.0 || totalRange < 1.0e-6) {
    grid.push_back(startEnergy);
    return grid;
  }

  return GenerateGrid(startEnergy, endEnergy, IdentifyResonances(startEnergy, endEnergy, compound));
}

std::vector<double> AdaptiveIntegrationGrid::GenerateGrid(double startEnergy, double endEnergy,
                                                          const std::vector<ResonanceInfo> &resonances) {
  std::vector<double> grid;
  if (startEnergy < endEnergy) std::swap(startEnergy, endEnergy);
  double totalRange = startEnergy - endEnergy;
  if (totalRange <= 0.0 || totalRange < 1.0e-6) {
    grid.push_back(startEnergy);
    return grid;
  }

  static const bool debugGrid = (std::getenv("AZR_DEBUG_GRID") != nullptr);
  if (debugGrid) {
    fprintf(stderr, "[AZR_DEBUG_GRID] window [%.10f, %.10f] MeV, %zu resonances found:\n",
            endEnergy, startEnergy, resonances.size());
    for (const ResonanceInfo &r : resonances) {
      fprintf(stderr, "[AZR_DEBUG_GRID]   res E=%.10f MeV particleWidth=%.6e MeV totalWidth=%.6e MeV\n",
              r.energy, r.particleWidth, r.totalWidth);
    }
  }

  // Use smooth adaptive stepping - step size varies continuously based on distance to resonances
  // This avoids sharp boundaries that cause integration artifacts.
  //
  // Near a resonance, though, points are snapped onto a FIXED lattice anchored
  // at that resonance's own energy (spacing = particleWidth / pointsPerWidth)
  // rather than accumulated by marching from the window's own startEnergy.
  // Two overlapping integration windows (e.g. the same resonance seen by
  // neighboring data points, whose windows both extend back past it) walk the
  // fine region from different starting points; free marching gives each one
  // an essentially arbitrary phase relative to the resonance, so a subinterval
  // can straddle the peak very differently window to window even when both
  // are individually well converged. For a resonance many orders of magnitude
  // narrower than the window (the extreme, but not the only, case this
  // matters for), that phase drift makes the quadrature result swing sharply
  // between adjacent evaluation energies instead of varying smoothly. Anchoring
  // the fine lattice to the resonance itself removes that drift: every window
  // that reaches the same resonance samples it at the same fixed energies.

  //
  // Outside the lattice the tail still matters for a narrow resonance: beyond
  // resonanceWidthMultiplier widths a Lorentzian holds about 1/(pi * multiplier)
  // of its area (1.6 % at the default 20).  Stepping straight from the lattice
  // edge to the coarse step put that tail on a single trapezoid whose upper node
  // was fixed by the window, not by the resonance (DeBoer 2026-10-22, 17O+a: a
  // 0.8 eV level in a 10 keV target window, 50 sub-points: the 16 meV lattice
  // met 200 eV steps 16 eV from the peak).  The integral was up to 8 % of the
  // resonance's area off, and since the grids are laid out again at every
  // evaluation, the error changed as a fit moved the level: chi2 oscillated in
  // the level energy with the period of the coarse step and the analytic
  // gradient disagreed with finite differences.  So the anchored set goes on
  // past the lattice with GRADED nodes, each step tailStepFraction of the
  // distance to the resonance, until the step reaches baseEnergyStep.  The whole
  // zone (lattice plus graded nodes) moves with the resonance, and the coarse
  // grid meets it where the tail is resolved by the coarse step anyway.
  std::vector<Zone> zones = AnchoredZones(resonances);

  grid.push_back(startEnergy);
  double currentEnergy = startEnergy;

  while (currentEnergy > endEnergy) {
    // Inside the zone of one or more anchored resonances the next point is the
    // closest node of their sets below the current one (the union of the sets).
    // Only a resonance whose lattice is FINER than the smooth step has a zone. A
    // broad pole's pitch can exceed the whole window -- a 10 MeV background level
    // at pointsPerWidth 50 gives 200 keV, against a 100 keV thick target -- and
    // anchoring to it stepped straight from the window's start to its end, past
    // everything in between, leaving the grid [start, E_R, end]. Such a resonance
    // is left to CalculateAdaptiveStep, which is capped at baseEnergyStep.
    double nextEnergy = -std::numeric_limits<double>::infinity();
    const Zone *anchorZone = nullptr;
    for (const Zone &z : zones) {
      if (std::abs(currentEnergy - z.energy) > z.half + 1.0e-12) continue;  // the edge itself counts (see below)
      double e = NextNodeBelow(z, currentEnergy);
      if (e > nextEnergy) {
        nextEnergy = e;
        anchorZone = &z;
      }
    }
    bool anchored = (anchorZone != nullptr && std::isfinite(nextEnergy));
    if (!anchored) {
      anchorZone = nullptr;
      double stepSize = CalculateAdaptiveStep(currentEnergy, resonances);
      nextEnergy = currentEnergy - stepSize;
      // Safety check to avoid infinite loops in the smooth/background regime;
      // the anchored branch above always advances to a node strictly below,
      // so it cannot stall and does not need this guard.
      if (stepSize < 1.0e-10) {
        if (nextEnergy < endEnergy) nextEnergy = endEnergy;
        grid.push_back(nextEnergy);
        break;
      }
    }

    // Never step INTO a zone from outside it: stop at its upper edge, so the
    // zone is sampled on its own nodes from there down. The smooth step is up to
    // baseEnergyStep long, which can land deep inside a narrow region (or clear
    // across it) -- a 4.5 keV step from 264.9 keV to 260.4 keV skipped the whole
    // upper wing of the 0.97 keV-wide 259.7 keV resonance in 14N(p,g), and the
    // trapezoid then spread sigma at 0.7 keV from the peak over the full step.
    for (const Zone &z : zones) {
      double upperEdge = z.energy + z.half;
      if (upperEdge < currentEnergy - 1.0e-12 && upperEdge > nextEnergy) nextEnergy = upperEdge;
    }

    if (nextEnergy < endEnergy) {
      nextEnergy = endEnergy;
    }

    if (debugGrid) {
      fprintf(stderr, "[AZR_DEBUG_GRID]   step: current=%.10f -> next=%.10f anchored=%d anchorE=%.10f anchorPitch=%.6e\n",
              currentEnergy, nextEnergy, anchored ? 1 : 0,
              anchorZone ? anchorZone->energy : -1.0, anchorZone ? anchorZone->pitch : -1.0);
    }

    grid.push_back(nextEnergy);
    currentEnergy = nextEnergy;
  }

  // Ensure endpoint is included
  if (std::abs(grid.back() - endEnergy) > 1.0e-10) {
    grid.push_back(endEnergy);
  }

  // Also ensure resonance centers are included as grid points for accuracy
  for (const ResonanceInfo &res : resonances) {
    if (res.energy >= endEnergy && res.energy <= startEnergy) {
      grid.push_back(res.energy);
    }
  }

  // Sort grid in descending order (high to low energy)
  std::sort(grid.begin(), grid.end(), std::greater<double>());

  // Remove duplicates (keep points that are sufficiently different). The
  // default 0.01 eV tolerance is fine for ordinary (>~keV-wide) resonances,
  // but a resonance narrow enough to need a sub-eV lattice pitch would have
  // most of its own fine points collapsed back together by a fixed 0.01 eV
  // tolerance, silently undoing the anchoring above. Scale the tolerance down
  // to track the finest lattice pitch actually in use, floored well below any
  // physically meaningful width so exact duplicates still collapse.
  double tolerance = 1.0e-8;  // 0.01 eV, unchanged default for the no-resonance case
  for (const ResonanceInfo &res : resonances) {
    if (res.particleWidth > 0.0 && config_.pointsPerWidth > 0.0) {
      double pitch = res.particleWidth / config_.pointsPerWidth;
      double candidateTol = pitch * 1.0e-3;
      if (candidateTol > 0.0 && candidateTol < tolerance) tolerance = candidateTol;
    }
  }
  if (tolerance < 1.0e-14) tolerance = 1.0e-14;

  std::vector<double> uniqueGrid;
  if (!grid.empty()) {
    uniqueGrid.push_back(grid[0]);
    for (size_t i = 1; i < grid.size(); i++) {
      if (std::abs(grid[i] - uniqueGrid.back()) > tolerance) {
        uniqueGrid.push_back(grid[i]);
      }
    }
  }

  return uniqueGrid;
}

namespace {
// Is a resonance of this particle width finer than the coarse step, so that it gets a
// lattice of its own (see GenerateGrid)?
bool HasLattice(double particleWidth, const AdaptiveIntegrationGrid::GridConfig &c) {
  return particleWidth > 0.0 && c.pointsPerWidth > 0.0 && particleWidth / c.pointsPerWidth < c.baseEnergyStep;
}
}  // namespace

std::vector<double> AdaptiveIntegrationGrid::GradedOffsets(double pitch, double latticeTop) const {
  std::vector<double> graded;
  const double f = config_.tailStepFraction;
  if (!(f > 0.0) || !(pitch > 0.0)) return graded;
  double x = latticeTop;
  for (;;) {
    const double step = std::max(pitch, f * x);
    if (step >= config_.baseEnergyStep) break;
    x += step;
    graded.push_back(x);
  }
  return graded;
}

double AdaptiveIntegrationGrid::ZoneHalfWidth(double particleWidth) const {
  const double margin = particleWidth * config_.resonanceWidthMultiplier;
  if (!HasLattice(particleWidth, config_)) return margin;
  const double pitch = particleWidth / config_.pointsPerWidth;
  const double latticeTop = std::ceil(margin / pitch - 1.0e-9) * pitch;
  std::vector<double> graded = GradedOffsets(pitch, latticeTop);
  return graded.empty() ? latticeTop : graded.back();
}

std::vector<AdaptiveIntegrationGrid::Zone>
AdaptiveIntegrationGrid::AnchoredZones(const std::vector<ResonanceInfo> &resonances) const {
  std::vector<Zone> zones;
  for (const ResonanceInfo &res : resonances) {
    if (!HasLattice(res.particleWidth, config_)) continue;
    Zone z;
    z.energy = res.energy;
    z.pitch = res.particleWidth / config_.pointsPerWidth;
    z.numLattice = (int)std::ceil(res.particleWidth * config_.resonanceWidthMultiplier / z.pitch - 1.0e-9);
    z.graded = GradedOffsets(z.pitch, z.numLattice * z.pitch);
    z.half = z.graded.empty() ? z.numLattice * z.pitch : z.graded.back();
    zones.push_back(std::move(z));
  }
  return zones;
}

double AdaptiveIntegrationGrid::NextNodeBelow(const Zone &z, double energy) {
  // Node energies are compared directly: (energy - E_R) / pitch is not a reliable
  // index, since energies of order 10 MeV carry rounding errors of 1e-7 of a
  // 16 neV pitch.
  const double eps = 1.0e-12;
  const double latticeTop = z.numLattice * z.pitch;
  const double t = energy - z.energy;
  if (t > latticeTop + eps) {
    // above the lattice: the largest graded offset below t, else the lattice top
    auto it = std::lower_bound(z.graded.begin(), z.graded.end(), t - eps);
    if (it != z.graded.begin()) return z.energy + *(it - 1);
    return z.energy + latticeTop;
  }
  if (t >= -latticeTop - eps) {
    double n = std::floor(t / z.pitch + 0.5);
    while (z.energy + n * z.pitch >= energy - eps) n -= 1.0;
    if (n >= -z.numLattice) return z.energy + n * z.pitch;
  }
  // below the lattice: the next graded offset down
  auto it = std::upper_bound(z.graded.begin(), z.graded.end(), std::max(-t, latticeTop) + eps);
  if (it != z.graded.end()) return z.energy - *it;
  return -std::numeric_limits<double>::infinity();
}

/*!
 * \brief Generate grid with detailed information about each point
 */
std::vector<AdaptiveIntegrationGrid::GridPoint>
AdaptiveIntegrationGrid::GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc *compound) {
  // First generate the regular grid
  std::vector<double> energyGrid = GenerateGrid(startEnergy, endEnergy, compound);

  // Then annotate with resonance information
  std::vector<ResonanceInfo> resonances = IdentifyResonances(startEnergy, endEnergy, compound);

  std::vector<GridPoint> gridInfo;
  for (double energy : energyGrid) {
    GridPoint point;
    point.energy = energy;
    point.isResonant = IsInResonantRegion(energy, resonances);
    gridInfo.push_back(point);
  }

  return gridInfo;
}

/*!
 * \brief Get expected point count
 */
int AdaptiveIntegrationGrid::GetExpectedPointCount(double startEnergy, double endEnergy, CNuc *compound) {
  std::vector<double> grid = GenerateGrid(startEnergy, endEnergy, compound);
  return grid.size();
}

/*!
 * \brief Set grid configuration
 */
void AdaptiveIntegrationGrid::SetConfig(const GridConfig &config) {
  config_ = config;
}

/*!
 * \brief Get current configuration
 */
AdaptiveIntegrationGrid::GridConfig AdaptiveIntegrationGrid::GetConfig() const {
  return config_;
}

/*!
 * \brief Identify resonances within the energy range with their total widths
 *
 * Loops through all J-groups and levels in the compound nucleus to find
 * resonances that fall within the integration energy range. For each resonance,
 * calculates the total width by summing all partial widths (Γ_total = Σ Γ_i).
 * Converts level energies from compound excitation energy to CM energy using
 * the separation energy of the entrance channel.
 */
std::vector<AdaptiveIntegrationGrid::ResonanceInfo>
AdaptiveIntegrationGrid::IdentifyResonances(double startEnergy, double endEnergy, CNuc *compound) {
  std::vector<ResonanceInfo> resonances;

  if (!compound) return resonances;

  // Get the entrance pair for separation energy conversion
  PPair *entrancePair = nullptr;
  if (config_.entranceKey > 0 && compound->IsPairKey(config_.entranceKey)) {
    int pairNum = compound->GetPairNumFromKey(config_.entranceKey);
    entrancePair = compound->GetPair(pairNum);
  }

  // If we don't have entrance pair info, we can't convert energies accurately
  if (!entrancePair) {
    return resonances;
  }

  double separationEnergy = entrancePair->GetSepE();
  double excitationEnergy = entrancePair->GetExE();

  // Loop through all J-groups
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jgroup = compound->GetJGroup(j);
    if (!jgroup || !jgroup->IsInRMatrix()) continue;
    // Levels of another compound nucleus are not resonances of this entrance pair.
    if (jgroup->GetNucleus() != entrancePair->GetNucleus()) continue;

    int numChannels = jgroup->NumChannels();

    // Loop through all levels in this J-group
    for (int l = 1; l <= jgroup->NumLevels(); l++) {
      ALevel *level = jgroup->GetLevel(l);
      if (!level || !level->IsInRMatrix()) continue;

      // Get level energy (in compound excitation energy)
      double levelExcitationEnergy = level->GetE();

      // Check if level has widths (Γ > 0) to be considered a resonance
      bool hasWidth = false;
      for (int ch = 1; ch <= numChannels; ch++) {
        if (std::abs(level->GetGamma(ch)) > 1.0e-6) {
          hasWidth = true;
          break;
        }
      }
      if (!hasWidth) continue;  // skip levels with no widths

      // Convert to CM energy: E_cm = E_excitation - S - E_ex
      double levelCMEnergy = levelExcitationEnergy - separationEnergy - excitationEnergy - config_.resonanceOffset;

      // Total width Gamma_total = Gamma_particle + Gamma_radiative, in MeV.
      // The Breit-Wigner cross section carries the total width in its
      // denominator -- every open decay channel broadens the resonance -- so
      // this is the width that shapes sigma(E) and that the grid must resolve.
      // For radiative capture Gamma_gamma can dominate (Gamma_gamma >>
      // Gamma_p); leaving it out makes the resonance look far narrower than it
      // is in the extrapolated cross section.
      //
      // GetGamma() returns the reduced width amplitude gamma (MeV^1/2) after
      // CNuc::TransformIn, so it is used directly -- NOT divided by 1e6.  That
      // 1e6 belongs to the pre-transformation convention (input widths in eV);
      // applying it here made every width 1e12 too small.  The particle sum is
      // level-shift normalised by (1 + sum gamma^2 dS/dE) so the result is the
      // *observed* width, matching CNuc::TransformOut / parameters.out.
      double totalWidth = 0.0;
      double normSum = 0.0;
      // Before the input transformation has run (which is the case whenever
      // the data structures -- and so these grids -- are filled), GetGamma()
      // still holds the values as read: the observed partial width in eV for
      // an open particle channel, an ANC for a closed one, Gamma_gamma in eV
      // for a photon channel.  Those are the widths directly; running them
      // through 2 gamma^2 P as if they were amplitudes gave every resonance a
      // width of order 2P/(dS/dE), independent of the level (a 0.6 keV 2+ in
      // 12C+alpha came out 0.8 MeV wide and was never resolved).
      if (config_.inputWidthsArePhysical) {
        double particleWidthEV = 0.0;
        double gammaWidthEV = 0.0;
        for (int ch = 1; ch <= numChannels; ch++) {
          AChannel *channel = jgroup->GetChannel(ch);
          PPair *chPair = compound->GetPair(channel->GetPairNum());
          double gamma = std::abs(level->GetGamma(ch));
          if (gamma <= 0.0) continue;
          if (channel->GetRadType() == 'P') {
            double localEnergy = level->GetE() - chPair->GetExE() - chPair->GetSepE();
            if (localEnergy <= 0.0) continue;  // closed channel: the value is an ANC
            particleWidthEV += gamma;
          } else if (channel->GetRadType() == 'M' || channel->GetRadType() == 'E') {
            gammaWidthEV += gamma;
          }
        }
        double particleWidth = particleWidthEV * 1.0e-6;
        double totalWidthMeV = (particleWidthEV + gammaWidthEV) * 1.0e-6;
        if (DebugGridEnabled()) {
          fprintf(stderr, "[AZR_DEBUG_GRID]     level E=%.6f: input widths particle=%.6e eV gamma=%.6e eV\n",
                  level->GetE(), particleWidthEV, gammaWidthEV);
        }
        double margin = ZoneHalfWidth(particleWidth);
        if (levelCMEnergy >= endEnergy - margin && levelCMEnergy <= startEnergy + margin) {
          ResonanceInfo resInfo;
          resInfo.energy = levelCMEnergy;
          resInfo.totalWidth = totalWidthMeV;
          resInfo.particleWidth = particleWidth;
          resonances.push_back(resInfo);
        }
        continue;
      }
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *channel = jgroup->GetChannel(ch);
        if (channel->GetRadType() != 'P') continue;  // particle channels
        PPair *chPair = compound->GetPair(channel->GetPairNum());
        double localEnergy = level->GetE() - chPair->GetExE() - chPair->GetSepE();
        if (localEnergy <= 0.0) continue;  // sub-threshold channel
        double gamma = std::abs(level->GetGamma(ch));
        if (gamma <= 0.0) continue;
        double radius = chPair->GetChRad();
        CoulFunc coulFunc(chPair, false);
        double pene = coulFunc.Penetrability(channel->GetL(), radius, localEnergy);
        double dSdE = coulFunc.PEShift_dE(channel->GetL(), radius, localEnergy);
        totalWidth += 2.0 * gamma * gamma * pene;
        normSum += dSdE * gamma * gamma;
        if (DebugGridEnabled()) {
          fprintf(stderr, "[AZR_DEBUG_GRID]     ch=%d l=%d Elocal=%.6e radius=%.4f gamma=%.6e "
                  "pene=%.6e dSdE=%.6e chanW=%.6e(eV)\n",
                  ch, channel->GetL(), localEnergy, radius, gamma, pene, dSdE,
                  2.0 * gamma * gamma * pene * 1.0e6);
        }
      }
      if (1.0 + normSum > 0.0) totalWidth /= (1.0 + normSum);

      double particleWidth = totalWidth;  // before radiative channels

      // Radiative channels (M/E): Gamma = 2 gamma^2 P_rad, with the radiation
      // penetrability following CNuc::TransformOut (the same special cases for a
      // ground-state transition and the RMC formalism).
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *channel = jgroup->GetChannel(ch);
        char radType = channel->GetRadType();
        if (radType != 'M' && radType != 'E') continue;
        PPair *chPair = compound->GetPair(channel->GetPairNum());
        double gamma = std::abs(level->GetGamma(ch));
        if (gamma <= 0.0) continue;
        double localEnergy = level->GetE() - chPair->GetExE() - chPair->GetSepE();
        double pene;
        if (std::abs(level->GetE() - chPair->GetExE()) < 1.0e-3 &&
            jgroup->GetJ() == chPair->GetJ(2) &&
            jgroup->GetPi() == chPair->GetPi(2)) {
          double jValue = jgroup->GetJ();
          pene = 1.0e-10;
          if (radType == 'M' && channel->GetL() == 1)
            pene = 3.0 * jValue / 4.0 / (jValue + 1.) / nuclearMagneton / nuclearMagneton;
          else if (radType == 'E' && channel->GetL() == 2)
            pene = 60.0 * jValue * (2. * jValue - 1.) / (jValue + 1.) / (2. * jValue + 3.);
        } else {
          pene = pow(std::abs(localEnergy) / hbarc, 2.0 * channel->GetL() + 1.0);
        }
        totalWidth += 2.0 * gamma * gamma * std::abs(pene);
      }

      // The resonance in sigma(E) is shaped by the particle width, so the grid
      // is sized by it (see CalculateAdaptiveStep / IsInResonantRegion).
      // A resonance outside the window whose zone (lattice and graded tail) reaches
      // into it counts: its nodes there move with it, as they do once it is inside.
      double margin = ZoneHalfWidth(particleWidth);
      if (levelCMEnergy >= endEnergy - margin && levelCMEnergy <= startEnergy + margin) {
        ResonanceInfo resInfo;
        resInfo.energy = levelCMEnergy;
        resInfo.totalWidth = totalWidth;
        resInfo.particleWidth = particleWidth;
        resonances.push_back(resInfo);
      }
    }
  }

  // Sort resonances by energy in ascending order for efficient searching
  std::sort(resonances.begin(), resonances.end(),
            [](const ResonanceInfo &a, const ResonanceInfo &b) {
              return a.energy < b.energy;
            });

  return resonances;
}

/*!
 * \brief Check if an energy is within a resonant region
 */
bool AdaptiveIntegrationGrid::IsInResonantRegion(double energy, const std::vector<ResonanceInfo> &resonances) const {
  for (const ResonanceInfo &resonance : resonances) {
    double resonantRegionWidth = resonance.particleWidth * config_.resonanceWidthMultiplier;
    if (std::abs(energy - resonance.energy) <= resonantRegionWidth) {
      return true;
    }
  }
  return false;
}

/*!
 * \brief Find the nearest resonance to a given energy
 */
const AdaptiveIntegrationGrid::ResonanceInfo *
AdaptiveIntegrationGrid::FindNearestResonance(double energy,
                                              const std::vector<ResonanceInfo> &resonances,
                                              double &distance) const {
  if (resonances.empty()) {
    distance = 1e10;
    return nullptr;
  }

  const ResonanceInfo *nearest = nullptr;
  double minDistance = 1e10;

  for (const ResonanceInfo &resonance : resonances) {
    double dist = std::abs(energy - resonance.energy);
    if (dist < minDistance) {
      minDistance = dist;
      nearest = &resonance;
    }
  }

  distance = minDistance;
  return nearest;
}

/*!
 * \brief Calculate adaptive step size based on proximity to resonances
 *
 * The step size is determined by the width of the nearest resonance:
 * - Near resonances: stepSize ~ Γ / pointsPerWidth
 * - Far from resonances: stepSize = baseEnergyStep
 * - Smooth transition between the two regimes
 *
 * This ensures that narrow resonances get fine grids and broad resonances
 * get appropriately coarser grids, all based on the actual physics.
 */
double AdaptiveIntegrationGrid::CalculateAdaptiveStep(double energy, const std::vector<ResonanceInfo> &resonances) const {
  // Find nearest resonance and distance to it
  double distance;
  const ResonanceInfo *nearestRes = FindNearestResonance(energy, resonances, distance);

  if (!nearestRes) {
    // No resonances - use base step
    return config_.baseEnergyStep;
  }

  // Fine step at resonance center (particle width shapes sigma(E)).
  double fineStep = nearestRes->particleWidth / config_.pointsPerWidth;

  double sigma = nearestRes->particleWidth * config_.resonanceWidthMultiplier / 2.0;

  // A zero width (e.g. a level whose channels are all sub-threshold) would make
  // the falloff below 0/0.  Such a level has no resonant structure to resolve,
  // so fall back to the base step.
  if (!(sigma > 0.0)) return config_.baseEnergyStep;

  // Gaussian falloff factor: 1 at center, approaches 0 far away
  double gaussianFactor = std::exp(-distance * distance / (2.0 * sigma * sigma));

  // Step size: fine at center, coarse far away
  double stepSize = fineStep * gaussianFactor + config_.baseEnergyStep * (1.0 - gaussianFactor);

  // For very narrow resonances, ensure we don't go below a reasonable limit
  if (stepSize < 1.0e-5) {
    stepSize = 1.0e-5;  // 0.01 keV absolute minimum
  }

  // Also cap at baseEnergyStep to avoid overly large steps
  if (stepSize > config_.baseEnergyStep) {
    stepSize = config_.baseEnergyStep;
  }

  return stepSize;
}
