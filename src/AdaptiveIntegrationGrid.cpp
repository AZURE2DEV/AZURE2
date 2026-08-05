#include "AdaptiveIntegrationGrid.h"
#include "JGroup.h"
#include "ALevel.h"
#include "PPair.h"
#include "CoulFunc.h"
#include <cmath>
#include <algorithm>

/*!
 * \brief Constructor
 */
AdaptiveIntegrationGrid::AdaptiveIntegrationGrid(const GridConfig& config)
  : config_(config) {
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
std::vector<double> AdaptiveIntegrationGrid::GenerateGrid(double startEnergy, double endEnergy, CNuc* compound) {
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

  // Identify resonances in the energy range with their widths
  std::vector<ResonanceInfo> resonances = IdentifyResonances(startEnergy, endEnergy, compound);

  // Use smooth adaptive stepping - step size varies continuously based on distance to resonances
  // This avoids sharp boundaries that cause integration artifacts

  grid.push_back(startEnergy);
  double currentEnergy = startEnergy;

  while (currentEnergy > endEnergy) {
    // Calculate adaptive step size at current energy
    double stepSize = CalculateAdaptiveStep(currentEnergy, resonances);

    // Ensure we don't overshoot the end
    double nextEnergy = currentEnergy - stepSize;
    if (nextEnergy < endEnergy) {
      nextEnergy = endEnergy;
    }

    // Add the next point
    grid.push_back(nextEnergy);
    currentEnergy = nextEnergy;

    // Safety check to avoid infinite loops
    if (stepSize < 1.0e-10) {
      break;
    }
  }

  // Ensure endpoint is included
  if (std::abs(grid.back() - endEnergy) > 1.0e-10) {
    grid.push_back(endEnergy);
  }

  // Also ensure resonance centers are included as grid points for accuracy
  for (const ResonanceInfo& res : resonances) {
    if (res.energy >= endEnergy && res.energy <= startEnergy) {
      grid.push_back(res.energy);
    }
  }

  // Sort grid in descending order (high to low energy)
  std::sort(grid.begin(), grid.end(), std::greater<double>());

  // Remove duplicates (keep points that are sufficiently different)
  std::vector<double> uniqueGrid;
  double tolerance = 1.0e-8; // 0.01 eV tolerance

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

/*!
 * \brief Generate grid with detailed information about each point
 */
std::vector<AdaptiveIntegrationGrid::GridPoint>
AdaptiveIntegrationGrid::GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc* compound) {
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
int AdaptiveIntegrationGrid::GetExpectedPointCount(double startEnergy, double endEnergy, CNuc* compound) {
  std::vector<double> grid = GenerateGrid(startEnergy, endEnergy, compound);
  return grid.size();
}

/*!
 * \brief Set grid configuration
 */
void AdaptiveIntegrationGrid::SetConfig(const GridConfig& config) {
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
AdaptiveIntegrationGrid::IdentifyResonances(double startEnergy, double endEnergy, CNuc* compound) {
  std::vector<ResonanceInfo> resonances;

  if (!compound) return resonances;

  // Get the entrance pair for separation energy conversion
  PPair* entrancePair = nullptr;
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
    JGroup* jgroup = compound->GetJGroup(j);
    if (!jgroup || !jgroup->IsInRMatrix()) continue;

    int numChannels = jgroup->NumChannels();

    // Loop through all levels in this J-group
    for (int l = 1; l <= jgroup->NumLevels(); l++) {
      ALevel* level = jgroup->GetLevel(l);
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
      if (!hasWidth) continue; // skip levels with no widths

      // Convert to CM energy: E_cm = E_excitation - S - E_ex
      double levelCMEnergy = levelExcitationEnergy - separationEnergy - excitationEnergy;

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
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel* channel = jgroup->GetChannel(ch);
        if (channel->GetRadType() != 'P') continue;   // particle channels
        PPair* chPair = compound->GetPair(channel->GetPairNum());
        double localEnergy = level->GetE() - chPair->GetExE() - chPair->GetSepE();
        if (localEnergy <= 0.0) continue;             // sub-threshold channel
        double gamma = std::abs(level->GetGamma(ch));
        if (gamma <= 0.0) continue;
        double radius = chPair->GetChRad();
        CoulFunc coulFunc(chPair, false);
        double pene = coulFunc.Penetrability(channel->GetL(), radius, localEnergy);
        totalWidth += 2.0 * gamma * gamma * pene;
        normSum += coulFunc.PEShift_dE(channel->GetL(), radius, localEnergy)
                   * gamma * gamma;
      }
      if (1.0 + normSum > 0.0) totalWidth /= (1.0 + normSum);

      double particleWidth = totalWidth;   // before radiative channels

      // Radiative channels (M/E): Gamma = 2 gamma^2 P_rad, with the radiation
      // penetrability following CNuc::TransformOut (the same special cases for a
      // ground-state transition and the RMC formalism).
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel* channel = jgroup->GetChannel(ch);
        char radType = channel->GetRadType();
        if (radType != 'M' && radType != 'E') continue;
        PPair* chPair = compound->GetPair(channel->GetPairNum());
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
            pene = 3.0*jValue/4.0/(jValue+1.)/nuclearMagneton/nuclearMagneton;
          else if (radType == 'E' && channel->GetL() == 2)
            pene = 60.0*jValue*(2.*jValue-1.)/(jValue+1.)/(2.*jValue+3.);
        } else {
          pene = pow(std::abs(localEnergy)/hbarc, 2.0*channel->GetL()+1.0);
        }
        totalWidth += 2.0 * gamma * gamma * std::abs(pene);
      }

      // The resonance in sigma(E) is shaped by the particle width, so the grid
      // is sized by it (see CalculateAdaptiveStep / IsInResonantRegion).
      double margin = particleWidth * config_.resonanceWidthMultiplier;
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
            [](const ResonanceInfo& a, const ResonanceInfo& b) {
              return a.energy < b.energy;
            });

  return resonances;
}

/*!
 * \brief Check if an energy is within a resonant region
 */
bool AdaptiveIntegrationGrid::IsInResonantRegion(double energy, const std::vector<ResonanceInfo>& resonances) const {
  for (const ResonanceInfo& resonance : resonances) {
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
const AdaptiveIntegrationGrid::ResonanceInfo*
AdaptiveIntegrationGrid::FindNearestResonance(double energy,
                                               const std::vector<ResonanceInfo>& resonances,
                                               double& distance) const {
  if (resonances.empty()) {
    distance = 1e10;
    return nullptr;
  }

  const ResonanceInfo* nearest = nullptr;
  double minDistance = 1e10;

  for (const ResonanceInfo& resonance : resonances) {
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
double AdaptiveIntegrationGrid::CalculateAdaptiveStep(double energy, const std::vector<ResonanceInfo>& resonances) const {
  // Find nearest resonance and distance to it
  double distance;
  const ResonanceInfo* nearestRes = FindNearestResonance(energy, resonances, distance);

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
    stepSize = 1.0e-5; // 0.01 keV absolute minimum
  }

  // Also cap at baseEnergyStep to avoid overly large steps
  if (stepSize > config_.baseEnergyStep) {
    stepSize = config_.baseEnergyStep;
  }

  return stepSize;
}
