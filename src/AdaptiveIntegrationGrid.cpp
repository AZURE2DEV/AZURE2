#include "AdaptiveIntegrationGrid.h"
#include "JGroup.h"
#include "ALevel.h"
#include "PPair.h"
#include <cmath>
#include <algorithm>
#include <iostream>

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

  // PASS 1: Add fine symmetric grids around each resonance FIRST
  // Create a list of resonance regions to exclude from coarse grid
  std::vector<std::pair<double, double>> resonanceRegions;

  for (const ResonanceInfo& res : resonances) {
    double resonanceEnergy = res.energy;
    double resonanceWidth = res.totalWidth;

    // Integration region: ±(multiplier × Γ) around resonance
    double regionHalfWidth = config_.resonanceWidthMultiplier * resonanceWidth / 2;

    // Calculate ideal symmetric grid limits
    double idealLowerLimit = resonanceEnergy - regionHalfWidth;
    double idealUpperLimit = resonanceEnergy + regionHalfWidth;

    // Fine step size based on desired points per width
    double fineStep = (idealUpperLimit - idealLowerLimit) / (config_.pointsPerWidth);

    // Clip to actual integration range
    double actualLowerLimit = std::max(idealLowerLimit, endEnergy);
    double actualUpperLimit = std::min(idealUpperLimit, startEnergy);

    // Store this region to exclude from coarse grid
    if (actualUpperLimit >= actualLowerLimit) {
      resonanceRegions.push_back(std::make_pair(actualLowerLimit, actualUpperLimit));
    }

    // Check if resonance is within range
    if (resonanceEnergy < endEnergy || resonanceEnergy > startEnergy) {
      // Resonance center is outside range - still add partial coverage
      if (actualUpperLimit >= actualLowerLimit) {
        // Some part of resonance region overlaps with integration range
        double currentEnergy = actualLowerLimit;
        while (currentEnergy <= actualUpperLimit) {
          grid.push_back(currentEnergy);
          currentEnergy += fineStep;
        }
      }
    } else {
      // Resonance center is within range - add symmetric grid

      // Add point at resonance center
      grid.push_back(resonanceEnergy);

      // Add symmetric points on both sides
      int maxSteps = static_cast<int>(regionHalfWidth / fineStep);

      for (int i = 1; i <= maxSteps; i++) {
        double offset = i * fineStep;

        // Point above resonance
        double upperPoint = resonanceEnergy + offset;
        if (upperPoint <= actualUpperLimit) {
          grid.push_back(upperPoint);
        }

        // Point below resonance (symmetric)
        double lowerPoint = resonanceEnergy - offset;
        if (lowerPoint >= actualLowerLimit) {
          grid.push_back(lowerPoint);
        }

        // Stop if we've reached the boundaries on both sides
        if (upperPoint > actualUpperLimit && lowerPoint < actualLowerLimit) {
          break;
        }
      }
    }
  }

  // PASS 2: Add coarse grid points ONLY in smooth regions (outside resonance regions)
  // This respects maxPoints limit
  double coarseStep = config_.baseEnergyStep;

  // Helper function to check if energy is in any resonance region
  auto isInResonanceRegion = [&resonanceRegions](double energy) -> bool {
    for (const auto& region : resonanceRegions) {
      if (energy >= region.first && energy <= region.second) {
        return true;
      }
    }
    return false;
  };

  // Add coarse points, skipping those in resonance regions
  double currentEnergy = startEnergy;
  while (currentEnergy >= endEnergy) {
    if (!isInResonanceRegion(currentEnergy)) {
      grid.push_back(currentEnergy);
    }
    currentEnergy -= coarseStep;
  }

  // Ensure endpoints
  grid.push_back(startEnergy);
  grid.push_back(endEnergy);

  // Sort grid in descending order (high to low energy)
  std::sort(grid.begin(), grid.end(), std::greater<double>());

  // Remove duplicates (keep points that are sufficiently different)
  std::vector<double> uniqueGrid;
  double tolerance = 1.0e-6; // 1 eV tolerance

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

      // Convert to CM energy: E_cm = E_excitation - S - E_ex
      double levelCMEnergy = levelExcitationEnergy - separationEnergy - excitationEnergy;

      // Calculate total width by summing all partial widths
      double totalWidth = 0.0;
      for (int ch = 1; ch <= numChannels; ch++) {
        double gamma = level->GetGamma(ch) / 1e6; // Convert eV to MeV
        totalWidth = sqrt(totalWidth * totalWidth + gamma * gamma);
      }

      // If width is zero or very small, use a minimum value
      if (totalWidth < 1.0e-6) {
        totalWidth = 0.001; // 1 keV default for very narrow resonances
      }

      // Check if this resonance falls within our integration range
      // Account for margin based on the resonance width
      double margin = totalWidth * config_.resonanceWidthMultiplier;
      if (levelCMEnergy >= endEnergy - margin && levelCMEnergy <= startEnergy + margin) {
        ResonanceInfo resInfo;
        resInfo.energy = levelCMEnergy;
        resInfo.totalWidth = totalWidth;
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
    double resonantRegionWidth = resonance.totalWidth * config_.resonanceWidthMultiplier;
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

  // Determine the resonant region extent for this particular resonance
  double resonantRegionWidth = nearestRes->totalWidth * config_.resonanceWidthMultiplier;

  double stepSize;

  // Debug output
  //std::cout << "DEBUG: Energy = " << energy << " MeV,"
  //          << " Nearest Resonance = " << nearestRes->energy << " MeV,"
  //          << " Distance = " << distance << " MeV,"
  //          << " Resonant Region Width = " << resonantRegionWidth << " MeV" << std::endl;


  if (distance <= resonantRegionWidth) {
    // Inside resonant region - step size based on resonance width
    // Fine step = Γ / pointsPerWidth
    double fineStep = nearestRes->totalWidth / config_.pointsPerWidth;

    // Smooth transition from fine step at resonance center to base step at edge
    // blendFactor = 0 at resonance center, = 1 at edge of resonant region
    double blendFactor = distance / resonantRegionWidth;

    // Cubic interpolation for smoother transition
    double smoothBlend = blendFactor * blendFactor * (3.0 - 2.0 * blendFactor);

    stepSize = fineStep + (config_.baseEnergyStep - fineStep) * smoothBlend;

    // For very narrow resonances, ensure we don't go below a reasonable limit
    // (avoid numerical issues with extremely small steps)
    if (stepSize < 1.0e-5) {
      stepSize = 1.0e-5; // 0.01 keV absolute minimum
    }

  } else {
    // Outside resonant regions - use base step
    stepSize = config_.baseEnergyStep;
  }

  return stepSize;
}
