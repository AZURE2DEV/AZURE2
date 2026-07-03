#ifndef ADAPTIVEINTEGRATIONGRID_H
#define ADAPTIVEINTEGRATIONGRID_H

#include <vector>
#include <algorithm>
#include "CNuc.h"
#include "Config.h"

/*!
 * \brief Adaptive integration grid generator for target integration
 *
 * This class generates intelligent energy grids for target integration that
 * adaptively place more points near resonances and fewer points in smooth regions.
 * The grid generation analyzes the compound nucleus level scheme to identify
 * resonance locations and allocates integration points accordingly.
 *
 * Key features:
 * - Resonance-aware: Places finer grids near level energies
 * - Adaptive: Uses fewer points in smooth regions
 * - Bounded: Respects a maximum point limit
 * - Compatible: Designed as a drop-in replacement for uniform grids
 */
class AdaptiveIntegrationGrid {
public:
  /*!
   * \brief Configuration parameters for adaptive grid generation
   */
  struct GridConfig {
    int maxPoints;                    ///< Maximum number of integration points (upper limit)
    double baseEnergyStep;            ///< Base energy step in smooth regions (MeV)
    double resonanceWidthMultiplier;  ///< Multiplier for resonance width region (e.g., 3.0 = 3×Γ)
    double pointsPerWidth;            ///< Target number of points per resonance width
    int entranceKey;                  ///< Entrance channel key for separation energy calculation

    // Default constructor with sensible defaults
    GridConfig()
      : maxPoints(1000),
        baseEnergyStep(0.001),
        resonanceWidthMultiplier(5.0),
        pointsPerWidth(50.0),
        entranceKey(0) {}
  };

  /*!
   * \brief Grid point information
   */
  struct GridPoint {
    double energy;      ///< Energy of the point (CM frame)
    bool isResonant;    ///< Flag indicating if point is in resonant region
  };

  /*!
   * \brief Resonance information extracted from CNuc
   */
  struct ResonanceInfo {
    double energy;      ///< Resonance energy (CM frame)
    double totalWidth;  ///< Total width Γ_total = Σ Γ_i (MeV)
  };

  /*!
   * \brief Constructor
   * \param config Grid generation configuration
   */
  AdaptiveIntegrationGrid(const GridConfig& config);

  /*!
   * \brief Generate adaptive grid for a given energy range
   *
   * Analyzes the compound nucleus level scheme and creates an adaptive
   * energy grid with more points near resonances and fewer in smooth regions.
   *
   * \param startEnergy Starting energy (highest energy, CM frame)
   * \param endEnergy Ending energy (lowest energy, CM frame)
   * \param compound Pointer to compound nucleus for level information
   * \return Vector of grid points (energies in descending order)
   */
  std::vector<double> GenerateGrid(double startEnergy, double endEnergy, CNuc* compound);

  /*!
   * \brief Generate grid with detailed point information
   * \param startEnergy Starting energy (highest energy, CM frame)
   * \param endEnergy Ending energy (lowest energy, CM frame)
   * \param compound Pointer to compound nucleus for level information
   * \return Vector of GridPoint structures with energy and resonance flags
   */
  std::vector<GridPoint> GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc* compound);

  /*!
   * \brief Get the actual number of points that would be generated
   * \param startEnergy Starting energy (CM frame)
   * \param endEnergy Ending energy (CM frame)
   * \param compound Pointer to compound nucleus
   * \return Number of points in the adaptive grid
   */
  int GetExpectedPointCount(double startEnergy, double endEnergy, CNuc* compound);

  /*!
   * \brief Get the resonances (CM energy + total width) within an energy range.
   *
   * Thin public accessor around the internal resonance finder so that other
   * integrators (e.g. the reaction-rate integration) can place integration
   * breakpoints at the resonance energies.
   *
   * \param startEnergy Highest energy of the range (CM frame)
   * \param endEnergy   Lowest energy of the range (CM frame)
   * \param compound    Pointer to compound nucleus for level information
   * \return Vector of ResonanceInfo (energy, totalWidth), sorted by energy
   */
  std::vector<ResonanceInfo> GetResonances(double startEnergy, double endEnergy, CNuc* compound) {
    return IdentifyResonances(startEnergy, endEnergy, compound);
  }

  /*!
   * \brief Set grid configuration
   * \param config New configuration parameters
   */
  void SetConfig(const GridConfig& config);

  /*!
   * \brief Get current grid configuration
   * \return Current configuration
   */
  GridConfig GetConfig() const;

private:
  GridConfig config_;  ///< Grid generation configuration

  /*!
   * \brief Identify resonances within the energy range with their widths
   * \param startEnergy Starting energy (CM frame)
   * \param endEnergy Ending energy (CM frame)
   * \param compound Pointer to compound nucleus
   * \return Vector of ResonanceInfo structures with energies and widths
   */
  std::vector<ResonanceInfo> IdentifyResonances(double startEnergy, double endEnergy, CNuc* compound);

  /*!
   * \brief Check if an energy is within a resonant region
   * \param energy Energy to check
   * \param resonances Vector of resonance information
   * \return True if energy is near a resonance
   */
  bool IsInResonantRegion(double energy, const std::vector<ResonanceInfo>& resonances) const;

  /*!
   * \brief Calculate appropriate step size at a given energy
   * \param energy Current energy
   * \param resonances Vector of resonance information
   * \return Adaptive step size based on local resonance width
   */
  double CalculateAdaptiveStep(double energy, const std::vector<ResonanceInfo>& resonances) const;

  /*!
   * \brief Find the nearest resonance to a given energy
   * \param energy Energy to check
   * \param resonances Vector of resonance information
   * \param distance Output parameter for distance to nearest resonance
   * \return Pointer to nearest ResonanceInfo, or nullptr if no resonances
   */
  const ResonanceInfo* FindNearestResonance(double energy,
                                            const std::vector<ResonanceInfo>& resonances,
                                            double& distance) const;
};

#endif // ADAPTIVEINTEGRATIONGRID_H
