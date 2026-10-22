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
    /// The level gammas are still the input values (partial widths in eV for
    /// open particle channels, ANCs for closed ones, gamma widths in eV): the
    /// parameter transformation is enabled but has not run yet.  This is the
    /// state during EData::Fill in every CLI and API flow.
    bool inputWidthsArePhysical;
    /// Subtracted from every resonance's c.m. energy.  The sub-points of a segment
    /// with an energy shift are laid out at the unshifted energies and moved with
    /// the data; offsetting the resonances by the segment's c.m. shift puts the
    /// dense part of the grid where the resonance is after the move
    /// (EData::ReplaceSubPointGrids).
    double resonanceOffset;
    /// Beyond the fine lattice (resonanceWidthMultiplier widths) the nodes of an
    /// anchored resonance go on with steps of this fraction of the distance to it,
    /// until the step reaches baseEnergyStep (see GenerateGrid).
    double tailStepFraction;

    // Default constructor with sensible defaults
    GridConfig() :
      maxPoints(1000),
      baseEnergyStep(0.001),
      // 20, not 5: the fine lattice has to span enough energy for the sub-point
      // grid to resolve the resonance. 5 widths was only ever adequate because
      // the width estimate was wrong and enormous, which over-gridded
      // everything; with correct widths the tests/13N and tests/hybrid_potential
      // integrals are 4 % to 2x off the converged answer at 5, and flat from 20
      // out to 200.
      resonanceWidthMultiplier(20.0),
      pointsPerWidth(50.0),
      entranceKey(0),
      inputWidthsArePhysical(false),
      resonanceOffset(0.0),
      tailStepFraction(0.1) {}
  };

  /*!
   * \brief Grid point information
   */
  struct GridPoint {
    double energy;    ///< Energy of the point (CM frame)
    bool isResonant;  ///< Flag indicating if point is in resonant region
  };

  /*!
   * \brief Resonance information extracted from CNuc
   */
  struct ResonanceInfo {
    double energy;         ///< Resonance energy (CM frame, MeV)
    double totalWidth;     ///< Gamma_particle + Gamma_radiative (MeV)
    double particleWidth;  ///< Gamma_particle only (MeV); shapes sigma(E) in non-RMC
  };

  /*!
   * \brief Constructor
   * \param config Grid generation configuration
   */
  AdaptiveIntegrationGrid(const GridConfig &config);

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
  std::vector<double> GenerateGrid(double startEnergy, double endEnergy, CNuc *compound);

  /*!
   * \brief The same grid for a given list of resonances (c.m. energies and widths, MeV),
   * as IdentifyResonances would return them; used by the CNuc overload and by
   * tests/reference/adaptive_grid_reference_test.
   */
  std::vector<double> GenerateGrid(double startEnergy, double endEnergy, const std::vector<ResonanceInfo> &resonances);

  /*!
   * \brief Generate grid with detailed point information
   * \param startEnergy Starting energy (highest energy, CM frame)
   * \param endEnergy Ending energy (lowest energy, CM frame)
   * \param compound Pointer to compound nucleus for level information
   * \return Vector of GridPoint structures with energy and resonance flags
   */
  std::vector<GridPoint> GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc *compound);

  /*!
   * \brief Get the actual number of points that would be generated
   * \param startEnergy Starting energy (CM frame)
   * \param endEnergy Ending energy (CM frame)
   * \param compound Pointer to compound nucleus
   * \return Number of points in the adaptive grid
   */
  int GetExpectedPointCount(double startEnergy, double endEnergy, CNuc *compound);

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
  std::vector<ResonanceInfo> GetResonances(double startEnergy, double endEnergy, CNuc *compound) {
    return IdentifyResonances(startEnergy, endEnergy, compound);
  }

  /*!
   * \brief Set grid configuration
   * \param config New configuration parameters
   */
  void SetConfig(const GridConfig &config);

  /*!
   * \brief Get current grid configuration
   * \return Current configuration
   */
  GridConfig GetConfig() const;

 private:
  GridConfig config_;  ///< Grid generation configuration

  /// The nodes anchored to one resonance whose lattice is finer than the coarse step:
  /// E_R + n pitch for |n| <= numLattice, then E_R +- graded[k], all within E_R +- half.
  struct Zone {
    double energy;                ///< E_R (c.m., MeV)
    double pitch;                 ///< particleWidth / pointsPerWidth
    int numLattice;               ///< lattice nodes on each side
    std::vector<double> graded;   ///< offsets beyond the lattice, ascending
    double half;                  ///< the outermost offset
  };
  /// Graded offsets beyond a lattice ending at \p latticeTop, for \p pitch.
  std::vector<double> GradedOffsets(double pitch, double latticeTop) const;
  /// Half-width of the zone of a resonance of particle width \p particleWidth (its
  /// lattice and graded nodes); resonanceWidthMultiplier * width for one without a zone.
  double ZoneHalfWidth(double particleWidth) const;
  std::vector<Zone> AnchoredZones(const std::vector<ResonanceInfo> &resonances) const;
  /// The largest node of \p z strictly below \p energy; -infinity if there is none.
  static double NextNodeBelow(const Zone &z, double energy);

  /*!
   * \brief Identify resonances within the energy range with their widths
   * \param startEnergy Starting energy (CM frame)
   * \param endEnergy Ending energy (CM frame)
   * \param compound Pointer to compound nucleus
   * \return Vector of ResonanceInfo structures with energies and widths
   */
  std::vector<ResonanceInfo> IdentifyResonances(double startEnergy, double endEnergy, CNuc *compound);

  /*!
   * \brief Check if an energy is within a resonant region
   * \param energy Energy to check
   * \param resonances Vector of resonance information
   * \return True if energy is near a resonance
   */
  bool IsInResonantRegion(double energy, const std::vector<ResonanceInfo> &resonances) const;

  /*!
   * \brief Calculate appropriate step size at a given energy
   * \param energy Current energy
   * \param resonances Vector of resonance information
   * \return Adaptive step size based on local resonance width
   */
  double CalculateAdaptiveStep(double energy, const std::vector<ResonanceInfo> &resonances) const;

  /*!
   * \brief Find the nearest resonance to a given energy
   * \param energy Energy to check
   * \param resonances Vector of resonance information
   * \param distance Output parameter for distance to nearest resonance
   * \return Pointer to nearest ResonanceInfo, or nullptr if no resonances
   */
  const ResonanceInfo *FindNearestResonance(double energy,
                                            const std::vector<ResonanceInfo> &resonances,
                                            double &distance) const;
};

#endif  // ADAPTIVEINTEGRATIONGRID_H
