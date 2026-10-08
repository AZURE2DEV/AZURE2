#ifndef ADAPTIVEINTEGRATIONGRID_H
#define ADAPTIVEINTEGRATIONGRID_H

#include <vector>
#include <algorithm>
#include "CNuc.h"
#include "Config.h"

/*!
 * \brief Adaptive integration grid generator for target integration
 *
 * Generates the sub-point energies of a target effect (Gaussian convolution,
 * target integration, beam profile).  The window is covered by a uniform
 * smooth lattice (step = baseEnergyStep) and, for every resonance narrower
 * than that step, by a lattice anchored at the resonance itself:
 *
 *  - a uniform core of pitch Gamma/pointsPerWidth out to
 *    +-resonanceWidthMultiplier*Gamma, then
 *  - geometric tails whose step grows in proportion to the distance to the
 *    resonance (see TailFactor), until it reaches the smooth step.
 *
 * The tails keep the relative spacing constant in the 1/(E-E_R)^2 wings of a
 * narrow peak, so the piecewise-linear quadrature over the sub-points
 * integrates each isolated resonance's area to a fixed relative accuracy
 * however narrow it is and however far the window extends.  A lattice that
 * stopped at a fixed number of widths and bridged the wings with the smooth
 * step instead over-counted a 33 eV THM resonance under a 21 keV Gaussian by up
 * to 16 % (the chord over a convex 1/x^2 wing lies above it).
 *
 * The anchors are quantized (see Quantize) so the lattice is a deterministic,
 * piecewise-constant function of the level parameters: the grid can then be
 * rebuilt when a fit moves a narrow level (EPoint::RefreshSubPointGrid)
 * without each evaluation seeing a slightly different grid.
 */
class AdaptiveIntegrationGrid {
 public:
  /*!
   * \brief Configuration parameters for adaptive grid generation
   */
  struct GridConfig {
    int maxPoints;                    ///< Number of smooth intervals across the window (targetInt field 3)
    double baseEnergyStep;            ///< Smooth step (MeV); the caller sets window / maxPoints
    double resonanceWidthMultiplier;  ///< Half-extent of the uniform core, in widths
    double pointsPerWidth;            ///< Core points per width; also sets the tail growth
    int entranceKey;                  ///< Entrance channel key for separation energy calculation
    /// The level gammas are still the input values (partial widths in eV for
    /// open particle channels, ANCs for closed ones, gamma widths in eV): the
    /// parameter transformation is enabled but has not run yet.  This is the
    /// state during EData::Fill in every CLI and API flow.
    bool inputWidthsArePhysical;
    /// Read the current fit values (ALevel::GetFitE / GetFitGamma: observed
    /// energies and Brune amplitudes) instead of the stored level parameters.
    bool useFitParameters;
    /// With useFitParameters: the fit values are formal R-matrix parameters
    /// (E_lambda, gamma; Brune formalism off), not observed energies and Brune
    /// amplitudes.  Each level is then anchored at the Thomas estimate of its
    /// observed energy, E_lambda - sum_c gamma_c^2 (S_c - B_c) / (1 + sum_c
    /// gamma_c^2 dS_c/dE), with the widths taken there.
    bool formalParameters;
    /// With useFitParameters: the fit amplitudes are Park's (--use-park),
    /// gamma_Park = gamma_Brune sqrt(J), J = 1 - sum_c gamma_c^2 dS_c/dE.  The
    /// widths are then those of the same Brune parameters, so that both
    /// parametrizations anchor the grid at the same places.
    bool parkAmplitudes;

    GridConfig() :
      maxPoints(1000),
      baseEnergyStep(0.001),
      // 20, not 5: before the lattice had geometric tails it had to span
      // enough energy for the sub-point grid to resolve the resonance; with
      // correct widths the tests/13N and tests/hybrid_potential integrals were
      // 4 % to 2x off the converged answer at 5.  A Gaussian convolution needs
      // far less (see convolutionCoreWidths).
      resonanceWidthMultiplier(20.0),
      pointsPerWidth(50.0),
      entranceKey(0),
      inputWidthsArePhysical(false),
      useFitParameters(false),
      formalParameters(false),
      parkAmplitudes(false) {}
  };

  /// Tail step / distance for a convolutionCoreWidths core: tailRatio /
  /// pointsPerWidth (0.05 at the default 50); see TailFactor.
  static constexpr double tailRatio = 2.5;

  /// Default core half-extent (widths) of a pure Gaussian convolution (no
  /// target integration, no beam profile) whose targetInt line leaves it
  /// out.  The kernel is smooth on the scale of the sub-point step, so the
  /// core only has to hold the peak: with the geometric tails, 2 widths at 50
  /// points per width fold an isolated resonance to ~3e-4 (tests/17O,
  /// tests/thm_narrow_fold) with a fifth of the sub-points of 20 widths.  A
  /// thick-target integral keeps 20: its integrand can be steep away from any
  /// resonance (tests/hybrid_potential integrates Rutherford scattering down to
  /// 1 keV), and the wide core is what samples it there.
  static constexpr double convolutionCoreWidths = 2.0;

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
    int id;                ///< 1000 * J-group + level: identifies the level
  };

  AdaptiveIntegrationGrid(const GridConfig &config);

  /*!
   * \brief Generate the grid for a window from the compound's levels.
   * \return Energies in descending order (CM frame)
   */
  std::vector<double> GenerateGrid(double startEnergy, double endEnergy, CNuc *compound);

  /*!
   * \brief Generate the grid for a window from already quantized anchors
   * (the output of Anchors()).  Only anchors whose lattice reaches the window
   * and is finer than the smooth step contribute.
   */
  std::vector<double> GenerateGrid(double startEnergy, double endEnergy,
                                   const std::vector<ResonanceInfo> &anchors) const;

  /*!
   * \brief The quantized resonance anchors of every level of the compound.
   *
   * \param frameShift  c.m. energy (MeV) subtracted from every level energy
   * before quantizing: the anchors of a segment whose energies are shifted by
   * frameShift, in the frame of its unshifted data (and so of its sub-point
   * grid), quantized there.
   */
  std::vector<ResonanceInfo> Anchors(CNuc *compound, double frameShift = 0.0);

  /*!
   * \brief The anchors that shape the grid of a window (see GenerateGrid).
   */
  std::vector<ResonanceInfo> AnchorsInReach(double startEnergy, double endEnergy,
                                            const std::vector<ResonanceInfo> &anchors) const;

  /*!
   * \brief Round a resonance to its anchor: the width to a power of 1.25 and
   * the energy to a quarter of that width.  The lattice then only changes
   * when a level moves by a quarter width or its width by 25 %, and a peak
   * never sits more than Gamma/8 from the centre of its uniform core.
   */
  static ResonanceInfo Quantize(const ResonanceInfo &resonance);

  std::vector<GridPoint> GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc *compound);

  int GetExpectedPointCount(double startEnergy, double endEnergy, CNuc *compound);

  /*!
   * \brief Get the resonances (CM energy + total width) within an energy range.
   *
   * Thin public accessor around the internal resonance finder so that other
   * integrators (e.g. the reaction-rate integration) can place integration
   * breakpoints at the resonance energies.  Not quantized.
   */
  std::vector<ResonanceInfo> GetResonances(double startEnergy, double endEnergy, CNuc *compound) {
    return IdentifyResonances(startEnergy, endEnergy, compound);
  }

  void SetConfig(const GridConfig &config);
  GridConfig GetConfig() const;

 private:
  GridConfig config_;  ///< Grid generation configuration

  /// Every level with a width, unquantized, sorted by energy.
  std::vector<ResonanceInfo> LevelResonances(CNuc *compound);

  /// The levels within resonanceWidthMultiplier widths of the range.
  std::vector<ResonanceInfo> IdentifyResonances(double startEnergy, double endEnergy, CNuc *compound);

  /// Distance out to which a resonance's lattice is finer than the smooth step.
  double Reach(const ResonanceInfo &resonance) const;

  /// Tail step per unit distance from the resonance.
  double TailFactor() const;
};

#endif  // ADAPTIVEINTEGRATIONGRID_H
