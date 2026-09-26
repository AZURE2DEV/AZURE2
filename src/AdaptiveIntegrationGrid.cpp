#include "AdaptiveIntegrationGrid.h"
#include "JGroup.h"
#include "ALevel.h"
#include "PPair.h"
#include "ChannelFunc.h"
#include "CoulFunc.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

namespace {
bool DebugGridEnabled() {
  static const bool enabled = (std::getenv("AZR_DEBUG_GRID") != nullptr);
  return enabled;
}
}  // namespace

AdaptiveIntegrationGrid::AdaptiveIntegrationGrid(const GridConfig &config) :
  config_(config) {
}

/*!
 * \brief Round a resonance to its lattice anchor.
 *
 * The width goes to the nearest power of 1.25 and the energy to the nearest
 * multiple of a quarter of that width.  The grid built from the anchors is
 * then piecewise constant in the level parameters: it changes only when a
 * level moves by about a quarter of its width or its width changes by about
 * 25 %, so a fit that rebuilds the grid as a narrow level moves (see
 * EPoint::RefreshSubPointGrid) gets the same grid at the same parameters in
 * every thread, and a finite-difference step almost never straddles a change.
 * The peak stays within Gamma/8 of the centre of its uniform core.
 */
AdaptiveIntegrationGrid::ResonanceInfo AdaptiveIntegrationGrid::Quantize(const ResonanceInfo &resonance) {
  ResonanceInfo q = resonance;
  if (!(resonance.particleWidth > 0.0)) return q;
  const double ratio = 1.25;
  double width = std::pow(ratio, std::round(std::log(resonance.particleWidth) / std::log(ratio)));
  double quantum = 0.25 * width;
  q.particleWidth = width;
  q.totalWidth = resonance.totalWidth * width / resonance.particleWidth;
  q.energy = quantum * std::round(resonance.energy / quantum);
  return q;
}

/*!
 * \brief The quantized anchors of every level of the compound.
 */
std::vector<AdaptiveIntegrationGrid::ResonanceInfo> AdaptiveIntegrationGrid::Anchors(CNuc *compound) {
  std::vector<ResonanceInfo> anchors = LevelResonances(compound);
  for (ResonanceInfo &r : anchors) r = Quantize(r);
  return anchors;
}

/*!
 * \brief Growth of the tail step with distance: step = TailFactor() x distance.
 *
 * A chord over a 1/x^2 wing whose step is a fraction r of the distance
 * over-counts that stretch by about r^2/2; summed over the wing beyond the
 * core (which holds ~Gamma/(pi core) of the peak's area on each side) the peak
 * area comes out high by about r^2 Gamma / (2 pi core).  So the ratio can grow
 * as sqrt(core) for the same bias: tailRatio/pointsPerWidth (0.05 at 50) for a
 * 2-width core, three times that for 20 widths.
 */
double AdaptiveIntegrationGrid::TailFactor() const {
  double coreWidths = std::max(convolutionCoreWidths, config_.resonanceWidthMultiplier);
  return tailRatio / config_.pointsPerWidth * std::sqrt(coreWidths / convolutionCoreWidths);
}

/*!
 * \brief Distance out to which a resonance's lattice is finer than the smooth
 * step: the uniform core, then the geometric tail up to where its step
 * (TailFactor() times the distance) reaches baseEnergyStep.
 * Zero for a resonance whose core pitch is already coarser than the smooth
 * step (a broad or background pole): the smooth lattice covers it.
 */
double AdaptiveIntegrationGrid::Reach(const ResonanceInfo &resonance) const {
  double width = resonance.particleWidth;
  if (!(width > 0.0) || !(config_.pointsPerWidth > 0.0)) return 0.0;
  double pitch = width / config_.pointsPerWidth;
  if (!(pitch < config_.baseEnergyStep)) return 0.0;
  double core = std::max(0.0, config_.resonanceWidthMultiplier) * width;
  double tail = config_.baseEnergyStep / TailFactor();
  return std::max(core, tail);
}

std::vector<AdaptiveIntegrationGrid::ResonanceInfo>
AdaptiveIntegrationGrid::AnchorsInReach(double startEnergy, double endEnergy,
                                        const std::vector<ResonanceInfo> &anchors) const {
  if (startEnergy < endEnergy) std::swap(startEnergy, endEnergy);
  std::vector<ResonanceInfo> inReach;
  for (const ResonanceInfo &r : anchors) {
    double reach = Reach(r);
    if (!(reach > 0.0)) continue;
    if (r.energy - reach > startEnergy || r.energy + reach < endEnergy) continue;
    inReach.push_back(r);
  }
  return inReach;
}

std::vector<double> AdaptiveIntegrationGrid::GenerateGrid(double startEnergy, double endEnergy, CNuc *compound) {
  return GenerateGrid(startEnergy, endEnergy, Anchors(compound));
}

/*!
 * \brief Generate the sub-point grid of a window [endEnergy, startEnergy].
 *
 * 1. A uniform smooth lattice with the step closest to baseEnergyStep that
 *    divides the window.
 * 2. For every anchor finer than that step: its own lattice, a fixed function
 *    of (E_R, Gamma) and never of the window -- so overlapping windows sample
 *    a resonance at the same energies and the convolved curve varies smoothly
 *    from one data point to the next -- made of a uniform core of pitch
 *    Gamma/pointsPerWidth out to resonanceWidthMultiplier widths and
 *    geometric tails (step = TailFactor() x distance) out to the
 *    smooth step.  Where lattices overlap, each energy keeps only the points
 *    of the finest one there; the smooth lattice fills in where none is finer.
 * 3. Sort (descending) and remove near-duplicates.
 */
std::vector<double> AdaptiveIntegrationGrid::GenerateGrid(double startEnergy, double endEnergy,
                                                          const std::vector<ResonanceInfo> &anchors) const {
  std::vector<double> grid;
  if (startEnergy < endEnergy) std::swap(startEnergy, endEnergy);
  double totalRange = startEnergy - endEnergy;
  if (totalRange <= 0.0 || totalRange < 1.0e-6) {
    grid.push_back(startEnergy);
    return grid;
  }

  double base = config_.baseEnergyStep;
  if (!(base > 0.0) || base > totalRange) base = totalRange;
  int nSmooth = (int)std::ceil(totalRange / base - 1.0e-9);
  if (nSmooth < 1) nSmooth = 1;

  std::vector<ResonanceInfo> resonances = AnchorsInReach(startEnergy, endEnergy, anchors);
  const double tailFactor = TailFactor();

  // Spacing of resonance r's lattice at energy e (base where it has none).
  // Where lattices overlap, only the finest one is kept: interleaving two
  // lattices adds points without resolving anything either does not (three
  // 35-55 keV poles with 20-width cores tripled the sub-points of
  // tests/hybrid_potential).
  auto localStep = [&](const ResonanceInfo &r, double e) {
    double pitch = r.particleWidth / config_.pointsPerWidth;
    double d = std::abs(e - r.energy);
    double step = (d <= std::max(0.0, config_.resonanceWidthMultiplier) * r.particleWidth * (1.0 + 1.0e-9))
        ? pitch : std::max(pitch, tailFactor * d);
    return std::min(step, base);
  };
  auto finestStep = [&](double e) {
    double s = base;
    for (const ResonanceInfo &r : resonances) s = std::min(s, localStep(r, e));
    return s;
  };

  // The smooth lattice, where no resonance lattice is finer (the edges always).
  for (int i = 0; i <= nSmooth; i++) {
    double e = startEnergy - totalRange * i / nSmooth;
    if (i == 0 || i == nSmooth || !(finestStep(e) < 0.5 * base)) grid.push_back(e);
  }

  static const bool debugGrid = DebugGridEnabled();
  if (debugGrid) {
    fprintf(stderr, "[AZR_DEBUG_GRID] window [%.10f, %.10f] MeV, smooth step %.6e MeV, %zu anchors in reach:\n",
            endEnergy, startEnergy, totalRange / nSmooth, resonances.size());
    for (const ResonanceInfo &r : resonances)
      fprintf(stderr, "[AZR_DEBUG_GRID]   anchor id=%d E=%.10f MeV particleWidth=%.6e MeV\n",
              r.id, r.energy, r.particleWidth);
  }

  double finestPitch = base;
  for (const ResonanceInfo &res : resonances) {
    const double width = res.particleWidth;
    const double pitch = width / config_.pointsPerWidth;
    const double core = std::max(0.0, config_.resonanceWidthMultiplier) * width;
    finestPitch = std::min(finestPitch, pitch);
    // Farthest window edge from the resonance: nothing beyond it is kept.
    const double farthest = std::max(startEnergy - res.energy, res.energy - endEnergy);
    // Nearest window edge: offsets below it fall outside on both sides.
    const double nearest = (res.energy > startEnergy) ? res.energy - startEnergy
                         : (res.energy < endEnergy) ? endEnergy - res.energy : 0.0;
    auto keep = [&](double e) {
      if (!(e >= endEnergy && e <= startEnergy)) return;
      double own = localStep(res, e);
      for (const ResonanceInfo &other : resonances) {
        if (&other == &res) continue;
        double s = localStep(other, e);
        if (s < own * (1.0 - 1.0e-9) || (s <= own * (1.0 + 1.0e-9) && &other < &res)) return;
      }
      grid.push_back(e);
    };
    double offset = 0.0;
    long nCore = (long)std::floor(core / pitch + 1.0e-9);
    // Core: offset = k * pitch.  Start at the first index that can reach the
    // window (an anchor outside it contributes only its inner flank).
    long k0 = (long)std::floor(nearest / pitch);
    if (k0 > nCore) k0 = nCore + 1;
    for (long k = k0; k <= nCore; k++) {
      offset = k * pitch;
      if (offset > farthest) break;
      keep(res.energy + offset);
      if (k > 0) keep(res.energy - offset);
    }
    offset = nCore * pitch;
    // Geometric tails.
    while (offset <= farthest) {
      double step = std::max(pitch, tailFactor * offset);
      if (!(step < base)) break;
      offset += step;
      if (offset >= nearest) {
        keep(res.energy + offset);
        keep(res.energy - offset);
      }
    }
  }

  std::sort(grid.begin(), grid.end(), std::greater<double>());

  // Remove near-duplicates.  The tolerance tracks the finest pitch in use, so
  // a sub-eV lattice is not collapsed; exact duplicates always go.
  double tolerance = std::min(1.0e-8, finestPitch * 1.0e-3);
  if (tolerance < 1.0e-14) tolerance = 1.0e-14;
  std::vector<double> uniqueGrid;
  for (double e : grid) {
    if (uniqueGrid.empty() || std::abs(e - uniqueGrid.back()) > tolerance) uniqueGrid.push_back(e);
  }
  // The window edges are exact grid points (the quadrature clips to them).
  uniqueGrid.front() = startEnergy;
  uniqueGrid.back() = endEnergy;
  if (debugGrid) fprintf(stderr, "[AZR_DEBUG_GRID]   %zu grid points\n", uniqueGrid.size());
  return uniqueGrid;
}

std::vector<AdaptiveIntegrationGrid::GridPoint>
AdaptiveIntegrationGrid::GenerateGridWithInfo(double startEnergy, double endEnergy, CNuc *compound) {
  std::vector<ResonanceInfo> anchors = Anchors(compound);
  std::vector<double> energyGrid = GenerateGrid(startEnergy, endEnergy, anchors);
  std::vector<ResonanceInfo> inReach = AnchorsInReach(startEnergy, endEnergy, anchors);
  std::vector<GridPoint> gridInfo;
  for (double energy : energyGrid) {
    GridPoint point;
    point.energy = energy;
    point.isResonant = false;
    for (const ResonanceInfo &r : inReach)
      if (std::abs(energy - r.energy) <= Reach(r)) point.isResonant = true;
    gridInfo.push_back(point);
  }
  return gridInfo;
}

int AdaptiveIntegrationGrid::GetExpectedPointCount(double startEnergy, double endEnergy, CNuc *compound) {
  return (int)GenerateGrid(startEnergy, endEnergy, compound).size();
}

void AdaptiveIntegrationGrid::SetConfig(const GridConfig &config) {
  config_ = config;
}

AdaptiveIntegrationGrid::GridConfig AdaptiveIntegrationGrid::GetConfig() const {
  return config_;
}

/*!
 * \brief The levels within 20 widths of an energy range, unquantized (for the
 * reaction-rate breakpoints).
 */
std::vector<AdaptiveIntegrationGrid::ResonanceInfo>
AdaptiveIntegrationGrid::IdentifyResonances(double startEnergy, double endEnergy, CNuc *compound) {
  if (startEnergy < endEnergy) std::swap(startEnergy, endEnergy);
  std::vector<ResonanceInfo> resonances;
  for (const ResonanceInfo &r : LevelResonances(compound)) {
    double margin = 20.0 * r.particleWidth;
    if (r.energy >= endEnergy - margin && r.energy <= startEnergy + margin) resonances.push_back(r);
  }
  return resonances;
}

/*!
 * \brief Every level of the compound as a resonance, with its widths
 *
 * Loops through all J-groups and levels in the compound nucleus. For each resonance,
 * calculates the total width by summing all partial widths (Γ_total = Σ Γ_i).
 * Converts level energies from compound excitation energy to CM energy using
 * the separation energy of the entrance channel.
 */
std::vector<AdaptiveIntegrationGrid::ResonanceInfo>
AdaptiveIntegrationGrid::LevelResonances(CNuc *compound) {
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

  // The level parameters to read: the stored ones, or the current fit values
  // (observed energies and Brune amplitudes during a fit).
  const bool fit = config_.useFitParameters;
  auto levelE = [fit](ALevel *lv) { return fit ? lv->GetFitE() : lv->GetE(); };
  auto levelGamma = [fit](ALevel *lv, int ch) { return fit ? lv->GetFitGamma(ch) : lv->GetGamma(ch); };

  // Loop through all J-groups
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jgroup = compound->GetJGroup(j);
    if (!jgroup || !jgroup->IsInRMatrix()) continue;

    int numChannels = jgroup->NumChannels();

    // Loop through all levels in this J-group
    for (int l = 1; l <= jgroup->NumLevels(); l++) {
      ALevel *level = jgroup->GetLevel(l);
      if (!level || !level->IsInRMatrix()) continue;

      // Get level energy (in compound excitation energy)
      double levelExcitationEnergy = levelE(level);

      // Check if level has widths (Γ > 0) to be considered a resonance
      bool hasWidth = false;
      for (int ch = 1; ch <= numChannels; ch++) {
        if (std::abs(levelGamma(level, ch)) > 1.0e-6) {
          hasWidth = true;
          break;
        }
      }
      if (!hasWidth) continue;  // skip levels with no widths

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
        // A channel flagged gammaIsRWA (the optional 33rd field of a level
        // line, used by THM projects) holds a reduced-width amplitude in
        // MeV^1/2, not a width in eV.  Summing it as eV made such a resonance
        // look orders of magnitude narrower than it is: the grid then packed
        // its points into a sliver around the peak and bridged the real flanks
        // with keV-wide linear interpolation, over-counting the resonance area
        // tenfold in every convolution or target integration across it.
        // Convert those channels exactly as the amplitude branch below does.
        double rwaParticleMeV = 0.0;
        double rwaNormSum = 0.0;
        double rwaGammaMeV = 0.0;
        for (int ch = 1; ch <= numChannels; ch++) {
          AChannel *channel = jgroup->GetChannel(ch);
          PPair *chPair = compound->GetPair(channel->GetPairNum());
          double gamma = std::abs(levelGamma(level, ch));
          if (gamma <= 0.0) continue;
          const bool isRWA = level->GammaIsRWA(ch);
          if (channel->GetRadType() == 'P') {
            double localEnergy = levelE(level) - chPair->GetExE() - chPair->GetSepE();
            // Closed channels -- below or at threshold (the value is an ANC
            // unless flagged), or just above a Coulomb threshold where P = 0
            // in double precision (ChannelFunc) -- are left out of this
            // estimate, on both sides of threshold alike, so it does not jump
            // when a level crosses it.  (Their gamma^2 dS/dE would belong in
            // the normalization, as in CNuc::TransformOut; the estimate only
            // sizes the integration grid, and adding it moves the grid of
            // existing projects -- tests/17O by 1.4e-3 in chi2.)
            ChannelFunc channelFunc(chPair, false);
            if (channelFunc.IsClosed(localEnergy)) continue;
            if (isRWA) {
              double pene = channelFunc.Penetrability(channel->GetL(), localEnergy);
              double dSdE = channelFunc.ShiftDerivative(channel->GetL(), localEnergy);
              rwaParticleMeV += 2.0 * gamma * gamma * pene;
              rwaNormSum += dSdE * gamma * gamma;
            } else {
              particleWidthEV += gamma;
            }
          } else if (channel->GetRadType() == 'M' || channel->GetRadType() == 'E') {
            if (isRWA) {
              double localEnergy = levelE(level) - chPair->GetExE() - chPair->GetSepE();
              rwaGammaMeV += 2.0 * gamma * gamma *
                             pow(std::abs(localEnergy) / hbarc, 2.0 * channel->GetL() + 1.0);
            } else {
              gammaWidthEV += gamma;
            }
          }
        }
        if (1.0 + rwaNormSum > 0.0) rwaParticleMeV /= (1.0 + rwaNormSum);
        double particleWidth = particleWidthEV * 1.0e-6 + rwaParticleMeV;
        double totalWidthMeV = particleWidth + gammaWidthEV * 1.0e-6 + rwaGammaMeV;
        if (DebugGridEnabled()) {
          fprintf(stderr, "[AZR_DEBUG_GRID]     level E=%.6f: input widths particle=%.6e eV gamma=%.6e eV "
                  "(total from RWA-flagged channels %.6e eV)\n",
                  levelE(level), particleWidthEV, gammaWidthEV, (rwaParticleMeV + rwaGammaMeV) * 1.0e6);
        }
        ResonanceInfo resInfo;
        resInfo.energy = levelCMEnergy;
        resInfo.totalWidth = totalWidthMeV;
        resInfo.particleWidth = particleWidth;
        resInfo.id = 1000 * j + l;
        resonances.push_back(resInfo);
        continue;
      }
      for (int ch = 1; ch <= numChannels; ch++) {
        AChannel *channel = jgroup->GetChannel(ch);
        if (channel->GetRadType() != 'P') continue;  // particle channels
        PPair *chPair = compound->GetPair(channel->GetPairNum());
        double localEnergy = levelE(level) - chPair->GetExE() - chPair->GetSepE();
        // Closed channels are left out, as above (ChannelFunc: P = 0 also just
        // above a Coulomb threshold).
        ChannelFunc channelFunc(chPair, false);
        if (channelFunc.IsClosed(localEnergy)) continue;
        double gamma = std::abs(levelGamma(level, ch));
        if (gamma <= 0.0) continue;
        double radius = chPair->GetChRad();
        double pene = channelFunc.Penetrability(channel->GetL(), localEnergy);
        double dSdE = channelFunc.ShiftDerivative(channel->GetL(), localEnergy);
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
        double gamma = std::abs(levelGamma(level, ch));
        if (gamma <= 0.0) continue;
        double localEnergy = levelE(level) - chPair->GetExE() - chPair->GetSepE();
        double pene;
        if (std::abs(levelE(level) - chPair->GetExE()) < 1.0e-3 &&
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
      // is sized by it.
      ResonanceInfo resInfo;
      resInfo.energy = levelCMEnergy;
      resInfo.totalWidth = totalWidth;
      resInfo.particleWidth = particleWidth;
      resInfo.id = 1000 * j + l;
      resonances.push_back(resInfo);
    }
  }

  // Sort resonances by energy in ascending order for efficient searching
  std::sort(resonances.begin(), resonances.end(),
            [](const ResonanceInfo &a, const ResonanceInfo &b) {
              return a.energy < b.energy;
            });

  return resonances;
}

