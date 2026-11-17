#ifndef AZURE2_EXPERIMENT_ISOSPIN_H
#define AZURE2_EXPERIMENT_ISOSPIN_H

#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * The records of <isospin> (docs/upgrade/dev_features.md section 5): isospin-labelled levels across compound nuclei,
 * Hale's charge-independent R-matrix.  A level whose <levels> lines span several compound nuclei under one level ID
 * is one isospin level of isospin T (`isospin_level`); in every nucleus its reduced-width amplitudes are the
 * Clebsch-Gordan projections of one amplitude per isospin channel, and its energy is the reference nucleus's
 * (`isospin_reference`) moved by the threshold difference and the nucleus's internal Coulomb shift (`isospin_shift`).
 * `isospin_mix` gives an isospin-forbidden channel an amplitude of its own, `isospin_factor` replaces a computed
 * factor, `isospin_pair` gives the particles of a pair their isospins where |N - Z|/2 is not it.
 *
 * Plain data in engine units (MeV, MeV^1/2), resolved against the pairs of <levels>; which level lines, channels and
 * parameters they name is the set-up's (the engine's IsospinSet, which knows the compound system).
 */

/// A number a fit may vary: `<value> fixed|fit [prior <centre> <width>]`.
struct IsospinNumber {
  double value = 0.0;
  bool fit = false;
  bool hasPrior = false;
  double priorCentre = 0.0, priorWidth = 0.0;  ///< the width in the engine unit (a percentage already applied)
};

/// `isospin_reference` (no shift) or `isospin_shift`: a compound nucleus of the isospin levels, the particle pair its
/// level energies are measured from, and the internal Coulomb shift of its isospin levels.
struct IsospinNucleusEntry {
  int Z = 0, A = 0;
  int pairKey = 0;          ///< the threshold pair, a particle pair of this nucleus
  bool reference = false;   ///< the reference nucleus: no shift, no parameter
  IsospinNumber shift;      ///< MeV (isospin_shift)
  std::string record;
  int line = 0;
};

/// `isospin_level`: the isospin of the level with this ID; with `soft` its derived parameters are free and tied by a
/// penalty ((value - tie)/sigma)^2 each.
struct IsospinLevelEntry {
  int levelId = 0;
  int twoT = 0;
  bool soft = false;
  double sigmaE = 0.0;      ///< MeV
  double sigmaGamma = 0.0;  ///< MeV^1/2
  std::string record;
  int line = 0;
};

/// `isospin_mix`: an amplitude of its own (MeV^1/2) for an isospin-forbidden channel (pair, l, s) of a level.
struct IsospinMixEntry {
  int levelId = 0, pairKey = 0, l = 0, twoS = 0;
  IsospinNumber amplitude;
  std::string record;
  int line = 0;
};

/// `isospin_factor`: the factor of the channels of a pair in a level, in place of the Clebsch-Gordan coefficient.
struct IsospinFactorEntry {
  int levelId = 0, pairKey = 0;
  double value = 0.0;
  std::string record;
  int line = 0;
};

/// `isospin_pair`: the isospins of the two particles of a pair (default |N - Z|/2 each).
struct IsospinPairEntry {
  int pairKey = 0;
  int twoT1 = 0, twoT2 = 0;  ///< light, heavy
  std::string record;
  int line = 0;
};

struct IsospinSpec {
  std::vector<IsospinNucleusEntry> nuclei;  ///< the reference and the shifts, in record order
  std::vector<IsospinLevelEntry> levels;
  std::vector<IsospinMixEntry> mixes;
  std::vector<IsospinFactorEntry> factors;
  std::vector<IsospinPairEntry> pairs;
  bool Empty() const { return nuclei.empty() && levels.empty() && mixes.empty() && factors.empty() && pairs.empty(); }
  /// Does an isospin_level record declare this level ID?  (Its <levels> lines may then span compound nuclei.)
  bool Declares(int levelId) const {
    for (const IsospinLevelEntry &e : levels)
      if (e.levelId == levelId) return true;
    return false;
  }
};

}  // namespace experiment
}  // namespace azure2

#endif
