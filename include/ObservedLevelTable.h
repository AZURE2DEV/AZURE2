#ifndef OBSERVED_LEVEL_TABLE_H
#define OBSERVED_LEVEL_TABLE_H

#include <vector>

#include "Constants.h"

/*!
 * The Brune (Park) and observed values of every level at one set of fit values: what CNuc::TransformOut computes for
 * the reports (parameters.out, the model data, the warnings, the observed parameter vector).  A value: the
 * transformation returns it and changes no level; its caller holds it as long as it reads it.
 *
 * Indexed like the compound's levels: J-group j and level la, both 1-based, and channels of the group 1-based.
 */
struct ObservedLevel {
  double E = 0.0;              ///< Brune (Park) energy, MeV (E_x)
  int iterations = 0;          ///< iterations of the formal -> Brune conversion (0 in the Brune and Park bases)
  double sqrtNF = 1.0;         ///< N_f^-1/2 formed with the Brune (Park) amplitudes; external capture reads it
  vector_r gamma;              ///< Brune (Park) amplitude per channel, MeV^1/2
  vector_c external;           ///< external part of a photon channel's amplitude, MeV^1/2 (0 for other channels)
  vector_r W;                  ///< observed partial width per channel, signed, MeV (physics/ observed_levels)
  vector_r value;              ///< the observed value as a <levels> line holds it (Gamma in eV, C, mu, Q, B)
  std::vector<unsigned char> kind;  ///< azure2::physics::ObservedKind of `value`
};

class ObservedLevelTable {
 public:
  ObservedLevelTable() = default;
  /// One empty entry per level: numLevels[j - 1] levels of numChannels[j - 1] channels in J-group j.
  ObservedLevelTable(const std::vector<int> &numLevels, const std::vector<int> &numChannels);
  const ObservedLevel &Level(int j, int la) const { return levels_.at(j - 1).at(la - 1); }
  ObservedLevel &Level(int j, int la) { return levels_.at(j - 1).at(la - 1); }

 private:
  std::vector<std::vector<ObservedLevel>> levels_;
};

#endif
