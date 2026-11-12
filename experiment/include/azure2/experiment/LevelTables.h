#ifndef AZURE2_EXPERIMENT_LEVELTABLES_H
#define AZURE2_EXPERIMENT_LEVELTABLES_H

#include <string>
#include <vector>

#include "azure2/experiment/Values.h"

namespace azure2 {
namespace experiment {

/*!
 * The pairs and levels of a project's <levels> section as plain tables, and the one rule that finds the pair a label
 * names (docs/upgrade/project_format3.md section 3.5).  The tables are made once from the <levels> rows
 * (LevelTablesOfRows, ProjectReader.h; tests write them by hand) and read by Schema::Build (data sets,
 * extrapolations) and BuildModelData (model data); the engine's model data use the same tables; the converter writes
 * its labels with PairLabel.  No engine headers.
 */

/// One pair of <levels> (the columns of its first row; the engine's PPair).
struct PairInfo {
  int key = 0;   ///< the pair key of <levels>
  int kind = 0;  ///< 0 particle pair, 10 photon (capture to a level of the compound nucleus), 20 beta
  int lightZ = 0, lightA = 0;  ///< the light particle (a photon: 0, 0; a beta: Z = -1 or +1, A = 0)
  int heavyZ = 0, heavyA = 0;  ///< the heavy particle (photon pair: the compound nucleus)
  double lightExcitation = 0.0;  ///< MeV
  double heavyExcitation = 0.0;  ///< MeV (photon pair: E_x of the final level)
  /// The two excitations as written in <levels> ("" = not known: the label formats the value)
  std::string lightExcitationText, heavyExcitationText;
  double threshold = 0.0;  ///< MeV, the E_x of the compound nucleus at which the pair opens (S + e1 + e2)
  double lightMass = 0.0, heavyMass = 0.0;  ///< u, the masses of <levels> (PPair's), for the frame of a c.m. grid
  double separationEnergy = 0.0;  ///< MeV, S of <levels>
  double radius = 0.0;            ///< fm, the channel radius of its first row
  int twoJLight = 0, twoJHeavy = 0;  ///< 2j of the two particles
  int heavyParity = 1;               ///< parity of the heavy particle (photon pair: of the final level)
  /// The external-capture multipole mask of the pair (column 31 of its first active row that has one; 0 = none)
  unsigned int externalCapture = 0;
  bool identical = false;            ///< the two particles are identical (PPair::IsIdentical)
};

/// One channel of a level's J^pi group.
struct LevelChannelInfo {
  int pairKey = 0;
  int l = 0;
  int twoS = 0;
};

/// A level of <levels> (the engine's ALevel and its J^pi group).
struct LevelInfo {
  std::vector<int> ids;  ///< the level IDs (column 9) of the <levels> lines read into this level
  int twoJ = 0;
  int parity = 0;
  double energy = 0.0;  ///< MeV, E_x as written
  /// MeV, the total width from the widths written in <levels> (the start of a fit), for the warnings; < 0: unknown
  double width = -1.0;
  int group = 0, level = 0;  ///< the engine's indices of the J^pi group and of the level in it (1-based)
  std::vector<LevelChannelInfo> channels;  ///< the channels of its J^pi group
  int nucleusZ = 0, nucleusA = 0;  ///< the compound nucleus of its J^pi group
};

/// Where one row of <levels> belongs: the engine's indices (1-based) of its J^pi group, of its level in the group and
/// of its channel in the group's channel list; all 0 for an inactive row (the engine does not read it).
struct LevelRowSlot {
  int group = 0, level = 0, channel = 0;
  double value = 0.0;  ///< column 12 as written: a width (eV), an ANC (fm^-1/2) or an amplitude (MeV^1/2)
};

/// The tables of one <levels> section.
struct LevelTables {
  std::vector<PairInfo> pairs;    ///< one per pair key, from its first row (active or not), in row order
  std::vector<LevelInfo> levels;  ///< the levels of the active rows, by group and level number
  /// The keys of the pairs with an external-capture mask on an active row, in the engine's pair order (the order of
  /// their first active rows)
  std::vector<int> externalCapturePairs;
  std::vector<LevelRowSlot> rowSlots;  ///< one per row, in row order
};

/// "15N", "4He": mass number and element symbol.
std::string NuclideLabel(int Z, int A);
/// The light particle as a label writes it: g (photon pair), b- / b+ (beta pair), n p d t h a, else the nuclide.
std::string LightParticleLabel(int kind, int Z, int A);
/// The pair as a label writes it: "p+15N", "a[3.0]+12C[4.439]" (the excitations as written when known).
std::string PairLabel(const PairInfo &pair);
/// The label of the reaction from pair `in` to pair `out`, "15N(p,a)12C" (states as PairLabel writes them), or with
/// `out` null the total capture of `in`, "12C(p,g)13N[total]".  The converter writes its reactions with it.
std::string ReactionLabel(const PairInfo &in, const PairInfo *out);
/// "key: label" for every pair, the list a refusal shows.
std::vector<std::string> PairListing(const std::vector<PairInfo> &pairs);

const PairInfo *PairByKey(int key, const std::vector<PairInfo> &pairs);

/// The keys of the pairs whose light particle is `light`, heavy particle `heavy` in state `state`: particles by Z
/// and A (g, b-, b+ by the pair kind), states to the digits written (a written 4.439 matches 4.43891 MeV of <levels>;
/// no state written is the ground state, excitation 0).
std::vector<int> PairsNamed(const Particle &light, const Particle &heavy, const ResidualState &state,
                            const std::vector<PairInfo> &pairs);

/// The messages of E_UNKNOWN_PAIR and E_AMBIGUOUS_PAIR (one wording for data sets and model data): `what` is the
/// side of the reaction or the pair as written ("the entrance of 15N(p,g)16O"); `remedy` ends the second.
std::string UnknownPairMessage(const std::string &what);
std::string AmbiguousPairMessage(const std::string &what, const std::vector<int> &keys, const std::string &remedy);

}  // namespace experiment
}  // namespace azure2

#endif
