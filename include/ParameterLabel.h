#ifndef PARAMETERLABEL_H
#define PARAMETERLABEL_H

#include <string>

class CNuc;
class EData;
class JGroup;
class ALevel;

/*!
 * Human-readable names for levels, channels and fit parameters.
 *
 * Diagnostics that identify a parameter only by an internal index -- "j=2 la=1
 * ch=3", or the Minuit name "width_1_2" -- leave the user to work out which
 * resonance and which channel is meant.  These helpers turn the same objects
 * into something that can be found in the GUI: spin-parity, level energy,
 * particle pair, L and S.
 *
 * Every function keeps the raw indices in its output as well, since those are
 * what appear in param.sav, parameters.out and the check files.
 */
namespace AZURELabel {

/*! Spin as a fraction: 1.5 -> "3/2", 2.0 -> "2". */
std::string Spin(double j);

/*! Spin-parity: (1.5, +1) -> "3/2+". */
std::string JPi(double j, int parity);

/*! Nuclide from charge and mass: (0,1) -> "n", (1,2) -> "d", (2,4) -> "4He". */
std::string Nuclide(int z, double mass);

/*! Particle pair by 1-based pair number: "d+t".  Photon pairs give "gamma". */
std::string Pair(CNuc *compound, int pairNum);

/*! e.g. "J^pi=3/2+, E=16.9285 MeV (J-group 1, level 1)". */
std::string Level(JGroup *jgroup, ALevel *level, int jgroupIndex, int levelIndex);

/*! e.g. "channel 2: d+t, L=0, S=3/2", or "channel 5: E1 gamma to n+4He". */
std::string Channel(CNuc *compound, JGroup *jgroup, int channelIndex);

/*! Level and channel together, for a width. */
std::string LevelAndChannel(CNuc *compound, JGroup *jgroup, ALevel *level,
                            int jgroupIndex, int levelIndex, int channelIndex);

/*!
 * Describe the fit parameter at \p minuitIndex (0-based, in the order
 * CNuc::FillMnParams then EData::FillMnParams add them).  \p data may be null,
 * in which case only the R-matrix parameters are described.
 *
 * Returns an empty string if the index is out of range, so callers can fall
 * back to the raw Minuit name.
 */
std::string Parameter(CNuc *compound, EData *data, int minuitIndex);

/*!
 * Total number of parameters the walk covers; matches the size of the Minuit
 * parameter list built from the same compound nucleus and data.
 */
int NumParameters(CNuc *compound, EData *data);

}  // namespace AZURELabel

#endif  // PARAMETERLABEL_H
