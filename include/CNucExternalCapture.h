#ifndef CNUCEXTERNALCAPTURE_H
#define CNUCEXTERNALCAPTURE_H

#include <vector>

#include "Constants.h"

class AChannel;
class ALevel;
class CNuc;
class Config;
class JGroup;

/*!
 * The engine's side of physics/ external_capture (docs/upgrade/external_capture.md section 2): it reads a level, its
 * group and the final state from CNuc and calls the physics functions.  Nothing else of the engine computes N_f or
 * gamma^ext.
 */

/// Is level `level` of J-group `group` already the final state of another photon pair than `pairNum` (pair numbers
/// 1-based)?  Then an ERROR naming both pair keys and the level is written: a photon pair has its own final state (a
/// second pair took the level's capture parameters, and capture into the first ran without its external part).
bool FinalStateOfAnotherPair(CNuc &compound, const Config &configure, int group, int level, int pairNum);

/*!
 * N_f^(-1/2) of level `level` of `group` (Eq. 1): n0 = 1 (Brune) or Park's J, gamma the level's amplitude per group
 * channel (1-based channel ch at gamma[ch - 1]); the channels with an N_f integral (ch <= NumNFIntegrals) contribute.
 * 1 / sqrt(N_f) as the engine has always formed it (no check of N_f > 0: a Park norm at or below 0 is reported where
 * it is computed).
 */
double SqrtNFFactor(CNuc &compound, JGroup *group, const ALevel *level, double n0, const std::vector<double> &gamma);

/*!
 * gamma^ext (MeV^1/2, complex) of the photon channel `channel` of level `level` in `group` (Eq. 4): the final state is
 * the first external-capture level whose pair is the channel's; the level's energy and amplitudes are its input
 * values (initial) or its transform slots (Brune values of the fit parameters), the final state's likewise, with its
 * stored N_f^(-1/2) and conversion factor K_c.  Zero unless E1, E2 or M1 with the final state's multipolarity flag.
 */
complex ExternalReducedWidth(CNuc &compound, JGroup *group, ALevel *level, AChannel *channel, bool initial,
                             const Config &configure);

#endif
