#ifndef CNUCEXTERNALCAPTURE_H
#define CNUCEXTERNALCAPTURE_H

#include <utility>
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

/// The ERROR of a photon pair (number `pairNum`, 1-based) whose final state, level `level` of its J^pi group, is the
/// final state of the earlier pair `otherPairKey` (experiment::ExternalCaptureFinalStates decides it): a photon pair
/// has its own final state (else the second pair took the level's capture parameters, and capture into the first ran
/// without its external part).
void FinalStateOfAnotherPair(CNuc &compound, const Config &configure, int otherPairKey, int level, int pairNum);

/*!
 * N_f^(-1/2) of level `level` of `group` (Eq. 1): n0 = 1 (Brune) or Park's J, gamma the level's amplitude per group
 * channel (1-based channel ch at gamma[ch - 1]); the channels with an N_f integral (ch <= NumNFIntegrals) contribute.
 * 1 / sqrt(N_f) as the engine has always formed it (no check of N_f > 0: a Park norm at or below 0 is reported where
 * it is computed).
 */
double SqrtNFFactor(CNuc &compound, JGroup *group, const ALevel *level, double n0, const std::vector<double> &gamma);

/*!
 * The final state of external capture through the photon channel `channel`: the first external-capture level whose
 * pair is the channel's, when its multipolarity flag (column 31) covers the channel's E1, M1 or E2; (group, level),
 * 1-based, or (0, 0) when the channel has no external part.
 */
std::pair<int, int> ExternalCaptureFinalState(CNuc &compound, AChannel *channel);

/// The numbers gamma^ext of one photon channel depends on (Eqs. 2 and 4), as values.
struct ExternalCaptureInputs {
  double levelEnergy = 0.0;        ///< E_lambda (MeV, E_x)
  std::vector<double> levelGamma;  ///< the level's amplitudes (MeV^1/2) per channel of its group, gamma[ch - 1]
  double finalEnergy = 0.0;        ///< E_f (MeV, E_x)
  std::vector<double> finalGamma;  ///< the final state's amplitudes (MeV^1/2) per channel of its group
  double finalSqrtNF = 1.0;        ///< N_f^(-1/2) of the final state
};

/// The partial derivatives of gamma^ext, exact since Eq. (4) is bilinear: d gamma^ext / d gamma_lambda c per channel
/// of the level's group (0 for a channel that is no particle channel) and d gamma^ext / d s_fc per channel of the
/// final state's group (MeV^1/2 per unit of s_fc).
struct ExternalCaptureSlopes {
  std::vector<complex> level;
  std::vector<complex> finalScale;
};

/*!
 * gamma^ext (MeV^1/2, complex) of the photon channel `channel` of a level of `group` into `finalState`
 * (ExternalCaptureFinalState) at the values `in` (Eq. 4, with s_fc = gamma~_fc N_f^(-1/2) K_c, Eq. 2): the radial
 * integrals are ECIntegral's direct ones with the run's per-thread scratch.  With `slopes`, also its derivatives,
 * from the same integrals.  Zero unless E1, E2 or M1 into a final state.
 */
complex ExternalReducedWidthAt(CNuc &compound, JGroup *group, AChannel *channel, std::pair<int, int> finalState,
                               const ExternalCaptureInputs &in, const Config &configure,
                               ExternalCaptureSlopes *slopes = nullptr);

/*!
 * gamma^ext of the photon channel `channel` of level `level` in `group` (ExternalReducedWidthAt): the level's energy
 * and amplitudes are its input values (initial) or its transform slots (Brune values of the fit parameters), the
 * final state's likewise, with its stored N_f^(-1/2).
 */
complex ExternalReducedWidth(CNuc &compound, JGroup *group, ALevel *level, AChannel *channel, bool initial,
                             const Config &configure);

#endif
