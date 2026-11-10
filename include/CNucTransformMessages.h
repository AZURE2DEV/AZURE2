#ifndef CNUC_TRANSFORM_MESSAGES_H
#define CNUC_TRANSFORM_MESSAGES_H

#include <cstddef>

#include "azure2/physics/LevelParameters.h"

class CNuc;
class Config;

/*!
 * The messages of the level-parameter transformations (CNuc::TransformIn, CNuc::TransformToBrune): what physics/
 * level_parameters reports about a level, worded for the user with the level's and channel's labels.  Each function
 * reads the nucleus through const access and writes to the configuration's output stream only.  j, la and ch are the
 * engine's 1-based J-group, level and channel numbers; lb and nc index the block arrays of the physics/ result.
 */

/// An input column that level_parameters could not convert as given (amplitude flag on a non-particle channel, a
/// width where P = 0, an ANC exactly at threshold); `input` is the file's value.
void WarnInputColumn(const CNuc &compound, const Config &configure, int j, int la, int ch,
                     azure2::physics::ColumnNote note, double input);
/// Observed widths that no reduced width reproduces: Park's J not positive (Park basis), or the Brune denominator
/// 1 - sum_c Gamma_c (dS_c/dE)/(2 P_c) negative, with its largest term and the likely cause.  Silent otherwise.
void WarnObservedOutOfReach(const CNuc &compound, const Config &configure, int j, int la,
                            const azure2::physics::BasisAmplitudes &amp, const azure2::physics::ProjectLevels &in,
                            std::size_t lb, std::size_t nc);
/// The imaginary part of the external-capture amplitude exceeds the total photon amplitude of (j, la, ch).
void WarnExternalExceedsTotal(const CNuc &compound, const Config &configure, int j, int la, int ch);
/// Park's J of level (j, la), as the level holds it, is not positive.
void WarnParkNormNotPositive(const CNuc &compound, const Config &configure, int j, int la);
/// A level whose formal -> Brune conversion failed, was bracketed or did not converge (brune_from_formal's status).
void WarnFormalToBrune(const CNuc &compound, const Config &configure, int j, int la,
                       const azure2::physics::LevelSet &out, std::size_t lb, std::size_t nc,
                       const azure2::physics::FormalToBruneOptions &options);

#endif
