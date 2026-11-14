#ifndef CNUC_REPORT_H
#define CNUC_REPORT_H

class CNuc;
class LevelInputs;
class LevelState;
class ObservedLevelTable;
class Config;

/*!
 * What the compound nucleus writes: the check files (checks/*.chk, or the screen, by the check masks of the
 * configuration) and the level parameters after a calculation (parameters.out, relativistic.out).  Each function
 * reads the nucleus through const access and writes text only; the values are those the engine has set (TransformOut
 * for the observed parameters).
 */

/// checks/compoundnucleus.chk: pairs, J-groups, levels and channels as read (before Initialize).
void WriteCompoundNucleusCheck(const CNuc &compound, const LevelInputs &inputs, const Config &configure);
/// checks/pathways.chk: the internal and external reaction pathways of the pathway table.
void WritePathwaysCheck(const CNuc &compound, const Config &configure);
/// checks/boundaryconditions.chk: B_c of every channel of the R-matrix groups.
void WriteBoundaryConditionsCheck(const CNuc &compound, const Config &configure);
/// checks/angulardistributions.chk: the interference terms and their coefficients, from the pathway table.
void WriteAngularDistributionsCheck(const CNuc &compound, const Config &configure);

/// parameters.out: every level's Brune (Park) energy and, per channel, the observed width, ANC, moment or B with the
/// internal and external amplitudes (the table CNuc::TransformOut returns); relativistic.out too with relativistic=1.
void WriteParametersOut(const CNuc &compound, const LevelState &state, const ObservedLevelTable &table,
                        const Config &configure);
/// relativistic.out: the parameters.out levels in the E(s) convention of the LANL code EDA (ENDF KRL=1).
void WriteRelativisticOut(const CNuc &compound, const ObservedLevelTable &table, const Config &configure);

#endif
