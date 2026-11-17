#ifndef LEVEL_REPORT_H
#define LEVEL_REPORT_H

class CompoundSystem;
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
void WriteCompoundNucleusCheck(const CompoundSystem &compound, const LevelInputs &inputs, const Config &configure);
/// checks/pathways.chk: the internal and external reaction pathways of the pathway table.
void WritePathwaysCheck(const CompoundSystem &compound, const Config &configure);
/// checks/boundaryconditions.chk: B_c of every channel of the R-matrix groups.
void WriteBoundaryConditionsCheck(const CompoundSystem &compound, const Config &configure);
/// checks/angulardistributions.chk: the interference terms and their coefficients, from the pathway table.
void WriteAngularDistributionsCheck(const CompoundSystem &compound, const Config &configure);

/// An observed width as parameters.out prints it: the value in keV, eV or meV, the largest unit in which it prints
/// (six decimals) as at least 1.000000.
struct ReportedWidth {
  double value;
  const char *unit;
};
ReportedWidth ReportWidth(double W_MeV);

/// parameters.out: every level's Brune (Park) energy and, per channel, the observed width, ANC, moment or B with the
/// internal and external amplitudes (the table CompoundSystem::TransformOut returns); relativistic.out too with relativistic=1.
void WriteParametersOut(const CompoundSystem &compound, const LevelState &state, const ObservedLevelTable &table,
                        const Config &configure);
/// relativistic.out: the parameters.out levels in the E(s) convention of the LANL code EDA (ENDF KRL=1).
void WriteRelativisticOut(const CompoundSystem &compound, const ObservedLevelTable &table, const Config &configure);

#endif
