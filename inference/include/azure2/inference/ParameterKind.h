#ifndef AZURE2_INFERENCE_PARAMETERKIND_H
#define AZURE2_INFERENCE_PARAMETERKIND_H

#include <string>
#include <utility>
#include <vector>

/*!
 * What a fit parameter is, and what it belongs to.
 *
 * Every entry of the parameter vector has one ParameterRecord, written where the
 * parameter is created (CompoundSystem::FillMnParams, WorkingCopy::FillMnParams, through
 * AZUREParams::Add) and held beside the Minuit parameters in AZUREParams.  It is
 * the one answer to "is this a normalization?", "which segment?", "does it carry
 * its own prior?": no code asks the parameter's name, which is a label for
 * people and files (docs/upgrade/parameter_kinds.md).
 *
 * Plain data and the standard library only: a part of the inference/ library
 * (docs/upgrade/inference.md), which describes the parameter vector and does not
 * hold its values (they are the engine's AZUREParams, on Minuit2's parameters).
 */

namespace azure2 {
namespace inference {

/// Kinds of fit parameter, in their order in the vector.  The values never change; a new kind takes the next value,
/// a row in ParamKindTable() and the fields that say what it belongs to in the API's record table
/// (api/src/pybind_records.cpp).
enum class ParamKind {
  LevelEnergy = 0,      ///< level energy E_lambda of a (J group, level)
  Gamma = 1,            ///< reduced-width amplitude gamma_lambda,c of a (J group, level, channel)
  Norm = 2,             ///< normalization of a data segment
  EnergyShift = 3,      ///< constant energy shift of a data segment
  EnergyShiftSqrt = 4,  ///< sqrt(E) energy-shift coefficient of a data segment
  EnergyShiftFn = 5,    ///< parameter a_j of a data segment's shift function ("shiftfn")
  Abundance = 6,        ///< isotopic abundance of a compound nucleus in a sample (<abundances>)
  Screening = 7,        ///< electron-screening potential U_e of a <targetInt> line
  Density = 8,          ///< fitted areal density of a <targetInt> line or density group
  Doppler = 9,          ///< fitted free-gas Doppler temperature of a <targetInt> line
  Thickness = 10,       ///< fitted sample thickness of a transmission / yield segment
  ThmCoherent = 11,     ///< Re or Im of a THM coherent-background coefficient (cbackground=)
  Apparatus = 12,       ///< a number of an apparatus record (ApparatusParameters.h: a stopping scale, a kernel's width)
  RadiativeDecay = 13,  ///< Gamma_gamma, a branching or a mixing ratio of a gamma_decay record, in a photon channel's slot
  IsospinShift = 14,    ///< internal Coulomb shift of the isospin levels of a compound nucleus (an isospin_shift record)
  IsospinMix = 15       ///< amplitude of an isospin-forbidden channel of an isospin level (an isospin_mix record)
};

/// A parameter in a level's slot of the vector: its energy, a channel's amplitude, or a number of a gamma_decay record
/// that the map of GammaDecayMap.h turns into a photon amplitude.  The R-matrix gradient covers them all.
inline bool IsLevelSlot(ParamKind kind) {
  return kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::RadiativeDecay;
}

/// A parameter the R-matrix adjoint differentiates: a level slot, or an isospin shift or mix, which derived level
/// slots follow (IsospinTies.h; the chain rule of GradAccum::Scatter).
inline bool IsRMatrixParameter(ParamKind kind) {
  return IsLevelSlot(kind) || kind == ParamKind::IsospinShift || kind == ParamKind::IsospinMix;
}

/// One row per kind: what a maintainer and the API need to know about it.
struct ParamKindRow {
  ParamKind kind;
  const char *key;    ///< the API's and pyazr's name of the kind (Parameter.kind)
  const char *unit;   ///< unit of the parameter's value in the vector
  const char *owner;  ///< which ParameterRecord fields say what it belongs to
  /// The parameter carries its own prior (segment error, <abundances>, target-effect or transmission block): a
  /// <parameterSettings> nuisance flag adds nothing to it.
  bool ownPrior;
  /// It belongs to a data segment: a run without data (extrapolation, reaction rate) does not take it from a
  /// parameter file, whose segment k is a data segment, not the run's test grid k.
  bool segmentField;
  /// Whether it varies is its project line's vary flag: a parameter file's zero or nonzero value does not fix or
  /// free it (a level parameter and a normalization are freed by a nonzero value and fixed by a zero).
  bool varyFromProject;
  const char *doc;
};

/// The kinds, in the order of their codes.  A constant table, built once.
inline const std::vector<ParamKindRow> &ParamKindTable() {
  static const std::vector<ParamKindRow> table = {
      {ParamKind::LevelEnergy, "energy", "MeV", "jGroup, level, levelNumber", false, false, false,
       "level energy E_lambda (excitation energy in the compound nucleus) in the run's basis"},
      {ParamKind::Gamma, "width", "MeV^1/2", "jGroup, level, channel, levelNumber", false, false, false,
       "reduced-width amplitude gamma_lambda,c in the run's basis (the .azr gives a partial width or an ANC)"},
      {ParamKind::Norm, "norm", "1", "segment, segmentKey", true, true, false,
       "normalization factor of a data segment (a profiled THM norm is not a parameter)"},
      {ParamKind::EnergyShift, "shift", "MeV", "segment, segmentKey", true, true, true,
       "constant lab-energy shift Delta E of a data segment"},
      {ParamKind::EnergyShiftSqrt, "shift_sqrt", "MeV^1/2", "segment, segmentKey", true, true, true,
       "coefficient b of the lab-energy shift b sqrt(E/MeV) of a data segment"},
      {ParamKind::EnergyShiftFn, "shift_fn", "as the shift function defines", "segment, segmentKey, fnTerm", true, true,
       true, "parameter a_j of a data segment's shift function (a shift in MeV lab)"},
      {ParamKind::Abundance, "abundance", "1", "sample, sampleName, nucleus, nucleusZ, nucleusA", true, false, false,
       "isotopic abundance of a compound nucleus in a sample of the <abundances> block"},
      {ParamKind::Screening, "screening", "eV", "effect", true, false, true,
       "electron-screening potential U_e of a <targetInt> line with a screening block"},
      {ParamKind::Density, "density", "atoms/cm^2", "effect", true, false, true,
       "fitted areal target density of a <targetInt> line (densityfit) or density group (its first line)"},
      {ParamKind::Doppler, "doppler", "K", "effect", true, false, true,
       "fitted free-gas Doppler temperature of a <targetInt> line ('doppler T vary sigma')"},
      {ParamKind::Thickness, "thickness", "atoms/barn", "segment, segmentKey", true, true, true,
       "fitted sample thickness of a transmission or yield segment ('transmission n vary sigma')"},
      {ParamKind::ThmCoherent, "cbkg", "HOES amplitude units (c1: per MeV)",
       "thmExperiment, jGroup, channel (entrance), exitChannel, thmTerm", false, false, false,
       "Re or Im of c0 or c1 of a THM coherent background c(E) = c0 + c1 E (cbackground=)"},
      {ParamKind::Apparatus, "apparatus", "the row's (apparatusUnit)", "effect, apparatus, owner", true, false, true,
       "a number of an apparatus record a fit may vary, one per sharing key (docs/upgrade/apparatus_parameters.md): "
       "the stopping scale of a (material, projectile), a table kernel's parameter (format 3: written as a nuisance)"},
      {ParamKind::RadiativeDecay, "radiative", "eV (total_width), 1 (branching, mixing)",
       "jGroup, level, channel, levelNumber, decayRecord, decayPart, finalState", true, false, true,
       "a number of a level's gamma_decay record (docs/upgrade/gamma_parametrisation.md): the total radiative width, the "
       "branching to a final state, or its mixing ratio (Krane-Steffen), at the position of one of the level's photon "
       "amplitudes"},
      {ParamKind::IsospinShift, "isospin_shift", "MeV", "isospinRecord, nucleus, nucleusZ, nucleusA", true, false, true,
       "internal Coulomb energy shift of the isospin levels of a compound nucleus (an isospin_shift record of <isospin>, "
       "named <record>.shift): "
       "the members of the nucleus sit at the reference member's energy above threshold plus this shift"},
      {ParamKind::IsospinMix, "isospin_mix", "MeV^1/2", "isospinRecord, jGroup, level, channel", true, false, true,
       "amplitude of an isospin-forbidden channel of an isospin level (an isospin_mix record of <isospin>, named "
       "<record>.amplitude), which the "
       "channel's slot of the level follows"},
  };
  return table;
}

/// The row of a kind.
inline const ParamKindRow &KindRow(ParamKind kind) { return ParamKindTable()[(int)kind]; }

/*!
 * What one parameter is and what it belongs to.  Fields that do not apply to
 * its kind are -1 (or empty); the owner column of its ParamKindTable() row lists
 * the ones that do.  Indices are 1-based, as the engine's.
 */
struct ParameterRecord {
  ParamKind kind = ParamKind::LevelEnergy;
  int ordinal = 0;  ///< place among the parameters of its kind, from 1 (set by AZUREParams::Add)

  int jGroup = -1;       ///< J group (LevelEnergy, Gamma, ThmCoherent)
  int level = -1;        ///< level within the J group (LevelEnergy, Gamma)
  int channel = -1;      ///< channel within the J group (Gamma; ThmCoherent: the entrance channel)
  int levelNumber = -1;  ///< level counted over all J groups (LevelEnergy, Gamma), as the GUI and param files number it

  int segment = -1;     ///< data segment, position in the segment list (Norm, EnergyShift*, Thickness)
  int segmentKey = -1;  ///< that segment's key (its <segmentsData> line)
  int fnTerm = -1;      ///< j of a_j (EnergyShiftFn)

  int effect = -1;  ///< the <targetInt> line's effect number (Screening, Density, Doppler; Apparatus: its first)

  int sample = -1;          ///< sample of the <abundances> block (Abundance)
  std::string sampleName;   ///< its name, empty for the unnamed sample (Abundance)
  int nucleus = -1;         ///< compound nucleus of the block (Abundance)
  int nucleusZ = -1, nucleusA = -1;

  std::string apparatus;      ///< the row of ApparatusRows() (Apparatus): "stopping_scale", "kernel"
  std::string apparatusUnit;  ///< the value's unit (Apparatus)
  /// What it belongs to (Apparatus), in the order of its name, whose parts these values are joined by '.':
  /// (quantity, stopping_scale), (material, TiN), (projectile, 1H) for stopping_scale.TiN.1H; (record, det1),
  /// (field, energy_resolution),
  /// (parameter, width) for det1.energy_resolution.width.
  std::vector<std::pair<std::string, std::string>> owner;

  std::string thmExperiment;  ///< THM experiment (ThmCoherent)
  int exitChannel = -1;       ///< exit channel within the J group (ThmCoherent)
  int thmTerm = -1;           ///< 0 Re c0, 1 Im c0, 2 Re c1, 3 Im c1 (ThmCoherent)

  std::string decayRecord;  ///< the gamma_decay record (RadiativeDecay)
  std::string decayPart;    ///< total_width, branching or mixing (RadiativeDecay)
  std::string finalState;   ///< the final state's photon pair as the record writes it (branching, mixing)

  std::string isospinRecord;  ///< the isospin_shift or isospin_mix record (IsospinShift, IsospinMix)

  const ParamKindRow &Row() const { return KindRow(kind); }

  static ParameterRecord LevelEnergy(int jGroup, int level, int levelNumber) {
    ParameterRecord r;
    r.kind = ParamKind::LevelEnergy;
    r.jGroup = jGroup;
    r.level = level;
    r.levelNumber = levelNumber;
    return r;
  }
  static ParameterRecord Width(int jGroup, int level, int channel, int levelNumber) {
    ParameterRecord r = LevelEnergy(jGroup, level, levelNumber);
    r.kind = ParamKind::Gamma;
    r.channel = channel;
    return r;
  }
  /// A number of a gamma_decay record in the slot of channel `channel` of (jGroup, level).
  static ParameterRecord Radiative(int jGroup, int level, int channel, int levelNumber, const std::string &record,
                                   const std::string &part, const std::string &finalState) {
    ParameterRecord r = Width(jGroup, level, channel, levelNumber);
    r.kind = ParamKind::RadiativeDecay;
    r.decayRecord = record;
    r.decayPart = part;
    r.finalState = finalState;
    return r;
  }
  /// The shift of an isospin_shift record: compound nucleus `nucleus` (1-based) with its Z and A.
  static ParameterRecord IsospinShift(const std::string &record, int nucleus, int z, int a) {
    ParameterRecord r;
    r.kind = ParamKind::IsospinShift;
    r.isospinRecord = record;
    r.nucleus = nucleus;
    r.nucleusZ = z;
    r.nucleusA = a;
    return r;
  }
  /// The amplitude of an isospin_mix record, which channel `channel` of (jGroup, level) follows.
  static ParameterRecord IsospinMix(const std::string &record, int jGroup, int level, int channel) {
    ParameterRecord r;
    r.kind = ParamKind::IsospinMix;
    r.isospinRecord = record;
    r.jGroup = jGroup;
    r.level = level;
    r.channel = channel;
    return r;
  }
  /// A parameter of data segment `segment` (Norm, EnergyShift, EnergyShiftSqrt, Thickness).
  static ParameterRecord OfSegment(ParamKind kind, int segment, int segmentKey) {
    ParameterRecord r;
    r.kind = kind;
    r.segment = segment;
    r.segmentKey = segmentKey;
    return r;
  }
  static ParameterRecord ShiftFnTerm(int segment, int segmentKey, int j) {
    ParameterRecord r = OfSegment(ParamKind::EnergyShiftFn, segment, segmentKey);
    r.fnTerm = j;
    return r;
  }
  /// A parameter of target effect `effect` (Screening, Density, Doppler).
  static ParameterRecord OfEffect(ParamKind kind, int effect) {
    ParameterRecord r;
    r.kind = kind;
    r.effect = effect;
    return r;
  }
  /// A parameter of an apparatus record (ApparatusParameters.h); \p effect its first target effect.
  static ParameterRecord Apparatus(int effect, const std::string &row, const std::string &unit,
                                   const std::vector<std::pair<std::string, std::string>> &owner) {
    ParameterRecord r = OfEffect(ParamKind::Apparatus, effect);
    r.apparatus = row;
    r.apparatusUnit = unit;
    r.owner = owner;
    return r;
  }
  static ParameterRecord Abundance(int sample, const std::string &sampleName, int nucleus, int z, int a) {
    ParameterRecord r;
    r.kind = ParamKind::Abundance;
    r.sample = sample;
    r.sampleName = sampleName;
    r.nucleus = nucleus;
    r.nucleusZ = z;
    r.nucleusA = a;
    return r;
  }
  static ParameterRecord ThmCoherent(const std::string &experiment, int jGroup, int entrance, int exit, int term) {
    ParameterRecord r;
    r.kind = ParamKind::ThmCoherent;
    r.thmExperiment = experiment;
    r.jGroup = jGroup;
    r.channel = entrance;
    r.exitChannel = exit;
    r.thmTerm = term;
    return r;
  }
};

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_PARAMETERKIND_H
