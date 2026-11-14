// The records pyazr builds its parameter table and its pairs from (pyazr/parameters.py): one dict per fit parameter
// (Session.parameter_records) and one per particle pair (Session.pair_records), with named fields.  Each field is one
// row of a table below: its name, its unit, which kinds of parameter it describes, whether the parameter table shows
// it as the parameter's owner, a line of doc, and the function that reads it from the engine.  The names are typed
// here only: pyazr reads them (parameter_record_fields, pair_record_fields; pyazr.kinds()), and a new field or a new
// kind of parameter is a row here (tests/pyazr/parameter_records_test.py fails for a kind no row describes).
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/LevelEdits.h"
#include "azure2/io/Document.h"
#include "AChannel.h"
#include "ALevel.h"
#include "AZUREAPI.h"
#include "ApparatusParameters.h"
#include "CNuc.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
#include "ParameterKind.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "ThmExperiment.h"
#include "WorkingCopy.h"

namespace py = pybind11;

namespace {

// A set of parameter kinds, one bit per ParamKind.
using KindSet = unsigned;
constexpr KindSet kEveryKind = ~0u;
KindSet Kinds(std::initializer_list<ParamKind> kinds) {
  KindSet set = 0;
  for (ParamKind k : kinds) set |= 1u << (unsigned)k;
  return set;
}
bool InSet(KindSet set, ParamKind kind) { return (set >> (unsigned)kind) & 1u; }

// What a field's reader sees of one parameter: the session, the parameter's place and record, and the engine objects
// it belongs to, looked up once (null where they do not apply).
struct ParameterView {
  const AZUREAPI *api;
  int index;  // among all parameters
  const ParameterRecord *record;
  CNuc *nuc;
  const DataSet *data;
  const JGroup *group = nullptr;          // level kinds, cbkg
  const ALevel *level = nullptr;          // level kinds
  const AChannel *channel = nullptr;      // width, radiative; cbkg: the entrance channel
  const AChannel *exitChannel = nullptr;  // cbkg
  const Segment *segment = nullptr;          // segment kinds
  const TargetEffect *effect = nullptr;      // screening, density, doppler
  const ApparatusParam *apparatus = nullptr;  // apparatus
  const SampleAbundance *abundance = nullptr;  // abundance
  int freeIndex = -1;                          // position in the fit vector, -1 when fixed
  const vector_r *start;                       // every parameter's start value, fit basis (params_all_rwa)
  const vector_r *value;                       // and as the project states it (params_all)
  const std::map<std::pair<int, int>, std::string> *levelLabels;
};

// The declared prior of a parameter, (centre, width) in its unit, the width 0 when it has none: the numbers the
// objective's diagonal priors are made of, read where they read them (include/Objective.h, ForEachDiagonalPrior;
// docs/upgrade/parameter_kinds.md section 4), here whether or not the parameter is free.
std::pair<double, double> PriorOf(const ParameterView &v) {
  const ParameterRecord &r = *v.record;
  switch (r.kind) {
    case ParamKind::Norm: {
      const double n0 = v.segment->GetNominalNorm();
      return {n0, n0 / 100.0 * v.segment->GetNormError()};
    }
    case ParamKind::EnergyShift: return {v.segment->GetNominalEnergyShift(), v.segment->GetEnergyShiftError()};
    case ParamKind::EnergyShiftSqrt:
      return {v.segment->GetEnergyShiftSqrtCentre(), v.segment->GetEnergyShiftSqrtError()};
    case ParamKind::EnergyShiftFn:
      return {v.segment->GetShiftFnCentres()[r.fnTerm], v.segment->GetShiftFnErrors()[r.fnTerm]};
    case ParamKind::Abundance: return {v.abundance->nominal, v.abundance->sigma};
    case ParamKind::Screening: return {v.effect->GetScreeningNominal(), v.effect->GetScreeningSigma()};
    case ParamKind::Density: return {v.effect->GetDensityCentre(), v.effect->GetDensitySigma()};
    case ParamKind::Doppler: return {v.effect->GetDopplerCentre(), v.effect->GetDopplerSigma()};
    case ParamKind::Thickness:
      return {v.segment->GetThicknessCentre(),
              v.segment->IsThicknessGroupFollower() ? 0.0 : v.segment->GetThicknessSigma()};
    case ParamKind::Apparatus: {
      const ApparatusDeclaration &d = v.apparatus->declaration;
      return {d.centre, d.vary ? d.sigma : 0.0};
    }
    default: return {0.0, 0.0};
  }
}

// The largest value the sub-point grids of a fitted density, Doppler temperature or apparatus number are laid out for.
double LayoutMaxOf(const ParameterView &v) {
  switch (v.record->kind) {
    case ParamKind::Density: return v.effect->GetLayoutDensity();
    case ParamKind::Doppler: return v.effect->GetDopplerLayoutTemperature();
    default: return v.apparatus->declaration.layout;
  }
}

py::object OrNone(int i) { return i < 0 ? py::object(py::none()) : py::object(py::int_(i)); }

struct ParameterField {
  const char *name;
  const char *unit;  // "unit": the parameter's own unit (its `unit` field); "" for a number without a unit or a text
  KindSet kinds;     // the kinds it describes
  bool owner;        // shown as what the parameter belongs to (Parameters.table(), pyazr.kinds()'s owner fields)
  const char *doc;
  py::object (*read)(const ParameterView &);
};

const KindSet kLevel = Kinds({ParamKind::LevelEnergy, ParamKind::Gamma, ParamKind::RadiativeDecay});
const KindSet kChannel = Kinds({ParamKind::Gamma, ParamKind::RadiativeDecay});
const KindSet kSegment = Kinds({ParamKind::Norm, ParamKind::EnergyShift, ParamKind::EnergyShiftSqrt,
                                ParamKind::EnergyShiftFn, ParamKind::Thickness});
const KindSet kPrior = Kinds({ParamKind::Norm, ParamKind::EnergyShift, ParamKind::EnergyShiftSqrt,
                              ParamKind::EnergyShiftFn, ParamKind::Abundance, ParamKind::Screening, ParamKind::Density,
                              ParamKind::Doppler, ParamKind::Thickness, ParamKind::Apparatus});

const std::vector<ParameterField> &ParameterFields() {
  static const std::vector<ParameterField> fields = {
      // -- every parameter
      {"index", "", kEveryKind, false, "position among all parameters (fixed ones included; params_all_rwa order)",
       [](const ParameterView &v) -> py::object { return py::int_(v.index); }},
      {"name", "", kEveryKind, false, "the parameter's name, as param.sav and the project's records give it",
       [](const ParameterView &v) -> py::object { return py::str(v.api->params_names(v.index)); }},
      {"kind", "", kEveryKind, false, "the key of its kind (pyazr.kinds()['parameters'])",
       [](const ParameterView &v) -> py::object { return py::str(v.record->Row().key); }},
      {"ordinal", "", kEveryKind, false, "its place among the parameters of its kind, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->ordinal); }},
      {"unit", "", kEveryKind, false, "the unit of start, bounds and the prior: the kind's, an apparatus row's own",
       [](const ParameterView &v) -> py::object {
         return py::str(v.record->kind == ParamKind::Apparatus ? v.record->apparatusUnit : v.record->Row().unit);
       }},
      {"start", "unit", kEveryKind, false, "the start value in the fit vector's basis (its entry of x0 when free)",
       [](const ParameterView &v) -> py::object { return py::float_((*v.start)[v.index]); }},
      {"value", "", kEveryKind, false,
       "the start value as the project states it: a level energy in MeV, a width as a partial width in eV or an ANC "
       "in fm^-1/2 (an amplitude in MeV^1/2 when the project gives one); the other kinds as start",
       [](const ParameterView &v) -> py::object { return py::float_((*v.value)[v.index]); }},
      {"fixed", "", kEveryKind, false, "held fixed: not an entry of the fit vector",
       [](const ParameterView &v) -> py::object { return py::bool_(v.freeIndex < 0); }},
      {"free_index", "", kEveryKind, false, "position in the fit vector x, None when fixed",
       [](const ParameterView &v) -> py::object { return OrNone(v.freeIndex); }},
      {"bounds", "unit", kEveryKind, false,
       "(lower, upper) limits of the fit, -inf / inf where it has none: those the parameter is created with, then "
       "those of a <parameterSettings> block; outside them Objective.log_prob is -inf, as the program's sampler has it",
       [](const ParameterView &v) -> py::object {
         const std::pair<double, double> b = v.api->GetParameterBounds(v.index);
         return py::make_tuple(b.first, b.second);
       }},
      // -- a level's parameters
      {"jgroup", "", kLevel | Kinds({ParamKind::ThmCoherent}), false, "the J group, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->jGroup); }},
      {"J", "", kLevel | Kinds({ParamKind::ThmCoherent}), false, "the J group's total spin",
       [](const ParameterView &v) -> py::object { return py::float_(v.group->GetJ()); }},
      {"parity", "", kLevel | Kinds({ParamKind::ThmCoherent}), false, "the J group's parity, +1 or -1",
       [](const ParameterView &v) -> py::object { return py::int_(v.group->GetPi()); }},
      {"jpi", "", kLevel | Kinds({ParamKind::ThmCoherent}), false,
       "the J group's J^pi as the engine writes it (experiment::JpiText): \"1/2+\", \"2-\"",
       [](const ParameterView &v) -> py::object {
         return py::str(azure2::experiment::JpiText((int)std::lround(2.0 * v.group->GetJ()), v.group->GetPi()));
       }},
      {"level", "", kLevel, false, "the level within its J group, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->level); }},
      {"level_number", "", kLevel, false, "the level counted over all J groups, from 1 (param files, energy_<n>)",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->levelNumber); }},
      {"level_energy", "MeV", kLevel, false, "the level's energy E_lambda as the session holds it (excitation energy)",
       [](const ParameterView &v) -> py::object {
         return py::float_(v.nuc->Inputs().E(v.record->jGroup, v.record->level));
       }},
      {"level_label", "", kLevel, true,
       "the level's label in the project, '<J^pi> <E_x as written>' (Project.levels())",
       [](const ParameterView &v) -> py::object {
         auto it = v.levelLabels->find({v.record->jGroup, v.record->level});
         return it == v.levelLabels->end() ? py::object(py::none()) : py::object(py::str(it->second));
       }},
      // -- a channel's amplitude (and the THM background's entrance channel)
      {"channel", "", kChannel | Kinds({ParamKind::ThmCoherent}), false,
       "the channel within the J group, from 1 (cbkg: the entrance channel)",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->channel); }},
      {"L", "", kChannel | Kinds({ParamKind::ThmCoherent}), true,
       "the channel's orbital angular momentum l (a photon channel: its multipolarity)",
       [](const ParameterView &v) -> py::object { return py::int_(v.channel->GetL()); }},
      {"S", "", kChannel | Kinds({ParamKind::ThmCoherent}), true, "the channel spin s",
       [](const ParameterView &v) -> py::object { return py::float_(v.channel->GetS()); }},
      {"pair", "", kChannel | Kinds({ParamKind::ThmCoherent}), false,
       "the channel's pair number (Pair.number; cbkg: the exit channel's)",
       [](const ParameterView &v) -> py::object {
         return py::int_((v.exitChannel ? v.exitChannel : v.channel)->GetPairNum());
       }},
      {"pair_key", "", kChannel, true, "the channel's pair key in <levels> (Pair.key)",
       [](const ParameterView &v) -> py::object {
         return py::int_(v.nuc->GetPair(v.channel->GetPairNum())->GetPairKey());
       }},
      {"radiation_type", "", kChannel, true, "'P' particle, 'E' or 'M' photon, 'F' or 'G' beta decay",
       [](const ParameterView &v) -> py::object { return py::str(std::string(1, v.channel->GetRadType())); }},
      {"wigner_limit", "MeV", kChannel, false,
       "the Wigner limit gamma^2_W = 3 hbar^2 / (2 mu a^2) of the channel's reduced width (not of the amplitude)",
       [](const ParameterView &v) -> py::object { return py::float_(v.channel->GetWignerLimit()); }},
      {"width_is_amplitude", "", kChannel, false,
       "the project gives the width as a reduced-width amplitude in MeV^1/2, not as a partial width or ANC",
       [](const ParameterView &v) -> py::object { return py::bool_(v.level->GammaIsRWA(v.record->channel)); }},
      // -- a number of a gamma_decay record
      {"decay_record", "", Kinds({ParamKind::RadiativeDecay}), true, "the gamma_decay record",
       [](const ParameterView &v) -> py::object { return py::str(v.record->decayRecord); }},
      {"decay_part", "", Kinds({ParamKind::RadiativeDecay}), true, "total_width (eV), branching or mixing",
       [](const ParameterView &v) -> py::object { return py::str(v.record->decayPart); }},
      {"final_state", "", Kinds({ParamKind::RadiativeDecay}), true,
       "the final state's photon pair as the record writes it, None for the total width",
       [](const ParameterView &v) -> py::object {
         return v.record->finalState.empty() ? py::object(py::none()) : py::object(py::str(v.record->finalState));
       }},
      // -- a data segment's parameters
      {"segment", "", kSegment, false, "the data segment, from 1 (the engine's order)",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->segment); }},
      {"segment_key", "", kSegment, true, "the segment's key: the data set's number among the data sets, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->segmentKey); }},
      {"fn_term", "", Kinds({ParamKind::EnergyShiftFn}), true, "j of the shift function's parameter a_j, from 0",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->fnTerm); }},
      {"sample_kind", "", Kinds({ParamKind::Thickness}), true, "'transmission' or 'yield'",
       [](const ParameterView &v) -> py::object { return py::str(v.segment->IsYield() ? "yield" : "transmission"); }},
      {"thickness_leader", "", Kinds({ParamKind::Thickness}), false,
       "the segment key of the thickness group's leader for a follower (whose value is the leader's), else None",
       [](const ParameterView &v) -> py::object {
         if (!v.segment->IsThicknessGroupFollower()) return py::none();
         return py::int_(v.data->GetSegments()[v.segment->GetThicknessLeaderPos()].GetSegmentKey());
       }},
      // -- a target effect's parameters
      {"effect", "", Kinds({ParamKind::Screening, ParamKind::Density, ParamKind::Doppler, ParamKind::Apparatus}),
       true, "the target effect, from 1 (an apparatus number: its first)",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->effect); }},
      {"apparatus", "", Kinds({ParamKind::Apparatus}), true, "the apparatus row: stopping_scale, kernel",
       [](const ParameterView &v) -> py::object { return py::str(v.record->apparatus); }},
      {"owner", "", Kinds({ParamKind::Apparatus}), false,
       "what the apparatus number belongs to, (what, value) pairs in the order of its name",
       [](const ParameterView &v) -> py::object { return py::cast(v.record->owner); }},
      // -- an isotopic abundance
      {"sample_index", "", Kinds({ParamKind::Abundance}), false, "the sample of the <abundances> block, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->sample); }},
      {"sample_name", "", Kinds({ParamKind::Abundance}), true, "the sample's name, '' for the unnamed sample",
       [](const ParameterView &v) -> py::object { return py::str(v.record->sampleName); }},
      {"nucleus_index", "", Kinds({ParamKind::Abundance}), false, "the compound nucleus of the block, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->nucleus); }},
      {"nucleus_Z", "", Kinds({ParamKind::Abundance}), true, "the compound nucleus's charge number",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->nucleusZ); }},
      {"nucleus_A", "", Kinds({ParamKind::Abundance}), true, "the compound nucleus's mass number",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->nucleusA); }},
      // -- a THM coherent background
      {"thm_experiment", "", Kinds({ParamKind::ThmCoherent}), true, "the THM experiment",
       [](const ParameterView &v) -> py::object { return py::str(v.record->thmExperiment); }},
      {"exit_channel", "", Kinds({ParamKind::ThmCoherent}), false, "the exit channel within the J group, from 1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->exitChannel); }},
      {"thm_jpi", "", Kinds({ParamKind::ThmCoherent}), true, "J^pi of the combination as cbackground= writes it",
       [](const ParameterView &v) -> py::object {
         return py::str(ThmSpinText(v.group->GetJ()) + (v.group->GetPi() > 0 ? "+" : "-"));
       }},
      {"thm_exit_key", "", Kinds({ParamKind::ThmCoherent}), true, "the exit pair's key",
       [](const ParameterView &v) -> py::object {
         return py::int_(v.nuc->GetPair(v.exitChannel->GetPairNum())->GetPairKey());
       }},
      {"thm_channels", "", Kinds({ParamKind::ThmCoherent}), true, "\"s,l,s',l'\" of the entrance and exit channels",
       [](const ParameterView &v) -> py::object {
         return py::str(ThmSpinText(v.channel->GetS()) + "," + std::to_string(v.channel->GetL()) + "," +
                        ThmSpinText(v.exitChannel->GetS()) + "," + std::to_string(v.exitChannel->GetL()));
       }},
      {"thm_term", "", Kinds({ParamKind::ThmCoherent}), true, "0 Re c0, 1 Im c0, 2 Re c1, 3 Im c1",
       [](const ParameterView &v) -> py::object { return py::int_(v.record->thmTerm); }},
      // -- a nuisance's prior and grids
      {"prior_centre", "unit", kPrior, false,
       "the centre of the prior the project declares (the penalty is ((p - centre) / prior_sigma)^2)",
       [](const ParameterView &v) -> py::object { return py::float_(PriorOf(v).first); }},
      {"prior_sigma", "unit", kPrior, false,
       "the width of that prior, 0 when it has none; paid only while the parameter is free",
       [](const ParameterView &v) -> py::object { return py::float_(PriorOf(v).second); }},
      {"layout_max", "unit", Kinds({ParamKind::Density, ParamKind::Doppler, ParamKind::Apparatus}), false,
       "the largest value the sub-point grids are laid out for: beyond it the engine uses this value and warns",
       [](const ParameterView &v) -> py::object { return py::float_(LayoutMaxOf(v)); }},
  };
  return fields;
}

// The engine objects a parameter belongs to, looked up from its record.  A record that names an object the session
// does not have is a bug of the engine's bookkeeping: refused, never read past.
void Resolve(ParameterView &v) {
  const ParameterRecord &r = *v.record;
  auto refuse = [&](const std::string &what) {
    throw std::logic_error("parameter_records: parameter " + std::to_string(v.index) + " (" +
                           v.api->params_names(v.index) + ", kind " + r.Row().key + ") " + what);
  };
  if (r.jGroup > 0) {
    if (r.jGroup > v.nuc->NumJGroups()) refuse("names a J group the session does not have");
    v.group = v.nuc->GetJGroup(r.jGroup);
    if (r.level > 0) v.level = v.group->GetLevel(r.level);
    if (r.channel > 0) v.channel = v.group->GetChannel(r.channel);
    if (r.exitChannel > 0) v.exitChannel = v.group->GetChannel(r.exitChannel);
  }
  if (r.segment > 0) v.segment = v.data->GetSegment(r.segment);
  switch (r.kind) {
    case ParamKind::Screening: v.effect = v.data->GetScreeningTargetEffect(r.ordinal); break;
    case ParamKind::Density: v.effect = v.data->GetDensityTargetEffect(r.ordinal); break;
    case ParamKind::Doppler: v.effect = v.data->GetDopplerTargetEffect(r.ordinal); break;
    case ParamKind::Apparatus: v.apparatus = &v.data->GetApparatusParam(r.ordinal); break;
    case ParamKind::Abundance: v.abundance = &v.data->Abundances().Sample(r.sample).entries[r.nucleus - 1]; break;
    case ParamKind::Thickness:
      if (v.data->GetThicknessSegment(r.ordinal) != r.segment) refuse("is not the thickness of its segment");
      break;
    default: break;
  }
}

// (J group, level) -> the level's label in the <levels> section of the session's project file
// (experiment::ReadLevelRowTable on that section alone): the names Project.levels() gives.  Empty when the file
// cannot be read now.
std::map<std::pair<int, int>, std::string> LevelLabels(const Config &config) {
  std::map<std::pair<int, int>, std::string> labels;
  std::unique_ptr<std::istream> file = config.SectionStream("levels");  // its tags and lines
  if (!file || !*file) return labels;
  std::ostringstream text;
  text << file->rdbuf();
  azure2::experiment::LevelRowTable table;
  std::string why;
  if (!azure2::experiment::ReadLevelRowTable(azure2::io::Document::Parse(text.str()), table, why)) return labels;
  for (const azure2::experiment::LevelEntry &l : table.levels)
    if (l.group > 0) labels[{l.group, l.level}] = l.label;
  return labels;
}

struct PairField {
  const char *name;
  const char *unit;
  const char *doc;
  py::object (*read)(const CNuc &, const WorkingCopy &, const PPair &, int number);
};

const std::vector<PairField> &PairFields() {
  static const std::vector<PairField> fields = {
      {"number", "", "the engine's pair number, from 1 (Parameter.pair)",
       [](const CNuc &, const WorkingCopy &, const PPair &, int p) -> py::object { return py::int_(p); }},
      {"key", "", "the pair's key in <levels> (engine messages and isospin.out name pairs by key)",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetPairKey()); }},
      {"pair_kind", "", "particle, photon (gamma + the final level) or beta (a beta decay)",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::str(q.GetPType() == 0 ? "particle" : q.GetPType() == 10 ? "photon" : "beta");
       }},
      {"ptype", "", "the pair type as <levels> writes it: 0 particle, 10 photon, 20 beta",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetPType()); }},
      {"is_entrance", "", "an entrance pair of the calculation: of an active data set (of an active extrapolation in a "
                          "session of the extrapolations), set when the session is built",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::bool_(q.IsEntrance()); }},
      {"J1", "", "the light particle's spin", [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetJ(1));
       }},
      {"parity1", "", "its parity, +1 or -1",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetPi(1)); }},
      {"Z1", "", "its charge number",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetZ(1)); }},
      {"M1", "u", "its mass",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetM(1));
       }},
      {"J2", "", "the heavy particle's spin", [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetJ(2));
       }},
      {"parity2", "", "its parity, +1 or -1",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetPi(2)); }},
      {"Z2", "", "its charge number",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetZ(2)); }},
      {"M2", "u", "its mass",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetM(2));
       }},
      {"sep_energy", "MeV", "the separation energy: the threshold is sep_energy + excitation",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::float_(q.GetSepE()); }},
      {"excitation", "MeV", "the pair's excitation energy, light plus heavy particle",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::float_(q.GetExE()); }},
      {"light_excitation", "MeV", "the light particle's excitation energy (column 32 of <levels>; 0 for gamma, beta)",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetLightExE());
       }},
      {"channel_radius", "fm", "the channel radius",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::float_(q.GetChRad()); }},
      {"i1i2factor", "", "1 / ((2 J1 + 1)(2 J2 + 1)), the statistical factor's denominator of an entrance pair",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetI1I2Factor());
       }},
      {"nucleus_Z", "", "the compound nucleus the pair forms: its charge number",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object { return py::int_(q.GetNucleusZ()); }},
      {"nucleus_A", "", "its mass number", [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::int_(q.GetNucleusA());
       }},
      {"abundance", "", "that nucleus's abundance in the default sample (1 without <abundances>)",
       [](const CNuc &, const WorkingCopy &data, const PPair &q, int) -> py::object {
         return py::float_(data.Abundance(data.Data().Abundances().Default(), q.GetNucleus()));
       }},
      {"binding_energy", "MeV", "the THM binding energy of the pair's composite projectile, 0 for other pairs",
       [](const CNuc &, const WorkingCopy &, const PPair &q, int) -> py::object {
         return py::float_(q.GetBindingEnergy());
       }},
  };
  return fields;
}

}  // namespace

py::list ParameterRecordDicts(const AZUREAPI &api) {
  const std::vector<ParameterRecord> &records = api.GetParameterRecords();
  const std::vector<bool> fixed = api.params_fixed();
  const std::map<std::pair<int, int>, std::string> labels = LevelLabels(api.configure());
  const vector_r start = api.params_all_rwa(), value = api.params_all();
  py::list out;
  int free = 0;
  for (int i = 0; i < (int)records.size(); ++i) {
    ParameterView v{&api, i, &records[i], api.compound(), &api.data()->Data()};
    v.levelLabels = &labels;
    v.start = &start;
    v.value = &value;
    v.freeIndex = fixed[i] ? -1 : free++;
    Resolve(v);
    py::dict d;
    for (const ParameterField &f : ParameterFields())
      if (InSet(f.kinds, records[i].kind)) d[f.name] = f.read(v);
    out.append(d);
  }
  return out;
}

py::list PairRecordDicts(const AZUREAPI &api) {
  const CNuc &nuc = *api.compound();
  const WorkingCopy &data = *api.data();
  py::list out;
  for (int p = 1; p <= nuc.NumPairs(); ++p) {
    py::dict d;
    for (const PairField &f : PairFields()) d[f.name] = f.read(nuc, data, *nuc.GetPair(p), p);
    out.append(d);
  }
  return out;
}

/// The parameter fields of a kind that say what it belongs to (pyazr.kinds()['parameters'][key]['owner']).
std::vector<std::string> OwnerFieldsOf(ParamKind kind) {
  std::vector<std::string> names;
  for (const ParameterField &f : ParameterFields())
    if (f.kinds != kEveryKind && InSet(f.kinds, kind)) names.push_back(f.name);
  return names;
}

void BindRecordFields(py::module_ &m) {
  m.def(
      "parameter_record_fields",
      []() {
        py::list out;
        for (const ParameterField &f : ParameterFields()) {
          py::dict d;
          d["name"] = f.name;
          d["unit"] = f.unit;
          py::list kinds;
          for (const ParamKindRow &row : ParamKindTable())
            if (InSet(f.kinds, row.kind)) kinds.append(row.key);
          d["kinds"] = kinds;
          d["every_kind"] = f.kinds == kEveryKind;
          d["owner"] = f.owner;
          d["doc"] = f.doc;
          out.append(d);
        }
        return out;
      },
      "The fields of a parameter record (Session.parameter_records), one dict each: name, unit ('unit': the "
      "parameter's own), kinds (the keys of the kinds it describes), every_kind, owner (shown as what the parameter "
      "belongs to), doc.");
  m.def(
      "pair_record_fields",
      []() {
        py::list out;
        for (const PairField &f : PairFields()) {
          py::dict d;
          d["name"] = f.name;
          d["unit"] = f.unit;
          d["doc"] = f.doc;
          out.append(d);
        }
        return out;
      },
      "The fields of a pair record (Session.pair_records), one dict each: name, unit, doc.");
}
