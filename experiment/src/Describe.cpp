// What Schema::Build made of a format-3 project, as text and JSON (Describe.h).
#include "azure2/experiment/Describe.h"

#include <filesystem>
#include <sstream>

#include "azure2/experiment/Json.h"

namespace azure2 {
namespace experiment {

namespace {

const char *ObservableName(Observable o) {
  switch (o) {
    case Observable::CrossSection: return "angle-integrated cross section";
    case Observable::CrossSectionE1: return "angle-integrated cross section, E1 part";
    case Observable::CrossSectionE2: return "angle-integrated cross section, E2 part";
    case Observable::TotalCapture: return "capture summed over every capture pair";
    case Observable::DifferentialLab: return "differential cross section, lab angle";
    case Observable::DifferentialCM: return "differential cross section, c.m. angle";
    case Observable::AnalyzingPower: return "analyzing power, c.m. angle";
    case Observable::PhaseShift: return "phase shift";
    case Observable::AngularCoefficients: return "Legendre coefficients (data)";
    case Observable::AngularDistribution: return "Legendre coefficients (computed)";
  }
  return "?";
}

/// For a ratio to a Coulomb cross section, what the observable is divided by ("" otherwise).
std::string ReferenceText(const ReactionSpec &r) {
  switch (r.coulombReference) {
    case CoulombReference::None: return "";
    case CoulombReference::Rutherford:
      return "  divided by the Rutherford cross section at the same c.m. energy and angle, folded alike; the value is "
             "the same in both frames (not converted at load)\n";
    case CoulombReference::Mott:
      return "  divided by the Mott cross section (identical particles, the particles' spin) at the same c.m. energy and "
             "angle, folded alike; the value is the same in both frames (not converted at load)\n";
  }
  return "";
}

const char *ReferenceName(CoulombReference c) {
  return c == CoulombReference::Rutherford ? "rutherford" : c == CoulombReference::Mott ? "mott" : "";
}

bool AngleObservable(Observable o) {
  return o == Observable::DifferentialLab || o == Observable::DifferentialCM || o == Observable::AnalyzingPower;
}

std::string Pair(const Format3Project &p, int key) {
  if (key == -1) return "every capture pair";
  const PairInfo *pair = PairByKey(key, p.levels.pairs);
  return std::to_string(key) + " (" + (pair ? PairLabel(*pair) : std::string("?")) + ")";
}

struct Term {
  std::string text;
  int entrance = 0, exit = 0;
  double scale = 1.0;
  bool hasAngle = false;
  double angle = 0.0;
};

std::vector<Term> Terms(const ValidatedRecord &r, const ReactionSpec &spec) {
  std::vector<Term> terms;
  const Value *v = r.Get("reaction");
  if (!v) return terms;
  for (std::size_t i = 0; i < v->reaction.terms.size(); i++) {
    Term t;
    t.text = v->reaction.terms[i].reaction.text;
    if (i == 0) {
      t.entrance = spec.entranceKey;
      t.exit = spec.exitKey;
    } else if (i - 1 < spec.components.size()) {
      const ComponentSpec &c = spec.components[i - 1];
      t.entrance = c.entranceKey;
      t.exit = c.exitKey;
      t.scale = c.scaling;
      t.hasAngle = c.hasAngle;
      t.angle = c.angle;
    }
    terms.push_back(t);
  }
  return terms;
}

struct Parameter {
  std::string name, kind, status;
};

std::vector<Parameter> DatasetParameters(const DatasetSpec &d, bool profiled) {
  std::vector<Parameter> out;
  const ParameterNames &n = d.names;
  if (profiled) out.push_back({"", "norm", "profiled (THM experiment, no parameter)"});
  else if (d.varyNorm) out.push_back({n.norm, "norm", n.norm.find('[') != std::string::npos ? "fixed (the group's is its first member's)" : "fit"});
  out.push_back({n.energyShift, "energy shift",
                 n.energyShift.find('[') != std::string::npos ? "fixed (the group's is its first member's)"
                                                              : d.varyEnergyShift ? "fit" : "fixed"});
  out.push_back({n.energyShiftSqrt, "energy shift sqrt", d.varyEnergyShiftSqrt ? "fit" : "fixed"});
  for (std::size_t j = 0; j < n.shiftFunction.size(); j++)
    out.push_back({n.shiftFunction[j], "energy shift function", j < d.shiftFunctionVary.size() && d.shiftFunctionVary[j] ? "fit" : "fixed"});
  if (d.thicknessFit)
    out.push_back({n.thickness, "thickness", n.thickness.find('[') != std::string::npos ? "fixed (the group's is its first member's)"
                                             : d.thicknessVary ? "fit" : "fixed"});
  return out;
}

std::vector<std::string> Effects(const SpecList<TargetEffectSpec> &effects, int key) {
  std::vector<std::string> out;
  for (const auto &item : effects.items)
    for (int k : item.spec.datasetKeys)
      if (k == key) out.push_back(item.spec.source);
  return out;
}

std::string ResolvedPath(const std::string &projectDir, const std::string &file) {
  std::error_code ec;
  const std::filesystem::path path = std::filesystem::path(projectDir.empty() ? "." : projectDir) / file;
  const std::filesystem::path absolute = std::filesystem::absolute(path, ec).lexically_normal();
  return ec ? path.string() : absolute.string();
}

/// The records of one kind in input order, the n-th of which is the spec with key n.
std::vector<const ValidatedRecord *> RecordsOf(const Format3Project &p, const std::string &kind) {
  std::vector<const ValidatedRecord *> out;
  for (const ValidatedRecord &r : p.validated.records)
    if (r.kind == kind) out.push_back(&r);
  return out;
}

bool Profiled(const Format3Project &p, int key) {
  for (const auto &x : p.built.thm.experiments)
    for (int k : x.datasetKeys)
      if (k == key) return true;
  return false;
}

std::string Fmt(double x) {
  std::ostringstream s;
  s << x;
  return s.str();
}

}  // namespace

std::string DescribeText(const Format3Project &p) {
  std::ostringstream out;
  const std::vector<const ValidatedRecord *> data = RecordsOf(p, "dataset"), extrapolations = RecordsOf(p, "extrapolation");
  for (std::size_t i = 0; i < data.size() && i < p.built.data.datasets.items.size(); i++) {
    const ValidatedRecord &r = *data[i];
    const DatasetSpec &d = p.built.data.datasets.items[i].spec;
    const Value *reaction = r.Get("reaction");
    out << "dataset " << r.name << " (data set " << d.key << ", line " << r.line << ")" << (d.active ? "" : ": inactive")
        << "\n";
    if (!d.active) continue;
    out << "  " << (r.quantity ? r.quantity->name : "?") << ": " << ObservableName(d.reaction.observable)
        << (reaction ? " of " + reaction->reaction.Tree() : "") << "\n";
    out << ReferenceText(d.reaction);
    for (const Term &t : Terms(r, d.reaction))
      out << "    " << t.text << ": entrance pair " << Pair(p, t.entrance) << ", exit pair " << Pair(p, t.exit)
          << (t.scale != 1.0 ? ", scale " + Fmt(t.scale) : "") << (t.hasAngle ? ", at " + Fmt(t.angle) + " deg lab" : "") << "\n";
    out << "  window E_lab " << Fmt(d.minE) << " - " << Fmt(d.maxE) << " MeV";
    if (AngleObservable(d.reaction.observable))
      out << ", theta_" << (d.reaction.observable == Observable::DifferentialLab ? "lab " : "lab (data converted to c.m.) ")
          << Fmt(d.minA) << " - " << Fmt(d.maxA) << " deg";
    out << "\n  data file " << ResolvedPath(p.projectDir, d.dataFile) << "\n";
    for (const std::string &e : Effects(p.built.data.targetEffects, d.key)) out << "  effect: " << e << "\n";
    if (d.reaction.sample != "0") out << "  sample: " << d.reaction.sample << "\n";
    out << "  parameters:";
    for (const Parameter &x : DatasetParameters(d, Profiled(p, d.key)))
      out << "\n    " << (x.name.empty() ? "-" : x.name) << " (" << x.kind << ", " << x.status << ")";
    out << "\n";
  }
  for (std::size_t i = 0; i < extrapolations.size() && i < p.built.extrapolations.extrapolations.items.size(); i++) {
    const ValidatedRecord &r = *extrapolations[i];
    const ExtrapolationSpec &x = p.built.extrapolations.extrapolations.items[i].spec;
    const Value *reaction = r.Get("reaction");
    out << "extrapolation " << r.name << " (extrapolation " << x.key << ", line " << r.line << ")"
        << (x.active ? "" : ": inactive") << "\n";
    if (!x.active) continue;
    out << "  " << (r.quantity ? r.quantity->name : "?") << ": " << ObservableName(x.reaction.observable)
        << (reaction ? " of " + reaction->reaction.Tree() : "") << "\n";
    out << ReferenceText(x.reaction);
    for (const Term &t : Terms(r, x.reaction))
      out << "    " << t.text << ": entrance pair " << Pair(p, t.entrance) << ", exit pair " << Pair(p, t.exit) << "\n";
    out << "  grid: " << x.energies.size() << (x.adaptive.on ? " start energies" : " energies") << " E_lab " << Fmt(x.minE)
        << " - " << Fmt(x.maxE) << " MeV";
    if (!x.angles.empty()) out << " at " << x.angles.size() << " angle(s) " << Fmt(x.minA) << " - " << Fmt(x.maxA) << " deg";
    out << "\n";
    if (!x.gridFile.empty()) out << "  grid file: " << x.gridFile << " (" << x.energies.size() << " energies)\n";
    if (x.adaptive.on)
      out << "  adaptive: refined in the extrapolation run to a relative tolerance " << Fmt(x.adaptive.tolerance)
          << ", at most " << x.adaptive.maxPoints << " energies\n";
    for (const std::string &e : Effects(p.built.extrapolations.targetEffects, x.key)) out << "  effect: " << e << "\n";
    out << "  parameters: " << x.names.energyShift << ", " << x.names.energyShiftSqrt << " (energy shift, fixed)\n";
  }
  for (const ModelDatum &d : p.built.modelData) {
    out << "datum " << d.name << " (line " << d.line << ")\n  " << d.quantityName << " = " << Fmt(d.value / d.unitFactor)
        << " +- " << Fmt(d.error / d.unitFactor) << " " << d.unit;
    if (d.levelId) out << "; level ID " << d.levelId << " (" << d.jpi << ", group " << d.group << ", level " << d.level << ")";
    if (d.entrancePair) out << "; entrance pair " << Pair(p, d.entrancePair);
    if (d.exitPair) out << "; exit pair " << Pair(p, d.exitPair);
    if (d.pair) out << "; pair " << Pair(p, d.pair);
    out << "\n";
  }
  return out.str();
}

std::string DescribeJson(const Format3Project &p) {
  JsonWriter j(1);
  j.BeginObject().Key("datasets").BeginArray();
  const std::vector<const ValidatedRecord *> data = RecordsOf(p, "dataset"), extrapolations = RecordsOf(p, "extrapolation");
  auto reactionTree = [&](const ValidatedRecord &r, const ReactionSpec &spec) {
    const Value *reaction = r.Get("reaction");
    j.Key("observable").String(ObservableName(spec.observable));
    j.Key("tree").String(reaction ? reaction->reaction.Tree() : "");
    j.Key("coulomb_reference").String(ReferenceName(spec.coulombReference));
    j.Key("terms").BeginArray();
    for (const Term &t : Terms(r, spec)) {
      j.BeginObject().Key("reaction").String(t.text).Key("entrance_pair").Integer(t.entrance).Key("exit_pair").Integer(t.exit);
      j.Key("scale").Number(t.scale);
      if (t.hasAngle) j.Key("angle_lab_deg").Number(t.angle);
      j.EndObject();
    }
    j.EndArray();
  };
  for (std::size_t i = 0; i < data.size() && i < p.built.data.datasets.items.size(); i++) {
    const ValidatedRecord &r = *data[i];
    const DatasetSpec &d = p.built.data.datasets.items[i].spec;
    j.BeginObject().Key("name").String(r.name).Key("key").Integer(d.key).Key("line").Integer(r.line).Key("active").Bool(d.active);
    j.Key("quantity").String(r.quantity ? r.quantity->name : "");
    reactionTree(r, d.reaction);
    j.Key("E_lab_window_MeV").BeginArray().Number(d.minE).Number(d.maxE).EndArray();
    if (AngleObservable(d.reaction.observable)) j.Key("angle_window_deg").BeginArray().Number(d.minA).Number(d.maxA).EndArray();
    j.Key("data_file").String(d.dataFile).Key("data_file_resolved").String(ResolvedPath(p.projectDir, d.dataFile));
    j.Key("effects").StringArray(Effects(p.built.data.targetEffects, d.key));
    j.Key("sample").String(d.reaction.sample == "0" ? "" : d.reaction.sample);
    j.Key("parameters").BeginArray();
    for (const Parameter &x : DatasetParameters(d, Profiled(p, d.key)))
      j.BeginObject().Key("name").String(x.name).Key("kind").String(x.kind).Key("status").String(x.status).EndObject();
    j.EndArray().EndObject();
  }
  j.EndArray().Key("extrapolations").BeginArray();
  for (std::size_t i = 0; i < extrapolations.size() && i < p.built.extrapolations.extrapolations.items.size(); i++) {
    const ValidatedRecord &r = *extrapolations[i];
    const ExtrapolationSpec &x = p.built.extrapolations.extrapolations.items[i].spec;
    j.BeginObject().Key("name").String(r.name).Key("key").Integer(x.key).Key("line").Integer(r.line).Key("active").Bool(x.active);
    j.Key("quantity").String(r.quantity ? r.quantity->name : "");
    reactionTree(r, x.reaction);
    j.Key("energies").Integer((long)x.energies.size()).Key("angles").Integer((long)x.angles.size());
    j.Key("grid_file").String(x.gridFile);
    j.Key("adaptive").Bool(x.adaptive.on);
    if (x.adaptive.on) j.Key("tolerance").Number(x.adaptive.tolerance).Key("max_points").Integer(x.adaptive.maxPoints);
    j.Key("effects").StringArray(Effects(p.built.extrapolations.targetEffects, x.key));
    j.EndObject();
  }
  j.EndArray().Key("model_data").BeginArray();
  for (const ModelDatum &d : p.built.modelData) {
    j.BeginObject().Key("name").String(d.name).Key("line").Integer(d.line).Key("quantity").String(d.quantityName);
    j.Key("value").Number(d.value).Key("error").Number(d.error).Key("unit").String(d.unit);
    j.Key("level").Integer(d.levelId).Key("entrance_pair").Integer(d.entrancePair).Key("exit_pair").Integer(d.exitPair);
    j.Key("pair").Integer(d.pair).EndObject();
  }
  j.EndArray().EndObject();
  return j.Str();
}

}  // namespace experiment
}  // namespace azure2
