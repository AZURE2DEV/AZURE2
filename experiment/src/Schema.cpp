// The record kinds of format 3 and their fields (Schema.h; docs/upgrade/project_format3.md sections 3.2-3.10).
// One function per record kind fills its RecordKind member by member; its fields are rows made by the helpers below,
// whose names say what each argument is (no positional flags).  A new field is one row in its kind's list, e.g.
//
//   RequiredWhen("charge_state", "Charge state", IntAtLeast(1), Condition::Dimension({"per_charge"}),
//                "the beam's charge state q; per-charge units need it (F8)"),
//
// its value type one of Values.h's, its structure (if any) a FormSpec of StructuredValues().
#include "azure2/experiment/Schema.h"

#include <utility>

namespace azure2 {
namespace experiment {

namespace {

// ---------------------------------------------------------------- value specs
ValueSpec Spec(ValueType t, const std::string &dimension = "") { return ValueSpec(t, dimension); }
/// Any real number.
ValueSpec Num() { return Spec(ValueType::Number); }
/// A real number > 0.
ValueSpec PositiveNum() {
  ValueSpec s = Spec(ValueType::Number);
  s.positive = true;
  return s;
}
ValueSpec IntAtLeast(long min) {
  ValueSpec s = Spec(ValueType::Integer);
  s.min = min;
  return s;
}
/// A number and a unit of `dimension`.
ValueSpec WithUnit(const std::string &dimension) { return Spec(ValueType::NumberWithUnit, dimension); }
/// A number > 0 and a unit of `dimension`.
ValueSpec PositiveWithUnit(const std::string &dimension) {
  ValueSpec s = Spec(ValueType::NumberWithUnit, dimension);
  s.positive = true;
  return s;
}
ValueSpec Nuis(const std::string &dimension) { return Spec(ValueType::Nuisance, dimension); }
ValueSpec Words(std::vector<std::string> values) {
  ValueSpec s = Spec(ValueType::Enum);
  s.values = std::move(values);
  return s;
}
/// The name of a record of one of `kinds`.
ValueSpec Ref(std::vector<std::string> kinds) {
  ValueSpec s = Spec(ValueType::Reference);
  s.values = std::move(kinds);
  return s;
}
/// The name of a record of one of `kinds`, optionally followed by @<epoch>.
ValueSpec RefWithEpoch(std::vector<std::string> kinds) {
  ValueSpec s = Ref(std::move(kinds));
  s.epoch = true;
  return s;
}
/// A grid of `dimension` (section 3.8).
ValueSpec GridOf(const std::string &dimension) { return Spec(ValueType::Grid, dimension); }
/// A grid of `dimension`, or the word 'none'.
ValueSpec GridOrNone(const std::string &dimension) {
  ValueSpec s = GridOf(dimension);
  s.none = true;
  return s;
}
ValueSpec Structured(const FormSpec &form) {
  ValueSpec s = Spec(ValueType::Form);
  s.form = &form;
  return s;
}
ValueSpec QuantityEnum(const std::string &recordKind) {
  ValueSpec s = Words(QuantityNames(recordKind));
  for (const Quantity &q : Quantities())
    for (const std::string &a : q.aliases)
      for (const std::string &k : q.recordKinds)
        if (k == recordKind) s.aliases.push_back({a, q.name});
  return s;
}

// ---------------------------------------------------------------- parts and alternatives of structured values
/// A positional part (positional parts come first, in order).
PartSpec Pos(const std::string &name, ValueSpec spec, const std::string &doc = "") {
  PartSpec p;
  p.name = name;
  p.spec = std::move(spec);
  p.doc = doc;
  return p;
}
/// A required part written "<name> <value>".
PartSpec Kw(const std::string &name, ValueSpec spec, const std::string &doc = "") {
  PartSpec p = Pos(name, std::move(spec), doc);
  p.keyword = true;
  return p;
}
/// An optional part written "<name> <value>"; `defaultText` is its value when absent ("" for none).
PartSpec OptionalKw(const std::string &name, ValueSpec spec, const std::string &doc, const std::string &defaultText = "") {
  PartSpec p = Kw(name, std::move(spec), doc);
  p.optional = true;
  p.defaultText = defaultText;
  return p;
}
/// The alternative led by the word `keyword`.
Alternative Keyword(const std::string &keyword, std::vector<PartSpec> parts, const std::string &doc) {
  Alternative a;
  a.keyword = keyword;
  a.parts = std::move(parts);
  a.doc = doc;
  return a;
}
/// The alternative chosen when no keyword of the form matches (or the form's only one).
Alternative Unkeyed(std::vector<PartSpec> parts, const std::string &doc) { return Keyword("", std::move(parts), doc); }

// ---------------------------------------------------------------- fields
const Condition kNever = Condition::Of(Condition::Always, true);

FieldDef Field(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  FieldDef f;
  f.key = key;
  f.title = title;
  f.spec = std::move(spec);
  f.doc = doc;
  return f;
}
/// Required in every record of the kind.
FieldDef Req(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  return Field(key, title, std::move(spec), doc);
}
/// Optional in every record of the kind.
FieldDef Opt(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  FieldDef f = Field(key, title, std::move(spec), doc);
  f.required = kNever;
  return f;
}
/// Where it applies and whether it is required is decided by the record's quantity row (Quantity::fields).
FieldDef ByQ(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  FieldDef f = Opt(key, title, std::move(spec), doc);
  f.byQuantity = true;
  return f;
}
/// Applies (and is optional) where `applies` holds; refused elsewhere.
FieldDef OptionalWhen(const std::string &key, const std::string &title, ValueSpec spec, Condition applies,
                      const std::string &doc) {
  FieldDef f = Opt(key, title, std::move(spec), doc);
  f.applies = std::move(applies);
  return f;
}
/// Applies and is required where `applies` holds; refused elsewhere.
FieldDef RequiredWhen(const std::string &key, const std::string &title, ValueSpec spec, Condition applies,
                      const std::string &doc) {
  FieldDef f = Field(key, title, std::move(spec), doc);
  f.applies = std::move(applies);
  return f;
}
/// The field, required (among the records it applies to) where `required` holds.
FieldDef RequiredIf(FieldDef f, Condition required) {
  f.required = std::move(required);
  return f;
}
FieldDef Repeatable(FieldDef f) {
  f.repeatable = true;
  return f;
}
FieldDef Defaulted(FieldDef f, const std::string &defaultText) {
  f.defaultText = defaultText;
  return f;
}

std::vector<FieldDef> Concat(std::vector<std::vector<FieldDef>> lists) {
  std::vector<FieldDef> all;
  for (std::vector<FieldDef> &l : lists) all.insert(all.end(), l.begin(), l.end());
  return all;
}

// ---------------------------------------------------------------- structured values (sections 3.3-3.10)
struct Forms {
  FormSpec shift, quadrature, foldWindows, secondaryGamma, orders, function, param, epoch, stopping, straggling,
      isotope, isotopeSum, doppler, spread, profile, resolution, udr, channelSpin;
};

Forms MakeForms() {
  Forms f;
  f.shift = {{Keyword("offset", {Pos("value", Nuis("energy"))}, "E' = E + offset"),
              Keyword("sqrt", {Pos("value", Nuis("sqrt_energy"))}, "E' = E + b sqrt(E), b in MeV^1/2 with E in MeV"),
              Unkeyed({Pos("scale", Ref({"energy_scale"}))}, "the terms of a shared energy_scale record")}};
  f.quadrature = {{Unkeyed({Kw("points", IntAtLeast(1), "sub-points of the fold; no default (every <targetInt> line states it)"),
                            OptionalKw("width_multiplier", PositiveNum(), "resonance-width multiplier of the adaptive grid", "20"),
                            OptionalKw("per_width", PositiveNum(), "points per resonance width", "50")},
                           "")}};
  f.foldWindows = {{Unkeyed({Pos("windows", Spec(ValueType::Windows, "energy"), "lab-energy windows, or auto[:tolerance[:margin]]"),
                             OptionalKw("blend", WithUnit("energy"), "blend width at a window edge"),
                             OptionalKw("auto_tolerance", PositiveNum(), "tolerance of the automatic windows")},
                            "")}};
  f.secondaryGamma = {{Unkeyed({Kw("L", IntAtLeast(1), "multipolarity"), Kw("Ic", Spec(ValueType::Spin), "spin of the intermediate state"),
                                Kw("delta", Num(), "mixing ratio")},
                               "unobserved primary, observed secondary photon (UPOS)")}};
  f.orders = {{Unkeyed({Pos("first", IntAtLeast(0)), Pos("last", IntAtLeast(0))}, "the Legendre orders first..last")}};
  f.function = {{Unkeyed({Pos("expression", Spec(ValueType::Formula), "x = E_lab of the point in MeV"),
                          Pos("unit", Spec(ValueType::Unit, "energy"))},
                         "a shift term function(x) in the stated unit; its parameters are 'param' fields")}};
  f.param = {{Unkeyed({Pos("name", Spec(ValueType::Name)), Pos("value", Nuis(""))}, "")}};
  f.epoch = {{Unkeyed({Pos("name", Spec(ValueType::Name)), Kw("areal_density", Nuis("areal_density"))},
                      "the areal density of the target in one epoch of a campaign (data sets name target@epoch)")}};
  f.stopping = {{Keyword("formula", {Pos("expression", Spec(ValueType::Formula), "x = E_lab/MeV"), Pos("unit", Spec(ValueType::Unit, "stopping"))},
                         "stopping cross section per active atom"),
                 Keyword("evaluation", {Kw("scale", Nuis(""))}, "the evaluated stopping of the composition (Phase 5), scaled")}};
  f.straggling = {{Keyword("none", {}, "no straggling"),
                   Keyword("bohr", {}, "Bohr straggling from the entrance charges and the stopping"),
                   Keyword("coefficient", {Pos("value", PositiveWithUnit("sqrt_energy")), Pos("frame", Words({"cm"}))},
                           "sigma = c sqrt(dE) with c and dE in the c.m. frame")}};
  f.isotope = {{Unkeyed({Pos("nucleus", Spec(ValueType::Nucleus)), Pos("abundance", Nuis(""))}, "abundance weight of a target nucleus")}};
  f.isotopeSum = {{Keyword("none", {}, "no constraint on the sum of the abundances"),
                   Keyword("prior", {Pos("centre", Num()), Pos("width", PositiveNum())}, "a Gaussian prior on the sum of the abundances")}};
  f.doppler = {{Keyword("free_gas", {Pos("temperature", Nuis("temperature"))}, "free-gas Doppler broadening of the target atoms")}};
  f.spread = {{Unkeyed({Pos("sigma", WithUnit("energy")), Pos("frame", Words({"lab", "cm"}))}, "a constant Gaussian sigma"),
               Keyword("formula",
                       {Pos("expression", Spec(ValueType::Formula), "x = E_lab/MeV"), Pos("unit", Spec(ValueType::Unit, "energy")),
                        Pos("frame", Words({"lab", "cm"}))},
                       "sigma(E_lab); 'cm' is the convention before 2026-10")}};
  f.profile = {{Keyword("skewed_gaussian",
                        {Kw("xi", WithUnit("energy")), Kw("omega", PositiveWithUnit("energy")), Kw("alpha", Num()), Kw("weight", PositiveNum())},
                        "one skewed Gaussian of the beam energy distribution")}};
  f.resolution = {{Unkeyed({Pos("sigma", WithUnit("energy")), Pos("frame", Words({"lab", "cm"}))}, "a Gaussian energy resolution")}};
  f.udr = {{Unkeyed({Kw("file", Spec(ValueType::Path)), Kw("flight_path", PositiveWithUnit("length")), Kw("burst", WithUnit("time")),
                     Kw("channel", WithUnit("time")), Kw("centred", Words({"yes", "no"}))},
                    "SAMMY's time-of-flight resolution (UDR)")}};
  f.channelSpin = {{Keyword("coherent", {}, "the coherent scattering length"), Unkeyed({Pos("spin", Spec(ValueType::Spin))}, "the channel spin")}};
  return f;
}

/// The structured values, built once (ValueSpec::form points into this table).
const Forms &StructuredValues() {
  static const Forms forms = MakeForms();
  return forms;
}

// ---------------------------------------------------------------- fields shared by several kinds
std::vector<FieldDef> ApparatusFields() {
  const Forms &f = StructuredValues();
  const Condition fold = Condition::Of(Condition::Fold);
  return {ByQ("target", "Target", RefWithEpoch({"target"}), "the target record (and epoch) the beam passed through"),
          ByQ("beam", "Beam", Ref({"beam"}), "the beam record: energy spread or profile"),
          ByQ("detector", "Detector", Ref({"detector"}), "the detector record: Q_k, resolution, UDR"),
          RequiredWhen("quadrature", "Quadrature", Structured(f.quadrature), fold, "sub-point grid of the fold"),
          OptionalWhen("fold_windows", "Fold windows", Structured(f.foldWindows), fold, "lab-energy windows of the fold")};
}

std::vector<FieldDef> QuantityFields() {
  const Forms &f = StructuredValues();
  return {ByQ("J", "J", Spec(ValueType::Spin), "total spin of the phase shift"),
          ByQ("l", "l", IntAtLeast(0), "orbital angular momentum of the phase shift"),
          ByQ("s", "s", Spec(ValueType::Spin), "channel spin of the phase shift"),
          ByQ("expansion", "Legendre expansion", Words({"a_over_a0", "absolute", "rsl", "2l2", "4pi"}),
              "the normalisation of the Legendre coefficients (data_kinds.md section 5.2)"),
          ByQ("multipole", "Multipole", Words({"E1", "E2"}), "the photon multipole the cross section is filtered to (M8)"),
          ByQ("secondary_gamma", "Secondary gamma", Structured(f.secondaryGamma), "UPOS: L, Ic, delta"),
          ByQ("charge_state", "Charge state", IntAtLeast(1), "the beam's charge state q; per-charge units need it (F8)"),
          ByQ("efficiency", "Efficiency x branching", Nuis(""), "detection efficiency times branching of counted photons (F9)"),
          ByQ("thm", "THM experiment", Spec(ValueType::Name), "the experiment[<name>] of <thm>, or global"),
          OptionalWhen("thm_weight", "THM weight file", Spec(ValueType::Path), Condition::Present("thm"),
                       "the weight table of <thm> weight[k]")};
}

std::vector<FieldDef> TextFields() {
  return {Opt("source", "Source", Spec(ValueType::Text), "where the numbers come from; never read by the physics"),
          Opt("exfor", "EXFOR", Spec(ValueType::Text), "the EXFOR entry and subentry"),
          Opt("doi", "DOI", Spec(ValueType::Text), "the DOI of the publication")};
}

// ---------------------------------------------------------------- record kinds (section 3.2)
RecordKind DatasetKind() {
  const Forms &f = StructuredValues();
  ValueSpec norm = Nuis("");
  norm.profiled = true;
  norm.refKinds = {"normalization"};
  RecordKind k;
  k.kind = "dataset";
  k.sections = {"data"};
  k.title = "Data set";
  k.doc = "Measured values of one quantity in one data file (section 3.4).";
  k.fields = Concat({{Req("reaction", "Reaction", Spec(ValueType::Reaction), "the reaction expression (section 3.5)"),
                      Req("quantity", "Quantity", QuantityEnum("dataset"), "what the numbers are (section 7.2)"),
                      Req("unit", "Unit", Spec(ValueType::Unit, "@quantity"), "unit of the value and error columns"),
                      Req("file", "Data file", Spec(ValueType::Path), "the data file, relative to the project"),
                      Req("columns", "Columns", Spec(ValueType::Columns), "role of each column (section 3.4.1)"),
                      Opt("energy_window", "Energy window", Spec(ValueType::Range, "energy"),
                          "rows outside are not used; the frame is the energy column's"),
                      OptionalWhen("angle_window", "Angle window", Spec(ValueType::Range, "angle"),
                                   Condition::Of(Condition::AngleAxis), "rows outside are not used"),
                      Req("norm", "Normalization", norm, "the normalization nuisance, or the name of a shared normalization record"),
                      Repeatable(Opt("shift", "Energy shift", Structured(f.shift), "a term of the energy shift; terms add")),
                      Opt("covariance", "Covariance file", Spec(ValueType::Path), "data_formats.rst, covariance files"),
                      Defaulted(Opt("active", "Active", Words({"yes", "no"}), "an inactive data set is not used"), "yes")},
                     ApparatusFields(), QuantityFields(), TextFields()});
  return k;
}

RecordKind NormalizationKind() {
  RecordKind k;
  k.kind = "normalization";
  k.sections = {"data"};
  k.title = "Shared normalization";
  k.doc = "One normalization shared by the data sets that name it (F6).";
  k.fields = {Req("value", "Value", Nuis(""), "the normalization nuisance")};
  return k;
}

RecordKind EnergyScaleKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "energy_scale";
  k.sections = {"data"};
  k.title = "Shared energy scale";
  k.doc = "Energy-shift terms shared by the data sets that name it (F6); summed in the order offset, sqrt, function.";
  k.fields = {Opt("offset", "Offset", Nuis("energy"), "E' = E + offset"),
              Opt("sqrt", "sqrt(E) term", Nuis("sqrt_energy"), "E' = E + b sqrt(E)"),
              Opt("function", "Function", Structured(f.function), "E' = E + f(x), x = E_lab/MeV"),
              Repeatable(OptionalWhen("param", "Parameter", Structured(f.param), Condition::Present("function"),
                                      "a parameter of the function"))};
  k.oneOf = {{"offset", "sqrt", "function"}};
  return k;
}

RecordKind ExtrapolationKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "extrapolation";
  k.sections = {"extrapolations"};
  k.title = "Extrapolation";
  k.doc = "Model values on a grid, without data (section 3.8, F20).";
  k.fields = Concat({{Req("reaction", "Reaction", Spec(ValueType::Reaction), "the reaction expression (section 3.5)"),
                      Req("quantity", "Quantity", QuantityEnum("extrapolation"), "what is computed (section 7.2)"),
                      Req("unit", "Unit", Spec(ValueType::Unit, "@quantity"), "unit of the output values"),
                      Req("energies", "Energies", GridOf("energy"), "the energy grid (F10: E_k = lo + k step)"),
                      Opt("angles", "Angles", GridOrNone("angle"), "the angle grid; required for quantities with an angle"),
                      Defaulted(Opt("active", "Active", Words({"yes", "no"}), "an inactive extrapolation is not computed"), "yes"),
                      ByQ("orders", "Legendre orders", Structured(f.orders), "the Legendre orders computed")},
                     ApparatusFields(), QuantityFields()});
  k.refused = {{"file", "an extrapolation has no data file (F20)"},
               {"columns", "an extrapolation has no data file (F20)"},
               {"norm", "an extrapolation has no normalization (F20)"},
               {"shift", "an extrapolation has no energy shift (F20)"},
               {"covariance", "an extrapolation has no data (F20)"},
               {"energy_window", "the grid is the extrapolation's energy range (F20)"},
               {"angle_window", "the grid is the extrapolation's angle range (F20)"}};
  return k;
}

RecordKind TargetKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "target";
  k.sections = {"targets"};
  k.title = "Target";
  k.doc = "The target the beam passes through: thickness, stopping, straggling, isotopes, Doppler, screening (F7).";
  k.fields = {Req("kind", "Kind", Words({"layer", "infinite", "sample", "none"}),
                  "layer: finite areal density; infinite: thick target; sample: neutron sample (atoms/b); none: isotopes only"),
              RequiredWhen("active", "Active nucleus", Spec(ValueType::Nucleus), Condition::Is("kind", {"layer", "infinite", "sample"}),
                           "the nucleus the areal density and the stopping count"),
              RequiredIf(OptionalWhen("composition", "Composition", Spec(ValueType::Composition), Condition::Is("kind", {"layer", "infinite"}),
                                      "element stoichiometry; needed by 'stopping evaluation'"),
                         Condition::Is("stopping", {"evaluation"})),
              RequiredWhen("areal_density", "Areal density", Nuis("areal_density"), Condition::Is("kind", {"layer", "sample"}),
                           "areal density of the active nucleus"),
              Repeatable(OptionalWhen("epoch", "Epoch", Structured(f.epoch), Condition::Present("areal_density"),
                                      "the areal density in one epoch (K5)")),
              RequiredWhen("stopping", "Stopping", Structured(f.stopping), Condition::Is("kind", {"layer", "infinite"}),
                           "stopping cross section"),
              RequiredWhen("straggling", "Straggling", Structured(f.straggling), Condition::Is("kind", {"layer", "infinite"}),
                           "energy straggling in the target"),
              RequiredWhen("phase", "Phase", Words({"solid", "gas"}), Condition::Is("stopping", {"evaluation"}),
                           "solid or gas stopping (with 'stopping evaluation')"),
              Repeatable(Opt("isotope", "Isotope", Structured(f.isotope), "abundance weight of a target nucleus")),
              OptionalWhen("isotope_sum", "Isotope sum", Structured(f.isotopeSum), Condition::Present("isotope"),
                           "constraint on the sum of the abundances"),
              Opt("doppler", "Doppler", Structured(f.doppler), "free-gas Doppler broadening"),
              Opt("screening", "Electron screening", Nuis("energy"), "the screening potential U_e")};
  return k;
}

RecordKind BeamKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "beam";
  k.sections = {"targets"};
  k.title = "Beam";
  k.doc = "The beam's energy distribution (F7).";
  k.fields = {Opt("spread", "Energy spread", Structured(f.spread), "Gaussian energy spread"),
              Repeatable(Opt("profile", "Profile", Structured(f.profile), "a component of the beam energy distribution")),
              RequiredWhen("truncation", "Truncation", Num(), Condition::Present("profile"), "truncation of the profile"),
              RequiredWhen("detailed_balance", "Detailed balance", Words({"yes", "no"}), Condition::Present("profile"),
                           "the profile kernel by detailed balance")};
  k.oneOf = {{"spread", "profile"}};
  return k;
}

RecordKind DetectorKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "detector";
  k.sections = {"targets"};
  k.title = "Detector";
  k.doc = "Finite solid angle, energy resolution, time-of-flight resolution (F7).";
  k.fields = {Opt("attenuation", "Attenuation Q_k", Spec(ValueType::NumberList), "Q_0 ... Q_k, lab"),
              Opt("energy_resolution", "Energy resolution", Structured(f.resolution), "resolution of a reconstructed energy"),
              Opt("udr", "UDR", Structured(f.udr), "time-of-flight resolution")};
  k.oneOf = {{"attenuation", "energy_resolution", "udr"}};
  return k;
}

RecordKind DatumKind() {
  const Forms &f = StructuredValues();
  RecordKind k;
  k.kind = "datum";
  k.sections = {"modelData"};
  k.title = "Model datum";
  k.doc = "A literature value of a level property (section 3.10, M10).";
  k.fields = Concat({{Req("quantity", "Quantity", QuantityEnum("datum"), "what the value is (section 7.2)"),
                      ByQ("level", "Level", IntAtLeast(0), "the level ID of <levels> (column 9)"),
                      ByQ("J", "J^pi", Spec(ValueType::SpinParity), "must equal the level's J^pi (F13)"),
                      ByQ("reaction", "Reaction", Spec(ValueType::Reaction), "entrance and exit"),
                      ByQ("frame", "Frame", Words({"E_x", "E_cm", "E_lab"}), "frame of a resonance energy"),
                      ByQ("pair", "Pair", Spec(ValueType::Pair), "the particle pair"),
                      ByQ("l", "l", IntAtLeast(0), "orbital angular momentum"),
                      ByQ("s", "s", Structured(f.channelSpin), "channel spin, or coherent"),
                      Req("value", "Value", WithUnit("@quantity"), "the value"),
                      Req("error", "Error", PositiveWithUnit("@quantity"), "its 1 sigma uncertainty"),
                      ByQ("imag", "Imaginary part", WithUnit("@quantity"), "Im of a scattering length"),
                      ByQ("imag_error", "Imaginary part error", PositiveWithUnit("@quantity"), "its 1 sigma uncertainty")},
                     TextFields()});
  return k;
}

}  // namespace

std::string Grammar(const ValueSpec &spec) {
  auto join = [](const std::vector<std::string> &v, const char *sep) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); i++) s += (i ? sep : "") + v[i];
    return s;
  };
  const std::string unit = spec.dimension.empty()            ? ""
                           : spec.dimension == "@quantity" ? " <unit of the quantity>"
                           : spec.dimension == "*"         ? " <unit>"
                                                           : " <" + spec.dimension + " unit>";
  switch (spec.type) {
    case ValueType::Number: return spec.positive ? "<number > 0>" : "<number>";
    case ValueType::Integer: return spec.min == LONG_MIN ? "<integer>" : "<integer >= " + std::to_string(spec.min) + ">";
    case ValueType::Unit: return unit.empty() ? "<unit>" : unit.substr(1);
    case ValueType::NumberWithUnit: return "<number>" + unit;
    case ValueType::Range: return "<lo> <hi>" + unit;
    case ValueType::Nuisance: {
      std::string g = "<value>" + (spec.dimension == "*" ? std::string(" [<unit>]") : unit) +
                      " fixed|fit [prior <centre> <width>[%]] [grid_max <value>]";
      if (spec.profiled) g += " | profiled";
      if (!spec.refKinds.empty()) g += " | <name of a " + join(spec.refKinds, " or ") + ">";
      return g;
    }
    case ValueType::Enum: return spec.values.empty() ? "<one word of the field's list>" : join(spec.values, "|");
    case ValueType::Reference: return "<name of a " + join(spec.values, " or ") + ">" + (spec.epoch ? "[@<epoch>]" : "");
    case ValueType::Name: return "<name>";
    case ValueType::Path: return "<path>";
    case ValueType::Text: return "<text>";
    case ValueType::Formula: return "\"<expression in x>\" [<name>=<number> ...]";
    case ValueType::NumberList: return "<number> ...";
    case ValueType::Spin: return "<J>";
    case ValueType::SpinParity: return "<J><parity>";
    case ValueType::Nucleus: return "<A><Symbol>";
    case ValueType::Pair: return "<particle>+<nucleus>[<E_x>]";
    case ValueType::Composition: return "<Symbol> <atoms> ...";
    case ValueType::Reaction: return "<reaction> { + <reaction> } | <reaction> / <reaction>";
    case ValueType::Columns: return "<role> <role> ... | from_file";
    case ValueType::Grid:
      return "range <lo> <hi>" + unit + " lab|cm step <s> | range <lo> <hi>" + unit + " lab|cm points <n> [log] | list <v> ..." +
             unit + " lab|cm | file <path>" + unit + " lab|cm | adaptive <lo> <hi>" + unit +
             " lab|cm tolerance <t> max_points <n> start_step <s>" + (spec.none ? " | none" : "");
    case ValueType::Windows: return "<lo>-<hi> ..." + unit + " | auto[:<tolerance>[:<margin>]]";
    case ValueType::Form: {
      std::vector<std::string> alts;
      for (const Alternative &a : spec.form->alternatives) {
        std::string s = a.keyword;
        for (const PartSpec &p : a.parts) {
          std::string part = p.keyword ? p.name + " " + Grammar(p.spec) : Grammar(p.spec);
          if (p.optional) part = "[" + part + "]";
          s += (s.empty() ? "" : " ") + part;
        }
        alts.push_back(s);
      }
      return join(alts, " | ");
    }
  }
  return "";
}

const char *ConditionKindName(Condition::Kind k) {
  switch (k) {
    case Condition::Always: return "always";
    case Condition::FieldIs: return "field_is";
    case Condition::FieldPresent: return "field_present";
    case Condition::UnitDimension: return "unit_dimension";
    case Condition::RatioReaction: return "ratio_reaction";
    case Condition::AngleAxis: return "angle_axis";
    case Condition::GammaExit: return "gamma_exit";
    case Condition::ThmExperiment: return "thm_experiment";
    case Condition::Fold: return "fold";
  }
  return "?";
}

std::string Describe(const Condition &c) {
  std::string values;
  for (std::size_t i = 0; i < c.values.size(); i++) values += (i ? " or " : "") + c.values[i];
  const char *no = c.negate ? " not" : "";
  switch (c.kind) {
    case Condition::Always: return c.negate ? "never" : "always";
    case Condition::FieldIs: return c.field + " is" + no + " " + values;
    case Condition::FieldPresent: return c.field + " is" + (c.negate ? " not" : "") + " given";
    case Condition::UnitDimension: return "the unit is" + std::string(no) + " " + values;
    case Condition::RatioReaction: return "the reaction is" + std::string(no) + " a ratio";
    case Condition::AngleAxis: return "the record has" + std::string(c.negate ? " no" : " an") + " angle axis (a theta column or an angle grid)";
    case Condition::GammaExit: return "the reaction has" + std::string(c.negate ? " no" : " a") + " photon exit";
    case Condition::ThmExperiment: return "thm" + std::string(c.negate ? " does not name" : " names") + " an experiment of <thm>";
    case Condition::Fold:
      return std::string("the record is") + no + " folded (a target layer or infinite target, a Doppler-broadening target, a beam, "
             "or a detector resolution)";
  }
  return "";
}

const FieldDef *RecordKind::Find(const std::string &key) const {
  for (const FieldDef &f : fields)
    if (f.key == key) return &f;
  return nullptr;
}

const std::vector<std::pair<std::string, std::vector<std::string>>> &Schema::Sections() {
  static const std::vector<std::pair<std::string, std::vector<std::string>>> sections = {
      {"data", {"dataset", "normalization", "energy_scale"}},
      {"targets", {"target", "beam", "detector"}},
      {"extrapolations", {"extrapolation"}},
      {"modelData", {"datum"}},
  };
  return sections;
}

std::vector<std::string> Schema::KindsOfSection(const std::string &section) {
  for (const auto &s : Sections())
    if (s.first == section) return s.second;
  return {};
}

const RecordKind *Schema::FindKind(const std::string &kind) {
  for (const RecordKind &k : Kinds())
    if (k.kind == kind) return &k;
  return nullptr;
}

const std::vector<RecordKind> &Schema::Kinds() {
  static const std::vector<RecordKind> kinds = {DatasetKind(),       NormalizationKind(), EnergyScaleKind(),
                                                ExtrapolationKind(), TargetKind(),        BeamKind(),
                                                DetectorKind(),      DatumKind()};
  return kinds;
}

bool ValidatedRecord::Has(const std::string &key) const { return Get(key) != nullptr; }

const Value *ValidatedRecord::Get(const std::string &key) const {
  for (const ValidatedField &f : fields)
    if (f.key == key) return &f.value;
  return nullptr;
}

std::vector<const Value *> ValidatedRecord::GetAll(const std::string &key) const {
  std::vector<const Value *> out;
  for (const ValidatedField &f : fields)
    if (f.key == key) out.push_back(&f.value);
  return out;
}

bool ValidationResult::Ok() const {
  for (const Diagnostic &d : diagnostics)
    if (d.severity == Severity::Error) return false;
  return true;
}

const ValidatedRecord *ValidationResult::Find(const std::string &name) const {
  for (const ValidatedRecord &r : records)
    if (r.name == name) return &r;
  return nullptr;
}

}  // namespace experiment
}  // namespace azure2
