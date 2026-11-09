// The record kinds of format 3 and their fields (Schema.h; docs/upgrade/project_format3.md sections 3.2-3.10).
// A new field is one row of its kind below; its value type is one of Values.h's, its structure (if any) a FormSpec.
#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace experiment {

namespace {

// ---------------------------------------------------------------- value specs
ValueSpec Spec(ValueType t, const std::string &dimension = "") { return ValueSpec(t, dimension); }
ValueSpec Num(bool positive = false) {
  ValueSpec s = Spec(ValueType::Number);
  s.positive = positive;
  return s;
}
ValueSpec Int(long min) {
  ValueSpec s = Spec(ValueType::Integer);
  s.min = min;
  return s;
}
ValueSpec WithUnit(const std::string &dimension, bool positive = false) {
  ValueSpec s = Spec(ValueType::NumberWithUnit, dimension);
  s.positive = positive;
  return s;
}
ValueSpec Nuis(const std::string &dimension) { return Spec(ValueType::Nuisance, dimension); }
ValueSpec Words(std::vector<std::string> values) {
  ValueSpec s = Spec(ValueType::Enum);
  s.values = std::move(values);
  return s;
}
ValueSpec Ref(std::vector<std::string> kinds, bool epoch = false) {
  ValueSpec s = Spec(ValueType::Reference);
  s.values = std::move(kinds);
  s.epoch = epoch;
  return s;
}
ValueSpec GridOf(const std::string &dimension, bool none) {
  ValueSpec s = Spec(ValueType::Grid, dimension);
  s.none = none;
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

PartSpec Pos(const std::string &name, ValueSpec spec, const std::string &doc = "") {
  return {name, false, std::move(spec), false, "", doc};
}
PartSpec Kw(const std::string &name, ValueSpec spec, const std::string &doc = "", bool optional = false,
            const std::string &defaultText = "") {
  return {name, true, std::move(spec), optional, defaultText, doc};
}

// ---------------------------------------------------------------- fields
const Condition kNever = Condition::Of(Condition::Always, true);

FieldDef Field(const std::string &key, const std::string &title, ValueSpec spec, Condition applies, Condition required,
               const std::string &doc) {
  FieldDef f;
  f.key = key;
  f.title = title;
  f.spec = std::move(spec);
  f.applies = std::move(applies);
  f.required = std::move(required);
  f.doc = doc;
  return f;
}
FieldDef Req(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  return Field(key, title, std::move(spec), {}, {}, doc);
}
FieldDef Opt(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  return Field(key, title, std::move(spec), {}, kNever, doc);
}
FieldDef ByQ(const std::string &key, const std::string &title, ValueSpec spec, const std::string &doc) {
  FieldDef f = Field(key, title, std::move(spec), {}, kNever, doc);
  f.byQuantity = true;
  return f;
}
FieldDef When(const std::string &key, const std::string &title, ValueSpec spec, Condition applies, bool required,
              const std::string &doc) {
  return Field(key, title, std::move(spec), applies, required ? Condition() : kNever, doc);
}
FieldDef Repeatable(FieldDef f) {
  f.repeatable = true;
  return f;
}
FieldDef Defaulted(FieldDef f, const std::string &defaultText) {
  f.defaultText = defaultText;
  return f;
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
      return std::string("the record is") + no + " folded (a target layer, infinite target or sample, a Doppler-broadening target, a "
             "beam, or a detector resolution)";
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
  // ---------------------------------------------------------------- structured values (section 3.3-3.10)
  static const FormSpec kShift = {{
      {"offset", {Pos("value", Nuis("energy"))}, "E' = E + offset"},
      {"sqrt", {Pos("value", Nuis("sqrt_energy"))}, "E' = E + b sqrt(E), b in MeV^1/2 with E in MeV"},
      {"", {Pos("scale", Ref({"energy_scale"}))}, "the terms of a shared energy_scale record"},
  }};
  static const FormSpec kQuadrature = {{
      {"",
       {Kw("points", Int(1), "sub-points of the fold; no default (every <targetInt> line states it)"),
        Kw("width_multiplier", Num(true), "resonance-width multiplier of the adaptive grid", true, "20"),
        Kw("per_width", Num(true), "points per resonance width", true, "50")},
       ""},
  }};
  static const FormSpec kFoldWindows = {{
      {"",
       {Pos("windows", Spec(ValueType::Windows, "energy"), "lab-energy windows, or auto[:tolerance[:margin]]"),
        Kw("blend", WithUnit("energy"), "blend width at a window edge", true),
        Kw("auto_tolerance", Num(true), "tolerance of the automatic windows", true)},
       ""},
  }};
  static const FormSpec kSecondaryGamma = {{
      {"", {Kw("L", Int(1), "multipolarity"), Kw("Ic", Spec(ValueType::Spin), "spin of the intermediate state"),
            Kw("delta", Num(), "mixing ratio")},
       "unobserved primary, observed secondary photon (UPOS)"},
  }};
  static const FormSpec kOrders = {{
      {"", {Pos("first", Int(0)), Pos("last", Int(0))}, "the Legendre orders first..last"},
  }};
  static const FormSpec kFunction = {{
      {"", {Pos("expression", Spec(ValueType::Formula), "x = E_lab of the point in MeV"), Pos("unit", Spec(ValueType::Unit, "energy"))},
       "a shift term function(x) in the stated unit; its parameters are 'param' fields"},
  }};
  static const FormSpec kParam = {{{"", {Pos("name", Spec(ValueType::Name)), Pos("value", Nuis(""))}, ""}}};
  static const FormSpec kEpoch = {{
      {"", {Pos("name", Spec(ValueType::Name)), Kw("areal_density", Nuis("areal_density"))},
       "the areal density of the target in one epoch of a campaign (data sets name target@epoch)"},
  }};
  static const FormSpec kStopping = {{
      {"formula", {Pos("expression", Spec(ValueType::Formula), "x = E_lab/MeV"), Pos("unit", Spec(ValueType::Unit, "stopping"))},
       "stopping cross section per active atom"},
      {"evaluation", {Kw("scale", Nuis(""))}, "the evaluated stopping of the composition (Phase 5), scaled"},
  }};
  static const FormSpec kStraggling = {{
      {"none", {}, "no straggling"},
      {"bohr", {}, "Bohr straggling from the entrance charges and the stopping"},
      {"coefficient", {Pos("value", WithUnit("sqrt_energy", true)), Pos("frame", Words({"cm"}))},
       "sigma = c sqrt(dE) with c and dE in the c.m. frame"},
  }};
  static const FormSpec kIsotope = {{
      {"", {Pos("nucleus", Spec(ValueType::Nucleus)), Pos("abundance", Nuis(""))}, "abundance weight of a target nucleus"},
  }};
  static const FormSpec kIsotopeSum = {{
      {"none", {}, "no constraint on the sum of the abundances"},
      {"prior", {Pos("centre", Num()), Pos("width", Num(true))}, "a Gaussian prior on the sum of the abundances"},
  }};
  static const FormSpec kDoppler = {{
      {"free_gas", {Pos("temperature", Nuis("temperature"))}, "free-gas Doppler broadening of the target atoms"},
  }};
  static const FormSpec kSpread = {{
      {"", {Pos("sigma", WithUnit("energy")), Pos("frame", Words({"lab", "cm"}))}, "a constant Gaussian sigma"},
      {"formula",
       {Pos("expression", Spec(ValueType::Formula), "x = E_lab/MeV"), Pos("unit", Spec(ValueType::Unit, "energy")),
        Pos("frame", Words({"lab", "cm"}))},
       "sigma(E_lab); 'cm' is the convention before 2026-10"},
  }};
  static const FormSpec kProfile = {{
      {"skewed_gaussian",
       {Kw("xi", WithUnit("energy")), Kw("omega", WithUnit("energy", true)), Kw("alpha", Num()), Kw("weight", Num(true))},
       "one skewed Gaussian of the beam energy distribution"},
  }};
  static const FormSpec kResolution = {{
      {"", {Pos("sigma", WithUnit("energy")), Pos("frame", Words({"lab", "cm"}))}, "a Gaussian energy resolution"},
  }};
  static const FormSpec kUdr = {{
      {"",
       {Kw("file", Spec(ValueType::Path)), Kw("flight_path", WithUnit("length", true)), Kw("burst", WithUnit("time")),
        Kw("channel", WithUnit("time")), Kw("centred", Words({"yes", "no"}))},
       "SAMMY's time-of-flight resolution (UDR)"},
  }};
  static const FormSpec kChannelSpin = {{
      {"coherent", {}, "the coherent scattering length"},
      {"", {Pos("spin", Spec(ValueType::Spin))}, "the channel spin"},
  }};

  // ---------------------------------------------------------------- shared fields
  ValueSpec norm = Nuis("");
  norm.profiled = true;
  norm.refKinds = {"normalization"};
  const Condition kFold = Condition::Of(Condition::Fold);
  const Condition kThm = Condition::Present("thm");
  const std::vector<FieldDef> apparatus = {
      ByQ("target", "Target", Ref({"target"}, true), "the target record (and epoch) the beam passed through"),
      ByQ("beam", "Beam", Ref({"beam"}), "the beam record: energy spread or profile"),
      ByQ("detector", "Detector", Ref({"detector"}), "the detector record: Q_k, resolution, UDR"),
      When("quadrature", "Quadrature", Structured(kQuadrature), kFold, true, "sub-point grid of the fold"),
      When("fold_windows", "Fold windows", Structured(kFoldWindows), kFold, false, "lab-energy windows of the fold"),
  };
  const std::vector<FieldDef> quantityFields = {
      ByQ("J", "J", Spec(ValueType::Spin), "total spin of the phase shift"),
      ByQ("l", "l", Int(0), "orbital angular momentum of the phase shift"),
      ByQ("s", "s", Spec(ValueType::Spin), "channel spin of the phase shift"),
      ByQ("expansion", "Legendre expansion", Words({"a_over_a0", "absolute", "rsl", "2l2", "4pi"}),
          "the normalisation of the Legendre coefficients (data_kinds.md section 5.2)"),
      ByQ("multipole", "Multipole", Words({"E1", "E2"}), "the photon multipole the cross section is filtered to (M8)"),
      ByQ("secondary_gamma", "Secondary gamma", Structured(kSecondaryGamma), "UPOS: L, Ic, delta"),
      ByQ("charge_state", "Charge state", Int(1), "the beam's charge state q; per-charge units need it (F8)"),
      ByQ("efficiency", "Efficiency x branching", Nuis(""), "detection efficiency times branching of counted photons (F9)"),
      ByQ("thm", "THM experiment", Spec(ValueType::Name), "the experiment[<name>] of <thm>, or global"),
      When("thm_weight", "THM weight file", Spec(ValueType::Path), kThm, false, "the weight table of <thm> weight[k]"),
  };
  auto concat = [](std::vector<FieldDef> a, const std::vector<FieldDef> &b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
  };
  const std::vector<FieldDef> text = {
      Opt("source", "Source", Spec(ValueType::Text), "where the numbers come from; never read by the physics"),
      Opt("exfor", "EXFOR", Spec(ValueType::Text), "the EXFOR entry and subentry"),
      Opt("doi", "DOI", Spec(ValueType::Text), "the DOI of the publication"),
  };

  static const std::vector<RecordKind> kinds = {
      {"dataset", {"data"}, "Data set", "Measured values of one quantity in one data file (section 3.4).",
       concat(concat(concat(std::vector<FieldDef>{Req("reaction", "Reaction", Spec(ValueType::Reaction), "the reaction expression (section 3.5)"),
                             Req("quantity", "Quantity", QuantityEnum("dataset"), "what the numbers are (section 7.2)"),
                             Req("unit", "Unit", Spec(ValueType::Unit, "@quantity"), "unit of the value and error columns"),
                             Req("file", "Data file", Spec(ValueType::Path), "the data file, relative to the project"),
                             Req("columns", "Columns", Spec(ValueType::Columns), "role of each column (section 3.4.1)"),
                             Opt("energy_window", "Energy window", Spec(ValueType::Range, "energy"),
                                 "rows outside are not used; the frame is the energy column's"),
                             When("angle_window", "Angle window", Spec(ValueType::Range, "angle"),
                                  Condition::Of(Condition::AngleAxis), false, "rows outside are not used"),
                             Req("norm", "Normalization", norm, "the normalization nuisance, or the name of a shared normalization record"),
                             Repeatable(Opt("shift", "Energy shift", Structured(kShift), "a term of the energy shift; terms add")),
                             Opt("covariance", "Covariance file", Spec(ValueType::Path), "data_formats.rst, covariance files"),
                             Defaulted(Opt("active", "Active", Words({"yes", "no"}), "an inactive data set is not used"), "yes")},
                            apparatus),
                     quantityFields),
              text),
       {},
       {}},
      {"normalization", {"data"}, "Shared normalization", "One normalization shared by the data sets that name it (F6).",
       {Req("value", "Value", Nuis(""), "the normalization nuisance")}, {}, {}},
      {"energy_scale", {"data"}, "Shared energy scale",
       "Energy-shift terms shared by the data sets that name it (F6); summed in the order offset, sqrt, function.",
       {Opt("offset", "Offset", Nuis("energy"), "E' = E + offset"),
        Opt("sqrt", "sqrt(E) term", Nuis("sqrt_energy"), "E' = E + b sqrt(E)"),
        Opt("function", "Function", Structured(kFunction), "E' = E + f(x), x = E_lab/MeV"),
        Repeatable(When("param", "Parameter", Structured(kParam), Condition::Present("function"), false,
                        "a parameter of the function"))},
       {{"offset", "sqrt", "function"}},
       {}},
      {"extrapolation", {"extrapolations"}, "Extrapolation", "Model values on a grid, without data (section 3.8, F20).",
       concat(concat(std::vector<FieldDef>{Req("reaction", "Reaction", Spec(ValueType::Reaction), "the reaction expression (section 3.5)"),
                      Req("quantity", "Quantity", QuantityEnum("extrapolation"), "what is computed (section 7.2)"),
                      Req("unit", "Unit", Spec(ValueType::Unit, "@quantity"), "unit of the output values"),
                      Req("energies", "Energies", GridOf("energy", false), "the energy grid (F10: E_k = lo + k step)"),
                      Opt("angles", "Angles", GridOf("angle", true), "the angle grid; required for quantities with an angle"),
                      Defaulted(Opt("active", "Active", Words({"yes", "no"}), "an inactive extrapolation is not computed"), "yes"),
                      ByQ("orders", "Legendre orders", Structured(kOrders), "the Legendre orders computed")},
                     apparatus),
              quantityFields),
       {},
       {{"file", "an extrapolation has no data file (F20)"},
        {"columns", "an extrapolation has no data file (F20)"},
        {"norm", "an extrapolation has no normalization (F20)"},
        {"shift", "an extrapolation has no energy shift (F20)"},
        {"covariance", "an extrapolation has no data (F20)"},
        {"energy_window", "the grid is the extrapolation's energy range (F20)"},
        {"angle_window", "the grid is the extrapolation's angle range (F20)"}}},
      {"target", {"targets"}, "Target",
       "The target the beam passes through: thickness, stopping, straggling, isotopes, Doppler, screening (F7).",
       {Req("kind", "Kind", Words({"layer", "infinite", "sample", "none"}),
            "layer: finite areal density; infinite: thick target; sample: neutron sample (atoms/b); none: isotopes only"),
        When("active", "Active nucleus", Spec(ValueType::Nucleus), Condition::Is("kind", {"layer", "infinite", "sample"}),
             true, "the nucleus the areal density and the stopping count"),
        Field("composition", "Composition", Spec(ValueType::Composition), Condition::Is("kind", {"layer", "infinite"}),
              Condition::Is("stopping", {"evaluation"}), "element stoichiometry; needed by 'stopping evaluation'"),
        When("areal_density", "Areal density", Nuis("areal_density"), Condition::Is("kind", {"layer", "sample"}), true,
             "areal density of the active nucleus"),
        Repeatable(When("epoch", "Epoch", Structured(kEpoch), Condition::Present("areal_density"), false,
                        "the areal density in one epoch (K5)")),
        When("stopping", "Stopping", Structured(kStopping), Condition::Is("kind", {"layer", "infinite"}), true,
             "stopping cross section"),
        When("straggling", "Straggling", Structured(kStraggling), Condition::Is("kind", {"layer", "infinite"}), true,
             "energy straggling in the target"),
        When("phase", "Phase", Words({"solid", "gas"}), Condition::Is("stopping", {"evaluation"}), true,
             "solid or gas stopping (with 'stopping evaluation')"),
        Repeatable(Opt("isotope", "Isotope", Structured(kIsotope), "abundance weight of a target nucleus")),
        When("isotope_sum", "Isotope sum", Structured(kIsotopeSum), Condition::Present("isotope"), false,
             "constraint on the sum of the abundances"),
        Opt("doppler", "Doppler", Structured(kDoppler), "free-gas Doppler broadening"),
        Opt("screening", "Electron screening", Nuis("energy"), "the screening potential U_e")},
       {},
       {}},
      {"beam", {"targets"}, "Beam", "The beam's energy distribution (F7).",
       {Opt("spread", "Energy spread", Structured(kSpread), "Gaussian energy spread"),
        Repeatable(Opt("profile", "Profile", Structured(kProfile), "a component of the beam energy distribution")),
        When("truncation", "Truncation", Num(), Condition::Present("profile"), true, "truncation of the profile"),
        When("detailed_balance", "Detailed balance", Words({"yes", "no"}), Condition::Present("profile"), true,
             "the profile kernel by detailed balance")},
       {{"spread", "profile"}},
       {}},
      {"detector", {"targets"}, "Detector", "Finite solid angle, energy resolution, time-of-flight resolution (F7).",
       {Opt("attenuation", "Attenuation Q_k", Spec(ValueType::NumberList), "Q_0 ... Q_k, lab"),
        Opt("energy_resolution", "Energy resolution", Structured(kResolution), "resolution of a reconstructed energy"),
        Opt("udr", "UDR", Structured(kUdr), "time-of-flight resolution")},
       {{"attenuation", "energy_resolution", "udr"}},
       {}},
      {"datum", {"modelData"}, "Model datum", "A literature value of a level property (section 3.10, M10).",
       concat(std::vector<FieldDef>{Req("quantity", "Quantity", QuantityEnum("datum"), "what the value is (section 7.2)"),
               ByQ("level", "Level", Int(0), "the level ID of <levels> (column 9)"),
               ByQ("J", "J^pi", Spec(ValueType::SpinParity), "must equal the level's J^pi (F13)"),
               ByQ("reaction", "Reaction", Spec(ValueType::Reaction), "entrance and exit"),
               ByQ("frame", "Frame", Words({"E_x", "E_cm", "E_lab"}), "frame of a resonance energy"),
               ByQ("pair", "Pair", Spec(ValueType::Pair), "the particle pair"),
               ByQ("l", "l", Int(0), "orbital angular momentum"),
               ByQ("s", "s", Structured(kChannelSpin), "channel spin, or coherent"),
               Req("value", "Value", WithUnit("@quantity"), "the value"),
               Req("error", "Error", WithUnit("@quantity", true), "its 1 sigma uncertainty"),
               ByQ("imag", "Imaginary part", WithUnit("@quantity"), "Im of a scattering length"),
               ByQ("imag_error", "Imaginary part error", WithUnit("@quantity", true), "its 1 sigma uncertainty")},
              text),
       {},
       {}},
  };
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
