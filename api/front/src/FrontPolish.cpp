// The facade's part of round 8 (docs/upgrade/qt_gui.md section 15): the fields a quantity uses (experiment's quantity
// table), the kernels a record of <targets> selects drawn from the engine's own kernels (the set-up call of
// docs/upgrade/physics_tab.md section 6.2: the engine's build of the project and its TargetEffect), a reaction
// expression split into its terms by the engine's parser, and the channels a level may have (the engine's
// experiment::AllowedChannels on the GUI's unsaved pairs).
#include "azure2/front/Engine.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

#include "Config.h"
#include "TargetEffect.h"
#include "azure2/experiment/Conditions.h"
#include "azure2/experiment/DataRows.h"
#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LevelEdits.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/LinearForm.h"
#include "azure2/experiment/Parsers.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Quantities.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/ValueText.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string Joined(const std::vector<std::string> &tokens) {
  std::string out;
  for (const std::string &t : tokens) {
    if (t.empty()) continue;
    if (!out.empty()) out += ' ';
    out += t;
  }
  return out;
}

/// The full width at half maximum of a curve sampled at x (any order of x, ascending or not), by linear
/// interpolation at the two outermost crossings of half the maximum; NaN when it does not fall below half on both
/// sides.
double FullWidthHalfMaximum(const std::vector<double> &x, const std::vector<double> &y) {
  if (x.size() < 3 || x.size() != y.size()) return kNaN;
  std::size_t top = 0;
  for (std::size_t i = 1; i < y.size(); i++)
    if (y[i] > y[top]) top = i;
  const double half = 0.5 * y[top];
  if (!(half > 0.0)) return kNaN;
  auto cross = [&](std::size_t a, std::size_t b) { return x[a] + (half - y[a]) * (x[b] - x[a]) / (y[b] - y[a]); };
  double left = kNaN, right = kNaN;
  for (std::size_t i = top; i > 0; i--)
    if (y[i - 1] < half) {
      left = cross(i - 1, i);
      break;
    }
  for (std::size_t i = top; i + 1 < y.size(); i++)
    if (y[i + 1] < half) {
      right = cross(i, i + 1);
      break;
    }
  return std::fabs(right - left);
}

/// The kernel's weight at each node of a uniform mesh, as a density: the weight over the node's share of the mesh
/// (h, h/2 at the two ends), the trapezoid rule's, so that the density's trapezoid sum is the weights' sum.
void Densities(const std::vector<double> &weights, double h, std::vector<double> &y) {
  y.assign(weights.size(), 0.0);
  for (std::size_t j = 0; j < weights.size(); j++)
    y[j] = weights[j] / ((j == 0 || j + 1 == weights.size()) ? 0.5 * h : h);
}

/// A curve of c.m. node energies (descending) and c.m. densities as the lab curve of the beam: x = E_cm / f
/// (ascending), y per lab MeV = f times the density per c.m. MeV.
void ToLab(const std::vector<double> &eCM, const std::vector<double> &densityCM, double f, KernelCurve &c) {
  c.x.clear();
  c.y.clear();
  for (std::size_t i = eCM.size(); i-- > 0;) {
    c.x.push_back(eCM[i] / f);
    c.y.push_back(densityCM[i] * f);
  }
  c.fwhm_keV = 1e3 * FullWidthHalfMaximum(c.x, c.y);
  c.fwhmCm_keV = f * c.fwhm_keV;
  c.xLabel = "E_lab (MeV)";
}

/// The first error diagnostic of a build, in words.
std::string FirstError(const ex::Format3Project &built) {
  for (const ex::Diagnostic &d : built.diagnostics)
    if (d.severity == ex::Severity::Error) return d.ToText();
  return "the project does not build";
}

/// The keys of the records of `kind` of `section` (1, 2, ... in file order), the key of a spec (DatasetSpec::key).
int KeyOf(const Project &project, const std::string &section, const std::string &kind, const std::string &name) {
  int key = 0;
  for (const Record &r : project.Records(section)) {
    if (r.kind != kind) continue;
    key++;
    if (r.name == name) return key;
  }
  return 0;
}

/// The energies of a data set's data file (lab MeV, the engine's conversion of each row): the column of its energy
/// role in the record's `columns`.  Empty when the file or the roles cannot be read.
std::vector<double> DataEnergies(const Project &project, const Record &dataset, const ex::ColumnConversion &conversion) {
  std::vector<double> out;
  ValueParts columns;
  std::string why;
  if (!ParseFieldValue(dataset.kind, "columns", dataset.Value("columns"), dataset.Value("quantity"), columns, why) ||
      columns.fromFile)
    return out;
  int energyColumn = -1;
  for (std::size_t i = 0; i < columns.columns.size(); i++)
    if (energyColumn < 0 && columns.columns[i].first.rfind("E_", 0) == 0) energyColumn = (int)i;
  ValueParts file;
  if (energyColumn < 0 || !ParseFieldValue(dataset.kind, "data_file", dataset.Value("data_file"), "", file, why)) return out;
  std::string path = file.text;
  if (!path.empty() && path[0] != '/') path = project.Directory() + "/" + path;
  std::ifstream in(path);
  for (std::string line; std::getline(in, line);) {
    std::istringstream row(line);
    std::vector<double> numbers;
    for (double v; row >> v;) numbers.push_back(v);
    if ((int)numbers.size() <= energyColumn) continue;
    double energy = numbers[energyColumn], angle = 0, value = 0, error = 0;
    std::vector<double> extra;
    ex::RowToEngine(conversion, energy, angle, value, error, extra);
    if (std::isfinite(energy)) out.push_back(energy);
  }
  return out;
}

/// One kernel of an effect made as a run makes it: an engine TargetEffect of the reduced spec, its kernels and the
/// parameters they read.  Its shape at c.m. energy eCM on `points` nodes.
bool Shape(const ex::TargetEffectSpec &spec, const std::string &directory, double cmFactor, int z1, int z2, double eCM,
           int points, KernelCurve &curve, std::string &why) {
  std::ostringstream notes;
  Config config(notes);
  config.baseDirectory = directory;
  const TargetEffect effect(spec, config);
  const ex::EffectKernels &kernels = effect.Kernels();
  const ex::KernelParameters p = effect.KernelParametersFor(config);
  const ex::Kernel *k = kernels.broadening ? kernels.broadening : kernels.resolution;
  if (!k) {
    why = "the engine selects no kernel for it";
    return false;
  }
  curve.kernel = k->name;
  if (k->name == "beam_profile") {
    // a lab distribution by itself: the profile's weight over its support
    double low = 0, high = 0;
    ex::BeamProfileSupport(p, low, high);
    if (!(high > low)) {
      why = "the profile has no width";
      return false;
    }
    for (int i = 0; i < points; i++) {
      const double e = low + (high - low) * i / (points - 1);
      curve.x.push_back(e);
      curve.y.push_back(ex::BeamProfileWeight(p, e));
    }
    double area = 0.0;
    for (std::size_t i = 1; i < curve.x.size(); i++)
      area += 0.5 * (curve.y[i] + curve.y[i - 1]) * (curve.x[i] - curve.x[i - 1]);
    if (area > 0.0)
      for (double &y : curve.y) y /= area;
    curve.xLabel = "E_lab (MeV)";
    curve.fwhm_keV = 1e3 * FullWidthHalfMaximum(curve.x, curve.y);
    curve.fwhmCm_keV = cmFactor * curve.fwhm_keV;
    return true;
  }
  ex::KernelPoint point;
  point.energy = eCM;
  point.labEnergy = point.nominalLabEnergy = eCM / cmFactor;
  point.cmFactor = cmFactor;
  point.z1 = z1;
  point.z2 = z2;
  point.binLow = point.binHigh = kNaN;
  point.temperature = spec.dopplerTemperature;
  double low = eCM, high = eCM;
  if (kernels.broadening) {
    ex::DopplerWindow(p, eCM, cmFactor, low, high);
  } else {
    ex::SupportPoint s;
    s.energy = eCM;
    s.labEnergy = eCM / cmFactor;
    s.cmFactor = cmFactor;
    s.z1 = z1;
    s.z2 = z2;
    s.binLow = s.binHigh = kNaN;
    ex::EffectKernels only;
    only.resolution = k;
    ex::EffectSupport(only, p, s, high, low);
  }
  if (!(high > low) || points < 3) {
    why = "the kernel has no width at this energy";
    return false;
  }
  // a uniform mesh, descending, a little wider than the support so that the curve reaches zero on both sides
  const double pad = 0.05 * (high - low);
  const double top = high + pad, bottom = std::max(low - pad, 0.0);
  const double h = (top - bottom) / (points - 1);
  std::vector<double> energy(points);
  for (int j = 0; j < points; j++) energy[j] = top - h * j;
  ex::KernelNodes nodes;
  nodes.n = points;
  nodes.energy = energy.data();
  std::vector<double> weights(points, 0.0);
  if (kernels.broadening) {
    std::vector<std::pair<int, double>> w;
    k->broadeningWeights(p, point, nodes, eCM, w, nullptr);
    for (const auto &t : w) weights[t.first - 1] += t.second;
  } else {
    std::vector<ex::LinearForm> unit;
    for (int j = 1; j <= points; j++) unit.push_back(ex::LinearForm::Unit(j));
    const ex::LinearForm form = k->weights(p, point, nodes, unit.data());
    for (const auto &t : form.Terms()) weights[t.first - 1] += t.second;
  }
  std::vector<double> density;
  Densities(weights, h, density);
  ToLab(energy, density, cmFactor, curve);
  return true;
}

}  // namespace

bool QuantityFieldRules(const std::string &quantity, std::vector<QuantityFieldRule> &out) {
  out.clear();
  const ex::Quantity *q = ex::FindQuantity(quantity);
  if (!q) return false;
  for (const ex::QuantityField &f : q->fields) {
    QuantityFieldRule r;
    r.field = f.field;
    r.required = f.required;
    r.conditional = f.when.kind != ex::Condition::Always || f.when.negate;
    if (r.conditional) r.when = ex::Describe(f.when);
    out.push_back(r);
  }
  return true;
}

bool KernelPreviewOf(const Project &project, const std::string &section, const Record &record, const std::string &user,
                     double labEnergy, int points, KernelPreview &out, std::string &why) {
  out = KernelPreview();
  if (points < 3) points = 3;
  Project copy = Project::FromText(project.Text(), project.Path());
  if (!copy.SetRecord(section, record, why)) return false;
  // the records that name it
  std::vector<KindInfo> kinds = RecordKinds();
  auto names = [&](const Record &r) {
    for (const KindInfo &k : kinds) {
      if (k.kind != r.kind) continue;
      for (const FieldInfo &f : k.fields) {
        if (f.value.type != "reference" ||
            std::find(f.value.values.begin(), f.value.values.end(), record.kind) == f.value.values.end())
          continue;
        const std::string v = r.Value(f.key);
        if (!v.empty() && v.substr(0, v.find('@')) == record.name) return true;
      }
    }
    return false;
  };
  struct User {
    std::string section;
    Record record;
  };
  std::vector<User> users;
  for (const char *s : {"data", "extrapolations"})
    for (const Record &r : copy.Records(s))
      if ((r.kind == "dataset" || r.kind == "extrapolation") && names(r)) {
        users.push_back({s, r});
        out.users.push_back(std::string(s) + ": " + r.name);
      }
  const User *chosen = nullptr;
  for (const User &u : users)
    if (!chosen && (user.empty() || user == u.section + ": " + u.record.name)) chosen = &u;
  if (!user.empty() && !chosen) {
    why = user + " does not name " + record.name;
    return false;
  }

  if (!chosen) {
    // no record names it: an extrapolation of the first particle pair's elastic scattering at E_0 that names it, in
    // the copy only, so that the engine builds the effect as for any user
    const io::Document doc = io::Document::Parse(copy.Text());
    bool rowsOk = false;
    const ex::LevelTables tables = ex::LevelTablesOfRows(ex::LevelRowsOf(doc, rowsOk));
    const ex::PairInfo *first = nullptr;
    for (const ex::PairInfo &p : tables.pairs)
      if (!first && p.kind == 0) first = &p;
    std::string field;
    for (const KindInfo &k : kinds)
      if (k.kind == "extrapolation")
        for (const FieldInfo &f : k.fields)
          if (field.empty() && f.value.type == "reference" &&
              std::find(f.value.values.begin(), f.value.values.end(), record.kind) != f.value.values.end())
            field = f.key;
    if (!first || field.empty()) {
      why = "no particle pair to draw the kernels for (the Particle Pairs tab)";
      return false;
    }
    const double e = std::isfinite(labEnergy) && labEnergy > 0.0 ? labEnergy : 1.0;
    Record x;
    x.kind = "extrapolation";
    x.name = "kernel_preview_of_" + record.name;
    x.fields = {{"reaction", ex::ReactionLabel(*first, first), ""},
                {"quantity", "sigma_nuclear", ""},  // elastic: the nuclear part, defined for charged particles too
                {"unit", "b", ""},
                {"energies", "list " + ex::ShortestNumberText(e) + " MeV lab", ""},
                {field, record.name, ""}};
    if (!copy.SetRecord("extrapolations", x, why)) return false;
    users.push_back({"extrapolations", x});
    chosen = &users.back();
    out.notes.push_back("no data set or extrapolation names " + record.name + ": drawn for " + ex::PairLabel(*first) +
                        " at the energy given");
  }
  const bool listed = !out.users.empty();
  // what is built: the record, the other records of <targets>, <levels>, the shared records of <data> (norms,
  // scales) and the one user; the other data sets, extrapolations and model data are left out of the copy, so that an error of theirs does not stop the drawing.  A
  // user without its quadrature (required once the record folds it) gets the smallest one: the kernels do not read it.
  for (const char *s : {"data", "extrapolations", "modelData"})
    for (const Record &r : copy.Records(s))
      if ((r.kind == "dataset" || r.kind == "extrapolation" || std::string(s) == "modelData") &&
          !(s == chosen->section && r.name == chosen->record.name)) {
        std::string ignored;
        copy.RemoveRecord(s, r.name, ignored);
      }
  {
    Record u = chosen->record;
    bool quadrature = false;
    for (const Field &f : u.fields) quadrature = quadrature || f.key == "quadrature";
    bool takes = false;
    for (const KindInfo &k : kinds)
      if (k.kind == u.kind)
        for (const FieldInfo &f : k.fields) takes = takes || f.key == "quadrature";
    if (takes && !quadrature) {
      u.fields.push_back({"quadrature", "points 1", ""});
      if (!copy.SetRecord(chosen->section, u, why)) return false;
    }
  }

  ex::Format3ReadOptions read;
  read.fileName = copy.Path();
  read.projectDir = copy.Directory();
  read.readDataFiles = false;
  read.dataFiles = ex::DataFiles::NotUsed;
  read.setUpChecks = false;
  read.allRecordsActive = true;
  const ex::Format3Project built = ex::ReadFormat3Text(copy.Text(), read);
  if (!built.Ok()) {
    why = FirstError(built);
    return false;
  }
  // the effect of the user, its entrance pair and its energies
  const ex::TargetEffectSpec *effect = nullptr;
  int entrance = 0;
  {
    if (listed) out.user = chosen->section + ": " + chosen->record.name;
    const bool data = chosen->section == "data";
    const int key = KeyOf(copy, chosen->section, chosen->record.kind, chosen->record.name);
    const auto &effects = data ? built.built.data.targetEffects.items : built.built.extrapolations.targetEffects.items;
    for (const auto &i : effects)
      if (!effect && std::find(i.spec.datasetKeys.begin(), i.spec.datasetKeys.end(), key) != i.spec.datasetKeys.end())
        effect = &i.spec;
    std::vector<double> energies;
    double minE = 0, maxE = 0;
    if (data) {
      for (const auto &i : built.built.data.datasets.items)
        if (i.spec.key == key) {
          entrance = i.spec.reaction.entranceKey;
          minE = i.spec.minE;
          maxE = i.spec.maxE;
          energies = DataEnergies(copy, chosen->record, i.spec.columns);
        }
    } else {
      for (const auto &i : built.built.extrapolations.extrapolations.items)
        if (i.spec.key == key) {
          entrance = i.spec.reaction.entranceKey;
          energies = i.spec.energies;
          minE = i.spec.minE;
          maxE = i.spec.maxE;
        }
    }
    std::vector<double> inside;
    for (double e : energies)
      if (e >= minE && e <= maxE) inside.push_back(e);
    if (!inside.empty()) {
      out.lowEnergy = *std::min_element(inside.begin(), inside.end());
      out.highEnergy = *std::max_element(inside.begin(), inside.end());
    } else if (maxE > minE && maxE < 1e300) {
      out.lowEnergy = minE;
      out.highEnergy = maxE;
    }
  }
  const ex::PairInfo *pair = nullptr;
  for (const ex::PairInfo &p : built.levels.pairs)
    if ((entrance > 0 && p.key == entrance) || (entrance <= 0 && !pair && p.kind == 0)) pair = &p;
  if (!pair || !(pair->lightMass > 0.0) || !(pair->heavyMass > 0.0)) {
    why = "no particle pair with masses to draw the kernels for (the Particle Pairs tab)";
    return false;
  }
  out.pair = ex::PairLabel(*pair);
  ex::Partition partition;
  partition.m1 = pair->lightMass;
  partition.m2 = pair->heavyMass;
  out.cmFactor = ex::BeamCMFactor(partition);
  if (!std::isfinite(labEnergy) || !(labEnergy > 0.0)) {
    labEnergy = out.highEnergy > 0.0 ? 0.5 * (out.lowEnergy + out.highEnergy) : 1.0;
    if (!(out.highEnergy > 0.0) && listed) out.notes.push_back("no energy of the data known: drawn at E_lab = 1 MeV");
  }
  out.labEnergy = labEnergy;
  out.cmEnergy = labEnergy * out.cmFactor;

  // the spec of the effect as the engine builds it; without a user, the record's own effect read from a project
  // in which one extrapolation names it would be needed: the record's fields are then taken from any effect that
  // carries this record's name, else refused
  if (!effect) {
    why = record.name + " selects no kernel the engine folds with";
    return false;
  }
  // one reduced spec per kernel the record selects: what of the effect comes from this record's kind
  std::vector<std::pair<std::string, ex::TargetEffectSpec>> parts;
  ex::TargetEffectSpec base;
  base.active = true;
  base.targetName = effect->targetName;
  base.datasetKeys = effect->datasetKeys;
  if (record.kind == "target" && effect->doppler) {
    ex::TargetEffectSpec d = base;
    d.doppler = true;
    d.dopplerTemperature = effect->dopplerTemperature;
    d.dopplerFit = effect->dopplerFit;
    d.dopplerVaryText = effect->dopplerVaryText;
    d.dopplerSigmaText = effect->dopplerSigmaText;
    d.dopplerMaxText = effect->dopplerMaxText;
    d.dopplerCentre = effect->dopplerCentre;
    parts.push_back({"free-gas Doppler broadening, T = " + ex::ShortestNumberText(effect->dopplerTemperature) + " K", d});
  }
  if (record.kind == "beam" && (effect->convolution || effect->convolutionCoefficients)) {
    ex::TargetEffectSpec g = base;
    g.convolution = effect->convolution;
    g.sigma = effect->sigma;
    g.convolutionCoefficients = effect->convolutionCoefficients;
    g.convolutionEquation = effect->convolutionEquation;
    g.convolutionCoefficientValues = effect->convolutionCoefficientValues;
    g.sigmaLab = effect->sigmaLab;
    parts.push_back({effect->convolutionCoefficients ? "Gaussian beam spread, sigma(E) a formula"
                                                     : "Gaussian beam spread, sigma = " +
                                                           ex::ShortestNumberText(1e3 * effect->sigma) + " keV lab",
                     g});
  }
  if (record.kind == "beam" && effect->beamProfile) {
    ex::TargetEffectSpec b = base;
    b.beamProfile = true;
    b.beamProfileComponents = effect->beamProfileComponents;
    b.beamTpcSigma = effect->beamTpcSigma;
    b.beamTruncation = effect->beamTruncation;
    b.beamDetailedBalance = effect->beamDetailedBalance;
    parts.push_back({"beam energy profile (lab, independent of E_0)", b});
  }
  if (record.kind == "detector" && !effect->kernel.empty()) {
    ex::TargetEffectSpec t = base;
    t.kernel = effect->kernel;
    t.kernelParameters = effect->kernelParameters;
    t.apparatusParameters = effect->apparatusParameters;
    std::string title = effect->kernel + " energy resolution";
    if (const ex::Kernel *k = ex::KernelNamed(effect->kernel))
      for (std::size_t i = 0; i < k->parameters.size() && i < effect->kernelParameters.size(); i++)
        title += ", " + k->parameters[i].name + " " +
                 (k->parameters[i].dimension == "energy" ? ex::ShortestNumberText(1e3 * effect->kernelParameters[i]) + " keV"
                                                         : ex::ShortestNumberText(effect->kernelParameters[i]) + " " +
                                                               k->parameters[i].unit);
    parts.push_back({title, t});
  }
  if (record.kind == "detector" && (effect->udr)) out.notes.push_back("the UDR kernel is not drawn");
  for (const auto &part : parts) {
    KernelCurve c;
    std::string reason;
    if (Shape(part.second, copy.Directory(), out.cmFactor, pair->lightZ, pair->heavyZ, out.cmEnergy, points, c, reason)) {
      c.title = part.first;
      out.curves.push_back(c);
    } else {
      out.notes.push_back(part.first + ": " + reason);
    }
  }
  if (out.curves.empty()) {
    why = out.notes.empty() ? record.name + " selects no kernel the preview draws" : out.notes.back();
    return false;
  }
  return true;
}

bool ParseReactionExpression(const Project &project, const std::string &text, ReactionExpression &out, std::string &why) {
  out = ReactionExpression();
  std::string comment;
  std::vector<io::RecordDiagnostic> problems;
  const std::vector<io::Token> tokens = io::Tokenize(text, 1, comment, problems);
  if (!problems.empty() || tokens.empty()) {
    why = tokens.empty() ? "a reaction is missing" : "the text does not read as tokens";
    return false;
  }
  std::vector<ex::Diagnostic> diagnostics;
  ex::ParseContext ctx;
  ctx.out = &diagnostics;
  ex::TokenCursor cur;
  cur.tokens = &tokens;
  cur.line = 1;
  cur.endCol = (int)text.size() + 1;
  ex::ReactionExpr expr;
  if (!ex::ParseReaction(cur, ctx, expr) || !cur.AtEnd()) {
    why = diagnostics.empty() ? "the reaction expression does not read" : diagnostics.front().message;
    const std::size_t at = why.find(": ");
    if (why.rfind("dataset", 0) == 0 && at != std::string::npos) why = why.substr(at + 2);
    return false;
  }
  out.op = expr.op == ex::ReactionExpr::Sum ? "sum" : expr.op == ex::ReactionExpr::Ratio ? "ratio" : "single";
  const io::Document doc = io::Document::Parse(project.Text());
  bool ok = false;
  const ex::LevelTables tables = ex::LevelTablesOfRows(ex::LevelRowsOf(doc, ok));
  for (const ex::ReactionTerm &t : expr.terms) {
    ReactionTermInfo info;
    if (t.hasScale) info.scale = t.scale.text;
    info.reaction = t.reaction.text;
    const std::vector<int> in = ex::PairsNamed(t.reaction.projectile, t.reaction.target, t.reaction.targetState, tables.pairs);
    info.entrance = in.size() == 1 ? in[0] : 0;
    if (t.reaction.state.kind == ex::ResidualState::Total) {
      info.exit = -1;
    } else {
      const std::vector<int> o = ex::PairsNamed(t.reaction.ejectile, t.reaction.residual, t.reaction.state, tables.pairs);
      info.exit = o.size() == 1 ? o[0] : 0;
    }
    if (t.hasAngle) {
      info.angle = t.angle.text;
      info.angleUnit = t.angleUnit ? t.angleUnit->token : "";
      info.angleFrame = t.angleFrame;
    }
    out.terms.push_back(info);
  }
  return true;
}

std::string ReactionExpressionText(const ReactionExpression &expression) {
  std::string out;
  const std::string join = expression.op == "ratio" ? " / " : " + ";
  for (std::size_t i = 0; i < expression.terms.size(); i++) {
    const ReactionTermInfo &t = expression.terms[i];
    if (i > 0) out += join;
    if (!t.scale.empty()) out += t.scale + "*";
    out += t.reaction;
    if (!t.angle.empty()) out += " @ " + Joined({t.angle, t.angleUnit, t.angleFrame});
  }
  return out;
}

bool AllowedChannels(const std::vector<LevelRow> &pairRows, int twoJ, int parity, int lMax, int multipolarityMax,
                     int multipolaritiesMax, const std::vector<int> &pairKeys, std::vector<AllowedChannelInfo> &out,
                     std::string &why) {
  out.clear();
  std::string text = "<levels>\n";
  for (const LevelRow &r : pairRows) text += Joined(r.tokens) + "\n";
  text += "</levels>\n";
  const io::Document doc = io::Document::Parse(text);
  ex::ChannelLimits limits;
  limits.lMax = lMax;
  limits.multipolarityMax = multipolarityMax;
  limits.multipolaritiesMax = multipolaritiesMax;
  const ex::AllowedChannelList list = ex::AllowedChannels(doc, twoJ, parity, limits, pairKeys);
  if (!list.error.empty()) {
    why = list.error;
    return false;
  }
  bool ok = false;
  const ex::LevelTables tables = ex::LevelTablesOfRows(ex::LevelRowsOf(doc, ok));
  for (const ex::LevelChannelInfo &c : list.channels) {
    const ex::PairInfo *p = ex::PairByKey(c.pairKey, tables.pairs);
    if (!p) continue;
    // identical particles: the symmetry of the pair's wave function leaves l + s even only (the engine's set-up
    // refuses the others for a pair with spin and drops them for a spin-0 pair)
    if (p->kind == 0 && p->identical && (c.l + c.twoS / 2) % 2 != 0) continue;
    // a photon channel couples J to J_f with a whole multipolarity: J and J_f both whole or both half-whole
    if (p->kind == 10 && (twoJ + p->twoJHeavy) % 2 != 0) continue;
    out.push_back({c.pairKey, c.l, c.twoS});
  }
  return true;
}

}  // namespace front
}  // namespace azure2
