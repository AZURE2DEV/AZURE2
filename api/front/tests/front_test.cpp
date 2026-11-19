// The front ends' facade without Qt (docs/upgrade/qt_gui.md section 10): front_test <source dir> <AZURE2 binary>.
//  (a) the tables: schema, commands, options, fitters, <config> keys;
//  (b) the test projects T1, T2, T3 and 13N opened and written back byte for byte, also after every record, every
//      <levels> row and every <config> entry is set to what it holds; an edit changes only its lines and is undone;
//      a record added, renamed and removed; --validate's diagnostics, a broken field named by its code;
//  (c) an old project refused on opening and converted by the converter of --convert;
//  (d) a run of `calculate` through the facade writes the same output files, byte for byte, as the program on the same
//      project, and returns the chi^2 the program prints; T1 and 13N;
//  (e) a fit stopped through the facade ends at its best point with exit 0 and fitStatus.out saying "stopped";
//  (f) physics: the rows of the physics table and one evaluated; a session of T1: poles, channel functions against
//      the row, the 1/2+ block's phase shift against the program's extrapolation, |U| = 1, |det A^-1| least at the
//      level; T2's THM vertex.
// Each check prints PASS or FAIL; exit 1 on a failure.
#include "azure2/front/Engine.h"

#include <map>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

std::string Sci(double x) {
  char b[32];
  std::snprintf(b, sizeof(b), "%.2e", x);
  return b;
}

std::string Bytes(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

/// A copy of a project's directory (the project file, its data and side files) under `to`.
fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  for (const char *d : {"output", "checks"}) {
    fs::remove_all(to / d);
    fs::create_directories(to / d);
  }
  return to / projectFile.filename();
}

class Collect : public front::RunObserver {
 public:
  void Message(const std::string &line, bool transient) override {
    std::lock_guard<std::mutex> lock(mutex);
    if (!transient) lines.push_back(line);
  }
  void Progress(const front::RunProgress &p) override {
    std::lock_guard<std::mutex> lock(mutex);
    progress++;
    if (stopAfter > 0 && progress >= stopAfter && stop) stop->store(true);
    last = p;
  }
  std::mutex mutex;
  std::vector<std::string> lines;
  long progress = 0;
  long stopAfter = 0;
  std::atomic<bool> *stop = nullptr;
  front::RunProgress last;
};

/// The output files of two runs are the same files with the same bytes.
bool SameOutputs(const fs::path &a, const fs::path &b, std::string &why) {
  std::vector<std::string> names;
  for (const auto &e : fs::directory_iterator(a))
    if (e.is_regular_file()) names.push_back(e.path().filename().string());
  int nb = 0;
  for (const auto &e : fs::directory_iterator(b))
    if (e.is_regular_file()) nb++;
  if ((int)names.size() != nb) {
    why = std::to_string(names.size()) + " files against " + std::to_string(nb);
    return false;
  }
  if (names.empty()) {
    why = "no output file";
    return false;
  }
  for (const std::string &n : names)
    if (!fs::exists(b / n) || Bytes(a / n) != Bytes(b / n)) {
      why = n + " differs";
      return false;
    }
  why = std::to_string(names.size()) + " files identical";
  return true;
}

void RoundTrip(const fs::path &file) {
  const std::string name = file.filename().string();
  front::Project p;
  std::string why;
  Check(front::Project::Open(file.string(), p, why), name + ": opens (" + why + ")");
  const std::string original = Bytes(file);
  Check(p.Text() == original && !p.Modified(), name + ": text as read");
  // every record, row and entry set to what it holds: not one byte moves
  for (const std::string &section : p.Sections())
    for (const front::Record &r : p.Records(section))
      if (!p.SetRecord(section, r, why)) Check(false, name + ": SetRecord " + r.name + ": " + why);
  if (!p.SetLevels(p.Levels(), why)) Check(false, name + ": SetLevels: " + why);
  for (const front::ConfigEntry &e : p.Config())
    if (e.name != "format" && !p.SetConfig(e.name, e.value, why)) Check(false, name + ": SetConfig " + e.name + ": " + why);
  Check(p.Text() == original, name + ": every record, row and entry set as it stands: the same bytes");
  const fs::path out = fs::path("front_test_work") / ("roundtrip_" + name);
  fs::create_directories(out.parent_path());
  Check(p.Save(out.string(), why) && Bytes(out) == original, name + ": saved byte for byte");
}


/// The value text with one blank between tokens outside quotes (what a writer of format 3 writes).
std::string OneBlank(const std::string &text) {
  std::string out;
  bool quote = false, blank = false;
  for (char c : text) {
    if (c == '"') quote = !quote;
    if (!quote && (c == ' ' || c == '\t')) {
      blank = true;
      continue;
    }
    if (blank && !out.empty()) out += ' ';
    blank = false;
    out += c;
  }
  return out;
}

/// (f) the forms of values (docs/upgrade/qt_gui.md section 11): every field of every record of the test projects
/// split into its pieces by the engine's parser and written back from them gives the same text; the summary
/// columns, the include switch and the groups come from the schema; pairs, reactions and the move of a record.
void FormsChecks(const fs::path &source) {
  const std::vector<fs::path> files = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                       source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr",
                                       source / "experiment/tests/data/section_examples.azr"};
  for (const fs::path &file : files) {
    front::Project p;
    std::string why;
    if (!front::Project::Open(file.string(), p, why)) {
      Check(false, file.filename().string() + ": opens");
      continue;
    }
    int fields = 0, same = 0;
    std::vector<std::string> differ;
    for (const auto &section : front::RecordSections())
      for (const front::Record &r : p.Records(section.first))
        for (const front::Field &f : r.fields) {
          front::ValueParts parts;
          fields++;
          if (!front::ParseFieldValue(r.kind, f.key, f.value, r.Value("quantity"), parts, why)) {
            differ.push_back(r.name + "." + f.key + " not read: " + why);
            continue;
          }
          const std::string written = front::FieldValueText(r.kind, f.key, parts);
          if (written == OneBlank(f.value))
            same++;
          else
            differ.push_back(r.name + "." + f.key + ": '" + f.value + "' -> '" + written + "'");
        }
    for (std::size_t i = 0; i < differ.size() && i < 8; i++) std::cout << "    " << differ[i] << std::endl;
    Check(fields > 0 && same == fields, file.filename().string() + ": " + std::to_string(same) + " of " +
                                            std::to_string(fields) + " field values split and written back unchanged");
  }
  // a piece changed: only it moves
  {
    front::ValueParts n;
    std::string why;
    const bool read = front::ParseFieldValue("target", "areal_density", "2e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17", "", n, why);
    Check(read && n.fit && n.hasPrior && n.priorPercent && n.priorWidth == "10" && n.gridMax == "3e+17" && n.unit == "atoms/cm2",
          "a fitted areal density with a prior and grid_max split into its pieces");
    n.fit = false;
    n.hasPrior = false;
    Check(front::FieldValueText("target", "areal_density", n) == "2e+17 atoms/cm2 fixed grid_max 3e+17",
          "the prior dropped and the mode fixed: '" + front::FieldValueText("target", "areal_density", n) + "'");
    front::ValueParts s;
    front::ParseFieldValue("target", "stopping", "formula \"a0*x^a1\" a0=6.0e-20 a1=-0.30 MeV*cm2/atom", "", s, why);
    Check(s.alternative == "formula" && s.Part("expression") && s.Part("expression")->expression == "a0*x^a1" &&
              s.Part("expression")->params.size() == 2 && s.Part("unit") && s.Part("unit")->text == "MeV*cm2/atom",
          "a stopping formula split into the alternative, the expression, its parameters and the unit");
    front::ValueParts b;
    front::ParseFieldValue("beam", "profile", "skewed_gaussian xi 0.34 MeV omega 0.035 MeV alpha 1.5 weight 1.0", "", b, why);
    Check(b.alternative == "skewed_gaussian" && b.Part("omega") && b.Part("omega")->number == "0.035" &&
              b.Part("alpha") && b.Part("alpha")->text == "1.5",
          "a beam profile split into its named parameters");
    const bool unitRefused = !front::ParseFieldValue("dataset", "unit", "barn_per_furlong", "sigma", b, why);
    Check(unitRefused && !why.empty() && why.find("field '") == std::string::npos, "a wrong unit refused in plain words: " + why);
    const bool modeRefused = !front::ParseFieldValue("target", "areal_density", "2e17 atoms/cm2", "", b, why);
    Check(modeRefused, "a nuisance without its mode refused: " + why);
    Check(front::ShortValueText("dataset", "norm", "1 fit prior 1 5%") == "1 (fit)", "the short text of a fitted norm");
    Check(front::ShortValueText("detector", "attenuation", "1 0.995  0.98") == "1 0.995 0.98", "the short text of a number list");
    // the beam profile's preview: the kernel's weight; a skewed Gaussian of alpha > 0 peaks above xi, and its
    // integral is the component's weight (the curve spans xi +- 5 omega)
    front::Record beam;
    beam.kind = "beam";
    beam.name = "b";
    beam.fields = {{"profile", "skewed_gaussian xi 340 keV omega 0.035 MeV alpha 1.5 weight 2", ""}, {"truncation", "0", ""}};
    front::ShapeCurve shape;
    const bool drawn = front::PreviewShape(beam, 401, shape, why);
    double area = 0, peak = 0, at = 0;
    for (std::size_t i = 1; i < shape.x.size(); i++) area += 0.5 * (shape.y[i] + shape.y[i - 1]) * (shape.x[i] - shape.x[i - 1]);
    for (std::size_t i = 0; i < shape.x.size(); i++)
      if (shape.y[i] > peak) {
        peak = shape.y[i];
        at = shape.x[i];
      }
    Check(drawn && std::fabs(area - 2.0) < 2e-3 && at > 0.34 && at < 0.34 + 0.035 && std::fabs(shape.x.front() - (0.34 - 0.175)) < 1e-12,
          "the beam profile's preview: area " + std::to_string(area) + " (weight 2), peak at " + std::to_string(at) + " MeV");
    front::Record target;
    target.kind = "target";
    Check(!front::PreviewShape(target, 10, shape, why), "no preview for a record without a shape (" + why + ")");
  }
  // the shapes from the schema
  {
    const front::KindInfo *dataset = nullptr, *target = nullptr, *beam = nullptr;
    const std::vector<front::KindInfo> kinds = front::RecordKinds();
    for (const front::KindInfo &k : kinds) {
      if (k.kind == "dataset") dataset = &k;
      if (k.kind == "target") target = &k;
      if (k.kind == "beam") beam = &k;
    }
    auto has = [](const front::KindInfo *k, const std::string &key) {
      return k && std::find(k->summaryFields.begin(), k->summaryFields.end(), key) != k->summaryFields.end();
    };
    std::string listed;
    if (dataset) for (const std::string &x : dataset->summaryFields) listed += x + " ";
    Check(dataset && dataset->summaryFields.size() >= 5 && dataset->summaryFields[0] == "reaction" && has(dataset, "norm") &&
              has(dataset, "target") && has(dataset, "energy_window") && !has(dataset, "columns") && !has(dataset, "source") && !has(dataset, "quadrature"),
          "a data set's summary columns: " + listed);
    Check(dataset && dataset->includeField == "active" && dataset->includeOn == "yes", "a data set's include switch is 'active'");
    Check(has(target, "stopping") && has(target, "areal_density") && has(beam, "profile"), "a target's and a beam's summary columns");
    std::string phase, truncation;
    if (target) for (const front::FieldInfo &f : target->fields) if (f.key == "phase") phase = f.group;
    if (beam) for (const front::FieldInfo &f : beam->fields) if (f.key == "truncation") truncation = f.group;
    Check(phase == "stopping" && truncation == "profile", "a stopping's phase grouped with the stopping, a profile's truncation with the profile");
    bool profileShape = false;
    if (beam)
      for (const front::FieldInfo &f : beam->fields)
        if (f.key == "profile")
          profileShape = f.value.alternatives.size() == 1 && f.value.alternatives[0].parts.size() == 4 &&
                         !f.value.alternatives[0].parts[0].value.units.empty();
    Check(profileShape, "the beam profile's shape: one alternative, four named parts with their units");
    const auto columns = front::DefaultColumns("dsigma_domega");
    Check(columns.size() == 4 && !columns[0].second.empty() && !columns[1].second.empty() && columns[3].second.empty(),
          "the columns of a new differential cross section: energy, angle, value, error");
    Check(!front::Quantities("dataset").empty() && !front::Quantities("dataset")[0].title.empty(), "the quantities with their titles");
  }
  // every grid kind the facade offers is one the engine reads
  {
    std::string why;
    int read = 0;
    const std::vector<front::GridKind> kinds = front::GridKinds(true);
    for (const front::GridKind &k : kinds) {
      front::ValueParts g;
      g.type = "grid";
      g.gridKind = k.word;
      g.lo = "0.1";
      g.hi = "1";
      g.unit = "MeV";
      g.frame = front::Frames()[1];
      g.step = "0.01";
      g.numbers = {"0.1", "0.2"};
      g.text = "grid.txt";
      g.tolerance = "1e-3";
      g.maxPoints = "100";
      g.startStep = "0.01";
      g.unit = k.unit ? "deg" : "";
      front::ValueParts back;
      const std::string angles = front::FieldValueText("extrapolation", "angles", g);
      read += front::ParseFieldValue("extrapolation", "angles", angles, "", back, why) && back.gridKind == k.word;
    }
    Check(read == (int)kinds.size(), std::to_string(read) + " of " + std::to_string(kinds.size()) + " grid kinds written and read back");
  }
  // pairs, reactions, the move of a record
  {
    front::Project p;
    std::string why;
    front::Project::Open((source / "tests/upgrade/T1/T1.azr").string(), p, why);
    const std::vector<front::PairChoice> pairs = front::Pairs(p);
    Check(!pairs.empty() && pairs[0].key == 1 && !pairs[0].label.empty(), "T1: the pairs of <levels> (" + (pairs.empty() ? std::string() : pairs[0].label) + ")");
    const std::string reaction = p.Records("data")[0].Value("reaction");
    int in = 0, out = 0;
    const bool resolved = front::ReactionPairs(p, reaction, in, out);
    Check(resolved && in > 0 && out > 0 && front::ReactionLabel(p, in, out) == reaction,
          "T1: the reaction " + reaction + " is the pairs " + std::to_string(in) + " -> " + std::to_string(out) + " and back");
    const std::string original = p.Text();
    const std::vector<front::Record> before = p.Records("data");
    Check(p.MoveRecord("data", before[1].name, -1, why) && p.Records("data")[0].name == before[1].name &&
              p.Records("data")[1].name == before[0].name,
          "T1: a data set moved up");
    Check(p.MoveRecord("data", before[1].name, +1, why) && p.Text() == original, "T1: and down again: the same bytes");
    const bool refused = !p.MoveRecord("data", before[0].name, -1, why);
    Check(refused, "T1: the first record cannot move up (" + why + ")");
  }
}

}  // namespace

/// Round 8 (polish, docs/upgrade/qt_gui.md section 15): the fields a quantity uses, the kernels drawn from the
/// engine's own (against the closed forms of experiment/tests/kernels_test.cpp), a reaction expression's terms, and
/// the channel rule on unsaved pairs.
void PolishChecks(const fs::path &source) {
  std::string why;
  // (a) the fields a quantity uses
  std::vector<front::QuantityFieldRule> rules;
  bool phase = front::QuantityFieldRules("phase_shift", rules);
  int required = 0;
  for (const front::QuantityFieldRule &r : rules) required += r.required && (r.field == "J" || r.field == "l" || r.field == "s");
  Check(phase && required == 3, "polish: phase_shift requires J, l and s");
  bool sigmaJ = false, multipole = false;
  front::QuantityFieldRules("sigma", rules);
  for (const front::QuantityFieldRule &r : rules) {
    sigmaJ = sigmaJ || r.field == "J";
    if (r.field == "multipole") multipole = r.conditional && !r.when.empty();
  }
  Check(!sigmaJ && multipole, "polish: sigma does not use J; its multipole only under a condition in words");
  Check(!front::QuantityFieldRules("no_such_quantity", rules), "polish: an unknown quantity is refused");

  // (b) kernels: T1 (p + 12C) with a Doppler target, a Gaussian spread and two table resolutions, each named by an
  // extrapolation at E_lab 1.2 MeV
  front::Project p;
  Check(front::Project::Open((source / "tests/upgrade/T1/T1.azr").string(), p, why), "polish: T1 opens");
  auto record = [](const std::string &kind, const std::string &name, std::vector<front::Field> fields) {
    front::Record r;
    r.kind = kind;
    r.name = name;
    r.fields = std::move(fields);
    return r;
  };
  const front::Record td = record("target", "TD", {{"kind", "none", ""}, {"doppler", "free_gas 300 K fixed", ""}});
  const front::Record bg = record("beam", "BG", {{"spread", "1 keV lab", ""}});
  const front::Record dr = record("detector", "DR", {{"energy_resolution", "rectangular width 20 keV fixed", ""}});
  const front::Record dt = record("detector", "DT", {{"energy_resolution", "triangular width 20 keV fixed", ""}});
  bool set = p.SetRecord("targets", td, why) && p.SetRecord("targets", bg, why) && p.SetRecord("targets", dr, why) &&
             p.SetRecord("targets", dt, why);
  const std::vector<std::pair<std::string, std::string>> users = {{"xd", "target TD"}, {"xg", "beam BG"},
                                                                 {"xr", "detector DR"}, {"xt", "detector DT"}};
  for (const auto &u : users) {
    const std::string kind = u.second.substr(0, u.second.find(' ')), name = u.second.substr(u.second.find(' ') + 1);
    set = set && p.SetRecord("extrapolations",
                             record("extrapolation", u.first,
                                    {{"reaction", "12C(p,g)13N", ""}, {"quantity", "sigma", ""}, {"unit", "b", ""},
                                     {"energies", "list 1.1 1.2 1.3 MeV lab", ""}, {kind, name, ""},
                                     {"quadrature", "points 20", ""}}),
                             why);
  }
  Check(set, "polish: kernel records and their extrapolations set (" + why + ")");
  auto area = [](const front::KernelCurve &c) {
    double a = 0;
    for (std::size_t i = 1; i < c.x.size(); i++) a += 0.5 * (c.y[i] + c.y[i - 1]) * (c.x[i] - c.x[i - 1]);
    return a;
  };
  const double gaussianFwhm = 2.0 * std::sqrt(2.0 * std::log(2.0));
  front::KernelPreview k;
  why.clear();
  bool ok = front::KernelPreviewOf(p, "targets", dr, "", std::nan(""), 801, k, why);
  Check(ok && k.user == "extrapolations: xr" && std::fabs(k.labEnergy - 1.2) < 1e-12 && k.curves.size() == 1 &&
            k.curves[0].kernel == "rectangular" && std::fabs(k.curves[0].fwhm_keV - 20.0) < 0.05 &&
            std::fabs(area(k.curves[0]) - 1.0) < 1e-6,
        "polish: rectangular width 20 keV: FWHM " + (ok && !k.curves.empty() ? std::to_string(k.curves[0].fwhm_keV) : why) +
            " keV, area 1, at the middle of the grid (1.2 MeV)");
  why.clear();
  ok = front::KernelPreviewOf(p, "targets", dt, "", 1.2, 801, k, why);
  Check(ok && k.curves.size() == 1 && std::fabs(k.curves[0].fwhm_keV - 10.0) < 0.05 && std::fabs(area(k.curves[0]) - 1.0) < 1e-6,
        "polish: triangular base 20 keV: FWHM " + (ok && !k.curves.empty() ? std::to_string(k.curves[0].fwhm_keV) : why) + " keV");
  why.clear();
  ok = front::KernelPreviewOf(p, "targets", bg, "", 1.2, 801, k, why);
  Check(ok && k.curves.size() == 1 && k.curves[0].kernel == "gaussian" &&
            std::fabs(k.curves[0].fwhm_keV - gaussianFwhm) < 2e-3 * gaussianFwhm && std::fabs(area(k.curves[0]) - 1.0) < 1e-6,
        "polish: Gaussian sigma 1 keV lab: FWHM " + (ok && !k.curves.empty() ? std::to_string(k.curves[0].fwhm_keV) : why) +
            " keV (2 sqrt(2 ln 2) = " + std::to_string(gaussianFwhm) + ")");
  why.clear();
  ok = front::KernelPreviewOf(p, "targets", td, "", 1.2, 801, k, why);
  // the free gas: a Gaussian in velocity of standard deviation u / sqrt 2, u^2 = kT m/M; in the lab-equivalent energy
  // E = w^2 its standard deviation is sqrt(2 E kT m/M) (E >> kT), m/M = (1 - f)/f
  const double kT = 8.617333262e-11 * 300.0, f = k.cmFactor;
  const double doppler = 1e3 * gaussianFwhm * std::sqrt(2.0 * 1.2 * kT * (1.0 - f) / f);
  Check(ok && k.curves.size() == 1 && k.curves[0].kernel == "doppler" && std::fabs(k.curves[0].fwhm_keV - doppler) < 5e-3 * doppler &&
            std::fabs(area(k.curves[0]) - 1.0) < 1e-3,
        "polish: free-gas Doppler at 300 K, E_lab 1.2 MeV, p + 12C: FWHM " +
            (ok && !k.curves.empty() ? std::to_string(k.curves[0].fwhm_keV) : why) + " keV (closed form " +
            std::to_string(doppler) + " keV)");
  // the record as the form holds it, not as the project: 600 K widens it by sqrt 2
  front::Record hot = td;
  hot.fields[1].value = "free_gas 600 K fixed";
  why.clear();
  front::KernelPreview k2;
  ok = front::KernelPreviewOf(p, "targets", hot, "", 1.2, 801, k2, why) && !k.curves.empty() && !k2.curves.empty();
  Check(ok && std::fabs(k2.curves[0].fwhm_keV / k.curves[0].fwhm_keV - std::sqrt(2.0)) < 5e-3,
        "polish: the form's unapplied 600 K drawn: the width times sqrt 2");
  why.clear();
  const front::Record lonely = record("beam", "BL", {{"spread", "1 keV lab", ""}});
  ok = front::KernelPreviewOf(p, "targets", lonely, "", 1.2, 801, k, why);
  Check(ok && k.users.empty() && k.user.empty() && !k.notes.empty() && k.curves.size() == 1 &&
            std::fabs(k.curves[0].fwhm_keV - gaussianFwhm) < 2e-3 * gaussianFwhm,
        "polish: a record no data set names is drawn for the first pair, and says so: " +
            (k.notes.empty() ? why : k.notes[0]));
  why.clear();
  const front::Record plain = record("target", "TP", {{"kind", "none", ""}});
  Check(!front::KernelPreviewOf(p, "targets", plain, "", 1.2, 101, k, why) && !why.empty(),
        "polish: a target without a kernel: nothing drawn (" + why + ")");

  // the columns of a data file: the quantity's first ones, then columns the engine skips, or fewer
  const auto four = front::ColumnsForFile("sigma", 4), two = front::ColumnsForFile("sigma", 2);
  Check(four.size() == 4 && four[0] == front::DefaultColumns("sigma")[0] && four[3].first == "_" && four[3].second.empty() &&
            two.size() == 2,
        "polish: the columns of a 4-column file of sigma: " + four[0].first + " ... " + four[3].first);

  // (c) a reaction expression's terms
  front::Project t2;
  Check(front::Project::Open((source / "tests/upgrade/T2/T2.azr").string(), t2, why), "polish: T2 opens");
  front::ReactionExpression e;
  const std::string sum = "0.5*15N(p,a)12C + 15N(p,p)15N";
  ok = front::ParseReactionExpression(t2, sum, e, why);
  Check(ok && e.op == "sum" && e.terms.size() == 2 && e.terms[0].scale == "0.5" && e.terms[0].entrance > 0 &&
            e.terms[0].exit > 0 && e.terms[1].scale.empty() && front::ReactionExpressionText(e) == sum,
        "polish: a sum's terms with their weight and pairs, written back the same: " + front::ReactionExpressionText(e));
  const std::string ratio = "15N(p,a)12C @ 90 deg cm / 15N(p,p)15N @ 90 deg cm";
  ok = front::ParseReactionExpression(t2, ratio, e, why);
  Check(ok && e.op == "ratio" && e.terms.size() == 2 && e.terms[1].angle == "90" && e.terms[1].angleFrame == "cm" &&
            front::ReactionExpressionText(e) == ratio,
        "polish: a ratio at a fixed angle written back the same");
  why.clear();
  Check(!front::ParseReactionExpression(t2, "15N(p,a)12C +", e, why) && !why.empty(), "polish: a sum without its term refused: " + why);

  // (d) the channel rule on the rows of T1, and an identical pair made of its first row
  const std::vector<front::LevelRow> rows = p.Levels();
  std::vector<front::AllowedChannelInfo> channels;
  ok = front::AllowedChannels(rows, 1, +1, 2, 2, 2, {}, channels, why);
  std::string list;
  for (const front::AllowedChannelInfo &c : channels)
    list += " (" + std::to_string(c.pairKey) + "," + std::to_string(c.l) + "," + std::to_string(c.twoS) + ")";
  Check(ok && channels.size() == 4 && channels[0].pairKey == 1 && channels[0].l == 0 && channels[0].twoS == 1,
        "polish: T1, 1/2+: p + 12C s-wave and the photon channels:" + list);
  Check(front::AllowedChannels(rows, 2, +1, 2, 2, 2, {}, channels, why) && channels.empty(),
        "polish: T1, 1+: no channel (p + 12C has half-whole J only; no photon from J = 1 to a J_f = 1/2 or 3/2)");
  front::LevelRow alpha;
  for (const front::LevelRow &r : rows)
    if (alpha.tokens.empty() && !r.tokens.empty()) alpha = r;
  for (const char *c : {"z1", "z2"}) alpha.tokens[front::LevelColumnIndex(c)] = "2";
  for (const char *c : {"m1", "m2"}) alpha.tokens[front::LevelColumnIndex(c)] = "4.00260325";
  for (const char *c : {"j1", "j2"}) alpha.tokens[front::LevelColumnIndex(c)] = "0";
  for (const char *c : {"pi1", "pi2"}) alpha.tokens[front::LevelColumnIndex(c)] = "1";
  for (const char *c : {"g1", "g2", "e2", "e1"})
    if (front::LevelColumnIndex(c) < (int)alpha.tokens.size()) alpha.tokens[front::LevelColumnIndex(c)] = "0";
  alpha.tokens[front::LevelColumnIndex("s2")] = "0";
  why.clear();
  ok = front::AllowedChannels({alpha}, 2, -1, 2, 2, 2, {}, channels, why);
  Check(ok && channels.empty(), "polish: alpha + alpha, 1-: l = 1 is left out (identical bosons: l + s even) " + why);
  ok = front::AllowedChannels({alpha}, 4, +1, 2, 2, 2, {}, channels, why);
  Check(ok && channels.size() == 1 && channels[0].l == 2, "polish: alpha + alpha, 2+: l = 2");
  why.clear();
  Check(!front::AllowedChannels(rows, 1, 0, 2, 2, 2, {}, channels, why) && !why.empty(), "polish: parity 0 refused: " + why);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: front_test <source dir> <AZURE2 binary>" << std::endl;
    return 2;
  }
  const fs::path source = argv[1];
  const std::string binary = argv[2];
  fs::create_directories("front_test_work");

  // (a) the tables
  Check(front::SchemaJson().find("\"dataset\"") != std::string::npos, "schema JSON names the record kind dataset");
  bool calculate = false, fit = false;
  for (const front::Command &c : front::Commands()) {
    calculate = calculate || c.name == "calculate";
    fit = fit || c.name == "fit";
  }
  Check(calculate && fit, "the commands calculate and fit");
  Check(!front::Options().empty() && !front::Fitters().empty(), "the options and the fitters");
  Check(front::ConfigTableJson().find("\"use_brune\"") != std::string::npos, "the <config> table names use_brune");
  {
    // the options' editors (the GUI's Calculate tab): the unit the row's help names, when it applies, its other way
    std::map<std::string, front::Option> byName;
    for (const front::Option &o : front::Options()) byName[o.option] = o;
    Check(byName.count("--temperature-min=") && byName["--temperature-min="].unit == "GK" &&
              byName["--temperature-min="].instead == "--temperature-file" &&
              byName["--temperature-min="].when == "without --temperature-file",
          "--temperature-min: unit GK, the other way --temperature-file, asked without it: unit '" +
              byName["--temperature-min="].unit + "', instead '" + byName["--temperature-min="].instead + "'");
    Check(byName["--energy-spread="].unit == "keV" && byName["--walkers="].unit.empty() &&
              byName["--covariance-band"].unit.empty(),
          "units: --energy-spread keV, none for --walkers and a flag ('" + byName["--walkers="].unit + "')");
    Check(byName["--scale-covariance"].needs == "--covariance-band" && byName["--covariance-band"].needs.empty() &&
              byName["--scale-covariance"].when == "with --covariance-band",
          "--scale-covariance applies with --covariance-band");
    Check(byName["--delta-chi2="].lowOpen && byName["--delta-chi2="].low == 0 && !byName["--walkers="].lowOpen,
          "--delta-chi2: above 0 (an open range), --walkers at least its low");
  }
  {
    bool reaction = false;
    for (const front::KindInfo &k : front::RecordKinds())
      if (k.kind == "dataset" && !k.fields.empty() && k.fields[0].key == "reaction") reaction = k.fields[0].required;
    Check(reaction, "the record kinds in the schema's order: dataset's first field the required reaction");
    const std::vector<std::string> units = front::QuantityUnits("sfactor");
    Check(!units.empty() && units[0] == "MeV*b", "the units of the S factor");
  }

  // (b) round trips and edits
  const std::vector<fs::path> projects = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                          source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr"};
  for (const fs::path &f : projects) RoundTrip(f);
  {
    front::Project p;
    std::string why;
    front::Project::Open(projects[0].string(), p, why);
    const std::string original = p.Text();
    std::vector<front::Record> data = p.Records("data");
    Check(data.size() > 3 && data[0].kind == "dataset", "T1: the data sets");
    front::Record r = data[0];
    for (front::Field &f : r.fields)
      if (f.key == "norm") f.value = "1.05 fit";
    Check(p.SetRecord("data", r, why), "T1: a norm set to 1.05 fit");
    // exactly one line differs
    std::istringstream a(original), b(p.Text());
    std::string la, lb;
    int differ = 0, lines = 0;
    while (std::getline(a, la) && std::getline(b, lb)) {
      lines++;
      if (la != lb) differ++;
    }
    Check(differ == 1, "T1: the edit changed one line of " + std::to_string(lines));
    Check(p.Records("data")[0].Value("norm") == "1.05 fit", "T1: the record reads back the new value");
    r.fields.push_back({"energy_shift", "0 keV fixed", ""});
    Check(p.SetRecord("data", r, why) && p.Records("data")[0].Value("energy_shift") == "0 keV fixed",
          "T1: a field added (" + why + ")");
    front::Record copy = data[1];
    copy.name = "copy_of_" + data[1].name;
    Check(p.SetRecord("data", copy, why) && p.Records("data").size() == data.size() + 1, "T1: a record appended");
    Check(p.RenameRecord("data", copy.name, "renamed", why) && p.Records("data").back().name == "renamed",
          "T1: a record renamed");
    Check(p.RemoveRecord("data", "renamed", why) && p.Records("data").size() == data.size(), "T1: a record removed");
    Check(p.SetRecord("data", data[0], why) && p.Text() == original, "T1: the first record set back: the same bytes");
    front::Validation v = p.Validate();
    Check(v.ok && v.format3, "T1: valid (" + std::to_string(v.diagnostics.size()) + " diagnostics)");
    // a broken field: the engine's code, the record and the field
    front::Record bad = data[0];
    for (front::Field &f : bad.fields)
      if (f.key == "quantity") f.value = "sigmaa";
    p.SetRecord("data", bad, why);
    v = p.Validate();
    bool named = false;
    for (const front::Diagnostic &d : v.diagnostics)
      named = named || (d.severity == "error" && d.record == bad.name && d.field == "quantity" && !d.code.empty());
    Check(!v.ok && named, "T1: quantity sigmaa refused with a code naming the record and the field");
    // <levels>: one width changed through the row writer
    p = front::Project::FromText(original, projects[0].string());
    std::vector<front::LevelRow> rows = p.Levels();
    Check(!rows.empty() && rows[0].tokens.size() >= 31, "T1: the <levels> rows");
    rows[0].tokens[11] = "1.9";
    Check(p.SetLevels(rows, why) && p.Levels()[0].tokens[11] == "1.9", "T1: a level width written by the row writer");
    rows[0].tokens[11] = "x";
    Check(!p.SetLevels(rows, why), "T1: a row the engine cannot read is refused (" + why + ")");
    // <config>
    Check(p.SetConfig("use_brune", "0", why) && p.ConfigValue("use_brune") == "0", "T1: <config> use_brune=0");
    Check(!p.SetConfig("no_such_key", "1", why), "T1: an unknown <config> key refused");
  }
  {
    // an empty project: what New makes, valid once saved
    front::Project e = front::Project::Empty((fs::path("front_test_work") / "empty.azr").string());
    std::string why;
    Check(e.Save("", why), "an empty project saved");
    front::Project back;
    Check(front::Project::Open(e.Path(), back, why), "the empty project opens (" + why + ")");
    Check(back.ConfigValue("format") == "3", "the empty project is format 3");
  }

  // (c) an old project: refused, then converted
  {
    const fs::path old = source / "tests/convert/fixtures/tests/13N/13N.azr";
    front::Project p;
    std::string why;
    Check(!front::Project::Open(old.string(), p, why) && why.find("--convert") != std::string::npos,
          "an old project refused with the conversion command");
    const fs::path dir = fs::path("front_test_work") / "convert";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const fs::path to = dir / "13N.f3.azr";
    std::string report;
    Check(front::ConvertProject(old.string(), to.string(), report) && front::Project::Open(to.string(), p, why),
          "the old project converted and opened");
    Check(!front::ConvertProject(old.string(), to.string(), report), "the converter overwrites nothing");
  }

  // (d) calculate through the facade = the program
  for (const fs::path &f : {projects[0], projects[3]}) {
    const std::string name = f.stem().string();
    const fs::path a = CopyProject(f, fs::path("front_test_work") / ("cli_" + name));
    const fs::path b = CopyProject(f, fs::path("front_test_work") / ("front_" + name));
    std::string command = "\"" + binary + "\" calculate --threads=2 \"" + a.string() + "\" > \"" +
                          (a.parent_path() / "stdout.txt").string() + "\" 2>&1";
#ifdef _WIN32
    // cmd.exe drops the first and the last quote of a line that starts with one: an outer pair keeps the inner ones.
    command = "\"" + command + "\"";
#endif
    const int status = std::system(command.c_str());
    Collect observer;
    front::RunRequest request;
    request.command = "calculate";
    request.project = b.string();
    request.threads = 2;
    const front::RunResult result = front::Run(request, observer, nullptr);
    Check(status == 0 && result.exitCode == 0, name + ": both runs end with exit 0");
    std::string why;
    const bool same = SameOutputs(a.parent_path() / "output", b.parent_path() / "output", why);
    Check(same, name + ": output files of the program and of the facade: " + why);
    Check(!result.outputFiles.empty(), name + ": the run lists " + std::to_string(result.outputFiles.size()) + " files");
    // the chi^2 as data equals the program's line
    const std::string out = Bytes(a.parent_path() / "stdout.txt");
    const std::size_t at = out.find("Total Chi-Squared: ");
    const double printed = at == std::string::npos ? NAN : std::atof(out.c_str() + at + 19);
    Check(result.haveChiSquared && std::fabs(result.totalChiSquared - printed) <= 1e-11 * std::fabs(printed),
          name + ": total chi^2 as data " + std::to_string(result.totalChiSquared) + " = printed " + std::to_string(printed));
    Check(!result.dataSets.empty(), name + ": " + std::to_string(result.dataSets.size()) + " data sets' chi^2");
    bool sawTotal = false;
    for (const std::string &l : observer.lines) sawTotal = sawTotal || l.rfind("Total Chi-Squared: ", 0) == 0;
    Check(sawTotal, name + ": the messages carry the program's lines");

    // (f) the run's curves as data: every row of every AZUREOut_*.out file in one curve, each data set's points
    // as many as chiSquared.out's N
    front::Project opened;
    front::Project::Open(b.string(), opened, why);
    const front::Curves curves = front::ReadCurves(opened);
    std::size_t rows = 0, points = 0;
    for (const auto &e : fs::directory_iterator(b.parent_path() / "output")) {
      const std::string file = e.path().filename().string();
      if (file.rfind("AZUREOut_aa=", 0) != 0 || e.path().extension() != ".out") continue;
      std::ifstream in(e.path());
      std::string line;
      while (std::getline(in, line)) rows += line.find_first_not_of(" \t") != std::string::npos;
    }
    for (const front::Curve &c : curves.curves) points += c.points.size() + c.coefficients.size();
    std::size_t chiN = 0;
    for (const front::DataSetChiSquared &d : result.dataSets) chiN += d.points;
    std::string notes;
    for (const std::string &n : curves.notes) notes += " " + n;
    Check(!curves.curves.empty() && points == rows && points == chiN,
          name + ": the curves hold " + std::to_string(points) + " points, the files " + std::to_string(rows) +
              " rows, chiSquared.out " + std::to_string(chiN) + notes);
    bool named = true;
    for (const front::Curve &c : curves.curves) named = named && !c.record.empty() && c.section == "data";
    Check(named, name + ": every curve names its data set");
  }

  // (g) the fit's parameters, <parameterSettings> and <mcmc> rows through their grammars
  {
    const fs::path b = CopyProject(projects[0], fs::path("front_test_work") / "params_T1");
    std::vector<front::Parameter> params;
    std::string why;
    const bool read = front::ReadParameters(b.string(), params, why);
    Check(read && !params.empty(),
          "T1: " + std::to_string(params.size()) + " parameters of a session (" + why + ")");
    bool energy = false;
    for (const front::Parameter &p : params) energy = energy || (p.name == "energy_1" && p.kind == "energy");
    Check(energy, "T1: energy_1 among them, of kind energy");
    front::Project p;
    front::Project::Open(b.string(), p, why);
    const std::string original = p.Text();
    front::ParameterSetting row;
    row.name = "energy_1";
    row.nominal = 0.5;
    row.lower = 0.4;
    row.upper = 0.6;
    Check(front::SetParameterSetting(p, row, why) && front::ParameterSettings(p).size() == 1 &&
              front::ParameterSettings(p)[0].upper == 0.6,
          "T1: a <parameterSettings> row written and read back");
    Check(p.Validate().ok, "T1: valid with the row");
    Check(front::SetParameterSetting(p, row, why, true) && front::ParameterSettings(p).empty(), "T1: the row removed");
    front::SamplerPrior prior{"energy_1", 0.5, 0.5, 0.01, true};
    front::SetSamplerPriors(p, {prior});
    const std::vector<front::SamplerPrior> back = front::SamplerPriors(p);
    Check(back.size() == 1 && back[0].name == "energy_1" && back[0].sigma == 0.01 && back[0].use,
          "T1: an <mcmc> prior row written and read back");
    front::SetSamplerPriors(p, {});
    Check(front::SamplerPriors(p).empty(), "T1: the <mcmc> prior rows cleared");
    std::vector<front::ParameterValue> values;
    Check(!front::ReadParameterFile((b.parent_path() / "none.sav").string(), values, why), "a missing parameter file");
    (void)original;
  }

  // (e) a fit stopped through the facade
  {
    const fs::path b = CopyProject(source / "tests/upgrade/T1/T1_start.azr", fs::path("front_test_work") / "stop");
    std::atomic<bool> stop{false};
    Collect observer;
    observer.stop = &stop;
    observer.stopAfter = 25;
    front::RunRequest request;
    request.command = "fit";
    request.project = b.string();
    request.threads = 2;
    const front::RunResult result = front::Run(request, observer, &stop);
    const std::string status = Bytes(b.parent_path() / "output" / "fitStatus.out");
    Check(result.stopped && result.exitCode == 0, "T1_start: the stopped fit ends with exit 0");
    Check(status.find("stopped") != std::string::npos, "T1_start: fitStatus.out says stopped");
    Check(observer.progress >= 25 && observer.progress < 2000, "T1_start: " + std::to_string(observer.progress) +
                                                                   " progress records, the fit ended soon after the request");
    Check(result.haveChiSquared, "T1_start: the final evaluation at the best point");
  }

  // (i) the sampler's chain read back by the engine's reader (a samples.mcmc written here: 2 steps, 2 walkers, 2
  // parameters, rows out of order), a missing file and a broken one refused
  {
    const fs::path b = CopyProject(projects[0], fs::path("front_test_work") / "chain_T1");
    front::Project p;
    std::string why;
    front::Project::Open(b.string(), p, why);
    front::SamplerChain chain;
    const bool none = !front::ReadSamplerChain(p, chain, why) && why.find("samples.mcmc") != std::string::npos;
    Check(none, "chain: none before a sampler's run (" + why + ")");
    why.clear();
    {
      std::ofstream f(b.parent_path() / "output" / "samples.mcmc");
      f << "step,walker,logprob,loglikelihood,logprior,param0,param1\n"
        << "0,0,-1.5,-1,-0.5,1.25,10\n0,1,-2.5,-2,-0.5,1.5,11\n1,1,-4.5,-4,-0.5,2,13\n1,0,-3.5,-3,-0.5,1.75,12\n";
    }
    const bool read = front::ReadSamplerChain(p, chain, why);
    Check(read && chain.steps == 2 && chain.walkers == 2 && chain.parameters == 2 && chain.chain.size() == 8 &&
              chain.chain[2 * 2 + 0 * 2 + 1] == 12 && chain.chain[2 * 2 + 1 * 2 + 0] == 2 && chain.logProbability[3] == -4.5 &&
              chain.logLikelihood[1] == -2 && chain.logPrior[2] == -0.5,
          "chain: [step][walker][parameter] as the engine reads samples.mcmc (" + why + ")");
    {
      std::ofstream f(b.parent_path() / "output" / "samples.mcmc");
      f << "step,walker,logprob,loglikelihood,logprior,param0\n0,0,-1,-1,0,1\n0,0,-1,-1,0,1\n1,1,-1,-1,0,1\n"
        << "1,0,-1,-1,0,1\n";
    }
    const bool twice = !front::ReadSamplerChain(p, chain, why) && why.find("twice") != std::string::npos;
    Check(twice, "chain: a row written twice refused (" + why + ")");
  }

  // (h) the Levels and Pairs tabs' facade (qt_gui.md section 11): columns by name, a channel at its level energy, a
  //     channel's radius from <channelRadii>, what names a pair
  {
    front::Project p;
    std::string why;
    Check(front::Project::Open(projects[0].string(), p, why), "T1 opened for the level functions");
    Check(front::LevelColumnNames().size() == 34 && front::LevelColumnIndex("levelE") == 2 &&
              front::LevelColumnIndex("gamma") == 11 && front::LevelColumnIndex("e1") == 33 &&
              front::LevelColumnIndex("nonsense") == -1,
          "the <levels> columns by name (34, levelE 2, gamma 11, e1 33)");
    const std::vector<front::LevelRow> rows = p.Levels();
    auto row = [&rows](const std::string &e, const std::string &key, const std::string &l2) {
      for (const front::LevelRow &r : rows)
        if (r.tokens.size() > 11 && r.tokens[2] == e && r.tokens[5] == key && r.tokens[7] == l2) return r;
      return front::LevelRow();
    };
    // 1/2+ at 2.365 MeV, p + 12C, l = 0, Gamma = 33000 eV, a = 3.6 fm: E_c = 2.365 - 1.9434900024034978 MeV.  The
    // values by mpmath 1.4.1 (coulombf, coulombg, diff; 30 digits; masses 1.00782503223 u and 12 u, 931.49410242
    // MeV/u, hbar c 197.3269804 MeV fm, alpha 1/137.035999084): P = 0.0112659697575860, S = -0.802139976806155,
    // gamma_W^2 = 5.20376154799012 MeV, Gamma_sp = 0.117250840450691 MeV, gamma^2 = Gamma/(2P) = 1.46458763471201
    // MeV, theta^2 = 0.281447875965357.
    const front::ChannelQuantities open = front::ChannelAtLevel(row("2.365", "1", "0"));
    auto near = [](double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); };
    Check(open.error.empty() && open.kind == "open" && open.unit == "eV" && near(open.channelEnergy, 0.4215099975965022, 1e-13),
          "T1 2.365 MeV p+12C l=0: an open channel, Gamma in eV, E_c 0.42151 MeV (" + open.kind + ")");
    Check(near(open.P, 0.0112659697575860, 1e-8) && near(open.S, -0.802139976806155, 1e-8),
          "T1 2.365 MeV: P and S of the engine = mpmath's to 1e-8 (P " + std::to_string(open.P) + ")");
    Check(near(open.wignerLimit, 5.20376154799012, 1e-12) && near(open.singleParticleWidth, 0.117250840450691, 1e-8) &&
              near(open.reducedWidth, 1.46458763471201, 1e-8) && near(open.wignerRatio, 0.281447875965357, 1e-8),
          "T1 2.365 MeV: gamma_W^2, Gamma_sp, gamma^2 and theta^2 = mpmath's");
    Check(near(open.singleParticleWidth, front::SingleParticleWidth(1, 6, 1.00782503223 * 12. / 13.00782503223, 0, 3.6,
                                                                    open.channelEnergy), 1e-14),
          "T1 2.365 MeV: Gamma_sp = SingleParticleWidth");
    const front::ChannelQuantities closed = front::ChannelAtLevel(row("0.0", "1", "2"));
    Check(closed.kind == "closed" && closed.unit == "fm^-1/2" && closed.channelEnergy < 0 && closed.P > 0 &&
              std::isfinite(closed.reducedWidth) && std::isnan(closed.singleParticleWidth),
          "T1 1/2- at 0 MeV p+12C l=1: closed, an ANC in fm^-1/2, gamma^2 = C^2/(2P) (" + closed.kind + ")");
    // below threshold, l = 1: P = mu a / (hbar^2 W^2), W = W_{-eta_b, 3/2}(2 kappa a); mpmath's whitw (30 digits,
    // the constants above, B = 1.9434900024034978 MeV): P = 1.03417877501829 MeV^-1 fm^-1, gamma^2 = 1.85^2 / (2P) =
    // 1.65469456668141 MeV (l = 0 would give P = 3.40154840962196)
    Check(near(closed.P, 1.03417877501829, 1e-7) && near(closed.reducedWidth, 1.65469456668141, 1e-7),
          "T1 1/2- at 0 MeV: the bound factor of l = 1 and gamma^2 = mpmath's (P " + std::to_string(closed.P) + ")");
    const front::ChannelQuantities moment = front::ChannelAtLevel(row("1.0", "3", "2"));
    const front::ChannelQuantities gamma = front::ChannelAtLevel(row("2.365", "2", "2"));
    Check(moment.kind == "magnetic_moment" && moment.unit == "μ_N" && gamma.kind == "photon" && gamma.unit == "eV",
          "T1: the 3/2- level's M1 to its own state is a magnetic moment; the 1/2+ level's E1 a gamma width (" +
              moment.kind + ", " + gamma.kind + ")");
    const front::ChannelQuantities own = front::ChannelAtLevel(row("2.365", "1", "0"), 4.5);
    Check(own.radius == 4.5 && own.P > open.P, "T1 2.365 MeV: a radius of 4.5 fm given: P grows");
    front::LevelRow broken = row("2.365", "1", "0");
    broken.tokens[11] = "x";
    Check(!front::ChannelAtLevel(broken).error.empty(), "a row the engine cannot read: an error");
    // <channelRadii>: none in T1 (the pair's radius); a record for pair 1, l = 0
    front::ChannelRadiusInfo r0 = front::ChannelRadius(p, 1, 0, 1, 1, 3.6);
    Check(!r0.own && r0.radius == 3.6, "T1: no <channelRadii>: the pair's radius");
    front::Record radius;
    radius.kind = "channel_radius";
    radius.name = "p12C_s";
    radius.fields = {{"pair", "p+12C", ""}, {"l", "0", ""}, {"radius", "4.2 fm", ""}};
    const bool set = p.SetRecord("channelRadii", radius, why);
    r0 = front::ChannelRadius(p, 1, 0, 1, 1, 3.6);
    const front::ChannelRadiusInfo r1 = front::ChannelRadius(p, 1, 1, 1, 1, 3.6);
    Check(set && r0.own && r0.radius == 4.2 && r0.record == "p12C_s" && !r1.own && r1.radius == 3.6,
          "T1: a channel_radius record gives pair 1, l = 0 its radius (4.2 fm), l = 1 keeps the pair's (" + why + ")");
    // what names a pair: T1's data sets name pairs 1 and 2
    const std::vector<std::string> uses1 = front::PairUses(p, 1), uses9 = front::PairUses(p, 9);
    std::string list;
    for (const std::string &u : uses1) list += " " + u;
    Check(!uses1.empty() && uses9.empty() && uses1[0].rfind("data: ", 0) == 0,
          "T1: the data sets that name pair 1:" + list);
  }

  // (f) physics (docs/upgrade/physics_tab.md section 7)
  {
    // The rows: the physics table pyazr.physics reads (26 on 2026-11-19), every row with its axis among its arguments.
    const std::vector<front::PhysicsFunction> rows = front::PhysicsFunctions();
    bool axes = !rows.empty(), penetrability = false;
    for (const front::PhysicsFunction &f : rows) {
      bool found = false;
      for (const front::PhysicsArgument &a : f.arguments) found = found || (a.name == f.axis && a.hasSweep);
      axes = axes && found;
      penetrability = penetrability || f.name == "penetrability";
    }
    Check(rows.size() >= 26 && axes && penetrability,
          std::to_string(rows.size()) + " physics rows, each with its axis an argument with a sweep");
    std::string why;
    front::PhysicsValues v;
    // A broadcast column, a value outside its domain (e_cm = -1: argument 7) and its words.
    const bool evaluated = front::EvaluatePhysics(
        "penetrability", {{0}, {1}, {6}, {1.00782503}, {12.0}, {5.0}, {0.5, 1.0, -1.0}, {0}}, v, why);
    Check(evaluated && v.status.size() == 3 && v.status[0] == 0 && v.status[2] == 7 && std::isnan(v.values[0][2]) &&
              v.values[0][1] > v.values[0][0],
          "penetrability on three energies: computed, rising, and e_cm = -1 refused as argument 7 (" + why + ")");
    Check(front::PhysicsStatusText("penetrability", 7).find("e_cm") != std::string::npos,
          "the status words name e_cm: " + front::PhysicsStatusText("penetrability", 7));
    why.clear();
    const bool unknownRow = !front::EvaluatePhysics("no_such_row", {}, v, why);
    Check(unknownRow && !why.empty(), "an unknown row is refused: " + why);

    // The stopping table (SRIM-2013): protons in Ta2O5.  total = electronic + nuclear, and two rows invert each
    // other: the areal density between E_in and energy_after(E_in, N) is N again.
    const std::vector<front::PhysicsFunction> stopping = front::StoppingFunctions();
    Check(stopping.size() >= 5 && stopping[0].arguments[0].kind == "text" && !front::StoppingCitation().empty(),
          std::to_string(stopping.size()) + " stopping rows, the formula first, SRIM's citation");
    why.clear();
    const bool sp = front::EvaluateStopping("stopping_power", "Ta2O5", {{0}, {1}, {1}, {1.00782503}, {0.3, 1.0}}, v,
                                            why);
    Check(sp && v.status == std::vector<int>({0, 0}) &&
              std::fabs(v.values[2][1] - v.values[0][1] - v.values[1][1]) <= 1e-12 * v.values[2][1] &&
              v.values[2][0] > v.values[2][1],
          "stopping_power of p in Ta2O5: total = electronic + nuclear, falling from 0.3 to 1 MeV (" + why + ")");
    front::PhysicsValues after, between;
    front::EvaluateStopping("energy_after", "Ta2O5", {{0}, {1}, {1}, {1.00782503}, {1.0}, {1e18}}, after, why);
    front::EvaluateStopping("areal_density_between", "Ta2O5", {{0}, {1}, {1}, {1.00782503}, {1.0}, after.values[0]},
                            between, why);
    Check(after.status[0] == 0 && between.status[0] == 0 && std::fabs(between.values[0][0] / 1e18 - 1.0) < 1e-6,
          "energy_after and areal_density_between invert each other to " +
              Sci(between.values.empty() ? 1.0 : std::fabs(between.values[0][0] / 1e18 - 1.0)));
    why.clear();
    const bool badFormula = !front::EvaluateStopping("range", "Xq7", {{0}, {1}, {1}, {1.0}, {1.0}}, v, why);
    Check(badFormula && !why.empty(), "a formula stopping/ refuses is refused: " + why);

    // A session of T1 (Brune's basis): pairs, groups, poles.
    const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", fs::path("front_test_work") / "physics");
    front::PhysicsSession s;
    why.clear();
    const bool opened = front::PhysicsSession::Open(t1.string(), "", false, s, why);
    Check(opened && s.IsOpen(), "T1: a physics session (" + why + ")");
    Check(s.Parametrization() == "brune", "T1: parametrization " + s.Parametrization());
    const std::vector<front::PhysicsPair> pairs = s.Pairs();
    Check(pairs.size() >= 2 && pairs[0].kind == "particle" && pairs[0].z1 * pairs[0].z2 == 6,
          "T1: pair 1 is p + 12C");
    int group = 0;
    for (const front::PhysicsGroup &g : s.Groups())
      if (g.J == 0.5 && g.parity == 1) group = g.jgroup;
    Check(group > 0, "T1: the J group 1/2+ is " + std::to_string(group));
    // The poles: the 1/2+ level at 2.365 MeV with Gamma_p = 33000 eV (T1.azr's <levels>), back from the engine's
    // transformation (Brune's basis: the numbers given are the observed ones).
    bool pole = false;
    for (const front::PhysicsPole &p : s.Poles())
      if (p.jgroup == group && std::fabs(p.energy - 2.365) < 1e-9 && !p.widths.empty())
        pole = std::fabs(p.widths[0].value - 33000.0) < 1e-6 && p.widths[0].kind == "width" &&
               std::fabs(p.widths[0].width - 0.033) < 1e-12 && p.totalWidth >= 0.033;
    Check(pole, "T1: the pole 1/2+ at E_x = 2.365 MeV with Gamma_p = 33000 eV");

    // The session's channel functions against the session-free row of the same pair (the row computes F and G
    // exactly; a session with a Coulomb table interpolates them to about 1e-7).
    front::PhysicsCurves cf;
    const std::vector<double> energies = {0.3, 0.6, 1.2};
    why.clear();
    const bool functions = s.ChannelFunctions(1, 0, 0.0, -1, energies, cf, why);
    Check(functions && cf.Find("P") && cf.Find("open"), "T1: channel functions of pair 1, l = 0 (" + why + ")");
    double radius = 0;
    for (const front::PhysicsGroup &g : s.Groups())
      for (const front::PhysicsChannel &c : g.channels)
        if (c.pair == 1 && c.l == 0) radius = c.radius;
    front::EvaluatePhysics("penetrability",
                           {{0}, {(double)pairs[0].z1}, {(double)pairs[0].z2}, {pairs[0].m1}, {pairs[0].m2}, {radius},
                            energies, {0}},
                           v, why);
    double worst = 0;
    for (std::size_t i = 0; i < energies.size() && cf.Find("P"); i++)
      worst = std::max(worst, std::fabs(cf.Find("P")->values[i] / v.values[0][i] - 1.0));
    Check(cf.Find("P") && worst < 1e-6, "T1: session P = the row's P to " + Sci(worst));

    // The 1/2+ block on the E_x of T1's extrapolation x3 (12C(p,p) phase shift, J = 1/2, l = 0): the facade's delta
    // against the program's output of `extrapolate` (AZUREOut, read by ReadCurves), and |U_pp| = 1 (one open particle
    // channel: the collision matrix is unitary there).
    Collect observer;
    front::RunRequest request;
    request.command = "extrapolate";
    request.project = t1.string();
    request.threads = 2;
    const front::RunResult run = front::Run(request, observer, nullptr);
    front::Project project;
    front::Project::Open(t1.string(), project, why);
    const front::Curves curves = front::ReadCurves(project);
    const front::Curve *x3 = nullptr;
    for (const front::Curve &c : curves.curves)
      if (c.record == "x3") x3 = &c;
    Check(run.exitCode == 0 && x3 && x3->points.size() > 10, "T1: extrapolate wrote the phase shifts of x3");
    if (x3) {
      // E_x of each point from the record's lab energies (the file prints 7 digits; dDelta/dE reaches 5000 deg/MeV at
      // the level): the row lab_to_cm_energy, the engine's conversion of a lab point, plus the pair's threshold.
      std::vector<double> lab;
      for (const front::Record &r : project.Records("extrapolations"))
        if (r.name == "x3") {
          std::istringstream words(r.Value("energies"));
          std::string w;
          while (words >> w)
            if (w != "list" && w.find_first_not_of("0123456789.eE+-") == std::string::npos) lab.push_back(std::stod(w));
        }
      front::EvaluatePhysics("lab_to_cm_energy", {{pairs[0].m1}, {pairs[0].m2}, lab, {0}}, v, why);
      std::vector<double> ex;
      for (std::size_t i = 0; i < lab.size() && v.values.size() == 1; i++)
        ex.push_back(v.values[0][i] + pairs[0].separation + pairs[0].excitation);
      Check(ex.size() == x3->points.size(), "T1: " + std::to_string(ex.size()) + " energies of x3");
      front::PhysicsCurves block;
      why.clear();
      const bool blocked = s.Block(group, ex, block, why);
      Check(blocked, "T1: the 1/2+ block on " + std::to_string(ex.size()) + " energies (" + why + ")");
      const front::PhysicsSeries *delta = block.Find("delta", {{"channel", 1}, {"channel_out", 1}});
      const front::PhysicsSeries *eta = block.Find("eta", {{"channel", 1}, {"channel_out", 1}});
      double dmax = 0, umax = 0;
      for (std::size_t i = 0; delta && eta && i < ex.size() && i < x3->points.size(); i++) {
        // the file's 7 significant digits
        dmax = std::max(dmax, std::fabs(delta->values[i] - x3->points[i].model) / std::max(1.0, std::fabs(x3->points[i].model)));
        umax = std::max(umax, std::fabs(eta->values[i] - 1.0));
      }
      Check(delta && dmax < 1e-6, "T1: delta of the block = the program's phase shift of x3 to " +
                                      Sci(dmax) + " relative");
      Check(eta && umax < 1e-12, "T1: |U_pp| = 1 to " + Sci(umax));
      // log10 |det A^-1| on a grid across the 2.365 MeV level (Gamma = 33 keV): its minimum lies within Gamma.
      std::vector<double> grid;
      for (int i = 0; i <= 400; i++) grid.push_back(2.2 + 0.3 * i / 400.0);
      s.Block(group, grid, block, why);
      const front::PhysicsSeries *det = block.Find("log10_abs_det");
      std::size_t at = 0;
      for (std::size_t i = 0; det && i < grid.size(); i++)
        if (det->values[i] < det->values[at]) at = i;
      Check(det && std::fabs(grid[at] - 2.365) < 0.033,
            "T1: |det A^-1| is least at E_x = " + std::to_string(grid[at]) + " MeV, within Gamma of 2.365");
    }
    why.clear();
    const bool noGroup = !s.Block(99, {2.0}, cf, why);
    Check(noGroup && !why.empty(), "T1: a J group that is not there is refused: " + why);

    // THM: T2's experiment LC09, its vertex at three energies (pyazr's Session.thm_vertex reads the same table).
    const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", fs::path("front_test_work") / "physics_t2");
    front::PhysicsSession s2;
    why.clear();
    const bool opened2 = front::PhysicsSession::Open(t2.string(), "", false, s2, why);
    Check(opened2, "T2: a physics session (" + why + ")");
    const std::vector<std::string> experiments = s2.ThmExperiments();
    Check(experiments.size() == 1 && experiments[0] == "LC09", "T2: the THM experiment LC09");
    front::PhysicsCurves vertex;
    why.clear();
    bool positive = s2.ThmVertex("LC09", {0.3, 0.5, 0.8}, vertex, why) && vertex.Find("M2");
    for (const front::PhysicsSeries &series : vertex.series)
      if (series.name == "M2")
        for (double x : series.values) positive = positive && x > 0 && std::isfinite(x);
    Check(positive, "T2: <|M_l|^2> of LC09 finite and positive (" + why + ")");
    why.clear();
    const bool noExperiment = !s2.ThmVertex("none", {0.5}, vertex, why);
    Check(noExperiment && !why.empty(), "T2: an unknown experiment refused: " + why);
  }

  FormsChecks(source);
  PolishChecks(source);

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
