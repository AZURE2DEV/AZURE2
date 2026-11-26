// azure2::front's Overview calls (Engine.h; docs/upgrade/qt_gui.md section 15): the project as it stands
// (SummarizeProject), the last run as its files tell it (ReadLastRun) and the run's parameters against the project's
// set-up (ReadRunParameters).  Nothing is recomputed: the numbers are the files' (chiSquared.out, written by
// WorkingCopy::WriteChiSquaredFile; fitStatus.out, by WriteFitStatus in src/ProgramFit.cpp; param.sav / param.par, by
// AZUREParams::WriteUserParameters), matched to the project's records by the engine's own build of the project, as
// ReadCurves (FrontOutputs.cpp) matches the curves.
#include "azure2/front/Engine.h"
#include "azure2/support/ClassicNumeric.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <system_error>

#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/ConfigSection.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;
namespace fs = std::filesystem;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

/// The build of the project's text, as ReadCurves makes it (the reader and the validator; no set-up checks).
ex::Format3Project Built(const Project &project) {
  ex::Format3ReadOptions read;
  read.fileName = project.Path();
  read.projectDir = project.Directory();
  read.setUpChecks = false;
  return ex::ReadFormat3Text(project.Text(), read);
}

std::string SpinLabel(int twoJ, int parity) {
  const std::string j = twoJ % 2 ? std::to_string(twoJ) + "/2" : std::to_string(twoJ / 2);
  return j + (parity < 0 ? "-" : "+");
}

/// The records of a section's kind `kind`, in file order (a data set's key is its place among them, from 1).
std::vector<Record> RecordsOfKind(const Project &project, const std::string &section, const std::string &kind) {
  std::vector<Record> out;
  for (const Record &r : project.Records(section))
    if (r.kind == kind) out.push_back(r);
  return out;
}

/// The include switch of a kind (KindInfo::includeField): a record is active unless it says the off word.
bool RecordActive(const Record &r) {
  for (const KindInfo &k : RecordKinds())
    if (k.kind == r.kind && !k.includeField.empty()) return r.Value(k.includeField) != k.includeOff;
  return true;
}

/// A file's modification time (precise, for comparisons) and in seconds since 1970 (for a caller).
struct FileTime {
  bool exists = false;
  fs::file_time_type time;
  long long seconds = 0;
};
FileTime TimeOf(const fs::path &path) {
  FileTime t;
  std::error_code ec;
  if (!fs::is_regular_file(path, ec)) return t;
  t.time = fs::last_write_time(path, ec);
  if (ec) return t;
  t.exists = true;
  // the file clock's epoch is not the system clock's before C++20: the offset of the two clocks now
  const auto system = std::chrono::system_clock::now() +
                      std::chrono::duration_cast<std::chrono::system_clock::duration>(t.time -
                                                                                      fs::file_time_type::clock::now());
  t.seconds = (long long)std::chrono::duration_cast<std::chrono::seconds>(system.time_since_epoch()).count();
  return t;
}

double ToDouble(const std::string &s) {
  char *end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  return end && end != s.c_str() && *end == '\0' ? v : kNaN;
}

/// One data line of chiSquared.out: "key,chi2,N,norm,normChi2".
bool ReadChiLine(const std::string &line, DataSetResult &r) {
  double norm = 1, normChi = 0;
  const int n = std::sscanf(line.c_str(), "%d,%lf,%d,%lf,%lf", &r.key, &r.chiSquared, &r.points, &norm, &normChi);
  if (n < 3) return false;
  if (n >= 4) r.norm = norm;
  if (n >= 5) r.normChiSquared = normChi;
  r.chiSquaredPerPoint = r.points > 0 ? r.chiSquared / r.points : kNaN;
  return true;
}

/// The total line: "Label: value Label: value ...", every pair in order.
void ReadTotalLine(const std::string &line, std::vector<std::pair<std::string, double>> &totals) {
  std::istringstream s(line);
  std::string token, label;
  while (s >> token) {
    if (token.size() > 1 && token.back() == ':') {
      label = token.substr(0, token.size() - 1);
    } else if (!label.empty()) {
      totals.push_back({label, ToDouble(token)});
      label.clear();
    }
  }
}

}  // namespace

ProjectSummary SummarizeProject(const Project &project) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  ProjectSummary out;
  out.file = project.Path();
  out.pairs = Pairs(project);
  for (const PairChoice &p : out.pairs) {
    bool seen = false;
    for (const std::string &n : out.nuclei) seen = seen || n == p.nucleus;
    if (!seen && !p.nucleus.empty()) out.nuclei.push_back(p.nucleus);
  }
  const ex::Format3Project built = Built(project);
  out.built = built.Ok();
  for (const ex::LevelInfo &level : built.levels.levels) {
    const std::string nucleus = ex::NuclideLabel(level.nucleusZ, level.nucleusA);
    LevelGroupCount *group = nullptr;
    for (LevelGroupCount &g : out.levelGroups)
      if (g.nucleus == nucleus && g.twoJ == level.twoJ && g.parity == level.parity) group = &g;
    if (!group) {
      out.levelGroups.push_back({nucleus, level.twoJ, level.parity, SpinLabel(level.twoJ, level.parity), 0});
      group = &out.levelGroups.back();
    }
    group->levels++;
    out.levels++;
  }
  out.externalCapturePairs = built.levels.externalCapturePairs;

  // the data sets and extrapolations: active as the build says (the record's switch when the project does not build)
  std::vector<Record> activeData;
  for (int pass = 0; pass < 2; pass++) {
    const bool extrapolation = pass == 1;
    const std::vector<Record> records = RecordsOfKind(project, extrapolation ? "extrapolations" : "data",
                                                      extrapolation ? "extrapolation" : "dataset");
    std::vector<bool> active(records.size(), false);
    if (out.built) {
      auto mark = [&active](int key, bool on) {
        if (key >= 1 && key <= (int)active.size()) active[key - 1] = on;
      };
      if (extrapolation)
        for (const auto &i : built.built.extrapolations.extrapolations.items) mark(i.spec.key, i.spec.active);
      else
        for (const auto &i : built.built.data.datasets.items) mark(i.spec.key, i.spec.active);
    } else {
      for (std::size_t i = 0; i < records.size(); i++) active[i] = RecordActive(records[i]);
    }
    int on = 0;
    for (std::size_t i = 0; i < records.size(); i++)
      if (active[i]) {
        on++;
        if (!extrapolation) activeData.push_back(records[i]);
      }
    (extrapolation ? out.extrapolations : out.dataSets) = (int)records.size();
    (extrapolation ? out.activeExtrapolations : out.activeDataSets) = on;
  }
  // the effect records by kind; one is used when an active data set has a field whose value is its name
  std::set<std::string> named;
  for (const Record &r : activeData)
    for (const Field &f : r.fields) {
      std::istringstream s(f.value);
      std::string first;
      s >> first;
      named.insert(first);
    }
  for (const Record &r : project.Records("targets")) {
    RecordUse *use = nullptr;
    for (RecordUse &u : out.effects)
      if (u.kind == r.kind) use = &u;
    if (!use) {
      out.effects.push_back({"targets", r.kind, 0, 0});
      use = &out.effects.back();
    }
    use->records++;
    if (named.count(r.name)) use->used++;
  }

  // the options of the calculation: the key the project writes, else the row's default
  for (const io::ConfigKey &k : io::ConfigKeys()) {
    if (!k.option) continue;
    OptionState o;
    o.name = k.name;
    o.title = k.title;
    o.effect = k.effect;
    o.defaultOn = k.defaultOn;
    o.on = k.defaultOn;
    const std::string value = project.ConfigValue(k.name, "\x01");
    bool on = false;
    if (value != "\x01") {
      o.written = true;
      if (io::ConfigOptionValue(value, on)) o.on = on;
      else out.notes.push_back(std::string("<config> ") + k.name + "=" + value + " does not read as 0 or 1");
    }
    out.options.push_back(o);
  }
  auto option = [&out](const std::string &name) {
    for (const OptionState &o : out.options)
      if (o.name == name) return o.on;
    return false;
  };
  out.basis = !option("use_brune") ? "formal (Lane-Thomas)" : option("use_park") ? "Park" : "Brune";
  const std::string amatrix = project.ConfigValue("amatrix", "true");
  out.matrix = amatrix == "false" || amatrix == "0" ? "R-matrix" : "A-matrix";
  if (!out.built) out.notes.push_back("the project has errors (Validate lists them): the records' activity is the switch as written");
  return out;
}

RunSummary ReadLastRun(const Project &project) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  RunSummary out;
  out.outputDirectory = ex::PathToOpen(project.Directory(), project.ConfigValue("output_dir"));
  if (!out.outputDirectory.empty() && out.outputDirectory.back() != '/') out.outputDirectory += '/';
  const fs::path dir(out.outputDirectory);
  const FileTime projectTime = TimeOf(project.Path());
  out.projectTime = projectTime.seconds;

  std::error_code ec;
  FileTime newestExtrap;
  if (fs::is_directory(dir, ec)) {
    for (const auto &e : fs::directory_iterator(dir, ec)) {
      if (!e.is_regular_file()) continue;
      const std::string name = e.path().filename().string();
      const FileTime t = TimeOf(e.path());
      out.files.push_back({name, t.seconds});
      if (name.size() > 7 && name.compare(name.size() - 7, 7, ".extrap") == 0 && name.rfind("AZUREOut", 0) == 0 &&
          (!newestExtrap.exists || t.time > newestExtrap.time))
        newestExtrap = t;
    }
    std::sort(out.files.begin(), out.files.end(),
              [](const OutputFileTime &a, const OutputFileTime &b) { return a.name < b.name; });
  }

  const FileTime chiTime = TimeOf(dir / "chiSquared.out");
  const FileTime fitTime = TimeOf(dir / "fitStatus.out");
  const FileTime parTime = TimeOf(dir / "param.par");
  const FileTime savTime = TimeOf(dir / "param.sav");
  FileTime runTime;
  if (chiTime.exists) {
    runTime = chiTime;
    // the order a run writes (docs/upgrade/qt_gui.md section 15): param.par at the set-up, fitStatus.out after the
    // fitter, chiSquared.out at the final evaluation
    const bool fit = fitTime.exists && fitTime.time <= chiTime.time && (!parTime.exists || parTime.time <= fitTime.time);
    out.command = fit ? "fit" : "calculate";
    if (newestExtrap.exists && newestExtrap.time > chiTime.time)
      out.notes.push_back("an extrapolation ran after this run (its files are newer than chiSquared.out)");
  } else if (newestExtrap.exists) {
    runTime = newestExtrap;
    out.command = "extrapolate";
  }
  out.present = runTime.exists;
  if (!out.present) {
    out.notes.push_back("no run yet: the output directory holds no chiSquared.out and no extrapolation");
    return out;
  }
  out.runTime = runTime.seconds;
  out.projectNewer = projectTime.exists && projectTime.time > runTime.time;
  out.stale = out.projectNewer;
  if (out.command == "fit" && savTime.exists) out.parameterFile = (dir / "param.sav").string();
  else if (parTime.exists) out.parameterFile = (dir / "param.par").string();

  if (out.command == "fit") {
    std::ifstream in(dir / "fitStatus.out");
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream s(line);
      std::string key, value;
      if (!(s >> key)) continue;
      std::getline(s >> std::ws, value);
      while (!value.empty() && (value.back() == '\r' || value.back() == ' ')) value.pop_back();
      out.fitStatus.push_back({key, value});
      if (key == "status") out.status = value;
      if (key == "minimizer") out.minimizer = value;
      if (key == "calls") out.calls = std::atol(value.c_str());
    }
  }

  if (chiTime.exists) {
    std::ifstream in(dir / "chiSquared.out");
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
      if (header) {  // "Segment#, Chi-Squared, N, Norm, Norm-Chi-Squared"
        header = false;
        continue;
      }
      DataSetResult r;
      if (line.rfind("Total-", 0) == 0) ReadTotalLine(line, out.totals);
      else if (ReadChiLine(line, r)) out.dataSets.push_back(r);
    }
    out.haveChiSquared = true;
    for (const auto &t : out.totals) {
      if (t.first == "Total-Chi-Squared") out.totalChiSquared = t.second;
      if (t.first == "Total-N") out.points = (int)std::lround(t.second);
    }
  }

  // the records of the keys, the norms' priors, the active data sets the run left out (the build of the project)
  const ex::Format3Project built = Built(project);
  const std::vector<Record> records = RecordsOfKind(project, "data", "dataset");
  if (!built.Ok()) {
    out.notes.push_back("the project has errors (Validate lists them): the run's data sets are shown by key");
    for (DataSetResult &r : out.dataSets) r.normPull = kNaN;
    return out;
  }
  std::map<int, const ex::DatasetSpec *> specs;
  std::map<std::string, int> groupOwner;  // a norm group's prior belongs to its first active member
  for (const auto &i : built.built.data.datasets.items) {
    specs[i.spec.key] = &i.spec;
    if (i.spec.active && !i.spec.normGroup.empty() && !groupOwner.count(i.spec.normGroup))
      groupOwner[i.spec.normGroup] = i.spec.key;
  }
  std::set<int> used;
  for (DataSetResult &r : out.dataSets) {
    used.insert(r.key);
    r.normPull = kNaN;
    if (r.key >= 1 && r.key <= (int)records.size()) {
      const Record &rec = records[r.key - 1];
      r.record = rec.name;
      r.reaction = rec.Value("reaction");
      r.quantity = rec.Value("quantity");
    }
    const auto s = specs.find(r.key);
    if (s == specs.end()) continue;
    const ex::DatasetSpec &spec = *s->second;
    const bool owner = spec.normGroup.empty() || groupOwner[spec.normGroup] == spec.key;
    r.normPrior = spec.varyNorm && spec.normError > 0 && owner;
    if (r.normPrior) {
      const double centre = spec.normCentre.Or(spec.norm);
      r.normPull = (r.norm < centre ? -1.0 : 1.0) * std::sqrt(std::max(0.0, r.normChiSquared));
    }
  }
  if (out.haveChiSquared)
    for (const auto &i : built.built.data.datasets.items)
      if (i.spec.active && !used.count(i.spec.key) && i.spec.key >= 1 && i.spec.key <= (int)records.size())
        out.leftOut.push_back(records[i.spec.key - 1].name);
  // an input file written after the run: its record's results are stale too
  for (const InputFile &f : InputFiles(project)) {
    const FileTime t = TimeOf(f.path);
    bool listed = false;
    for (const std::string &n : out.changedInputs) listed = listed || n == f.record;
    if (t.exists && t.time > runTime.time && !listed) out.changedInputs.push_back(f.record);
  }
  if (!out.changedInputs.empty()) out.stale = true;
  return out;
}

std::vector<InputFile> InputFiles(const Project &project) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<InputFile> out;
  const ex::Format3Project built = Built(project);
  if (!built.Ok()) return out;
  auto add = [&out, &project](const std::string &file, const std::string &section, const std::string &record,
                              const char *role) {
    if (!file.empty()) out.push_back({ex::PathToOpen(project.Directory(), file), section, record, role});
  };
  const std::vector<Record> data = RecordsOfKind(project, "data", "dataset");
  for (const auto &i : built.built.data.datasets.items) {
    if (!i.spec.active || i.spec.key < 1 || i.spec.key > (int)data.size()) continue;
    const std::string &name = data[i.spec.key - 1].name;
    add(i.spec.dataFile, "data", name, "data");
    add(i.spec.covarianceFile, "data", name, "covariance");
  }
  const std::vector<Record> extrapolations = RecordsOfKind(project, "extrapolations", "extrapolation");
  for (const auto &i : built.built.extrapolations.extrapolations.items)
    if (i.spec.active && i.spec.key >= 1 && i.spec.key <= (int)extrapolations.size())
      add(i.spec.gridFile, "extrapolations", extrapolations[i.spec.key - 1].name, "grid");
  return out;
}

bool ReadRunParameters(const std::string &projectFile, const std::string &parameterFile, bool fit,
                       RunParameters &out, std::string &why) {
  out = RunParameters();
  // the set-up writes its external-capture integrals where its output directory is: a scratch directory, so that the
  // project's outputs are read, never written (an intEC.dat of the output directory is copied there to be reused)
  std::error_code ec;
  const fs::path scratch = fs::temp_directory_path(ec) /
                           ("azure2-overview-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
                            "-" + std::to_string(reinterpret_cast<std::uintptr_t>(&out)));
  fs::create_directories(scratch, ec);
  if (ec) {
    why = "no scratch directory for the set-up: " + ec.message();
    return false;
  }
  const fs::path outputs = !parameterFile.empty() ? fs::path(parameterFile).parent_path() : fs::path();
  for (const char *cache : {"intEC.dat", "intEC.dat.sig"})
    if (!outputs.empty() && fs::is_regular_file(outputs / cache, ec)) fs::copy_file(outputs / cache, scratch / cache, ec);
  std::vector<Parameter> parameters;
  const bool read = ReadParameters(projectFile, scratch.string(), parameters, why);
  fs::remove_all(scratch, ec);
  if (!read) return false;
  std::map<std::string, ParameterValue> values;
  if (!parameterFile.empty()) {
    std::vector<ParameterValue> rows;
    if (!ReadParameterFile(parameterFile, rows, why)) return false;
    for (const ParameterValue &v : rows) values[v.name] = v;
    out.file = parameterFile;
  }
  out.parameters = (int)parameters.size();
  for (const Parameter &p : parameters) {
    if (p.fixed) continue;
    out.free++;
    bool counted = false;
    for (auto &k : out.freeByKind)
      if (k.first == p.kind) k.second++, counted = true;
    if (!counted) out.freeByKind.push_back({p.kind, 1});
    const auto v = values.find(p.name);
    if (v == values.end()) continue;
    ParameterFlag f;
    f.name = p.name;
    f.kind = p.kind;
    f.unit = p.unit;
    f.value = v->second.value;
    f.error = v->second.error;
    f.lower = p.lower;
    f.upper = p.upper;
    auto near = [](double x, double limit) {
      return std::isfinite(limit) && std::fabs(x - limit) <= 1e-6 * std::max(1.0, std::fabs(limit));
    };
    // a parameter without limits has lower == upper (0 0 or -inf +inf): never at a limit
    f.atLimit = p.lower < p.upper && (near(f.value, p.lower) || near(f.value, p.upper));
    f.noError = fit && (!std::isfinite(f.error) || f.error == 0);
    if (f.atLimit || f.noError) out.flagged.push_back(f);
  }
  return true;
}

double ChiSquaredPerDegreeOfFreedom(const RunSummary &run, const RunParameters &parameters) {
  const int dof = run.points - parameters.free;
  return dof > 0 ? run.totalChiSquared / dof : kNaN;
}

}  // namespace front
}  // namespace azure2
