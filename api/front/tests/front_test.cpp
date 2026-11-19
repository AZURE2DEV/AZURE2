// The front ends' facade without Qt (docs/upgrade/qt_gui.md section 10): front_test <source dir> <AZURE2 binary>.
//  (a) the tables: schema, commands, options, fitters, <config> keys;
//  (b) the test projects T1, T2, T3 and 13N opened and written back byte for byte, also after every record, every
//      <levels> row and every <config> entry is set to what it holds; an edit changes only its lines and is undone;
//      a record added, renamed and removed; --validate's diagnostics, a broken field named by its code;
//  (c) an old project refused on opening and converted by the converter of --convert;
//  (d) a run of `calculate` through the facade writes the same output files, byte for byte, as the program on the same
//      project, and returns the chi^2 the program prints; T1 and 13N;
//  (e) a fit stopped through the facade ends at its best point with exit 0 and fitStatus.out saying "stopped".
// Each check prints PASS or FAIL; exit 1 on a failure.
#include "azure2/front/Engine.h"

#include <map>
#include <atomic>
#include <cmath>
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

}  // namespace

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

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
