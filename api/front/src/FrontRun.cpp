// azure2::front::Run and the tables a front end builds its forms from (Engine.h; docs/upgrade/qt_gui.md section 3):
// the program's run (RunProgram) with a stream that hands each line to the caller's observer, the stop request and the
// sinks of ProgramObserver, and the program's own tables of commands, options, fitters and <config> keys.
#include "azure2/front/Engine.h"
#include "azure2/support/ClassicNumeric.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <mutex>
#include <ostream>
#include <sstream>
#include <streambuf>

#include <gsl/gsl_errno.h>

#include "CommandLineValues.h"
#include "Config.h"
#include "Fitters.h"
#include "GSLException.h"
#include "ProgramMain.h"
#include "ProgramObserver.h"
#include "ProjectCommands.h"
#include "RunRequest.h"
#include "azure2/convert/Convert.h"
#include "azure2/experiment/ProjectPaths.h"
#include "azure2/experiment/Schema.h"
#include "azure2/io/ConfigSection.h"
#include "azure2/support/Environment.h"
#include "azure2/support/FloatingPointPrecision.h"

namespace azure2 {
namespace front {

RunObserver::~RunObserver() = default;

namespace {

/// A stream buffer that hands every line to the observer: a line ends at '\n', or at '\r' (a console progress line
/// that the next overwrites: passed as transient).  Written from the run's threads under the engine's own locks or from
/// one thread; the buffer locks anyway.
class LineBuffer : public std::streambuf {
 public:
  explicit LineBuffer(RunObserver &observer) : observer_(observer) {}
  ~LineBuffer() override { Flush(false); }
  void Flush(bool transient) {
    if (!line_.empty()) observer_.Message(line_, transient);
    line_.clear();
  }

 protected:
  int_type overflow(int_type c) override {
    if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
    std::lock_guard<std::mutex> lock(mutex_);
    const char ch = traits_type::to_char_type(c);
    if (ch == '\n') {
      observer_.Message(line_, false);
      line_.clear();
    } else if (ch == '\r') {
      if (!line_.empty()) observer_.Message(line_, true);
      line_.clear();
    } else {
      line_.push_back(ch);
    }
    return c;
  }
  std::streamsize xsputn(const char *s, std::streamsize n) override {
    for (std::streamsize i = 0; i < n; i++) overflow(traits_type::to_int_type(s[i]));
    return n;
  }

 private:
  RunObserver &observer_;
  std::string line_;
  std::mutex mutex_;
};

}  // namespace

RunResult Run(const RunRequest &request, RunObserver &observer, const std::atomic<bool> *stop) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  // The calling thread's x87 precision (FloatingPointPrecision.h), as main's; the regions of the run take their own.
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  gsl_set_error_handler(&GSLException::GSLErrorHandler);  // as main: the engine catches GSL's errors as exceptions

  RunResult result;
  const auto started = std::filesystem::file_time_type::clock::now();
  std::vector<std::string> arguments;
  arguments.push_back(request.command);
  for (const std::string &o : request.options) arguments.push_back(o);
  arguments.push_back(request.project);

  // The threads of this run, in this thread (OpenMP's thread count is the calling thread's): --threads=N as main
  // applies it, the environment record read from the environment alone (as a pyazr session's).
  azure2::ThreadsInForce threads;
  {
    std::vector<std::string> threadOption{"--threads=" + std::to_string(request.threads < 1 ? 1 : request.threads)};
    try {
      threads = azure2::ApplyThreadsOption(threadOption, azure2::support::ProcessEnvironment());
    } catch (const std::invalid_argument &e) {
      observer.Message(std::string("ERROR: ") + e.what(), false);
      result.exitCode = 2;
      return result;
    }
  }

  ProgramObserver programObserver;
  programObserver.stop = stop;
  std::mutex chiMutex;
  programObserver.progress = [&observer](const FitProgress &record) {
    RunProgress p;
    p.stage = record.stage;
    p.evaluation = record.kind == FitProgress::Kind::kIteration ? record.iteration : record.evaluation;
    p.calls = record.calls;
    p.objective = record.objective;
    p.best = record.best;
    observer.Progress(p);
  };
  programObserver.dataSetChiSquared = [&result, &chiMutex](int segment, double chi2, int points) {
    std::lock_guard<std::mutex> lock(chiMutex);
    result.dataSets.push_back({segment, chi2, points});
  };
  programObserver.samplerStep = [&observer](int step, int total, double logPosterior, double logLikelihood, double) {
    RunProgress p;
    p.stage = "mcmc";
    p.evaluation = step;
    p.calls = total;
    p.objective = logPosterior;
    p.best = logLikelihood;
    observer.Progress(p);
  };
  programObserver.totalChiSquared = [&result, &chiMutex](double chi2) {
    std::lock_guard<std::mutex> lock(chiMutex);
    result.haveChiSquared = true;
    result.totalChiSquared = chi2;
  };

  {
    LineBuffer buffer(observer);
    std::ostream stream(&buffer);
    Config configure(stream);
    configure.observer = &programObserver;
    try {
      result.exitCode = RunProgram(configure, arguments, threads);
    } catch (const std::exception &e) {
      stream << "ERROR: " << e.what() << "\nCalculation was aborted." << std::endl;
      result.exitCode = 1;
    }
    stream.flush();
    buffer.Flush(false);
    result.outputDirectory = configure.outputdir.empty() ? std::string() : configure.OpenPath(configure.outputdir);
  }
  result.stopped = stop && stop->load();

  // The files of the output directory the run wrote (modified since it started).
  std::error_code ec;
  if (!result.outputDirectory.empty() && std::filesystem::is_directory(result.outputDirectory, ec)) {
    for (const auto &entry : std::filesystem::directory_iterator(result.outputDirectory, ec)) {
      if (!entry.is_regular_file(ec)) continue;
      if (entry.last_write_time(ec) >= started) result.outputFiles.push_back(entry.path().filename().string());
    }
    std::sort(result.outputFiles.begin(), result.outputFiles.end());
  }
  return result;
}

std::string SchemaJson() { return azure2::experiment::Schema::ToJson(); }

namespace {

std::string Quoted(const std::string &s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    if (c == '\n') {
      out += "\\n";
      continue;
    }
    out += c;
  }
  return out + "\"";
}

}  // namespace

std::string ConfigTableJson() {
  std::ostringstream out;
  out << "[";
  bool first = true;
  auto words = [](const char *s) {
    std::istringstream in(s ? s : "");
    std::string w, list;
    while (in >> w) list += (list.empty() ? "" : ",") + Quoted(w);
    return "[" + list + "]";
  };
  for (const io::ConfigPositional &p : io::ConfigPositionals()) {
    out << (first ? "\n " : ",\n ") << "{\"name\":" << Quoted(p.name) << ",\"kind\":" << Quoted(p.token ? "word" : "path")
        << ",\"title\":" << Quoted(p.title) << ",\"group\":" << Quoted(p.group) << ",\"doc\":" << Quoted(p.doc)
        << ",\"values\":" << words(p.values) << ",\"positional\":true,\"option\":false,\"rare\":"
        << (p.rare ? "true" : "false") << "}";
    first = false;
  }
  for (const io::ConfigKey &k : io::ConfigKeys()) {
    const char *kind = k.kind == io::ConfigKey::Bool     ? "bool"
                       : k.kind == io::ConfigKey::Int    ? "int"
                       : k.kind == io::ConfigKey::Number ? "number"
                       : k.kind == io::ConfigKey::Path   ? "path"
                                                         : "version";
    out << ",\n {\"name\":" << Quoted(k.name) << ",\"kind\":" << Quoted(kind) << ",\"minimum\":" << k.minimum
        << ",\"title\":" << Quoted(k.title) << ",\"group\":" << Quoted(k.group) << ",\"doc\":" << Quoted(k.doc)
        << ",\"absent\":" << Quoted(k.absent) << ",\"positional\":false,\"option\":" << (k.option ? "true" : "false")
        << ",\"flag\":" << Quoted(k.flag) << ",\"negation\":" << Quoted(k.negation)
        << ",\"default\":" << (k.defaultOn ? "true" : "false") << ",\"effect\":" << Quoted(k.effect)
        << ",\"rare\":" << (k.rare ? "true" : "false") << "}";
  }
  out << "\n]";
  return out.str();
}

std::vector<Command> Commands() {
  std::vector<Command> out;
  for (const RunCommandRow &row : RunCommands()) out.push_back({row.name, row.help, row.built});
  return out;
}

namespace {

/// The unit a row's help names last: the token after its last ", " when that is one word of letters ("GK", "keV").
std::string UnitOfHelp(const std::string &help) {
  const std::size_t comma = help.rfind(", ");
  if (comma == std::string::npos) return "";
  const std::string tail = help.substr(comma + 2);
  if (tail.empty() || tail.size() > 6) return "";
  for (const char c : tail)
    if (!std::isalpha(static_cast<unsigned char>(c))) return "";
  return tail;
}

}  // namespace

std::vector<Option> Options() {
  std::vector<Option> out;
  for (const RunQuestion &q : RunQuestions()) {
    Option o;
    o.option = q.option;
    o.kind = q.kind == AnswerKind::kFlag ? "flag" : q.kind == AnswerKind::kWord ? "word"
             : q.kind == AnswerKind::kNumber ? "number" : q.kind == AnswerKind::kWholeNumber ? "whole" : "file";
    if (o.kind != "flag") o.option += "=";
    for (const RunCommandRow &row : RunCommands())
      if (q.commands & CommandBit(row.command)) o.commands.push_back(row.name);
    o.required = q.required;
    o.low = q.low;
    o.high = q.high;
    o.defaultValue = q.defaultValue;
    o.help = q.help;
    for (std::size_t from = 0; q.words && *q.words;) {  // "rwa|physical"
      const std::string words = q.words;
      const std::size_t bar = words.find('|', from);
      o.words.push_back(words.substr(from, bar == std::string::npos ? std::string::npos : bar - from));
      if (bar == std::string::npos) break;
      from = bar + 1;
    }
    o.lowOpen = q.lowOpen;
    o.unit = UnitOfHelp(o.help);
    o.instead = q.instead ? q.instead : "";
    switch (q.when) {
      case AskedWhen::kAlways: break;
      case AskedWhen::kBandChosen: o.when = "with --covariance-band", o.needs = "--covariance-band"; break;
      case AskedWhen::kTemperaturesFromFile: o.when = "temperatures from a file"; break;
      case AskedWhen::kTemperaturesInRange: o.when = "without --temperature-file"; break;
      case AskedWhen::kProjectNeedsCapture: o.when = "when the run has external capture"; break;
    }
    out.push_back(o);
  }
  Option maxCalls;
  maxCalls.option = "--max-calls=";
  maxCalls.kind = "whole";
  maxCalls.commands = {"fit", "minos"};
  maxCalls.high = 1e9;
  maxCalls.help = "the fit's call limit (0: none)";
  out.push_back(maxCalls);
  return out;
}

std::vector<Fitter> Fitters() {
  std::vector<Fitter> out;
  for (const FitterRow &row : ::Fitters()) out.push_back({row.name, row.option, row.doc});
  return out;
}

bool ConvertProject(const std::string &from, const std::string &to, std::string &report) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  ConvertLineProblem ignored;
  return ConvertProject(from, to, report, ignored);
}

bool WriteRepairedCopy(const std::string &from, int line, const std::string &text, const std::string &to, std::string &why) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::ifstream in(from, std::ios::binary);
  if (!in) {
    why = "cannot read " + from;
    return false;
  }
  std::ostringstream bytes;
  bytes << in.rdbuf();
  const std::string all = bytes.str();
  std::error_code ec;
  if (std::filesystem::exists(to, ec)) {
    why = to + " exists: nothing is overwritten";
    return false;
  }
  std::string out;
  std::size_t start = 0;
  int number = 0;
  bool replaced = false;
  while (start < all.size()) {
    std::size_t end = all.find('\n', start);
    const bool last = end == std::string::npos;
    end = last ? all.size() : end + 1;
    number++;
    if (number == line) {
      // the line's own end kept (CRLF or LF)
      const std::string eol = last ? "" : (end >= 2 && all[end - 2] == '\r' ? "\r\n" : "\n");
      out += text + eol;
      replaced = true;
    } else {
      out += all.substr(start, end - start);
    }
    start = end;
  }
  if (!replaced) {
    why = from + " has no line " + std::to_string(line);
    return false;
  }
  std::ofstream o(to, std::ios::binary);
  o << out;
  if (!o) {
    why = "cannot write " + to;
    return false;
  }
  return true;
}

bool ConvertProject(const std::string &from, const std::string &to, std::string &report, ConvertLineProblem &problem) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  problem = ConvertLineProblem();
  // AZURE2 --convert OLD NEW: the same converter, the same refusals of an existing output
  std::ifstream in(from, std::ios::binary);
  if (!in) {
    report = "cannot read " + from;
    return false;
  }
  std::ostringstream bytes;
  bytes << in.rdbuf();
  if (from == to) {
    report = "the input and the output are the same file";
    return false;
  }
  azure2::convert::ConvertOptions o;
  o.projectDir = azure2::experiment::ProjectDirectoryOf(from);
  const azure2::convert::ConvertResult r = azure2::convert::ConvertToFormat3(bytes.str(), o);
  std::ostringstream text;
  for (const azure2::convert::ConvertNote &n : r.notes)
    text << "note: " << (n.record.empty() ? "" : n.record + ": ") << n.text << "\n";
  for (const azure2::convert::ConvertQuestion &q : r.questions)
    text << "question " << q.key << ": " << q.text << (q.proposal.empty() ? "" : " [" + q.proposal + "]") << "\n";
  if (!r.ok) {
    for (const azure2::experiment::Diagnostic &d : r.diagnostics) text << d.ToText() << "\n";
    text << from << " is not converted: " << r.refusal << "\n";
    report = text.str();
    if (r.lineProblem.Found()) {
      problem.line = r.lineProblem.fileLine;
      problem.section = r.lineProblem.section;
      problem.text = r.lineProblem.text;
      problem.problem = r.lineProblem.problem;
      problem.expected = r.lineProblem.expected;
      for (const azure2::convert::LineReading &g : r.lineProblem.readings) problem.readings.push_back({g.label, g.line});
    }
    return false;
  }
  const std::string outDir = azure2::experiment::ProjectDirectoryOf(to);
  std::vector<std::pair<std::string, std::string>> files(r.sideFiles.begin(), r.sideFiles.end());
  files.push_back({std::filesystem::path(to).filename().string(), r.text});
  std::error_code ec;
  for (const auto &f : files) {
    const std::filesystem::path p = std::filesystem::path(outDir) / f.first;
    if (std::filesystem::exists(p, ec)) {
      text << p.string() << " exists: nothing is overwritten\n";
      report = text.str();
      return false;
    }
  }
  for (const auto &f : files) {
    std::ofstream out((std::filesystem::path(outDir) / f.first).string(), std::ios::binary);
    out << f.second;
    if (!out) {
      text << "cannot write " << f.first << "\n";
      report = text.str();
      return false;
    }
  }
  text << "Wrote " << to << " (format 3): " << r.notes.size() << " note(s), " << r.questions.size() << " question(s).\n";
  report = text.str();
  return true;
}

std::string Version() {
#ifdef AZURE2_VERSION_STRING
  return AZURE2_VERSION_STRING;
#else
  return "unknown";
#endif
}

}  // namespace front
}  // namespace azure2
