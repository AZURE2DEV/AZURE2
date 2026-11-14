// The table of the program's commands and of the options of a run, and the command line's reader of a run request
// (include/RunRequest.h).
#include "RunRequest.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

constexpr unsigned kCalculate = CommandBit(RunCommand::kCalculate);
constexpr unsigned kFit = CommandBit(RunCommand::kFit);
constexpr unsigned kExtrapolate = CommandBit(RunCommand::kExtrapolate);
constexpr unsigned kMinos = CommandBit(RunCommand::kMinos);
constexpr unsigned kRate = CommandBit(RunCommand::kRate);
constexpr unsigned kMcmc = CommandBit(RunCommand::kMcmc);
constexpr unsigned kEvery = kCalculate | kFit | kExtrapolate | kMinos | kRate | kMcmc;

/// A row of the option table, written by name: the kind's constructor, then what differs from its defaults.
struct Row {
  RunQuestion q;
  Row(const char *option, AnswerKind kind, const char *value, const char *help) {
    q = RunQuestion{option, kind, 0u, AskedWhen::kAlways, false, "", -kInfinity, kInfinity, false, false, 0., "", value,
                    help, nullptr, nullptr, nullptr, nullptr};
  }
  /// The commands that take it.
  Row &Commands(unsigned taken) { q.commands = taken; return *this; }
  Row &When(AskedWhen when) { q.when = when; return *this; }
  Row &Required() { q.required = true; return *this; }
  /// Required, or the other way `other` (help names both).
  Row &RequiredOr(const char *other) { q.required = true; q.instead = other; return *this; }
  Row &Default(double value) { q.defaultValue = value; return *this; }
  Row &AtLeast(double low) { q.low = low; return *this; }
  Row &Above(double low) { q.low = low; q.lowOpen = true; return *this; }
  Row &AtMost(double high) { q.high = high; return *this; }
  Row &NotZero() { q.excludesZero = true; return *this; }
  Row &Words(const char *words) { q.words = words; return *this; }
  Row &Into(bool RunRequest::*member) { q.flag = member; return *this; }
  Row &Into(int RunRequest::*member) { q.whole = member; return *this; }
  Row &Into(double RunRequest::*member) { q.number = member; return *this; }
  Row &Into(std::string RunRequest::*member) { q.text = member; return *this; }
};

std::vector<RunQuestion> MakeQuestions() {
  using K = AnswerKind;
  using W = AskedWhen;
  const char *kRange = "--temperature-min, --temperature-max and --temperature-step";
  std::vector<Row> rows = {
      Row("--delta-chi2", K::kNumber, "X", "MINOS: the rise of chi-squared that bounds each interval (Minuit's Up; "
                                           "1 for one standard deviation of one parameter)")
          .Commands(kMinos).Required().Above(0.).Into(&RunRequest::deltaChi2),
      Row("--covariance-band", K::kFlag, "", "write the cross-section uncertainty band from the parameter covariance "
                                             "(this fit's, else the covariance.dat a fit saved)")
          .Commands(kCalculate | kFit | kExtrapolate | kMinos).Into(&RunRequest::covarianceBand),
      Row("--scale-covariance", K::kFlag, "", "with --covariance-band: inflate the covariance by the reduced "
                                              "chi-squared when it exceeds 1 (Birge ratio)")
          .Commands(kFit | kMinos).When(W::kBandChosen).Into(&RunRequest::scaleCovariance),
      Row("--parameter-file", K::kFile, "FILE", "start from the parameters of FILE (a param.sav), applied by name; "
                                                "without it, from the <levels> of the project")
          .Commands(kEvery).Into(&RunRequest::parameterFile),
      Row("--rate-entrance", K::kWholeNumber, "N", "the entrance pair key of the reaction rate")
          .Commands(kRate).Required().AtLeast(1.).Into(&RunRequest::rateEntrance),
      Row("--rate-exit", K::kWholeNumber, "N", "the exit pair key of the reaction rate; -1 the total radiative "
                                               "capture of the entrance pair")
          .Commands(kRate).Required().AtLeast(-1.).NotZero().Into(&RunRequest::rateExit),
      Row("--temperature-file", K::kFile, "FILE", "the rate's temperatures, GK, one per line")
          .Commands(kRate).When(W::kTemperaturesFromFile).RequiredOr(kRange).Into(&RunRequest::temperatureFile),
      Row("--temperature-min", K::kNumber, "T", "the rate's lowest temperature, GK")
          .Commands(kRate).When(W::kTemperaturesInRange).RequiredOr("--temperature-file").AtLeast(0.)
          .Into(&RunRequest::temperatureMin),
      Row("--temperature-max", K::kNumber, "T", "the rate's highest temperature, GK")
          .Commands(kRate).When(W::kTemperaturesInRange).RequiredOr("--temperature-file").AtLeast(0.)
          .Into(&RunRequest::temperatureMax),
      Row("--temperature-step", K::kNumber, "T", "the rate's temperature step, GK")
          .Commands(kRate).When(W::kTemperaturesInRange).RequiredOr("--temperature-file").AtLeast(0.)
          .Into(&RunRequest::temperatureStep),
      Row("--capture-amplitudes", K::kFile, "FILE", "read the external-capture amplitudes from FILE (an intEC.dat of "
                                                    "an earlier run, checked against its signature) instead of "
                                                    "computing them; used only when the run has external capture")
          .Commands(kEvery & ~kRate).When(W::kProjectNeedsCapture).Into(&RunRequest::captureFile),
      Row("--walkers", K::kWholeNumber, "N", "the sampler's walkers")
          .Commands(kMcmc).Required().AtLeast(2.).AtMost(1000.).Into(&RunRequest::walkers),
      Row("--steps", K::kWholeNumber, "N", "the sampler's steps")
          .Commands(kMcmc).Required().AtLeast(100.).AtMost(10000000.).Into(&RunRequest::steps),
      Row("--parameter-spread", K::kNumber, "P", "the walkers' initial spread, percent of each parameter")
          .Commands(kMcmc).Default(5.).Above(0.).AtMost(50.).Into(&RunRequest::parameterSpread),
      // Absolute, not a percentage: level energies scattered by a percentage of their value put walkers on unrelated
      // resonance structures, and the ensemble never contracts.
      Row("--energy-spread", K::kNumber, "E", "the walkers' initial spread of a level energy, keV")
          .Commands(kMcmc).Default(1.).Above(0.).AtMost(1.).Into(&RunRequest::energySpread),
      Row("--sampler-threads", K::kWholeNumber, "N", "the sampler's threads")
          .Commands(kMcmc).Required().AtLeast(1.).Into(&RunRequest::samplerThreads),
      Row("--basis", K::kWord, "rwa|physical", "the sampled parameters: reduced width amplitudes (rwa) or the "
                                               "project's physical parameters")
          .Commands(kMcmc).Required().Words("rwa|physical").Into(&RunRequest::reducedWidthAmplitudes),
      Row("--overwrite-samples", K::kFlag, "", "overwrite samples.mcmc; without it the sampler appends to the file")
          .Commands(kMcmc).Into(&RunRequest::overwriteSamples),
  };
  std::vector<RunQuestion> questions;
  for (const Row &r : rows) questions.push_back(r.q);
  return questions;
}

/// Refuses the command line (CommandRefusal's three parts).
[[noreturn]] void Refuse(const std::string &what, RunCommand command, const std::string &fix) {
  throw std::invalid_argument(CommandRefusal(what, command, fix));
}

bool BothPairsDiffer(const RunRequest &r) { return r.rateEntrance != r.rateExit; }

/// The option of `q` as written with its value ("--walkers=N", "--covariance-band").
std::string Written(const RunQuestion &q) {
  return std::string(q.option) + (q.kind != AnswerKind::kFlag ? std::string("=") + q.value : "");
}

bool TakesValue(const RunQuestion &q) { return q.kind != AnswerKind::kFlag; }

/// "--a, --b=N" of the options `command` takes.
std::string OptionsOf(RunCommand command) {
  std::string list;
  for (const RunQuestion &q : RunQuestions())
    if (q.commands & CommandBit(command)) list += (list.empty() ? "" : ", ") + Written(q);
  return list.empty() ? std::string("no options") : list;
}

/// The commands that take `q`, by name.
std::string CommandsOf(const RunQuestion &q) {
  std::string list;
  for (const RunCommandRow &c : RunCommands())
    if (c.built && (q.commands & CommandBit(c.command))) list += (list.empty() ? "" : ", ") + std::string(c.name);
  return list;
}

std::string Number(double x) {
  std::ostringstream s;
  s << std::setprecision(12) << x;
  return s.str();
}

/// What a number question takes, in words ("a whole number from 2 to 1000").
std::string RangeText(const RunQuestion &q) {
  std::string text = q.kind == AnswerKind::kWholeNumber ? "a whole number" : "a number";
  const bool low = std::isfinite(q.low), high = std::isfinite(q.high);
  if (low && high && !q.lowOpen)
    text += " from " + Number(q.low) + " to " + Number(q.high);
  else if (low && high)
    text += " above " + Number(q.low) + " and up to " + Number(q.high);
  else if (low)
    text += (q.lowOpen ? " above " : " of at least ") + Number(q.low);
  if (q.excludesZero) text += " other than 0";
  return text;
}

/// The value of `argument` ("--name=value") as the answer of `q` into `request`, or the refusal.
void ReadValue(const RunQuestion &q, const std::string &argument, RunCommand command, RunRequest &request) {
  const std::size_t eq = argument.find('=');
  const bool hasValue = eq != std::string::npos;
  const std::string value = hasValue ? argument.substr(eq + 1) : "";
  if (!TakesValue(q)) {
    if (hasValue)
      Refuse(argument + ": " + q.option + " takes no value.", command,
                                          std::string("Give it alone: ") + q.option + ".");
    request.*q.flag = true;
    return;
  }
  if (value.empty())
    Refuse(std::string(q.option) + " needs a value: " + Written(q) + ".", command,
                                        "Give it as " + Written(q) + " (" + q.help + ").");
  switch (q.kind) {
    case AnswerKind::kWord: {
      const std::string words = q.words;
      const std::string yes = words.substr(0, words.find('|')), no = words.substr(words.find('|') + 1);
      if (value != yes && value != no)
        Refuse(argument + ": " + q.option + " takes " + yes + " or " + no + ".", command,
                                            "Give " + std::string(q.option) + "=" + yes + " or " + q.option + "=" +
                                                no + ".");
      request.*q.flag = value == yes;
      return;
    }
    case AnswerKind::kNumber:
    case AnswerKind::kWholeNumber: {
      // a number as operator>> reads it, and nothing after it
      std::istringstream in(value);
      double x = 0.;
      int n = 0;
      const bool read = q.kind == AnswerKind::kWholeNumber ? static_cast<bool>(in >> n) : static_cast<bool>(in >> x);
      if (q.kind == AnswerKind::kWholeNumber) x = n;
      std::string rest;
      if (!read || (in >> rest) || !InRange(q, x))
        Refuse(argument + ": " + q.option + " takes " + RangeText(q) + ".", command,
                                            "Give " + Written(q) + " with " + q.value + " " + RangeText(q) + ".");
      if (q.kind == AnswerKind::kWholeNumber)
        request.*q.whole = n;
      else
        request.*q.number = x;
      return;
    }
    case AnswerKind::kFile:
      request.*q.text = value;  // read when the project is read: relative to its directory (RequestFileRefusal)
      return;
    case AnswerKind::kFlag:
      return;
  }
}

std::string JsonText(const std::string &text) {
  std::string out = "\"";
  for (const char c : text) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + "\"";
}

std::string JsonNumber(double x) { return std::isfinite(x) ? Number(x) : "null"; }

const char *KindName(AnswerKind kind) {
  switch (kind) {
    case AnswerKind::kFlag: return "flag";
    case AnswerKind::kWord: return "word";
    case AnswerKind::kNumber: return "number";
    case AnswerKind::kWholeNumber: return "whole number";
    case AnswerKind::kFile: return "file";
  }
  return "";
}

const char *WhenName(AskedWhen when) {
  switch (when) {
    case AskedWhen::kAlways: return "always";
    case AskedWhen::kBandChosen: return "with --covariance-band";
    case AskedWhen::kTemperaturesFromFile: return "temperatures from a file";
    case AskedWhen::kTemperaturesInRange: return "without --temperature-file";
    case AskedWhen::kProjectNeedsCapture: return "when the run has external capture";
  }
  return "";
}

}  // namespace

std::string CommandRefusal(const std::string &what, RunCommand command, const std::string &fix) {
  const RunCommandRow &c = RunCommandOf(command);
  return what + "\n       " + c.name + " takes: " + OptionsOf(command) + ".\n       " + fix;
}

const std::vector<RunCommandRow> &RunCommands() {
#ifdef USE_MCMC
  constexpr bool kSampler = true;
#else
  constexpr bool kSampler = false;
#endif
  static const std::vector<RunCommandRow> commands = {
      {RunCommand::kCalculate, "calculate",
       "calculate the data sets of the project at their data points from the parameters (no fit)", true},
      {RunCommand::kFit, "fit",
       "fit the parameters to the data sets, then calculate them (param.sav, covariance.dat)", true},
      {RunCommand::kExtrapolate, "extrapolate",
       "calculate the extrapolations of the project (no data)", true},
      {RunCommand::kMinos, "minos",
       "fit, then the MINOS interval of every free parameter at --delta-chi2", true},
      {RunCommand::kRate, "rate",
       "the reaction rate of one entrance and exit pair over a grid of temperatures", true},
      {RunCommand::kMcmc, "mcmc",
       "sample the posterior of the parameters with the ensemble sampler (samples.mcmc)", kSampler},
  };
  return commands;
}

const RunCommandRow &RunCommandOf(RunCommand command) {
  for (const RunCommandRow &c : RunCommands())
    if (c.command == command) return c;
  return RunCommands().front();
}

const std::vector<RunQuestion> &RunQuestions() {
  static const std::vector<RunQuestion> questions = MakeQuestions();
  return questions;
}

std::string OptionName(const std::string &argument) { return argument.substr(0, argument.find('=')); }

const RunQuestion *FindRunQuestion(const std::string &option) {
  for (const RunQuestion &q : RunQuestions())
    if (option == q.option) return &q;
  return nullptr;
}

const std::vector<RunRule> &RunRules() {
  static const std::vector<RunRule> rules = {
      {"--rate-entrance", "--rate-exit", &BothPairsDiffer, "Cannot calculate rate for elastic scattering.",
       "Give two different pair keys: --rate-entrance=N --rate-exit=M (-1: the total radiative capture)."},
  };
  return rules;
}

bool QuestionApplies(const RunQuestion &question, const RunRequest &request) {
  switch (question.when) {
    case AskedWhen::kAlways: return true;
    case AskedWhen::kBandChosen: return request.covarianceBand;
    case AskedWhen::kTemperaturesFromFile: return request.temperaturesFromFile;
    case AskedWhen::kTemperaturesInRange: return !request.temperaturesFromFile;
    case AskedWhen::kProjectNeedsCapture: return true;
  }
  return true;
}

bool InRange(const RunQuestion &question, double value) {
  if (question.excludesZero && value == 0.) return false;
  if (question.lowOpen ? !(value > question.low) : !(value >= question.low)) return false;
  return value <= question.high;
}

bool ReadRunCommand(const std::vector<std::string> &arguments, RunRequest &request, std::vector<std::string> &rest) {
  rest.clear();
  std::vector<std::size_t> words;  // the arguments that are not options
  for (std::size_t i = 0; i < arguments.size(); i++)
    if (arguments[i].compare(0, 2, "--") != 0) words.push_back(i);
  const RunCommandRow *row = nullptr;
  if (!words.empty())
    for (const RunCommandRow &c : RunCommands())
      if (arguments[words.front()] == c.name) row = &c;
  std::string names;
  for (const RunCommandRow &c : RunCommands())
    if (c.built) names += (names.empty() ? "" : ", ") + std::string(c.name);
  if (!row) {
    if (words.size() > 1)
      throw std::invalid_argument(arguments[words.front()] + " is not a command of AZURE2 (and a run takes one " +
                                  "project file).\n       The commands: " + names + ".\n       Give one of them " +
                                  "before the project file: AZURE2 <command> [options] FILE (AZURE2 --help lists " +
                                  "them).");
    rest = arguments;
    return false;
  }
  if (!row->built)
    throw std::invalid_argument(std::string(row->name) + " is not a command of this build of AZURE2 (it has no " +
                                "sampler).\n       The commands of this build: " + names + ".\n       Build AZURE2 " +
                                "with the sampler: cmake -DUSE_MCMC=ON.");
  const RunCommand command = row->command;
  request = RunRequest();
  request.command = command;
  std::string project;
  std::size_t projectAt = 0;  // its place among the arguments
  std::set<const RunQuestion *> given;
  for (std::size_t i = 0; i < arguments.size(); i++) {
    const std::string &argument = arguments[i];
    if (i == words.front()) continue;
    if (argument.compare(0, 2, "--") != 0) {
      if (!project.empty()) {
        // "--max-calls 50 FILE": the first word is an option's value written after a space
        const std::string &before = projectAt > 0 ? arguments[projectAt - 1] : std::string();
        const bool value = before.compare(0, 2, "--") == 0 && before.find('=') == std::string::npos;
        Refuse("Two project files: " + project + " and " + argument + "." +
                   (value ? " " + before + " needs a value: " + before + "=" + project + " (with '=')." : ""),
               command, std::string("Give one: AZURE2 ") + row->name + " [options] FILE.");
      }
      project = argument;
      projectAt = i;
      continue;
    }
    const RunQuestion *q = FindRunQuestion(OptionName(argument));
    if (!q) {
      rest.push_back(argument);  // an option of the program (parseOptions)
      continue;
    }
    if (!(q->commands & CommandBit(command)))
      Refuse(std::string(q->option) + " is not an option of the command " + row->name + ".", command,
             "Leave it out, or give it to a command that takes it: " + CommandsOf(*q) + ".");
    if (!given.insert(q).second)
      Refuse(std::string(q->option) + " is given twice.", command, "Give it once.");
    ReadValue(*q, argument, command, request);
  }
  if (project.empty())
    Refuse(std::string("The command ") + row->name + " needs the project file.", command,
                                        std::string("Give it last: AZURE2 ") + row->name + " [options] FILE.");
  // the temperatures come from a file when --temperature-file is given, else from the minimum, maximum and step
  request.temperaturesFromFile = given.count(FindRunQuestion("--temperature-file")) != 0;
  std::string missing;
  for (const RunQuestion &q : RunQuestions()) {
    if (!(q.commands & CommandBit(command))) continue;
    const bool applies = QuestionApplies(q, request);
    if (given.count(&q) && !applies)
      Refuse(std::string(q.option) + " is used only " + WhenName(q.when) + ".", command,
                                          std::string("Leave ") + q.option + " out, or give what it goes with.");
    if (given.count(&q) || !applies) continue;
    if (q.required)
      missing += (missing.empty() ? "" : ", ") + Written(q);
    else if (q.number)
      request.*q.number = q.defaultValue;
  }
  if (!missing.empty())
    Refuse(std::string("The command ") + row->name + " needs " + missing + ".", command,
                                        std::string("Add ") + (missing.find(',') == std::string::npos ? "it" : "them") +
                                            ": AZURE2 " + row->name + " " + missing + " [options] FILE" +
                                            (command == RunCommand::kRate
                                                 ? " (the temperatures: --temperature-file=FILE, or "
                                                   "--temperature-min, --temperature-max and --temperature-step)."
                                                 : "."));
  for (const RunRule &rule : RunRules()) {
    const RunQuestion *first = FindRunQuestion(rule.first), *last = FindRunQuestion(rule.last);
    if (!(first->commands & last->commands & CommandBit(command))) continue;
    if (!rule.holds(request))
      Refuse(std::string(rule.first) + " and " + rule.last + ": " + rule.message,
                                          command, rule.fix);
  }
  rest.push_back(project);
  return true;
}

std::string NoCommandRefusal() {
  std::string names;
  for (const RunCommandRow &c : RunCommands())
    if (c.built) names += (names.empty() ? "" : ", ") + std::string(c.name);
  return "No command: AZURE2 runs a calculation only when a command names it (it reads nothing from standard "
         "input).\n       The commands: " + names + ".\n       Give one before the project file: AZURE2 <command> "
         "[options] FILE (AZURE2 --help lists the commands and their options; tools/menu_to_command.py turns the "
         "answers of an old menu script into the command).";
}

std::string CommandOptionWithoutCommand(const RunQuestion &question) {
  return std::string(question.option) + " is an option of a command, and no command is given.\n       The commands "
         "that take it: " + CommandsOf(question) +
         ".\n       Give the command before the project file: AZURE2 <command> " +
         Written(question) + " FILE (AZURE2 --help lists the commands).";
}

int CommandError(const std::string &message) {
  std::cout << "ERROR: " << message << std::endl
            << "\tSyntax: AZURE2 <command> [options] FILE  (AZURE2 --help lists the commands and their options)"
            << std::endl;
  return 2;
}

std::string RunCommandsHelp() {
  std::ostringstream out;
  out << "Syntax: AZURE2 <command> [options] FILE" << std::endl
      << std::endl
      << "Commands (each runs one calculation of the project FILE and reads nothing from standard input):"
      << std::endl;
  for (const RunCommandRow &c : RunCommands()) {
    if (!c.built) continue;
    std::string head = std::string("\t") + c.name + ":";
    out << head << std::string(head.size() < 25 ? 25 - head.size() : 1, ' ') << c.help << "." << std::endl
        << std::string(25, ' ') << "Options: " << OptionsOf(c.command) << "." << std::endl;
  }
  out << "The options of the commands (each is refused by a command that does not take it):" << std::endl;
  for (const RunQuestion &q : RunQuestions()) {
    if (CommandsOf(q).empty()) continue;
    std::string text = std::string(TakesValue(q) ? std::string("=") + q.value + ": " : "") + q.help;
    if (q.kind == AnswerKind::kNumber || q.kind == AnswerKind::kWholeNumber) text += "; " + RangeText(q);
    if (q.required)
      text += *q.instead ? std::string("; required, or instead ") + q.instead : std::string("; required");
    else if (q.number)
      text += "; default " + Number(q.defaultValue);
    text += " (" + CommandsOf(q) + ").";
    out << std::setw(25) << std::left << ("\t" + std::string(q.option) + ":") << std::setw(0) << text << std::endl;
  }
  return out.str();
}

std::string RunCommandsJson() {
  std::ostringstream out;
  out << "[";
  bool firstCommand = true;
  for (const RunCommandRow &c : RunCommands()) {
    if (!c.built) continue;
    out << (firstCommand ? "\n  " : ",\n  ") << "{\"name\": " << JsonText(c.name)
        << ", \"help\": " << JsonText(c.help) << ", \"options\": [";
    firstCommand = false;
    bool firstOption = true;
    for (const RunQuestion &q : RunQuestions()) {
      if (!(q.commands & CommandBit(c.command))) continue;
      const bool number = q.kind == AnswerKind::kNumber || q.kind == AnswerKind::kWholeNumber;
      out << (firstOption ? "\n   " : ",\n   ") << "{\"option\": " << JsonText(q.option)
          << ", \"written\": " << JsonText(Written(q))
          << ", \"kind\": " << JsonText(TakesValue(q) ? KindName(q.kind) : "flag")
          << ", \"required\": " << (q.required ? "true" : "false")
          << ", \"instead\": " << (*q.instead ? JsonText(q.instead) : "null")
          << ", \"default\": " << (q.number && !q.required ? JsonNumber(q.defaultValue) : "null")
          << ", \"low\": " << (number ? JsonNumber(q.low) : "null")
          << ", \"low_included\": " << (number && !q.lowOpen ? "true" : "false")
          << ", \"high\": " << (number ? JsonNumber(q.high) : "null")
          << ", \"excludes_zero\": " << (q.excludesZero ? "true" : "false")
          << ", \"when\": " << JsonText(WhenName(q.when)) << ", \"help\": " << JsonText(q.help) << "}";
      firstOption = false;
    }
    out << "]}";
  }
  out << "\n ]";
  return out.str();
}
