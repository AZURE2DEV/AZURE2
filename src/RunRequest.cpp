// The table of the program's commands and of the questions of a run, and the command line's reader of a run request
// (include/RunRequest.h).  The menu (src/RunMenu.cpp) reads the same rows.
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

/// A row of the question table, written by name: the kind's constructor, then what differs from its defaults.
struct Row {
  RunQuestion q;
  Row(const char *option, AnswerKind kind, const char *value, const char *help) {
    q = RunQuestion{option, kind, 0u, 0u, AskedWhen::kAlways, false, -kInfinity, kInfinity, false, false, 0.,
                    "", value, help, "", 0, false, "", "", nullptr, nullptr, nullptr, nullptr};
  }
  /// The commands that take it on the command line and the commands whose menu asks it.
  Row &Commands(unsigned taken, unsigned asked) { q.commands = taken; q.menuAsks = asked; return *this; }
  Row &When(AskedWhen when) { q.when = when; return *this; }
  Row &Required() { q.required = true; return *this; }
  Row &Default(double value) { q.defaultValue = value; return *this; }
  Row &AtLeast(double low) { q.low = low; return *this; }
  Row &Above(double low) { q.low = low; q.lowOpen = true; return *this; }
  Row &AtMost(double high) { q.high = high; return *this; }
  Row &NotZero() { q.excludesZero = true; return *this; }
  Row &Words(const char *words) { q.words = words; return *this; }
  /// The menu's prompt as printed: its text, field width, a blank line before it, the name an end of input gives.
  Row &Prompt(const char *prompt, int width, const char *endName) {
    q.prompt = prompt; q.promptWidth = width; q.endName = endName; return *this;
  }
  Row &BlankLineBefore() { q.blankLineBefore = true; return *this; }
  Row &Retry(const char *retry) { q.retry = retry; return *this; }
  Row &Into(bool RunRequest::*member) { q.flag = member; return *this; }
  Row &Into(int RunRequest::*member) { q.whole = member; return *this; }
  Row &Into(double RunRequest::*member) { q.number = member; return *this; }
  Row &Into(std::string RunRequest::*member) { q.text = member; return *this; }
};

std::vector<RunQuestion> MakeQuestions() {
  using K = AnswerKind;
  using W = AskedWhen;
  const char *kPositive = "Please enter a positive number.";
  std::vector<Row> rows = {
      Row("--delta-chi2", K::kNumber, "X", "MINOS: the rise of chi-squared that bounds each interval (Minuit's Up; "
                                           "1 for one standard deviation of one parameter)")
          .Commands(kMinos, kMinos).Required().AtLeast(0.)
          .Prompt("Allowed Chi-Squared Variance: ", 30, "Allowed Chi-Squared Variance").Retry(kPositive)
          .Into(&RunRequest::deltaChi2),
      Row("--covariance-band", K::kFlag, "", "write the cross-section uncertainty band from the parameter covariance "
                                             "(this fit's, else the covariance.dat a fit saved)")
          .Commands(kCalculate | kFit | kExtrapolate | kMinos, kFit | kExtrapolate | kMinos)
          .Prompt("Calculate cross-section uncertainty band?", 0, "").Into(&RunRequest::covarianceBand),
      Row("--scale-covariance", K::kFlag, "", "with --covariance-band: inflate the covariance by the reduced "
                                              "chi-squared when it exceeds 1 (Birge ratio)")
          .Commands(kFit | kMinos, kFit | kMinos).When(W::kBandChosen)
          .Prompt("Scale covariance to reduced chi-squared = 1?", 0, "").Into(&RunRequest::scaleCovariance),
      Row("--parameter-file", K::kFile, "FILE", "start from the parameters of FILE (a param.sav), applied by name; "
                                                "without it, from the <levels> of the project")
          .Commands(kEvery, kEvery).BlankLineBefore()
          .Prompt("External Parameter File (leave blank for new file): ", 0, "").Into(&RunRequest::parameterFile),
      Row("--rate-entrance", K::kWholeNumber, "N", "the entrance pair key of the reaction rate")
          .Commands(kRate, kRate).Required().NotZero()
          .Prompt("Reaction Rate Entrance Pair: ", 38, "Reaction Rate Entrance Pair")
          .Retry("Please enter an integer greater than zero.").Into(&RunRequest::rateEntrance),
      Row("--rate-exit", K::kWholeNumber, "N", "the exit pair key of the reaction rate; -1 the total radiative "
                                               "capture of the entrance pair")
          .Commands(kRate, kRate).Required().AtLeast(-1.).NotZero()
          .Prompt("Reaction Rate Exit Pair: ", 38, "Reaction Rate Exit Pair")
          .Retry("Please enter an integer greater than zero, or -1 for the total radiative capture.")
          .Into(&RunRequest::rateExit),
      Row("", K::kWord, "", "the rate's temperatures from a file (--temperature-file) or from a minimum, maximum and "
                            "step (--temperature-min, --temperature-max, --temperature-step)")
          .Commands(kRate, kRate).Required()
          .Prompt("Use temperatures from file (yes/no): ", 38, "Use temperatures from file")
          .Retry("Please type 'yes' or 'no'.").Into(&RunRequest::temperaturesFromFile),
      Row("--temperature-file", K::kFile, "FILE", "the rate's temperatures, GK, one per line")
          .Commands(kRate, kRate).When(W::kTemperaturesFromFile).Required()
          .Prompt("Temperature File Name: ", 38, "Temperature File Name").Into(&RunRequest::temperatureFile),
      Row("--temperature-min", K::kNumber, "T", "the rate's lowest temperature, GK (without --temperature-file)")
          .Commands(kRate, kRate).When(W::kTemperaturesInRange).Required().AtLeast(0.)
          .Prompt("Reaction Rate Min Temp [GK]: ", 38, "Reaction Rate Min Temp").Retry(kPositive)
          .Into(&RunRequest::temperatureMin),
      Row("--temperature-max", K::kNumber, "T", "the rate's highest temperature, GK (without --temperature-file)")
          .Commands(kRate, kRate).When(W::kTemperaturesInRange).Required().AtLeast(0.)
          .Prompt("Reaction Rate Max Temp [GK]: ", 38, "Reaction Rate Max Temp").Retry(kPositive)
          .Into(&RunRequest::temperatureMax),
      Row("--temperature-step", K::kNumber, "T", "the rate's temperature step, GK (without --temperature-file)")
          .Commands(kRate, kRate).When(W::kTemperaturesInRange).Required().AtLeast(0.)
          .Prompt("Reaction Rate Temp Step [GK]: ", 38, "Reaction Rate Temp Step").Retry(kPositive)
          .Into(&RunRequest::temperatureStep),
      Row("--capture-amplitudes", K::kFile, "FILE", "read the external-capture amplitudes from FILE (an intEC.dat of "
                                                    "an earlier run, checked against its signature) instead of "
                                                    "computing them; used only when the run has external capture")
          .Commands(kEvery & ~kRate, kEvery & ~kRate).When(W::kProjectNeedsCapture).BlankLineBefore()
          .Prompt("External Capture Amplitude File (leave blank for new file): ", 0, "")
          .Into(&RunRequest::captureFile),
      Row("--walkers", K::kWholeNumber, "N", "the sampler's walkers, 2 to 1000")
          .Commands(kMcmc, kMcmc).Required().AtLeast(2.).AtMost(1000.)
          .Prompt("Number of Walkers (2-1000): ", 38, "Number of Walkers")
          .Retry("Please enter an integer between 2 and 1000.").Into(&RunRequest::walkers),
      Row("--steps", K::kWholeNumber, "N", "the sampler's steps, 100 to 10000000")
          .Commands(kMcmc, kMcmc).Required().AtLeast(100.).AtMost(10000000.)
          .Prompt("Number of Steps (100-10000000): ", 38, "Number of Steps")
          .Retry("Please enter an integer between 100 and 10000000.").Into(&RunRequest::steps),
      Row("--parameter-spread", K::kNumber, "P", "the walkers' initial spread, percent of each parameter, above 0 "
                                                 "and up to 50")
          .Commands(kMcmc, kMcmc).Default(5.).Above(0.).AtMost(50.)
          .Prompt("Initial Parameter Spread % (0.001-50) [default=5]: ", 38, "")
          .Retry("Please enter a number between 0.001 and 50, or press Enter for default (5%).")
          .Into(&RunRequest::parameterSpread),
      // Absolute, not a percentage: level energies scattered by a percentage of their value put walkers on unrelated
      // resonance structures, and the ensemble never contracts.
      Row("--energy-spread", K::kNumber, "E", "the walkers' initial spread of a level energy, keV, above 0 and up to 1")
          .Commands(kMcmc, kMcmc).Default(1.).Above(0.).AtMost(1.)
          .Prompt("Level Energy Spread keV (0.001-1) [default=1]: ", 38, "")
          .Retry("Please enter a number between 0.001 and 1, or press Enter for default (1 keV).")
          .Into(&RunRequest::energySpread),
      Row("--sampler-threads", K::kWholeNumber, "N", "the sampler's threads (1 or more)")
          .Commands(kMcmc, kMcmc).Required().AtLeast(1.)
          .Prompt("Number of Threads (1 or more): ", 38, "Number of Threads")
          .Retry("Please enter a positive integer.").Into(&RunRequest::samplerThreads),
      Row("--basis", K::kWord, "rwa|physical", "the sampled parameters: reduced width amplitudes (rwa) or the "
                                               "project's physical parameters")
          .Commands(kMcmc, kMcmc).Required().Words("rwa|physical")
          .Prompt("Use Reduced Width Amplitudes (yes/no): ", 38, "Use Reduced Width Amplitudes")
          .Retry("Please type 'yes' or 'no'.").Into(&RunRequest::reducedWidthAmplitudes),
      Row("--overwrite-samples", K::kWord, "", "overwrite samples.mcmc; without it the sampler appends to the file")
          .Commands(kMcmc, kMcmc)
          .Prompt("Overwrite existing samples.mcmc? (yes/no) [default=no]: ", 38, "")
          .Retry("Please type 'yes' or 'no', or press Enter for default (no).").Into(&RunRequest::overwriteSamples),
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
  // a kWord row without words is a flag on the command line (its default is "no")
  const bool takesValue = q.kind != AnswerKind::kFlag && !(q.kind == AnswerKind::kWord && !*q.words);
  return std::string(q.option) + (takesValue ? std::string("=") + q.value : "");
}

bool TakesValue(const RunQuestion &q) { return Written(q) != q.option; }

/// "--a, --b=N" of the options `command` takes.
std::string OptionsOf(RunCommand command) {
  std::string list;
  for (const RunQuestion &q : RunQuestions())
    if (*q.option && (q.commands & CommandBit(command))) list += (list.empty() ? "" : ", ") + Written(q);
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
      // the menu's reading of a number (operator>>), and nothing after it
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
      {RunCommand::kCalculate, "calculate", 1, "Calculate Segments From Data",
       "calculate the data sets of the project at their data points from the parameters (no fit)", true},
      {RunCommand::kFit, "fit", 2, "Fit Segments From Data",
       "fit the parameters to the data sets, then calculate them (param.sav, covariance.dat)", true},
      {RunCommand::kExtrapolate, "extrapolate", 3, "Calculate Segments Without Data",
       "calculate the extrapolations of the project (no data)", true},
      {RunCommand::kMinos, "minos", 4, "Perform MINOS Error Analysis",
       "fit, then the MINOS interval of every free parameter at --delta-chi2", true},
      {RunCommand::kRate, "rate", 5, "Calculate Reaction Rate",
       "the reaction rate of one entrance and exit pair over a grid of temperatures", true},
      {RunCommand::kMcmc, "mcmc", 6, "Perform MCMC Bayesian Sampling",
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
    if (*q.option && option == q.option) return &q;
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
                                  "them); without a command, give the project file alone for the menu.");
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
  std::set<const RunQuestion *> given;
  for (std::size_t i = 0; i < arguments.size(); i++) {
    const std::string &argument = arguments[i];
    if (i == words.front()) continue;
    if (argument.compare(0, 2, "--") != 0) {
      if (!project.empty())
        Refuse("Two project files: " + project + " and " + argument + ".", command,
                                            std::string("Give one: AZURE2 ") + row->name + " [options] FILE.");
      project = argument;
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
  // the menu's "temperatures from a file?" is answered by --temperature-file
  request.temperaturesFromFile = given.count(FindRunQuestion("--temperature-file")) != 0;
  std::string missing;
  for (const RunQuestion &q : RunQuestions()) {
    if (!*q.option || !(q.commands & CommandBit(command))) continue;
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
      Refuse(std::string(rule.first) + " and " + rule.last + ": " + rule.menuMessage,
                                          command, rule.fix);
  }
  rest.push_back(project);
  return true;
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
  out << "Commands: AZURE2 <command> [options] FILE runs one calculation and reads nothing from standard input"
      << std::endl
      << "(without a command, AZURE2 FILE asks for the command and its answers in a menu):" << std::endl;
  for (const RunCommandRow &c : RunCommands()) {
    if (!c.built) continue;
    std::string head = std::string("\t") + c.name + ":";
    out << head << std::string(head.size() < 25 ? 25 - head.size() : 1, ' ') << c.help << " (menu " << c.menuNumber
        << ")." << std::endl
        << std::string(25, ' ') << "Options: " << OptionsOf(c.command) << "." << std::endl;
  }
  out << "The options of the commands (each is refused by a command that does not take it):" << std::endl;
  for (const RunQuestion &q : RunQuestions()) {
    if (!*q.option || CommandsOf(q).empty()) continue;
    std::string text = std::string(TakesValue(q) ? std::string("=") + q.value + ": " : "") + q.help;
    if (q.kind == AnswerKind::kNumber || q.kind == AnswerKind::kWholeNumber) text += "; " + RangeText(q);
    if (q.required)
      text += "; required";
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
    out << (firstCommand ? "\n  " : ",\n  ") << "{\"name\": " << JsonText(c.name) << ", \"menu\": " << c.menuNumber
        << ", \"help\": " << JsonText(c.help) << ", \"options\": [";
    firstCommand = false;
    bool firstOption = true;
    for (const RunQuestion &q : RunQuestions()) {
      if (!*q.option || !(q.commands & CommandBit(c.command))) continue;
      const bool number = q.kind == AnswerKind::kNumber || q.kind == AnswerKind::kWholeNumber;
      out << (firstOption ? "\n   " : ",\n   ") << "{\"option\": " << JsonText(q.option)
          << ", \"written\": " << JsonText(Written(q))
          << ", \"kind\": " << JsonText(TakesValue(q) ? KindName(q.kind) : "flag")
          << ", \"required\": " << (q.required ? "true" : "false")
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
