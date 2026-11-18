#ifndef RUNREQUEST_H
#define RUNREQUEST_H

// What a run of the program is asked to do, and the one table of its options (docs/source/reference/
// command_line.rst, "Commands"): the command line gives a request as `AZURE2 <command> [options] FILE`
// (ReadRunCommand), and PrepareRun (src/RunPreparation.cpp) turns it into the run's Config.  The program reads nothing
// from standard input.  A new command is a row of RunCommands() and its case in SelectCalculation
// (src/RunPreparation.cpp); a new option is a member of RunRequest and a row of RunQuestions().

#include <ostream>
#include <string>
#include <vector>

/// The calculation a run does: the commands of the command line.
enum class RunCommand { kCalculate, kFit, kExtrapolate, kMinos, kRate, kMcmc };

/// Every answer of a run.  Paths are as given (relative to the project file's directory, Config::OpenPath).
struct RunRequest {
  RunCommand command = RunCommand::kCalculate;
  bool covarianceBand = false;      ///< write the cross-section uncertainty band
  bool scaleCovariance = false;     ///< inflate the band's covariance to a reduced chi^2 of 1 (fit, MINOS)
  double deltaChi2 = 0.;            ///< MINOS: the rise of chi^2 that bounds an interval (Minuit's Up)
  std::string parameterFile;        ///< the parameters to start from; empty: the <levels> of the project
  std::string captureFile;          ///< external-capture amplitudes to read; empty: compute them (intEC.dat)
  int rateEntrance = 0;             ///< the rate's entrance pair key
  int rateExit = 0;                 ///< the rate's exit pair key; -1 the total radiative capture
  bool temperaturesFromFile = false;  ///< --temperature-file given (else the minimum, maximum and step)
  std::string temperatureFile;      ///< one temperature per line, GK
  double temperatureMin = 0.;       ///< GK
  double temperatureMax = 0.;       ///< GK
  double temperatureStep = 0.;      ///< GK
  int walkers = 0;                  ///< the sampler's walkers
  int steps = 0;                    ///< the sampler's steps
  double parameterSpread = 0.;      ///< the walkers' initial spread, percent of each parameter
  double energySpread = 0.;         ///< the walkers' initial spread of a level energy, keV
  int samplerThreads = 0;           ///< the sampler's threads
  bool reducedWidthAmplitudes = false;  ///< sample reduced width amplitudes (else the physical parameters)
  bool overwriteSamples = false;    ///< overwrite samples.mcmc (else the sampler continues it)
};

/// One command: its name on the command line and what it computes.
struct RunCommandRow {
  RunCommand command;
  const char *name;       ///< `AZURE2 <name> FILE`
  const char *help;       ///< what it computes, for --help and --describe
  bool built;             ///< false: this build lacks it (mcmc without the sampler)
};

/// The commands, in the order --help lists them.
const std::vector<RunCommandRow> &RunCommands();
/// The row of `command`.
const RunCommandRow &RunCommandOf(RunCommand command);

/// The kind of an answer: how the command line gives it.
enum class AnswerKind {
  kFlag,         ///< a flag without a value: given, yes; left out, no
  kWord,         ///< one of two words (the row's words): the first is yes, the second no
  kNumber,       ///< a real number in the row's range
  kWholeNumber,  ///< a whole number in the row's range
  kFile,         ///< a file that can be read
};

/// When an option of its command applies.
enum class AskedWhen {
  kAlways,
  kBandChosen,             ///< with the uncertainty band
  kTemperaturesFromFile,   ///< the rate's temperatures from a file
  kTemperaturesInRange,    ///< the rate's temperatures from a minimum, maximum and step
  kProjectNeedsCapture,    ///< the project's run has external capture (decided when the project is read)
};

/// One option of a run: its name, what it takes, which commands take it, and where the answer goes (exactly one of
/// the four members is set, the one of the row's kind).
struct RunQuestion {
  const char *option;      ///< "--parameter-file"
  AnswerKind kind;
  unsigned commands;       ///< the commands that take it (bits CommandBit)
  AskedWhen when;
  bool required;           ///< no default: a command line without it (where it applies) is refused
  const char *instead;     ///< for a required option that is one of two ways: the other way ("" for none), for help
  double low;              ///< range of a number (kNumber, kWholeNumber)
  double high;
  bool lowOpen;            ///< the range excludes `low`
  bool excludesZero;       ///< 0 is no answer (a pair key)
  double defaultValue;     ///< a number's value when the answer is left blank (not `required`)
  const char *words;       ///< kWord on the command line: "yes-word|no-word"
  const char *value;       ///< the value's name in --help ("FILE", "N", "X")
  const char *help;        ///< what it is; its range and default are the row's (printed from them)
  // the answer
  bool RunRequest::*flag;
  int RunRequest::*whole;
  double RunRequest::*number;
  std::string RunRequest::*text;
};

/// The bit of `command` in RunQuestion::commands.
constexpr unsigned CommandBit(RunCommand command) { return 1u << static_cast<unsigned>(command); }

/// The options, in the order --help lists them.
const std::vector<RunQuestion> &RunQuestions();
/// The option's name: "--name" of "--name=VALUE".
std::string OptionName(const std::string &argument);
/// The question whose option is `option` ("--name", without a value), or null.
const RunQuestion *FindRunQuestion(const std::string &option);

/// A rule between two answers: a command line that breaks it is refused.
struct RunRule {
  const char *first;       ///< the options the rule relates
  const char *last;
  bool (*holds)(const RunRequest &request);
  const char *message;     ///< what is wrong
  const char *fix;         ///< how to fix the command line
};
const std::vector<RunRule> &RunRules();

/// Whether `question` applies to `request` (its condition; kProjectNeedsCapture counts as applying: the project
/// decides, PrepareRun).
bool QuestionApplies(const RunQuestion &question, const RunRequest &request);
/// Whether `value` is an answer of the number question `question` (its range).
bool InRange(const RunQuestion &question, double value);

/// The command line's reader: whether `arguments` (argv without the program's name) give a command; if so, `request`
/// holds its answers (defaults where an optional answer is not given) and `rest` the arguments that are not the
/// command's (the program's other options and the project file); if not, `rest` is `arguments` (the caller refuses a
/// run without a command, NoCommandRefusal).  Refuses (std::invalid_argument, the message in
/// three parts: what, what the command takes, how to fix) an unknown or unbuilt command, an option of another command,
/// a value that does not parse or is out of range, a missing required answer, a broken rule and a missing project.
bool ReadRunCommand(const std::vector<std::string> &arguments, RunRequest &request, std::vector<std::string> &rest);
/// The three lines of a refusal of a command line: what, what `command` takes, how to fix it.
std::string CommandRefusal(const std::string &what, RunCommand command, const std::string &fix);
/// The refusal of a run without a command, in three parts.
std::string NoCommandRefusal();
/// The refusal of a command option given without a command (for parseOptions).
std::string CommandOptionWithoutCommand(const RunQuestion &question);
/// A refused command line: "ERROR: <message>" and the commands' usage line on standard output; returns 2.
int CommandError(const std::string &message);
/// The same on `out` (the run's stream of a front end in its own process, ProgramMain.h).
int CommandError(std::ostream &out, const std::string &message);

/// The commands and their options for --help.
std::string RunCommandsHelp();
/// The commands and their options as JSON (the member "commands" of --describe --json).
std::string RunCommandsJson();

#endif
