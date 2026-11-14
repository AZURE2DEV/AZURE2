#ifndef RUNREQUEST_H
#define RUNREQUEST_H

// What a run of the program is asked to do, and the one table of its questions (docs/source/reference/
// command_line.rst, "Commands"): the command line reads a request from `AZURE2 <command> [options] FILE`
// (ReadRunCommand), the text menu fills one by asking (src/RunMenu.cpp), and PrepareRun (src/RunPreparation.cpp) turns
// it into the run's Config.  A new command is a row of RunCommands() and its case in PrepareRun; a new question is a
// member of RunRequest and a row of RunQuestions().

#include <string>
#include <vector>

/// The calculation a run does: the commands of the command line, the modes 1 to 6 of the menu.
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
  bool temperaturesFromFile = false;
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

/// One command: its name on the command line, its number and line in the menu.
struct RunCommandRow {
  RunCommand command;
  const char *name;       ///< `AZURE2 <name> FILE`
  int menuNumber;         ///< its number in the menu (1 to 6)
  const char *menuText;   ///< its line in the menu, as printed
  const char *help;       ///< what it computes, for --help and --describe
  bool built;             ///< false: this build lacks it (mcmc without the sampler)
};

/// The commands, in the menu's order.
const std::vector<RunCommandRow> &RunCommands();
/// The row of `command`.
const RunCommandRow &RunCommandOf(RunCommand command);

/// The kind of an answer: how the menu asks it and how the command line gives it.
enum class AnswerKind {
  kFlag,         ///< a yes/no question whose blank answer is no; on the command line a flag without a value
  kWord,         ///< a question answered yes or no in full; on the command line one of two words (the row's words)
  kNumber,       ///< a real number in the row's range
  kWholeNumber,  ///< a whole number in the row's range
  kFile,         ///< a file that can be read
};

/// When a question of its command is asked.
enum class AskedWhen {
  kAlways,
  kBandChosen,             ///< with the uncertainty band
  kTemperaturesFromFile,   ///< the rate's temperatures from a file
  kTemperaturesInRange,    ///< the rate's temperatures from a minimum, maximum and step
  kProjectNeedsCapture,    ///< the project's run has external capture (decided when the project is read)
};

/// One question of a run: its option, what it takes, which commands ask it, the menu's texts, and where the answer
/// goes (exactly one of the four members is set, the one of the row's kind).
struct RunQuestion {
  const char *option;      ///< "--parameter-file"; "" for the menu question the command line answers by other options
  AnswerKind kind;
  unsigned commands;       ///< the commands that take it (bits CommandBit), on the command line
  unsigned menuAsks;       ///< the commands whose menu asks it
  AskedWhen when;
  bool required;           ///< no default: a command line without it is refused, the menu asks until answered
  double low;              ///< range of a number (kNumber, kWholeNumber)
  double high;
  bool lowOpen;            ///< the range excludes `low`
  bool excludesZero;       ///< 0 is no answer (a pair key)
  double defaultValue;     ///< a number's value when the answer is left blank (not `required`)
  const char *words;       ///< kWord on the command line: "yes-word|no-word"
  const char *value;       ///< the value's name in --help ("FILE", "N", "X")
  const char *help;
  // the menu
  const char *prompt;      ///< as printed
  int promptWidth;         ///< the prompt's field width (std::setw), 0 none
  bool blankLineBefore;    ///< a blank line before the prompt
  const char *endName;     ///< the prompt named when the input ends before a required answer
  const char *retry;       ///< printed after an answer out of range (kFile: after a file that cannot be read)
  // the answer
  bool RunRequest::*flag;
  int RunRequest::*whole;
  double RunRequest::*number;
  std::string RunRequest::*text;
};

/// The bit of `command` in RunQuestion::commands and menuAsks.
constexpr unsigned CommandBit(RunCommand command) { return 1u << static_cast<unsigned>(command); }

/// The questions, in the order the menu asks them.
const std::vector<RunQuestion> &RunQuestions();
/// The option's name: "--name" of "--name=VALUE".
std::string OptionName(const std::string &argument);
/// The question whose option is `option` ("--name", without a value), or null.
const RunQuestion *FindRunQuestion(const std::string &option);

/// A rule between answers: the menu asks again from `first` when it fails, the command line refuses.
struct RunRule {
  const char *first;       ///< the option of the question the menu asks again from
  const char *last;        ///< the option of the question after which the rule is checked
  bool (*holds)(const RunRequest &request);
  const char *menuMessage; ///< as the menu prints it
  const char *fix;         ///< how to fix a command line
};
const std::vector<RunRule> &RunRules();

/// Whether `question` is asked in `request` as answered so far (its command and its condition; kProjectNeedsCapture
/// counts as asked: the project decides, PrepareRun).
bool QuestionApplies(const RunQuestion &question, const RunRequest &request);
/// Whether `value` is an answer of the number question `question` (its range).
bool InRange(const RunQuestion &question, double value);

/// The command line's reader: whether `arguments` (argv without the program's name) give a command; if so, `request`
/// holds its answers (defaults where an optional answer is not given) and `rest` the arguments that are not the
/// command's (the program's other options and the project file).  Refuses (std::invalid_argument, the message in
/// three parts: what, what the command takes, how to fix) an unknown or unbuilt command, an option of another command,
/// a value that does not parse or is out of range, a missing required answer, a broken rule and a missing project.
bool ReadRunCommand(const std::vector<std::string> &arguments, RunRequest &request, std::vector<std::string> &rest);
/// The three lines of a refusal of a command line: what, what `command` takes, how to fix it.
std::string CommandRefusal(const std::string &what, RunCommand command, const std::string &fix);
/// The refusal of a command option given without a command (for parseOptions).
std::string CommandOptionWithoutCommand(const RunQuestion &question);
/// A refused command line: "ERROR: <message>" and the commands' usage line on standard output; returns 2.
int CommandError(const std::string &message);

/// The commands and their options for --help.
std::string RunCommandsHelp();
/// The commands and their options as JSON (the member "commands" of --describe --json).
std::string RunCommandsJson();

#endif
