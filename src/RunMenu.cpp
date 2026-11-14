// The text menu, the reader of a run request that asks (include/RunMenu.h): the prompts, their order and texts are
// those of the menu that drove every run before the commands existed, character for character (tests compare logs).
#include "RunMenu.h"

#include "Config.h"
#include "RunDiagnostics.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#ifdef _WIN32
#include <io.h>
#define AZURE2_ISATTY(fd) _isatty(fd)
#else
#include <unistd.h>
#define AZURE2_ISATTY(fd) isatty(fd)
#endif
#ifndef NO_READLINE
#include <readline/history.h>
#include <readline/readline.h>
#endif

namespace {

/*!
 * Ends the run when the input is exhausted at a prompt that needs an answer.
 *
 * A prompt loop re-asks until it gets a valid answer; at end of file (a piped
 * script one line short, a closed stdin, Ctrl-D) getline() fails at once on
 * every pass, so such a loop used to print its prompt forever -- one run wrote
 * 10 GB of "Number of Threads" prompts.  Prompts whose blank answer is a
 * documented default (the parameter and integral file names, the band and
 * overwrite questions, the MCMC spreads) take that default at end of file
 * instead, as before.
 */

[[noreturn]] void inputEnded(const Config &configure, const char *prompt) {
  configure.outStream << std::endl
                      << "ERROR: the input ended (end of file) at the prompt \"" << prompt
                      << "\", which needs an answer.  AZURE2 stops." << std::endl;
  configure.outStream.flush();
  ExitRun(configure, 1);
}

/*!
 * One answer from std::cin for the prompt \p prompt (named in the error);
 * ends the run at end of file (inputEnded).
 */

std::string readAnswer(const Config &configure, const char *prompt) {
  std::string answer;
  if (!getline(std::cin, answer)) inputEnded(configure, prompt);
  return answer;
}

/*!
 * One answer for a prompt whose blank answer means its default: at end of
 * file the answer is blank, i.e. the default.
 */

std::string readAnswerOrDefault() {
  std::string answer;
  if (!getline(std::cin, answer)) answer.clear();
  return answer;
}

/// Whether `answer` is blank (spaces, tabs, newlines or nothing).
bool Blank(const std::string &answer) { return answer.find_first_not_of(" \t\n") == std::string::npos; }

/*!
 * A y/n question of the text menu.  Answers (any case): y, yes, 1, t, true; n, no, 0, f, false; a blank line is "no"
 * (1 and t were "yes" before, and job scripts answer 1).  Anything else is refused.  From a pipe or a file the run stops with exit status 1: asking again would hand the NEXT scripted line to
 * this question, and that is how the trap worked -- a script that left out the band question of modes 2 and 3
 * ('3\nparam.sav\n\n7') gave the parameter-file name to it, which was read as "no", and the blank line after it to
 * the parameter-file prompt, so the run used the <levels> of the .azr without a word.  At a terminal the question is
 * asked again.  Before 2026-11-01 any answer not starting with y, 1 or t was "no".
 */
bool askYesNo(Config &configure, const std::string &prompt) {
  const bool terminal = AZURE2_ISATTY(0) != 0;
  while (true) {
    configure.outStream << prompt << " (y/n): ";
    std::string in;
    if (!getline(std::cin, in)) return false;  // end of input: the next prompt stops the run
    const size_t a = in.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return false;
    const size_t b = in.find_last_not_of(" \t\r\n");
    const std::string answer = in.substr(a, b - a + 1);
    std::string lower = answer;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    if (lower == "y" || lower == "yes" || lower == "1" || lower == "t" || lower == "true") return true;
    if (lower == "n" || lower == "no" || lower == "0" || lower == "f" || lower == "false") return false;
    configure.outStream << std::endl
                        << "ERROR: \"" << prompt << "\" takes y or n, not \"" << answer << "\"." << std::endl;
    if (terminal) continue;
    configure.outStream << "Modes 2 and 3 ask about the uncertainty band BEFORE the parameter file, so a script answers"
                        << std::endl
                        << "it first, e.g. printf '3\\nn\\nparam.sav\\n\\n7\\n' | AZURE2 --no-readline project.azr"
                        << std::endl
                        << "(--covariance-band answers it yes without a prompt).  Nothing was calculated." << std::endl;
    ExitRun(configure, EXIT_FAILURE);
  }
}

/*!
 * Reads one line of user input, either through readline or from std::cin.
 *
 * Returns false once the input stream is exhausted (end of a piped script, or
 * Ctrl-D at the terminal) and leaves \p inFile empty.  Callers must honour that:
 * readline() hands back a null pointer at EOF, and assigning it to a std::string
 * used to dereference null, so any scripted run that fell one line short of the
 * prompts crashed instead of stopping.  The plain std::cin branch had the mirror
 * problem -- getline() simply kept failing, so a prompt that insists on a
 * non-empty answer spun forever.
 */

bool getInputLine(bool useReadline, const char *prompt, std::string &inFile) {
  inFile.clear();
#ifndef NO_READLINE
  if (useReadline) {
    char *line = readline(prompt);
    if (!line) return false;
    inFile = line;
    size_t endpos = inFile.find_last_not_of(" \t");
    if (std::string::npos != endpos)
      inFile = inFile.substr(0, endpos + 1);
    else
      inFile.clear();
    if (*line) add_history(line);
    free(line);
    return true;
  }
#else
  (void)prompt;
  (void)useReadline;
#endif
  if (!getline(std::cin, inFile)) {
    inFile.clear();
    return false;
  }
  return true;
}

/// The prompt right-aligned in its field, as std::setw prints it (readline's prompt).
std::string Padded(const RunQuestion &q) {
  const std::string prompt = q.prompt;
  const std::size_t width = static_cast<std::size_t>(q.promptWidth);
  return prompt.size() < width ? std::string(width - prompt.size(), ' ') + prompt : prompt;
}

}  // namespace

std::optional<RunCommand> RunMenu::AskCommand() {
  int command = 0;
  int last = 0;
  configure_.outStream << "Please select from the following options: " << std::endl;
  for (const RunCommandRow &c : RunCommands())
    if (c.built) configure_.outStream << "\t" << (last = c.menuNumber) << ". " << c.menuText << std::endl;
  const int exitNumber = last + 1;
  configure_.outStream << "\t" << exitNumber << ". Exit" << std::endl;
  while (command < 1 || command > exitNumber) {
    configure_.outStream << "azure2: ";
    std::string inString = readAnswer(configure_, "azure2: (menu option)");
    if (inString.empty()) continue;
    std::istringstream in;
    in.str(inString);
    if (!(in >> command))
      configure_.outStream << "Please enter an integer." << std::endl;
    else if (command < 1 || command > exitNumber)
      configure_.outStream << "Invalid option.  Please try again."
                           << std::endl;
  }
  for (const RunCommandRow &c : RunCommands())
    if (c.built && c.menuNumber == command) return c.command;
  return std::nullopt;
}

void RunMenu::Ask(const RunQuestion &q, RunRequest &request) {
  std::ostream &out = configure_.outStream;
  switch (q.kind) {
    case AnswerKind::kFlag: {
      // An answer the command line gave (--covariance-band, --scale-covariance) is not asked.  With the band given
      // there, the scaling is not asked either: it is --scale-covariance's alone.
      const unsigned int preset = q.flag == &RunRequest::covarianceBand ? Config::CALCULATE_COVARIANCE_BAND
                                                                          : Config::SCALE_COVARIANCE_BY_CHI2;
      if ((configure_.paramMask & preset) ||
          (q.when == AskedWhen::kBandChosen && (configure_.paramMask & Config::CALCULATE_COVARIANCE_BAND))) {
        request.*q.flag = (configure_.paramMask & preset) != 0;
        return;
      }
      request.*q.flag = askYesNo(configure_, q.prompt);
      return;
    }
    case AnswerKind::kWord:
      while (true) {
        out << std::setw(q.promptWidth) << q.prompt;
        const std::string answer = q.required ? readAnswer(configure_, q.endName) : readAnswerOrDefault();
        if (!q.required && Blank(answer)) {  // the default: no
          request.*q.flag = false;
          return;
        }
        std::string trimmed;
        for (const char c : answer)
          if (c != ' ' && c != '\t' && c != '\n') trimmed += c;
        if (trimmed == "yes" || trimmed == "no") {
          request.*q.flag = trimmed == "yes";
          return;
        }
        out << q.retry << std::endl;
      }
    case AnswerKind::kNumber:
    case AnswerKind::kWholeNumber:
      while (true) {
        out << std::setw(q.promptWidth) << q.prompt;
        const std::string answer = q.required ? readAnswer(configure_, q.endName) : readAnswerOrDefault();
        if (!q.required && Blank(answer)) {
          request.*q.number = q.defaultValue;
          return;
        }
        std::istringstream in(answer);
        if (q.kind == AnswerKind::kWholeNumber) {
          int n = 0;
          if ((in >> n) && InRange(q, n)) {
            request.*q.whole = n;
            return;
          }
        } else {
          double x = 0.;
          if ((in >> x) && InRange(q, x)) {
            request.*q.number = x;
            return;
          }
        }
        out << q.retry << std::endl;
      }
    case AnswerKind::kFile: {
      if (q.blankLineBefore) out << std::endl;
      const std::string padded = Padded(q);
      if (!useReadline_) out << std::setw(q.promptWidth) << q.prompt;
      while (true) {
        std::string name;
        // At the end of the input an optional file takes its default (none); a required one stops the run.
        if (!getInputLine(useReadline_, padded.c_str(), name)) {
          if (q.required) inputEnded(configure_, q.endName);
          return;
        }
        if (name.empty() && !q.required) return;
        if (!name.empty()) {
          std::ifstream in(configure_.OpenPath(name).c_str());  // a name typed at a prompt is relative to the project
          if (in) {
            request.*q.text = name;
            return;
          }
        }
        if (name.empty())
          out << "Please enter a file name." << std::endl;
        else
          out << "Cannot Read From " << name << ". Please reenter file." << std::endl;
        if (!useReadline_) out << std::setw(q.promptWidth) << q.prompt;
      }
    }
  }
}

void RunMenu::Refuse(const RunRule &rule) { configure_.outStream << rule.menuMessage << std::endl; }
