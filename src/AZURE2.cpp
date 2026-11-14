/*
 * AZURE2.cpp
 * Ethan Uberseder, University of Notre Dame, 2011
 *
 * This file contains the main() function, as well as other C functions used to
 * initialize the program.  The functions defined here are (mostly) proceedural
 * C, and used to initialize the AZUREMain object from command line input.  All
 * actual calculations are managed by the AZUREMain object, and (almost) all
 * subsequent routines/classes follow object-oriented C++ programming conventions.
 */

#include "azure2/io/Document.h"
#include "azure2/experiment/Schema.h"
#include "AZUREMain.h"
#include "MeshPlanner.h"
#include "Profile.h"
#include "azure2/support/Environment.h"
#include "ProjectCommands.h"
#include "ProjectSetup.h"
#include "CaptureAyRequests.h"
#include "CommandLineValues.h"
#include "Config.h"
#include "RunDiagnostics.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/ConfigSection.h"
#include "GSLException.h"
#include "ECAmplitudeCache.h"
#include "CoulombTable.h"
#include <atomic>
#include <chrono>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <gsl/gsl_errno.h>
#ifdef _WIN32
#include <io.h>
#define AZURE2_ISATTY(fd) _isatty(fd)
#else
#include <unistd.h>
#define AZURE2_ISATTY(fd) isatty(fd)
#endif
#ifndef NO_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif
#include <string.h>
#include <cstdlib>

struct SegPairs {
  int firstPair;
  int secondPair;
};

/*!
 * This function displays the welcome banner.
 */

void welcomeMessage(const Config &configure) {
#ifndef AZURE2_VERSION_STRING
#define AZURE2_VERSION_STRING "unknown"
#endif
  std::string versionColumn = std::string("Version ") + AZURE2_VERSION_STRING;
  if (versionColumn.size() > 25) versionColumn.resize(25);
  configure.outStream << std::endl
                      << "O--------------------------O-------------------------------O" << std::endl
                      << "| #### #### #  # ###  ####  ##  | " << std::left << std::setw(25) << versionColumn << std::right << "|" << std::endl
                      << "| #  #    # #  # #  # #    #  # |                          |" << std::endl
                      << "| ####   #  #  # ###  ##     #  |                          |" << std::endl
                      << "| #  #  #   #  # # #  #     #   | Joint Institute for      |" << std::endl
                      << "| #  # ####  ##  #  # #### #### | Nuclear Astrophysics     |" << std::endl
                      << "O--------------------------O-------------------------------O" << std::endl
                      << std::endl;
}

/*!
 * This function prints a message upon successful termination of the program.
 */

void exitMessage(const Config &configure) {
  configure.outStream << std::endl
                      << "Thanks for using AZURE2." << std::endl;
}

/*!
 * The help lines of the options of the calculation, from their rows of io::ConfigKeys() (the flags the parser takes from
 * the same rows, ApplyOptionFlag): each flag says the option, the value it gives and what it implies.
 */
static void printOptionFlags() {
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    if (!k.option) continue;
    for (const bool on : {true, false}) {
      const char *flag = on ? k.flag : k.negation;
      if (!*flag) continue;
      std::string text = std::string(k.name) + "=" + (on ? "1" : "0") + " for this run" +
                         (on == k.defaultOn ? " (the default)" : "") + ", over the project's <config> key";
      for (const azure2::io::ConfigOptionImplication &i : azure2::io::ConfigOptionImplications())
        if (k.name == std::string(i.name) && on == i.value) text += std::string(" and ") + i.implied + "=" + (i.impliedValue ? "1" : "0");
      if (on) text += std::string(": ") + k.doc;
      std::cout << std::setw(25) << std::left << ("\t" + std::string(flag) + ":") << std::setw(0) << text << "." << std::endl;
    }
  }
}

/*!
 * This function prints the response to the --help command which consists of
 * available runtime options.
 */

void printHelp() {
  // Every flag parseOptions (and main) accepts, and only those
  // (tests/cli_help/check.sh compares the two).
  std::cout << "Syntax: AZURE2 <options> configfile" << std::endl
            << std::endl
            << "Options:" << std::endl
            << std::setw(25) << std::left << "\t--help:" << std::setw(0) << "Print this list and exit." << std::endl
            << std::setw(25) << std::left << "\t--no-gui:" << std::setw(0) << "Accepted and ignored, so that older scripts keep working;" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "this build has no graphical setup utility." << std::endl
            << std::setw(25) << std::left << "\t--no-readline:" << std::setw(0) << "Do not use the readline package (no effect in a build without it)." << std::endl;
  printOptionFlags();
  std::cout << std::setw(25) << std::left << "\t--threads=N:" << std::setw(0)
            << "The OpenMP threads of the run (N >= 1), over OMP_NUM_THREADS; printed by the run and by --describe."
            << std::endl
            << "The switches that can move a number (each also its environment variable; the command line wins over "
               "the options file, the file over the variable):"
            << std::endl
            << azure2::support::SwitchOptionsHelp();
  std::cout
            << std::setw(25) << std::left << "\t--convert:" << std::setw(0) << "--convert OLD NEW [--json-report R] [--answers A] [--exact-grids] [--strict] [--allow-questions]: write the old project OLD (format 1 or 2)" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "in format 3 to NEW with its side files <f>.f3<ext>, and exit (OLD kept, nothing overwritten; this AZURE2 reads format 3 only). R the JSON report" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "(notes, questions, parameter names old -> new), A the answers {\"<record>.<field>\": value}, --exact-grids the grids the old engine" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "accumulated (F10), --strict ask where the old meaning was implicit, --allow-questions '?' for unanswered questions." << std::endl
            << std::setw(25) << std::left << "\t--convert-params:" << std::setw(0) << "--convert-params OLD.sav NEW.sav --project OLD.azr: a param.sav/param.par of the old project OLD.azr with the names of its conversion, and exit." << std::endl
            << std::setw(25) << std::left << "\t--validate:" << std::setw(0) << "--validate FILE [--json]: check a format-3 project against the table and print every problem" << std::endl
            << std::setw(25) << std::left << "\t--describe:" << std::setw(0) << "--describe FILE [--json] [--record NAME]: print what the engine makes of a format-3 project, or of its record NAME (data sets, reactions and pair keys, data files, effects, parameters)" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "(file:line: record.field: message, or a JSON array); exit 0 without errors, 2 with errors, 3 for a file that is not format 3." << std::endl
            << std::setw(25) << std::left << "\t--print-predictions:" << std::setw(0) << "--print-predictions FILE: print the expression tree of every data set and extrapolation of FILE, one line each, and exit." << std::endl
            << std::setw(25) << std::left << "\t--schema:" << std::setw(0) << "Print the table of project format 3 (record kinds, fields, quantities, units, diagnostic codes) as JSON and exit." << std::endl
            << std::setw(25) << std::left << "\t--build-coulomb-table[=FILE]:" << std::setw(0) << "Build the Coulomb table (about 5 min on 4 threads, 88 MB; FILE as for --coulomb-table) and exit." << std::endl
#ifdef USE_NLOPT
            << std::setw(25) << std::left << "\t--use-nlopt:" << std::setw(0) << "Use NLopt minimizer (BOBYQA) instead of Minuit2." << std::endl
#endif
            << std::setw(25) << std::left << "\t--use-gradient:" << std::setw(0) << "Use Minuit2 (MIGRAD) with the analytic gradient (default: numerical)." << std::endl
            << std::setw(25) << std::left << "\t--refreeze-covariance[=N]:" << std::setw(0) << "Refit up to N times (default 5) with relvector covariance terms re-frozen at the fitted parameters." << std::endl
            << std::setw(25) << std::left << "\t--max-calls=N:" << std::setw(0) << "Stop a MIGRAD or NLopt fit after N function calls (default: no limit; the fit ends at convergence)." << std::endl
            << std::setw(25) << std::left << "\t--use-lm:" << std::setw(0) << "Use the Levenberg-Marquardt minimizer (analytic Jacobian; falls back to MIGRAD)." << std::endl
            << std::setw(25) << std::left << "\t--use-gsl-lm:" << std::setw(0) << "Use the GSL trust-region least-squares minimizer with geodesic acceleration (analytic Jacobian; falls back to MIGRAD)." << std::endl
            << std::setw(25) << std::left << "\t--covariance-band:" << std::setw(0) << "Write the cross-section uncertainty band from the parameter covariance, without the interactive prompt." << std::endl
            << std::setw(25) << std::left << "\t--scale-covariance:" << std::setw(0) << "With --covariance-band after a fit: inflate the covariance by the reduced chi-squared when it exceeds 1 (Birge ratio)." << std::endl
            << std::setw(25) << std::left << "\t--use-api:" << std::setw(0) << "Removed (the socket API is gone): exits with a pointer to pyazr." << std::endl;
}


/*!
 * --build-coulomb-table[=FILE]: build the universal Coulomb table with the default
 * specification (CoulombTable::BuildSpec) on all OpenMP threads, write it to FILE (default
 * CoulombTable::DefaultPath()) and report progress.  The table is not distributed, so this
 * is how a user of a ready-built AZURE2 makes one.
 */
int buildCoulombTable(const std::string &path) {
  const std::string file = path.empty() ? CoulombTable::DefaultPath() : path;
  CoulombTable::BuildSpec spec;
#ifdef _OPENMP
  spec.threads = omp_get_max_threads();
#endif
  std::atomic<int> done{0};
  std::atomic<bool> finished{false};
  spec.rowsDone = &done;
  const int rows = CoulombTable::BuildRows(spec);
  std::cout << "Building the Coulomb table " << file << " on " << spec.threads
            << " thread(s): about 5 min on 4 threads, 88 MB." << std::endl;
  const auto t0 = std::chrono::steady_clock::now();
  bool ok = false;
  std::string error;
  std::thread worker([&] {
    ok = CoulombTable::BuildFile(file, spec, error);
    finished = true;
  });
  int reported = 0;
  while (!finished) {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    const int percent = int(100.0 * done.load() / rows);
    if (percent >= reported + 10) {
      reported = percent / 10 * 10;
      std::cout << "  " << reported << " %" << std::endl;
    }
  }
  worker.join();
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (!ok) {
    std::cerr << "ERROR: Coulomb table not built: " << error << std::endl;
    return 1;
  }
  std::cout << "Coulomb table written to " << file << " in " << std::fixed << std::setprecision(0) << seconds
            << " s." << std::endl;
  return 0;
}

/*!
 * This function parses the command line options given, and sets the appropriate variables in the Config
 * structure.  It also reads and parses the configuration file if the appropriate environment
 * variable is set.  An unknown option, or a value-taking flag without a value that parses, throws
 * std::invalid_argument naming it (include/CommandLineValues.h).
 */

bool parseOptions(int argc, char *argv[], Config &configure) {
#ifndef NO_READLINE
  bool useReadline = true;
#else
  bool useReadline = false;
#endif
  configure.configfile = "";
  std::vector<std::string> options;
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (!arg.empty() && arg.substr(0, 2) == "--")
      options.push_back(arg);
    else
      configure.configfile = arg;
  }
  // The options file (--options-file=FILE, else AZURE_OPTIONS_FILE): its lines, read once with the record
  const azure2::support::Environment &environment = azure2::support::ProcessEnvironment();
  if (environment.OptionsFileUnreadable())
    configure.outStream << "The options file " << environment.Text(azure2::support::Switch::kOptionsFile)
                        << " is named but not readable." << std::endl;
  options.insert(options.end(), environment.OptionsFileLines().begin(), environment.OptionsFileLines().end());
  configure.refineQuadrature = azure2::RefineQuadratureFromEnvironment();
  for (std::vector<std::string>::iterator it = options.begin(); it < options.end(); it++) {
    if (*it == "--no-readline")  // accepted in every build (scripts pass it); without readline it changes nothing
      useReadline = false;
    else if (ApplyOptionFlag(configure, *it))  // an option of the calculation: its row of io::ConfigKeys()
      continue;
    else if (it->compare(0, 16, "--coulomb-table=") == 0) {  // the option and the table file of this machine
      GiveOption(configure, "use_coulomb_table", true, OptionSource::kCommandLine);
      configure.coulombTablePath = azure2::OptionText(*it, "--coulomb-table=");
    }
    else if (std::string value; azure2::support::ParseSwitchOption(*it, value))  // in the record (StartProcessEnvironment)
      continue;
    else if (azure2::OptionIs(*it, "--threads="))  // set by main before anything ran (ApplyThreadsOption)
      continue;
    else if (*it == "--use-api")
      configure.paramMask |= Config::USE_API;
#ifdef USE_NLOPT
    else if (*it == "--use-nlopt")
      configure.paramMask |= Config::USE_NLOPT_MINIMIZER;
#endif
    else if (*it == "--use-lm")
      configure.paramMask |= Config::USE_LM_MINIMIZER;
    else if (*it == "--use-gsl-lm")
      configure.paramMask |= Config::USE_GSL_LM_MINIMIZER;
    else if (*it == "--use-gradient")
      configure.paramMask |= Config::USE_ANALYTIC_GRADIENT;
    else if (*it == "--refreeze-covariance")
      configure.refreezeCovariance = 5;
    else if (it->compare(0, 22, "--refreeze-covariance=") == 0)
      configure.refreezeCovariance = (int)azure2::OptionWholeNumber(*it, "--refreeze-covariance=", 0, 1000);
    else if (azure2::OptionIs(*it, "--max-calls="))  // also given bare: an error, never the default
      configure.maxCalls = (unsigned int)azure2::OptionWholeNumber(*it, "--max-calls=", 0, 1000000000L);
    else if (*it == "--covariance-band")
      configure.paramMask |= Config::CALCULATE_COVARIANCE_BAND;
    else if (*it == "--scale-covariance")
      configure.paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
    else if (*it == "--no-gui")
      continue;
    else  // a mistyped flag must not run another calculation: exit 2 with the usage line (main, OptionError)
      throw std::invalid_argument("Unknown option " + *it + ".");
  }
  return useReadline;
}

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

[[noreturn]] static void inputEnded(const Config &configure, const char *prompt) {
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

static std::string readAnswer(const Config &configure, const char *prompt) {
  std::string answer;
  if (!getline(std::cin, answer)) inputEnded(configure, prompt);
  return answer;
}

/*!
 * One answer for a prompt whose blank answer means its default: at end of
 * file the answer is blank, i.e. the default.
 */

static std::string readAnswerOrDefault() {
  std::string answer;
  if (!getline(std::cin, answer)) answer.clear();
  return answer;
}

/*!
 * This function handles the command shell in AZURE2.  The function will not terminate until the user
 * enters a valid integer option, or the input ends (an error).  Upon successful entry, the integer
 * option is returned.
 */

int commandShell(const Config &configure) {
  int command = 0;

  configure.outStream << "Please select from the following options: " << std::endl
                      << "\t1. Calculate Segments From Data" << std::endl
                      << "\t2. Fit Segments From Data" << std::endl
                      << "\t3. Calculate Segments Without Data" << std::endl
                      << "\t4. Perform MINOS Error Analysis" << std::endl
                      << "\t5. Calculate Reaction Rate" << std::endl
#ifdef USE_MCMC
                      << "\t6. Perform MCMC Bayesian Sampling" << std::endl
                      << "\t7. Exit" << std::endl;
#else
                      << "\t6. Exit" << std::endl;
#endif

#ifdef USE_MCMC
  while (command < 1 || command > 7) {
#else
  while (command < 1 || command > 6) {
#endif
    configure.outStream << "azure2: ";
    std::string inString = readAnswer(configure, "azure2: (menu option)");
    if (inString.empty()) continue;
    std::istringstream in;
    in.str(inString);
    if (!(in >> command))
      configure.outStream << "Please enter an integer." << std::endl;
#ifdef USE_MCMC
    else if (command < 1 || command > 7)
#else
    else if (command < 1 || command > 6)
#endif
      configure.outStream << "Invalid option.  Please try again."
                          << std::endl;
  }
  return command;
}

/*!
 * This function takes the returned option from the commandShell function and
 * sets the appropriate flags in the Config structure.
 */

// Read a yes/no answer from stdin; empty/unrecognized input counts as "no".
/*!
 * A y/n question of the text menu.  Answers (any case): y, yes, 1, t, true; n, no, 0, f, false; a blank line is "no"
 * (1 and t were "yes" before, and job scripts answer 1).  Anything else is refused.  From a pipe or a file the run stops with exit status 1: asking again would hand the NEXT scripted line to
 * this question, and that is how the trap worked -- a script that left out the band question of modes 2 and 3
 * ('3\nparam.sav\n\n7') gave the parameter-file name to it, which was read as "no", and the blank line after it to
 * the parameter-file prompt, so the run used the <levels> of the .azr without a word.  At a terminal the question is
 * asked again.  Before 2026-11-01 any answer not starting with y, 1 or t was "no".
 */
static bool askYesNo(Config &configure, const std::string &prompt) {
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

// Offer the cross-section uncertainty band (and, for fit modes, the reduced-chi^2
// covariance scaling).  If --covariance-band was passed the choice is already
// made non-interactively, so no prompt appears -- this keeps scripted/headless
// runs deterministic while interactive runs are asked.
static void promptBandOptions(Config &configure, bool fitMode) {
  if (configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) return;  // set via flag
  if (!askYesNo(configure, "Calculate cross-section uncertainty band?")) return;
  configure.paramMask |= Config::CALCULATE_COVARIANCE_BAND;
  if (fitMode && !(configure.paramMask & Config::SCALE_COVARIANCE_BY_CHI2)) {
    if (askYesNo(configure, "Scale covariance to reduced chi-squared = 1?"))
      configure.paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  }
}

void processCommand(int command, Config &configure) {
  if (command == 2) {
    configure.paramMask |= Config::PERFORM_FIT;
    promptBandOptions(configure, true);
  } else if (command == 3) {
    configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
    promptBandOptions(configure, false);  // extrapolation reuses a saved covariance
  } else if (command == 4) {
    bool goodAnswer = false;
    while (!goodAnswer) {
      configure.outStream << std::setw(30) << "Allowed Chi-Squared Variance: ";
      std::string inString = readAnswer(configure, "Allowed Chi-Squared Variance");
      std::istringstream stm;
      stm.str(inString);
      if (!(stm >> configure.chiVariance) || configure.chiVariance < 0.)
        configure.outStream << "Please enter a positive number." << std::endl;
      else
        goodAnswer = true;
    }
    configure.paramMask |= Config::PERFORM_FIT;
    configure.paramMask |= Config::PERFORM_ERROR_ANALYSIS;
    promptBandOptions(configure, true);
  } else if (command == 5) {
    configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
    configure.paramMask |= Config::CALCULATE_REACTION_RATE;
#ifdef USE_MCMC
  } else if (command == 6) {
    configure.paramMask |= Config::PERFORM_MCMC;
  } else if (command == 7) {
    exitMessage(configure);
    ExitRun(configure, 0);
  }
#else
  } else if (command == 6) {
    exitMessage(configure);
    ExitRun(configure, 0);
  }
#endif
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

static bool getInputLine(bool useReadline, const char *prompt, std::string &inFile) {
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

/*!
 * This function prompts for a parameter file and sets the corresponding configure
 * flags and variables based on the user response.
 */

void getParameterFile(bool useReadline, Config &configure) {
  bool validInfile = false;
  configure.outStream << std::endl;
  if (!useReadline) configure.outStream << "External Parameter File (leave blank for new file): ";
  while (!validInfile) {
    std::string inFile;
    // At EOF fall back to the "leave blank" behaviour: build a new parameter file.
    if (!getInputLine(useReadline, "External Parameter File (leave blank for new file): ", inFile)) return;
    if (!inFile.empty()) {
      std::ifstream in;
      in.open(configure.OpenPath(inFile).c_str());  // a name typed at a prompt is relative to the project, as its paths
      if (in) {
        validInfile = true;
        configure.paramMask |= Config::USE_PREVIOUS_PARAMETERS;
        configure.paramfile = inFile;
        in.close();
      }
      in.clear();
    } else
      validInfile = true;
    if (!validInfile) {
      configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
      if (!useReadline) configure.outStream << "External Parameter File (leave blank for new file): ";
    }
  }
}

/*!
 * The section <name> of the project as a stream over its text (Config::SectionStream: io::Document, whole-line
 * tags; docs/upgrade/project.md), or a stream in the fail state when there is no such section.
 */
static std::istringstream SectionText(const Config &configure, const std::string &name) {
  std::istringstream in;
  if (std::unique_ptr<std::istream> section = configure.SectionStream(name)) {
    std::ostringstream text;
    text << section->rdbuf();
    in.str(text.str());
  } else {
    in.setstate(std::ios::failbit);
  }
  return in;
}

/*!
 * This function reads the segment file, and stores the active entrance and exit pair keys
 * for cross reference with the External capture file.  Only if an active external capture segment
 * is required is the user prompted for an external integrals file.
 */

bool readSegmentFile(const Config &configure, std::vector<SegPairs> &segPairs) {
  // the active data sets (or extrapolations, for a run without data) from their specs (ProjectSetup.h)
  for (const std::pair<int, int> &p : ProjectReactionPairs(configure, (configure.paramMask & Config::CALCULATE_WITH_DATA) != 0))
    segPairs.push_back({p.first, p.second});
  return true;
}

/*!
 * If reaction rate is desired, the user may specify a file containing temperatures for the calculation.
 * This function prompts for that file name, checks for access, and stores path.
 */

void getTemperatureFile(bool useReadline, Config &configure) {
  bool validInfile = false;
  if (!useReadline) configure.outStream << std::setw(38) << "Temperature File Name: ";
  while (!validInfile) {
    std::string inFile;
    // No file name can be obtained at EOF, and the rate needs one.
    if (!getInputLine(useReadline, "               Temperature File Name: ", inFile))
      inputEnded(configure, "Temperature File Name");
    if (!inFile.empty()) {
      std::ifstream in;
      in.open(configure.OpenPath(inFile).c_str());  // a name typed at a prompt is relative to the project, as its paths
      if (in) {
        validInfile = true;
        configure.rateParams.temperatureFile = inFile;
        in.close();
      }
      in.clear();
    }
    if (!validInfile) {
      if (inFile.empty())
        configure.outStream << "Please enter a file name." << std::endl;
      else
        configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
      if (!useReadline) configure.outStream << "               Temperature File Name: ";
    }
  }
}

/*!
 * This function prompts the user for MCMC sampling parameters.
 */

#ifdef USE_MCMC
// Upper bound on the initial spread of a level-energy parameter, in keV.
// Scattering level energies by a percentage of their value puts walkers on
// unrelated resonance structures and the ensemble never contracts.
static const double kMaxEnergySpreadKeV = 1.0;

struct MCMCParams {
  int nwalkers;
  int nsteps;
  double chainSpread;
  double energySpreadKeV;
  int nthreads;
  bool useRWA;
  bool overwriteSamples;
};

void getMCMCParams(Config &configure, MCMCParams &mcmcParams) {
  mcmcParams.nwalkers = 0;
  mcmcParams.nsteps = 0;
  mcmcParams.chainSpread = -1.0;      // Will be prompted
  mcmcParams.energySpreadKeV = -1.0;  // Will be prompted
  mcmcParams.nthreads = 1;
  mcmcParams.useRWA = false;

  // Get number of walkers
  while (mcmcParams.nwalkers < 2 || mcmcParams.nwalkers > 1000) {
    configure.outStream << std::setw(38) << "Number of Walkers (2-1000): ";
    std::string inString = readAnswer(configure, "Number of Walkers");
    std::istringstream stm;
    stm.str(inString);
    if (!(stm >> mcmcParams.nwalkers) || mcmcParams.nwalkers < 2 || mcmcParams.nwalkers > 1000)
      configure.outStream << "Please enter an integer between 2 and 1000." << std::endl;
  }

  // Get number of steps
  while (mcmcParams.nsteps < 100 || mcmcParams.nsteps > 10000000) {
    configure.outStream << std::setw(38) << "Number of Steps (100-10000000): ";
    std::string inString = readAnswer(configure, "Number of Steps");
    std::istringstream stm;
    stm.str(inString);
    if (!(stm >> mcmcParams.nsteps) || mcmcParams.nsteps < 100 || mcmcParams.nsteps > 10000000)
      configure.outStream << "Please enter an integer between 100 and 10000000." << std::endl;
  }

  // Get initial parameter spread percentage
  while (mcmcParams.chainSpread <= 0.0 || mcmcParams.chainSpread > 50.0) {
    configure.outStream << std::setw(38) << "Initial Parameter Spread % (0.001-50) [default=5]: ";
    std::string inString = readAnswerOrDefault();

    // Allow empty input for default
    if (inString.empty() || inString.find_first_not_of(" \t\n") == std::string::npos) {
      mcmcParams.chainSpread = 5.0;
      break;
    }

    std::istringstream stm;
    stm.str(inString);
    if (!(stm >> mcmcParams.chainSpread) || mcmcParams.chainSpread <= 0.0 || mcmcParams.chainSpread > 50.0)
      configure.outStream << "Please enter a number between 0.001 and 50, or press Enter for default (5%)." << std::endl;
  }

  // Get the initial spread of the level energies, which is absolute rather than
  // a percentage of the level energy.
  while (mcmcParams.energySpreadKeV <= 0.0 || mcmcParams.energySpreadKeV > kMaxEnergySpreadKeV) {
    configure.outStream << std::setw(38) << "Level Energy Spread keV (0.001-1) [default=1]: ";
    std::string inString = readAnswerOrDefault();

    // Allow empty input for default
    if (inString.empty() || inString.find_first_not_of(" \t\n") == std::string::npos) {
      mcmcParams.energySpreadKeV = kMaxEnergySpreadKeV;
      break;
    }

    std::istringstream stm;
    stm.str(inString);
    if (!(stm >> mcmcParams.energySpreadKeV) || mcmcParams.energySpreadKeV <= 0.0 ||
        mcmcParams.energySpreadKeV > kMaxEnergySpreadKeV)
      configure.outStream << "Please enter a number between 0.001 and " << kMaxEnergySpreadKeV
                          << ", or press Enter for default (1 keV)." << std::endl;
  }

  // Get number of threads
  mcmcParams.nthreads = 0;
  while (mcmcParams.nthreads < 1) {
    configure.outStream << std::setw(38) << "Number of Threads (1 or more): ";
    std::string inString = readAnswer(configure, "Number of Threads");
    std::istringstream stm;
    stm.str(inString);
    if (!(stm >> mcmcParams.nthreads) || mcmcParams.nthreads < 1)
      configure.outStream << "Please enter a positive integer." << std::endl;
  }

  // Ask if using RWA parameters
  bool goodAnswer = false;
  std::string rwaAnswer = "";
  while (!goodAnswer) {
    configure.outStream << std::setw(38) << "Use Reduced Width Amplitudes (yes/no): ";
    rwaAnswer = readAnswer(configure, "Use Reduced Width Amplitudes");
    std::string trimmedAnswer;
    for (int i = 0; i < rwaAnswer.length(); i++)
      if (rwaAnswer[i] != ' ' && rwaAnswer[i] != '\t' && rwaAnswer[i] != '\n')
        trimmedAnswer += rwaAnswer[i];
    if (trimmedAnswer != "yes" && trimmedAnswer != "no")
      configure.outStream << "Please type 'yes' or 'no'." << std::endl;
    else {
      goodAnswer = true;
      rwaAnswer = trimmedAnswer;
    }
  }
  mcmcParams.useRWA = (rwaAnswer == "yes");

  // Ask if overwriting existing samples.mcmc file
  goodAnswer = false;
  std::string overwriteAnswer = "";
  while (!goodAnswer) {
    configure.outStream << std::setw(38) << "Overwrite existing samples.mcmc? (yes/no) [default=no]: ";
    overwriteAnswer = readAnswerOrDefault();

    // Allow empty input for default (no)
    if (overwriteAnswer.empty() || overwriteAnswer.find_first_not_of(" \t\n") == std::string::npos) {
      mcmcParams.overwriteSamples = false;
      break;
    }

    std::string trimmedAnswer;
    for (int i = 0; i < overwriteAnswer.length(); i++)
      if (overwriteAnswer[i] != ' ' && overwriteAnswer[i] != '\t' && overwriteAnswer[i] != '\n')
        trimmedAnswer += overwriteAnswer[i];
    if (trimmedAnswer != "yes" && trimmedAnswer != "no")
      configure.outStream << "Please type 'yes' or 'no', or press Enter for default (no)." << std::endl;
    else {
      goodAnswer = true;
      mcmcParams.overwriteSamples = (trimmedAnswer == "yes");
    }
  }
}
#endif

/*!
 * This function  prompts the user for the required parameters if reaction rate is
 * to be calculated.
 */

void getRateParams(Config &configure, std::vector<SegPairs> &segPairs, bool useReadline) {
  configure.rateParams.entrancePair = 0;
  configure.rateParams.exitPair = 0;
  configure.rateParams.minTemp = -1.;
  configure.rateParams.maxTemp = -1.;
  configure.rateParams.tempStep = -1.;
  while (configure.rateParams.entrancePair == configure.rateParams.exitPair) {
    while (!configure.rateParams.entrancePair) {
      configure.outStream << std::setw(38) << "Reaction Rate Entrance Pair: ";
      std::string inString = readAnswer(configure, "Reaction Rate Entrance Pair");
      std::istringstream stm;
      stm.str(inString);
      if (!(stm >> configure.rateParams.entrancePair) || configure.rateParams.entrancePair == 0)
        configure.outStream << "Please enter an integer greater than zero." << std::endl;
    }
    while (!configure.rateParams.exitPair) {
      configure.outStream << std::setw(38) << "Reaction Rate Exit Pair: ";
      std::string inString = readAnswer(configure, "Reaction Rate Exit Pair");
      std::istringstream stm;
      stm.str(inString);
      // -1: the total radiative capture, summed over every capture pair of the entrance nucleus (the exit
      // key a total-capture segment line carries).
      if (!(stm >> configure.rateParams.exitPair) || configure.rateParams.exitPair == 0 ||
          configure.rateParams.exitPair < -1) {
        configure.outStream << "Please enter an integer greater than zero, or -1 for the total radiative capture."
                            << std::endl;
        configure.rateParams.exitPair = 0;
      }
    }
    if (configure.rateParams.entrancePair == configure.rateParams.exitPair) {
      configure.outStream << "Cannot calculate rate for elastic scattering." << std::endl;
      configure.rateParams.entrancePair = 0;
      configure.rateParams.exitPair = 0;
    }
  }
  SegPairs tempSet = {configure.rateParams.entrancePair, configure.rateParams.exitPair};
  segPairs.push_back(tempSet);
  bool goodAnswer = false;
  std::string fileAnswer = "";
  while (!goodAnswer) {
    configure.outStream << std::setw(38) << "Use temperatures from file (yes/no): ";
    fileAnswer = readAnswer(configure, "Use temperatures from file");
    std::string trimmedAnswer;
    for (int i = 0; i < fileAnswer.length(); i++)
      if (fileAnswer[i] != ' ' && fileAnswer[i] != '\t' && fileAnswer[i] != '\n')
        trimmedAnswer += fileAnswer[i];
    if (trimmedAnswer != "yes" && trimmedAnswer != "no")
      configure.outStream << "Please type 'yes' or 'no'." << std::endl;
    else {
      goodAnswer = true;
      fileAnswer = trimmedAnswer;
    }
  }
  configure.rateParams.useFile = (fileAnswer == "yes") ? (true) : (false);
  if (configure.rateParams.useFile)
    getTemperatureFile(useReadline, configure);
  else {
    while (configure.rateParams.minTemp < 0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Min Temp [GK]: ";
      std::string inString = readAnswer(configure, "Reaction Rate Min Temp");
      std::istringstream stm;
      stm.str(inString);
      if (!(stm >> configure.rateParams.minTemp) || configure.rateParams.minTemp < 0.)
        configure.outStream << "Please enter a positive number." << std::endl;
    }
    while (configure.rateParams.maxTemp < 0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Max Temp [GK]: ";
      std::string inString = readAnswer(configure, "Reaction Rate Max Temp");
      std::istringstream stm;
      stm.str(inString);
      if (!(stm >> configure.rateParams.maxTemp) || configure.rateParams.maxTemp < 0.)
        configure.outStream << "Please enter a positive number." << std::endl;
    }
    while (configure.rateParams.tempStep < 0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Temp Step [GK]: ";
      std::string inString = readAnswer(configure, "Reaction Rate Temp Step");
      std::istringstream stm;
      stm.str(inString);
      if (!(stm >> configure.rateParams.tempStep) || configure.rateParams.tempStep < 0.)
        configure.outStream << "Please enter a positive number." << std::endl;
    }
  }
}

/*!
 * This function checks the external capture file against a vector of segment key
 * pairs. Only if the calculation includes external capture segments is the user
 * prompted for an integrals file.  The appropriate configure flag is set here.
 */

bool checkExternalCapture(Config &configure, const std::vector<SegPairs> &segPairs) {
  std::istringstream in = SectionText(configure, "levels");
  if (!in) {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;
  }
  std::string line = "";
  while (line != "<levels>" && !in.eof()) getline(in, line);
  if (line != "<levels>") {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;
  }
  line = "";
  while (line != "</levels>" && !in.eof() &&
         !(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (line != "</levels>" && !in.eof()) {
      std::istringstream stm;
      stm.str(line);
      azure2::io::LevelRow tempNucLine(stm);
      if ((stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit))) {
        configure.outStream << "Problem reading nuclear information. Check configuration file.  The line:\n  "
                            << line << std::endl;
        return false;
      }
      if (tempNucLine.ecMultMask() != 0) {
        for (int i = 0; i < segPairs.size(); i++) {
          if (tempNucLine.ir() == segPairs[i].secondPair || segPairs[i].secondPair == -1) {
            configure.paramMask |= Config::USE_EXTERNAL_CAPTURE;
            break;
          }
        }
      }
    }
  }
  in.clear();
  if ((configure.paramMask & Config::USE_EXTERNAL_CAPTURE) && (configure.paramMask & Config::USE_RMC_FORMALISM)) {
    configure.outStream << "WARNING: External capture is not compatible with Reich-Moore.  Ignoring external capture." << std::endl;
    configure.paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  }
  return true;
}

/*!
 * This function promps the user for an external capture integrals file, checks for it's validity,
 * and stores the path.
 */

void getExternalCaptureFile(bool useReadline, Config &configure) {
  if ((configure.paramMask & Config::USE_EXTERNAL_CAPTURE) &&
      !(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
    configure.outStream << std::endl;
    if (!useReadline) configure.outStream << "External Capture Amplitude File (leave blank for new file): ";
    bool validInfile = false;
    while (!validInfile) {
      std::string inFile;
      // At EOF fall back to the "leave blank" behaviour: build a new integrals file.
      if (!getInputLine(useReadline, "External Capture Amplitude File (leave blank for new file): ", inFile)) return;
      if (!inFile.empty()) {
        std::ifstream in;
        in.open(configure.OpenPath(inFile).c_str());
        if (in) {
          validInfile = true;
          configure.paramMask |= Config::USE_PREVIOUS_INTEGRALS;
          configure.integralsfile = inFile;
          in.close();
        }
        in.clear();
      } else
        validInfile = true;
      if (!validInfile) {
        configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
        if (!useReadline) configure.outStream << "External Capture Amplitude File (leave blank for new file): ";
      }
    }
  }
}

/*!
 * This function prints a brief start message describing the type of calculation that will be performed.
 */

void startMessage(const Config &configure) {
  if (configure.paramMask & Config::PERFORM_ERROR_ANALYSIS)
    configure.outStream << "Calling AZURE2 for MINOS Error Analysis..." << std::endl;
#ifdef USE_MCMC
  else if (configure.paramMask & Config::PERFORM_MCMC)
    configure.outStream << "Calling AZURE2 for MCMC Bayesian sampling..." << std::endl;
#endif
  else if (configure.paramMask & Config::CALCULATE_REACTION_RATE)
    configure.outStream << "Calling AZURE2 for reaction rate calculation..." << std::endl;
  else if (!(configure.paramMask & Config::CALCULATE_WITH_DATA))
    configure.outStream << "Calling AZURE2 for calculation of segments without data..." << std::endl;
  else if (configure.paramMask & Config::PERFORM_FIT)
    configure.outStream << "Calling AZURE2 for fitting of segments from data..." << std::endl;
  else
    configure.outStream << "Calling AZURE2 for calculation of segments from data..." << std::endl;
}

/*!
 * This function runs the MCMC sampling procedure.
 */

#ifdef USE_MCMC
#include "CNuc.h"
#include "WorkingCopy.h"
#include "AZUREParams.h"
#include "AZURECalcMCMC.h"
#include <iomanip>
#include <ctime>

// The CLI's progress line of an MCMC run: printed once per step, with an estimate of the time left from the steps
// done since `start`.
static void PrintMcmcProgress(time_t start, int current_step, int total_steps, double logprob, double loglikelihood,
                              double logprior) {
  // Calculate time estimates
  time_t current_time = time(NULL);
  double elapsed_seconds = difftime(current_time, start);
  double steps_per_second = (elapsed_seconds > 0) ? current_step / elapsed_seconds : 0;
  double remaining_steps = total_steps - current_step;
  double estimated_remaining_seconds = (steps_per_second > 0) ? remaining_steps / steps_per_second : 0;

  // Format time
  int hours = (int)(estimated_remaining_seconds / 3600);
  int minutes = (int)((estimated_remaining_seconds - hours * 3600) / 60);
  int seconds = (int)(estimated_remaining_seconds - hours * 3600 - minutes * 60);

  std::cout << "\rStep " << std::setw(6) << current_step << "/" << total_steps
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * current_step / total_steps) << "%)"
            << " | logP=" << std::setprecision(2) << logprob
            << " | logL=" << loglikelihood
            << " | logPrior=" << logprior
            << " | ETA: " << std::setfill('0') << std::setw(2) << hours << ":"
            << std::setw(2) << minutes << ":" << std::setw(2) << seconds
            << std::setfill(' ') << std::flush;
}

// The name the MCMC table prints for a parameter (the GUI's MCMC tab writes the same), from its record.
static std::string PrettyParameterName(const ParameterRecord &r, const std::string &name) {
  const std::string seg = "Segment " + std::to_string(r.segmentKey);
  switch (r.kind) {
    case ParamKind::LevelEnergy: return "Level " + std::to_string(r.levelNumber) + " Energy (MeV)";
    case ParamKind::Gamma:
      return "Level " + std::to_string(r.levelNumber) + " Channel " + std::to_string(r.channel) + " Width (eV)";
    case ParamKind::Norm: return seg + " Normalization";
    case ParamKind::Thickness: return seg + " Sample Thickness (atoms/b)";
    case ParamKind::EnergyShift: return seg + " Energy Shift (keV)";
    case ParamKind::EnergyShiftSqrt: return seg + " sqrt(E) Energy-Shift Coefficient (MeV^1/2)";
    case ParamKind::EnergyShiftFn: return seg + " Energy-Shift Function a" + std::to_string(r.fnTerm);
    case ParamKind::Screening: return "Screening U_e (eV), effect " + std::to_string(r.effect);
    case ParamKind::Density: return "Target density (atoms/cm^2), effect " + std::to_string(r.effect);
    case ParamKind::Doppler: return "Doppler temperature (K), effect " + std::to_string(r.effect);
    case ParamKind::Abundance:
      return "Abundance " + (r.sampleName.empty() ? std::string() : "sample " + r.sampleName + " ") +
             "Z=" + std::to_string(r.nucleusZ) + " A=" + std::to_string(r.nucleusA);
    case ParamKind::ThmCoherent: case ParamKind::Apparatus: case ParamKind::RadiativeDecay: break;  // their names say what they are
  }
  return name;
}

/*!
 * The parameter a row of the <mcmc> <parameters> block names, or -1.  The row gives the parameter's own name (whole,
 * not a prefix), or one of the GUI's short forms E<n> (energy of level n), G<n><c> (width of level n, one digit, in
 * channel c), N<k> (normalization of segment k), S<k> (constant energy shift of segment k), which are resolved through
 * the parameters' records.
 */
static int FindMcmcPriorParameter(const std::string &token, const AZUREParams &params) {
  const int exact = params.IndexOfName(token);
  if (exact >= 0) return exact;
  if (token.size() < 2 || !std::isdigit((unsigned char)token[1])) return -1;
  ParamKind kind = ParamKind::LevelEnergy;
  int first = 0, second = -1;
  try {
    if (token[0] == 'E' || token[0] == 'N' || token[0] == 'S') {
      kind = token[0] == 'E' ? ParamKind::LevelEnergy : token[0] == 'N' ? ParamKind::Norm : ParamKind::EnergyShift;
      first = std::stoi(token.substr(1));
    } else if (token[0] == 'G' && token.size() > 2) {
      kind = ParamKind::Gamma;
      first = std::stoi(token.substr(1, 1));
      second = std::stoi(token.substr(2));
      if (second <= 0) return -1;
    } else {
      return -1;
    }
  } catch (...) {
    return -1;
  }
  if (first <= 0) return -1;
  for (int i = 0; i < (int)params.Records().size(); i++) {
    const ParameterRecord &r = params.Records()[i];
    if (r.kind != kind) continue;
    if ((kind == ParamKind::LevelEnergy && r.levelNumber == first) ||
        (kind == ParamKind::Gamma && r.levelNumber == first && r.channel == second) ||
        ((kind == ParamKind::Norm || kind == ParamKind::EnergyShift) && r.segmentKey == first))
      return i;
  }
  return -1;
}

int runMCMC(Config &configure, const MCMCParams &mcmcParams) {
  try {
    CNuc *compound = new CNuc();
    WorkingCopy *data = new WorkingCopy();

    configure.outStream << "Filling Compound Nucleus..." << std::endl;
    if (compound->Fill(configure) == -1) {
      configure.outStream << "Could not fill compound nucleus from file." << std::endl;
      delete compound;
      delete data;
      return -1;
    } else if (compound->NumPairs() == 0 || compound->NumJGroups() == 0) {
      configure.outStream << "No nuclear data exists. Calculation not possible." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Fill the data object from segments and data file
    configure.outStream << "Filling Data Structures..." << std::endl;
    if (configure.paramMask & Config::CALCULATE_WITH_DATA) {
      if (const int filled = data->Fill(configure, compound); filled != 0) {
        if (filled != WorkingCopy::kFillRefused) configure.outStream << "Could not fill data object from file." << std::endl;
        delete compound;
        delete data;
        return -1;
      } else if (data->Data().NumSegments() == 0 && compound->GetModelData().NumResiduals() == 0) {  // model data alone run
        configure.outStream << "There is no data provided." << std::endl;
        delete compound;
        delete data;
        return -1;
      }
      RequestDataEntrances(*compound, configure);
    } else {
      if (const int made = data->MakePoints(configure, compound); made != 0) {
        if (made != WorkingCopy::kFillRefused) configure.outStream << "Could not fill data object from file." << std::endl;
        delete compound;
        delete data;
        return -1;
      } else if (data->Data().NumSegments() == 0) {
        configure.outStream << "Extrapolation segments produce no data." << std::endl;
        delete compound;
        delete data;
        return -1;
      }
      RequestExtrapolationEntrances(*compound, configure);
    }

    // Initialize compound nucleus, with the capture-A_y terms the data sets need
    RequestCaptureAyTerms(*compound, *data);
    try {
      compound->Initialize(configure);
    } catch (GSLException e) {
      configure.outStream << e.what() << std::endl;
      configure.outStream << "Calculation was aborted." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Initialize EC amplitude cache if using external capture
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      InitializeECAmplitudeCache(configure);
      configure.outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
    }

    const int initialized = data->Initialize(compound, configure);  // AZURE_PROFILE: the set-up's stages on their own
    AZUREProfile::Report(configure.profile, std::cerr, "set-up");
    if (initialized == -1) {
      configure.outStream << "Calculation was aborted." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Create parameter structure
    AZUREParams params;
    compound->FillMnParams(params, &configure);
    data->FillMnParams(params);

    // Read existing parameters if provided
    if (configure.paramMask & Config::USE_PREVIOUS_PARAMETERS) {
      configure.outStream << "Reading User Parameter File..." << std::endl;
      if (!params.ReadUserParameters(configure) || !params.ReconcileBasis(compound, configure)) return -1;
    } else {
      configure.outStream << "Creating New param.par File..." << std::endl;
      params.WriteUserParameters(configure, false);
    }

    // Setup compound nucleus from parameters
    compound->FillCompoundFromParams(params.GetMinuitParams().Params());
    if (configure.paramMask & Config::USE_BRUNE_FORMALISM) {
      compound->CalcShiftFunctions(configure);
    }

    // Extract initial parameter values (only non-fixed parameters)
    std::vector<double> initialParams;
    if (mcmcParams.useRWA) {
      // Use RWA parameters directly
      for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
        if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
          initialParams.push_back(params.GetMinuitParams().Parameter(i).Value());
        }
      }
    } else {
      // Use physical parameters
      compound->TransformOut(configure);
      vector_r physicalParams = compound->GetTransformParams(configure);

      // Extend with data parameters if needed
      for (size_t i = physicalParams.size(); i < params.GetMinuitParams().Params().size(); i++) {
        physicalParams.push_back(params.GetMinuitParams().Parameter(i).Value());
      }

      // Extract non-fixed parameters
      for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
        if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
          initialParams.push_back(physicalParams[i]);
        }
      }
    }

    if (initialParams.empty()) {
      configure.outStream << "Error: No free parameters found for MCMC sampling." << std::endl;
      configure.outStream << "Please ensure at least one parameter is not fixed." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Read priors from AZURE2 file if they exist
    std::map<int, std::tuple<double, double, bool>> priorMap;  // parameter index -> (priorMean, priorStd, usePrior)
    std::istringstream priorFile = SectionText(configure, "mcmc");
    if (priorFile.good()) {
      std::string line;
      bool inMCMCSection = false;
      bool inParametersSection = false;

      while (std::getline(priorFile, line)) {
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start != std::string::npos) {
          line = line.substr(start);
          size_t end = line.find_last_not_of(" \t\r\n");
          if (end != std::string::npos) line = line.substr(0, end + 1);
        }

        if (line == "<mcmc>") {
          inMCMCSection = true;
          continue;
        } else if (line == "</mcmc>") {
          inMCMCSection = false;
          break;
        }

        if (inMCMCSection) {
          if (line == "<parameters>") {
            inParametersSection = true;
            continue;
          } else if (line == "</parameters>") {
            inParametersSection = false;
            continue;
          }

          if (inParametersSection && !line.empty()) {
            // Parse: "Parameter (raw_param_name) currentValue priorMean priorStd useGaussianPrior"
            std::istringstream iss(line);
            std::vector<std::string> parts;
            std::string part;
            while (iss >> part) parts.push_back(part);

            if (parts.size() >= 5) {
              // Last 4 parts are: currentValue, priorMean, priorStd, useGaussianPrior
              // Everything before that is the parameter name (possibly "Parameter (raw_name)")
              std::string fileParamName;
              for (size_t i = 0; i < parts.size() - 4; i++) {
                if (i > 0) fileParamName += " ";
                fileParamName += parts[i];
              }

              // Extract raw parameter name from "Parameter (raw_name)" format
              std::string rawParamName;
              size_t openParen = fileParamName.find('(');
              size_t closeParen = fileParamName.find(')');
              if (openParen != std::string::npos && closeParen != std::string::npos && closeParen > openParen) {
                // Extract the raw name from inside parentheses
                rawParamName = fileParamName.substr(openParen + 1, closeParen - openParen - 1);
              } else {
                // No parentheses - use the whole name
                rawParamName = fileParamName;
              }

              const int paramOfRow = FindMcmcPriorParameter(rawParamName, params);

              double priorMean = std::atof(parts[parts.size() - 3].c_str());
              double priorStd = std::atof(parts[parts.size() - 2].c_str());
              bool usePrior = (std::atoi(parts[parts.size() - 1].c_str()) == 1);

              if (paramOfRow < 0) {
                // A prior that is silently not applied is a wrong chain that looks right.
                configure.outStream << "ERROR: <mcmc> <parameters> row \"" << line << "\": \"" << rawParamName
                                    << "\" names no parameter of this model.  A row names a parameter by its whole "
                                    << "name as in param.sav (energy_1, width_1_2, segment_1_norm, ...), optionally "
                                    << "as \"label (name)\", or by the short forms E<n>, G<n><c>, N<k>, S<k>, "
                                    << "followed by: value, prior mean, prior sigma, use (0 or 1)." << std::endl;
                delete compound;
                delete data;
                return -1;
              }
              priorMap[paramOfRow] = std::make_tuple(priorMean, priorStd, usePrior);
            } else {
              configure.outStream << "ERROR: <mcmc> <parameters> row \"" << line << "\": expected a parameter name "
                                  << "followed by: value, prior mean, prior sigma, use (0 or 1)." << std::endl;
              delete compound;
              delete data;
              return -1;
            }
          }
        }
      }
    }

    // Build prior vectors for MCMC calculator (matching GUI approach)
    std::vector<double> priorMeans;
    std::vector<double> priorStds;
    std::vector<bool> usePriors;

    // Print parameter information
    configure.outStream << std::endl;
    configure.outStream << "Free Parameters for MCMC:" << std::endl;
    configure.outStream << std::setw(40) << std::left << "Parameter Name"
                        << std::setw(15) << "Current Value"
                        << std::setw(15) << "Prior Mean"
                        << std::setw(15) << "Prior Std"
                        << std::setw(20) << "Use Prior" << std::endl;
    configure.outStream << std::string(105, '-') << std::endl;

    int paramIndex = 0;
    int priorsFound = 0;
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
        std::string paramName = PrettyParameterName(params.Record(i), params.GetMinuitParams().GetName(i));
        double currentValue = initialParams[paramIndex];

        // Check if prior exists for this parameter
        std::string priorMeanStr = "N/A";
        std::string priorStdStr = "N/A";
        std::string usePriorStr = "No (uniform)";

        double priorMean = 0.0;
        double priorStd = 0.0;
        bool usePrior = false;

        if (priorMap.find(i) != priorMap.end()) {
          priorMean = std::get<0>(priorMap[i]);
          priorStd = std::get<1>(priorMap[i]);
          usePrior = std::get<2>(priorMap[i]);

          // Only mark as found if usePrior is true AND std is positive
          if (usePrior && priorStd > 0.0) {
            std::ostringstream meanStream, stdStream;
            meanStream << std::scientific << std::setprecision(4) << priorMean;
            stdStream << std::scientific << std::setprecision(4) << priorStd;
            priorMeanStr = meanStream.str();
            priorStdStr = stdStream.str();
            usePriorStr = "Yes (Gaussian)";
            priorsFound++;
          } else {
            // Found in map but usePrior is false or std is invalid - reset to defaults
            priorMean = 0.0;
            priorStd = 0.0;
            usePrior = false;
          }
        }

        // Add to prior vectors - these must be added for ALL parameters to maintain index alignment
        priorMeans.push_back(priorMean);
        priorStds.push_back(priorStd);
        usePriors.push_back(usePrior);

        configure.outStream << std::setw(40) << std::left << paramName
                            << std::setw(15) << std::scientific << std::setprecision(6) << currentValue
                            << std::setw(15) << priorMeanStr
                            << std::setw(15) << priorStdStr
                            << std::setw(20) << usePriorStr << std::endl;
        paramIndex++;
      }
    }
    configure.outStream << std::endl;
    if (priorsFound > 0) {
      configure.outStream << "Found " << priorsFound << " Gaussian prior(s) in AZURE2 file." << std::endl;
      configure.outStream << "Priors will be applied during MCMC sampling." << std::endl;
    } else {
      configure.outStream << "No priors found in AZURE2 file. Using uniform priors for all parameters." << std::endl;
      configure.outStream << "To use Gaussian priors, use the GUI to configure them and save the .azr file." << std::endl;
    }
    configure.outStream << std::endl;

    // Handle overwriting existing samples.mcmc file if requested
    if (mcmcParams.overwriteSamples) {
      std::string samplesFile = configure.outputdir + "samples.mcmc";
      std::ifstream checkFile(configure.OpenPath(samplesFile).c_str());
      if (checkFile.good()) {
        checkFile.close();
        if (std::remove(configure.OpenPath(samplesFile).c_str()) == 0) {
          configure.outStream << "Removed existing samples.mcmc file." << std::endl;
        } else {
          configure.outStream << "Warning: Could not remove existing samples.mcmc file." << std::endl;
        }
      }
    }

    configure.outStream << std::endl;
    configure.outStream << "MCMC Configuration:" << std::endl;
    configure.outStream << "  Number of walkers: " << mcmcParams.nwalkers << std::endl;
    configure.outStream << "  Number of steps: " << mcmcParams.nsteps << std::endl;
    configure.outStream << "  Initial parameter spread: " << mcmcParams.chainSpread << "%" << std::endl;
    configure.outStream << "  Level energy spread: " << mcmcParams.energySpreadKeV << " keV" << std::endl;
    configure.outStream << "  Number of threads: " << mcmcParams.nthreads << std::endl;
    configure.outStream << "  Parameter mode: " << (mcmcParams.useRWA ? "RWA" : "Physical") << std::endl;
    configure.outStream << "  Number of free parameters: " << initialParams.size() << std::endl;
    configure.outStream << std::endl;

    // Create MCMC calculator
    AZURECalcMCMC mcmcCalc(data, compound, configure);

    // Set priors in MCMC calculator (same as GUI does in AZUREMCMCThread.cpp)
    if (priorsFound > 0) {
      mcmcCalc.SetPriors(priorMeans, priorStds, usePriors);
      configure.outStream << "Applied " << priorsFound << " Gaussian prior(s) to MCMC calculator." << std::endl;
      configure.outStream << std::endl;
    }

    // The console's progress line, one per step (the sampler calls the progress callback only when an iteration
    // callback is set).
    const time_t mcmcStart = time(NULL);
    int lastPrintedStep = -1;
    mcmcCalc.SetIterationCallback([](int, int) {});
    mcmcCalc.SetProgressCallback([mcmcStart, &lastPrintedStep](int step, int total, double logprob, double loglik,
                                                               double logprior) {
      if (step == lastPrintedStep) return;
      lastPrintedStep = step;
      PrintMcmcProgress(mcmcStart, step, total, logprob, loglik, logprior);
    });

    // Run MCMC sampling
    std::vector<std::vector<double>> samples;
    const bool sampled = mcmcCalc.RunMCMCSampling(mcmcParams.nwalkers, mcmcParams.nsteps, initialParams,
                                                  samples, mcmcParams.chainSpread, mcmcParams.nthreads,
                                                  mcmcParams.useRWA, mcmcParams.energySpreadKeV);

    if (!sampled) {
      if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) CleanupECAmplitudeCache(configure);
      delete compound;
      delete data;
      return -1;
    }

    configure.outStream << std::endl
                        << std::endl;
    configure.outStream << "MCMC sampling completed successfully!" << std::endl;
    configure.outStream << "Results saved to: " << configure.outputdir << "samples.mcmc" << std::endl;
    configure.outStream << "Total samples generated: " << samples.size() << std::endl;

    // Cleanup EC amplitude cache if it was used
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      if (RunECAmplitudes(configure)) configure.outStream << "EC amplitude cache statistics:" << std::endl;
      CleanupECAmplitudeCache(configure);
    }

    // Cleanup
    delete compound;
    delete data;

    return 0;
  } catch (const std::exception &e) {
    configure.outStream << "Error during MCMC sampling: " << e.what() << std::endl;
    return -1;
  } catch (...) {
    configure.outStream << "Unknown error during MCMC sampling." << std::endl;
    return -1;
  }
}
#endif

/*!
 * This is the main function. All the above initialization functions are called from here.
 * The AZUREMain oject is created and called.
 */

int main(int argc, char *argv[]) {
  // Check for --help option first.  If set, print help and exit.
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--help") == 0) {
      printHelp();
      return 0;
    }
  // The process's environment switches with the command line's switch options, made once, before anything reads them
  // (azure2/support/Environment.h); then the threads (--threads=N) before any parallel region.  A switch option without
  // its value, a checked switch with a bad value, or a bad --threads: exit 2.
  const std::vector<std::string> arguments(argv + 1, argv + argc);
  azure2::ThreadsInForce threads;
  try {
    const azure2::support::Environment &environment = azure2::support::StartProcessEnvironment(arguments);
    azure2::support::CheckSwitchValues(environment);
    threads = azure2::ApplyThreadsOption(arguments, environment);
  } catch (const std::invalid_argument &e) {
    return azure2::OptionError(e.what());
  }

  // Set GSL Error Handler
  gsl_set_error_handler(&GSLException::GSLErrorHandler);

  // The project-file commands (src/ProjectCommands.cpp: --convert, --validate, --schema, ...), each exits
  if (int projectCommandExit = 0; RunProjectCommand(argc, argv, threads, projectCommandExit)) return projectCommandExit;

  // --build-coulomb-table[=FILE] builds the table and exits
  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "--build-coulomb-table") return buildCoulombTable("");
    if (arg == "--build-coulomb-table=") return azure2::OptionError("--build-coulomb-table= needs a FILE");
    if (arg.compare(0, 22, "--build-coulomb-table=") == 0) return buildCoulombTable(arg.substr(22));
  }

  // Create new configuration structure, and parse the command line parameters
  Config configure(std::cout);
  InitializeECAmplitudeCache(configure);
  // An unknown option, a value-taking flag without its value or with one that does not parse, or a bad
  // AZURE_REFINE_QUADRATURE: exit 2
  bool useReadline = false;
  try {
    useReadline = parseOptions(argc, argv, configure);
  } catch (const std::invalid_argument &e) {
    return azure2::OptionError(e.what());
  }

  // The TCP/IP API server is gone; Python access is the in-process _azure2
  // module.  Refuse rather than fall through to the interactive shell, which
  // an old script's closed stdin would spin in forever.
  if (configure.paramMask & Config::USE_API) {
    configure.outStream << "ERROR: --use-api has been removed.  AZURE2 no longer serves a socket API;" << std::endl
                        << "       drive it from Python instead:  from pyazr import azure2" << std::endl;
    return 1;
  }

  // Read the parameters from the runtime configuration file
  if (configure.configfile.empty()) {
    configure.outStream << "A valid configuration file must be specified." << std::endl
                        << "\tSyntax: AZURE2 <options> configfile" << std::endl;
    return -1;
  }
  int configStatus = configure.ReadConfigFile();
  if (configStatus == Config::kConfigUnreadable) {
    configure.outStream << "Could not open " << configure.configfile << ".  Check that file exists."
                        << std::endl;
    return -1;
  } else if (configStatus == Config::kConfigNotFormat3) {
    return 3;  // not a format-3 project (docs/upgrade/project_format3.md section 8.4): the refusal is the last line
  } else if (configStatus != 0) {
    return 1;  // a refused format-3 project or a malformed block, like a refused set-up: its diagnostics are printed
  }
#ifndef NO_STAT
  else if (configure.CheckForInputFiles() == -1)
    return -1;
#endif
  // Joint THM+direct fits couple an arbitrary HOES scale (~1e-4) to widths in
  // eV; MIGRAD's own numerical differences stall far from the minimum on that
  // scaling (7Li(p,a) Tumino+Rolfs benchmark: default MIGRAD stops at
  // chi2 ~1057 with the direct norm dragged to 0.18, while the hybrid
  // analytic/THM-finite-difference gradient reaches the physical minimum,
  // chi2 47.9 with the direct norm at 1.000).  A THM segment is an active
  // data set with a THM observable (its spec), so enable the analytic
  // gradient for such projects.
  if (!(configure.paramMask & Config::USE_ANALYTIC_GRADIENT)) {
    const bool hasTHMSegment = ProjectHasThmData(configure);  // an active data set with a THM observable
    if (hasTHMSegment) {
      configure.paramMask |= Config::USE_ANALYTIC_GRADIENT;
      configure.outStream << "THM segments present: enabling the analytic gradient "
                             "(numerical MIGRAD stalls on joint THM+direct fits)."
                          << std::endl;
    }
  }
  // Before 2026-11-01 Brune OR ignore-externals forced the A-matrix ("--ignore-externals and --use-brune options
  // require A-Matrix"), and as ignore-externals is on by default with no flag to clear it the CLI could never use
  // the R-matrix route.  Brune has been implemented on the R-matrix route (RMatrixFunc; equal to the A-matrix route
  // to ~1e-13); ignore-externals still is not.
  configure.ResolveMatrixRoute();

  // Exit status: 1 when the calculation was aborted or failed ("Calculation was aborted.", a refused project --
  // e.g. a <thm> weight that does not fit the data --, a failed MCMC run), so that scripts and job chains can
  // tell; before 2026-11-01 it was 0 whatever happened on dev, and -1 (255) on thm.
  int returnValue = 0;
  {
    // Print welcome message
    welcomeMessage(configure);
    configure.outStream << "Threads: " << threads.count << " (" << threads.source << ")" << std::endl;

    // Read and process command, setting appropriate configuration flags
    processCommand(commandShell(configure), configure);

    // Open history file for readline
#ifndef NO_READLINE
    if (useReadline) read_history("./.azure_history");
#endif

    // Read the external parameter file to be used, if any
    getParameterFile(useReadline, configure);

    // Parse the segment files for entrance,exit pairs
    std::vector<SegPairs> segPairs;
#ifdef USE_MCMC
    // For MCMC, we still need to check for external capture but use segments from config file
    if (configure.paramMask & Config::PERFORM_MCMC) {
      if (!readSegmentFile(configure, segPairs)) ExitRun(configure, 1);
      if (!checkExternalCapture(configure, segPairs)) ExitRun(configure, 1);
      getExternalCaptureFile(useReadline, configure);
    } else {
#endif
      if (!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
        if (!readSegmentFile(configure, segPairs)) ExitRun(configure, 1);
      } else
        getRateParams(configure, segPairs, useReadline);

      // Check if the entrance,exit pairs are in the external capture file
      //  If so, external capture will be needed
      if (!checkExternalCapture(configure, segPairs)) ExitRun(configure, 1);

      // Read the external capture file name to be used, if any
      getExternalCaptureFile(useReadline, configure);
#ifdef USE_MCMC
    }
#endif

#ifdef USE_MCMC
    // Handle MCMC separately
    if (configure.paramMask & Config::PERFORM_MCMC) {
      MCMCParams mcmcParams;
      getMCMCParams(configure, mcmcParams);

      configure.outStream << std::endl;
      startMessage(configure);
      returnValue = runMCMC(configure, mcmcParams);

      if (returnValue != 0) {
        configure.outStream << "MCMC sampling failed." << std::endl;
      }
    } else {
#endif
      // Create instance of main AZURE function, print start message,
      //  and execute.  On the heap and deliberately never deleted: the process
      //  ends right after this, and freeing the millions of Unit objects one by
      //  one took 2.3-3.7 s of a natural_neon calculation (more with more threads)
      //  after every output file had been written.  No engine object holds a file
      //  or another resource its destructor would release; the OS reclaims the
      //  memory.
      AZUREMain *azureMain = new AZUREMain(configure);
      configure.outStream << std::endl;
      startMessage(configure);
      try {
        returnValue = (*azureMain)();
      } catch (const azure2::MeshPlanBudgetError &e) {  // a plan too large for memory (MeshPlanner.h): refused
        returnValue = (configure.outStream << "ERROR: " << e.what() << "\nCalculation was aborted." << std::endl, 1);
      }
#ifdef USE_MCMC
    }
#endif

    // Print exit message
    exitMessage(configure);
  }

  // Write readline history file
#ifndef NO_READLINE
  if (useReadline) write_history("./.azure_history");
#endif

  // Cleanup caches
  CleanupECAmplitudeCache(configure);

  return returnValue != 0 ? 1 : 0;
}
