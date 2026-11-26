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
#include "azure2/support/FloatingPointPrecision.h"
#include "ProjectCommands.h"
#include "ProjectSetup.h"
#include "CommandLineValues.h"
#include "Config.h"
#include "Fitters.h"
#include "RunDiagnostics.h"
#include "azure2/io/ConfigSection.h"
#include "GSLException.h"
#include "ECAmplitudeCache.h"
#include "CoulombTable.h"
#include "RunPreparation.h"
#include "RunRequest.h"
#include "ProgramMain.h"
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
#include <optional>
#include <tuple>
#include <gsl/gsl_errno.h>
#include <string.h>
#include <cstdlib>

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
  std::cout << RunCommandsHelp() << std::endl
            << "The options of the program (with any command; the project-file actions --convert ... --schema run "
               "without one):" << std::endl
            << std::setw(25) << std::left << "\t--help:" << std::setw(0) << "Print this list and exit." << std::endl
            << std::setw(25) << std::left << "\t--no-gui:" << std::setw(0) << "Accepted and ignored, so that older scripts keep working;" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "the Qt setup program is AZURESetup, a program of its own." << std::endl;
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
            << std::setw(25) << std::left << "\t--schema-out=FILE:" << std::setw(0) << "As --schema, written to FILE: to a temporary file beside it, renamed over FILE when whole; on an error FILE is unchanged (exit 1)." << std::endl
            << std::setw(25) << std::left << "\t--build-coulomb-table[=FILE]:" << std::setw(0) << "Build the Coulomb table (about 5 min on 4 threads, 88 MB; FILE as for --coulomb-table) and exit." << std::endl
            << std::setw(25) << std::left << "\t--refreeze-covariance[=N]:" << std::setw(0) << "Refit up to N times (default 5) with relvector covariance terms re-frozen at the fitted parameters." << std::endl
            << FitterSwitchesHelp()
            << std::setw(25) << std::left << "\t--covariance-band:" << std::setw(0) << "An option of a command (above); in an options file it applies to every command." << std::endl
            << std::setw(25) << std::left << "\t--scale-covariance:" << std::setw(0) << "An option of a command (above); in an options file it applies to every command." << std::endl
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
  spec.threads = CoulombTable::DefaultBuildThreads();
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
 * This is the main function. All the above initialization functions are called from here.
 * The AZUREMain oject is created and called.
 */

int main(int argc, char *argv[]) {
  // The main thread's x87 precision at 64 bits, whatever the C runtime's start-up set (the OpenMP regions take the
  // same scope in each thread: azure2/support/FloatingPointPrecision.h).
  const azure2::support::ExtendedPrecisionScope fpPrecision;
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

  // The run of the command (src/ProgramMain.cpp): the same function a front end in its own process calls
  // (azure2::front::Run, docs/upgrade/qt_gui.md), with standard output as the run's stream.
  Config configure(std::cout);
  return RunProgram(configure, arguments, threads);
}
