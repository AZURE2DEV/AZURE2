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
#include "Fitters.h"
#include "RunDiagnostics.h"
#include "azure2/io/ConfigSection.h"
#include "GSLException.h"
#include "ECAmplitudeCache.h"
#include "CoulombTable.h"
#include "RunPreparation.h"
#include "RunRequest.h"
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
  std::cout << RunCommandsHelp() << std::endl
            << "The options of the program (with any command; the project-file actions --convert ... --schema run "
               "without one):" << std::endl
            << std::setw(25) << std::left << "\t--help:" << std::setw(0) << "Print this list and exit." << std::endl
            << std::setw(25) << std::left << "\t--no-gui:" << std::setw(0) << "Accepted and ignored, so that older scripts keep working;" << std::endl
            << std::setw(25) << std::left << "\t" << std::setw(0) << "this build has no graphical setup utility." << std::endl;
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

void parseOptions(const std::vector<std::string> &arguments, Config &configure) {
  configure.configfile = "";
  std::vector<std::string> options;
  for (const std::string &arg : arguments) {
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
    if (*it == "--no-readline")  // a menu script's option: the menu is gone
      throw std::invalid_argument("Unknown option --no-readline.  AZURE2 has no menu: it takes a command and reads "
                                  "nothing from standard input.\n       Drop --no-readline and give the command: "
                                  "AZURE2 <command> [options] FILE (AZURE2 --help lists them; "
                                  "tools/menu_to_command.py turns an old script's answers into the command).");
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
    else if (ApplyFitterSwitch(configure, *it))  // a fitter's switch: its row of Fitters(); one not built is refused
      continue;
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
    else if (const RunQuestion *question = FindRunQuestion(OptionName(*it)))
      throw std::invalid_argument(CommandOptionWithoutCommand(*question));
    else  // a mistyped flag must not run another calculation: exit 2 with the usage line (main, OptionError)
      throw std::invalid_argument("Unknown option " + *it + ".");
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

int runMCMC(Config &configure, const RunRequest &request) {
  try {
    CNuc *compound = new CNuc();
    WorkingCopy *data = new WorkingCopy();
    // The run's levels (made by the set-up below) and the idle levels of the sampler's walkers.
    LevelState state;
    LevelInputs inputs;
    StatePool<WalkerLevels> walkerPool;

    configure.outStream << "Filling Compound Nucleus..." << std::endl;
    if (compound->Fill(configure, state, inputs) == -1) {
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
      if (const int filled = data->Fill(configure, compound, inputs); filled != 0) {
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
      if (const int made = data->MakePoints(configure, compound, inputs); made != 0) {
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
      compound->Initialize(configure, state, inputs);
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
    compound->FillMnParams(params, state, inputs, &configure);
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
    compound->FillCompoundFromParams(params.GetMinuitParams().Params(), state);
    if (configure.paramMask & Config::USE_BRUNE_FORMALISM) {
      compound->CalcShiftFunctions(configure, state);
    }

    // Extract initial parameter values (only non-fixed parameters)
    std::vector<double> initialParams;
    if (request.reducedWidthAmplitudes) {
      // Use RWA parameters directly
      for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
        if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
          initialParams.push_back(params.GetMinuitParams().Parameter(i).Value());
        }
      }
    } else {
      // Use physical parameters
      vector_r physicalParams = compound->GetTransformParams(state, compound->TransformOut(configure, state, inputs));
      data->Data().Abundances().AppendTo(params.GetMinuitParams().Params(), physicalParams);

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
    if (request.overwriteSamples) {
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
    configure.outStream << "  Number of walkers: " << request.walkers << std::endl;
    configure.outStream << "  Number of steps: " << request.steps << std::endl;
    configure.outStream << "  Initial parameter spread: " << request.parameterSpread << "%" << std::endl;
    configure.outStream << "  Level energy spread: " << request.energySpread << " keV" << std::endl;
    configure.outStream << "  Number of threads: " << request.samplerThreads << std::endl;
    configure.outStream << "  Parameter mode: " << (request.reducedWidthAmplitudes ? "RWA" : "Physical") << std::endl;
    configure.outStream << "  Number of free parameters: " << initialParams.size() << std::endl;
    configure.outStream << std::endl;

    // Create MCMC calculator
    AZURECalcMCMC mcmcCalc(data, compound, state, inputs, walkerPool, configure);

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
    const bool sampled = mcmcCalc.RunMCMCSampling(request.walkers, request.steps, initialParams,
                                                  samples, request.parameterSpread, request.samplerThreads,
                                                  request.reducedWidthAmplitudes, request.energySpread);

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
  // A command (AZURE2 <command> [options] FILE) and its options: a refused command line exits 2 in three parts
  // (src/RunRequest.cpp); the rest of the arguments are the program's options and the project file.
  RunRequest request;
  bool byCommand = false;
  std::vector<std::string> programArguments;
  try {
    byCommand = ReadRunCommand(arguments, request, programArguments);
  } catch (const std::invalid_argument &e) {
    return CommandError(e.what());
  }
  try {
    parseOptions(programArguments, configure);
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
  // A run is a command: without one, refused (exit 2) before the project is read; nothing is read from the input.
  if (!byCommand) return CommandError(NoCommandRefusal());

  // Read the parameters from the runtime configuration file (ReadRunCommand has refused a command without one)
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
  // A file a command's option names, relative to the project's directory: refused before anything is calculated
  if (const std::string refusal = RequestFileRefusal(configure, request); !refusal.empty())
    return CommandError(refusal);
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

    // The command and its answers into the Config
    PrepareRun(configure, request);

#ifdef USE_MCMC
    // Handle MCMC separately
    if (configure.paramMask & Config::PERFORM_MCMC) {
      configure.outStream << std::endl;
      startMessage(configure);
      returnValue = runMCMC(configure, request);

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

  // Cleanup caches
  CleanupECAmplitudeCache(configure);

  return returnValue != 0 ? 1 : 0;
}
