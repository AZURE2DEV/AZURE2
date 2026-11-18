// The program's run of one command on a project (docs/upgrade/qt_gui.md section 3): what main does after the process's
// start-up (the environment record, the threads, GSL's error handler, the project-file commands).  The command line's
// main calls it with standard output; a front end in its own process (azure2::front::Run, api/front/) with its own
// stream, so that the two run the same code and write the same files.
#include "ProgramMain.h"

#include "AZUREMain.h"
#include "CommandLineValues.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "Fitters.h"
#include "MeshPlanner.h"
#include "Profile.h"
#include "ProjectSetup.h"
#include "RunPreparation.h"
#include "RunRequest.h"
#include "azure2/io/ConfigSection.h"
#include "azure2/support/Environment.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

/*!
 * This function displays the welcome banner.
 */

static void welcomeMessage(const Config &configure) {
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

static void exitMessage(const Config &configure) {
  configure.outStream << std::endl
                      << "Thanks for using AZURE2." << std::endl;
}

/*!
 * This function parses the command line options given, and sets the appropriate variables in the Config
 * structure.  It also reads and parses the configuration file if the appropriate environment
 * variable is set.  An unknown option, or a value-taking flag without a value that parses, throws
 * std::invalid_argument naming it (include/CommandLineValues.h).
 */

static void parseOptions(const std::vector<std::string> &arguments, Config &configure) {
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
    else if (const azure2::io::RemovedConfigKey *removed = azure2::io::RemovedConfigKeyOfFlag(*it))
      throw std::invalid_argument("Unknown option " + *it + ".  " +  // the parser's words for a flag it does not take
                                  azure2::io::RemovedConfigKeyText(*removed, azure2::io::OptionSpelling::kFlag));
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

static void startMessage(const Config &configure) {
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
#include "CompoundSystem.h"
#include "CompoundSystemSetUp.h"
#include "ProgramRun.h"
#include "WorkingCopy.h"
#include "AZUREParams.h"
#include "SamplerRun.h"
#include <iomanip>
#include <ctime>

// The CLI's progress line of an MCMC run: printed once per step, with an estimate of the time left from the steps
// done since `start`.
static void PrintMcmcProgress(std::ostream &out, time_t start, int current_step, int total_steps, double logprob, double loglikelihood,
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

  out << "\rStep " << std::setw(6) << current_step << "/" << total_steps
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * current_step / total_steps) << "%)"
            << " | logP=" << std::setprecision(2) << logprob
            << " | logL=" << loglikelihood
            << " | logPrior=" << logprior
            << " | ETA: " << std::setfill('0') << std::setw(2) << hours << ":"
            << std::setw(2) << minutes << ":" << std::setw(2) << seconds
            << std::setfill(' ') << std::flush;
}

static int runMCMC(Config &configure, const RunRequest &request) {
  try {
    WorkingCopy *data = new WorkingCopy();

    // The compound system from the <levels> rows, the data set up in the middle, then its last part: as the
    // program's run does it (SetUpProgramCompound), its check files and its words on a refusal included.
    CompoundSetUp setUp = SetUpProgramCompound(configure, *data, DataCheckAt::kDataSetUp);
    if (setUp.refusal != CompoundSetUp::Refusal::kNone) {
      delete data;
      return -1;
    }
    const std::unique_ptr<const CompoundSystem> compound = std::move(setUp.compound);
    LevelState &state = setUp.state;
    LevelInputs &inputs = setUp.inputs;

    // Initialize EC amplitude cache if using external capture
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      InitializeECAmplitudeCache(configure);
      configure.outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
    }

    const int initialized = data->Initialize(compound.get(), configure);  // AZURE_PROFILE: the set-up's stages on their own
    AZUREProfile::Report(configure.profile, std::cerr, "set-up");
    if (initialized == -1) {
      configure.outStream << "Calculation was aborted." << std::endl;
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
      if (!params.ReadUserParameters(configure) || !params.ReconcileBasis(compound.get(), configure)) return -1;
    } else {
      configure.outStream << "Creating New param.par File..." << std::endl;
      params.WriteUserParameters(configure, false);
    }

    // Setup compound nucleus from parameters
    compound->FillCompoundFromParams(params.GetMinuitParams().Params(), state);
    if (configure.paramMask & Config::USE_BRUNE_FORMALISM) {
      compound->CalcShiftFunctions(configure, state);
    }

    // From here on the run is the engine's (SamplerRun.h), the same for a session of pyazr: the start vector, the
    // priors of <mcmc>, the table of the free parameters, the sampler and its run.  The console's progress line, one
    // per step.
    const time_t mcmcStart = time(NULL);
    int lastPrintedStep = -1;
    SamplerHooks hooks;
    hooks.afterStep = [mcmcStart, &lastPrintedStep, &configure](int step, int total, double logprob, double loglik,
                                                                double logprior) {
      // A front end's stop request (ProgramObserver::stop) ends the run after this step; the command line has none.
      if (step == lastPrintedStep) return configure.StopRequested();
      lastPrintedStep = step;
      PrintMcmcProgress(configure.outStream, mcmcStart, step, total, logprob, loglik, logprior);
      return configure.StopRequested();
    };
    std::vector<std::vector<double>> samples;
    const SamplerOutcome outcome = RunSampler(configure, data, compound.get(), state, inputs, params, request, hooks,
                                              samples);
    if (outcome == SamplerOutcome::kRefused) {
      delete data;
      return -1;
    }
    const bool sampled = outcome == SamplerOutcome::kSampled;

    if (!sampled) {
      if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) CleanupECAmplitudeCache(configure);
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
 * The run of a command: the arguments (argv without the program's name), the threads in force.  Returns main's exit
 * value.
 */
int RunProgram(Config &configure, const std::vector<std::string> &arguments, const azure2::ThreadsInForce &threads) {
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
    return CommandError(configure.outStream, e.what());
  }
  try {
    parseOptions(programArguments, configure);
  } catch (const std::invalid_argument &e) {
    return azure2::OptionError(configure.outStream, e.what());
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
  if (!byCommand) return CommandError(configure.outStream, NoCommandRefusal());

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
    return CommandError(configure.outStream, refusal);
  // Joint THM+direct fits couple an arbitrary HOES scale (~1e-4) to widths in
  // eV; MIGRAD's own numerical differences stall far from the minimum on that
  // scaling (7Li(p,a) Tumino+Rolfs benchmark: default MIGRAD stops at
  // chi2 ~1057 with the direct norm dragged to 0.18, while the hybrid
  // analytic/THM-finite-difference gradient reaches the physical minimum,
  // chi2 47.9 with the direct norm at 1.000).  A THM segment is an active
  // data set with a THM observable (its spec), so enable the analytic
  // gradient for such projects: the row migrad-gradient, under the switches' precedence (a switch of an earlier row
  // still wins; docs/upgrade/fit.md section 5).
  if (std::string(FitterOfConfig(configure).name) != "migrad-gradient") {
    const bool hasTHMSegment = ProjectHasThmData(configure);  // an active data set with a THM observable
    if (hasTHMSegment) {
      ChooseFitter(configure, "migrad-gradient");
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
