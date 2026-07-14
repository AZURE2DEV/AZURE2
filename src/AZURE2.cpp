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

#include "AZUREMain.h"
#include "Config.h"
#include "NucLine.h"
#include "SegLine.h"
#include "ExtrapLine.h"
#include "GSLException.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cctype>
#include <vector>
#include <map>
#include <tuple>
#include <gsl/gsl_errno.h>
#ifndef NO_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif
#include <string.h>

// Global Config pointer for use across the codebase
Config* g_config = nullptr;

#ifdef GUI_BUILD
extern int start_gui(int argc, char *argv[]);
#endif
struct SegPairs {int firstPair; int secondPair;};

extern int start_api(int port, Config& configure);

/*!
 * This function displays the welcome banner.
 */

void welcomeMessage(const Config& configure) {
  configure.outStream << std::endl
	    << "O--------------------------O-------------------------------O" << std::endl
	    << "| #### #### #  # ###  ####  ##  | Version 1.0              |" << std::endl
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

void exitMessage(const Config& configure) {
  configure.outStream << std::endl
	    << "Thanks for using AZURE2." << std::endl;
}

/*!
 * This function prints the response to the --help command which consists of
 * available runtime options.
 */

void printHelp() {
  std::cout  << "Syntax: AZURE2 <options> configfile" << std::endl << std::endl
	     << "Options:" << std::endl
	     << std::setw(25) << std::left << "\t--no-gui:" << std::setw(0) << "Do not use graphical setup utility (if built)." << std::endl 
	     << std::setw(25) << std::left << "\t" << std::setw(0) << "If this flag is not set all other options are ignored," << std::endl 
	     << std::setw(25) << std::left << "\t" << std::setw(0) << "and configuration occurs within the setup utility." << std::endl 
#ifndef NO_READLINE
	     << std::setw(25) << std::left << "\t--no-readline:" << std::setw(0) << "Do not use readline package." <<  std::endl
#endif
	     << std::setw(25) << std::left << "\t--no-transform:" << std::setw(0) << "Do not perform initial parameter transformations." << std::endl
	     << std::setw(25) << std::left << "\t--no-long-wavelenth:" << std::setw(0) << "Do not use long wavelength approximation for EL capture." << std::endl
	     << std::setw(25) << std::left << "\t--use-brune:" << std::setw(0) << "Use the alternative level matrix of C.R. Brune." << std::endl
	     << std::setw(25) << std::left << "\t--ignore-externals:" << std::setw(0) << "Ignore external resonant capture amplitude if internal width is zero." << std::endl
	     << std::setw(25) << std::left << "\t--use-rmc:" << std::setw(0) << "Use Reich-Moore approximation for capture (neutron capture only)." << std::endl
	     << std::setw(25) << std::left << "\t--gsl-coul:" << std::setw(0) << "Use GSL Coulomb functions (faster, but less accurate)." << std::endl
#ifdef USE_NLOPT
	     << std::setw(25) << std::left << "\t--use-nlopt:" << std::setw(0) << "Use NLopt minimizer (BOBYQA) instead of Minuit2." << std::endl
#endif
	     << std::setw(25) << std::left << "\t--use-gradient:" << std::setw(0) << "Use Minuit2 (MIGRAD) with the analytic gradient (default: numerical)." << std::endl
	     << std::setw(25) << std::left << "\t--use-lm:" << std::setw(0) << "Use the Levenberg-Marquardt minimizer (analytic Jacobian; falls back to MIGRAD)." << std::endl
	     ;
}

/*!
 * This function parses the command line options given, and sets the appropriate variables in the Config
 * structure.  It also reads and parses the configuration file if the appropriate environment
 * variable is set.
 */

bool parseOptions(int argc, char *argv[], Config& configure) {
#ifndef NO_READLINE
  bool useReadline=true;
#else
  bool useReadline = false;
#endif
  configure.configfile="";
  std::vector<std::string> options;
  for(int i=1;i<argc;i++) {
    std::string arg=argv[i];
    if(!arg.empty()&&arg.substr(0,2)=="--") options.push_back(arg);
    else configure.configfile=arg;
  }
  char* optionsFile = getenv("AZURE_OPTIONS_FILE");
  if(optionsFile) {
    std::ifstream in(optionsFile);
    if(in) {
      while(!in.eof()) {
	std::string newOption;
	getline(in,newOption);
	if(!in.eof()) {
	  std::string trimmedOption="";
	  for(std::string::iterator it = newOption.begin();it<newOption.end();it++) 
	    if(*it!=' '&&*it!='\t'&&*it!='\n')
	      trimmedOption+=*it;
	  if(!trimmedOption.empty()&&trimmedOption.substr(0,2)=="--") options.push_back(trimmedOption);
	}
      }
      in.close();
    } else configure.outStream << "AZURE_OPTIONS_FILE variable set, but file not readable." << std::endl;
  }
  for(std::vector<std::string>::iterator it = options.begin();it<options.end();it++) {
#ifndef NO_READLINE
    if(*it=="--no-readline") useReadline=false;
    else if(*it=="--no-transform") configure.paramMask &= ~Config::TRANSFORM_PARAMETERS;
#else
    if (*it == "--no-transform") configure.paramMask &= ~Config::TRANSFORM_PARAMETERS;
#endif
    else if(*it=="--no-long-wavelength") configure.paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
    else if(*it=="--use-brune") configure.paramMask |= Config::USE_BRUNE_FORMALISM;
    else if(*it=="--gsl-coul") configure.paramMask |= Config::USE_GSL_COULOMB_FUNC;
    else if(*it=="--ignore-externals") configure.paramMask |= Config::IGNORE_ZERO_WIDTHS;
    else if(*it=="--use-rmc") configure.paramMask |= Config::USE_RMC_FORMALISM;
    else if(*it=="--use-api") configure.paramMask |= Config::USE_API;
#ifdef USE_NLOPT
    else if(*it=="--use-nlopt") configure.paramMask |= Config::USE_NLOPT_MINIMIZER;
#endif
    else if(*it=="--use-lm") configure.paramMask |= Config::USE_LM_MINIMIZER;
    else if(*it=="--use-gradient") configure.paramMask |= Config::USE_ANALYTIC_GRADIENT;
    else if(*it=="--covariance-band") configure.paramMask |= Config::CALCULATE_COVARIANCE_BAND;
    else if(*it=="--scale-covariance") configure.paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
    else if(*it=="--no-gui") continue;
    else configure.outStream << "WARNING: Unknown option " << *it << '.' << std::endl;
  }
  return useReadline;
}

/*!
 * This function handles the command shell in AZURE2.  The function will not terminate until the user
 * enters a valid integer option.  Upon successful entry, the integer option is returned.
 */

int commandShell(const Config& configure) {
  int command=0;

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
  while(command<1||command>7) {
#else
  while(command<1||command>6) {
#endif
    configure.outStream << "azure2: ";
    std::string inString;
    getline(std::cin,inString);
    if(inString.empty()) continue;
    std::istringstream in;
    in.str(inString);
    if(!(in>>command))
      configure.outStream << "Please enter an integer." << std::endl;
#ifdef USE_MCMC
    else if(command<1||command>7)
#else
    else if(command<1||command>6)
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
static bool askYesNo(Config& configure, const std::string& prompt) {
  configure.outStream << prompt << " (y/n): ";
  std::string in;
  getline(std::cin, in);
  size_t a = in.find_first_not_of(" \t\r\n");
  if(a == std::string::npos) return false;
  char c = (char)std::tolower(in[a]);
  return (c=='y' || c=='1' || c=='t');
}

// Offer the cross-section uncertainty band (and, for fit modes, the reduced-chi^2
// covariance scaling).  If --covariance-band was passed the choice is already
// made non-interactively, so no prompt appears -- this keeps scripted/headless
// runs deterministic while interactive runs are asked.
static void promptBandOptions(Config& configure, bool fitMode) {
  if(configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) return;   // set via flag
  if(!askYesNo(configure, "Calculate cross-section uncertainty band?")) return;
  configure.paramMask |= Config::CALCULATE_COVARIANCE_BAND;
  if(fitMode && !(configure.paramMask & Config::SCALE_COVARIANCE_BY_CHI2)) {
    if(askYesNo(configure, "Scale covariance to reduced chi-squared = 1?"))
      configure.paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  }
}

void processCommand(int command, Config& configure) {
  if(command==2) {
    configure.paramMask |= Config::PERFORM_FIT;
    promptBandOptions(configure, true);
  }
  else if(command==3) {
    configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
    promptBandOptions(configure, false);   // extrapolation reuses a saved covariance
  }
  else if(command==4) {
    bool goodAnswer=false;
    while (!goodAnswer) {
      configure.outStream << std::setw(30) << "Allowed Chi-Squared Variance: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm>>configure.chiVariance) || configure.chiVariance<0.)
	configure.outStream << "Please enter a positive number." << std::endl;
      else goodAnswer=true;
    }
    configure.paramMask |= Config::PERFORM_FIT;
    configure.paramMask |= Config::PERFORM_ERROR_ANALYSIS;
    promptBandOptions(configure, true);
  } else if(command==5) {
    configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
    configure.paramMask |= Config::CALCULATE_REACTION_RATE;
#ifdef USE_MCMC
  } else if(command==6) {
    configure.paramMask |= Config::PERFORM_MCMC;
  } else if(command==7) {
    exitMessage(configure);
    exit(0);
  }
#else
  } else if(command==6) {
    exitMessage(configure);
    exit(0);
  }
#endif
}

/*!
 * This function prompts for a parameter file and sets the corresponding configure
 * flags and variables based on the user response.
 */

void getParameterFile(bool useReadline, Config& configure) {
  bool validInfile=false;
  configure.outStream << std::endl;
  if(!useReadline) configure.outStream << "External Parameter File (leave blank for new file): ";
  while(!validInfile) {
    std::string inFile;
    if(!useReadline) getline(std::cin,inFile);
#ifndef NO_READLINE
    else {
      char *line = readline("External Parameter File (leave blank for new file): ");
      inFile=line;
      size_t endpos = inFile.find_last_not_of(" \t");
      if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
      if(line && *line) add_history(line);
      free(line);
    }
#endif
    if(!inFile.empty()) {
      std::ifstream in;
      in.open(inFile.c_str());
      if(in) {
	validInfile=true;
	configure.paramMask |= Config::USE_PREVIOUS_PARAMETERS;
	configure.paramfile=inFile;
	in.close();
      }
      in.clear();
    } else validInfile=true;
    if(!validInfile) {
      configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
      if(!useReadline) configure.outStream << "External Parameter File (leave blank for new file): ";
    }
  }
}

/*!
 * This function reads the segment file, and stores the active entrance and exit pair keys
 * for cross reference with the External capture file.  Only if an active external capture segment
 * is required is the user prompted for an external integrals file.
 */

bool readSegmentFile(const Config& configure,std::vector<SegPairs>& segPairs) {
  std::ifstream in;
  std::string startTag,stopTag;
  if(configure.paramMask & Config::CALCULATE_WITH_DATA) {
    startTag="<segmentsData>";
    stopTag="</segmentsData>";
  } else {
    startTag="<segmentsTest>";
    stopTag="</segmentsTest>";    
  }
  in.open(configure.configfile.c_str());
  if(in) {
    std::string line="";
    while(line!=startTag&&!in.eof()) getline(in,line);
    if(line==startTag) {
      line="";
      while(line!=stopTag&&!in.eof()) {
	getline(in,line);
	bool empty=true;
	for(unsigned int i=0;i<line.size();++i) 
	  if(line[i]!=' '&&line[i]!='\t') {
	    empty=false;
	    break;
	  }
	if(empty==true) continue;
	if(line!=stopTag&&!in.eof()) {
	  std::istringstream stm;
	  stm.str(line);
	  if(configure.paramMask & Config::CALCULATE_WITH_DATA) {
	    SegLine segment(stm);
	    if(!(stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit))&&segment.isActive()==1) {
	      SegPairs tempSet={segment.entranceKey(),segment.exitKey()};
	      segPairs.push_back(tempSet);
	    }
	  } else {
	    ExtrapLine segment(stm);
	    if(!(stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit))&&segment.isActive()==1) {
	      SegPairs tempSet={segment.entranceKey(),segment.exitKey()};
	      segPairs.push_back(tempSet);
	    }
	  }
	}
      }
      if(line!=stopTag) {
	configure.outStream << "Problem reading segments. Check configuration file." << std::endl;
	return false;
      }
    } else {
      configure.outStream << "Problem reading segments. Check configuration file." << std::endl;
      return false;
    }
    in.close();
  } else {
    configure.outStream << "Cannot read segments. Check configuration file." << std::endl;
    return false;
  }
  in.clear();
  return true;
}

/*!
 * If reaction rate is desired, the user may specify a file containing temperatures for the calculation.
 * This function prompts for that file name, checks for access, and stores path.
 */

void getTemperatureFile(bool useReadline, Config& configure) {
  bool validInfile=false;
  if(!useReadline) configure.outStream << std::setw(38) << "Temperature File Name: ";
  while(!validInfile) {
    std::string inFile;
    if(!useReadline) getline(std::cin,inFile);
#ifndef NO_READLINE
    else {
      char *line = readline("               Temperature File Name: ");
      inFile=line;
      size_t endpos = inFile.find_last_not_of(" \t");
      if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
      if(line && *line) add_history(line);
      free(line);
    }
#endif
    if(!inFile.empty()) {
      std::ifstream in;
      in.open(inFile.c_str());
      if(in) {
	validInfile=true;
	configure.rateParams.temperatureFile=inFile;
	in.close();
      }
      in.clear();
    }
    if(!validInfile) {
      if(inFile.empty()) configure.outStream << "Please enter a file name." << std::endl;
      else configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
      if(!useReadline) configure.outStream << "               Temperature File Name: ";
    }
  }
}

/*!
 * This function prompts the user for MCMC sampling parameters.
 */

#ifdef USE_MCMC
struct MCMCParams {
  int nwalkers;
  int nsteps;
  double chainSpread;
  int nthreads;
  bool useRWA;
  bool overwriteSamples;
};

void getMCMCParams(Config& configure, MCMCParams& mcmcParams) {
  mcmcParams.nwalkers = 0;
  mcmcParams.nsteps = 0;
  mcmcParams.chainSpread = -1.0;  // Will be prompted
  mcmcParams.nthreads = 1;
  mcmcParams.useRWA = false;

  // Get number of walkers
  while(mcmcParams.nwalkers < 2 || mcmcParams.nwalkers > 1000) {
    configure.outStream << std::setw(38) << "Number of Walkers (2-1000): ";
    std::string inString;
    getline(std::cin, inString);
    std::istringstream stm;
    stm.str(inString);
    if(!(stm >> mcmcParams.nwalkers) || mcmcParams.nwalkers < 2 || mcmcParams.nwalkers > 1000)
      configure.outStream << "Please enter an integer between 2 and 1000." << std::endl;
  }

  // Get number of steps
  while(mcmcParams.nsteps < 100 || mcmcParams.nsteps > 10000000) {
    configure.outStream << std::setw(38) << "Number of Steps (100-10000000): ";
    std::string inString;
    getline(std::cin, inString);
    std::istringstream stm;
    stm.str(inString);
    if(!(stm >> mcmcParams.nsteps) || mcmcParams.nsteps < 100 || mcmcParams.nsteps > 10000000)
      configure.outStream << "Please enter an integer between 100 and 10000000." << std::endl;
  }

  // Get initial parameter spread percentage
  while(mcmcParams.chainSpread <= 0.0 || mcmcParams.chainSpread > 50.0) {
    configure.outStream << std::setw(38) << "Initial Parameter Spread % (0.001-50) [default=5]: ";
    std::string inString;
    getline(std::cin, inString);

    // Allow empty input for default
    if(inString.empty() || inString.find_first_not_of(" \t\n") == std::string::npos) {
      mcmcParams.chainSpread = 5.0;
      break;
    }

    std::istringstream stm;
    stm.str(inString);
    if(!(stm >> mcmcParams.chainSpread) || mcmcParams.chainSpread <= 0.0 || mcmcParams.chainSpread > 50.0)
      configure.outStream << "Please enter a number between 0.001 and 50, or press Enter for default (5%)." << std::endl;
  }

  // Get number of threads
  while(mcmcParams.nthreads < 1) {
    configure.outStream << std::setw(38) << "Number of Threads (1 or more): ";
    std::string inString;
    getline(std::cin, inString);
    std::istringstream stm;
    stm.str(inString);
    if(!(stm >> mcmcParams.nthreads) || mcmcParams.nthreads < 1)
      configure.outStream << "Please enter a positive integer." << std::endl;
  }

  // Ask if using RWA parameters
  bool goodAnswer = false;
  std::string rwaAnswer = "";
  while(!goodAnswer) {
    configure.outStream << std::setw(38) << "Use Reduced Width Amplitudes (yes/no): ";
    getline(std::cin, rwaAnswer);
    std::string trimmedAnswer;
    for(int i = 0; i < rwaAnswer.length(); i++)
      if(rwaAnswer[i] != ' ' && rwaAnswer[i] != '\t' && rwaAnswer[i] != '\n')
        trimmedAnswer += rwaAnswer[i];
    if(trimmedAnswer != "yes" && trimmedAnswer != "no")
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
  while(!goodAnswer) {
    configure.outStream << std::setw(38) << "Overwrite existing samples.mcmc? (yes/no) [default=no]: ";
    getline(std::cin, overwriteAnswer);

    // Allow empty input for default (no)
    if(overwriteAnswer.empty() || overwriteAnswer.find_first_not_of(" \t\n") == std::string::npos) {
      mcmcParams.overwriteSamples = false;
      break;
    }

    std::string trimmedAnswer;
    for(int i = 0; i < overwriteAnswer.length(); i++)
      if(overwriteAnswer[i] != ' ' && overwriteAnswer[i] != '\t' && overwriteAnswer[i] != '\n')
        trimmedAnswer += overwriteAnswer[i];
    if(trimmedAnswer != "yes" && trimmedAnswer != "no")
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

void getRateParams(Config& configure, std::vector<SegPairs>& segPairs,bool useReadline) {
  configure.rateParams.entrancePair=0;
  configure.rateParams.exitPair=0;
  configure.rateParams.minTemp=-1.;
  configure.rateParams.maxTemp=-1.;
  configure.rateParams.tempStep=-1.;
  while(configure.rateParams.entrancePair==configure.rateParams.exitPair){
    while(!configure.rateParams.entrancePair) {
      configure.outStream << std::setw(38) << "Reaction Rate Entrance Pair: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm >> configure.rateParams.entrancePair) || configure.rateParams.entrancePair==0) 
	configure.outStream << "Please enter an integer greater than zero." << std::endl ;
    }
    while(!configure.rateParams.exitPair) {
      configure.outStream << std::setw(38) << "Reaction Rate Exit Pair: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm >> configure.rateParams.exitPair) || configure.rateParams.exitPair==0) 
	configure.outStream << "Please enter an integer greater than zero." << std::endl;
    }
    if(configure.rateParams.entrancePair==configure.rateParams.exitPair) {
      configure.outStream << "Cannot calculate rate for elastic scattering." << std::endl;
      configure.rateParams.entrancePair=0;configure.rateParams.exitPair=0;
    }
  }
  SegPairs tempSet={configure.rateParams.entrancePair,configure.rateParams.exitPair};
  segPairs.push_back(tempSet);
  bool goodAnswer = false;
  std::string fileAnswer="";
  while(!goodAnswer) {   
    configure.outStream << std::setw(38) << "Use temperatures from file (yes/no): ";
    getline(std::cin,fileAnswer);
    std::string trimmedAnswer;
    for(int i =0;i<fileAnswer.length();i++)
      if(fileAnswer[i]!=' '&&fileAnswer[i]!='\t'&&fileAnswer[i]!='\n')
	trimmedAnswer+=fileAnswer[i];
    if(trimmedAnswer!="yes"&&trimmedAnswer!="no") 
     configure.outStream << "Please type 'yes' or 'no'." <<std::endl;
    else {
      goodAnswer = true;
      fileAnswer=trimmedAnswer;
    }
  }
  configure.rateParams.useFile = (fileAnswer=="yes") ? (true) : (false); 
  if(configure.rateParams.useFile) getTemperatureFile(useReadline,configure);
  else {
    while(configure.rateParams.minTemp<0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Min Temp [GK]: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm >> configure.rateParams.minTemp) || configure.rateParams.minTemp<0.) 
	configure.outStream << "Please enter a positive number." << std::endl;
    }
    while(configure.rateParams.maxTemp<0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Max Temp [GK]: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm >> configure.rateParams.maxTemp) || configure.rateParams.maxTemp<0.) 
	configure.outStream << "Please enter a positive number." << std::endl;
    }
    while(configure.rateParams.tempStep<0.) {
      configure.outStream << std::setw(38) << "Reaction Rate Temp Step [GK]: ";
      std::string inString;
      getline(std::cin,inString);
      std::istringstream stm;
      stm.str(inString);
      if(!(stm >> configure.rateParams.tempStep) || configure.rateParams.tempStep<0.) 
	configure.outStream << "Please enter a positive number." << std::endl;
    }
  }
}

/*!
 * This function checks the external capture file against a vector of segment key 
 * pairs. Only if the calculation includes external capture segments is the user
 * prompted for an integrals file.  The appropriate configure flag is set here.
 */

bool checkExternalCapture(Config& configure, const std::vector<SegPairs>& segPairs) {
  std::ifstream in;
  in.open(configure.configfile.c_str());
  if(!in) {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;
  }
  std::string line="";
  while(line!="<levels>"&&!in.eof()) getline(in,line);
  if(line!="<levels>") {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;    
  }
  line="";
  while(line!="</levels>"&&!in.eof()&&
	!(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) {
    getline(in,line);
    bool empty=true;
    for(unsigned int i=0;i<line.size();++i) 
      if(line[i]!=' '&&line[i]!='\t') {
	      empty=false;
	      break;
      }
    if(empty==true) continue;
    if(line!="</levels>"&&!in.eof()) {
      std::istringstream stm;
      stm.str(line);
      NucLine tempNucLine(stm);
      if((stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit))) {
	      configure.outStream << "Problem reading nuclear information. Check configuration file." << std::endl;
	      return false;
      }
      if(tempNucLine.ecMultMask()!=0) {
	      for(int i=0;i<segPairs.size();i++) {
	        if(tempNucLine.ir()==segPairs[i].secondPair || segPairs[i].secondPair==-1) {
	          configure.paramMask |= Config::USE_EXTERNAL_CAPTURE;
	          break;
	        }
	      }
      }
    }
  }
  in.close();
  in.clear();
  if((configure.paramMask & Config::USE_EXTERNAL_CAPTURE)&& (configure.paramMask & Config::USE_RMC_FORMALISM)) {
    configure.outStream << "WARNING: External capture is not compatible with Reich-Moore.  Ignoring external capture." << std::endl;
    configure.paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  }
  return true;
}

/*!
 * This function promps the user for an external capture integrals file, checks for it's validity, 
 * and stores the path. 
 */

void getExternalCaptureFile(bool useReadline, Config& configure) {
  if((configure.paramMask & Config::USE_EXTERNAL_CAPTURE)&&
     !(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
    configure.outStream << std::endl;
    if(!useReadline) configure.outStream << "External Capture Amplitude File (leave blank for new file): ";
    bool validInfile=false;
    while(!validInfile) {
      std::string inFile;
      if(!useReadline) getline(std::cin,inFile);
#ifndef NO_READLINE
      else {
		char *line = readline("External Capture Amplitude File (leave blank for new file): ");
		inFile=line;
		size_t endpos = inFile.find_last_not_of(" \t");	
		if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
		if(line && *line) add_history(line);
		free(line);
      }
#endif
      if(!inFile.empty()) {
	std::ifstream in;
	in.open(inFile.c_str());
	if(in) {
	  validInfile=true;
	  configure.paramMask |= Config::USE_PREVIOUS_INTEGRALS;
	  configure.integralsfile=inFile;
	  in.close();
	}
	in.clear();
      } else validInfile=true;
      if(!validInfile) {
	configure.outStream << "Cannot Read From " << inFile << ". Please reenter file." << std::endl;
	if(!useReadline) configure.outStream << "External Capture Amplitude File (leave blank for new file): ";
      }
    }
  } 
}

/*!
 * This function prints a brief start message describing the type of calculation that will be performed.
 */

void startMessage(const Config& configure) {
  if(configure.paramMask & Config::PERFORM_ERROR_ANALYSIS)
    configure.outStream << "Calling AZURE2 for MINOS Error Analysis..." << std::endl;
#ifdef USE_MCMC
  else if(configure.paramMask & Config::PERFORM_MCMC)
    configure.outStream << "Calling AZURE2 for MCMC Bayesian sampling..." << std::endl;
#endif
  else if(configure.paramMask & Config::CALCULATE_REACTION_RATE)
    configure.outStream << "Calling AZURE2 for reaction rate calculation..." << std::endl;
  else if(!(configure.paramMask & Config::CALCULATE_WITH_DATA))
    configure.outStream << "Calling AZURE2 for calculation of segments without data..." << std::endl;
  else if(configure.paramMask & Config::PERFORM_FIT)
    configure.outStream << "Calling AZURE2 for fitting of segments from data..." << std::endl;
  else  configure.outStream << "Calling AZURE2 for calculation of segments from data..." << std::endl;
}

/*!
 * This function runs the MCMC sampling procedure.
 */

#ifdef USE_MCMC
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#include "AZURECalcMCMC.h"
#include <iomanip>
#include <ctime>

// CLI progress callback - prints progress to stdout
static int g_last_printed_step = -1;
static time_t g_start_time = 0;

void cli_progress_callback(int current_step, int total_steps, double logprob, double loglikelihood, double logprior) {
  // Print progress every step (but the iteration callback controls frequency)
  if(current_step != g_last_printed_step) {
    g_last_printed_step = current_step;

    // Calculate time estimates
    time_t current_time = time(NULL);
    double elapsed_seconds = difftime(current_time, g_start_time);
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
}

void cli_iteration_callback(int current_iteration, int total_iterations) {
  // This gets called every iteration - we only print every 10 steps to avoid too much output
  if(current_iteration % 10 == 0 || current_iteration == total_iterations) {
    // The progress callback will handle the actual printing
  }
}

// Helper function to create pretty parameter names (same logic as MCMCTab::createPrettyParameterName)
std::string createPrettyParameterName(const std::string& paramName) {
  // Handle different parameter naming conventions from AZURE2

  // Convention 1: Short form (from .sav files)
  // Energy parameters: E1, E2, etc. -> "Level X Energy (MeV)"
  if(paramName.length() > 1 && paramName[0] == 'E' && std::isdigit(paramName[1])) {
    try {
      int levelNum = std::stoi(paramName.substr(1));
      if(levelNum > 0) {
        return "Level " + std::to_string(levelNum) + " Energy (MeV)";
      }
    } catch(...) {}
  }

  // Width parameters: G11, G12, G21, etc. -> "Level X Channel Y Width (eV)"
  if(paramName.length() > 2 && paramName[0] == 'G' && std::isdigit(paramName[1])) {
    try {
      std::string numPart = paramName.substr(1);
      if(numPart.length() >= 2) {
        int levelNum = std::stoi(numPart.substr(0, 1));
        int channelNum = std::stoi(numPart.substr(1));
        if(levelNum > 0 && channelNum > 0) {
          return "Level " + std::to_string(levelNum) + " Channel " + std::to_string(channelNum) + " Width (eV)";
        }
      }
    } catch(...) {}
  }

  // Normalization: N1, N2, etc. -> "Segment X Normalization"
  if(paramName.length() > 1 && paramName[0] == 'N' && std::isdigit(paramName[1])) {
    try {
      int segmentNum = std::stoi(paramName.substr(1));
      if(segmentNum > 0) {
        return "Segment " + std::to_string(segmentNum) + " Normalization";
      }
    } catch(...) {}
  }

  // Energy shift: S1, S2, etc. -> "Segment X Energy Shift (keV)"
  if(paramName.length() > 1 && paramName[0] == 'S' && std::isdigit(paramName[1])) {
    try {
      int segmentNum = std::stoi(paramName.substr(1));
      if(segmentNum > 0) {
        return "Segment " + std::to_string(segmentNum) + " Energy Shift (keV)";
      }
    } catch(...) {}
  }

  // Convention 2: Descriptive form (from Minuit parameter names)
  // Width parameters: "width_X_Y" -> "Level X Channel Y Width (eV)"
  if(paramName.find("width_") == 0) {
    try {
      size_t firstUnderscore = paramName.find('_');
      size_t secondUnderscore = paramName.find('_', firstUnderscore + 1);
      if(secondUnderscore != std::string::npos) {
        int levelNum = std::stoi(paramName.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1));
        int channelNum = std::stoi(paramName.substr(secondUnderscore + 1));
        return "Level " + std::to_string(levelNum) + " Channel " + std::to_string(channelNum) + " Width (eV)";
      }
    } catch(...) {}
  }

  // Energy parameters: "energy_X" -> "Level X Energy (MeV)"
  if(paramName.find("energy_") == 0) {
    try {
      int levelNum = std::stoi(paramName.substr(7)); // Skip "energy_"
      return "Level " + std::to_string(levelNum) + " Energy (MeV)";
    } catch(...) {}
  }

  // Segment normalization: "segment_X_norm" -> "Segment X Normalization"
  if(paramName.find("segment_") == 0 && paramName.find("_norm") != std::string::npos) {
    try {
      size_t firstUnderscore = paramName.find('_');
      size_t secondUnderscore = paramName.find('_', firstUnderscore + 1);
      if(secondUnderscore != std::string::npos) {
        int segmentNum = std::stoi(paramName.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1));
        return "Segment " + std::to_string(segmentNum) + " Normalization";
      }
    } catch(...) {}
  }

  // Segment energy shift: "segment_X_energy_shift" -> "Segment X Energy Shift (keV)"
  if(paramName.find("segment_") == 0 && paramName.find("_energy_shift") != std::string::npos) {
    try {
      size_t firstUnderscore = paramName.find('_');
      size_t secondUnderscore = paramName.find('_', firstUnderscore + 1);
      if(secondUnderscore != std::string::npos) {
        int segmentNum = std::stoi(paramName.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1));
        return "Segment " + std::to_string(segmentNum) + " Energy Shift (keV)";
      }
    } catch(...) {}
  }

  // Generic fallback - keep original name
  return paramName;
}

int runMCMC(Config& configure, const MCMCParams& mcmcParams) {
  try {
    // Create compound nucleus and experimental data objects
    CNuc* compound = new CNuc();
    EData* data = new EData();

    // Fill compound nucleus from configuration file
    configure.outStream << "Filling Compound Nucleus..." << std::endl;
    if(compound->Fill(configure) == -1) {
      configure.outStream << "Could not fill compound nucleus from file." << std::endl;
      delete compound;
      delete data;
      return -1;
    } else if(compound->NumPairs() == 0 || compound->NumJGroups() == 0) {
      configure.outStream << "No nuclear data exists. Calculation not possible." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Fill the data object from segments and data file
    configure.outStream << "Filling Data Structures..." << std::endl;
    if(configure.paramMask & Config::CALCULATE_WITH_DATA) {
      if(data->Fill(configure, compound) == -1) {
        configure.outStream << "Could not fill data object from file." << std::endl;
        delete compound;
        delete data;
        return -1;
      } else if(data->NumSegments() == 0) {
        configure.outStream << "There is no data provided." << std::endl;
        delete compound;
        delete data;
        return -1;
      }
    } else {
      if(data->MakePoints(configure, compound) == -1) {
        configure.outStream << "Could not fill data object from file." << std::endl;
        delete compound;
        delete data;
        return -1;
      } else if(data->NumSegments() == 0) {
        configure.outStream << "Extrapolation segments produce no data." << std::endl;
        delete compound;
        delete data;
        return -1;
      }
    }

    // Initialize compound nucleus
    try {
      compound->Initialize(configure);
    } catch(GSLException e) {
      configure.outStream << e.what() << std::endl;
      configure.outStream << "Calculation was aborted." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Initialize EC amplitude cache if using external capture
    if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      InitializeECAmplitudeCache();
      // Clear cache to ensure fresh start for this calculation
      if(g_ecAmplitudeCache) {
        g_ecAmplitudeCache->Clear();
      }
      configure.outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
    }

    // Initialize data object
    if(data->Initialize(compound, configure) == -1) {
      configure.outStream << "Calculation was aborted." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Create parameter structure
    AZUREParams params;
    compound->FillMnParams(params.GetMinuitParams(), &configure);
    data->FillMnParams(params.GetMinuitParams());

    // Read existing parameters if provided
    if(configure.paramMask & Config::USE_PREVIOUS_PARAMETERS) {
      configure.outStream << "Reading User Parameter File..." << std::endl;
      params.ReadUserParameters(configure);
    } else {
      configure.outStream << "Creating New param.par File..." << std::endl;
      params.WriteUserParameters(configure, false);
    }

    // Setup compound nucleus from parameters
    compound->FillCompoundFromParams(params.GetMinuitParams().Params());
    if(configure.paramMask & Config::USE_BRUNE_FORMALISM) {
      compound->CalcShiftFunctions(configure);
    }

    // Extract initial parameter values (only non-fixed parameters)
    std::vector<double> initialParams;
    if(mcmcParams.useRWA) {
      // Use RWA parameters directly
      for(int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
        if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
          initialParams.push_back(params.GetMinuitParams().Parameter(i).Value());
        }
      }
    } else {
      // Use physical parameters
      compound->TransformOut(configure);
      vector_r physicalParams = compound->GetTransformParams(configure);

      // Extend with data parameters if needed
      for(size_t i = physicalParams.size(); i < params.GetMinuitParams().Params().size(); i++) {
        physicalParams.push_back(params.GetMinuitParams().Parameter(i).Value());
      }

      // Extract non-fixed parameters
      for(int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
        if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
          initialParams.push_back(physicalParams[i]);
        }
      }
    }

    if(initialParams.empty()) {
      configure.outStream << "Error: No free parameters found for MCMC sampling." << std::endl;
      configure.outStream << "Please ensure at least one parameter is not fixed." << std::endl;
      delete compound;
      delete data;
      return -1;
    }

    // Read priors from AZURE2 file if they exist
    std::map<std::string, std::tuple<double, double, bool>> priorMap; // paramName -> (priorMean, priorStd, usePrior)
    std::ifstream priorFile(configure.configfile.c_str());
    if(priorFile.good()) {
      std::string line;
      bool inMCMCSection = false;
      bool inParametersSection = false;

      while(std::getline(priorFile, line)) {
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if(start != std::string::npos) {
          line = line.substr(start);
          size_t end = line.find_last_not_of(" \t\r\n");
          if(end != std::string::npos) line = line.substr(0, end + 1);
        }

        if(line == "<mcmc>") {
          inMCMCSection = true;
          continue;
        } else if(line == "</mcmc>") {
          inMCMCSection = false;
          break;
        }

        if(inMCMCSection) {
          if(line == "<parameters>") {
            inParametersSection = true;
            continue;
          } else if(line == "</parameters>") {
            inParametersSection = false;
            continue;
          }

          if(inParametersSection && !line.empty()) {
            // Parse: "Parameter (raw_param_name) currentValue priorMean priorStd useGaussianPrior"
            std::istringstream iss(line);
            std::vector<std::string> parts;
            std::string part;
            while(iss >> part) parts.push_back(part);

            if(parts.size() >= 5) {
              // Last 4 parts are: currentValue, priorMean, priorStd, useGaussianPrior
              // Everything before that is the parameter name (possibly "Parameter (raw_name)")
              std::string fileParamName;
              for(size_t i = 0; i < parts.size() - 4; i++) {
                if(i > 0) fileParamName += " ";
                fileParamName += parts[i];
              }

              // Extract raw parameter name from "Parameter (raw_name)" format
              std::string rawParamName;
              size_t openParen = fileParamName.find('(');
              size_t closeParen = fileParamName.find(')');
              if(openParen != std::string::npos && closeParen != std::string::npos && closeParen > openParen) {
                // Extract the raw name from inside parentheses
                rawParamName = fileParamName.substr(openParen + 1, closeParen - openParen - 1);
              } else {
                // No parentheses - use the whole name
                rawParamName = fileParamName;
              }

              // Transform raw name to pretty name
              std::string prettyParamName = createPrettyParameterName(rawParamName);

              double priorMean = std::atof(parts[parts.size() - 3].c_str());
              double priorStd = std::atof(parts[parts.size() - 2].c_str());
              bool usePrior = (std::atoi(parts[parts.size() - 1].c_str()) == 1);

              // Store using the pretty parameter name for lookup
              priorMap[prettyParamName] = std::make_tuple(priorMean, priorStd, usePrior);
            }
          }
        }
      }
      priorFile.close();
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
    for(int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
        std::string rawParamName = params.GetMinuitParams().GetName(i);
        std::string paramName = createPrettyParameterName(rawParamName);
        double currentValue = initialParams[paramIndex];

        // Check if prior exists for this parameter
        std::string priorMeanStr = "N/A";
        std::string priorStdStr = "N/A";
        std::string usePriorStr = "No (uniform)";

        double priorMean = 0.0;
        double priorStd = 0.0;
        bool usePrior = false;

        if(priorMap.find(paramName) != priorMap.end()) {
          priorMean = std::get<0>(priorMap[paramName]);
          priorStd = std::get<1>(priorMap[paramName]);
          usePrior = std::get<2>(priorMap[paramName]);

          // Only mark as found if usePrior is true AND std is positive
          if(usePrior && priorStd > 0.0) {
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
    if(priorsFound > 0) {
      configure.outStream << "Found " << priorsFound << " Gaussian prior(s) in AZURE2 file." << std::endl;
      configure.outStream << "Priors will be applied during MCMC sampling." << std::endl;
    } else {
      configure.outStream << "No priors found in AZURE2 file. Using uniform priors for all parameters." << std::endl;
      configure.outStream << "To use Gaussian priors, use the GUI to configure them and save the .azr file." << std::endl;
    }
    configure.outStream << std::endl;

    // Handle overwriting existing samples.mcmc file if requested
    if(mcmcParams.overwriteSamples) {
      std::string samplesFile = configure.outputdir + "samples.mcmc";
      std::ifstream checkFile(samplesFile.c_str());
      if(checkFile.good()) {
        checkFile.close();
        if(std::remove(samplesFile.c_str()) == 0) {
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
    configure.outStream << "  Number of threads: " << mcmcParams.nthreads << std::endl;
    configure.outStream << "  Parameter mode: " << (mcmcParams.useRWA ? "RWA" : "Physical") << std::endl;
    configure.outStream << "  Number of free parameters: " << initialParams.size() << std::endl;
    configure.outStream << std::endl;

    // Create MCMC calculator
    AZURECalcMCMC mcmcCalc(data, compound, configure);

    // Set priors in MCMC calculator (same as GUI does in AZUREMCMCThread.cpp)
    if(priorsFound > 0) {
      mcmcCalc.SetPriors(priorMeans, priorStds, usePriors);
      configure.outStream << "Applied " << priorsFound << " Gaussian prior(s) to MCMC calculator." << std::endl;
      configure.outStream << std::endl;
    }

    // Register CLI progress callbacks for console output
    g_start_time = time(NULL);
    g_last_printed_step = -1;
    AZURECalcMCMC::SetGUIProgressCallback(cli_progress_callback);
    AZURECalcMCMC::SetGUIIterationCallback(cli_iteration_callback);

    // Run MCMC sampling
    std::vector<std::vector<double>> samples;
    mcmcCalc.RunMCMCSampling(mcmcParams.nwalkers, mcmcParams.nsteps, initialParams,
                             samples, mcmcParams.chainSpread, mcmcParams.nthreads, mcmcParams.useRWA);

    // Clear callbacks
    AZURECalcMCMC::SetGUIProgressCallback(nullptr);
    AZURECalcMCMC::SetGUIIterationCallback(nullptr);

    configure.outStream << std::endl << std::endl;
    configure.outStream << "MCMC sampling completed successfully!" << std::endl;
    configure.outStream << "Results saved to: " << configure.outputdir << "samples.mcmc" << std::endl;
    configure.outStream << "Total samples generated: " << samples.size() << std::endl;

    // Cleanup EC amplitude cache if it was used
    if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      if(g_ecAmplitudeCache) {
        configure.outStream << "EC amplitude cache statistics:" << std::endl;
        // Optionally print cache statistics if available
        // g_ecAmplitudeCache->PrintStats();
      }
      CleanupECAmplitudeCache();
    }

    // Cleanup
    delete compound;
    delete data;

    return 0;
  } catch(const std::exception& e) {
    configure.outStream << "Error during MCMC sampling: " << e.what() << std::endl;
    return -1;
  } catch(...) {
    configure.outStream << "Unknown error during MCMC sampling." << std::endl;
    return -1;
  }
}
#endif

/*!
 * This is the main function. All the above initialization functions are called from here.
 * The AZUREMain oject is created and called.
 */

int main(int argc,char *argv[]){

  //Check for --help option first.  If set, print help and exit.
  for(int i=1;i<argc;i++)
    if(strcmp(argv[i],"--help")==0) {
      printHelp();
      return 0;
    }

  //Set GSL Error Handler
  gsl_set_error_handler (&GSLException::GSLErrorHandler);

  //Initialize caches for performance optimization
  InitializeCoulFuncCache();
  InitializeECAmplitudeCache();
  
  //If GUI is built, look for --no-gui option.  If not set, hand control to GUI.
#ifdef GUI_BUILD
  bool useGUI=true;
  for(int i=1;i<argc;i++) 
    if(strcmp(argv[i],"--no-gui")==0 || strcmp(argv[i],"--use-api")==0) useGUI=false;
  if(useGUI) return start_gui(argc,argv);
#endif

  //Create new configuration structure, and parse the command line parameters
  Config configure(std::cout);
  g_config = &configure; // Set global Config pointer for use across codebase
  bool useReadline = parseOptions(argc,argv,configure);

  //Read the parameters from the runtime configuration file
  if(configure.configfile.empty()) {
    configure.outStream << "A valid configuration file must be specified." << std::endl
	      << "\tSyntax: AZURE2 <options> configfile" << std::endl;
    return -1;
  }
  if(configure.ReadConfigFile()==-1) {
    configure.outStream << "Could not open " << configure.configfile << ".  Check that file exists." 
	      << std::endl;
    return -1;
  }
#ifndef NO_STAT
  else if(configure.CheckForInputFiles() == -1) return -1;
#endif
  if((configure.paramMask & Config::USE_RMC_FORMALISM) && (configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    configure.outStream << "WARNING: --use-brune is incompatible with --use-rmc. Ignoring --use-brune." << std::endl;
    configure.paramMask &= ~Config::USE_BRUNE_FORMALISM;
  }
  if((configure.paramMask & Config::USE_BRUNE_FORMALISM)||(configure.paramMask & Config::IGNORE_ZERO_WIDTHS)) {
    if(!(configure.paramMask & Config::USE_AMATRIX)) 
      configure.outStream << "WARNING: R-Matrix specified but --ignore-externals and --use-brune options require A-Matrix.  A-Matrix will be used."
		<< std::endl;
    configure.paramMask |= Config::USE_AMATRIX;
  }

  int port;
  bool useAPI=false;
  for(int i=1;i<argc;i++){ 
    if(strcmp(argv[i],"--use-api")==0){ useAPI=true;
      port = atoi(argv[i+1]);
    }
  }
  if(useAPI) start_api(port, configure);

  else{

    //Print welcome message
    welcomeMessage(configure);

    //Read and process command, setting appropriate configuration flags
    processCommand(commandShell(configure),configure);
  
    //Open history file for readline
  #ifndef NO_READLINE
    if(useReadline) read_history("./.azure_history");
  #endif
  
    //Read the external parameter file to be used, if any
    getParameterFile(useReadline,configure);
  
    //Parse the segment files for entrance,exit pairs
    std::vector<SegPairs> segPairs;
#ifdef USE_MCMC
    // For MCMC, we still need to check for external capture but use segments from config file
    if(configure.paramMask & Config::PERFORM_MCMC) {
      if(!readSegmentFile(configure,segPairs)) exit(1);
      if(!checkExternalCapture(configure,segPairs)) exit(1);
      getExternalCaptureFile(useReadline,configure);
    } else {
#endif
      if(!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
        if(!readSegmentFile(configure,segPairs)) exit(1);
      } else getRateParams(configure,segPairs,useReadline);

      //Check if the entrance,exit pairs are in the external capture file
      // If so, external capture will be needed
      if(!checkExternalCapture(configure,segPairs)) exit(1);

      //Read the external capture file name to be used, if any
      getExternalCaptureFile(useReadline,configure);
#ifdef USE_MCMC
    }
#endif

#ifdef USE_MCMC
    // Handle MCMC separately
    if(configure.paramMask & Config::PERFORM_MCMC) {
      MCMCParams mcmcParams;
      getMCMCParams(configure, mcmcParams);

      configure.outStream << std::endl;
      startMessage(configure);
      int returnValue = runMCMC(configure, mcmcParams);

      if(returnValue != 0) {
        configure.outStream << "MCMC sampling failed." << std::endl;
      }
    } else {
#endif
      //Create instance of main AZURE function, print start message,
      // and execute
      AZUREMain azureMain(configure);
      configure.outStream << std::endl; startMessage(configure);
      int returnValue = azureMain();
#ifdef USE_MCMC
    }
#endif

    //Print exit message
    exitMessage(configure);

  }

  //Write readline history file
#ifndef NO_READLINE
  if(useReadline) write_history("./.azure_history");
#endif

  //Cleanup caches
  CleanupCoulFuncCache();
  CleanupECAmplitudeCache();

  return 0;
}
