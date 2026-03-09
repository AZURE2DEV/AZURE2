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
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <string.h>

#ifdef GUI_BUILD
extern int start_gui(int argc, char *argv[]);
#endif
struct SegPairs {int firstPair; int secondPair;};

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
	     << std::setw(25) << std::left << "\t--no-readline:" << std::setw(0) << "Do not use readline package." <<  std::endl
	     << std::setw(25) << std::left << "\t--no-transform:" << std::setw(0) << "Do not perform initial parameter transformations." << std::endl
	     << std::setw(25) << std::left << "\t--no-long-wavelenth:" << std::setw(0) << "Do not use long wavelength approximation for EL capture." << std::endl
	     << std::setw(25) << std::left << "\t--use-brune:" << std::setw(0) << "Use the alternative level matrix of C.R. Brune." << std::endl
	     << std::setw(25) << std::left << "\t--ignore-externals:" << std::setw(0) << "Ignore external resonant capture amplitude if internal width is zero." << std::endl
	     << std::setw(25) << std::left << "\t--use-rmc:" << std::setw(0) << "Use Reich-Moore approximation for capture (neutron capture only)." << std::endl
	     << std::setw(25) << std::left << "\t--gsl-coul:" << std::setw(0) << "Use GSL Coulomb functions (faster, but less accurate)." << std::endl;
}

/*!
 * This function parses the command line options given, and sets the appropriate variables in the Config
 * structure.  It also reads and parses the configuration file if the appropriate environment
 * variable is set.
 */

bool parseOptions(int argc, char *argv[], Config& configure) {
  bool useReadline=true;
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
    if(*it=="--no-readline") useReadline=false;
    else if(*it=="--no-transform") configure.paramMask &= ~Config::TRANSFORM_PARAMETERS;
    else if(*it=="--no-long-wavelength") configure.paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
    else if(*it=="--use-brune") configure.paramMask |= Config::USE_BRUNE_FORMALISM;
    else if(*it=="--gsl-coul") configure.paramMask |= Config::USE_GSL_COULOMB_FUNC;
    else if(*it=="--ignore-externals") configure.paramMask |= Config::IGNORE_ZERO_WIDTHS;
    else if(*it=="--use-rmc") configure.paramMask |= Config::USE_RMC_FORMALISM;
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
	    << "\t6. Exit" << std::endl;

  while(command<1||command>6) {
    configure.outStream << "azure2: ";
    std::string inString;
    getline(std::cin,inString);
    if(inString.empty()) continue;
    std::istringstream in;
    in.str(inString);
    if(!(in>>command)) 
      configure.outStream << "Please enter an integer." << std::endl;
    else if(command<1||command>6) 
      configure.outStream << "Invalid option.  Please try again."
		<< std::endl;
  } 
  return command;
}

/*!
 * This function takes the returned option from the commandShell function and
 * sets the appropriate flags in the Config structure.
 */

void processCommand(int command, Config& configure) {
  if(command==2) configure.paramMask |= Config::PERFORM_FIT;
  else if(command==3) configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
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
  } else if(command==5) {
    configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
    configure.paramMask |= Config::CALCULATE_REACTION_RATE;
  } else if(command==6) {
    exitMessage(configure);
    exit(0);
  }
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
    else {
      char *line = readline("External Parameter File (leave blank for new file): ");
      inFile=line;
      size_t endpos = inFile.find_last_not_of(" \t");
      if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
      if(line && *line) add_history(line);
      free(line);
    }
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
    else {
      char *line = readline("               Temperature File Name: ");
      inFile=line;
      size_t endpos = inFile.find_last_not_of(" \t");
      if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
      if(line && *line) add_history(line);
      free(line);
    }
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
	  if(tempNucLine.ir()==segPairs[i].secondPair||
	     segPairs[i].secondPair==-1) {
	    configure.paramMask |= Config::USE_EXTERNAL_CAPTURE;
	    break;
	  }
	}
      }
    }
  }
  in.close();
  in.clear();
  if((configure.paramMask & Config::USE_EXTERNAL_CAPTURE)&&
     (configure.paramMask & Config::USE_RMC_FORMALISM)) {
    configure.outStream << "WARNING: External capture is not compatible with Reich-Moore.  Ignoring external capture." 
			<< std::endl;
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
      else {
	char *line = readline("External Capture Amplitude File (leave blank for new file): ");
	inFile=line;
	size_t endpos = inFile.find_last_not_of(" \t");	
	if( std::string::npos != endpos ) inFile = inFile.substr( 0, endpos+1 );
	if(line && *line) add_history(line);
	free(line);
      }
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
  else if(configure.paramMask & Config::CALCULATE_REACTION_RATE) 
    configure.outStream << "Calling AZURE2 for reaction rate calculation..." << std::endl;  
  else if(!(configure.paramMask & Config::CALCULATE_WITH_DATA)) 
    configure.outStream << "Calling AZURE2 for calculation of segments without data..." << std::endl; 
  else if(configure.paramMask & Config::PERFORM_FIT) 
    configure.outStream << "Calling AZURE2 for fitting of segments from data..." << std::endl;
  else  configure.outStream << "Calling AZURE2 for calculation of segments from data..." << std::endl;
}

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
  
  //If GUI is built, look for --no-gui option.  If not set, hand control to GUI.
#ifdef GUI_BUILD
  bool useGUI=true;
  for(int i=1;i<argc;i++) 
    if(strcmp(argv[i],"--no-gui")==0) useGUI=false;
  if(useGUI) return start_gui(argc,argv);
#endif

  //Create new configuration structure, and parse the command line parameters
  Config configure(std::cout);
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

  //Print welcome message
  welcomeMessage(configure);

  //Read and process command, setting appropriate configuration flags
  processCommand(commandShell(configure),configure);
  
  //Open history file for readline
  if(useReadline) read_history("./.azure_history");
  
  //Read the external parameter file to be used, if any
  getParameterFile(useReadline,configure);
  
  //Parse the segment files for entrance,exit pairs
  std::vector<SegPairs> segPairs;
  if(!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
    if(!readSegmentFile(configure,segPairs)) exit(1);
  } else getRateParams(configure,segPairs,useReadline);

  //Check if the entrance,exit pairs are in the external capture file
  // If so, external capture will be needed
  if(!checkExternalCapture(configure,segPairs)) exit(1);
  
  //Read the external capture file name to be used, if any
  getExternalCaptureFile(useReadline,configure);
   
  //Create instance of main AZURE function, print start message,
  // and execute
  AZUREMain azureMain(configure);
  configure.outStream << std::endl; startMessage(configure);
  int returnValue = azureMain();
  
  //Print exit message
  exitMessage(configure);

  //Write readline history file
  if(useReadline) write_history("./.azure_history");
  
  return returnValue;
}
