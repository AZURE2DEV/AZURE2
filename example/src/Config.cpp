#include "Config.h"
#ifndef NO_STAT
#include <sys/stat.h>
#endif
#include <iostream>

/*!
 * The constructor of the Config class sets defaults and the 
 * stream reference for output.
 */

Config::Config(std::ostream& stream) : outStream(stream) {
  Reset();
}

/*!
 * This function resets Config structure.
 */

void Config::Reset() {
  chiVariance=1.0;
  screenCheckMask=0;
  fileCheckMask=0;
  paramMask=0;
  paramMask |= (USE_AMATRIX|TRANSFORM_PARAMETERS|CALCULATE_WITH_DATA|USE_LONGWAVELENGTH_APPROX);
  stopFlag=false;
  outputdir="";
  checkdir="";
}

/*!
 * This funciton reads the configuration file and parses various options.
 */

int Config::ReadConfigFile() {
  std::string dummy; std::string temp;
  std::ifstream in(configfile.c_str());
  if(!in) return -1;
  std::string line="";
  while(line!="<config>"&&!in.eof()) getline(in,line);
  if(line!="<config>") return -1;
  in >> temp;getline(in,dummy);
  if(temp=="true") paramMask |= USE_AMATRIX;
  else paramMask &= ~USE_AMATRIX;
  getline(in,dummy);
  int poundSignPos=dummy.find_last_of('#');
  if(poundSignPos==std::string::npos) temp=dummy;
  else temp=dummy.substr(0,poundSignPos);
  int p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {  
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    outputdir=temp.substr(p1,(p2-p1)+1);
  } else outputdir=std::string();  
  getline(in,dummy);
  poundSignPos=dummy.find_last_of('#');
  if(poundSignPos==std::string::npos) temp=dummy;
  else temp=dummy.substr(0,poundSignPos);
  p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {  
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    checkdir=temp.substr(p1,(p2-p1)+1);
  } else checkdir=std::string();   
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_COMPOUND_NUCLEUS;
  else if(temp=="file") fileCheckMask |= CHECK_COMPOUND_NUCLEUS;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  else if(temp=="file") fileCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_DATA;
  else if(temp=="file") fileCheckMask |= CHECK_DATA;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_ENERGY_DEP;
  else if(temp=="file") fileCheckMask |= CHECK_ENERGY_DEP;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_LEGENDRE;
  else if(temp=="file") fileCheckMask |= CHECK_LEGENDRE;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_COUL_AMPLITUDES;
  else if(temp=="file") fileCheckMask |= CHECK_COUL_AMPLITUDES;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_PATHWAYS;
  else if(temp=="file") fileCheckMask |= CHECK_PATHWAYS;
  in >> temp;getline(in,dummy);
  if(temp=="screen") screenCheckMask |= CHECK_ANGULAR_DISTS;
  else if(temp=="file") fileCheckMask |= CHECK_ANGULAR_DISTS;
  line="";
  while(line!="</config>"&&!in.eof()) getline(in,line);
  if(line!="</config>") return -1;
  in.close();
  return 0;
}

/*!
 * If stat() is enabled, this function checks for the output and checks 
 * directories at runtime.
 */

#ifndef NO_STAT
int Config::CheckForInputFiles() {
  struct stat buffer;
  if(stat(outputdir.c_str(),&buffer) != 0) {
    outStream << "Could not find output directory: " << outputdir << ". Check that it exists." << std::endl;
    return -1;
  }
  if(stat(checkdir.c_str(),&buffer) != 0) {
    outStream << "Could not find checks directory: " << checkdir << ". Check that it exists." << std::endl;
    return -1;
  }
  return 0;
}
#endif
