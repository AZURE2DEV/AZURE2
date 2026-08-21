#include "Config.h"
#include "NuclearPotentialManager.h"
#ifndef NO_STAT
#include <sys/stat.h>
#endif
#include <iostream>
#include <sstream>
#include <stdexcept>

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
  paramMask |= (USE_AMATRIX|USE_BRUNE_FORMALISM|IGNORE_ZERO_WIDTHS|TRANSFORM_PARAMETERS|CALCULATE_WITH_DATA|USE_LONGWAVELENGTH_APPROX);
  stopFlag=false;
  outputdir="";
  checkdir="";
  nloptAlgorithm=0; // Default to SBPLX
  useHybridMethod=false; // Default to disabled
  useAdaptiveGrid=true; // Default to adaptive grid
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
  return this->ReadPotentialBlock();
}

/*!
 * Reads the <potential> block of the configuration file and configures the
 * hybrid Coulomb method accordingly.
 *
 * The block is optional -- a file without one, or one with
 * useHybridPotential=0, leaves the defaults set by Reset() untouched, so this
 * is a no-op for every existing project.  The format is the one written by the
 * setup utility (gui/src/AZURESetup.cpp) and read back by
 * NuclearPotentialTab::readPotentialSettings:
 *
 *   <potential>
 *   useHybridPotential=1
 *   useAdaptiveGrid=1
 *   potentialType=0        # 0 = Woods-Saxon, 1 = Gaussian
 *   V0=80                  # depth, MeV
 *   R=3.6                  # radius, fm      (Woods-Saxon)
 *   a=0.6                  # diffuseness, fm (Woods-Saxon)
 *   r0=5.0                 # width, fm       (Gaussian)
 *   pair=2                 # everything below applies to pair 2 alone
 *   potentialType=1
 *   V0=90
 *   r0=4.0
 *   pair=5
 *   useHybridPotential=0   # pair 5 opts out
 *   </potential>
 *
 * A nuclear potential belongs to a particle pair -- it bends the radial wave
 * functions of that channel and no other -- so the keys are read per pair.
 * Everything before the first pair= is the default, which stands in for every
 * pair that does not name itself; a pair= line opens a section that starts
 * from that default and only has to state what differs.  A file with no pair=
 * line is therefore read exactly as it was when the model was global.
 *
 * useAdaptiveGrid is a property of the target-effect integration rather than
 * of any one channel, so it stays global wherever it appears.
 *
 * Parsing it here rather than in the setup utility is what makes the hybrid
 * model reachable from --no-gui and from pyazr: this function is on the path
 * both of them take.
 *
 * Returns 0 on success (including "no block present") and -1 if the block is
 * present but malformed.
 */

int Config::ReadPotentialBlock() {
  std::ifstream in(configfile.c_str());
  if(!in) return -1;

  std::string line="";
  while(line!="<potential>"&&!in.eof()) getline(in,line);
  if(line!="<potential>") return 0;               // optional block, absent

  NuclearPotentialManager& manager=NuclearPotentialManager::instance();
  manager.resetToDefault();

  NuclearPotentialSetting current=manager.getDefaultSetting();
  int currentPair=0;                              // 0 = the default section
  bool currentHasType=false, defaultHasType=false, closed=false;

  // Install whatever the section just read describes.  A section that asks for
  // the hybrid model without naming a potentialType, or names one that does not
  // exist, is refused rather than silently given the default shape.
  auto flush=[&]() {
    if(current.enabled&&!currentHasType) {
      outStream << "WARNING: <potential> requests the hybrid method";
      if(currentPair) outStream << " for pair " << currentPair;
      outStream << " but gives no potentialType; it is disabled." << std::endl;
      current.enabled=false;
    }
    try {
      if(currentPair) manager.setSetting(currentPair,current);
      else manager.setDefaultSetting(current);
    } catch(const std::exception& e) {
      outStream << "WARNING: invalid potential parameters (" << e.what()
                << "); the hybrid method is disabled";
      if(currentPair) outStream << " for pair " << currentPair;
      outStream << "." << std::endl;
      current.enabled=false;
      if(currentPair) manager.clearPairSetting(currentPair);
      else manager.setDefaultEnabled(false);
    }
  };

  while(!in.eof()) {
    getline(in,line);
    size_t b=line.find_first_not_of(" \t\r\n");
    if(b==std::string::npos) continue;
    size_t e=line.find_last_not_of(" \t\r\n");
    std::string trimmed=line.substr(b,e-b+1);
    if(trimmed=="</potential>") { closed=true; break; }

    size_t eq=trimmed.find('=');
    if(eq==std::string::npos) continue;
    std::string key=trimmed.substr(0,eq);
    std::istringstream value(trimmed.substr(eq+1));

    // A pair= line closes the section being read and opens the next one.  Every
    // section starts from the default, so a pair only has to state what differs.
    if(key=="pair") {
      int nextPair=0;
      value >> nextPair;
      if(nextPair<=0) {
        outStream << "WARNING: <potential> has pair=" << nextPair
                  << ", which is not a pair key; the section is ignored." << std::endl;
        continue;
      }
      flush();
      currentPair=nextPair;
      current=manager.getDefaultSetting();
      currentHasType=defaultHasType;              // inherited with the shape
      continue;
    }

    if(key=="useHybridPotential") { int v=0; value >> v; current.enabled=(v==1); }
    else if(key=="useAdaptiveGrid") { int v=1; value >> v; useAdaptiveGrid=(v==1); }
    else if(key=="potentialType") {
      int typeCode=0;
      value >> typeCode;
      if(typeCode==0) { current.type="WoodsSaxon"; currentHasType=true; }
      else if(typeCode==1) { current.type="Gaussian"; currentHasType=true; }
      else {
        outStream << "WARNING: unknown potentialType " << typeCode
                  << " in <potential>";
        if(currentPair) outStream << " for pair " << currentPair;
        outStream << "; the hybrid method is disabled there." << std::endl;
        current.enabled=false;
        currentHasType=true;                      // refused, not merely missing
      }
      if(!currentPair) defaultHasType=currentHasType;
    }
    else if(key=="V0") value >> current.V0;
    else if(key=="R") value >> current.R;
    else if(key=="a") value >> current.a;
    else if(key=="r0") value >> current.r0;
  }
  in.close();

  if(!closed) return -1;                          // unterminated block
  flush();

  // The global switch stays the master: it gates the paramMask bit and the
  // GUI tab, while the manager decides which pairs the model applies to.
  useHybridMethod=manager.isAnyEnabled();
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
