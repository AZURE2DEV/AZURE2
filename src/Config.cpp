#include "Config.h"

#include "azure2/convert/PotentialBlock.h"
#include "azure2/experiment/ProjectReader.h"
#include "InputFiles.h"
#include "azure2/io/Convert.h"
#include "azure2/io/Document.h"
#include <cstdlib>
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

Config::Config(std::ostream &stream) :
  outStream(stream) {
  Reset();
}

/*!
 * This function resets Config structure.
 */

void Config::Reset() {
  chiVariance = 1.0;
  screenCheckMask = 0;
  fileCheckMask = 0;
  paramMask = 0;
  paramMask |= (USE_AMATRIX | USE_BRUNE_FORMALISM | IGNORE_ZERO_WIDTHS | TRANSFORM_PARAMETERS | CALCULATE_WITH_DATA | USE_LONGWAVELENGTH_APPROX);
  stopFlag = false;
  outputdir = "";
  checkdir = "";
  nloptAlgorithm = 0;       // Default to SBPLX
  useHybridMethod = false;  // Default to disabled
  potentials = azure2::experiment::HybridPotentialSettings();
  refreezeCovariance = 0;
  maxCalls = 0;
  useAdaptiveGrid = true;   // Default to adaptive grid
  thm = ThmOptions();
  useCoulombTable = true;   // the universal Coulomb table where it covers a channel (--no-coulomb-table)
  replaceSubPointGrids = true;  // re-place the adaptive sub-point grids (--no-grid-replacement)
  coulombTablePath.clear();
  useECTable = true;            // external-capture energy tables (--no-ec-table)
  ecConvectionCurrent = true;   // convection current in exact-form EC (--no-ec-convection-current)
  ecTableNodesPerDecade = 50;   // ECTable node density; $AZURE2_EC_TABLE_NODES overrides
  includeAllJPi = false;    // Off unless the <config> block says so: older projects compute as they always did
  hardSphereMaxL = -1;
  abundanceSumSigma = 0.0;
  relativistic = false;      // Off unless the <config> block says so: non-relativistic kinematics as before
  relativisticReferencePair = 0;
  priorCovariance.clear();   // no prior parameter covariance unless the <config> block names a file
}

void Config::Note(const std::string &line) const {
  if (noteLog) {
    std::lock_guard<std::mutex> lock(noteLog->mutex);
    if (!noteLog->printed.insert(line).second) return;
  }
  outStream << line << std::endl;
}

/*!
 * This funciton reads the configuration file and parses various options.
 */

bool Config::ResolveMatrixRoute() {
  if (paramMask & USE_AMATRIX) return false;
  const bool ignore = (paramMask & IGNORE_ZERO_WIDTHS), park = (paramMask & USE_PARK_FORMALISM);
  if (!ignore && !park) return false;
  outStream << "WARNING: the R-matrix route (<config> first line false) does not implement ";
  if (ignore)
    outStream << "ignore-externals, which is on (CLI: the default, --no-ignore-externals turns it off; pyazr:"
              << " ignore_externals=False)";
  if (ignore && park) outStream << ", nor ";
  if (park) outStream << "Park's parametrization";
  outStream << ".  The A-matrix route will be used." << std::endl;
  paramMask |= USE_AMATRIX;
  return true;
}

std::unique_ptr<std::istream> Config::SectionStream(const std::string &name) const {
  std::unique_ptr<std::istream> file = InputFiles::Open(configfile, cacheInputFiles);
  if (!file || !*file) return nullptr;
  std::ostringstream bytes;
  bytes << file->rdbuf();
  azure2::io::Document doc = azure2::io::Document::Parse(bytes.str());
  // A format-1 file is read as its format-2 conversion (docs/upgrade/project.md section 4); ReadConfigFile notes it.
  azure2::io::ConvertToFormat2(doc);
  const azure2::io::Document::Section *s = doc.Find(name);
  if (!s) return nullptr;
  std::string text;
  const std::size_t end =
      (s->close == azure2::io::Document::Section::npos) ? doc.NumLines() : s->close + 1;
  for (std::size_t i = s->open; i < end; i++) {
    // The tags as written by the readers' tests ("<name>", "</name>"): trimmed (P2).
    bool closing = false;
    const std::string tag = azure2::io::TagOf(doc.Line(i), closing);
    text += (i == s->open || (i == s->close && !tag.empty())) ? (closing ? "</" : "<") + tag + ">" : doc.Line(i);
    text += "\n";
  }
  return std::unique_ptr<std::istream>(new std::istringstream(text));
}

/*!
 * A format-3 project, read once through the one reading path (docs/upgrade/project_format3.md section 7.3): records ->
 * Validate -> <levels> tables -> Schema::Build.  Every diagnostic is printed (file:line: record.field: message); an
 * error refuses the project.  Data files are read relative to the working directory, as the set-up opens them; one
 * that cannot be read is reported by the run that needs it (ProjectSetup.cpp), not here.
 */
int Config::ReadFormat3Project(const azure2::io::Document &doc) {
  azure2::experiment::Format3ReadOptions read;
  read.fileName = configfile;
  // Which run follows is chosen later (menu, pyazr call): a data file that cannot be read is reported by the set-up
  // that would open it (src/ProjectSetup.cpp: an error for a data run, a warning for an extrapolation)
  read.dataFiles = azure2::experiment::DataFiles::NotUsed;
  read.setUpChecks = false;  // the set-up makes them itself (CNuc::ParseExternalCapture), once, with its message
  auto project = std::make_shared<azure2::experiment::Format3Project>(azure2::experiment::ReadFormat3Project(doc, read));
  for (const azure2::experiment::Diagnostic &d : project->diagnostics)
    if (d.code != azure2::experiment::codes::kDataFileNotRead) outStream << (d.severity == azure2::experiment::Severity::Error ? "ERROR: " : "WARNING: ") << d.ToText() << std::endl;
  if (!project->Ok() || !project->levelsRead) {
    if (project->Ok()) outStream << "ERROR: " << configfile << ": the <levels> section is missing or malformed." << std::endl;
    outStream << "ERROR: " << configfile << " (format 3) is refused; AZURE2 --validate " << configfile
              << " lists the same diagnostics." << std::endl;
    return -1;
  }
  format3 = project;
  return 0;
}

int Config::ReadConfigFile() {
  std::string dummy;
  std::string temp;
  // The <config> section of the project (io::Document; a tag is a whole trimmed line).
  std::unique_ptr<std::istream> section = SectionStream("config");
  if (!section) return -1;
  // The file's structural problems (an unclosed or repeated section, crossing tags), once per read.
  {
    std::unique_ptr<std::istream> file = InputFiles::Open(configfile, cacheInputFiles);
    std::ostringstream bytes;
    if (file && *file) bytes << file->rdbuf();
    azure2::io::Document doc = azure2::io::Document::Parse(bytes.str());  // kept alive for the loops
    for (const std::string &problem : doc.Problems())
      outStream << "WARNING: " << configfile << ": " << problem << std::endl;
    format3.reset();
    if (azure2::io::FormatOf(doc) == 3 && ReadFormat3Project(doc) != 0) return -1;
    if (azure2::io::FormatOf(doc) == 1) {
      const std::vector<std::string> changes = azure2::io::ConvertToFormat2(doc);
      std::string list;
      for (const std::string &c : changes) list += (list.empty() ? "" : "; ") + c;
      Note("Note: " + configfile + " is in the old format (format 1); it is read as its format-2 conversion (" + list +
           "). AZURE2 --convert writes the converted file.");
    }
  }
  useAdaptiveGrid = true;  // the default; a <config> key (format 2) or <potential> (format 1, converted) sets it
  std::istream &in = *section;
  std::string line = "";
  getline(in, line);  // <config>
  in >> temp;
  getline(in, dummy);
  if (temp == "true")
    paramMask |= USE_AMATRIX;
  else
    paramMask &= ~USE_AMATRIX;
  getline(in, dummy);
  int poundSignPos = dummy.find_last_of('#');
  if (poundSignPos == std::string::npos)
    temp = dummy;
  else
    temp = dummy.substr(0, poundSignPos);
  int p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    outputdir = temp.substr(p1, (p2 - p1) + 1);
  } else
    outputdir = std::string();
  getline(in, dummy);
  poundSignPos = dummy.find_last_of('#');
  if (poundSignPos == std::string::npos)
    temp = dummy;
  else
    temp = dummy.substr(0, poundSignPos);
  p2 = temp.find_last_not_of(" \n\t\r");
  if (p2 != std::string::npos) {
    int p1 = temp.find_first_not_of(" \n\t\r");
    if (p1 == std::string::npos) p1 = 0;
    checkdir = temp.substr(p1, (p2 - p1) + 1);
  } else
    checkdir = std::string();
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_COMPOUND_NUCLEUS;
  else if (temp == "file")
    fileCheckMask |= CHECK_COMPOUND_NUCLEUS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  else if (temp == "file")
    fileCheckMask |= CHECK_BOUNDARY_CONDITIONS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_DATA;
  else if (temp == "file")
    fileCheckMask |= CHECK_DATA;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_ENERGY_DEP;
  else if (temp == "file")
    fileCheckMask |= CHECK_ENERGY_DEP;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_LEGENDRE;
  else if (temp == "file")
    fileCheckMask |= CHECK_LEGENDRE;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_COUL_AMPLITUDES;
  else if (temp == "file")
    fileCheckMask |= CHECK_COUL_AMPLITUDES;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_PATHWAYS;
  else if (temp == "file")
    fileCheckMask |= CHECK_PATHWAYS;
  in >> temp;
  getline(in, dummy);
  if (temp == "screen")
    screenCheckMask |= CHECK_ANGULAR_DISTS;
  else if (temp == "file")
    fileCheckMask |= CHECK_ANGULAR_DISTS;
  // Optional key=value lines after the fixed ones. A reader that does not know a key skips it, so a
  // project carrying them still opens in older builds.
  line = "";
  while (line != "</config>" && !in.eof()) {
    getline(in, line);
    std::string entry = line;
    size_t hash = entry.find('#');
    if (hash != std::string::npos) entry = entry.substr(0, hash);
    size_t eq = entry.find('=');
    if (eq == std::string::npos) continue;
    std::string key = entry.substr(0, eq), value = entry.substr(eq + 1);
    key.erase(0, key.find_first_not_of(" \t\r"));
    key.erase(key.find_last_not_of(" \t\r") + 1);
    value.erase(0, value.find_first_not_of(" \t\r"));
    value.erase(value.find_last_not_of(" \t\r") + 1);
    if (key == "includeAllJPi")
      includeAllJPi = (value == "1" || value == "true");
    else if (key == "maxOrbitalMomentum")
      hardSphereMaxL = atoi(value.c_str());
    else if (key == "abundanceSumSigma")
      abundanceSumSigma = atof(value.c_str());
    else if (key == "relativistic")
      relativistic = (value == "1" || value == "true");
    else if (key == "relativisticReferencePair")
      relativisticReferencePair = atoi(value.c_str());
    else if (key == "priorCovariance")
      priorCovariance = value;
    else if (key == "useAdaptiveGrid")  // format 2 (format 1: in <potential>)
      useAdaptiveGrid = (value == "1" || value == "true");
    else if (key == "format")  // the file format (docs/upgrade/project.md section 4); read by io::FormatOf
      continue;
  }
  if (line != "</config>") return -1;
  // A malformed optional block has already been reported: -2, so the caller
  // does not add a misleading "could not open" on top.
  if (this->ReadPotentialBlock() != 0) return -2;
  return this->ReadThmBlock() != 0 ? -2 : 0;
}

/*!
 * The hybrid potentials of the project, per pair.  A format-3 project states them as `potential` records of
 * <potentials>, which Schema::Build has made into settings (one per pair, no default); the defaults of the table it
 * used are noted.  An older project states them in the <potential> block, read with the old block's one parser
 * (convert/ PotentialBlock.h, which the converter to format 3 uses too) and its warnings written here; without the
 * block no pair has a potential.  useAdaptiveGrid in the block (format 1) is a property of the target-effect
 * integration and stays global.
 *
 * Returns 0 on success (including "no block present") and -1 if the block is present but not closed.
 */
int Config::ReadPotentialBlock() {
  potentials = azure2::experiment::HybridPotentialSettings();
  useHybridMethod = false;
  if (format3) {
    potentials = format3->built.potentials;
    for (const std::string &note : format3->built.potentialNotes) Note("Note: " + configfile + ": " + note);
    useHybridMethod = potentials.AnyEnabled();
    return 0;
  }
  std::unique_ptr<std::istream> section = SectionStream("potential");
  if (!section) return 0;  // optional block, absent: no hybrid potential (also after a project that had one)
  const azure2::convert::PotentialBlock block = azure2::convert::ReadPotentialBlock(*section);
  for (const std::string &warning : block.warnings) outStream << warning << std::endl;
  if (block.hasAdaptiveGrid) useAdaptiveGrid = block.adaptiveGrid;
  potentials = block.settings;
  if (!block.closed) return -1;  // unterminated block
  // On when any pair has the model; CNuc::Fill gives each pair its potential.
  useHybridMethod = potentials.AnyEnabled();
  return 0;
}

/*!
 * If stat() is enabled, this function checks for the output and checks
 * directories at runtime.
 */

#ifndef NO_STAT
int Config::CheckForInputFiles() {
  struct stat buffer;
  if (stat(outputdir.c_str(), &buffer) != 0) {
    outStream << "Could not find output directory: " << outputdir << ". Check that it exists." << std::endl;
    return -1;
  }
  if (stat(checkdir.c_str(), &buffer) != 0) {
    outStream << "Could not find checks directory: " << checkdir << ". Check that it exists." << std::endl;
    return -1;
  }
  return 0;
}
#endif
