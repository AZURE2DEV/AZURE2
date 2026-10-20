#include "AZUREParams.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "Config.h"
#include "CNuc.h"

namespace {
// segment_<key>_energy_shift, ..._energy_shift_sqrt, ..._energy_shift_fn<j>
bool IsEnergyShiftName(const std::string &name) {
  return name.compare(0, 8, "segment_") == 0 && name.find("_energy_shift") != std::string::npos;
}
}  // namespace

/*!
 * This function returns the MnUserParameters object used by Minuit to store
 * the fit parameters.
 */

ROOT::Minuit2::MnUserParameters &AZUREParams::GetMinuitParams() {
  return params_;
}

/*!
 * This function reads the user specified parameters from a given file.
 * These parameters are formal R-matrix parameters, and overwrite any
 * initial parameters determined from the nuclear input file.
 */

void AZUREParams::ReadUserParameters(const Config &configure) {
  ApplyUserParameterFile(configure.paramfile, configure.outStream);
}

void AZUREParams::ReadUserParameters(const std::string &filename) {
  ApplyUserParameterFile(filename, std::cout);
}

/*!
 * Reads a parameter file -- one "name value error [fixed]" line per parameter --
 * and applies it to the parameters of the same names.
 *
 * The lines are parsed one at a time and the numbers with strtod, which accepts
 * "nan" and "inf".  The file used to be read with operator>>, which does not: a
 * failed extraction left the stream in fail() but never eof(), so the loop never
 * ended and grew the vectors until std::bad_alloc.  A fit whose covariance is
 * invalid writes "nan" in the error column, so its own .sav could not be read
 * back (2026-10-19, A=4 session).  A non-finite error now leaves that
 * parameter's step size as it was; a line without a readable value is skipped
 * with a warning.
 */
void AZUREParams::ApplyUserParameterFile(const std::string &filename, std::ostream &out) {
  std::vector<std::string> names;
  vector_r values;
  vector_r errors;
  std::vector<bool> fixed;

  std::ifstream in(filename.c_str());
  if (!in) {
    out << "Could not read user parameter file." << std::endl;
    return;
  }
  auto number = [](const std::string &token, double &v) {
    if (token.empty()) return false;
    char *end = nullptr;
    v = std::strtod(token.c_str(), &end);
    return end && *end == '\0';
  };
  std::string line;
  int lineNumber = 0, skipped = 0;
  while (std::getline(in, line)) {
    lineNumber++;
    std::istringstream fields(line);
    std::string name, valueToken, errorToken, rest, word;
    if (!(fields >> name)) continue;  // blank line
    fields >> valueToken >> errorToken;
    while (fields >> word) rest += word;
    double value, error = 0.0;
    if (!number(valueToken, value)) {
      if (skipped++ < 5)
        out << "WARNING: " << filename << " line " << lineNumber << ": no number for " << name << "; line skipped."
            << std::endl;
      continue;
    }
    if (name == "#parametrization" || name == "parametrization") {
      basisTag_ = (int)value;
      continue;
    }
    if (!number(errorToken, error)) error = std::nan("");
    names.push_back(name);
    values.push_back(value);
    errors.push_back(error);
    fixed.push_back(rest == "fixed");
  }

  for (int i = 0; i < GetMinuitParams().Params().size(); i++) {
    for (int ii = 0; ii < names.size(); ii++) {
      if (GetMinuitParams().GetName(i) == names[ii]) {
        // A level parameter fixed at zero is freed by a nonzero .sav value, and a
        // free one fixed by a zero.  An energy shift is not: it follows the vary
        // flag of its .azr segment line, because its penalty does (a freed shift
        // would float with no penalty).  The .sav still sets its value.
        const bool isShift = IsEnergyShiftName(GetMinuitParams().GetName(i));
        if (!isShift && GetMinuitParams().Value(i) == 0.0 && values[ii] != 0.0 &&
            GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Release(i);
        if (!isShift && GetMinuitParams().Value(i) != 0.0 && values[ii] == 0.0 &&
            !GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Fix(i);
        GetMinuitParams().SetValue(i, values[ii]);
        if (std::isfinite(errors[ii])) GetMinuitParams().SetError(i, errors[ii]);
        if (fixed[ii] && !GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Fix(i);
      }
    }
  }
}

/*!
 * This function writes the formal R-matrix parameters to a file.
 */

void AZUREParams::WriteUserParameters(const Config &configure, bool fitParameters) {
  char filename[256];
  if (fitParameters)
    snprintf(filename, sizeof(filename), "%sparam.sav", configure.outputdir.c_str());
  else
    snprintf(filename, sizeof(filename), "%sparam.par", configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if (out) {
    out.precision(7);
    // First line: which amplitudes the file holds (standard, Brune or Park
    // reduced widths differ by level-dependent factors).  Written as a
    // '#'-prefixed pseudo-parameter: numpy.loadtxt skips it as a comment, so
    // scripts that read the file by row position are unaffected; readers that
    // key on the first column see an unknown name; older AZURE2 binaries
    // parse it as a parameter no model has and ignore it.
    out << std::setw(20) << "#parametrization"
        << std::scientific << std::setw(20) << (double)BasisForConfig(configure)
        << std::scientific << std::setw(20) << 0.0 << std::endl;
    for (int i = 0; i < GetMinuitParams().Params().size(); i++) {
      out << std::setw(20) << GetMinuitParams().GetName(i)
          << std::scientific << std::setw(20) << GetMinuitParams().Value(i)
          << std::scientific << std::setw(20) << GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure.outStream << "Could not save param.par file." << std::endl;
}

int AZUREParams::BasisForConfig(const Config &configure) {
  if (configure.paramMask & Config::USE_PARK_FORMALISM) return kBasisPark;
  if (configure.paramMask & Config::USE_BRUNE_FORMALISM) return kBasisBrune;
  return kBasisStandard;
}

const char *AZUREParams::BasisName(int basis) {
  switch (basis) {
    case kBasisStandard: return "standard (constant boundary condition)";
    case kBasisBrune: return "Brune";
    case kBasisPark: return "Park (observed)";
    default: return "untagged";
  }
}

/*!
 * This function writes the parameter errors to a file if Minos has been invoked.
 */

void AZUREParams::WriteParameterErrors(const std::vector<std::pair<double, double>> &errors, const Config &configure) {
  char filename[256];
  snprintf(filename, sizeof(filename), "%sparam.errors", configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if (out) {
    out.precision(7);
    for (int i = 0; i < GetMinuitParams().Params().size(); i++) {
      out << std::setw(20) << GetMinuitParams().GetName(i)
          << std::scientific << std::setw(20) << GetMinuitParams().Value(i)
          << std::scientific << std::setw(20) << fabs(errors[i].first)
          << std::scientific << std::setw(20) << fabs(errors[i].second)
          << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure.outStream << "Could not save param.errors file." << std::endl;
}

bool AZUREParams::ReconcileBasis(CNuc *compound, const Config &configure) {
  const int mode = BasisForConfig(configure);
  int tag = basisTag_;
  if (tag == kBasisUnknown) {
    // Files written before the tag existed hold the amplitudes of the mode
    // that wrote them, which could not have been Park's.
    if (mode != kBasisPark) return true;
    configure.outStream << "The parameter file carries no parametrization tag; taking its amplitudes as Brune's."
                        << std::endl;
    tag = kBasisBrune;
  }
  if (tag == mode) return true;
  if (compound->ConvertAmplitudeBasis(params_, tag, mode, configure)) {
    configure.outStream << "Converted the parameter file from " << BasisName(tag) << " to "
                        << BasisName(mode) << " amplitudes." << std::endl;
    basisTag_ = mode;
    return true;
  }
  configure.outStream << "**ERROR: the parameter file holds " << BasisName(tag) << " amplitudes and this run uses "
                      << BasisName(mode) << " ones; they cannot be converted level by level." << std::endl;
  return false;
}
