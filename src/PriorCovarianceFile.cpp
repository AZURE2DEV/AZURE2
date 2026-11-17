// The prior covariance of a run: path, cached open and messages (include/PriorCovarianceFile.h); the grammar and the
// arithmetic are azure2/inference/PriorCovariance.h.
#include "PriorCovarianceFile.h"

#include <istream>

#include "Config.h"
#include "InputFiles.h"
#include "azure2/experiment/ProjectPaths.h"

std::string ResolvePriorCovariancePath(const Config &configure) {
  const std::string &p = configure.priorCovariance;
  if (p.empty() || azure2::experiment::IsAbsolutePath(p)) return p;
  const size_t slash = configure.configfile.find_last_of("/\\");
  if (slash == std::string::npos) return p;
  return configure.configfile.substr(0, slash + 1) + p;
}

bool ReadPriorCovarianceFile(const std::string &path, PriorCovariance::File &out, std::string &err, bool cached) {
  out = PriorCovariance::File();
  std::unique_ptr<std::istream> inStream = InputFiles::Open(path, cached);
  std::istream &in = *inStream;
  if (!in) {
    err = "cannot open the prior covariance file " + path;
    return false;
  }
  return PriorCovariance::Parse(in, path, out, err);
}

int LoadPriorCovariance(const Config &configure, const std::vector<std::string> &paramNames, const std::vector<bool> &fixed,
                        int runBasis, std::shared_ptr<const PriorCovariance> &out) {
  out.reset();
  if (configure.priorCovariance.empty()) return 0;
  auto prior = std::make_shared<PriorCovariance>();
  const std::string path = ResolvePriorCovariancePath(configure);
  PriorCovariance::File file;
  std::string err;
  std::vector<std::string> notes;
  if (!ReadPriorCovarianceFile(path, file, err, configure.cacheInputFiles) ||
      !prior->Setup(file, paramNames, fixed, runBasis, err, notes)) {
    configure.outStream << "**ERROR: priorCovariance: " << err << "." << std::endl;
    return -1;
  }
  for (const std::string &note : notes) configure.Note(note);
  configure.Note("Note: prior covariance " + prior->Path() + ": a multivariate Gaussian prior on " +
                 std::to_string(prior->Size()) + " parameter" + (prior->Size() == 1 ? "" : "s") +
                 " (Total-Prior-Covariance-Chi-Squared in chiSquared.out).");
  out = prior;
  return 1;
}
