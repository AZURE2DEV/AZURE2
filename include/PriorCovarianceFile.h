#ifndef PRIORCOVARIANCEFILE_H
#define PRIORCOVARIANCEFILE_H

#include <memory>
#include <string>
#include <vector>

#include "InferenceTypes.h"

class Config;

/*!
 * The prior covariance of a run (<config> priorCovariance=<path>): the engine's side of
 * azure2/inference/PriorCovariance.h, which holds the grammar and the arithmetic.  Here: the path relative to the
 * project, the file opened through the cache of input files, the messages to the run's stream.
 */

/// The key's path resolved against the .azr's directory (an absolute path is kept); empty without the key.
std::string ResolvePriorCovariancePath(const Config &configure);
/// Open and parse the file (PriorCovariance::Parse).  False with err set when it cannot be opened or parsed.
/// cached: read through InputFiles' cache (Config::cacheInputFiles, a pyazr session's rebuilds).
bool ReadPriorCovarianceFile(const std::string &path, PriorCovariance::File &out, std::string &err, bool cached = false);
/*!
 * Read and set up from the run's <config> key, printing the Notes and, on a refusal, the ERROR line to
 * configure.outStream.  paramNames / fixed: the run's whole parameter vector and its fixed flags; runBasis:
 * AZUREParams::BasisForConfig.  Returns 0 without the key (out reset), 1 with a prior set up (out), -1 on a refusal.
 */
int LoadPriorCovariance(const Config &configure, const std::vector<std::string> &paramNames, const std::vector<bool> &fixed,
                        int runBasis, std::shared_ptr<const PriorCovariance> &out);

#endif  // PRIORCOVARIANCEFILE_H
