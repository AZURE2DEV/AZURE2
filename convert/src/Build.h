#ifndef AZURE2_CONVERT_BUILD_H
#define AZURE2_CONVERT_BUILD_H

// The format-3 records of an old project (Build.cpp): section 8.1 of docs/upgrade/project_format3.md, construct by
// construct.  Internal to the converter.
#include <map>
#include <string>
#include <vector>

#include "OldProject.h"
#include "azure2/convert/Convert.h"
#include "azure2/io/Record.h"

namespace azure2 {
namespace convert {

/// The records of the four format-3 sections, the new texts of the kept sections, and what was noted and asked.
struct Built {
  io::RecordSet data, targets, extrapolations, modelData, potentials;
  /// Replacement bodies (lines between the tags) of sections that stay: thm, parameterSettings; "" key absent = keep.
  std::map<std::string, std::vector<std::string>> bodies;
  bool hasThm = false;
  /// <config> edits: key -> new value ("" removes the line).
  std::map<std::string, std::string> configEdits;
};

/// Builds the records.  Notes, questions, counts, the parameter map and side files go into `result`.  Returns "" or
/// why the project cannot be converted.
std::string BuildFormat3(const io::Document &doc, const OldProject &old, const ConvertOptions &options, Built &out,
                         ConvertResult &result);

}  // namespace convert
}  // namespace azure2

#endif
