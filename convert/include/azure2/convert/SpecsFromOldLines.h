#ifndef AZURE2_CONVERT_SPECSFROMOLDLINES_H
#define AZURE2_CONVERT_SPECSFROMOLDLINES_H

// TEMPORARY: the old project sections as the engine's set-up specs (experiment/'s SetupSpecs.h).  The engine reads
// a format-1/2 project through these functions until step 3d of project_format3.md deletes them with the old
// sections; a format-3 project gives the same specs through Schema::Build.  Each function reads one section's text
// with the old grammar's one parser (SegLine.h, ExtrapLine.h, TargetIntLine.h, AbundanceLine.h, PriorCentreRow.h)
// and makes the specs and the messages the old engine printed for it, in its order and with its words.

#include <istream>
#include <string>

#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace convert {

/// The <segmentsData> section (its whole text, tags included).  `relativistic`: the project's relativistic flag (the
/// old engine refused a THM line with it while reading the line).
experiment::SpecList<experiment::DatasetSpec> DatasetSpecsFromOldLines(std::istream &section, bool relativistic);

/// The <segmentsTest> section.
experiment::SpecList<experiment::ExtrapolationSpec> ExtrapolationSpecsFromOldLines(std::istream &section,
                                                                                   bool relativistic);

/// The <targetInt> section.  Every line is an item; `readFailed` marks a line the old reader could not finish.
experiment::SpecList<experiment::TargetEffectSpec> TargetEffectSpecsFromOldLines(std::istream &section);

/// One <targetInt> line (the old engine's TargetEffect read one line from a stream).  `readFailed` as above.
experiment::TargetEffectSpec TargetEffectSpecFromOldLine(std::istream &line, bool &readFailed);

/// The <abundances> section; `section` null when the project has none.
experiment::AbundanceSpecs AbundanceSpecsFromOldLines(std::istream *section);

/// The prior-centre rows of <parameterSettings> (`section` null when the project has none), written into the data
/// sets they name: a row names a parameter by its old name, segment_<N>_norm or segment_<N>_energy_shift, N the key
/// of a <segmentsData> line (inactive lines counted); the last row for a parameter wins.  Returns the message of
/// the first row it refuses ("" = none): a name that matches no line, a centre that is not a finite number, a
/// normalization centre that is not positive.  `datasets` holds every line of the section.
std::string PriorCentresFromOldRows(std::istream *section, experiment::SpecList<experiment::DatasetSpec> &datasets);

/// The old name of a data set's normalization (isNorm) or energy-shift parameter: the names prior_centre rows of
/// an old project use (the adaptor writes the same names into DatasetSpec::names, which the engine gives its parameters).
std::string OldSegmentParameterName(int key, bool isNorm);

}  // namespace convert
}  // namespace azure2

#endif
