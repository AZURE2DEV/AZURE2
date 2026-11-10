#ifndef AZURE2_EXPERIMENT_DESCRIBE_H
#define AZURE2_EXPERIMENT_DESCRIBE_H

#include <string>

#include "azure2/experiment/ProjectReader.h"

namespace azure2 {
namespace experiment {

/*!
 * `AZURE2 --describe FILE [--json]` and pyazr.describe (docs/upgrade/project_format3.md section 6.3): what
 * Schema::Build made of a format-3 project, the place an agent looks to see how its file was understood.  For every
 * data set and extrapolation: the quantity, the observable the engine computes, the reaction (its tree) with the
 * <levels> pair keys each term resolved to, the frames of its windows, the data file as the engine opens it (the
 * resolved path: data files are read relative to the working directory), the target, beam and detector effects
 * linked to it, and its fit parameters with their names and kinds (fit, fixed, shared, profiled); then the model
 * data with their levels and pairs.  Library azure2_experiment_reader.  `project` must have been built (Ok()).
 */
std::string DescribeText(const Format3Project &project);
std::string DescribeJson(const Format3Project &project);

}  // namespace experiment
}  // namespace azure2

#endif
