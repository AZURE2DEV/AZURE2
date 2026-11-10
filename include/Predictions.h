#ifndef AZURE2_PREDICTIONS_H
#define AZURE2_PREDICTIONS_H

#include <iosfwd>
#include <string>
#include <vector>

#include "azure2/experiment/Prediction.h"

class CNuc;
class ESegment;

/*!
 * The expression trees of the engine's data sets and extrapolations (docs/upgrade/data_model.md sections 2.3 and 3,
 * step S10).  A segment's tree is built once, after EData::Fill or MakePoints has attached its components, resolved
 * its pairs and its sample (EData::BuildPredictions), from what the segment then holds; the model pass
 * (src/ModelPass.cpp) evaluates it at every point.
 */
namespace azure2 {

/// The setup of a segment's tree: its reaction as read, and its terms as the set-up made them (term 0 the segment's own
/// points, then its leaves in order), with the compound nucleus of each entrance pair.
experiment::PredictionSetup PredictionSetupOf(const ESegment &segment, CNuc &compound);

/// One line per data set and extrapolation of a project, "data <label>: <tree>" then "extrapolation <label>: <tree>",
/// in the order of the file (AZURE2 --print-predictions, pyazr.predictions); the label is the record's name in a
/// format-3 project, the key (1, 2, ...) in an old one.  The project is set up as a calculation
/// would set it up (data files read relative to the working directory); a set-up that fails gives one line
/// "data: set-up failed" or "extrapolation: set-up failed", and its messages go to `messages`.
std::vector<std::string> ProjectPredictions(const std::string &projectFile, std::ostream &messages);

}  // namespace azure2

#endif
