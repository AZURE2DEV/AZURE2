#ifndef AZURE2_SEGMENTSTATE_H
#define AZURE2_SEGMENTSTATE_H

#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "ChannelTable.h"
#include "Equation.h"
#include "azure2/experiment/Frames.h"

/*!
 * What a working copy holds of one term of a data set (docs/upgrade/data_structure.md section 1.4): where its points
 * lie in the copy's unit vector, its rows in the c.m. frame at the present energy shift, and the channel tables of
 * its last shift states.  The term's settings and its frame at load are the data set's (Segment, immutable).
 */
struct TermState {
  struct ShiftTable {
    std::vector<double> state;  ///< energy shift (MeV), sqrt shift (MeV^1/2), shift-function values
    std::shared_ptr<const ChannelTable> table;
  };
  int first = 0;                                ///< the points: units [first, first + rows) of the working copy
  azure2::experiment::MeasurementFrame frame;   ///< the rows in the c.m. frame at the present shift
  std::vector<ShiftTable> shiftTables;          ///< most recent first (Segment::kShiftTables)
  std::vector<std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>>> shiftRows;  ///< per point, in a refill
};

/*!
 * What a working copy holds of one data set (a segment of the project): the values the hand-over of an evaluation
 * writes (normalization, energy shifts, shift function, sample thickness) and its terms' states.  Plain values: a
 * copy of the working copy copies them.  The settings and the start values are the data set's (Segment).
 */
struct SegmentState {
  double norm = 1.0;                  ///< normalization applied to the data
  double energyShift = 0.0;           ///< constant lab-energy shift (MeV)
  double lastEnergyShift = std::numeric_limits<double>::quiet_NaN();      ///< the shift the points carry (NaN: none yet)
  double energyShiftSqrt = 0.0;       ///< coefficient b of b sqrt(E/MeV) (MeV^1/2)
  double lastEnergyShiftSqrt = std::numeric_limits<double>::quiet_NaN();  ///< likewise
  /// The applied shift function with its present parameter values (Segment::HasShiftFn), and those values as the
  /// points carry them.
  Equation shiftFn;
  std::vector<double> shiftFnValues, lastShiftFnValues;
  double thickness = 0.0;             ///< sample thickness of a transmission or yield segment (atoms/barn)
  std::vector<TermState> terms;       ///< term 0 the segment's own, k its k-th leaf
};

/*!
 * The present values of one target effect's fit parameters (docs/upgrade/data_structure.md section 1.4, "the
 * apparatus values that depend on the parameters"); the effect itself (lines, kernels, priors) is the data set's.
 */
struct TargetEffectValues {
  double screeningUe = 0.0;         ///< electron-screening potential U_e (eV)
  double density = 0.0;             ///< areal density of the target (atoms/cm^2)
  double dopplerTemperature = 0.0;  ///< free-gas temperature (K)
  double stoppingScale = 1.0;       ///< stopping scale s (TargetEffect::StoppingAt), dimensionless
};

#endif
