#ifndef AZURE2_REFERENCETABLELEAF_H
#define AZURE2_REFERENCETABLELEAF_H

#include "azure2/experiment/ReferenceTable.h"

class Config;
class WorkingCopy;
class Unit;

namespace azure2 {

/*!
 * The value of the reference leaf of a ratio to a tabulated cross section (experiment::LeafPart::Reference; quantity
 * ratio_to_reference, data_kinds.md section 4.9) at a point, b: experiment::ReferenceAt at the point's present lab
 * energy (after an energy shift that re-energised it), folded as the point's observable is (FoldReferenceAtPoint:
 * the point's own sub-points and kernel, never a second quadrature), so the model is <sigma> / <sigma_ref>.  NaN where
 * an energy is outside the table (the set-up refuses rows and grid points there; a shift or a fold's sub-point can
 * still reach outside, and the point's model is then NaN).  `data` holds the point (its target effect).
 */
double ReferenceTableValue(const WorkingCopy &data, Unit *point, const experiment::ReferenceTable &table,
                           const Config &configure);

}  // namespace azure2

#endif
