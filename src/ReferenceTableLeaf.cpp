// The reference leaf of a ratio to a tabulated cross section (ReferenceTableLeaf.h).
#include "ReferenceTableLeaf.h"

#include "CoulombReference.h"
#include "Unit.h"

namespace azure2 {

double ReferenceTableValue(const WorkingCopy &data, Unit *point, const experiment::ReferenceTable &table,
                           const Config &configure) {
  return FoldReferenceAtPoint(data, point, configure,
                              [&](const Site *p) { return experiment::ReferenceAt(table, p->GetLabEnergy()); });
}

}  // namespace azure2
