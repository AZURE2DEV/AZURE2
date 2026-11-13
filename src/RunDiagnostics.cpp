// What the run's diagnostics count (RunDiagnostics.h).
#include "RunDiagnostics.h"

#include <iostream>

namespace azure2 {

RunDiagnostics::~RunDiagnostics() {
  const long h = coulombTable.lookups.load(), m = coulombTable.outside.load();
  if (h + m > 0)
    std::cerr << "Coulomb table: " << h << " lookups, " << m << " outside the table (exact routine used)" << std::endl;
}

}  // namespace azure2
