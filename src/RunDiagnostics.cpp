// What the run's diagnostics count (RunDiagnostics.h).
#include "RunDiagnostics.h"

#include <cstdlib>
#include <iostream>

#include "Config.h"

namespace azure2 {

RunDiagnostics::~RunDiagnostics() { Finish(); }

void RunDiagnostics::Finish() {
  if (finished_.exchange(true)) return;
  const long h = coulombTable.lookups.load(), m = coulombTable.outside.load();
  if (h + m > 0)
    std::cerr << "Coulomb table: " << h << " lookups, " << m << " outside the table (exact routine used)" << std::endl;
}

}  // namespace azure2

void ExitRun(const Config &configure, int status) {
  configure.diagnostics->Finish();
  std::exit(status);
}
