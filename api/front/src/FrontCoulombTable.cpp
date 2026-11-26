// The Coulomb table of this computer for a front end (Engine.h, docs/upgrade/qt_gui.md section 19.12): the state a run
// finds and the command line's build (CoulombTable::BuildFile, which AZURE2 --build-coulomb-table calls).
#include "azure2/front/Engine.h"
#include "azure2/support/ClassicNumeric.h"

#include <algorithm>
#include <filesystem>
#include <system_error>

#include "CoulombTable.h"

namespace azure2 {
namespace front {

namespace {

/// A front end's build runs beside the window and its runs: never more than 4 threads.
constexpr int kMaxFrontThreads = 4;

CoulombTable::BuildSpec SpecOf(const CoulombTableBuild &b) {
  CoulombTable::BuildSpec spec;
  spec.lMax = b.lMax;
  spec.etaMin = b.etaMin;
  spec.etaMax = b.etaMax;
  spec.rhoMin = b.rhoMin;
  spec.rhoMax = b.rhoMax;
  spec.perDecade = b.perDecade;
  spec.etaRhoMax = b.etaRhoMax;
  spec.threads = std::max(1, std::min(b.threads > 0 ? b.threads : CoulombTable::DefaultBuildThreads(), kMaxFrontThreads));
  return spec;
}

}  // namespace

CoulombTableState CoulombTableStatus() {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  CoulombTableState s;
  s.path = CoulombTable::DefaultPath();
  std::error_code ec;
  s.exists = std::filesystem::is_regular_file(s.path, ec);
  if (!s.exists) return s;
  const auto size = std::filesystem::file_size(s.path, ec);
  s.bytes = ec ? 0 : static_cast<long long>(size);
  // the loader the runs call: a usable table is read once and kept for them.  A file it could not read (or that was
  // not there when it looked) is not tried again until ForgetMissing, so its reason is asked of a load of our own; a
  // file that loads now (put there since) is handed to the loader again.
  const CoulombTable *t = CoulombTable::Shared("", s.path);
  if (!t) {
    CoulombTable probe;
    if (!probe.Load(s.path, s.problem)) {
      if (s.problem.empty()) s.problem = "cannot be read";
      return s;
    }
    CoulombTable::ForgetMissing();
    t = CoulombTable::Shared("", s.path);
    if (!t) {
      s.problem = "cannot be read";
      return s;
    }
  }
  s.usable = true;
  s.lMax = t->LMax();
  s.nodesPerDecade = t->EtaAxis().dx > 0 ? 1.0 / t->EtaAxis().dx : 0.0;
  s.etaRhoMax = t->EtaRhoMax();
  return s;
}

CoulombTableBuild DefaultCoulombTableBuild() {
  const CoulombTable::BuildSpec spec;
  CoulombTableBuild b;
  b.lMax = spec.lMax;
  b.etaMin = spec.etaMin;
  b.etaMax = spec.etaMax;
  b.rhoMin = spec.rhoMin;
  b.rhoMax = spec.rhoMax;
  b.perDecade = spec.perDecade;
  b.etaRhoMax = spec.etaRhoMax;
  b.threads = CoulombTable::DefaultBuildThreads();
  return b;
}

int CoulombTableBuildRows(const CoulombTableBuild &b) { return CoulombTable::BuildRows(SpecOf(b)); }

bool BuildCoulombTable(const CoulombTableBuild &b, std::atomic<int> *rowsDone, const std::atomic<bool> *cancel,
                       std::string &error) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  CoulombTable::BuildSpec spec = SpecOf(b);
  spec.rowsDone = rowsDone;
  spec.cancel = cancel;
  const std::string path = b.path.empty() ? CoulombTable::DefaultPath() : b.path;
  if (!CoulombTable::BuildFile(path, spec, error)) return false;
  CoulombTable::ForgetMissing();
  return true;
}

}  // namespace front
}  // namespace azure2
