// CoulombTable::BuildFile, the build behind AZURE2 --build-coulomb-table and the GUI's
// Configure > Build Coulomb Table...: it writes a loadable table (creating the directory),
// counts its rows for the progress display, leaves nothing behind when cancelled, and
// ForgetMissing() lets a running process pick up a table built after it first looked.
#include "CoulombTable.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;
static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}

int main() {
  const fs::path dir = fs::temp_directory_path() / ("azure2_ctable_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::remove_all(dir);
  const std::string file = (dir / "sub" / "table.bin").string();

  CoulombTable::BuildSpec spec;  // a small table: seconds, not minutes
  spec.lMax = 1;
  spec.etaMin = 0.1, spec.etaMax = 10.0;
  spec.rhoMin = 0.1, spec.rhoMax = 10.0;
  spec.perDecade = 4.0;
  spec.threads = 2;

  // cancelled before it starts: false, "cancelled", no file and no .new
  std::atomic<bool> cancel{true};
  spec.cancel = &cancel;
  std::string error;
  check(!CoulombTable::BuildFile(file, spec, error) && error == "cancelled", "a cancelled build returns false");
  check(!fs::exists(file) && !fs::exists(file + ".new"), "a cancelled build leaves no file");

  // the running process looks for the table before it exists
  check(CoulombTable::Shared(file) == nullptr, "no table before the build");

  // the build: rows counted, directory created, file loadable, no .new left
  cancel = false;
  std::atomic<int> done{0};
  spec.rowsDone = &done;
  check(CoulombTable::BuildFile(file, spec, error), "the build succeeds");
  check(done.load() == CoulombTable::BuildRows(spec), "rowsDone reaches BuildRows");
  check(!fs::exists(file + ".new"), "no .new file left");
  CoulombTable t;
  check(t.Load(file, error) && t.LMax() == 1, "the file loads");

  // Shared() remembers the miss until ForgetMissing()
  check(CoulombTable::Shared(file) == nullptr, "Shared still remembers the miss");
  CoulombTable::ForgetMissing();
  const CoulombTable *shared = CoulombTable::Shared(file);
  check(shared != nullptr && shared->LMax() == 1, "Shared loads it after ForgetMissing");

  fs::remove_all(dir);
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
