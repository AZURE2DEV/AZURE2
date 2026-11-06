// InputFiles::Open(path, cached = true) serves a file from memory while its stamp (size, modification time, and
// on POSIX inode/device/ctime) is unchanged.  A stamp cannot tell two versions apart when the second is written
// within the clock's resolution: on Windows the stamp is size + modification time only, and two same-size writes
// in quick succession share it (tests/pyazr/natural_target_test failed on the Windows CI from at least 2026-11-04:
// a rewritten data.dat was served with the previous contents, chi2 x96).  Since 2026-11-05 a file whose
// modification time lies within InputFiles' racy window of the moment it was read is re-read and compared rather
// than trusted (git's "racily clean" rule).
//
// (a) a same-size rewrite with the modification time put back to the old value, within one second (so ctime, whole
//     seconds on POSIX, is unchanged too): the second open returns the NEW contents;
// (b) a file last modified long before it is read is still served from memory (no second disk read).
#include "InputFiles.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>

static int failures = 0;
static void check(bool ok, const char *what) {
  std::cout << (ok ? "  ok    " : "  FAIL  ") << what << std::endl;
  if (!ok) failures++;
}

static std::string readAll(const std::string &p) {
  auto s = InputFiles::Open(p, true);
  std::ostringstream o;
  o << s->rdbuf();
  return o.str();
}

int main() {
  namespace fs = std::filesystem;
  const std::string p = (fs::temp_directory_path() / "azure2_input_files_cache_test.dat").string();
  InputFiles::Clear();

  // (a) start just after a second boundary so the whole sequence stays within one second
  auto now = std::chrono::system_clock::now().time_since_epoch();
  std::this_thread::sleep_for(std::chrono::seconds(1) - (now % std::chrono::seconds(1)) + std::chrono::milliseconds(5));
  { std::ofstream f(p); f << "1.000000 0 1.000000000000e+00 1.000000000000e-03\n"; }
  const auto t = fs::last_write_time(p);
  const std::string a = readAll(p);
  { std::ofstream f(p); f << "1.000000 0 1.000000000000e+00 1.000000000000e-02\n"; }  // same size, other error
  fs::last_write_time(p, t);                                                              // the same stamp
  const std::string b = readAll(p);
  check(a.find("e-03") != std::string::npos, "(a) first open reads the first version");
  check(b.find("e-02") != std::string::npos, "(a) a same-size rewrite with the old stamp is read, not served stale");

  // (b) an old file: modification time an hour back, read twice, the second from memory
  InputFiles::Clear();
  { std::ofstream f(p); f << "old file\n"; }
  fs::last_write_time(p, fs::file_time_type::clock::now() - std::chrono::hours(1));
  const long reads0 = InputFiles::GetStats().reads;
  const std::string c1 = readAll(p), c2 = readAll(p);
  const InputFiles::Stats st = InputFiles::GetStats();
  check(c1 == "old file\n" && c2 == c1, "(b) contents of an old file");
  check(st.reads - reads0 == 1, "(b) an old, unchanged file is read once and then served from memory");

  std::remove(p.c_str());
  if (failures) {
    std::cout << "FAIL: " << failures << " check(s)" << std::endl;
    return 1;
  }
  std::cout << "PASS: input file cache" << std::endl;
  return 0;
}
