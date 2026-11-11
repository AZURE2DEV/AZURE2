#include "InputFiles.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <sstream>
#include <system_error>
#ifndef _WIN32
#include <sys/stat.h>
#endif

namespace InputFiles {

namespace {

// What identifies one version of a file on disk.
struct Stamp {
  std::uintmax_t size = 0;
  long long mtime = 0;  // file_clock ticks (ns on Linux)
  unsigned long long ino = 0, dev = 0;
  long long ctime = 0;
  bool operator==(const Stamp &o) const {
    return size == o.size && mtime == o.mtime && ino == o.ino && dev == o.dev && ctime == o.ctime;
  }
};

bool StampOf(const std::string &path, Stamp &s) {
  std::error_code ec;
  const std::filesystem::path p(path);
  if (!std::filesystem::is_regular_file(p, ec) || ec) return false;
  s.size = std::filesystem::file_size(p, ec);
  if (ec) return false;
  s.mtime = static_cast<long long>(std::filesystem::last_write_time(p, ec).time_since_epoch().count());
  if (ec) return false;
#ifndef _WIN32
  struct stat st;
  if (::stat(path.c_str(), &st) != 0) return false;
  s.ino = static_cast<unsigned long long>(st.st_ino);
  s.dev = static_cast<unsigned long long>(st.st_dev);
  s.ctime = static_cast<long long>(st.st_ctime);
#endif
  return true;
}

struct Entry {
  Stamp stamp;
  std::shared_ptr<const std::string> text;
  long long readAt = 0;  // file_clock ticks when the text was read
};

// A file modified within this window of the moment it was read is "racily clean" (git's term): a later write in the
// same clock tick can leave every part of the stamp unchanged -- on Windows the stamp is size + modification time
// only, and two same-size writes in quick succession share it -- so its stamp is not trusted and the file is read
// again.  3 s covers file systems with 2 s timestamps (FAT).  Before 2026-11-05 a matching stamp was always trusted:
// pyazr served a rewritten data file with its previous contents (Windows CI, pyazr_natural_target).
long long RacyWindowTicks() {
  return std::chrono::duration_cast<std::filesystem::file_time_type::duration>(std::chrono::seconds(3)).count();
}
long long FileClockNow() {
  return static_cast<long long>(std::filesystem::file_time_type::clock::now().time_since_epoch().count());
}
bool Trusted(const Entry &e) { return e.stamp.mtime + RacyWindowTicks() < e.readAt; }

// Keyed by the absolute path (relative paths resolved against the working directory at the time of
// the open, so sessions in different directories do not evict each other's "data/x.dat"); the stamp
// (inode, device) also tells two versions of one name apart.
std::mutex g_mutex;
std::map<std::string, Entry> g_files;
std::size_t g_bytes = 0;
Stats g_stats;
const std::size_t kMaxBytes = std::size_t(512) << 20;

}  // namespace

std::unique_ptr<std::istream> Open(const std::string &path, bool cached) {
  if (!cached) return std::unique_ptr<std::istream>(new std::ifstream(path.c_str()));
  std::error_code aec;
  const std::string key = std::filesystem::absolute(std::filesystem::path(path), aec).lexically_normal().string();
  if (aec) return std::unique_ptr<std::istream>(new std::ifstream(path.c_str()));
  Stamp now;
  if (!StampOf(path, now)) {
    // Missing or not a regular file: the caller sees the same failed stream as before.
    std::lock_guard<std::mutex> lock(g_mutex);
    auto it = g_files.find(key);
    if (it != g_files.end()) {
      g_bytes -= it->second.text->size();
      g_files.erase(it);
    }
    return std::unique_ptr<std::istream>(new std::ifstream(path.c_str()));
  }
  {
    std::lock_guard<std::mutex> lock(g_mutex);
    auto it = g_files.find(key);
    if (it != g_files.end() && it->second.stamp == now && Trusted(it->second)) {
      g_stats.served++;
      return std::unique_ptr<std::istream>(new std::istringstream(*it->second.text));
    }
  }
  // Text mode, as the plain ifstream reads it (line endings on Windows).
  const long long readAt = FileClockNow();
  std::ifstream in(path.c_str());
  if (!in) return std::unique_ptr<std::istream>(new std::ifstream(path.c_str()));
  std::ostringstream buf;
  buf << in.rdbuf();
  auto text = std::make_shared<const std::string>(buf.str());
  // A file that changed while it was read is not kept: the next open reads it again.
  Stamp after;
  const bool stable = StampOf(path, after) && after == now;
  {
    std::lock_guard<std::mutex> lock(g_mutex);
    g_stats.reads++;
    auto it = g_files.find(key);
    if (it != g_files.end()) {
      g_bytes -= it->second.text->size();
      g_files.erase(it);
    }
    if (stable && text->size() <= kMaxBytes) {
      if (g_bytes + text->size() > kMaxBytes) {
        g_files.clear();
        g_bytes = 0;
      }
      g_files[key] = Entry{now, text, readAt};
      g_bytes += text->size();
    }
  }
  return std::unique_ptr<std::istream>(new std::istringstream(*text));
}

Stats GetStats() {
  std::lock_guard<std::mutex> lock(g_mutex);
  Stats s = g_stats;
  s.files = static_cast<long>(g_files.size());
  return s;
}

void Clear() {
  std::lock_guard<std::mutex> lock(g_mutex);
  g_files.clear();
  g_bytes = 0;
}

}  // namespace InputFiles
