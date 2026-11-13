// The standard output of one pyazr session: the stream its Config writes to (Config::outStream).
//
// The engine prints its set-up and evaluation messages to Config::outStream; the CLI gives it std::cout.  A session
// gives it a NoticeBuffer instead, which passes every line on to std::cout as before, except the lines of the
// table kHeldNotices (pybind_notices.cpp): those it keeps, and pyazr turns them into warnings
// (pyazr.EngineWarning) that name the file concerned and the consequence (Session::take_notices).
#ifndef AZURE2_PYBIND_NOTICES_H
#define AZURE2_PYBIND_NOTICES_H

#include <mutex>
#include <ostream>
#include <streambuf>
#include <string>
#include <vector>

namespace azure2_binding {

/// One notice of the engine that a session turns into a warning: the line as the engine prints it, a key that says
/// which file it concerns (the session resolves it to a path), and the consequence in words.
struct HeldNotice {
  const char *text;
  const char *file;
  const char *consequence;
};

/// The notices a session turns into warnings (the line must be the engine's whole line).
const std::vector<HeldNotice> &HeldNotices();

/// A line buffer in front of `out`: a whole line that is one of HeldNotices() is kept, every other line is written to
/// `out` when its newline arrives (or at a flush).  Thread-safe: the engine may write from its OpenMP threads.
class NoticeBuffer : public std::streambuf {
 public:
  explicit NoticeBuffer(std::ostream &out) : out_(out) {}
  ~NoticeBuffer() override;
  NoticeBuffer(const NoticeBuffer &) = delete;
  NoticeBuffer &operator=(const NoticeBuffer &) = delete;

  /// The kept notices since the last call, in order (index into HeldNotices()); clears them.
  std::vector<size_t> Take();

 protected:
  int overflow(int c) override;
  std::streamsize xsputn(const char *s, std::streamsize n) override;
  int sync() override;

 private:
  void Put(char c);           // caller holds mutex_
  void EndLine();             // caller holds mutex_
  void FlushPartial();        // caller holds mutex_

  std::ostream &out_;
  std::mutex mutex_;
  std::string line_;
  std::vector<size_t> held_;
};

}  // namespace azure2_binding

#endif
