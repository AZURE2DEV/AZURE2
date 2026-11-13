// The standard output of one pyazr session (pybind_notices.h).
#include "pybind_notices.h"

namespace azure2_binding {

const std::vector<HeldNotice> &HeldNotices() {
  // A row per notice: the engine's line, the file it concerns, what follows from it.  The file keys are resolved by
  // Session::take_notices (pybind_api.cpp).
  static const std::vector<HeldNotice> table = {
      {"Could not write to EC Amplitude File.", "ec_amplitudes",
       "The external-capture amplitudes of this session were computed and are used (its results are right), but they "
       "are not saved, so the next session computes them again.  Create the output directory before the session "
       "(write_output_files creates it, too late for this file), or name another one in the second line of "
       "<config>."},
      {"Could not append to EC Amplitude File.", "ec_amplitudes",
       "The external-capture amplitudes of the component segments were computed and are used (the results are "
       "right), but not saved, so the next session computes them again.  Create the output directory before the "
       "session, or name another one in the second line of <config>."},
      {"Could not write normalization file.", "normalizations",
       "The other output files may have been written; this one is missing.  Check that the output directory can be "
       "written, then call write_output_files() again."},
      {"Could not write energy shifts file.", "shifts",
       "The other output files may have been written; this one is missing.  Check that the output directory can be "
       "written, then call write_output_files() again."},
  };
  return table;
}

NoticeBuffer::~NoticeBuffer() {
  std::lock_guard<std::mutex> lock(mutex_);
  FlushPartial();
}

std::vector<size_t> NoticeBuffer::Take() {
  std::lock_guard<std::mutex> lock(mutex_);
  std::vector<size_t> out;
  out.swap(held_);
  return out;
}

int NoticeBuffer::overflow(int c) {
  if (c == traits_type::eof()) return traits_type::not_eof(c);
  std::lock_guard<std::mutex> lock(mutex_);
  Put(static_cast<char>(c));
  return c;
}

std::streamsize NoticeBuffer::xsputn(const char *s, std::streamsize n) {
  std::lock_guard<std::mutex> lock(mutex_);
  for (std::streamsize i = 0; i < n; ++i) Put(s[i]);
  return n;
}

int NoticeBuffer::sync() {
  std::lock_guard<std::mutex> lock(mutex_);
  // A notice is a whole line: a part line that is the start of one waits for its newline.
  bool startOfNotice = false;
  for (const HeldNotice &n : HeldNotices())
    if (!line_.empty() && std::string(n.text).compare(0, line_.size(), line_) == 0) startOfNotice = true;
  if (!startOfNotice) FlushPartial();
  out_.flush();
  return 0;
}

void NoticeBuffer::Put(char c) {
  if (c == '\n')
    EndLine();
  else
    line_.push_back(c);
}

void NoticeBuffer::EndLine() {
  const std::vector<HeldNotice> &notices = HeldNotices();
  for (size_t i = 0; i < notices.size(); ++i) {
    if (line_ == notices[i].text) {
      held_.push_back(i);
      line_.clear();
      return;
    }
  }
  out_ << line_ << '\n';
  line_.clear();
}

void NoticeBuffer::FlushPartial() {
  if (line_.empty()) return;
  out_ << line_;
  line_.clear();
}

}  // namespace azure2_binding
