#ifndef AZURE2_SUPPORT_CLASSICNUMERIC_H
#define AZURE2_SUPPORT_CLASSICNUMERIC_H

#include <clocale>
#include <string>

#if defined(_WIN32)
#include <locale.h>
#else
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif

namespace azure2 {
namespace support {

/*!
 * The "C" numeric locale for the calling thread while the object lives (docs/upgrade/qt_gui.md section 19.10).  The
 * project's and the data files' numbers are written with a decimal point, and the readers parse them with strtod,
 * std::stod and printf-style writers, which follow LC_NUMERIC: a program that sets the user's locale (a Qt application
 * does, from the environment) and runs with a decimal-comma language (de_DE, it_IT, pl_PL, ...) would read "5.441" as
 * 5 and write "5,441".  The command-line program never sets the locale and is not affected; the front ends' facade
 * (api/front) takes this scope at each call that reads, checks, edits, converts, saves or runs a project.
 *
 * Only the calling thread changes (uselocale; on Windows a per-thread locale), and only LC_NUMERIC: other threads,
 * the process's global locale and every other category stay as they were.  In the "C" locale it changes nothing.
 * Threads the call starts (the run's OpenMP threads) keep the process's locale.
 */
class ClassicNumericScope {
 public:
#if defined(_WIN32)
  ClassicNumericScope() {
    mode_ = _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
    const char *now = std::setlocale(LC_NUMERIC, nullptr);
    old_ = now ? now : "C";
    std::setlocale(LC_NUMERIC, "C");
  }
  ~ClassicNumericScope() {
    std::setlocale(LC_NUMERIC, old_.c_str());
    _configthreadlocale(mode_);
  }

 private:
  int mode_ = 0;
  std::string old_;
#else
  ClassicNumericScope() {
    // the thread's locale as it is, with LC_NUMERIC "C"
    locale_t base = duplocale(uselocale((locale_t)0));
    classic_ = base ? newlocale(LC_NUMERIC_MASK, "C", base) : (locale_t)0;
    if (!classic_ && base) freelocale(base);
    if (classic_) old_ = uselocale(classic_);
  }
  ~ClassicNumericScope() {
    if (!classic_) return;
    uselocale(old_);
    freelocale(classic_);
  }

 private:
  locale_t classic_ = (locale_t)0;
  locale_t old_ = (locale_t)0;
#endif

 public:
  ClassicNumericScope(const ClassicNumericScope &) = delete;
  ClassicNumericScope &operator=(const ClassicNumericScope &) = delete;
};

}  // namespace support
}  // namespace azure2

#endif
