#ifndef AZURE2_SUPPORT_FLOATING_POINT_PRECISION_H
#define AZURE2_SUPPORT_FLOATING_POINT_PRECISION_H

/*!
 * The x87 precision of a thread for the length of a scope: the engine computes `long double` with a 64-bit mantissa in
 * every thread that runs it, whoever made the thread.
 *
 * On x86 and x86-64 `long double` is the x87 80-bit format, and how many mantissa bits an x87 operation keeps is the
 * precision control (bits 8-9) of the thread's control word: 64 bits (0x37F, Linux's and macOS's default in every
 * thread) or 53 (0x27F).  On Windows a thread's control word is not inherited: MinGW-w64's start-up (_fpreset) gives
 * the program's main thread 64 bits, and a thread made by another runtime (a Python threading.Thread, made by the
 * Windows C runtime) starts with 53.  The Windows CI saw it: the same pyazr call gave other bits in a Python thread
 * than in the main thread (tests/pyazr/session_lock_test.py (b), objective_object_test.py (d); on Linux the same
 * difference, 27 of 442 residuals of T1 by up to 1.4e-12 relative, appears when a thread's precision is lowered to
 * 53 bits: tests/pyazr/x87_precision_test.py).  MinGW-w64's own mathematical functions use the x87 too.
 *
 * ExtendedPrecisionScope sets the precision control of the calling thread to 64 bits and puts the control word back
 * as it found it when the scope ends; the other bits (rounding, exception masks) are not touched.  It holds the word
 * it found and nothing else (no static state).  Where it is taken: every entry of the Python binding into the engine
 * (api/src/pybind_engine_call.h, EngineCallScope) and the first statement of the body of every OpenMP parallel region
 * of the engine, so that a worker thread of the OpenMP runtime computes as its master does whatever the runtime gave
 * it at its start [repo_consistency check 60].
 *
 * It is compiled on x86 and x86-64 with GCC or clang (Linux, macOS, MinGW), where it changes nothing in a thread
 * already at 64 bits (no number of a Linux run moves), and lets a test lower a thread's precision on any of them;
 * elsewhere (no x87: arm64) it is empty.  X87ControlWord() reads the calling thread's word (0 where there is none),
 * for the tests and the diagnostics.
 */
namespace azure2 {
namespace support {

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__i386__) || defined(__x86_64__))
#define AZURE2_X87_PRECISION_SCOPE 1

/// The x87 control word of the calling thread.
inline unsigned short X87ControlWord() noexcept {
  unsigned short word = 0;
  __asm__ __volatile__("fnstcw %0" : "=m"(word) : : "memory");
  return word;
}

/// Load `word` into the x87 control word of the calling thread.
inline void SetX87ControlWord(unsigned short word) noexcept {
  __asm__ __volatile__("fldcw %0" : : "m"(word) : "memory");
}

/// The x87 precision control at 64 bits for the length of the scope (see the head of this file).
class ExtendedPrecisionScope {
 public:
  ExtendedPrecisionScope() noexcept : found_(X87ControlWord()) {
    if ((found_ & kPrecisionBits) != kPrecisionBits) SetX87ControlWord(found_ | kPrecisionBits);
  }
  ~ExtendedPrecisionScope() {
    if ((found_ & kPrecisionBits) != kPrecisionBits) SetX87ControlWord(found_);
  }
  ExtendedPrecisionScope(const ExtendedPrecisionScope &) = delete;
  ExtendedPrecisionScope &operator=(const ExtendedPrecisionScope &) = delete;
  /// The control word the scope found (the diagnostics of api/src/pybind_precision.cpp).
  unsigned short Found() const noexcept { return found_; }

 private:
  static constexpr unsigned short kPrecisionBits = 0x0300;  // precision control 11: a 64-bit mantissa
  unsigned short found_;  // the calling thread's control word when the scope began
};

#else

inline unsigned short X87ControlWord() noexcept { return 0; }
inline void SetX87ControlWord(unsigned short) noexcept {}

class ExtendedPrecisionScope {
 public:
  ExtendedPrecisionScope() noexcept {}
  ~ExtendedPrecisionScope() {}
  ExtendedPrecisionScope(const ExtendedPrecisionScope &) = delete;
  ExtendedPrecisionScope &operator=(const ExtendedPrecisionScope &) = delete;
  unsigned short Found() const noexcept { return 0; }
};

#endif

}  // namespace support
}  // namespace azure2

#endif
