#ifndef GSLEXCEPTION_H
#define GSLEXCEPTION_H

#include <iostream>
#include <exception>
#include <mutex>
#include <sstream>
#include <gsl/gsl_errno.h>

/*!
 *  The GSLException class is an exception class thrown by the CoulFunc class.
 *  It should not be used directly.
 */

class GSLException : public std::exception {
 public:
  GSLException(std::string message, std::string line = "", std::string file = "") {
    std::ostringstream stm;
    if (line != "" && file != "") {
      stm << "Exception thrown from line " << line << " of file " << file << " with message: " << std::endl
          << message;
    } else {
      stm << "Exception thrown with message: " << std::endl
          << message;
    }
    messageString_ = stm.str();
    message_ = messageString_.c_str();
  };
  ~GSLException() throw() {
  };
  virtual const char *what() const throw() {
    return message_;
  };
  /// GSL error handler that throws a GSLException instead of calling abort().
  static void GSLErrorHandler(const char *, const char *, int, int);

 private:
  std::string messageString_;
  const char *message_;
};

/*!
 * Turns the GSL error handler off for its lifetime, around GSL `_e` calls whose
 * status the caller checks itself.
 *
 * The handler is one process-wide pointer.  A bare save / off / restore run from
 * OpenMP threads interleaves: thread B saves A's "off" and restores it after A
 * has restored the throwing handler, which is then off for good (later Coulomb
 * function errors pass silently), or A restores the throwing handler in the
 * middle of B's quiet call.  Here the first scope in turns it off and the last
 * one out restores what was there, under a lock.
 */
class GslQuiet {
 public:
  GslQuiet() {
    std::lock_guard<std::mutex> lock(Mutex());
    if (Count()++ == 0) Saved() = gsl_set_error_handler_off();
  }
  ~GslQuiet() {
    std::lock_guard<std::mutex> lock(Mutex());
    if (--Count() == 0) gsl_set_error_handler(Saved());
  }
  GslQuiet(const GslQuiet &) = delete;
  GslQuiet &operator=(const GslQuiet &) = delete;

 private:
  static std::mutex &Mutex() {
    static std::mutex m;
    return m;
  }
  static int &Count() {
    static int n = 0;
    return n;
  }
  static gsl_error_handler_t *&Saved() {
    static gsl_error_handler_t *h = nullptr;
    return h;
  }
};

#endif
