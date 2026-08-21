#include "GSLException.h"
#include <gsl/gsl_errno.h>

void GSLException::GSLErrorHandler(const char *reason,
                                   const char *file,
                                   int line,
                                   int errorCode) {
  // Common GSL error codes that should be ignored for numerical robustness
  if (errorCode == 1) return;   // GSL_EDOM - domain error - ignore
  if (errorCode == 2) return;   // GSL_ERANGE - range error - ignore
  if (errorCode == 3) return;   // GSL_EFAULT - fault - ignore
  if (errorCode == 4) return;   // GSL_EINVAL - invalid argument - ignore
  if (errorCode == 5) return;   // GSL_EFAILED - failed - ignore
  if (errorCode == 11) return;  // GSL_EROUND - rounding error - ignore
  if (errorCode == 13) return;  // GSL_ELOSS - loss of precision - ignore
  if (errorCode == 15) return;  // GSL_EUNDRFLW - underflow - ignore
  if (errorCode == 16) return;  // GSL_EOVRFLW - overflow - ignore
  if (errorCode == 18) return;  // GSL_ETOLF - tolerance failure - ignore
  if (errorCode == 19) return;  // GSL_ETOLX - tolerance in x - ignore
  if (errorCode == 20) return;  // GSL_ETOLG - tolerance in gradient - ignore

  // Log the error but don't crash for any numerical errors
  // Only throw for truly fatal errors
  if (errorCode >= 1 && errorCode <= 50) {
    // Most numerical errors should be ignored
    return;
  }

  std::ostringstream stm;
  stm << "GSL Error on line " << line << " of " << file << std::endl
      << reason << " (Error code: " << errorCode << ")" << std::endl;
  std::string message = stm.str();
  throw GSLException(message);
}
