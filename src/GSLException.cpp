#include "GSLException.h"
#include <gsl/gsl_errno.h>

void GSLException::GSLErrorHandler(const char* reason, 
				   const char* file,
				   int line,
				   int errorCode) {
  if(errorCode==GSL_EOVRFLW) return;  // overflow - ignore
  if(errorCode==GSL_EUNDRFLW) return; // underflow - ignore
  if(errorCode==GSL_ERANGE) return;   // range error - ignore
  if(errorCode==GSL_ELOSS) return;    // loss of precision - ignore
  std::ostringstream stm;
  stm << "GSL Error on line " << line << " of " << file << std::endl
      << reason << std::endl;
  std::string message = stm.str();
  throw GSLException(message);
}
