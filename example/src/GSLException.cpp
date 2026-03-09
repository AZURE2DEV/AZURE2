#include "GSLException.h"
#include <gsl/gsl_errno.h>

void GSLException::GSLErrorHandler(const char* reason, 
				   const char* file,
				   int line,
				   int errorCode) {
  if(errorCode==GSL_EOVRFLW) return;
  std::ostringstream stm;
  stm << "GSL Error on line " << line << " of " << file << std::endl
      << reason << std::endl;
  std::string message = stm.str();
  throw GSLException(message);
}
