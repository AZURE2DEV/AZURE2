#ifndef GSLEXCEPTION_H
#define GSLEXCEPTION_H

#include <iostream>
#include <exception>
#include <sstream>

/*!
 *  The GSLException class is an exception class thrown by the CoulFunc class.
 *  It should not be used directly.
 */

class GSLException : public std::exception {
 public:
  GSLException(std::string message, std::string line = "", std::string file = "") {
    std::ostringstream stm;
    if(line != "" && file != "") {
      stm << "Exception thrown from line " << line << " of file " << file << " with message: " << std::endl
	  << message;
    } else {
      stm << "Exception thrown with message: " << std::endl
	  << message;
    }
    messageString_=stm.str();
    message_=messageString_.c_str();
  }; 
  ~GSLException() throw() {
  };
  virtual const char* what() const throw() {
    return message_;
  };
  static void GSLErrorHandler(const char*, const char*, int, int);
 private:
  std::string messageString_;
  const char* message_;
};

#endif
