#ifndef AZUREMAIN_H
#define AZUREMAIN_H

#include "CompoundSystem.h"
#include "WorkingCopy.h"

/// The top-level AZURE function class

/*!
 * The AZUREMain function class is the top level function class in the AZURE package.  It is called directly from them main()
 * using the configuration parameters read from the runtime file as well as from the command shell prompt.
 */

class AZUREMain {
 public:
  /*!
   * The AZUREMain function class is created using a Config structure.  A new WorkingCopy object is created at
   * initialization of an AZUREMain object; the compound system is the run's (ProgramRun), made by its set-up.
   */
  AZUREMain(const Config &configure) :
    configure_(configure) {
    data_ = new WorkingCopy;
  };
  /*!
   * The WorkingCopy object is destroyed with the AZUREMain instance.
   */
  ~AZUREMain() {
    delete data_;
  };
  /*!
   * The parenthesis operator is defined so the instance of AZUREMain can be called as a function.  This executes
   * AZURE against the configuration parameters.
   */
  int operator()();
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const { return configure_; };
  /*!
   * Returns a pointer to the WorkingCopy object.
   */
  WorkingCopy *data() const { return data_; };

 private:
  const Config &configure_;
  WorkingCopy *data_;
};

#endif