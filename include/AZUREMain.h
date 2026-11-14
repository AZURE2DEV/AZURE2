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
   * The AZUREMain function class is created using a Config structure.  New CompoundSystem and WorkingCopy objects are created at initialization
   * of an AZUREMain object.
   */
  AZUREMain(const Config &configure) :
    configure_(configure) {
    compound_ = new CompoundSystem;
    data_ = new WorkingCopy;
  };
  /*!
   * The CompoundSystem and WorkingCopy objects are destroyed with the AZUREMain instance.
   */
  ~AZUREMain() {
    delete compound_;
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
   * Returns a pointer to the CompoundSystem object.
   */
  const CompoundSystem *compound() const { return compound_; };
  /*!
   * Returns a pointer to the WorkingCopy object.
   */
  WorkingCopy *data() const { return data_; };

 private:
  const Config &configure_;
  CompoundSystem *compound_;
  WorkingCopy *data_;
};

#endif