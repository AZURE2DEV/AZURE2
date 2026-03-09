#ifndef AZUREMAIN_H
#define AZUREMAIN_H

#include "CNuc.h"
#include "EData.h"

///The top-level AZURE function class

/*!
 * The AZUREMain function class is the top level function class in the AZURE package.  It is called directly from them main()
 * using the configuration parameters read from the runtime file as well as from the command shell prompt. 
 */

class AZUREMain {
 public:
  /*!
   * The AZUREMain function class is created using a Config structure.  New CNuc and EData objects are created at initialization 
   * of an AZUREMain object.
   */
  AZUREMain(const Config &configure) : configure_(configure) {
    compound_ = new CNuc;
    data_ = new EData;
  };
    /*! 
     * The CNuc and EData objects are destroyed with the AZUREMain instance.
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
  const Config &configure() const {return configure_;};
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a pointer to the EData object.
   */
  EData *data() const {return data_;};
 private:
  const Config &configure_;
  CNuc *compound_;
  EData *data_;
};

#endif
