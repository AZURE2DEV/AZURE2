#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>

/// A structure holding the reaction rate calculation configuration

/*!
 * The RateParams structure holds the configuration information for a reaction rate calculation.
 */

struct RateParams {
  ///False for looped temperatures, true for temperatures from file.
  bool useFile;
  ///String containing filename with temperatures to use
  std::string temperatureFile;
  ///The entrance pair number for the rate calculation
  int entrancePair;
  ///The exit pair number for the rate calculation
  int exitPair;
  ///The minimum temperature for the rate calculation
  double minTemp;
  ///The maximum temperature for the rate calculation
  double maxTemp;
  ///The temperature step for the rate calculation
  double tempStep;
};

/// A configuration structure for AZURE

/*!
 * The configuration structure is created from the runtime file passed to the 
 * AZURE executable, as well as the options specified in the command shell prompt.
 */

class Config {
 public:
  Config(std::ostream& stream);
  void Reset();
  /*!
   * Bit flags for various options in AZURE2.
   */
  enum ParameterFlags {
    USE_AMATRIX               =(1<<0),
    PERFORM_ERROR_ANALYSIS    =(1<<1),
    PERFORM_FIT               =(1<<2),
    CALCULATE_WITH_DATA       =(1<<3),
    USE_PREVIOUS_PARAMETERS   =(1<<4),
    USE_EXTERNAL_CAPTURE      =(1<<5),
    USE_PREVIOUS_INTEGRALS    =(1<<6),
    CALCULATE_REACTION_RATE   =(1<<7),
    TRANSFORM_PARAMETERS      =(1<<8),
    USE_BRUNE_FORMALISM       =(1<<9),
    IGNORE_ZERO_WIDTHS        =(1<<10),
    USE_RMC_FORMALISM         =(1<<11),
    USE_GSL_COULOMB_FUNC      =(1<<12),
    USE_LONGWAVELENGTH_APPROX =(1<<13)
  };
  /*!
   * Bit flags for check file control in AZURE2.
   */
  enum CheckFileFlags {
    CHECK_COMPOUND_NUCLEUS    =(1<<0),
    CHECK_PATHWAYS            =(1<<1),
    CHECK_DATA                =(1<<2),
    CHECK_ENERGY_DEP          =(1<<3),
    CHECK_LEGENDRE            =(1<<4),
    CHECK_BOUNDARY_CONDITIONS =(1<<5),
    CHECK_ANGULAR_DISTS       =(1<<6),
    CHECK_COUL_AMPLITUDES     =(1<<7)
  };
  /// Output stream
  std::ostream &outStream;
  ///The runtime configuration file name.
  std::string configfile;
  ///A control variable to stop AZURE calculation
  bool stopFlag;
  ///A bitmask for the encoding of configuration flags
  unsigned int paramMask;
  ///A bitmask storing which checks are printed to screen.
  unsigned int screenCheckMask;
  ///A bitmask storing which checks are printed to file.
  unsigned int fileCheckMask;
  ///If performError is true, sets the value of Up (the acceptable variance from the minimum chi-squared.
  double chiVariance;
  ///The path of the output directory
  std::string outputdir;
  ///The path of the check files directory.
  std::string checkdir;
  ///The name of the parameters file from which to read.
  std::string paramfile;
  ///The name of the external capture amplitudes file from which to read.
  std::string integralsfile;
  ///Parameters for calculating reaction rate.
  RateParams  rateParams;
  ///A constant indicating the maximum order of the Legendre polynomials to calculate.
  static const int maxLOrder=20;
  int ReadConfigFile();
#ifndef NO_STAT
  int CheckForInputFiles();
#endif
};

#endif
