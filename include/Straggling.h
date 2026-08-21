#ifndef STRAGGLING_H
#define STRAGGLING_H

#include <vector>
#include <cmath>

/**
 * @class Straggling
 * @brief Handles energy straggling calculations using ERYA-Profiler methodology
 *
 * This class implements energy straggling calculations for nuclear reaction analysis,
 * based on the ERYA-Profiler approach. It provides functionality for calculating
 * energy loss straggling functions and convolution with beam distributions.
 */
class Straggling {
 public:
  /**
   * @brief Constructor for Straggling class
   * @param element Atomic number of the target element
   * @param energy Initial beam energy in keV
   */
  Straggling(int element, double energy);

  /**
   * @brief Destructor
   */
  ~Straggling();

  /**
   * @brief Get stopping power for given energy
   * @param energy Energy in keV
   * @return Stopping power in MeV cm²/atoms
   */
  double getStoppingPower(double energy);

  /**
   * @brief Calculate straggling function value for given energy loss
   * @param deltaE Energy loss in keV
   * @return Straggling function value
   */
  double getStraggling(double deltaE);

  /**
   * @brief Calculate gamma function for straggling distribution
   * @param x Energy variable
   * @param a Shape parameter
   * @param b Rate parameter
   * @return Gamma function value
   */
  double gammaFunction(double x, double a, double b);

  /**
   * @brief Convolve straggling with Gaussian beam distribution
   * @param deltaE Energy loss in keV
   * @param sigma Beam width in keV
   * @param energyArray Output energy array
   * @param convolutionArray Output convolution array
   * @param arraySize Size of output arrays
   */
  void convolveWithGaussian(double deltaE, double sigma,
                            std::vector<double> &energyArray,
                            std::vector<double> &convolutionArray);

  /**
   * @brief Calculate double convolution for target effect
   * @param beamEnergy Beam energy in keV
   * @param beamSigma Beam width in keV
   * @param targetThickness Target thickness in keV
   * @param resultEnergy Output energy array
   * @param resultIntensity Output intensity array
   */
  void calculateDoubleConvolution(double beamEnergy, double beamSigma,
                                  double targetThickness,
                                  std::vector<double> &resultEnergy,
                                  std::vector<double> &resultIntensity);

  /**
   * @brief Set target element
   * @param element Atomic number
   */
  void setElement(int element) { element_ = element; }

  /**
   * @brief Set initial energy
   * @param energy Energy in keV
   */
  void setEnergy(double energy) { initialEnergy_ = energy; }

  /**
   * @brief Enable/disable oxide layer consideration
   * @param enable Enable oxide layer
   * @param oxideThickness Oxide thickness in keV (default 0.15)
   */
  void setOxideLayer(bool enable, double oxideThickness = 0.15);

 private:
  int element_;            ///< Target element atomic number
  double initialEnergy_;   ///< Initial beam energy in keV
  bool hasOxideLayer_;     ///< Whether to include oxide layer
  double oxideThickness_;  ///< Oxide layer thickness in keV

  // Straggling parameters from ERYA fitting
  double stragglingCoefficient_;  ///< Fitted straggling coefficient

  /**
   * @brief Initialize straggling coefficient based on element
   */
  void initializeStraggling();

  /**
   * @brief Calculate energy loss through target
   * @param energy Current energy in keV
   * @param element Target element
   * @return Energy loss per thickness step in keV
   */
  double calculateEnergyLoss(double energy, int element);

  /**
   * @brief Calculate straggling variance for given conditions
   * @param energy Energy in keV
   * @param deltaE Energy loss in keV
   * @return Straggling variance
   */
  double calculateStraggling(double energy, double deltaE);

  /**
   * @brief Gaussian kernel for convolution
   * @param x Energy point
   * @param center Center energy
   * @param sigma Width parameter
   * @return Gaussian value
   */
  double gaussianKernel(double x, double center, double sigma);

  /**
   * @brief Calculate gamma function using Stirling's approximation
   * @param z Input value
   * @return Gamma function value
   */
  double gammaApproximation(double z);
};

#endif  // STRAGGLING_H