#include "Straggling.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

Straggling::Straggling(int element, double energy) 
    : element_(element), initialEnergy_(energy), hasOxideLayer_(false), oxideThickness_(0.15) {
    initializeStraggling();
}

Straggling::~Straggling() {
    // Cleanup if needed
}

void Straggling::initializeStraggling() {
    // Initialize straggling coefficient based on ERYA methodology
    // This is a simplified approximation - in full implementation would use ERYA SRIM data
    // From the notebook example: straggling coefficient ~ 0.04 for Al at ~1 MeV
    stragglingCoefficient_ = 0.04; // keV^0.5 per keV^0.5 of energy loss
}

double Straggling::getStoppingPower(double energy) {
    // Simplified Ziegler stopping power formula
    // In full implementation, this would use ERYA SRIM tables
    // For Al-27: approximate stopping power in MeV cm²/mg
    
    double energyMeV = energy / 1000.0; // Convert keV to MeV
    
    // Simplified Bethe-Bloch approximation for protons in aluminum
    double mass = 27.0; // Al atomic mass
    double Z = element_; // Atomic number
    
    // Very simplified stopping power calculation
    // Real implementation would use SRIM tables from ERYA
    double stoppingPower = 0.1 * Z / energyMeV; // MeV cm²/atoms (very approximate)
    
    return stoppingPower;
}

double Straggling::getStraggling(double deltaE) {
    // Calculate straggling function value based on ERYA methodology
    // From notebook: straggling = coefficient * sqrt(deltaE)
    
    if (deltaE <= 0.0) return 0.0;
    
    double straggling = stragglingCoefficient_ * std::sqrt(deltaE);
    
    // Apply correction factor from notebook (1.2 factor)
    straggling /= 1.2;
    
    return straggling;
}

double Straggling::gammaFunction(double x, double a, double b) {
    // Calculate gamma function distribution as in ERYA notebook
    // Gamma(x; a, b) = (b^a / Gamma(a)) * x^(a-1) * exp(-b*x)
    
    if (x <= 0.0) return 0.0;
    
    double gammaA = gammaApproximation(a);
    double result = std::pow(b, a) * std::pow(x, a - 1.0) * std::exp(-b * x) / gammaA;
    
    return result;
}

double Straggling::gammaApproximation(double z) {
    // Stirling's approximation for gamma function
    // Gamma(z) ≈ sqrt(2π/z) * (z/e)^z for large z
    // For small z, use series expansion or lookup
    
    if (z < 1.0) {
        // For small z, use recurrence relation: Gamma(z) = Gamma(z+1)/z
        return gammaApproximation(z + 1.0) / z;
    } else if (z < 2.0) {
        // Linear interpolation between Gamma(1) = 1 and Gamma(2) = 1
        return 1.0;
    } else {
        // Stirling's approximation
        const double pi = 3.14159265358979323846;
        return std::sqrt(2.0 * pi / z) * std::pow(z / std::exp(1.0), z);
    }
}

void Straggling::convolveWithGaussian(double deltaE, double sigma,
                                     std::vector<double>& energyArray,
                                     std::vector<double>& convolutionArray) {
    // Implement convolution as shown in ERYA notebook
    // This creates a convolved distribution of straggling + Gaussian beam
    
    const int nPoints = 5000;
    const double energyRange = 5.0; // keV range
    
    energyArray.clear();
    convolutionArray.clear();
    energyArray.reserve(nPoints);
    convolutionArray.reserve(nPoints);
    
    // Create energy grid
    for (int i = 0; i < nPoints; ++i) {
        double energy = (i * energyRange / nPoints) + 0.01; // Start from 0.01 keV
        energyArray.push_back(energy);
    }
    
    if (deltaE > 0) {
        // Calculate gamma distribution parameters as in notebook
        double straggling = getStraggling(deltaE);
        double a = deltaE * deltaE / (straggling * straggling); // Shape parameter
        double b = deltaE / (straggling * straggling);          // Rate parameter
        
        // Calculate gamma distribution
        for (double energy : energyArray) {
            double gammaValue = gammaFunction(energy, a, b);
            convolutionArray.push_back(gammaValue);
        }
    } else {
        // Dirac delta function for zero energy loss
        for (int i = 0; i < nPoints; ++i) {
            convolutionArray.push_back(i == 0 ? 1e6 : 0.0);
        }
    }
    
    // Convolve with Gaussian (simplified FFT convolution would be used in practice)
    std::vector<double> tempResult(nPoints);
    const int kernelSize = 1000;
    const double dx = energyRange / nPoints;
    
    // Create Gaussian kernel
    std::vector<double> kernel(kernelSize);
    for (int i = 0; i < kernelSize; ++i) {
        double x = (i - kernelSize/2) * dx;
        kernel[i] = gaussianKernel(x, 0.0, sigma);
    }
    
    // Normalize kernel
    double kernelSum = 0.0;
    for (double k : kernel) kernelSum += k;
    for (double& k : kernel) k /= kernelSum;
    
    // Perform convolution (simplified implementation)
    for (int i = 0; i < nPoints; ++i) {
        double sum = 0.0;
        for (int j = 0; j < kernelSize; ++j) {
            int idx = i - kernelSize/2 + j;
            if (idx >= 0 && idx < nPoints) {
                sum += convolutionArray[idx] * kernel[j];
            }
        }
        tempResult[i] = sum;
    }
    
    convolutionArray = tempResult;
    
    // Normalize result
    double totalArea = 0.0;
    for (double value : convolutionArray) totalArea += value * dx;
    if (totalArea > 0.0) {
        for (double& value : convolutionArray) value /= totalArea;
    }
}

void Straggling::calculateDoubleConvolution(double beamEnergy, double beamSigma,
                                           double targetThickness,
                                           std::vector<double>& resultEnergy,
                                           std::vector<double>& resultIntensity) {
    // Calculate double convolution as in ERYA notebook
    // This combines beam convolution with straggling convolution
    
    // For now, this is a placeholder implementation
    // Full implementation would follow the notebook's double convolution approach
    
    resultEnergy.clear();
    resultIntensity.clear();
    
    const int nPoints = 100;
    const double energyRange = 2.0 * beamSigma; // ±2σ range
    
    for (int i = 0; i < nPoints; ++i) {
        double energy = beamEnergy - energyRange + (2.0 * energyRange * i / nPoints);
        double deltaE = beamEnergy - energy;
        
        // Get straggling contribution
        std::vector<double> stragglingEnergy, stragglingIntensity;
        convolveWithGaussian(deltaE, beamSigma, stragglingEnergy, stragglingIntensity);
        
        // For simplicity, take the peak value
        // Full implementation would properly integrate the double convolution
        auto maxIt = std::max_element(stragglingIntensity.begin(), stragglingIntensity.end());
        double intensity = (maxIt != stragglingIntensity.end()) ? *maxIt : 0.0;
        
        resultEnergy.push_back(energy);
        resultIntensity.push_back(intensity);
    }
}

double Straggling::gaussianKernel(double x, double center, double sigma) {
    // Standard Gaussian function
    const double pi = 3.14159265358979323846;
    double exponent = -(x - center) * (x - center) / (2.0 * sigma * sigma);
    return std::exp(exponent) / (sigma * std::sqrt(2.0 * pi));
}

double Straggling::calculateEnergyLoss(double energy, int element) {
    // Calculate energy loss per unit thickness
    // This would use SRIM data in full implementation
    return getStoppingPower(energy) * 0.1; // Simplified calculation
}

double Straggling::calculateStraggling(double energy, double deltaE) {
    // Calculate straggling variance
    // This implements the ERYA straggling model
    return getStraggling(deltaE);
}

void Straggling::setOxideLayer(bool enable, double oxideThickness) {
    hasOxideLayer_ = enable;
    oxideThickness_ = oxideThickness;
}