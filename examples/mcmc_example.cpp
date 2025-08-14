#ifdef USE_MCMC
// Example of how to use AZURE2 with MCMC sampling
#include "AZURECalcMCMC.h"
#include "AZUREMain.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"

int main() {
    // Initialize AZURE2 configuration
    Config config(std::cout);
    
    // Set up compound nucleus and experimental data
    CNuc* compound = new CNuc();
    EData* data = new EData();
    
    // Create MCMC calculator
    AZURECalcMCMC mcmcCalc(data, compound, config);
    
    // Set up initial parameters (these would come from your fit)
    std::vector<double> initialParams = {1.0, 2.0, 0.5}; // Example parameters
    
    // Configure MCMC sampling
    int nWalkers = 50;
    int nSteps = 10000;
    
    // Run MCMC sampling
    std::vector<std::vector<double>> samples;
    mcmcCalc.RunMCMCSampling(nWalkers, nSteps, initialParams, samples);
    
    // Analyze results
    if (!samples.empty()) {
        std::cout << "MCMC sampling completed successfully!" << std::endl;
        std::cout << "Generated " << samples.size() << " samples" << std::endl;
        
        // Calculate parameter statistics
        int nParams = initialParams.size();
        std::vector<double> means(nParams, 0.0);
        std::vector<double> stds(nParams, 0.0);
        
        // Calculate means
        for (const auto& sample : samples) {
            for (int i = 0; i < nParams; i++) {
                means[i] += sample[i];
            }
        }
        for (int i = 0; i < nParams; i++) {
            means[i] /= samples.size();
        }
        
        // Calculate standard deviations
        for (const auto& sample : samples) {
            for (int i = 0; i < nParams; i++) {
                stds[i] += (sample[i] - means[i]) * (sample[i] - means[i]);
            }
        }
        for (int i = 0; i < nParams; i++) {
            stds[i] = sqrt(stds[i] / (samples.size() - 1));
        }
        
        // Print results
        std::cout << "\nParameter Statistics:" << std::endl;
        for (int i = 0; i < nParams; i++) {
            std::cout << "Parameter " << i << ": " 
                     << means[i] << " +/- " << stds[i] << std::endl;
        }
    }
    
    // Cleanup
    delete compound;
    delete data;
    
    return 0;
}
#else
#include <iostream>
int main() {
    std::cout << "MCMC support not enabled. Rebuild with USE_MCMC=ON to use MCMC functionality." << std::endl;
    return 1;
}
#endif