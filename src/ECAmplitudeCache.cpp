#include "ECAmplitudeCache.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// Global cache instance
ECAmplitudeCache* g_ecAmplitudeCache = nullptr;

void ECAmplitudeCache::AddAmplitude(const AmplitudeKey& key, double energy, complex amplitude) {
    std::lock_guard<std::mutex> lock(cache_mutex_);

    // If amplitude 0 or NaN, do not cache
    if (amplitude == complex(0.0, 0.0) || std::isnan(amplitude.real()) || std::isnan(amplitude.imag())) {
        return;
    }

    auto& data = cache_[key];
    data.energies.push_back(energy);
    data.amplitudes.push_back(amplitude);
    // Ensure energies and amplitudes are sorted together
    auto combined = std::vector<std::pair<double, complex>>(data.energies.size());
    for (size_t i = 0; i < data.energies.size(); ++i) {
        combined[i] = std::make_pair(data.energies[i], data.amplitudes[i]);
    }
    std::sort(combined.begin(), combined.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    for (size_t i = 0; i < combined.size(); ++i) {
        data.energies[i] = combined[i].first;
        data.amplitudes[i] = combined[i].second;
    }

}

complex ECAmplitudeCache::GetInterpolatedAmplitude(const AmplitudeKey& key, double energy) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto it = cache_.find(key);
    if (it == cache_.end()) {
        return complex(0.0, 0.0);
    }
    
    const AmplitudeData& data = it->second;
    if (data.energies.empty()) {
        return complex(0.0, 0.0);
    }
    
    // If we only have one data point, return it
    if (data.energies.size() == 1) {
        complex result = data.amplitudes[0];
        
        return result;
    }
    
    // Find the closest energy points for interpolation
    // Since energies are sorted, we need to find the bracketing indices
    int lowerIndex = -1;
    int upperIndex = -1;
    
    // First check for exact match
    for (int i = 0; i < data.energies.size(); i++) {
        if (std::abs(data.energies[i] - energy) < 1e-10) {
            // Exact match found
            complex result = data.amplitudes[i];
            
            return result;
        }
    }
    
    // Find bracketing indices
    for (int i = 0; i < data.energies.size(); i++) {
        if (data.energies[i] < energy) {
            lowerIndex = i;  // Keep updating to get the highest index < energy
        } else if (data.energies[i] > energy && upperIndex == -1) {
            upperIndex = i;  // Take the first index > energy
            break;
        }
    }
    
    // Handle edge cases - simple extrapolation
    if (lowerIndex == -1) {
        // Energy is below all cached energies - use first value
        complex result = data.amplitudes[0];
        
        return result;
    }
    if (upperIndex == -1) {
        // Energy is above all cached energies - use last value
        complex result = data.amplitudes[data.energies.size() - 1];
        
        return result;
    }
    // Check if we have valid bracketing indices for interpolation
    if (lowerIndex == -1 || upperIndex == -1) {
        // We should have handled these cases above, but as a safety check
        complex result;
        std::string type;
        if (lowerIndex != -1) {
            result = data.amplitudes[lowerIndex];  // Use closest lower value
            type = "FALLBACK_LOWER";
        } else if (upperIndex != -1) {
            result = data.amplitudes[upperIndex];  // Use closest upper value
            type = "FALLBACK_UPPER";
        } else {
            result = data.amplitudes[0];  // Fallback
            type = "FALLBACK_FIRST";
        }
        
        return result;
    }
    
    // Linear interpolation between two points
    double e1 = data.energies[lowerIndex];
    double e2 = data.energies[upperIndex];
    complex amp1 = data.amplitudes[lowerIndex];
    complex amp2 = data.amplitudes[upperIndex];
    
    complex result = Interpolate(energy, e1, e2, amp1, amp2);

    return result;
}

bool ECAmplitudeCache::HasData(const AmplitudeKey& key) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto it = cache_.find(key);
    return it != cache_.end() && !it->second.energies.empty();
}

void ECAmplitudeCache::Clear() {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    cache_.clear();
}

void ECAmplitudeCache::PrintStats() const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    
    int totalPoints = 0;
    for (const auto& pair : cache_) {
        totalPoints += static_cast<int>(pair.second.energies.size());
    }

}

complex ECAmplitudeCache::Interpolate(double energy, double e1, double e2, complex amp1, complex amp2) const {
    if (std::abs(e2 - e1) < 1e-10) {
        // Energies are too close, return first amplitude
        return amp1;
    }
    
    double t = (energy - e1) / (e2 - e1);
    return amp1 + t * (amp2 - amp1);
}


// Global functions
void InitializeECAmplitudeCache() {
    if (g_ecAmplitudeCache) {
        delete g_ecAmplitudeCache;
    }
    g_ecAmplitudeCache = new ECAmplitudeCache();
    
    // Clear any existing cache files to ensure clean start (optional)
    // std::cout << "ECAmplitudeCache: Initialized fresh cache (no file reading)" << std::endl;
}

void CleanupECAmplitudeCache() {
    if (g_ecAmplitudeCache) {
        delete g_ecAmplitudeCache;
        g_ecAmplitudeCache = nullptr;
    }
}