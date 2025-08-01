#include "ECIntegralCache.h"
#include "ECIntegral.h"
#include "PPair.h"
#include "Config.h"
#include "EData.h"
#include "CNuc.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <set>

// Global cache instance
ECIntegralCache* g_ecIntegralCache = nullptr;

// Comparison operator for IntegralKey
bool ECIntegralCache::IntegralKey::operator<(const IntegralKey& other) const {
    if (liValue != other.liValue) return liValue < other.liValue;
    if (lfValue != other.lfValue) return lfValue < other.lfValue;
    if (std::abs(siValue - other.siValue) > 1e-10) return siValue < other.siValue;
    if (std::abs(sfValue - other.sfValue) > 1e-10) return sfValue < other.sfValue;
    if (std::abs(jInitial - other.jInitial) > 1e-10) return jInitial < other.jInitial;
    if (std::abs(jFinal - other.jFinal) > 1e-10) return jFinal < other.jFinal;
    if (multL != other.multL) return multL < other.multL;
    if (radType != other.radType) return radType < other.radType;
    if (std::abs(bindingEnergy - other.bindingEnergy) > 1e-10) return bindingEnergy < other.bindingEnergy;
    if (isChannelCapture != other.isChannelCapture) return isChannelCapture < other.isChannelCapture;
    return std::abs(separationEnergy - other.separationEnergy) > 1e-10 ? separationEnergy < other.separationEnergy : false;
}

bool ECIntegralCache::IntegralKey::operator==(const IntegralKey& other) const {
    return liValue == other.liValue &&
           lfValue == other.lfValue &&
           std::abs(siValue - other.siValue) < 1e-10 &&
           std::abs(sfValue - other.sfValue) < 1e-10 &&
           std::abs(jInitial - other.jInitial) < 1e-10 &&
           std::abs(jFinal - other.jFinal) < 1e-10 &&
           multL == other.multL &&
           radType == other.radType &&
           std::abs(bindingEnergy - other.bindingEnergy) < 1e-10 &&
           isChannelCapture == other.isChannelCapture &&
           std::abs(separationEnergy - other.separationEnergy) < 1e-10;
}

ECIntegralCache::ECIntegralCache(const std::string& cacheFile)
    : cacheFilename_(cacheFile), cacheInitialized_(false), energyTolerance_(1e-6) {
    if (!cacheFile.empty()) {
        LoadCache(cacheFile);
    }
}

ECIntegralCache::~ECIntegralCache() {
    if (!cacheFilename_.empty() && cacheInitialized_) {
        SaveCache();
    }
}

complex ECIntegralCache::GetIntegral(const IntegralKey& key, double energy, PPair* pair, const Config& configure, bool forceAdd) {
    // Create a separation-energy-aware cache key
    IntegralKey cacheKey = key;
    cacheKey.separationEnergy = pair->GetSepE() + pair->GetExE();
    
    auto it = cache_.find(cacheKey);
    
    // If not in cache, calculate and potentially cache it
    if (it == cache_.end()) {
        // Calculate directly
        ECIntegral ecIntegral(pair, configure);
        complex result = ecIntegral(key.liValue, key.lfValue, key.siValue, key.sfValue, 
                                  key.jInitial, key.jFinal, key.multL, key.radType, 
                                  energy, key.bindingEnergy, key.isChannelCapture, true);  // forceRecalc=true
        
        // If forceAdd is true, add this point to cache for future use
        if (forceAdd) {
            AddIntegralPoint(key, energy, result);
        }
        
        return result;
    }
    
    const CachedIntegrals& cached = it->second;
    
    // Check if energy is within cached range
    if (energy < cached.minEnergy || energy > cached.maxEnergy) {
        // Energy outside cached range, calculate directly
        ECIntegral ecIntegral(pair, configure);
        complex result = ecIntegral(key.liValue, key.lfValue, key.siValue, key.sfValue, 
                                  key.jInitial, key.jFinal, key.multL, key.radType, 
                                  energy, key.bindingEnergy, key.isChannelCapture, true);  // forceRecalc=true
        
        // If forceAdd is true, add this point to cache (extending the range)
        if (forceAdd) {
            AddIntegralPoint(key, energy, result);
        }
        
        return result;
    }
    
    // Use interpolation for cached data - this is now safe because separation energy is part of cache key
    return InterpolateIntegral(cached, energy);
}

void ECIntegralCache::PrecomputeIntegrals(const IntegralKey& key, double minEnergy, double maxEnergy, 
                                         double deltaEnergy, PPair* pair, const Config& configure) {
    CachedIntegrals cached;
    cached.minEnergy = minEnergy;
    cached.maxEnergy = maxEnergy;
    cached.deltaEnergy = deltaEnergy;
    
    // Calculate number of points
    int nPoints = static_cast<int>((maxEnergy - minEnergy) / deltaEnergy) + 1;
    cached.energies.reserve(nPoints);
    cached.integrals.reserve(nPoints);
    
    // Create ECIntegral instance for calculations
    ECIntegral ecIntegral(pair, configure);
    
    // Compute integrals over energy grid
    for (int i = 0; i < nPoints; ++i) {
        double energy = minEnergy + i * deltaEnergy;
        if (energy > maxEnergy) energy = maxEnergy;
        complex integral = ecIntegral(key.liValue, key.lfValue, key.siValue, key.sfValue, 
                                    key.jInitial, key.jFinal, key.multL, key.radType, 
                                    energy, key.bindingEnergy, key.isChannelCapture);
        
        cached.energies.push_back(energy);
        cached.integrals.push_back(integral);
        
        if (energy >= maxEnergy) break;
    }
    
    // Store in cache
    cache_[key] = cached;
    cacheInitialized_ = true;
}

complex ECIntegralCache::InterpolateIntegral(const CachedIntegrals& cached, double energy) const {
    // Find bounding indices
    auto indices = FindEnergyIndices(cached.energies, energy);
    int i1 = indices.first;
    int i2 = indices.second;
    
    // If exact match or at boundary, return stored value
    if (i1 == i2 || std::abs(cached.energies[i1] - energy) < energyTolerance_) {
        return cached.integrals[i1];
    }
    
    // Linear interpolation for complex values
    double e1 = cached.energies[i1];
    double e2 = cached.energies[i2];
    complex int1 = cached.integrals[i1];
    complex int2 = cached.integrals[i2];
    
    double t = (energy - e1) / (e2 - e1);
    return int1 + t * (int2 - int1);
}

std::pair<int, int> ECIntegralCache::FindEnergyIndices(const std::vector<double>& energies, double energy) const {
    // Binary search for energy bounds
    auto it = std::lower_bound(energies.begin(), energies.end(), energy);
    
    if (it == energies.end()) {
        // Energy beyond range
        int idx = static_cast<int>(energies.size()) - 1;
        return std::make_pair(idx, idx);
    }
    
    if (it == energies.begin()) {
        // Energy below range or exact match at start
        return std::make_pair(0, 0);
    }
    
    int i2 = static_cast<int>(it - energies.begin());
    int i1 = i2 - 1;
    
    return std::make_pair(i1, i2);
}

bool ECIntegralCache::LoadCache(const std::string& filename) {
    std::string file = filename.empty() ? cacheFilename_ : filename;
    if (file.empty()) return false;
    
    std::ifstream in(file.c_str());
    if (!in) return false;
    
    cache_.clear();
    std::string line;
    
    // Read header
    if (!std::getline(in, line) || line != "# AZURE2 ECIntegral Cache v1.0") {
        return false;
    }
    
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string keyStr;
        int nPoints;
        
        if (!(iss >> keyStr >> nPoints)) continue;
        
        IntegralKey key;
        if (!ParseCacheKey(keyStr, key)) continue;
        
        CachedIntegrals cached;
        cached.energies.reserve(nPoints);
        cached.integrals.reserve(nPoints);
        
        // Read energy-integral pairs
        for (int i = 0; i < nPoints; ++i) {
            if (!std::getline(in, line)) break;
            std::istringstream dataIss(line);
            double energy, realPart, imagPart;
            
            if (!(dataIss >> energy >> realPart >> imagPart)) break;
            
            cached.energies.push_back(energy);
            cached.integrals.push_back(complex(realPart, imagPart));
        }
        
        if (cached.energies.size() == nPoints && nPoints > 0) {
            cached.minEnergy = cached.energies.front();
            cached.maxEnergy = cached.energies.back();
            if (nPoints > 1) {
                cached.deltaEnergy = (cached.maxEnergy - cached.minEnergy) / (nPoints - 1);
            } else {
                cached.deltaEnergy = 0.0;
            }
            cache_[key] = cached;
        }
    }
    
    cacheInitialized_ = !cache_.empty();
    return cacheInitialized_;
}

bool ECIntegralCache::SaveCache(const std::string& filename) const {
    std::string file = filename.empty() ? cacheFilename_ : filename;
    if (file.empty()) return false;
    
    std::ofstream out(file.c_str());
    if (!out) return false;
    
    out << "# AZURE2 ECIntegral Cache v1.0" << std::endl;
    out << "# Format: key nPoints" << std::endl;
    out << "# Followed by nPoints lines of: energy real_part imag_part" << std::endl;
    out << std::endl;
    
    for (const auto& pair : cache_) {
        const IntegralKey& key = pair.first;
        const CachedIntegrals& cached = pair.second;
        
        std::string keyStr = GenerateCacheKey(key);
        out << keyStr << " " << cached.energies.size() << std::endl;
        
        for (size_t i = 0; i < cached.energies.size(); ++i) {
            out << std::fixed << std::setprecision(10)
                << cached.energies[i] << " "
                << cached.integrals[i].real() << " "
                << cached.integrals[i].imag() << std::endl;
        }
        out << std::endl;
    }
    
    return true;
}

std::string ECIntegralCache::GenerateCacheKey(const IntegralKey& key) const {
    std::ostringstream oss;
    oss << key.liValue << ":" << key.lfValue << ":"
        << std::fixed << std::setprecision(6)
        << key.siValue << ":" << key.sfValue << ":"
        << key.jInitial << ":" << key.jFinal << ":"
        << key.multL << ":" << key.radType << ":"
        << key.bindingEnergy << ":" << (key.isChannelCapture ? 1 : 0);
    return oss.str();
}

bool ECIntegralCache::ParseCacheKey(const std::string& keyStr, IntegralKey& key) const {
    std::istringstream iss(keyStr);
    std::string token;
    std::vector<std::string> tokens;
    
    while (std::getline(iss, token, ':')) {
        tokens.push_back(token);
    }
    
    if (tokens.size() != 10) return false;
    
    try {
        key.liValue = std::stoi(tokens[0]);
        key.lfValue = std::stoi(tokens[1]);
        key.siValue = std::stod(tokens[2]);
        key.sfValue = std::stod(tokens[3]);
        key.jInitial = std::stod(tokens[4]);
        key.jFinal = std::stod(tokens[5]);
        key.multL = std::stoi(tokens[6]);
        key.radType = tokens[7][0];
        key.bindingEnergy = std::stod(tokens[8]);
        key.isChannelCapture = (std::stoi(tokens[9]) != 0);
        return true;
    } catch (...) {
        return false;
    }
}

void ECIntegralCache::ClearCache() {
    cache_.clear();
    cacheInitialized_ = false;
}

bool ECIntegralCache::IsCached(const IntegralKey& key, double energy) const {
    auto it = cache_.find(key);
    if (it == cache_.end()) return false;
    
    const CachedIntegrals& cached = it->second;
    return energy >= cached.minEnergy && energy <= cached.maxEnergy;
}

bool ECIntegralCache::HasKey(const IntegralKey& key) const {
    return cache_.find(key) != cache_.end();
}

void ECIntegralCache::PrintCacheStats(std::ostream& out) const {
    out << "ECIntegral Cache Statistics:" << std::endl;
    out << "  Number of cached parameter sets: " << cache_.size() << std::endl;
    
    int totalPoints = 0;
    for (const auto& pair : cache_) {
        totalPoints += static_cast<int>(pair.second.energies.size());
    }
    out << "  Total cached energy points: " << totalPoints << std::endl;
    out << "  Cache file: " << (cacheFilename_.empty() ? "none" : cacheFilename_) << std::endl;
    out << "  Energy tolerance: " << energyTolerance_ << std::endl;
}

void ECIntegralCache::PopulateCacheForData(EData* data, CNuc* compound, const Config& configure) {
    if (!data || !compound) return;
    
    // Find energy range from all data points
    double minEnergy = 1e10;
    double maxEnergy = -1e10;
    
    for(EDataIterator dataIt = data->begin(); dataIt != data->end(); dataIt++) {
        double energy = dataIt.point()->GetCMEnergy();
        if(energy < minEnergy) minEnergy = energy;
        if(energy > maxEnergy) maxEnergy = energy;
    }
    
    // Expand energy range by 50% on each side to account for potential energy shifts
    double energyRange = maxEnergy - minEnergy;
    double expandedMin = minEnergy - 0.01 * energyRange;
    double expandedMax = maxEnergy + 0.01 * energyRange;
    
    // Use fine energy grid for good interpolation
    double deltaEnergy = energyRange / 1000.0;  // 1000 points across original range
    if(deltaEnergy < 0.001) deltaEnergy = 0.001;
    
    std::cout << "ECIntegralCache: Pre-populating cache from " << expandedMin 
              << " to " << expandedMax << " MeV with " << deltaEnergy << " MeV steps" << std::endl;
    
    // Collect all unique integral keys from entrance pair decays
    std::set<IntegralKey> uniqueKeys;
    
    // Loop through all pairs to find entrance pairs
    for(int pairNum = 1; pairNum <= compound->NumPairs(); pairNum++) {
        PPair* entrancePair = compound->GetPair(pairNum);
        if(!entrancePair->IsEntrance()) continue;
        
        // Loop through all decays of this entrance pair
        for(int decayNum = 1; decayNum <= entrancePair->NumDecays(); decayNum++) {
            Decay* theDecay = entrancePair->GetDecay(decayNum);
            
            // Loop through all KGroups in this decay
            for(int k = 1; k <= theDecay->NumKGroups(); k++) {
                KGroup* theKGroup = theDecay->GetKGroup(k);
                
                // Loop through all ECMGroups in this KGroup  
                for(int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
                    ECMGroup* theECMGroup = theKGroup->GetECMGroup(ecm);
                    
                    // Create integral key
                    IntegralKey key;
                    key.liValue = theECMGroup->GetL(); // Initial L from ECMGroup
                    
                    // Get final channel information
                    AChannel* theFinalChannel = compound->GetJGroup(theECMGroup->GetJGroupNum())
                                                         ->GetChannel(theECMGroup->GetFinalChannel());
                    key.lfValue = theFinalChannel->GetL();
                    key.siValue = theKGroup->GetS();
                    key.sfValue = theFinalChannel->GetS();
                    key.jInitial = theECMGroup->GetJ();
                    key.jFinal = compound->GetJGroup(theECMGroup->GetJGroupNum())->GetJ();
                    key.multL = theECMGroup->GetMult();
                    key.radType = theECMGroup->GetRadType();
                    
                    // Get level energy
                    ALevel* level = compound->GetJGroup(theECMGroup->GetJGroupNum())
                                             ->GetLevel(theECMGroup->GetLevelNum());
                    key.bindingEnergy = level->GetE();
                    key.isChannelCapture = theECMGroup->IsChannelCapture();
                    
                    uniqueKeys.insert(key);
                }
            }
        }
    }
    
    std::cout << "ECIntegralCache: Found " << uniqueKeys.size() << " unique integral keys to cache" << std::endl;
    
    // Pre-compute integrals for each unique key
    int keyCount = 0;
    for(const auto& key : uniqueKeys) {
        keyCount++;
        if(keyCount % 10 == 0) {
            std::cout << "ECIntegralCache: Processed " << keyCount << "/" << uniqueKeys.size() << " keys" << std::endl;
        }
        
        // Find the appropriate PPair for this key (final pair based on ExitKey)
        PPair* theFinalPair = nullptr;
        // This is a simplified approach - we'll use the first non-entrance pair
        // In a full implementation, you'd match based on the ECMGroup's exit pair
        for(int pairNum = 1; pairNum <= compound->NumPairs(); pairNum++) {
            PPair* pair = compound->GetPair(pairNum);
            if(!pair->IsEntrance()) {
                theFinalPair = pair;
                break;
            }
        }
        
        if(theFinalPair) {
            PrecomputeIntegrals(key, expandedMin, expandedMax, deltaEnergy, theFinalPair, configure);
        }
    }
    
    cacheInitialized_ = true;
    std::cout << "ECIntegralCache: Cache pre-population complete" << std::endl;
}

void ECIntegralCache::AddIntegralPoint(const IntegralKey& key, double energy, complex integral) {
    // Create a separation-energy-aware cache key
    IntegralKey cacheKey = key;
    // Note: separationEnergy should already be set in the key by the caller
    
    auto it = cache_.find(cacheKey);
    
    if (it == cache_.end()) {
        // Create a new cache entry with just this point
        CachedIntegrals cached;
        cached.energies.push_back(energy);
        cached.integrals.push_back(integral);
        cached.minEnergy = energy;
        cached.maxEnergy = energy;
        cached.deltaEnergy = 0.0;
        cache_[cacheKey] = cached;
    } else {
        // Add to existing cache entry and maintain sorted order
        CachedIntegrals& cached = it->second;
        
        // Find insertion point to maintain sorted order
        auto insertPos = std::lower_bound(cached.energies.begin(), cached.energies.end(), energy);
        int insertIndex = insertPos - cached.energies.begin();
        
        // Check if energy already exists (within tolerance)
        if (insertPos != cached.energies.end() && std::abs(*insertPos - energy) < energyTolerance_) {
            // Update existing value
            cached.integrals[insertIndex] = integral;
        } else {
            // Insert new point
            cached.energies.insert(insertPos, energy);
            cached.integrals.insert(cached.integrals.begin() + insertIndex, integral);
            
            // Update range
            cached.minEnergy = cached.energies.front();
            cached.maxEnergy = cached.energies.back();
            if (cached.energies.size() > 1) {
                cached.deltaEnergy = (cached.maxEnergy - cached.minEnergy) / (cached.energies.size() - 1);
            }
        }
    }
    
    cacheInitialized_ = true;
}

// Global functions
void InitializeECIntegralCache(const std::string& cacheFile) {
    if (g_ecIntegralCache) {
        delete g_ecIntegralCache;
    }
    g_ecIntegralCache = new ECIntegralCache(cacheFile);
}

void CleanupECIntegralCache() {
    if (g_ecIntegralCache) {
        delete g_ecIntegralCache;
        g_ecIntegralCache = nullptr;
    }
}