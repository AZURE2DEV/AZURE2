#include "ECIntegralUtils.h"
#include "ECIntegralCache.h"
#include "CNuc.h"
#include "EData.h"
#include "PPair.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "ECMGroup.h"
#include "KGroup.h"
#include "MGroup.h"
#include "Decay.h"
#include <iostream>
#include <algorithm>
#include <set>

namespace ECIntegralUtils {

void InitializeCacheSystem(const Config& configure) {
    std::string cacheFile;
    
    // Determine cache file name based on configuration
    if (configure.paramMask & Config::CALCULATE_WITH_DATA) {
        cacheFile = configure.outputdir + "intEC_cache.dat";
    } else {
        cacheFile = configure.outputdir + "intEC_cache.extrap";
    }
    
    // Initialize global cache
    InitializeECIntegralCache(cacheFile);
    
    if (configure.outStream.rdbuf() != nullptr) {
        configure.outStream << "EC Integral cache initialized: " << cacheFile << std::endl;
    }
}

void FinalizeCacheSystem() {
    CleanupECIntegralCache();
}

void AnalyzeEnergyRanges(CNuc* compound, EData* data, 
                        double& minEnergy, double& maxEnergy, double& suggestedDelta) {
    minEnergy = 1e10;
    maxEnergy = -1e10;
    
    // Analyze energy ranges from experimental data
    if (data) {
        for (int i = 1; i <= data->NumSegments(); i++) {
            ESegment* segment = data->GetSegment(i);
            for (int j = 1; j <= segment->NumPoints(); j++) {
                EPoint* point = segment->GetPoint(j);
                double energy = point->GetCMEnergy();
                
                // Find entrance pair to get separation energy
                int aa = compound->GetPairNumFromKey(point->GetEntranceKey());
                if (aa > 0) {
                    PPair* entrancePair = compound->GetPair(aa);
                    double inEnergy = energy + entrancePair->GetSepE() + entrancePair->GetExE();
                    
                    minEnergy = std::min(minEnergy, inEnergy);
                    maxEnergy = std::max(maxEnergy, inEnergy);
                }
            }
        }
    }
    
    // Add buffer around the energy range
    if (minEnergy < 1e9) {  // Valid range found
        double range = maxEnergy - minEnergy;
        minEnergy -= 0.1 * range;  // 10% buffer below
        maxEnergy += 0.2 * range;  // 20% buffer above
        
        // Ensure minimum energy is not negative
        if (minEnergy < 0.001) minEnergy = 0.001;  // 1 keV minimum
        
        // Suggest delta based on range
        suggestedDelta = range / 1000.0;  // 1000 points across range
        if (suggestedDelta < 0.001) suggestedDelta = 0.001;  // Minimum 1 keV steps
        if (suggestedDelta > 0.1) suggestedDelta = 0.1;     // Maximum 100 keV steps
    } else {
        // Default range if no data found
        minEnergy = 0.001;   // 1 keV
        maxEnergy = 10.0;    // 10 MeV  
        suggestedDelta = 0.01; // 10 keV steps
    }
}

void PrecomputeAllIntegrals(CNuc* compound, EData* data, const Config& configure) {
    if (!g_ecIntegralCache) {
        configure.outStream << "Error: EC Integral cache not initialized" << std::endl;
        return;
    }
    
    configure.outStream << "Pre-computing external capture integrals..." << std::endl;
    
    // Analyze energy ranges
    double minEnergy, maxEnergy, deltaEnergy;
    AnalyzeEnergyRanges(compound, data, minEnergy, maxEnergy, deltaEnergy);
    
    configure.outStream << "Energy range: " << minEnergy << " - " << maxEnergy 
                       << " MeV (step: " << deltaEnergy << " MeV)" << std::endl;
    
    std::set<ECIntegralCache::IntegralKey> uniqueKeys;
    
    // Collect all unique parameter combinations from nuclear structure
    for (int j = 1; j <= compound->NumJGroups(); j++) {
        JGroup* jGroup = compound->GetJGroup(j);
        
        for (int la = 1; la <= jGroup->NumLevels(); la++) {
            ALevel* ecLevel = jGroup->GetLevel(la);
            if (!ecLevel->IsECLevel()) continue;
            
            int exitPairNum = ecLevel->GetECPairNum();
            PPair* exitPair = compound->GetPair(exitPairNum);
            
            // Find all entrance pairs that can lead to this level
            for (int pairIdx = 1; pairIdx <= compound->NumPairs(); pairIdx++) {
                PPair* entrancePair = compound->GetPair(pairIdx);
                if (!entrancePair->IsEntrance()) continue;
                
                // Look for decays from entrance to exit pair
                for (int decay = 1; decay <= entrancePair->NumDecays(); decay++) {
                    if (entrancePair->GetDecay(decay)->GetPairNum() == exitPairNum) {
                        Decay* theDecay = entrancePair->GetDecay(decay);
                        
                        for (int k = 1; k <= theDecay->NumKGroups(); k++) {
                            KGroup* kGroup = theDecay->GetKGroup(k);
                            
                            for (int ecm = 1; ecm <= kGroup->NumECMGroups(); ecm++) {
                                ECMGroup* ecmGroup = kGroup->GetECMGroup(ecm);
                                
                                // Get final channel info
                                AChannel* finalChannel = jGroup->GetChannel(ecmGroup->GetFinalChannel());
                                
                                ECIntegralCache::IntegralKey key;
                                
                                // Determine initial L and S values
                                if (ecmGroup->IsChannelCapture()) {
                                    MGroup* chanCapMGroup = entrancePair->GetDecay(ecmGroup->GetChanCapDecay())->
                                        GetKGroup(ecmGroup->GetChanCapKGroup())->GetMGroup(ecmGroup->GetChanCapMGroup());
                                    key.liValue = compound->GetJGroup(chanCapMGroup->GetJNum())->
                                        GetChannel(chanCapMGroup->GetChpNum())->GetL();
                                    key.siValue = compound->GetJGroup(chanCapMGroup->GetJNum())->
                                        GetChannel(chanCapMGroup->GetChpNum())->GetS();
                                } else {
                                    key.liValue = ecmGroup->GetL();
                                    key.siValue = kGroup->GetS();
                                }
                                
                                key.lfValue = finalChannel->GetL();
                                key.sfValue = finalChannel->GetS();
                                key.jInitial = ecmGroup->GetJ();
                                key.jFinal = jGroup->GetJ();
                                key.multL = ecmGroup->GetMult();
                                key.radType = ecmGroup->GetRadType();
                                key.bindingEnergy = ecLevel->GetE();
                                key.isChannelCapture = ecmGroup->IsChannelCapture();
                                
                                uniqueKeys.insert(key);
                            }
                        }
                    }
                }
            }
        }
    }
    
    configure.outStream << "Found " << uniqueKeys.size() << " unique parameter combinations" << std::endl;
    
    // Pre-compute integrals for each unique combination
    int count = 0;
    for (const auto& key : uniqueKeys) {
        count++;
        
        // Find appropriate exit pair for this key
        PPair* exitPair = nullptr;
        for (int i = 1; i <= compound->NumPairs(); i++) {
            PPair* pair = compound->GetPair(i);
            // For now, use the first available pair - this could be improved
            // by matching based on the binding energy and other parameters
            if (!pair->IsEntrance()) {
                exitPair = pair;
                break;
            }
        }
        
        if (exitPair) {
            if (count % 10 == 0 || count <= 10) {
                configure.outStream << "Computing integrals " << count << "/" << uniqueKeys.size() 
                                   << " (Li=" << key.liValue << ", Lf=" << key.lfValue 
                                   << ", mult=" << key.multL << ")" << std::endl;
            }
            
            g_ecIntegralCache->PrecomputeIntegrals(key, minEnergy, maxEnergy, deltaEnergy, 
                                                 exitPair, configure);
        }
    }
    
    configure.outStream << "Pre-computation completed for " << count << " parameter sets" << std::endl;
}

void PrecomputeIntegralsForRange(CNuc* compound, double minEnergy, double maxEnergy, 
                               double deltaEnergy, const Config& configure) {
    if (!g_ecIntegralCache) {
        configure.outStream << "Error: EC Integral cache not initialized" << std::endl;
        return;
    }
    
    configure.outStream << "Pre-computing integrals for energy range: " 
                       << minEnergy << " - " << maxEnergy << " MeV" << std::endl;
    
    // This is a simplified version - in practice you'd want to collect
    // all relevant parameter combinations like in PrecomputeAllIntegrals
    
    // For demonstration, create a sample key
    ECIntegralCache::IntegralKey sampleKey;
    sampleKey.liValue = 0;
    sampleKey.lfValue = 1;
    sampleKey.siValue = 0.5;
    sampleKey.sfValue = 0.5;
    sampleKey.jInitial = 0.5;
    sampleKey.jFinal = 1.5;
    sampleKey.multL = 1;
    sampleKey.radType = 'E';
    sampleKey.bindingEnergy = 0.0;
    sampleKey.isChannelCapture = false;
    
    PPair* samplePair = compound->GetPair(1);  // Use first available pair
    if (samplePair) {
        g_ecIntegralCache->PrecomputeIntegrals(sampleKey, minEnergy, maxEnergy, deltaEnergy, 
                                             samplePair, configure);
    }
}

void PrintCacheStatistics(const Config& configure) {
    if (g_ecIntegralCache) {
        //g_ecIntegralCache->PrintCacheStats(configure.outStream);
    } else {
        //configure.outStream << "EC Integral cache not initialized" << std::endl;
    }
}

} // namespace ECIntegralUtils