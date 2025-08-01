#include "ECAmplitudeCache.h"
#include "ECIntegral.h"
#include "ECIntegralCache.h"
#include "CNuc.h"
#include "Config.h"
#include "CoulFunc.h"
#include <algorithm>
#include <iostream>
#include <cmath>

// Global cache instance
ECAmplitudeCache* g_ecAmplitudeCache = nullptr;

void ECAmplitudeCache::AddAmplitude(const AmplitudeKey& key, double energy, complex amplitude) {
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

complex ECAmplitudeCache::GetInterpolatedAmplitude(const AmplitudeKey& key, double energy, bool calculateIfMissing, 
                                                  CNuc* theCNuc, const Config* configure) const {
    auto it = cache_.find(key);
    if (it == cache_.end()) {
        if (calculateIfMissing && theCNuc && configure) {
            // Calculate and add the missing amplitude
            return CalculateAndAddAmplitude(key, energy, theCNuc, *configure);
        }
        return complex(0.0, 0.0);
    }
    
    const AmplitudeData& data = it->second;
    if (data.energies.empty()) {
        return complex(0.0, 0.0);
    }
    
    // If we only have one data point, return it
    if (data.energies.size() == 1) {
        return data.amplitudes[0];
    }
    
    // Find the closest energy points for interpolation
    int lowerIndex = -1;
    int upperIndex = -1;
    
    for (int i = 0; i < data.energies.size(); i++) {
        if (data.energies[i] <= energy) {
            lowerIndex = i;
        }
        if (data.energies[i] >= energy && upperIndex == -1) {
            upperIndex = i;
            break;
        }
    }
    
    // Handle edge cases
    if (lowerIndex == -1) {
        // Energy is below all cached energies
        if (calculateIfMissing && theCNuc && configure) {
            // Calculate and add the missing amplitude
            std::cout << "ECAmplitudeCache: Energy below cached data, calculating amplitude." << std::endl;
            return CalculateAndAddAmplitude(key, energy, theCNuc, *configure);
        }
        return data.amplitudes[0];
    }
    if (upperIndex == -1) {
        // Energy is above all cached energies
        if (calculateIfMissing && theCNuc && configure) {
            // Calculate and add the missing amplitude
            std::cout << "ECAmplitudeCache: Energy above cached data, calculating amplitude." << std::endl;
            return CalculateAndAddAmplitude(key, energy, theCNuc, *configure);
        }
        return data.amplitudes[data.energies.size() - 1];
    }
    if (lowerIndex == upperIndex) {
        // Exact match
        return data.amplitudes[lowerIndex];
    }
    
    // Linear interpolation between two points
    double e1 = data.energies[lowerIndex];
    double e2 = data.energies[upperIndex];
    complex amp1 = data.amplitudes[lowerIndex];
    complex amp2 = data.amplitudes[upperIndex];
    
    return Interpolate(energy, e1, e2, amp1, amp2);
}

bool ECAmplitudeCache::HasData(const AmplitudeKey& key) const {
    auto it = cache_.find(key);
    return it != cache_.end() && !it->second.energies.empty();
}

void ECAmplitudeCache::Clear() {
    cache_.clear();
}

void ECAmplitudeCache::PrintStats() const {
    std::cout << "ECAmplitudeCache Statistics:" << std::endl;
    std::cout << "  Number of cached keys: " << cache_.size() << std::endl;
    
    int totalPoints = 0;
    for (const auto& pair : cache_) {
        totalPoints += static_cast<int>(pair.second.energies.size());
    }
    std::cout << "  Total cached amplitude points: " << totalPoints << std::endl;
}

complex ECAmplitudeCache::Interpolate(double energy, double e1, double e2, complex amp1, complex amp2) const {
    if (std::abs(e2 - e1) < 1e-10) {
        // Energies are too close, return first amplitude
        return amp1;
    }
    
    double t = (energy - e1) / (e2 - e1);
    return amp1 + t * (amp2 - amp1);
}

complex ECAmplitudeCache::CalculateAndAddAmplitude(const AmplitudeKey& key, double energy, 
                                                 CNuc* theCNuc, const Config& configure) const {
    // This is a complex calculation that mirrors what's done in EPoint::CalculateECAmplitudes
    // We need to find the appropriate parameters and calculate the EC amplitude
    
    // Find entrance pair
    int aa = theCNuc->GetPairNumFromKey(key.entranceKey);
    if (!theCNuc->GetPair(aa)->IsEntrance()) {
        return complex(0.0, 0.0);
    }
    PPair* entrancePair = theCNuc->GetPair(aa);
    
    // Find corresponding EC pathway
    for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
        for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
            if (theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
                ALevel* ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
                int ir = theCNuc->GetPairNumFromKey(key.exitKey);
                if (ecLevel->GetECPairNum() == ir) {
                    // Found matching level, now find the specific pathway
                    for (int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
                        if (k != key.kGroupNum) continue;
                        
                        KGroup* theKGroup = entrancePair->GetDecay(ir)->GetKGroup(k);
                        for (int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
                            if (ecm != key.ecMGroupNum) continue;
                            
                            ECMGroup* theECMGroup = theKGroup->GetECMGroup(ecm);
                            
                            // Calculate entrance phase
                            CoulFunc theCoulombFunction(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
                            struct CoulWaves coul = theCoulombFunction(theECMGroup->GetL(), entrancePair->GetChRad(), energy);
                            
                            double eta = sqrt(uconv/2.) * fstruc * entrancePair->GetZ(1) * entrancePair->GetZ(2) *
                                       sqrt(entrancePair->GetRedMass() / energy);
                            
                            complex expCP(1.0, 0.0);
                            for (int ll = 1; ll <= theECMGroup->GetL(); ll++) {
                                expCP *= complex((double)ll / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)),
                                               eta / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)));
                            }
                            
                            complex expHSP(coul.G / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)),
                                         -coul.F / sqrt(pow(coul.F, 2.0) + pow(coul.G, 2.0)));
                            
                            double inEnergy = energy + entrancePair->GetSepE() + entrancePair->GetExE();
                            double levelEnergy = ecLevel->GetE();
                            double sqrtGammaPene = pow((inEnergy - levelEnergy) / hbarc, theECMGroup->GetMult() + 0.5);
                            
                            // Get integral parameters
                            AChannel* theFinalChannel = theCNuc->GetJGroup(j)->GetChannel(theECMGroup->GetFinalChannel());
                            PPair* theFinalPair = theCNuc->GetPair(theFinalChannel->GetPairNum());
                            
                            int theInitialLValue;
                            double theInitialSValue;
                            if (theECMGroup->IsChannelCapture()) {
                                MGroup* theChanCapMGroup = entrancePair->GetDecay(theECMGroup->GetChanCapDecay())->
                                    GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
                                theInitialLValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->
                                    GetChannel(theChanCapMGroup->GetChpNum())->GetL();
                                theInitialSValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->
                                    GetChannel(theChanCapMGroup->GetChpNum())->GetS();
                            } else {
                                theInitialLValue = theECMGroup->GetL();
                                theInitialSValue = theKGroup->GetS();
                            }
                            
                            // Calculate integral using ECIntegralCache (with forceAdd=true)
                            complex integrals;
                            if (g_ecIntegralCache) {
                                ECIntegralCache::IntegralKey integralKey;
                                integralKey.liValue = theInitialLValue;
                                integralKey.lfValue = theFinalChannel->GetL();
                                integralKey.siValue = theInitialSValue;
                                integralKey.sfValue = theFinalChannel->GetS();
                                integralKey.jInitial = theECMGroup->GetJ();
                                integralKey.jFinal = theCNuc->GetJGroup(j)->GetJ();
                                integralKey.multL = theECMGroup->GetMult();
                                integralKey.radType = theECMGroup->GetRadType();
                                integralKey.bindingEnergy = levelEnergy;
                                integralKey.isChannelCapture = theECMGroup->IsChannelCapture();
                                integralKey.separationEnergy = theFinalPair->GetSepE() + theFinalPair->GetExE();
                                
                                // Use forceAdd=true to add calculated values to integral cache
                                integrals = g_ecIntegralCache->GetIntegral(integralKey, inEnergy, theFinalPair, configure, true);
                            } else {
                                // Fall back to direct calculation
                                ECIntegral theECIntegral(theFinalPair, configure);
                                integrals = theECIntegral(theInitialLValue, theFinalChannel->GetL(),
                                                        theInitialSValue, theFinalChannel->GetS(),
                                                        theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
                                                        theECMGroup->GetMult(), theECMGroup->GetRadType(),
                                                        inEnergy, levelEnergy,
                                                        theECMGroup->IsChannelCapture());
                            }
                            
                            // Calculate final EC amplitude
                            complex ecAmplitude = expCP * expHSP * sqrtGammaPene * integrals;
                            
                            // Add to cache
                            const_cast<ECAmplitudeCache*>(this)->AddAmplitude(key, energy, ecAmplitude);
                            
                            return ecAmplitude;
                        }
                    }
                }
            }
        }
    }
    
    return complex(0.0, 0.0);
}

// Global functions
void InitializeECAmplitudeCache() {
    if (g_ecAmplitudeCache) {
        delete g_ecAmplitudeCache;
    }
    g_ecAmplitudeCache = new ECAmplitudeCache();
}

void CleanupECAmplitudeCache() {
    if (g_ecAmplitudeCache) {
        delete g_ecAmplitudeCache;
        g_ecAmplitudeCache = nullptr;
    }
}