// The queries of the compound system (CompoundSystem.h): its nuclei, pairs, J^pi groups and their numbering.
#include <cmath>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include "CompoundSystem.h"
#include "RunScratch.h"
#include "RunDiagnostics.h"
#include "azure2/physics/Constants.h"
#include <mutex>
#include <exception>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ParameterLabel.h"
#include "Config.h"
#include "ChannelFunc.h"
#include "EngineChannelFunctions.h"
#include "SystemOfCompound.h"
#include "CompoundSystemChecks.h"
#include "LevelReport.h"
#include "LevelInputMessages.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "GSLException.h"
#include "ProjectSetup.h"
#include "azure2/experiment/FinalStates.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/LevelRows.h"
#include "Minuit2/MnUserParameters.h"
#include "AZUREParams.h"
#include "NFIntegral.h"
#include "ExternalCaptureAmplitudes.h"
#include "ShftFunc.h"
#include "azure2/io/LevelGrouping.h"
#include "azure2/physics/ObservableCoupling.h"
#include "azure2/physics/Pathways.h"

/*!
 * Returns true if a specified pair key exists in the PPair vector, otherwise returns false.
 */

bool CompoundSystem::IsPairKey(int key) const {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumPairs()) {
    if (key == this->GetPair(c + 1)->GetPairKey()) b = true;
    c++;
  }
  return b;
}

std::vector<int> CompoundSystem::CapturePairKeys(int entranceKey) const {
  std::vector<int> keys;
  const PPair *entrance = GetPairFromKey(entranceKey);
  if (!entrance) return keys;
  for (int i = 1; i <= NumPairs(); i++)
    if (GetPair(i)->GetPType() == 10 && GetPair(i)->GetNucleus() == entrance->GetNucleus())
      keys.push_back(GetPair(i)->GetPairKey());
  return keys;
}

std::string CompoundSystem::RatePairsProblem(int entranceKey, int exitKey) const {
  const PPair *entrance = GetPairFromKey(entranceKey);
  if (!entrance) return "Reaction rate entrance pair " + std::to_string(entranceKey) + " does not exist in compound nucleus.";
  if (entrance->GetPType() != 0)
    return "Reaction rate entrance pair " + std::to_string(entranceKey) + " is not a particle pair.";
  if (exitKey == -1) {
    if (CapturePairKeys(entranceKey).empty())
      return "Total capture rate (exit key -1): the compound nucleus of entrance pair " + std::to_string(entranceKey) +
             " has no capture (gamma) pair.";
    return "";
  }
  if (!GetPairFromKey(exitKey)) return "Reaction rate exit pair " + std::to_string(exitKey) + " does not exist in compound nucleus.";
  if (exitKey == entranceKey) return "Cannot calculate rate for elastic scattering.";
  return "";
}

/*!
 * Returns the number of particle pairs in the PPair vector.
 */

int CompoundSystem::NumPairs() const {
  return pairs_.size();
}

/*!
 * Returns the number of \f$ J^\pi \f$ groups in the JGroup vector.
 */

int CompoundSystem::NumJGroups() const {
  return jgroups_.size();
}

/*!
 * Tests if a \f$ J^\pi \f$ group exists in the JGroup vector.  If the group exists, the position in the vector
 * is returned.  Otherwise, the function returns 0.
 */

int CompoundSystem::IsJGroup(JGroup jGroup) const {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumJGroups()) {
    // A group belongs to one compound nucleus: the same J^pi of another
    // nucleus in the file is a different group, with its own level matrix.
    if (jGroup.GetJ() == this->GetJGroup(c + 1)->GetJ() &&
        jGroup.GetPi() == this->GetJGroup(c + 1)->GetPi() &&
        jGroup.GetNucleus() == this->GetJGroup(c + 1)->GetNucleus()) b = true;
    c++;
  }
  if (b)
    return c;
  else
    return 0;
}

/*!
 * The compound nuclei of the project, one per distinct (Z, A) the pairs form.
 */

int CompoundSystem::NumNuclei() const {
  return nuclei_.size();
}

const CompoundNucleus *CompoundSystem::GetNucleus(int nucleusNum) const {
  return &nuclei_[nucleusNum - 1];
}

std::string CompoundSystem::NucleusLabel(int nucleusNum) const {
  std::ostringstream out;
  out << "Z=" << nuclei_[nucleusNum - 1].Z << " A=" << nuclei_[nucleusNum - 1].A;
  return out.str();
}

/*!
 * Returns the position of a particle pair in the PPair vector based on the pair key.
 * Pair keys are how particle pairs are specified in
 * the setup files, but may not correspond to the position of the particle pair in the PPair vector.
 * If the pair exists, the position in the vector is returned.  Otherwise, the function returns 0.
 */

int CompoundSystem::GetPairNumFromKey(int key) const {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumPairs()) {
    if (key == this->GetPair(c + 1)->GetPairKey()) b = true;
    c++;
  }
  if (b)
    return c;
  else
    return 0;
}

/*!
 * Returns a pointer to the particle pair specified by a position in the PPair vector.
 */

PPair *CompoundSystem::Pair(int pairNum) {
  return &pairs_[pairNum - 1];
}

const PPair *CompoundSystem::GetPair(int pairNum) const {
  return &pairs_[pairNum - 1];
}

const PPair *CompoundSystem::GetPairFromKey(int key) const {
  const int num = GetPairNumFromKey(key);
  return num > 0 ? GetPair(num) : nullptr;
}

/*!
 * Returns a pointer to the \f$ J^\pi \f$ group specified by a position in the JGroup vector.
 */

JGroup *CompoundSystem::Group(int jGroupNum) {
  return &jgroups_[jGroupNum - 1];
}

const JGroup *CompoundSystem::GetJGroup(int jGroupNum) const {
  return &jgroups_[jGroupNum - 1];
}

