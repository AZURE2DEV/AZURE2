// The build of the compound system (CompoundSystem.h): the structure from the <levels> rows, the hard-sphere groups,
// the boundary conditions, the input transformation, the pathway table and the external-capture set-up.
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

namespace {

//! 2j as an integer, for the physics/ coupling coefficients (docs/upgrade/cnuc_split.md section 9).
int Twice(double j) { return (int)std::lround(2.0 * j); }

}  // namespace

/*!
 * Tests if a particle pair exists in the PPair vector. If pair exists, the position in
 * the vector is returned.  Otherwise, the function returns 0.
 */

int CompoundSystem::IsPair(PPair pair) {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumPairs()) {
    if (pair.GetPairKey() == this->Pair(c + 1)->GetPairKey()) b = true;
    c++;
  }
  if (b)
    return c;
  else
    return 0;
}

int CompoundSystem::FindOrAddNucleus(int z, int a) {
  for (unsigned int i = 0; i < nuclei_.size(); i++)
    if (nuclei_[i].Z == z && nuclei_[i].A == a) return i + 1;
  CompoundNucleus nucleus;
  nucleus.Z = z;
  nucleus.A = a;
  nuclei_.push_back(nucleus);
  return nuclei_.size();
}

/*!
 * Fills the compound nucleus object, and all nested objects, from data specified in the nuclear and external capture
 * input files.  Returns -1 if the files could not be read, and 0 if the files were read successfully.
 */

int CompoundSystem::Fill(const Config &configure, LevelState &state, LevelInputs &inputs, const std::map<int, double> &radii) {
  inputs = LevelInputs();
  if (configure.relativistic && configure.useHybridMethod) {
    // The hybrid potential acts in a Schrodinger equation between a_c and the matching radius: no relativistic form.
    configure.outStream << "ERROR: relativistic=1 cannot be combined with the hybrid nuclear-potential "
                           "Coulomb functions (a <potential> block)." << std::endl;
    return -1;
  }
  // Every pair of the run computes its functions with the run's settings (PPair.h).
  PairFunctionSettings pairSettings;
  pairSettings.relativistic = configure.relativistic;
  pairSettings.useCoulombTable = configure.useCoulombTable;
  pairSettings.coulombTablePath = configure.coulombTablePath;
  pairSettings.coulombTableDefault = configure.environment.coulombTableDefault;
  pairSettings.boundShiftNumerical = configure.environment.boundShiftNumerical;
  pairSettings.noNeutralClosedForm = configure.environment.noNeutralClosedForm;
  pairSettings.shiftDerivative1eV = configure.environment.shiftDerivative1eV;
  pairSettings.tableUse = &configure.diagnostics->coulombTable;
  int PairNum, LevelNum, ChannelNum, JGroupNum;
  int maxLValue = 0;
  std::map<std::pair<int, int>, std::vector<int> > channelsListed;  // (group, level) -> group channel numbers, in line order
  azure2::io::LevelGrouper levelGrouper;  // the level of each line (io/LevelGrouping.h)
  // The rows of the <levels> section, read once (ProjectSetup.h: a format-3 project's by Config::ReadConfigFile).
  azure2::io::LevelSection section;
  if (!ProjectLevelSection(configure, section)) return -1;  // no <levels> section
  std::map<int, int> ecPairs;
  // Which compound nucleus each level (by its ID, column 9) belongs to: the
  // channel lines of one level must all come from pairs of one nucleus.
  std::map<int, int> levelNucleus;
  std::map<int, double> lightExcitation;  // pair key -> e1 of the first line naming it
  for (const azure2::io::LevelRow &Line : section.rows) {
    // The light particle's excitation energy, the optional 32nd column (2026-11-05).  A pair's columns are
    // read from the first line that names it, so a different e1 on a later line of the same pair would be
    // dropped without a word -- and it is a different threshold and a different identical-particle test, so
    // refuse it.  A photon or an electron has no excited state.
    if (Line.e1() < 0.) {
      configure.outStream << "ERROR: pair " << Line.ir() << " has a negative light-particle excitation energy (column 32), "
                          << Line.e1() << " MeV." << std::endl;
      return -1;
    }
    if (Line.e1() != 0. && (Line.pType() == 10 || Line.pType() == 20)) {
      configure.outStream << "ERROR: pair " << Line.ir() << " is a " << (Line.pType() == 10 ? "gamma" : "beta")
                          << " pair (type " << Line.pType() << ") but has a light-particle excitation energy (column 32) of "
                          << Line.e1() << " MeV; only a particle pair can have one." << std::endl;
      return -1;
    }
    {
      std::map<int, double>::const_iterator seen = lightExcitation.find(Line.ir());
      if (seen == lightExcitation.end()) {
        lightExcitation[Line.ir()] = Line.e1();
      } else if (fabs(seen->second - Line.e1()) > 1.0e-9) {
        configure.outStream << "ERROR: pair " << Line.ir() << " has two light-particle excitation energies (column 32), "
                            << seen->second << " and " << Line.e1() << " MeV, on different <levels> lines; a pair's "
                            << "columns must be the same on every line." << std::endl;
        return -1;
      }
    }
    if (Line.l() > maxLValue && Line.pType() == 0) maxLValue = Line.l();
    if (Line.isActive() != 1 && configure.includeAllJPi) {
      // A pair is a property of the system, not of the level that happens to mention it: with every
      // J^pi scattering through its hard sphere, a pair whose levels are all switched off still has
      // channels, so it is kept.
      PPair NewPair(Line, pairSettings);
      if (!this->IsPair(NewPair)) this->AddPair(NewPair);
    }
    if (Line.isActive() == 1) {
      PPair NewPair(Line, pairSettings);
      PairNum = this->IsPair(NewPair);
      if (!PairNum) {
        this->AddPair(NewPair);
        PairNum = this->IsPair(NewPair);
      }
      if (Line.ecMultMask() != 0) {
        std::map<int, int>::iterator it = ecPairs.find(PairNum);
        if (it == ecPairs.end()) ecPairs[PairNum] = Line.ecMultMask();
      }
      const int nucleusNum = this->Pair(PairNum)->GetNucleus();
      std::map<int, int>::iterator lvl = levelNucleus.find(Line.levelID());
      if (lvl == levelNucleus.end()) {
        levelNucleus[Line.levelID()] = nucleusNum;
      } else if (lvl->second != nucleusNum) {
        configure.outStream << "ERROR: level ID " << Line.levelID() << " (J = " << Line.levelJ()
                            << (Line.levelPi() == -1 ? "-" : "+") << ", E = " << Line.levelE()
                            << " MeV) has channels of two compound nuclei, "
                            << this->NucleusLabel(lvl->second) << " and " << this->NucleusLabel(nucleusNum)
                            << ": a level belongs to one nucleus." << std::endl;
        return -1;
      }
      JGroup NewJGroup(Line);
      NewJGroup.SetNucleus(nucleusNum, this->Pair(PairNum)->GetNucleusZ(), this->Pair(PairNum)->GetNucleusA());
      JGroupNum = this->IsJGroup(NewJGroup);
      if (!JGroupNum) {
        this->AddJGroup(NewJGroup);
        JGroupNum = this->IsJGroup(NewJGroup);
      }
      AChannel NewChannel(Line, PairNum);
      ChannelNum = this->Group(JGroupNum)->IsChannel(NewChannel);
      if (!ChannelNum) {
        this->Group(JGroupNum)->AddChannel(NewChannel);
        ChannelNum = this->Group(JGroupNum)->IsChannel(NewChannel);
        if (this->Group(JGroupNum)->GetChannel(ChannelNum)->GetL() > maxLValue &&
            this->Group(JGroupNum)->GetChannel(ChannelNum)->GetRadType() == 'P')
          maxLValue = this->Group(JGroupNum)->GetChannel(ChannelNum)->GetL();

        // Calculate and set Wigner Limit for the newly created channel
        PPair *channelPair = this->Pair(PairNum);
        this->Group(JGroupNum)->GetChannel(ChannelNum)->SetRadius(channelPair->GetPairRadius(), channelPair->GetRedMass());
      }
      // A line belongs to a level of its J-group that carries the same level ID within 1 keV, or failing that
      // to one at the identical energy.  Until 2026-10-24 any level within 1 keV was taken (JGroup::IsLevel),
      // whatever its ID, so a narrow level put next to a broad one of the same J^pi (17O+a: 3/2- at 8.2900 and
      // 8.2893 MeV, IDs 91 and 22) was silently merged into it: its widths overwrote the broad level's and the
      // fit lost three parameters.  Lines with other IDs within 1 keV now make their own level, with a warning.
      // The level of the line: io::LevelGrouper, the rule pyazr uses too (docs/upgrade/project.md, P2).
      ALevel NewLevel(Line);
      JGroup *group = this->Group(JGroupNum);
      const azure2::io::LevelGrouper::Slot slot = levelGrouper.Assign(JGroupNum, NewLevel.FileE(), Line.levelID());
      if (slot.level > group->NumLevels()) group->AddLevel(NewLevel);
      LevelNum = slot.level;
      group->GetLevel(LevelNum)->AddID(Line.levelID());
      if (slot.kind == azure2::io::LevelGrouper::Kind::SameEnergy)
        configure.outStream << "WARNING: level IDs " << slot.otherID << " and " << Line.levelID() << " (J = "
                            << group->GetJ() << (group->GetPi() == -1 ? "-" : "+") << ", E = " << NewLevel.FileE()
                            << " MeV) have the same J^pi and energy; their lines are read as one level." << std::endl;
      if (slot.kind == azure2::io::LevelGrouper::Kind::NewNearby)
        configure.outStream << "WARNING: level ID " << Line.levelID() << " (J = " << group->GetJ()
                            << (group->GetPi() == -1 ? "-" : "+") << ", E = " << NewLevel.FileE()
                            << " MeV) is within 1 keV of level ID " << slot.otherID << " (E = "
                            << group->GetLevel(slot.other)->FileE() << " MeV); they are kept as two levels."
                            << std::endl;
      // The width belongs to the GROUP's channel number: levels of one J-group may list their channel
      // lines in a different order, or not list a channel another level of the group carries.
      this->Group(JGroupNum)->GetLevel(LevelNum)->SetGammaFromLine(ChannelNum, Line);
      channelsListed[std::make_pair(JGroupNum, LevelNum)].push_back(ChannelNum);
    }
  }
  // A malformed row, or no </levels>: the rows before it are read (and refused for their own reasons first).
  if (section.malformed || !section.complete) return -1;

  // Every level of a J-group needs one width per group channel.  Until 2026-10-17 widths were appended in the
  // order of each level's own lines, so a level listing fewer channels than its group left a short width vector
  // (TransformIn wrote past its end: heap corruption, "free(): invalid next size" at exit), and a level listing
  // the same channels in another order had its widths attached to the wrong channels without a message
  // (reported by the A = 4 isospin session).  Widths now go to their group channel number; a channel a level
  // does not list gets a zero, fixed width, with a warning, and a different line order is reported.
  for (const auto &entry : channelsListed) {
    JGroup *group = this->Group(entry.first.first);
    ALevel *level = group->GetLevel(entry.first.second);
    const std::vector<int> &listed = entry.second;
    std::ostringstream what;
    what << "level J = " << group->GetJ() << (group->GetPi() == -1 ? "-" : "+") << ", E = " << level->FileE() << " MeV";
    std::set<int> seen;
    bool reordered = false;
    for (size_t k = 0; k < listed.size(); k++) {
      if (!seen.insert(listed[k]).second)
        configure.outStream << "WARNING: " << what.str() << " lists channel " << listed[k]
                            << " of its J-group twice; the later line's width is used." << std::endl;
      if (k > 0 && listed[k] < listed[k - 1]) reordered = true;
    }
    if (reordered)
      configure.outStream << "WARNING: " << what.str() << " lists its channels in a different order from its "
                          << "J-group; each width is matched to its channel (pair, l, s)." << std::endl;
    for (int ch = 1; ch <= group->NumChannels(); ch++) {
      if (seen.count(ch)) continue;
      AChannel *c = group->GetChannel(ch);
      configure.outStream << "WARNING: " << what.str() << " has no line for channel " << ch << " of its J-group (pair "
                          << this->Pair(c->GetPairNum())->GetPairKey() << ", l = " << c->GetL() << ", s = " << c->GetS()
                          << "); a zero, fixed width is used." << std::endl;
    }
    level->PadChannels(group->NumChannels());
  }

  // Channel-radius overrides (the API's SetRadius).  All of them, every time: before 2026-11-02 only the latest
  // call's pair was passed, so a second SetRadius, a potential change or a mode switch rebuilt the model with the
  // earlier pairs back at their .azr radius, silently.
  for (std::map<int, double>::const_iterator it = radii.begin(); it != radii.end(); ++it) {
    if (it->first < 1 || it->first > this->NumPairs()) {
      configure.outStream << "ERROR: channel-radius override for pair " << it->first << ", but the model has pairs 1-"
                          << this->NumPairs() << "." << std::endl;
      return -1;
    }
    PPair *changedPair = this->Pair(it->first);
    changedPair->SetChRad(it->second);
    // The channels' radii and Wigner limits above were set inside the parse loop, i.e. with
    // the radius as written in the .azr.  Redo them for the overridden pair,
    // or every channel quantity after a radius change is at the old radius.
    for (int j = 1; j <= this->NumJGroups(); j++) {
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        AChannel *theChannel = this->Group(j)->GetChannel(ch);
        if (theChannel->GetPairNum() == it->first)
          theChannel->SetRadius(changedPair->GetPairRadius(), changedPair->GetRedMass());
      }
    }
  }

  for (int p = 1; p <= this->NumPairs(); p++) this->Pair(p)->SetHybridPotential(configure.potentials.PotentialFor(this->Pair(p)->GetPairKey()));  // a value, none: plain Coulomb
  if (configure.includeAllJPi)
    this->AddHardSphereGroups(configure, configure.hardSphereMaxL >= 0 ? configure.hardSphereMaxL : maxLValue);
  if (RefuseForbiddenIdenticalChannels(*this, configure) == -1) return -1;
  if ((configure.paramMask & Config::USE_EXTERNAL_CAPTURE) && this->NumJGroups() > 0 && this->NumPairs() > 0)
    if (this->ParseExternalCapture(configure, ecPairs, section.rows) == -1) return -1;

  if (this->NumNuclei() > 1)
    for (int j = 1; j <= this->NumJGroups(); j++) this->Group(j)->SetNucleusLabelled(true);
  this->BuildSystems(configure);
  std::vector<int> numLevels, numChannels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    numLevels.push_back(this->Group(j)->NumLevels());
    numChannels.push_back(this->Group(j)->NumChannels());
  }
  state = LevelState(numLevels, numChannels);
  parameterLayout_ = std::make_shared<const LevelParameterLayout>(numLevels, numChannels);

  if (std::shared_ptr<const ModelDataSet> read = ModelData::Read(*this, configure, section.rows)) modelData_ = read; else return -1;
  if (std::shared_ptr<const GammaDecaySet> read = GammaDecaySet::Read(*this, configure, section.rows)) gammaDecays_ = read; else return -1;
  // The input values: the <levels> rows', a gamma_decay record's level with the partial widths of the record.
  inputs = LevelInputs(numLevels, numChannels);
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      const ALevel *level = this->Group(j)->GetLevel(la);
      inputs.SetE(j, la, level->FileE());
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        inputs.SetGamma(j, la, ch, level->FileGamma(ch));
        state.SetZeroInputWidth(j, la, ch, level->FileZeroWidth(ch));
      }
    }
  inputs.SetGammaDecayNumbers(gammaDecays_->FileNumbers());
  gammaDecays_->WriteInputWidths(inputs);
  // right after the pairs are known, before any data file is read (a project whose data fail still gets it)
  WarnMassSeparationMismatch(*this, configure);

  return 0;
}

/*!
 * Every J^pi group in the calculation scatters through the hard-sphere phases of its channels whether or
 * not it has levels (the elastic T-matrix diagonal is e^{2i omega} - e^{2i(omega - phi)}(1 + U), and U is
 * the level sum). The <levels> block, though, only creates a group when an active level carries it, so a
 * J^pi with no level was simply absent, and an elastic calculation lost its potential scattering in that
 * partial wave unless the user remembered a zero-width "dummy" level. With includeAllJPi set this adds
 * the missing groups: for every particle pair, channel spin s and l <= maxL, every (J, pi) that has no
 * group yet gets one with the particle channels of all pairs and no levels -- the same channel set the
 * setup utility builds for a level of that J^pi. Capture and beta-decay pairs add nothing here: with no
 * level there is nothing for a photon or a beta to connect.
 *
 * A level switched off in the file still leaves its J^pi scattering, which is what a decomposition into
 * resonant and potential parts needs; the old convention (the whole partial wave disappears with its
 * levels) is what includeAllJPi=0 gives.
 */

int CompoundSystem::AddHardSphereGroups(const Config &configure, int maxL) {
  int added = 0;
  for (int aa = 1; aa <= this->NumPairs(); aa++) {
    PPair *pair = this->Pair(aa);
    if (pair->GetPType() != 0) continue;
    double s1 = pair->GetJ(1), s2 = pair->GetJ(2);
    int sPi = pair->GetPi(1) * pair->GetPi(2);
    for (double chS = fabs(s1 - s2); chS <= s1 + s2 + 1.0e-6; chS += 1.0) {
      for (int chL = 0; chL <= maxL; chL++) {
        int chPi = sPi * (chL % 2 == 0 ? 1 : -1);
        for (double jValue = fabs(chS - chL); jValue <= chS + chL + 1.0e-6; jValue += 1.0) {
          JGroup newJGroup(jValue, chPi);
          newJGroup.SetNucleus(pair->GetNucleus(), pair->GetNucleusZ(), pair->GetNucleusA());
          if (this->IsJGroup(newJGroup)) continue;
          newJGroup.SetInRMatrix(true);
          this->AddJGroup(newJGroup);
          int jGroupNum = this->IsJGroup(newJGroup);
          JGroup *jg = this->Group(jGroupNum);
          for (int ir = 1; ir <= this->NumPairs(); ir++) {
            PPair *chPair = this->Pair(ir);
            if (chPair->GetPType() != 0) continue;
            if (chPair->GetNucleus() != pair->GetNucleus()) continue;   // one nucleus per group
            double t1 = chPair->GetJ(1), t2 = chPair->GetJ(2);
            int tPi = chPair->GetPi(1) * chPair->GetPi(2);
            bool identicalPair = chPair->IsIdentical();
            for (double cS = fabs(t1 - t2); cS <= t1 + t2 + 1.0e-6; cS += 1.0) {
              for (int cL = 0; cL <= maxL; cL++) {
                if (tPi * (cL % 2 == 0 ? 1 : -1) != chPi) continue;
                if (fabs(cS - cL) > jValue + 1.0e-6 || jValue > cS + cL + 1.0e-6) continue;
                // Bose/Fermi symmetry for an identical pair: only l+s even (RefuseForbiddenIdenticalChannels).
                if (identicalPair && (cL + (int)(cS + 0.5)) % 2 != 0) continue;
                AChannel newChannel(cL, cS, ir, 'P');
                if (jg->IsChannel(newChannel)) continue;
                jg->AddChannel(newChannel);
                jg->GetChannel(jg->NumChannels())->SetRadius(chPair->GetPairRadius(), chPair->GetRedMass());
              }
            }
          }
          added++;
          if (configure.screenCheckMask & Config::CHECK_COMPOUND_NUCLEUS)
            configure.outStream << "Hard-sphere group added: J = " << jValue << (chPi == -1 ? "-" : "+")
                                << " with " << jg->NumChannels() << " channels" << std::endl;
        }
      }
    }
  }
  return added;
}

/*!
 * Initializes the compound nucleus object.  This includes calculating the boundary conditions, transforming from
 * physical to formal parameters, creating and sorting all reaction pathways, and calculating angular interference
 * contributions and coefficients.  A CompoundSystem object can only be initialized for use AFTER it is filled.
 */

void CompoundSystem::Initialize(const Config &configure, LevelState &state, LevelInputs &inputs,
                                const CompoundSetUpRequests &requests) {
  WarnForbiddenSpinlessChannels(*this, configure);

  // Calculate Boundary Conditions
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Boundary Conditions..." << std::endl;
  this->CalcBoundaryConditions(configure, inputs);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_BOUNDARY_CONDITIONS)
    WriteBoundaryConditionsCheck(*this, configure);

  // Transform Input Parameters
  if (configure.paramMask & Config::TRANSFORM_PARAMETERS) {
    if (!(configure.paramMask & Config::USE_API))
      configure.outStream << "Performing Input Parameter Transformation..." << std::endl;
    // No standard R-matrix may correspond to the input (TransformIn): stop, as for any other GSL failure here.
    if (!this->TransformIn(configure, state, inputs)) throw GSLException(inputs.TransformError());
  }

  this->SetUpECFinalBlocks(configure, inputs);

  // Sort reaction pathways
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Sorting Reaction Pathways..." << std::endl;
  this->BuildPathways(configure, requests);
  this->SetUpFormalCaptureMasks(configure, state);
  // Everything the last part writes is written; what follows reads it.
  finished_ = true;
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_PATHWAYS)
    WritePathwaysCheck(*this, configure);

  // Calculate Angular Distribution Coefficients
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Angular Distribution Coefficients..." << std::endl;
  WarnTruncatedLegendreOrders(*this, configure);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_ANGULAR_DISTS)
    WriteAngularDistributionsCheck(*this, configure);
}

/*!
 * Adds a particle pair to the PPair vector.
 */

void CompoundSystem::AddPair(PPair pPair) {
  pPair.SetNucleus(this->FindOrAddNucleus(pPair.GetNucleusZ(), pPair.GetNucleusA()));
  pairs_.push_back(pPair);
}

/*!
 * Adds a \f$ J^\pi \f$ group to the JGroup vector.
 */

void CompoundSystem::AddJGroup(JGroup jGroup) {
  jgroups_.push_back(jGroup);
}

/*!
 * Calculates the boundary conditions: B_c = S_c at the lowest R-matrix level of each \f$ J^\pi \f$ group.
 */

void CompoundSystem::CalcBoundaryConditions(const Config &configure, const LevelInputs &inputs) {
  // B_c = S_c(E of the group's lowest R-matrix level) for every particle channel; a non-particle channel copies channel
  // 1's value (it never enters the level matrix); a group without R-matrix levels gets 0, where B_c drops out of the
  // collision matrix (physics/ boundary_values, docs/upgrade/level_parameters.md Eq. 3.1; HANDOFF decision 10).
  const bool formal = !(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  try {
    const SystemOfCompound &structure = this->SystemStructure();
    const azure2::physics::ChannelFunctions cf = EngineChannelFunctions(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      JGroup *theJGroup = this->Group(structure.map.groupOfBlock[b]);
      std::vector<double> energies;
      for (int la : structure.map.levelOfBlockLevel[b])
        energies.push_back(inputs.E(structure.map.groupOfBlock[b], la));
      const std::vector<double> B = azure2::physics::boundary_values(structure.system, b, energies, cf);
      const std::vector<int> &channels = structure.map.channelOfBlockChannel[b];
      for (std::size_t k = 0; k < channels.size(); k++) theJGroup->GetChannel(channels[k])->SetBoundaryCondition(B[k]);
      // Builds before 2026-11-07 took B_c at the first listed level: formal parameters they wrote for a group listed
      // out of energy order (param.sav, param.par, --no-transform input) belong to another B_c.
      if (formal && !energies.empty()) {
        const double lowest = *std::min_element(energies.begin(), energies.end());
        if (lowest < energies[0])
          configure.outStream << "Note: J = " << theJGroup->GetJ() << (theJGroup->GetPi() > 0 ? "+" : "-")
                              << ": the boundary condition is set at the lowest level, E = " << lowest
                              << " MeV, not at the first listed one, E = " << energies[0]
                              << " MeV (builds before 2026-11-07); formal parameters written by those builds for this "
                              << "group (param.sav, param.par, --no-transform input) do not apply." << std::endl;
      }
    }
  } catch (const GSLException &) {
    throw;
  } catch (const std::exception &e) {
    // The System refuses what it cannot represent (a spin that is not a multiple of 1/2, an unknown pair type):
    // stop as for any other failure of the initialization.
    throw GSLException(std::string("Boundary conditions: ") + e.what());
  }
}


/// The external-capture final states of the compound after Fill, in J-group and level order: what physics/
/// reaction_pathways needs of each (2J, parity, photon partition, multipoles, channels) and its J-group and level.
static std::vector<azure2::PathwayFinalState> ECFinalStates(const CompoundSystem &c) {
  std::vector<azure2::PathwayFinalState> out;
  for (int j = 1; j <= c.NumJGroups(); j++) {
    const JGroup *theFinalJGroup = c.GetJGroup(j);
    for (int la = 1; la <= theFinalJGroup->NumLevels(); la++) {
      const ALevel *theFinalLevel = theFinalJGroup->GetLevel(la);
      if (!theFinalLevel->IsECLevel()) continue;
      azure2::PathwayFinalState f;
      f.state.two_J = Twice(theFinalJGroup->GetJ());
      f.state.parity = theFinalJGroup->GetPi();
      f.state.photon_partition = theFinalLevel->GetECPairNum() - 1;
      f.state.multipoles = theFinalLevel->GetECMultMask();
      for (int chp = 1; chp <= theFinalJGroup->NumChannels(); chp++) {
        const AChannel *finalChannel = theFinalJGroup->GetChannel(chp);
        f.state.channels.push_back(azure2::physics::FinalStateChannel{
            finalChannel->GetPairNum() - 1, c.GetPair(finalChannel->GetPairNum())->GetPType() == 0, finalChannel->GetL(),
            Twice(finalChannel->GetS())});
      }
      f.group = j;
      f.level = la;
      out.push_back(f);
    }
  }
  return out;
}


const azure2::PathwayTable &CompoundSystem::Pathways() const {
  if (!pathways_) throw std::logic_error("CompoundSystem::Pathways: the compound nucleus has not been initialized");
  return *pathways_;
}

void CompoundSystem::RefuseUnfinished(const char *accessor) const {
  if (!finished_)
    throw std::logic_error(std::string(accessor) + ": read before the compound system's last part wrote it "
                                                   "(only SetUpCompoundSystem finishes a compound system)");
}

bool CompoundSystem::HasFormalECFinalBlocks() const {
  RefuseUnfinished("CompoundSystem::HasFormalECFinalBlocks");
  return ecFinalBlocks_ != nullptr;
}

/*!
 * Builds the pathway table from the System's structure, the external-capture final states and the configuration:
 * R-matrix capture, the largest Legendre order, the entrance pairs of the data's requests (flagged on their PPair
 * here) and their capture-A_y pairs.  Run by Initialize after the input transformation; nothing writes the table
 * afterwards.
 */
void CompoundSystem::BuildPathways(const Config &configure, const CompoundSetUpRequests &requests) {
  for (int pair : requests.entrancePairs) pairs_[pair - 1].SetEntrance();
  azure2::PathwayTableOptions options;
  options.rmatrix_capture = (configure.paramMask & Config::USE_RMC_FORMALISM) != 0;
  options.max_l_order = configure.maxLOrder;
  for (int aa = 1; aa <= this->NumPairs(); aa++)
    if (pairs_[aa - 1].IsEntrance()) options.entrances.push_back(aa);
  options.capture_ay = requests.captureAyPairs;
  pathways_ = azure2::BuildPathwayTable(this->SystemStructure(),
                                        ECFinalStates(*this), options);
}
