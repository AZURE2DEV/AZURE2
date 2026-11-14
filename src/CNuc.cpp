#include <cmath>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include "CNuc.h"
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
#include "CNucLevelParameters.h"
#include "CNucSystem.h"
#include "CNucChecks.h"
#include "CNucPathways.h"
#include "CNucReport.h"
#include "CNucTransformMessages.h"
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
#include "CNucExternalCapture.h"
#include "ShftFunc.h"
#include "azure2/io/LevelGrouping.h"
#include "azure2/physics/ObservableCoupling.h"
#include "azure2/physics/Pathways.h"

namespace {

//! 2j as an integer, for the physics/ coupling coefficients (docs/upgrade/cnuc_split.md section 9).
int Twice(double j) { return (int)std::lround(2.0 * j); }

}  // namespace

/*!
 * Returns true if a specified pair key exists in the PPair vector, otherwise returns false.
 */

bool CNuc::IsPairKey(int key) const {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumPairs()) {
    if (key == this->GetPair(c + 1)->GetPairKey()) b = true;
    c++;
  }
  return b;
}

std::vector<int> CNuc::CapturePairKeys(int entranceKey) const {
  std::vector<int> keys;
  const PPair *entrance = GetPairFromKey(entranceKey);
  if (!entrance) return keys;
  for (int i = 1; i <= NumPairs(); i++)
    if (GetPair(i)->GetPType() == 10 && GetPair(i)->GetNucleus() == entrance->GetNucleus())
      keys.push_back(GetPair(i)->GetPairKey());
  return keys;
}

std::string CNuc::RatePairsProblem(int entranceKey, int exitKey) const {
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

int CNuc::NumPairs() const {
  return pairs_.size();
}

/*!
 * Returns the number of \f$ J^\pi \f$ groups in the JGroup vector.
 */

int CNuc::NumJGroups() const {
  return jgroups_.size();
}

/*!
 * Tests if a particle pair exists in the PPair vector. If pair exists, the position in
 * the vector is returned.  Otherwise, the function returns 0.
 */

int CNuc::IsPair(PPair pair) {
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

/*!
 * Tests if a \f$ J^\pi \f$ group exists in the JGroup vector.  If the group exists, the position in the vector
 * is returned.  Otherwise, the function returns 0.
 */

int CNuc::IsJGroup(JGroup jGroup) {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumJGroups()) {
    // A group belongs to one compound nucleus: the same J^pi of another
    // nucleus in the file is a different group, with its own level matrix.
    if (jGroup.GetJ() == this->Group(c + 1)->GetJ() &&
        jGroup.GetPi() == this->Group(c + 1)->GetPi() &&
        jGroup.GetNucleus() == this->Group(c + 1)->GetNucleus()) b = true;
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

int CNuc::NumNuclei() const {
  return nuclei_.size();
}

const CompoundNucleus *CNuc::GetNucleus(int nucleusNum) const {
  return &nuclei_[nucleusNum - 1];
}

std::string CNuc::NucleusLabel(int nucleusNum) const {
  std::ostringstream out;
  out << "Z=" << nuclei_[nucleusNum - 1].Z << " A=" << nuclei_[nucleusNum - 1].A;
  return out.str();
}

int CNuc::FindOrAddNucleus(int z, int a) {
  for (unsigned int i = 0; i < nuclei_.size(); i++)
    if (nuclei_[i].Z == z && nuclei_[i].A == a) return i + 1;
  CompoundNucleus nucleus;
  nucleus.Z = z;
  nucleus.A = a;
  nuclei_.push_back(nucleus);
  return nuclei_.size();
}

/*!
 * Returns the position of a particle pair in the PPair vector based on the pair key.
 * Pair keys are how particle pairs are specified in
 * the setup files, but may not correspond to the position of the particle pair in the PPair vector.
 * If the pair exists, the position in the vector is returned.  Otherwise, the function returns 0.
 */

int CNuc::GetPairNumFromKey(int key) const {
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
 * Fills the compound nucleus object, and all nested objects, from data specified in the nuclear and external capture
 * input files.  Returns -1 if the files could not be read, and 0 if the files were read successfully.
 */

int CNuc::Fill(const Config &configure, const std::map<int, double> &radii) {
  levelInputs_ = LevelInputs();
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
        this->Group(JGroupNum)->GetChannel(ChannelNum)->SetWignerLimit(channelPair->GetRedMass(), channelPair->GetChRad());
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
    // The Wigner limits above were computed inside the parse loop, i.e. with
    // the radius as written in the .azr.  Redo them for the overridden pair,
    // or every theta^2 reported after a radius change is against the old limit.
    for (int j = 1; j <= this->NumJGroups(); j++) {
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        AChannel *theChannel = this->Group(j)->GetChannel(ch);
        if (theChannel->GetPairNum() == it->first)
          theChannel->SetWignerLimit(changedPair->GetRedMass(),
                                     changedPair->GetChRad());
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
  levelState_ = LevelState(numLevels, numChannels);
  parameterLayout_ = std::make_shared<const LevelParameterLayout>(numLevels, numChannels);

  if (std::shared_ptr<const ModelDataSet> read = ModelData::Read(*this, configure, section.rows)) modelData_ = read; else return -1;
  if (std::shared_ptr<const GammaDecaySet> read = GammaDecaySet::Read(*this, configure, section.rows)) gammaDecays_ = read; else return -1;
  // The input values: the <levels> rows', a gamma_decay record's level with the partial widths of the record.
  levelInputs_ = LevelInputs(numLevels, numChannels);
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      const ALevel *level = this->Group(j)->GetLevel(la);
      levelInputs_.SetE(j, la, level->FileE());
      levelState_.SetInputEnergy(j, la, level->FileE());
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        levelInputs_.SetGamma(j, la, ch, level->FileGamma(ch));
        levelState_.SetZeroInputWidth(j, la, ch, level->FileZeroWidth(ch));
      }
    }
  levelInputs_.SetGammaDecayNumbers(gammaDecays_->FileNumbers());
  gammaDecays_->WriteInputWidths(levelInputs_);
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

int CNuc::AddHardSphereGroups(const Config &configure, int maxL) {
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
                jg->GetChannel(jg->NumChannels())->SetWignerLimit(chPair->GetRedMass(), chPair->GetChRad());
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
 * Marks the final state of each photon pair with external capture (`ecPairs`: pair number -> multipole mask) and
 * stores its N_f integrals and Whittaker conversion factors; the final states are decided by
 * experiment::ExternalCaptureFinalStates on the <levels> `rows`, as the validation decides them.
 */

int CNuc::ParseExternalCapture(const Config &configure, std::map<int, int> &ecPairs,
                               const std::vector<azure2::io::LevelRow> &rows) {
  std::vector<int> ecKeys;
  for (const auto &ec : ecPairs) ecKeys.push_back(this->Pair(ec.first)->GetPairKey());
  const azure2::experiment::LevelTables tables = azure2::experiment::LevelTablesOfRows(rows);
  const std::vector<azure2::experiment::FinalStateOfPair> finals = azure2::experiment::ExternalCaptureFinalStates(tables, ecKeys);
  std::set<int> usedExitKeys = ProjectExitPairKeys(configure);  // the exit pairs of the active data sets (their specs)
  // A total-capture rate (exit key -1) exits into every capture pair of its entrance nucleus.
  if ((configure.paramMask & Config::CALCULATE_REACTION_RATE) && configure.rateParams.exitPair == -1)
    for (int key : CapturePairKeys(configure.rateParams.entrancePair)) usedExitKeys.insert(key);
  for (std::map<int, int>::iterator ec = ecPairs.begin(); ec != ecPairs.end(); ec++) {
    PPair *exitPair = this->Pair(ec->first);
    if (exitPair->GetPType() != 10) {
      configure.outStream << "Final state is not a capture pair." << std::endl;
      continue;
    }
    // The final state must be a level of the project (user guide, "Gamma-Ray Calculations"; it carries the ANC and the
    // gamma widths): refused for a pair an active data set (or the rate) exits into, else warned about and switched off,
    // so no level is ever made up.  tests/ec_final_state_listed.
    double jValue = exitPair->GetJ(2);
    int parity = exitPair->GetPi(2);
    JGroup newJGroup(jValue, parity);
    newJGroup.SetNucleus(exitPair->GetNucleus(), exitPair->GetNucleusZ(), exitPair->GetNucleusA());
    int jGroupNum = this->IsJGroup(newJGroup);
    const azure2::experiment::FinalStateOfPair &finalState = *azure2::experiment::FinalStateOf(finals, exitPair->GetPairKey());
    if (finalState.levels.empty()) {
      const bool used = usedExitKeys.count(exitPair->GetPairKey()) > 0;
      configure.outStream << (used ? "ERROR" : "WARNING") << ": external capture into pair " << exitPair->GetPairKey()
                          << " needs its final state, J^pi = " << jValue << (parity > 0 ? "+" : "-") << " at E_x = "
                          << exitPair->GetHeavyExE() << " MeV, as a level of the project (with its ANC and gamma widths); none is listed"
                          << (jGroupNum ? " within 1 keV in that J^pi group." : " and that J^pi group does not exist.")
                          << (used ? "" : "  No segment uses this pair: external capture into it is switched off for this run.")
                          << std::endl;
      if (used) return -1;
      continue;
    }
    // Two levels of one J^pi within 1 keV of each other (different level IDs): the final state is ambiguous.
    if (finalState.levels.size() > 1) {
      configure.outStream << "ERROR: external capture into pair " << exitPair->GetPairKey() << ": " << finalState.levels.size()
                          << " levels of J^pi = " << jValue << (parity > 0 ? "+" : "-") << " lie within 1 keV of its final-state energy E_x = "
                          << exitPair->GetHeavyExE() << " MeV, so the final state is ambiguous.  Move the extra level(s) by more than 1 keV"
                          << " or give the pair the energy of the one it decays to." << std::endl;
      return -1;
    }
    const azure2::experiment::LevelInfo &level = tables.levels[finalState.levels[0]];
    if (level.group != jGroupNum) throw std::logic_error("ParseExternalCapture: the level tables number J^pi groups unlike Fill");
    const int levelNum = level.level;
    if (finalState.sharedWithPairKey != 0) {
      FinalStateOfAnotherPair(*this, configure, finalState.sharedWithPairKey, levelNum, ec->first);
      return -1;
    }
    this->Group(jGroupNum)->GetLevel(levelNum)->SetECParams(ec->first, ec->second);
    for (int ch = 1; ch <= this->Group(jGroupNum)->NumChannels(); ch++) {
      PPair *theFinalPair = this->Pair(this->Group(jGroupNum)->GetChannel(ch)->GetPairNum());
      double nfIntegralValue = 0.;
      double ecConvert = 0.;
      if (theFinalPair->GetPType() == 0) {
        NFIntegral newNFIntegral(theFinalPair);
        nfIntegralValue = newNFIntegral(this->Group(jGroupNum)->GetChannel(ch)->GetL(), theFinalPair->SpinOrbitTwoJ(jValue), exitPair->GetHeavyExE());
        WhitFunc newWhitFunc(theFinalPair);
        double whitConv = newWhitFunc(this->Group(jGroupNum)->GetChannel(ch)->GetL(), theFinalPair->SpinOrbitTwoJ(jValue),
                                      theFinalPair->GetChRad(),
                                      fabs(exitPair->GetHeavyExE() - theFinalPair->GetSepE() - theFinalPair->GetExE()));
        ecConvert = sqrt(2.0 * theFinalPair->GetRedMass() * theFinalPair->GetChRad() * uconv / pow(hbarc, 2.0)) / whitConv;
      }
      this->Group(jGroupNum)->GetLevel(levelNum)->AddNFIntegral(nfIntegralValue);
      this->Group(jGroupNum)->GetLevel(levelNum)->AddECConversionFactor(ecConvert);
    }
  }
  return 0;
}

/*!
 * Initializes the compound nucleus object.  This includes calculating the boundary conditions, transforming from
 * physical to formal parameters, creating and sorting all reaction pathways, and calculating angular interference
 * contributions and coefficients.  A CNuc object can only be initialized for use AFTER it is filled.
 */

void CNuc::Initialize(const Config &configure) {
  WarnForbiddenSpinlessChannels(*this, configure);

  // Calculate Boundary Conditions
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Boundary Conditions..." << std::endl;
  this->CalcBoundaryConditions(configure);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_BOUNDARY_CONDITIONS)
    WriteBoundaryConditionsCheck(*this, configure);

  // Transform Input Parameters
  if (configure.paramMask & Config::TRANSFORM_PARAMETERS) {
    if (!(configure.paramMask & Config::USE_API))
      configure.outStream << "Performing Input Parameter Transformation..." << std::endl;
    // No standard R-matrix may correspond to the input (TransformIn): stop, as for any other GSL failure here.
    if (!this->TransformIn(configure)) throw GSLException(levelInputs_.TransformError());
  }

  this->SetUpECFinalBlocks(configure);

  // Sort reaction pathways
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Sorting Reaction Pathways..." << std::endl;
  this->BuildPathways(configure);
  this->SetUpFormalCaptureMasks(configure);
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

void CNuc::AddPair(PPair pPair) {
  pPair.SetNucleus(this->FindOrAddNucleus(pPair.GetNucleusZ(), pPair.GetNucleusA()));
  pairs_.push_back(pPair);
}

/*!
 * Adds a \f$ J^\pi \f$ group to the JGroup vector.
 */

void CNuc::AddJGroup(JGroup jGroup) {
  jgroups_.push_back(jGroup);
}

namespace phys = azure2::physics;

/*!
 * Performs the initial parameter transformations from physical to formal parameters: the observed widths, ANCs and
 * amplitudes of the levels become Brune's (or Park's) reduced width amplitudes, the external-capture amplitude is
 * subtracted from the photon channels, and, with the Brune parametrization off, the levels become formal ones.  The
 * arithmetic is physics/ level_parameters' (docs/upgrade/level_parameters.md sections 2-4), with the engine's channel
 * functions (ChannelFunctionsFromCNuc); the formal parameters take B_c at the lowest level of a group, and each formal
 * level the orientation of its Brune level (HANDOFF decision 10, cnuc_split.md section 7).  This function words the
 * warnings and keeps the bound-state normalization of external capture.
 */

bool CNuc::TransformIn(const Config &configure, LevelState &state, LevelInputs &inputs) const {
  inputs.StartTransformation();
  // Park (2021): the fit parameter is the observed reduced width amplitude, Gamma = 2 P gamma^2 (Eqs. 31-33), so the
  // factor that turns an observed width into Brune's reduced width is not applied; its J (Eq. 28) must stay positive.
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool formal = !(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  try {
    const CNucSystem &structure = this->SystemStructure();
    const phys::System &system = structure.system;
    const CNucSystemMap &map = structure.map;
    const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, map, configure);

    // 1. Observed widths, ANCs and amplitudes -> Brune or Park amplitudes (Eqs. 2.1-2.3).
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      const phys::ProjectLevels in = ProjectLevelsFromCNuc(*this, inputs, map, b);
      const phys::BasisAmplitudes amp = phys::amplitudes_from_observed(
          system, b, phys::internal_channels(system, b), in, park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        const ALevel *theLevel = theJGroup->GetLevel(la);
        for (std::size_t k = 0; k < nc; k++) {
          const int ch = (int)k + 1;
          const double input = in.value[lb * nc + k];
          // The --ignore-externals mask is a property of the input (external_capture.md section 3).
          state.SetZeroInputWidth(j, la, ch, input == 0.0);
          WarnInputColumn(*this, inputs, configure, j, la, ch, amp.notes[lb * nc + k], input);
        }
        WarnObservedOutOfReach(*this, inputs, configure, j, la, amp, in, lb, nc);
        const double parkNorm = amp.park_norm[lb];
        // The bound-state normalization of an external-capture final state is defined with Brune's reduced widths
        // gamma_obs^2 / J, so with Park's amplitudes it becomes 1 / sqrt(J + sum_c A_c gamma_obs^2) (see
        // CalcShiftFunctions).
        if (park) state.SetParkNorm(j, la, parkNorm);
        const std::vector<double> g(amp.gamma.begin() + lb * nc, amp.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, la, (int)k + 1, g[k]);
        state.SetSqrtNF(j, la, SqrtNFFactor(*this, theJGroup, theLevel, park ? parkNorm : 1.0, g));
      }
    }

    // 2. External capture subtracted from the photon amplitudes (Eq. 3.2) of every group, while every level still
    // holds its Brune (Park) amplitudes: the final state enters external capture with its own pole's amplitude in
    // every basis (cnuc_split.md section 8).  The subtraction changes photon amplitudes only, which no other
    // level's external capture reads.
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(inputs.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(inputs.Gamma(j, la, (int)k + 1));
      }
      const std::vector<phys::ExternalWidth> ext = ExternalWidthsFromCNuc(*this, state, inputs, map, b, configure);
      std::vector<std::string> warnings;
      std::vector<std::size_t> exceeded;
      phys::subtract_external_capture(system, b, gamma, ext, warnings, &exceeded);
      for (std::size_t i : exceeded) {
        WarnExternalExceedsTotal(*this, inputs, configure, j, levels[i / nc], (int)(i % nc) + 1);
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, levels[lb], (int)k + 1, gamma[lb * nc + k]);
    }
    // 3. Without the Brune parametrization, the formal parameters (Eq. 4.1), group by group.
    for (std::size_t b = 0; b < system.num_blocks() && formal; b++) {
      const int j = map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(inputs.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(inputs.Gamma(j, la, (int)k + 1));
      }
      // A hard-sphere group (includeAllJPi) has channels but no levels: nothing to transform.
      if (!E.empty()) {
        std::vector<double> B(nc);
        for (std::size_t k = 0; k < nc; k++) B[k] = theJGroup->GetChannel((int)k + 1)->GetBoundaryCondition();
        phys::LevelSet f;
        try {
          // Each level keeps its position and receives the eigenvalue of the rank of its Brune energy, so a level
          // listed first (a background pole, say) keeps its own energy, widths and fix flags.  TransformOut uses the
          // same rank mapping in the other direction.
          f = phys::formal_from_brune(system, b, E, gamma, B, cf);
        } catch (const std::domain_error &) {
          // N a = E M a has a real solution only for a positive-definite M (Brune, PRC 66, 044611 (2002)); M is
          // 1 - sum_c gamma gamma' dS_c/dE, so widths too large for the channel radii make it indefinite.
          std::ostringstream msg;
          msg << "No standard R-matrix (constant boundary conditions) corresponds to the levels of J^pi = "
              << theJGroup->GetJ() << (theJGroup->GetPi() > 0 ? "+" : "-")
              << ": the observed widths are too large for the channel radii (the matrix 1 - sum gamma gamma' dS/dE "
              << "of the transformation is not positive definite).  Use the Brune parametrization (the default; "
              << "drop --no-brune) or smaller widths for these levels:";
          for (int la : levels)
            msg << std::endl
                << "    " << AZURELabel::Level(theJGroup, inputs.E(j, la), j, la);
          inputs.SetTransformError(msg.str());
          return false;
        }
        E.swap(f.E);
        gamma.swap(f.gamma);
        for (std::size_t lb = 0; lb < levels.size(); lb++) {
          inputs.SetE(j, levels[lb], E[lb]);
          state.SetInputEnergy(j, levels[lb], E[lb]);
        }
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) inputs.SetGamma(j, levels[lb], (int)k + 1, gamma[lb * nc + k]);
    }
  } catch (const std::exception &e) {
    inputs.SetTransformError(std::string("Parameter transformation: ") + e.what());
    return false;
  }
  return true;
}

/*!
 * Calculates the boundary conditions: B_c = S_c at the lowest R-matrix level of each \f$ J^\pi \f$ group.
 */

void CNuc::CalcBoundaryConditions(const Config &configure) {
  // B_c = S_c(E of the group's lowest R-matrix level) for every particle channel; a non-particle channel copies channel
  // 1's value (it never enters the level matrix); a group without R-matrix levels gets 0, where B_c drops out of the
  // collision matrix (physics/ boundary_values, docs/upgrade/level_parameters.md Eq. 3.1; HANDOFF decision 10).
  const bool formal = !(configure.paramMask & Config::USE_BRUNE_FORMALISM);
  try {
    const CNucSystem &structure = this->SystemStructure();
    const azure2::physics::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      JGroup *theJGroup = this->Group(structure.map.groupOfBlock[b]);
      std::vector<double> energies;
      for (int la : structure.map.levelOfBlockLevel[b])
        energies.push_back(levelInputs_.E(structure.map.groupOfBlock[b], la));
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

/*!
 * Performs the final parameter transformations from formal to physical parameters.
 */

ObservedLevelTable CNuc::TransformOut(const Config &configure, LevelState &state, const LevelInputs &inputs) const {
  ObservedLevelTable table = this->EmptyObservedTable();
  this->TransformToBrune(configure, Messages::kPrint, state, inputs, table);
  this->ObservedWidths(configure, state, table);
  return table;
}

/*!
 * The Brune (Park) values of the fit parameters, E and gamma, with N_f^-1/2 formed from them, in `table`: in the
 * formal basis by physics/ brune_from_formal, otherwise a copy (N_f^-1/2 the fit's).  kPrint: print the warnings
 * for levels that do not convert.
 */

void CNuc::TransformToBrune(const Config &configure, Messages messages, LevelState &state, const LevelInputs &inputs,
                            ObservedLevelTable &table) const {
  // Park (2021): the fit parameters already are the observed reduced width
  // amplitudes, so Gamma = 2 P gamma^2 with no 1 + sum gamma^2 dS/dE factor.
  // Refresh J and the bound-state normalization for the current parameters.
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  if (park) {
    this->CalcShiftFunctions(configure, state);
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
        const ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
        if (messages == Messages::kPrint && theLevel->IsInRMatrix() && state.ParkNorm(j, la) <= 0.0) WarnParkNormNotPositive(*this, state, inputs, configure, j, la);
      }
  }
  // Every level first as a copy of its fit values; the formal basis then converts the levels of the R-matrix.
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      const ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
      ObservedLevel &out = table.Level(j, la);
      out.E = state.E(j, la);
      out.iterations = 0;
      out.sqrtNF = state.SqrtNF(j, la);
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) out.gamma[ch - 1] = state.Gamma(j, la, ch);
    }
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    // Formal -> Brune by the fixed-point iteration of physics/ brune_from_formal (level_parameters.md section 4), on
    // the fit values.  A level not in the R-matrix keeps its values.
    const phys::FormalToBruneOptions options;  // the engine's 1e-6 MeV and 1000 iterations
    const CNucSystem &structure = this->SystemStructure();
    const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      const int j = structure.map.groupOfBlock[b];
      const JGroup *theJGroup = this->GetJGroup(j);
      const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma, B(nc);
      for (int la : levels) {
        E.push_back(state.E(j, la));
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(state.Gamma(j, la, (int)k + 1));
      }
      for (std::size_t k = 0; k < nc; k++) B[k] = theJGroup->GetChannel((int)k + 1)->GetBoundaryCondition();
      const phys::LevelSet converted = phys::brune_from_formal(structure.system, b, E, gamma, B, cf, options);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        if (messages == Messages::kPrint) WarnFormalToBrune(*this, inputs, configure, j, la, converted, lb, nc, options);
        ObservedLevel &out = table.Level(j, la);
        out.E = converted.E[lb];
        out.iterations = converted.iterations[lb];
        const std::vector<double> g(converted.gamma.begin() + lb * nc, converted.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) out.gamma[k] = g[k];
        out.sqrtNF = SqrtNFFactor(*this, theJGroup, theJGroup->GetLevel(la), 1.0, g);
      }
    }
  }
}

/*!
 * The observed quantities of every level (cnuc_split.md section 8) from its Brune (Park) values in `table`
 * (TransformToBrune) and, with external capture, the external part of its photon amplitudes: physics/
 * observed_levels.  Sets the external amplitude, W and the reported value with its kind in `table`.
 */

void CNuc::ObservedWidths(const Config &configure, const LevelState &state, ObservedLevelTable &table) const {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool external = (configure.paramMask & Config::USE_EXTERNAL_CAPTURE);
  const bool ignoreZero = (configure.paramMask & Config::IGNORE_ZERO_WIDTHS);
  const CNucSystem &all = this->SystemAllGroups();
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, all.map, configure);
  for (std::size_t b = 0; b < all.system.num_blocks(); b++) {
    const int j = all.map.groupOfBlock[b];
    const JGroup *theJGroup = this->GetJGroup(j);
    const std::vector<int> &levels = all.map.levelOfBlockLevel[b];
    const std::size_t nc = theJGroup->NumChannels();
    std::vector<double> E, gamma;
    std::vector<unsigned char> rwa;
    std::vector<phys::ExternalWidth> ext(levels.size() * nc);
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      const ALevel *theLevel = theJGroup->GetLevel(levels[lb]);
      ObservedLevel &out = table.Level(j, levels[lb]);
      E.push_back(out.E);
      for (std::size_t k = 0; k < nc; k++) {
        const int ch = (int)k + 1;
        const AChannel *theChannel = theJGroup->GetChannel(ch);
        gamma.push_back(out.gamma[k]);
        rwa.push_back(theLevel->GammaIsRWA(ch) ? 1 : 0);
        complex externalWidth(0.0, 0.0);
        if ((theChannel->GetRadType() == 'M' || theChannel->GetRadType() == 'E') && theLevel->IsInRMatrix() &&
            external && !(ignoreZero && state.ZeroInputWidth(j, levels[lb], ch))) {
          externalWidth = ExternalReducedWidthOfTable(*this, table, j, levels[lb], theChannel, configure);
          ext[lb * nc + k].present = true;
          ext[lb * nc + k].value = externalWidth;
        }
        out.external[k] = externalWidth;
      }
    }
    const phys::ObservedLevels obs =
        phys::observed_levels(all.system, b, phys::internal_channels(all.system, b), E, gamma, ext, rwa,
                              park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      for (std::size_t k = 0; k < nc; k++) {
        ObservedLevel &out = table.Level(j, levels[lb]);
        out.W[k] = obs.W[lb * nc + k];
        out.value[k] = obs.value[lb * nc + k];
        out.kind[k] = (unsigned char)obs.kind[lb * nc + k];
      }
  }
}

/*!
 * The observed widths of a parameter vector, as parameters.out reports them (TransformOut without its messages):
 * particle and photon widths with the level's normalization and, with external capture, the external part of the
 * photon amplitude (physics/ observed_levels).  The fit state is set back to the vector afterwards, since
 * TransformToBrune in Park's basis recomputes the shift functions and the bound-state normalizations.
 */

std::optional<ObservedLevelTable> CNuc::ObservedWidthsAt(const Config &configure, const vector_r &params, LevelState &state,
                                                       const LevelInputs &inputs) const {
  this->FillCompoundFromParams(params, state);
  ObservedLevelTable table = this->EmptyObservedTable();
  try {
    this->TransformToBrune(configure, Messages::kQuiet, state, inputs, table);
    this->ObservedWidths(configure, state, table);
  } catch (const std::exception &) {
    this->FillCompoundFromParams(params, state);
    return std::nullopt;
  }
  this->FillCompoundFromParams(params, state);
  return table;
}

/*!
 * Returns a pointer to the particle pair specified by a position in the PPair vector.
 */

PPair *CNuc::Pair(int pairNum) {
  return &pairs_[pairNum - 1];
}

const PPair *CNuc::GetPair(int pairNum) const {
  return &pairs_[pairNum - 1];
}

const PPair *CNuc::GetPairFromKey(int key) const {
  const int num = GetPairNumFromKey(key);
  return num > 0 ? GetPair(num) : nullptr;
}

/*!
 * Returns a pointer to the \f$ J^\pi \f$ group specified by a position in the JGroup vector.
 */

JGroup *CNuc::Group(int jGroupNum) {
  return &jgroups_[jGroupNum - 1];
}

const JGroup *CNuc::GetJGroup(int jGroupNum) const {
  return &jgroups_[jGroupNum - 1];
}

/*!
 * Creates a new copy of the CNuc object in memory and returns a pointer to the new object.
 * Used in AZURECalc function class for thread safety.
 */

CNuc *CNuc::Clone() const {
  CNuc *localCompound = new CNuc(*this);
  return localCompound;
}
