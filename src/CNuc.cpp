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
  transformedIn_ = false;
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
      const azure2::io::LevelGrouper::Slot slot = levelGrouper.Assign(JGroupNum, NewLevel.GetE(), Line.levelID());
      if (slot.level > group->NumLevels()) group->AddLevel(NewLevel);
      LevelNum = slot.level;
      group->GetLevel(LevelNum)->AddID(Line.levelID());
      if (slot.kind == azure2::io::LevelGrouper::Kind::SameEnergy)
        configure.outStream << "WARNING: level IDs " << slot.otherID << " and " << Line.levelID() << " (J = "
                            << group->GetJ() << (group->GetPi() == -1 ? "-" : "+") << ", E = " << NewLevel.GetE()
                            << " MeV) have the same J^pi and energy; their lines are read as one level." << std::endl;
      if (slot.kind == azure2::io::LevelGrouper::Kind::NewNearby)
        configure.outStream << "WARNING: level ID " << Line.levelID() << " (J = " << group->GetJ()
                            << (group->GetPi() == -1 ? "-" : "+") << ", E = " << NewLevel.GetE()
                            << " MeV) is within 1 keV of level ID " << slot.otherID << " (E = "
                            << group->GetLevel(slot.other)->GetE() << " MeV); they are kept as two levels."
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
    what << "level J = " << group->GetJ() << (group->GetPi() == -1 ? "-" : "+") << ", E = " << level->GetE() << " MeV";
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

  if (std::shared_ptr<const ModelDataSet> read = ModelData::Read(*this, configure, section.rows)) modelData_ = read; else return -1;
  if (std::shared_ptr<const GammaDecaySet> read = GammaDecaySet::Read(*this, configure, section.rows)) gammaDecays_ = read; else return -1;
  gammaDecayInputs_ = gammaDecays_->FileNumbers();
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
    if (!this->TransformIn(configure)) throw GSLException(this->TransformError());
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

/*!
 * Fills the CNuc object from the parameter array.
 */

void CNuc::FillCompoundFromParamsPhysical(const vector_r &observed) {
  const vector_r &p = gammaDecays_->Empty() ? observed : gammaDecays_->WithPartialWidths(observed);
  if (!gammaDecays_->Empty()) gammaDecayInputs_ = gammaDecays_->NumbersAt(observed);
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      ALevel *level = this->Group(j)->GetLevel(la);
      level->SetE(p[i]);
      i++;
      std::vector<double> g;
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        level->SetGamma(ch, p[i]);
        g.push_back(p[i]);
        i++;
      }
      level->SetSqrtNFFactor(SqrtNFFactor(*this, this->Group(j), level, 1.0, g));
    }
  }
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

bool CNuc::TransformIn(const Config &configure) {
  transformedIn_ = true;
  transformError_.clear();
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
      JGroup *theJGroup = this->Group(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      const phys::ProjectLevels in = ProjectLevelsFromCNuc(*this, map, b);
      const phys::BasisAmplitudes amp = phys::amplitudes_from_observed(
          system, b, phys::internal_channels(system, b), in, park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        ALevel *theLevel = theJGroup->GetLevel(la);
        for (std::size_t k = 0; k < nc; k++) {
          const int ch = (int)k + 1;
          const double input = in.value[lb * nc + k];
          // The --ignore-externals mask is a property of the input (external_capture.md section 3).
          theLevel->SetZeroInputWidth(ch, input == 0.0);
          WarnInputColumn(*this, configure, j, la, ch, amp.notes[lb * nc + k], input);
        }
        WarnObservedOutOfReach(*this, configure, j, la, amp, in, lb, nc);
        const double parkNorm = amp.park_norm[lb];
        // The bound-state normalization of an external-capture final state is defined with Brune's reduced widths
        // gamma_obs^2 / J, so with Park's amplitudes it becomes 1 / sqrt(J + sum_c A_c gamma_obs^2) (see
        // CalcShiftFunctions).
        if (park) theLevel->SetParkNorm(parkNorm);
        const std::vector<double> g(amp.gamma.begin() + lb * nc, amp.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) theLevel->SetGamma((int)k + 1, g[k]);
        theLevel->SetSqrtNFFactor(SqrtNFFactor(*this, theJGroup, theLevel, park ? parkNorm : 1.0, g));
      }
    }

    // 2. External capture subtracted from the photon amplitudes (Eq. 3.2) of every group, while every level still
    // holds its Brune (Park) amplitudes: the final state enters external capture with its own pole's amplitude in
    // every basis (cnuc_split.md section 8).  The subtraction changes photon amplitudes only, which no other
    // level's external capture reads.
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      JGroup *theJGroup = this->Group(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(theJGroup->GetLevel(la)->GetE());
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(theJGroup->GetLevel(la)->GetGamma((int)k + 1));
      }
      const std::vector<phys::ExternalWidth> ext = ExternalWidthsFromCNuc(*this, map, b, configure);
      std::vector<std::string> warnings;
      std::vector<std::size_t> exceeded;
      phys::subtract_external_capture(system, b, gamma, ext, warnings, &exceeded);
      for (std::size_t i : exceeded) {
        WarnExternalExceedsTotal(*this, configure, j, levels[i / nc], (int)(i % nc) + 1);
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) theJGroup->GetLevel(levels[lb])->SetGamma((int)k + 1, gamma[lb * nc + k]);
    }
    // 3. Without the Brune parametrization, the formal parameters (Eq. 4.1), group by group.
    for (std::size_t b = 0; b < system.num_blocks() && formal; b++) {
      const int j = map.groupOfBlock[b];
      JGroup *theJGroup = this->Group(j);
      const std::vector<int> &levels = map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma;
      for (int la : levels) {
        E.push_back(theJGroup->GetLevel(la)->GetE());
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(theJGroup->GetLevel(la)->GetGamma((int)k + 1));
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
                << "    " << AZURELabel::Level(theJGroup, theJGroup->GetLevel(la), j, la);
          transformError_ = msg.str();
          return false;
        }
        E.swap(f.E);
        gamma.swap(f.gamma);
        for (std::size_t lb = 0; lb < levels.size(); lb++) theJGroup->GetLevel(levels[lb])->SetE(E[lb]);
      }
      for (std::size_t lb = 0; lb < levels.size(); lb++)
        for (std::size_t k = 0; k < nc; k++) theJGroup->GetLevel(levels[lb])->SetGamma((int)k + 1, gamma[lb * nc + k]);
    }
  } catch (const std::exception &e) {
    transformError_ = std::string("Parameter transformation: ") + e.what();
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
      for (int la : structure.map.levelOfBlockLevel[b]) energies.push_back(theJGroup->GetLevel(la)->GetE());
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
 * Fills the Minuit parameter array from initial values in the CNuc object.
 */

void CNuc::FillMnParams(AZUREParams &params, const Config *config) {
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  std::string varname;
  int energyIndex = 1;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      ALevel *level = this->Group(j)->GetLevel(la);
      varname = "energy_" + std::to_string(energyIndex);
      params.Add(varname, level->GetE(), 0.1 * level->GetE(), ParameterRecord::LevelEnergy(j, la, energyIndex));
      bool isUnbound = false;
      for (int ir = 1; ir <= this->NumPairs(); ir++) {
        PPair *pair = this->Pair(ir);
        if (pair->GetNucleus() != this->Group(j)->GetNucleus()) continue;   // thresholds of this nucleus only
        if (pair->GetPType() == 0 &&
            level->GetE() > (pair->GetSepE() + pair->GetExE())) isUnbound = true;
      }
      // Jakub's Fix: it can create problems with the fit with pyazr and brick
      if (!isUnbound) p.Fix(varname);
      if (level->EnergyFixed() && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);

      // Parameter settings will be applied by ParameterLimitsManager during fit
      int widthIndex = 1;
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        varname = "width_" + std::to_string(energyIndex) + "_" + std::to_string(widthIndex);
        widthIndex++;
        if (gammaDecays_->AddParameter(params, gammaDecayInputs_, j, la, ch, energyIndex)) continue;  // a gamma_decay record's number
        params.Add(varname, level->GetGamma(ch), 0.1 * level->GetGamma(ch), ParameterRecord::Width(j, la, ch, energyIndex));
        // A width given as 0 is not fitted, in every basis: the input decides, not the transformed amplitude, which
        // in the formal basis is 0 or rounding by the order of the sums (external_capture.md section 9).  Without a
        // transformation the input is the parameter.
        if (this->IsTransformedIn() ? level->ZeroInputWidth(ch) : level->GetGamma(ch) == 0.0) p.Fix(varname);
        if (level->ChannelFixed(ch) && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);
        // Wigner-limit bounds: GetWignerLimit is gamma^2_W = 3 hbar^2/(2 mu a^2) (MeV), the limit on the reduced width
        // squared; the fit parameter is the amplitude (MeV^1/2), so the bound is its square root (gamma^2_W itself
        // would let theta^2 = gamma^2/gamma^2_W run past 1 by a factor gamma^2_W).
        if (config && (config->paramMask & Config::USE_WIGNER_LIMITS)) {
          AChannel *channel = this->Group(j)->GetChannel(ch);
          double wignerLimit = channel->GetWignerLimit();
          if (channel->GetRadType() == 'P' && wignerLimit > 0.0 && std::isfinite(wignerLimit)) {
            double gammaLimit = sqrt(wignerLimit);
            p.SetLimits(varname, -gammaLimit, gammaLimit);
          }
        }
      }
      energyIndex++;
    }
  }
  // The isotopic abundances follow, written by WorkingCopy::FillMnParams (AbundanceSamples::AddParameters).
}

/*!
 * Fills the CNuc object from the Minuit parameter array.
 */

void CNuc::FillCompoundFromParams(const vector_r &params) {
  const vector_r &p = gammaDecays_->Empty() ? params : gammaDecays_->ToEngine(*this, params);  // GammaDecayMap.h
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      ALevel *level = this->Group(j)->GetLevel(la);
      level->SetFitE(p[i]);
      i++;
      std::vector<double> g;
      for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++) {
        level->SetFitGamma(ch, p[i]);
        level->SetECGamma(ch, p[i]);
        g.push_back(p[i]);
        i++;
      }
      level->SetSqrtNFFactor(SqrtNFFactor(*this, this->Group(j), level, 1.0, g));
    }
  }
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal();
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams();
}

/*!
 * Performs the final parameter transformations from formal to physical parameters.
 */

void CNuc::TransformOut(const Config &configure) {
  this->TransformToBrune(configure, true);
  this->ObservedWidths(configure);
}

/*!
 * The Brune (Park) values of the fit parameters, E and gamma, in the level's transform slots: in the formal basis by
 * physics/ brune_from_formal, otherwise a copy.  report: print the warnings for levels that do not convert.
 */

void CNuc::TransformToBrune(const Config &configure, bool report) {
  // Park (2021): the fit parameters already are the observed reduced width
  // amplitudes, so Gamma = 2 P gamma^2 with no 1 + sum gamma^2 dS/dE factor.
  // Refresh J and the bound-state normalization for the current parameters.
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  if (park) {
    this->CalcShiftFunctions(configure);
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
        ALevel *theLevel = this->Group(j)->GetLevel(la);
        if (report && theLevel->IsInRMatrix() && theLevel->GetParkNorm() <= 0.0) WarnParkNormNotPositive(*this, configure, j, la);
      }
  }
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    // Formal -> Brune by the fixed-point iteration of physics/ brune_from_formal (level_parameters.md section 4), on
    // the fit values.  A level not in the R-matrix keeps its values.
    const phys::FormalToBruneOptions options;  // the engine's 1e-6 MeV and 1000 iterations
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
        ALevel *theLevel = this->Group(j)->GetLevel(la);
        if (theLevel->IsInRMatrix() && this->Group(j)->IsInRMatrix()) continue;
        theLevel->SetTransformE(theLevel->GetFitE());
        theLevel->SetTransformIterations(0);
        for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++)
          theLevel->SetTransformGamma(ch, theLevel->GetFitGamma(ch));
      }
    const CNucSystem &structure = this->SystemStructure();
    const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      const int j = structure.map.groupOfBlock[b];
      JGroup *theJGroup = this->Group(j);
      const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
      const std::size_t nc = theJGroup->NumChannels();
      std::vector<double> E, gamma, B(nc);
      for (int la : levels) {
        E.push_back(theJGroup->GetLevel(la)->GetFitE());
        for (std::size_t k = 0; k < nc; k++) gamma.push_back(theJGroup->GetLevel(la)->GetFitGamma((int)k + 1));
      }
      for (std::size_t k = 0; k < nc; k++) B[k] = theJGroup->GetChannel((int)k + 1)->GetBoundaryCondition();
      const phys::LevelSet out = phys::brune_from_formal(structure.system, b, E, gamma, B, cf, options);
      for (std::size_t lb = 0; lb < levels.size(); lb++) {
        const int la = levels[lb];
        ALevel *theLevel = theJGroup->GetLevel(la);
        if (report) WarnFormalToBrune(*this, configure, j, la, out, lb, nc, options);
        theLevel->SetTransformE(out.E[lb]);
        theLevel->SetTransformIterations(out.iterations[lb]);
        const std::vector<double> g(out.gamma.begin() + lb * nc, out.gamma.begin() + (lb + 1) * nc);
        for (std::size_t k = 0; k < nc; k++) theLevel->SetTransformGamma((int)k + 1, g[k]);
        theLevel->SetSqrtNFFactor(SqrtNFFactor(*this, theJGroup, theLevel, 1.0, g));
      }
    }
  } else {
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
        this->Group(j)->GetLevel(la)->SetTransformIterations(0);
        this->Group(j)->GetLevel(la)->SetTransformE(this->Group(j)->GetLevel(la)->GetFitE());
        for (int ch = 1; ch <= this->Group(j)->NumChannels(); ch++)
          this->Group(j)->GetLevel(la)->SetTransformGamma(ch, this->Group(j)->GetLevel(la)->GetFitGamma(ch));
      }
  }

}

/*!
 * The observed quantities of every level (cnuc_split.md section 8) from its Brune (Park) values in the transform
 * slots (TransformToBrune) and, with external capture, the external part of its photon amplitudes: physics/
 * observed_levels.  Sets the external amplitude, W (SetBigGamma) and the reported value with its kind.
 */

void CNuc::ObservedWidths(const Config &configure) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool external = (configure.paramMask & Config::USE_EXTERNAL_CAPTURE);
  const bool ignoreZero = (configure.paramMask & Config::IGNORE_ZERO_WIDTHS);
  const CNucSystem &all = this->SystemAllGroups();
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, all.map, configure);
  for (std::size_t b = 0; b < all.system.num_blocks(); b++) {
    JGroup *theJGroup = this->Group(all.map.groupOfBlock[b]);
    const std::vector<int> &levels = all.map.levelOfBlockLevel[b];
    const std::size_t nc = theJGroup->NumChannels();
    std::vector<double> E, gamma;
    std::vector<unsigned char> rwa;
    std::vector<phys::ExternalWidth> ext(levels.size() * nc);
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      ALevel *theLevel = theJGroup->GetLevel(levels[lb]);
      E.push_back(theLevel->GetTransformE());
      for (std::size_t k = 0; k < nc; k++) {
        const int ch = (int)k + 1;
        AChannel *theChannel = theJGroup->GetChannel(ch);
        gamma.push_back(theLevel->GetTransformGamma(ch));
        rwa.push_back(theLevel->GammaIsRWA(ch) ? 1 : 0);
        complex externalWidth(0.0, 0.0);
        if ((theChannel->GetRadType() == 'M' || theChannel->GetRadType() == 'E') && theLevel->IsInRMatrix() &&
            external && !(ignoreZero && theLevel->ZeroInputWidth(ch))) {
          externalWidth = ExternalReducedWidth(*this, theJGroup, theLevel, theChannel, false, configure);
          ext[lb * nc + k].present = true;
          ext[lb * nc + k].value = externalWidth;
        }
        theLevel->SetExternalGamma(ch, externalWidth);
      }
    }
    const phys::ObservedLevels obs =
        phys::observed_levels(all.system, b, phys::internal_channels(all.system, b), E, gamma, ext, rwa,
                              park ? phys::LevelBasis::Park : phys::LevelBasis::Brune, cf);
    for (std::size_t lb = 0; lb < levels.size(); lb++)
      for (std::size_t k = 0; k < nc; k++) {
        ALevel *theLevel = theJGroup->GetLevel(levels[lb]);
        theLevel->SetBigGamma((int)k + 1, obs.W[lb * nc + k]);
        theLevel->SetObserved((int)k + 1, obs.value[lb * nc + k], (unsigned char)obs.kind[lb * nc + k]);
      }
  }
}

/*!
 * The observed widths of a parameter vector, as parameters.out reports them (TransformOut without its messages):
 * particle and photon widths with the level's normalization and, with external capture, the external part of the
 * photon amplitude (physics/ observed_levels), in the levels' observed slots.  The fit state is set back to the vector
 * afterwards, since TransformToBrune also sets the bound-state normalizations from the Brune values.
 */

bool CNuc::ObservedWidthsAt(const Config &configure, const vector_r &params) {
  this->FillCompoundFromParams(params);
  try {
    this->TransformToBrune(configure, false);
    this->ObservedWidths(configure);
  } catch (const std::exception &) {
    this->FillCompoundFromParams(params);
    return false;
  }
  this->FillCompoundFromParams(params);
  return true;
}

/* Gets the transformet parameters*/

vector_r CNuc::GetTransformParams(const Config &) {
  // The observed values TransformOut computed (physics/ observed_levels, cnuc_split.md section 8): Gamma (eV),
  // C (fm^-1/2), mu, Q, B and the Brune amplitude of a gammaIsRWA column, signed, as a <levels> line holds them.
  vector_r params;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    JGroup *theJGroup = this->Group(j);
    for (int la = 1; la <= theJGroup->NumLevels(); la++) {
      ALevel *theLevel = theJGroup->GetLevel(la);
      params.push_back(theLevel->GetTransformE());
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) params.push_back(gammaDecays_->ValueOr(*this, j, la, ch, theLevel->GetObservedValue(ch)));
    }
  }
  // The level block only: the caller appends the abundance block of the fit's vector (AbundanceSamples::AppendTo).
  return params;
}

/*!
 * This function is called for each iteration to calculate the shift
 * functions at new level energies when the Brune parametrization is used.
 *
 * In the Park formalism the first and second energy derivatives at the level
 * energies are numerical (ChannelFunc::ShiftDerivative, CoulFunc::PEShift_d2E,
 * ShftFunc), four
 * and five more shift functions per channel and level, most of a serial
 * evaluation at new parameters.  They are taken after the shift functions, for
 * all levels at once and in parallel; the Coulomb functions are computed at
 * every call (physics/coulomb, no memo), so the values do not depend on the
 * number of threads.  Channels of a level that differ only in
 * the channel spin share one derivative, and the working copies of a fit share
 * them through the run's exact-argument memo (RunScratch::shiftDerivatives).
 */

void CNuc::CalcShiftFunctions(const Config &configure, bool secondDerivatives) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  struct Derivative {
    ALevel *level;
    PPair *pair;
    int ch, l, two_j;  ///< two_j: the channel's 2J key (PPair::SpinOrbitTwoJ)
    double energy;
    bool first, second;
    double d1, d2;
    int source;  ///< an earlier entry of the same level, pair and l (another channel spin), or -1
  };
  std::vector<Derivative> derivatives;
  struct Moved {
    JGroup *group;
    ALevel *level;
    double energy;
    bool recompute, second;
  };
  std::vector<Moved> levels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    if (this->Group(j)->IsInRMatrix()) {
      JGroup *theJGroup = this->Group(j);
      for (int la = 1; la <= theJGroup->NumLevels(); la++) {
        ALevel *theLevel = theJGroup->GetLevel(la);
        if (theLevel->IsInRMatrix()) {
          const double levelEnergy = theLevel->GetFitE();
          // The shift functions depend on the level energy alone, and a
          // minimizer moves one parameter at a time, so most levels have not
          // moved since the last call.  Park's J below depends on the widths
          // too and is always refreshed.
          const bool recompute = !theLevel->ShiftFunctionsValidAt(levelEnergy);
          // d^2S/dE^2 (Park) enters only the gradient: the numerical second
          // derivative costs five shift functions per channel, so a chi-squared
          // alone leaves it out, and a later call that asks for it fills it in.
          const bool second = park && secondDerivatives && !theLevel->ShiftSecondDerivativesValidAt(levelEnergy);
          const int firstOfLevel = (int)derivatives.size();
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
            AChannel *theChannel = theJGroup->GetChannel(ch);
            PPair *thePair = this->Pair(theChannel->GetPairNum());
            if (thePair->GetPType() == 0) {
              if (!recompute && !second) continue;
              int lValue = theChannel->GetL();
              const bool parkP = park && theChannel->GetRadType() == 'P';
              double resonanceEnergy = levelEnergy - (thePair->GetSepE() + thePair->GetExE());
              if (recompute) {
                // S as everywhere else (ChannelFunc: continuous through threshold).
                ChannelFunc channelFunc(thePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
                theLevel->SetShiftFunction(ch, channelFunc.Shift(lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), resonanceEnergy));
              }
              if (parkP && (recompute || second)) {
                // Channels that differ only in the channel spin share the derivative.
                int source = -1;
                for (int k = firstOfLevel; k < (int)derivatives.size() && source < 0; k++)
                  if (derivatives[k].source < 0 && derivatives[k].pair == thePair && derivatives[k].l == lValue && derivatives[k].two_j == thePair->SpinOrbitTwoJ(theJGroup->GetJ())) source = k;
                derivatives.push_back({theLevel, thePair, ch, lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), levelEnergy, recompute, second, 0.0, 0.0, source});
              }
            } else {
              theLevel->SetShiftFunction(ch, theJGroup->GetLevel(1)->GetShiftFunction(1));
            }
          }
          levels.push_back({theJGroup, theLevel, levelEnergy, recompute, second});
        }
      }
    }
  }

  // Park: dS/dE and d^2S/dE^2, in parallel (see above); a GSL error is
  // carried out of the parallel region and rethrown.
  const bool gslCoulomb = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  const int numDerivatives = (int)derivatives.size();
  azure2::ShiftDerivativeMemo &memo = configure.scratch->shiftDerivatives;
  std::exception_ptr error;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) if (numDerivatives > 1 && !omp_in_parallel())
#endif
  for (int k = 0; k < numDerivatives; k++) {
    try {
      Derivative &d = derivatives[k];
      if (d.source >= 0) continue;
      const double threshold = d.pair->GetSepE() + d.pair->GetExE();
      const double resonanceEnergy = d.energy - threshold;
      CoulFunc theCoulombFunction(d.pair, gslCoulomb);
      const azure2::ShiftDerivativeMemo::Key key{d.pair->GetZ(1), d.pair->GetZ(2), d.l, d.two_j, gslCoulomb, azure2::physics::PotentialKey(theCoulombFunction.HybridPotential()),
                              theCoulombFunction.UsesTable(), d.pair->FunctionSettings().relativistic, d.pair->GetRedMass(),
                              d.pair->GetM(1) + d.pair->GetM(2), d.pair->GetChRad(), threshold, d.energy};
      bool need1 = d.first, need2 = d.second;
      {
        std::lock_guard<std::mutex> lock(memo.mutex);
        auto it = memo.values.find(key);
        if (it != memo.values.end()) {
          if (need1 && it->second.has1) {
            d.d1 = it->second.d1;
            need1 = false;
          }
          if (need2 && it->second.has2) {
            d.d2 = it->second.d2;
            need2 = false;
          }
        }
      }
      if (!need1 && !need2) continue;
      // dS/dE as everywhere else (ChannelFunc, also BuildShiftDerivTable's
      // under Brune); d^2S/dE^2 from the branch the level is on.
      if (need1) {
        ChannelFunc channelFunc(d.pair, gslCoulomb);
        d.d1 = channelFunc.ShiftDerivative(d.l, d.two_j, resonanceEnergy);
      }
      if (need2) {
        if (resonanceEnergy < 0.0) {
          ShftFunc theShiftFunction(d.pair);
          d.d2 = theShiftFunction.EnergySecondDerivative(d.l, d.two_j, d.energy);
        } else {
          d.d2 = theCoulombFunction.PEShift_d2E(d.l, d.two_j, d.pair->GetChRad(), resonanceEnergy);
        }
      }
      std::lock_guard<std::mutex> lock(memo.mutex);
      if (memo.values.size() >= azure2::ShiftDerivativeMemo::kMax) memo.values.clear();
      azure2::ShiftDerivativeMemo::Value &v = memo.values[key];
      if (need1) {
        v.has1 = true;
        v.d1 = d.d1;
      }
      if (need2) {
        v.has2 = true;
        v.d2 = d.d2;
      }
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(calc_shift_functions_error)
#endif
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
  for (Derivative &d : derivatives) {
    if (d.source >= 0) {
      d.d1 = derivatives[d.source].d1;
      d.d2 = derivatives[d.source].d2;
    }
    if (d.first) d.level->SetShiftDerivative(d.ch, d.d1);
    if (d.second) d.level->SetShiftSecondDerivative(d.ch, d.d2);
  }

  for (const Moved &m : levels) {
    ALevel *theLevel = m.level;
    JGroup *theJGroup = m.group;
    if (m.recompute) theLevel->SetShiftFunctionsEnergy(m.energy);
    if (m.second) theLevel->SetShiftSecondDerivativesEnergy(m.energy);
    if (park) {
      // Park (2021) Eq. (28): J = 1 - sum_c gamma_c^2 (dS_c/dE) at E_lambda.
      double parkNorm = 1.0;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        AChannel *theChannel = theJGroup->GetChannel(ch);
        if (theChannel->GetRadType() == 'P' && this->Pair(theChannel->GetPairNum())->GetPType() == 0)
          parkNorm -= pow(theLevel->GetFitGamma(ch), 2.0) * theLevel->GetShiftDerivative(ch);
      }
      theLevel->SetParkNorm(parkNorm);
      // N_f^{1/2} multiplies the level's reduced width amplitude wherever
      // it serves as an external-capture final state.  It is defined with
      // Brune's amplitudes, gamma_Park / sqrt(J); folding the 1/sqrt(J)
      // in here keeps gamma * N_f^{1/2} the same number in both
      // parametrizations.
      std::vector<double> g;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) g.push_back(theLevel->GetFitGamma(ch));
      theLevel->SetSqrtNFFactor(SqrtNFFactor(*this, theJGroup, theLevel, parkNorm, g));
    }
  }
}

/*!
 * Converts a parameter set read from a file written in the other alternative
 * parametrization.  Brune -> Park multiplies every amplitude of a level by
 * sqrt(J), J = 1/(1 + sum_c gamma~^2 dS_c/dE); Park -> Brune divides by
 * sqrt(J), J = 1 - sum_c gamma^2 dS_c/dE.  Photon and beta amplitudes rescale
 * with the level like the particle ones.  Levels outside the R matrix
 * (external-capture final states) are the same in both and are left alone.
 * The parameter order is that of FillMnParams.
 */

bool CNuc::ConvertAmplitudeBasis(ROOT::Minuit2::MnUserParameters &p, int from, int to, const Config &configure) {
  if (from == to) return true;
  const bool toPark = (to == 2);
  if (!((from == 1 && to == 2) || (from == 2 && to == 1))) return false;
  // Index of each level's energy in the parameter vector (E, then one amplitude per channel, group by group).
  std::map<std::pair<int, int>, int> energyIndex;
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++) {
      energyIndex[{j, la}] = i;
      i += 1 + this->Group(j)->NumChannels();
    }
  // physics/ park_from_brune, brune_from_park on each R-matrix block (level_parameters.md Eq. 2.3); the errors
  // scale with the amplitudes, by J^1/2 or J^-1/2.
  const CNucSystem &structure = this->SystemStructure();
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
  for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
    const int j = structure.map.groupOfBlock[b];
    const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
    const int nc = this->Group(j)->NumChannels();
    std::vector<double> E, gamma, J;
    for (int la : levels) {
      const int e = energyIndex[{j, la}];
      E.push_back(p.Value(e));
      for (int ch = 1; ch <= nc; ch++) gamma.push_back(p.Value(e + ch));
    }
    std::vector<double> out;
    try {
      out = toPark ? phys::park_from_brune(structure.system, b, E, gamma, cf, &J)
                   : phys::brune_from_park(structure.system, b, E, gamma, cf, &J);
    } catch (const std::domain_error &error) {
      configure.outStream << "**ERROR: " << error.what() << std::endl
                          << "  No Brune (or standard) parameters correspond to these widths." << std::endl;
      return false;
    }
    for (std::size_t lb = 0; lb < levels.size(); lb++) {
      const int e = energyIndex[{j, levels[lb]}];
      const double factor = toPark ? std::sqrt(J[lb]) : 1.0 / std::sqrt(J[lb]);
      for (int ch = 1; ch <= nc; ch++) {
        p.SetValue(e + ch, out[lb * nc + ch - 1]);
        p.SetError(e + ch, p.Error(e + ch) * factor);
      }
    }
  }
  return true;
}

/*!
 * Forgets the cached shift functions of every level and the run's shift-derivative memo.
 */

void CNuc::InvalidateShiftFunctions(const Config &configure) {
  {
    azure2::ShiftDerivativeMemo &memo = configure.scratch->shiftDerivatives;
    std::lock_guard<std::mutex> lock(memo.mutex);
    memo.values.clear();
  }
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      this->Group(j)->GetLevel(la)->InvalidateShiftFunctions();
}

/*!
 * Park formalism: sum over the R-matrix levels of (J/kParkNormScale)^2 for
 * J < 0.  CalcShiftFunctions must have run for the current parameters.  A
 * level with J <= 0 has no standard R-matrix counterpart (its Brune reduced
 * widths would be imaginary), and nothing in the level matrix itself keeps a
 * minimizer out of that region, so this is the wall.
 */

double CNuc::ParkNormPenalty() {
  double penalty = 0.0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    JGroup *theJGroup = this->Group(j);
    if (!theJGroup->IsInRMatrix()) continue;
    for (int la = 1; la <= theJGroup->NumLevels(); la++) {
      ALevel *theLevel = theJGroup->GetLevel(la);
      if (!theLevel->IsInRMatrix()) continue;
      double parkNorm = theLevel->GetParkNorm();
      if (parkNorm < 0.0) penalty += pow(parkNorm / kParkNormScale, 2.0);
    }
  }
  return penalty;
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

ALevel *CNuc::MutableLevel(int j, int la) {
  return jgroups_[j - 1].GetLevel(la);
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
