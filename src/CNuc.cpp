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

bool CNuc::IsPairKey(int key) {
  bool b = false;
  int c = 0;
  while (!b && c < this->NumPairs()) {
    if (key == this->GetPair(c + 1)->GetPairKey()) b = true;
    c++;
  }
  return b;
}

std::vector<int> CNuc::CapturePairKeys(int entranceKey) {
  std::vector<int> keys;
  PPair *entrance = GetPairFromKey(entranceKey);
  if (!entrance) return keys;
  for (int i = 1; i <= NumPairs(); i++)
    if (GetPair(i)->GetPType() == 10 && GetPair(i)->GetNucleus() == entrance->GetNucleus())
      keys.push_back(GetPair(i)->GetPairKey());
  return keys;
}

std::string CNuc::RatePairsProblem(int entranceKey, int exitKey) {
  PPair *entrance = GetPairFromKey(entranceKey);
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
    if (pair.GetPairKey() == this->GetPair(c + 1)->GetPairKey()) b = true;
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

int CNuc::NumNuclei() const {
  return nuclei_.size();
}

CompoundNucleus *CNuc::GetNucleus(int nucleusNum) {
  return &nuclei_[nucleusNum - 1];
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

int CNuc::GetPairNumFromKey(int key) {
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
      PPair NewPair(Line);
      if (!this->IsPair(NewPair)) this->AddPair(NewPair);
    }
    if (Line.isActive() == 1) {
      PPair NewPair(Line);
      PairNum = this->IsPair(NewPair);
      if (!PairNum) {
        this->AddPair(NewPair);
        PairNum = this->IsPair(NewPair);
      }
      if (Line.ecMultMask() != 0) {
        std::map<int, int>::iterator it = ecPairs.find(PairNum);
        if (it == ecPairs.end()) ecPairs[PairNum] = Line.ecMultMask();
      }
      const int nucleusNum = this->GetPair(PairNum)->GetNucleus();
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
      NewJGroup.SetNucleus(nucleusNum, this->GetPair(PairNum)->GetNucleusZ(), this->GetPair(PairNum)->GetNucleusA());
      JGroupNum = this->IsJGroup(NewJGroup);
      if (!JGroupNum) {
        this->AddJGroup(NewJGroup);
        JGroupNum = this->IsJGroup(NewJGroup);
      }
      AChannel NewChannel(Line, PairNum);
      ChannelNum = this->GetJGroup(JGroupNum)->IsChannel(NewChannel);
      if (!ChannelNum) {
        this->GetJGroup(JGroupNum)->AddChannel(NewChannel);
        ChannelNum = this->GetJGroup(JGroupNum)->IsChannel(NewChannel);
        if (this->GetJGroup(JGroupNum)->GetChannel(ChannelNum)->GetL() > maxLValue &&
            this->GetJGroup(JGroupNum)->GetChannel(ChannelNum)->GetRadType() == 'P')
          maxLValue = this->GetJGroup(JGroupNum)->GetChannel(ChannelNum)->GetL();

        // Calculate and set Wigner Limit for the newly created channel
        PPair *channelPair = this->GetPair(PairNum);
        this->GetJGroup(JGroupNum)->GetChannel(ChannelNum)->SetWignerLimit(channelPair->GetRedMass(), channelPair->GetChRad());
      }
      // A line belongs to a level of its J-group that carries the same level ID within 1 keV, or failing that
      // to one at the identical energy.  Until 2026-10-24 any level within 1 keV was taken (JGroup::IsLevel),
      // whatever its ID, so a narrow level put next to a broad one of the same J^pi (17O+a: 3/2- at 8.2900 and
      // 8.2893 MeV, IDs 91 and 22) was silently merged into it: its widths overwrote the broad level's and the
      // fit lost three parameters.  Lines with other IDs within 1 keV now make their own level, with a warning.
      // The level of the line: io::LevelGrouper, the rule pyazr uses too (docs/upgrade/project.md, P2).
      ALevel NewLevel(Line);
      JGroup *group = this->GetJGroup(JGroupNum);
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
      this->GetJGroup(JGroupNum)->GetLevel(LevelNum)->SetGammaFromLine(ChannelNum, Line);
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
    JGroup *group = this->GetJGroup(entry.first.first);
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
                          << this->GetPair(c->GetPairNum())->GetPairKey() << ", l = " << c->GetL() << ", s = " << c->GetS()
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
    PPair *changedPair = this->GetPair(it->first);
    changedPair->SetChRad(it->second);
    // The Wigner limits above were computed inside the parse loop, i.e. with
    // the radius as written in the .azr.  Redo them for the overridden pair,
    // or every theta^2 reported after a radius change is against the old limit.
    for (int j = 1; j <= this->NumJGroups(); j++) {
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) {
        AChannel *theChannel = this->GetJGroup(j)->GetChannel(ch);
        if (theChannel->GetPairNum() == it->first)
          theChannel->SetWignerLimit(changedPair->GetRedMass(),
                                     changedPair->GetChRad());
      }
    }
  }

  for (int p = 1; p <= this->NumPairs(); p++) this->GetPair(p)->SetHybridPotential(configure.potentials.PotentialFor(this->GetPair(p)->GetPairKey()));  // a value, none: plain Coulomb
  if (configure.includeAllJPi)
    this->AddHardSphereGroups(configure, configure.hardSphereMaxL >= 0 ? configure.hardSphereMaxL : maxLValue);
  if (RefuseForbiddenIdenticalChannels(*this, configure) == -1) return -1;
  if ((configure.paramMask & Config::USE_EXTERNAL_CAPTURE) && this->NumJGroups() > 0 && this->NumPairs() > 0)
    if (this->ParseExternalCapture(configure, ecPairs) == -1) return -1;

  if (this->NumNuclei() > 1)
    for (int j = 1; j <= this->NumJGroups(); j++) this->GetJGroup(j)->SetNucleusLabelled(true);

  if (this->ReadAbundances(configure) == -1) return -1;
  if (std::shared_ptr<const ModelDataSet> read = ModelData::Read(*this, configure, section.rows)) modelData_ = read; else return -1;
  // right after the pairs are known, before any data file is read (a project whose data fail still gets it)
  WarnMassSeparationMismatch(*this, configure);

  return 0;
}

/*!
 * The abundances of the samples (AbundancesSpec): one row per compound nucleus and sample, or a sample's sum-prior
 * width.  A nucleus without a row keeps 1.0, fixed.  A row for a (Z, A) the pairs do not form is reported and
 * ignored.  Without the block no abundance becomes a
 * Minuit parameter, so every project written before the block existed keeps
 * its parameter layout, param.par and <parameterSettings> indices.
 */

int CNuc::SetAbundances(const Config &configure, const azure2::experiment::AbundanceSpecs &specs) {
  hasAbundanceBlock_ = specs.abundances.present;
  samples_.clear();
  if (!hasAbundanceBlock_) return 0;
  // A sample is created the first time a row sets one of its values; the unnamed sample is "".
  auto sampleOf = [&](const std::string &name) -> AbundanceSample & {
    for (unsigned int s = 0; s < samples_.size(); s++)
      if (samples_[s].name == name) return samples_[s];
    AbundanceSample sample;
    sample.name = name;
    sample.entries.assign(nuclei_.size(), SampleAbundance());
    samples_.push_back(sample);
    return samples_.back();
  };
  for (const azure2::experiment::AbundanceRowSpec &row : specs.abundances.rows) {
    if (row.isSum) {
      sampleOf(row.sample).sumSigma = row.sumSigma;
      continue;
    }
    bool found = false;
    for (unsigned int i = 0; i < nuclei_.size(); i++) {
      if (nuclei_[i].Z != row.z || nuclei_[i].A != row.a) continue;
      SampleAbundance &e = sampleOf(row.sample).entries[i];
      e.abundance = e.nominal = row.value;
      e.vary = row.vary;
      e.sigma = row.sigma;
      e.isotope = row.isotope;
      found = true;
    }
    if (!found)
      configure.outStream << "WARNING: <abundances> names Z=" << row.z << " A=" << row.a
                          << ", which no pair of the project forms; ignored." << std::endl;
  }
  if (specs.failed) {
    if (!specs.error.empty()) configure.outStream << specs.error << std::endl;
    return -1;
  }
  // An empty block still means "abundances are parameters": one unnamed sample at the defaults.
  if (samples_.empty()) sampleOf("");
  return 0;
}

int CNuc::FindSample(const std::string &name) const {
  for (unsigned int s = 0; s < samples_.size(); s++)
    if (samples_[s].name == name) return (int)s + 1;
  return 0;
}

int CNuc::DefaultSample() const {
  int unnamed = this->FindSample("");
  if (unnamed > 0) return unnamed;
  return samples_.size() == 1 ? 1 : 0;
}

void CNuc::AbundanceOf(int ordinal, int &sample, int &nucleus) const {
  const int n = (int)nuclei_.size();
  sample = (ordinal - 1) / n + 1;
  nucleus = (ordinal - 1) % n + 1;
}

std::string CNuc::AbundanceParamName(int sample, int nucleus) const {
  const std::string &name = samples_[sample - 1].name;
  const std::string &isotope = samples_[sample - 1].entries[nucleus - 1].isotope;
  if (!isotope.empty()) return name + ".isotope." + isotope;  // format 3: <target>.isotope.<isotope>
  const std::string zA = std::to_string(nuclei_[nucleus - 1].Z) + "_" + std::to_string(nuclei_[nucleus - 1].A);
  return "abundance_" + (name.empty() ? std::string() : name + "_") + zA;
}

double CNuc::GetAbundance(int sample, int nucleus) const {
  if (sample < 1 || sample > (int)samples_.size()) return 1.0;
  return samples_[sample - 1].entries[nucleus - 1].abundance;
}

double CNuc::AbundanceSum(int sample) const {
  double sum = 0.0;
  for (unsigned int i = 0; i < nuclei_.size(); i++) sum += samples_[sample - 1].entries[i].abundance;
  return sum;
}

double CNuc::SampleSumSigma(int sample, double globalSumSigma) const {
  const double own = samples_[sample - 1].sumSigma;
  return own >= 0.0 ? own : globalSumSigma;
}

double CNuc::AbundancePenalty(double globalSumSigma) const {
  if (!hasAbundanceBlock_) return 0.0;
  double penalty = 0.0;
  for (unsigned int s = 0; s < samples_.size(); s++) {
    for (unsigned int i = 0; i < nuclei_.size(); i++) {
      const SampleAbundance &e = samples_[s].entries[i];
      if (e.sigma > 0.0 && !e.priorReplaced) penalty += pow((e.abundance - e.nominal) / e.sigma, 2.0);
    }
    const double sig = this->SampleSumSigma((int)s + 1, globalSumSigma);
    if (sig > 0.0) penalty += pow((this->AbundanceSum((int)s + 1) - 1.0) / sig, 2.0);
  }
  return penalty;
}

/*!
 * The sum-to-one part of AbundancePenalty alone: ((sum_i w_i) - 1)^2 / s^2 per sample with s > 0.  The MCMC
 * likelihood uses this, because AZURECalcMCMC::BuildAutoPriors already turns each abundance's (w0, sigma) into a
 * Gaussian prior (PARAM_ABUNDANCE); adding the full AbundancePenalty there counted that prior twice and made the
 * posterior width of an abundance 1/sqrt(2) too narrow (found 2026-10-17 by the isospin-levels session).
 */
double CNuc::AbundanceSumPenalty(double globalSumSigma) const {
  if (!hasAbundanceBlock_) return 0.0;
  double penalty = 0.0;
  for (unsigned int s = 0; s < samples_.size(); s++) {
    const double sig = this->SampleSumSigma((int)s + 1, globalSumSigma);
    if (sig > 0.0) penalty += pow((this->AbundanceSum((int)s + 1) - 1.0) / sig, 2.0);
  }
  return penalty;
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
    PPair *pair = this->GetPair(aa);
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
          JGroup *jg = this->GetJGroup(jGroupNum);
          for (int ir = 1; ir <= this->NumPairs(); ir++) {
            PPair *chPair = this->GetPair(ir);
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
 * Fills the ECLevel vector from information in the external capture file.  Also tests if the final state for external capture
 * exists from the nuclear file.  If not, the state is created.
 */

int CNuc::ParseExternalCapture(const Config &configure, std::map<int, int> &ecPairs) {
  std::set<int> usedExitKeys = ProjectExitPairKeys(configure);  // the exit pairs of the active data sets (their specs)
  // A total-capture rate (exit key -1) exits into every capture pair of its entrance nucleus.
  if ((configure.paramMask & Config::CALCULATE_REACTION_RATE) && configure.rateParams.exitPair == -1)
    for (int key : CapturePairKeys(configure.rateParams.entrancePair)) usedExitKeys.insert(key);
  for (std::map<int, int>::iterator ec = ecPairs.begin(); ec != ecPairs.end(); ec++) {
    PPair *exitPair = this->GetPair(ec->first);
    if (exitPair->GetPType() != 10) {
      configure.outStream << "Final state is not a capture pair." << std::endl;
      continue;
    }
    // The final state of an external-capture pair must be a level of the project (user guide, "Gamma-Ray
    // Calculations"): its ANC and gamma widths live on that level's lines.  Before 2026-10-26 a missing level
    // was silently created here with a fixed energy and FREE particle widths of 0.1 MeV^1/2 that were in no
    // .azr line and not in the GUI, and -- because it was appended to its J-group only when external capture
    // was on, i.e. only in modes with a capture segment on the CLI, always in the API -- every later level's
    // energy_N / width_N_c name moved by one between runs of the same project, so a param.sav landed on the
    // wrong levels without a message.  Now: a pair that an active segment (or the rate pair) exits into is
    // refused; one that nothing uses (a level-isolation or scratch model whose final-state lines are switched
    // off) is warned about and gets no external capture and no hidden level.  tests/ec_final_state_listed.
    double jValue = exitPair->GetJ(2);
    int parity = exitPair->GetPi(2);
    JGroup newJGroup(jValue, parity);
    newJGroup.SetNucleus(exitPair->GetNucleus(), exitPair->GetNucleusZ(), exitPair->GetNucleusA());
    int jGroupNum = this->IsJGroup(newJGroup);
    ALevel probeLevel(exitPair->GetHeavyExE());
    if (!jGroupNum || !this->GetJGroup(jGroupNum)->IsLevel(probeLevel)) {
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
    const int levelNum = this->GetJGroup(jGroupNum)->IsLevel(probeLevel);
    // IsLevel is an energy-only match within 1 keV and returns the first level that fits.  Since e239a12 two
    // levels of one J^pi may sit within 1 keV of each other (different level IDs), so the pair's E_x can match
    // more than one: refuse that instead of taking the first silently (audit 2026-10-25, checklist item 28).
    {
      int matches = 0;
      for (int lv = 1; lv <= this->GetJGroup(jGroupNum)->NumLevels(); lv++)
        if (fabs(this->GetJGroup(jGroupNum)->GetLevel(lv)->GetE() - exitPair->GetHeavyExE()) <= 1.0e-3) matches++;
      if (matches > 1) {
        configure.outStream << "ERROR: external capture into pair " << exitPair->GetPairKey() << ": " << matches
                            << " levels of J^pi = " << jValue << (parity > 0 ? "+" : "-") << " lie within 1 keV of its final-state energy E_x = "
                            << exitPair->GetHeavyExE() << " MeV, so the final state is ambiguous.  Move the extra level(s) by more than 1 keV"
                            << " or give the pair the energy of the one it decays to." << std::endl;
        return -1;
      }
    }
    if (FinalStateOfAnotherPair(*this, configure, jGroupNum, levelNum, ec->first)) return -1;
    this->GetJGroup(jGroupNum)->GetLevel(levelNum)->SetECParams(ec->first, ec->second);
    for (int ch = 1; ch <= this->GetJGroup(jGroupNum)->NumChannels(); ch++) {
      PPair *theFinalPair = this->GetPair(this->GetJGroup(jGroupNum)->GetChannel(ch)->GetPairNum());
      double nfIntegralValue = 0.;
      double ecConvert = 0.;
      if (theFinalPair->GetPType() == 0) {
        NFIntegral newNFIntegral(theFinalPair);
        nfIntegralValue = newNFIntegral(this->GetJGroup(jGroupNum)->GetChannel(ch)->GetL(), theFinalPair->SpinOrbitTwoJ(jValue), exitPair->GetHeavyExE());
        WhitFunc newWhitFunc(theFinalPair);
        double whitConv = newWhitFunc(this->GetJGroup(jGroupNum)->GetChannel(ch)->GetL(), theFinalPair->SpinOrbitTwoJ(jValue),
                                      theFinalPair->GetChRad(),
                                      fabs(exitPair->GetHeavyExE() - theFinalPair->GetSepE() - theFinalPair->GetExE()));
        ecConvert = sqrt(2.0 * theFinalPair->GetRedMass() * theFinalPair->GetChRad() * uconv / pow(hbarc, 2.0)) / whitConv;
      }
      this->GetJGroup(jGroupNum)->GetLevel(levelNum)->AddNFIntegral(nfIntegralValue);
      this->GetJGroup(jGroupNum)->GetLevel(levelNum)->AddECConversionFactor(ecConvert);
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
  this->SortPathways(configure);
  this->BuildPathways(configure);
  this->SetUpFormalCaptureMasks(configure);
  if ((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_PATHWAYS)
    WritePathwaysCheck(*this, configure);

  // Calculate Angular Distribution Coefficients
  if (!(configure.paramMask & Config::USE_API))
    configure.outStream << "Calculating Angular Distribution Coefficients..." << std::endl;
  WarnTruncatedLegendreOrders(*this, configure);
  this->CalcAngularDists(configure.maxLOrder);
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

void CNuc::FillCompoundFromParamsPhysical(const vector_r &p) {
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      ALevel *level = this->GetJGroup(j)->GetLevel(la);
      level->SetE(p[i]);
      i++;
      std::vector<double> g;
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) {
        level->SetGamma(ch, p[i]);
        g.push_back(p[i]);
        i++;
      }
      level->SetSqrtNFFactor(SqrtNFFactor(*this, this->GetJGroup(j), level, 1.0, g));
    }
  }
  if (hasAbundanceBlock_)
    for (unsigned int s = 0; s < samples_.size(); s++)
      for (unsigned int k = 0; k < nuclei_.size(); k++) samples_[s].entries[k].abundance = p[i++];
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
    const CNucSystem structure = BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure);
    const phys::System &system = structure.system;
    const CNucSystemMap &map = structure.map;
    const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, map, configure);

    // 1. Observed widths, ANCs and amplitudes -> Brune or Park amplitudes (Eqs. 2.1-2.3).
    for (std::size_t b = 0; b < system.num_blocks(); b++) {
      const int j = map.groupOfBlock[b];
      JGroup *theJGroup = this->GetJGroup(j);
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
      JGroup *theJGroup = this->GetJGroup(j);
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
      JGroup *theJGroup = this->GetJGroup(j);
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
 * Calculates internal and external reaction pathways.
 */

void CNuc::SortPathways(const Config &configure) {
  // The enumeration is physics/ reaction_pathways (docs/upgrade/cnuc_split.md section 10); the trees on PPair are
  // filled from it in its order, through the System map (block -> J-group, block channel -> group channel).
  namespace phys = azure2::physics;
  const CNucSystem structure = BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure);
  const CNucSystemMap &map = structure.map;
  // The external-capture final states, in J-group and level order.
  std::vector<phys::FinalState> finalStates;
  std::vector<std::pair<int, int> > finalLevel;  // (J-group, level) of each
  for (const azure2::PathwayFinalState &f : ECFinalStates(*this)) {
    finalStates.push_back(f.state);
    finalLevel.push_back(std::make_pair(f.group, f.level));
  }
  phys::PathwayOptions options;
  options.rmatrix_capture = (configure.paramMask & Config::USE_RMC_FORMALISM) != 0;

  for (int aa = 1; aa <= this->NumPairs(); aa++) {
    PPair *entrancePair = this->GetPair(aa);
    if (!entrancePair->IsEntrance()) continue;
    const std::vector<phys::ExitPathways> exits = phys::reaction_pathways(structure.system, aa - 1, finalStates, options);
    for (const phys::ExitPathways &exit : exits) {
      Decay NewDecay(map.pairOfPartition[exit.partition]);
      int decayNum = entrancePair->IsDecay(NewDecay);
      if (!decayNum) {
        entrancePair->AddDecay(NewDecay);
        decayNum = entrancePair->IsDecay(NewDecay);
      }
      Decay *theDecay = entrancePair->GetDecay(decayNum);
      for (const phys::SpinGroup &g : exit.groups) {
        KGroup NewKGroup = g.has_sp2 ? KGroup(g.two_s / 2., g.two_sp / 2., g.two_sp2 / 2.) : KGroup(g.two_s / 2., g.two_sp / 2.);
        int kGroupNum = theDecay->IsKGroup(NewKGroup, g.has_sp2);
        if (!kGroupNum) {
          theDecay->AddKGroup(NewKGroup);
          kGroupNum = theDecay->IsKGroup(NewKGroup, g.has_sp2);
        }
        KGroup *theKGroup = theDecay->GetKGroup(kGroupNum);
        for (const phys::InternalPathway &p : g.internal) {
          MGroup NewMGroup(map.groupOfBlock[p.block], map.channelOfBlockChannel[p.block][p.k_in],
                           map.channelOfBlockChannel[p.block][p.k_out]);
          int mGroupNum = theKGroup->IsMGroup(NewMGroup);
          if (!mGroupNum) {
            theKGroup->AddMGroup(NewMGroup);
            mGroupNum = theKGroup->IsMGroup(NewMGroup);
          }
          if (p.has_stat_spin_factor) theKGroup->GetMGroup(mGroupNum)->SetStatSpinFactor(p.stat_spin_factor);
        }
      }
    }
    // External pathways, after the internal ones of this entrance (they refer to them by position).
    for (std::size_t e = 0; e < exits.size(); e++) {
      Decay *theDecay = entrancePair->GetDecay(entrancePair->IsDecay(map.pairOfPartition[exits[e].partition]));
      for (std::size_t g = 0; g < exits[e].groups.size(); g++)
        for (const phys::ExternalPathway &x : exits[e].groups[g].external) {
          const char radType = (x.X == phys::Multipole::E) ? 'E' : 'M';
          const int chp = (int)x.final_channel + 1;
          const int j = finalLevel[x.final_state].first, la = finalLevel[x.final_state].second;
          KGroup *theKGroup = theDecay->GetKGroup((int)g + 1);
          if (!x.channel_capture) {
            theKGroup->AddECMGroup(ECMGroup(radType, x.L, x.l, x.two_J / 2., chp, j, la));
            continue;
          }
          const phys::InternalPathway &cc = exits[x.cc_exit].groups[x.cc_group].internal[x.cc_pathway];
          const int internalChannel = x.photon_channel < 0 ? 0 : map.channelOfBlockChannel[cc.block][x.photon_channel];
          theKGroup->AddECMGroup(ECMGroup(radType, x.L, x.l, x.two_J / 2., chp, j, la, (int)x.cc_exit + 1,
                                          (int)x.cc_group + 1, (int)x.cc_pathway + 1, internalChannel));
        }
    }
  }
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
    const CNucSystem structure = BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure);
    const azure2::physics::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      JGroup *theJGroup = this->GetJGroup(structure.map.groupOfBlock[b]);
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
 * Creates and sorts the KLGroup and Interference objects and calculates the appropriate coefficients.
 */

void CNuc::CalcAngularDists(int maxL) {
  for (int aa = 1; aa <= this->NumPairs(); aa++) {
    PPair *entrancePair = this->GetPair(aa);
    if (entrancePair->GetPType() == 20) continue;
    for (int ir = 1; ir <= this->GetPair(aa)->NumDecays(); ir++) {
      Decay *theDecay = this->GetPair(aa)->GetDecay(ir);
      for (int k = 1; k <= theDecay->NumKGroups(); k++) {
        for (int lOrder = 0; lOrder <= maxL; lOrder++) {
          for (int m1 = 1; m1 <= theDecay->GetKGroup(k)->NumMGroups() + theDecay->GetKGroup(k)->NumECMGroups(); m1++) {
            for (int m2 = 1; m2 <= theDecay->GetKGroup(k)->NumMGroups() + theDecay->GetKGroup(k)->NumECMGroups(); m2++) {
              std::string interferenceType;
              double j1, j2, l1, l1p, l2, l2p;
              double s1pCh = -1., s2pCh = -1.;  // exit channel spins of internal pathways
              int w1p, w2p, path1, path2;
              if (m1 > theDecay->GetKGroup(k)->NumMGroups()) {
                int m1_ec = m1 - theDecay->GetKGroup(k)->NumMGroups();
                ECMGroup *theECMGroup1 = theDecay->GetKGroup(k)->GetECMGroup(m1_ec);
                j1 = theECMGroup1->GetJ();
                l1 = (double)theECMGroup1->GetL();
                l1p = (double)theECMGroup1->GetMult();
                if (theECMGroup1->GetRadType() == 'M')
                  w1p = 0;
                else
                  w1p = 1;
                interferenceType = 'E';
                path1 = m1_ec;
              } else {
                JGroup *jgroup1 = this->GetJGroup(theDecay->GetKGroup(k)->GetMGroup(m1)->GetJNum());
                AChannel *channel1 = jgroup1->GetChannel(theDecay->GetKGroup(k)->GetMGroup(m1)->GetChNum());
                AChannel *channel1p = jgroup1->GetChannel(theDecay->GetKGroup(k)->GetMGroup(m1)->GetChpNum());
                j1 = jgroup1->GetJ();
                l1 = (double)channel1->GetL();
                l1p = (double)channel1p->GetL();
                s1pCh = channel1p->GetS();
                if (channel1p->GetRadType() == 'M' || channel1p->GetRadType() == 'P')
                  w1p = 0;
                else
                  w1p = 1;
                interferenceType = 'R';
                path1 = m1;
              }
              if (m2 > theDecay->GetKGroup(k)->NumMGroups()) {
                int m2_ec = m2 - theDecay->GetKGroup(k)->NumMGroups();
                ECMGroup *theECMGroup2 = theDecay->GetKGroup(k)->GetECMGroup(m2_ec);
                j2 = theECMGroup2->GetJ();
                l2 = (double)theECMGroup2->GetL();
                l2p = (double)theECMGroup2->GetMult();
                if (theECMGroup2->GetRadType() == 'M')
                  w2p = 0;
                else
                  w2p = 1;
                interferenceType += 'E';
                path2 = m2_ec;
              } else {
                JGroup *jgroup2 = this->GetJGroup(theDecay->GetKGroup(k)->GetMGroup(m2)->GetJNum());
                AChannel *channel2 = jgroup2->GetChannel(theDecay->GetKGroup(k)->GetMGroup(m2)->GetChNum());
                AChannel *channel2p = jgroup2->GetChannel(theDecay->GetKGroup(k)->GetMGroup(m2)->GetChpNum());
                j2 = jgroup2->GetJ();
                l2 = (double)channel2->GetL();
                l2p = (double)channel2p->GetL();
                s2pCh = channel2p->GetS();
                if (channel2p->GetRadType() == 'M' || channel2p->GetRadType() == 'P')
                  w2p = 0;
                else
                  w2p = 1;
                interferenceType += 'R';
                path2 = m2;
              }
              double s = theDecay->GetKGroup(k)->GetS();
              double sp = theDecay->GetKGroup(k)->GetSp();
              if ((int)(l1 + l2 + lOrder) % 2 == 0 && (int)(l1p + l2p + w1p + w2p + lOrder) % 2 == 0) {
                double z1z2 = 0.0;
                if (this->GetPair(theDecay->GetPairNum())->GetPType() == 0) {
                  z1z2 = azure2::physics::particle_interference((int)l1, (int)l1p, Twice(j1), (int)l2, (int)l2p, Twice(j2),
                                                                Twice(s), Twice(sp), lOrder);
                  // A cross-spin UPOS group mixes two exit channel spins; only the
                  // secondary-gamma coefficient below applies to it.
                  if (theDecay->GetKGroup(k)->IsCrossSpin()) z1z2 = 0.;
                } else if (this->GetPair(theDecay->GetPairNum())->GetPType() == 10) {
                  double jf = this->GetPair(theDecay->GetPairNum())->GetJ(2);
                  z1z2 = azure2::physics::photon_interference((int)l1, (int)l1p, Twice(j1), (int)l2, (int)l2p, Twice(j2),
                                                              Twice(s), Twice(jf), lOrder);
                }
                // Calculate z1z2_upos for Unobserved Primary, Observed Secondary (UPOS) case
                double z1z2_upos = 0.;
                // Eq. (22) of Brune and deBoer, PRC 102, 024628 (2020): the first
                // amplitude has exit channel spin s' = sp, the second s2' = sp2.  Eq. (22)
                // as printed corresponds to coupling the exit channel spin as
                // (I_B I_b | s'); AZURE2 couples the light particle first, (I_b I_B | s'),
                // as Eq. (17) of the paper and the analyzing-power code (PolarizationFunc)
                // do, which multiplies each amplitude by (-1)^(I_b+I_B-s') and the product
                // by (-1)^(sp-sp2).
                const double sp2 = theDecay->GetKGroup(k)->GetSp2();
                if (this->GetPair(theDecay->GetPairNum())->GetPType() == 0 && aa != theDecay->GetPairNum() &&
                    s1pCh == sp && s2pCh == sp2) {
                  double j1f = this->GetPair(theDecay->GetPairNum())->GetJ(1);
                  double j2f = this->GetPair(theDecay->GetPairNum())->GetJ(2);
                  z1z2_upos = azure2::physics::upos_interference((int)l1, (int)l1p, Twice(j1), (int)l2, Twice(j2), Twice(s),
                                                                 Twice(sp), Twice(sp2), Twice(j1f), Twice(j2f), lOrder);
                }
                if (fabs(z1z2) > 1e-10 || fabs(z1z2_upos) > 1e-10) {
                  KLGroup NewKLGroup(k, lOrder);
                  int KLGroupNum = theDecay->IsKLGroup(NewKLGroup);
                  if (!KLGroupNum) {
                    theDecay->AddKLGroup(NewKLGroup);
                    KLGroupNum = theDecay->IsKLGroup(NewKLGroup);
                  }
                  if (aa == theDecay->GetPairNum() || this->GetPair(theDecay->GetPairNum())->GetPType() == 10) {
                    // Standard interference (elastic or gamma capture)
                    Interference NewInterference(path1, path2, z1z2, interferenceType);
                    int InterNum = theDecay->GetKLGroup(KLGroupNum)->IsInterference(NewInterference);
                    if (!InterNum) {
                      theDecay->GetKLGroup(KLGroupNum)->AddInterference(NewInterference);
                      InterNum = theDecay->GetKLGroup(KLGroupNum)->IsInterference(NewInterference);
                    }
                  } else {
                    // UPOS interference (particle exit, different entrance/exit pair)
                    Interference NewInterference(path1, path2, z1z2, z1z2_upos, interferenceType);
                    int InterNum = theDecay->GetKLGroup(KLGroupNum)->IsInterference(NewInterference);
                    if (!InterNum) {
                      theDecay->GetKLGroup(KLGroupNum)->AddInterference(NewInterference);
                      InterNum = theDecay->GetKLGroup(KLGroupNum)->IsInterference(NewInterference);
                    }
                  }
                }
              }
            }
          }
        }
      }
      theDecay->BuildFlatInterferences();
      std::vector<DiagonalPathway> diag;
      for (int k = 1; k <= theDecay->NumKGroups(); k++)
        for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
          MGroup *mg = theDecay->GetKGroup(k)->GetMGroup(m);
          AChannel *en = this->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChNum());
          AChannel *ex = this->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChpNum());
          if (en == ex) diag.push_back(DiagonalPathway{k, m, en->GetS(), en->GetL(), mg->GetStatSpinFactor()});
        }
      theDecay->SetDiagonalPathways(diag);
      theDecay->SetTablesBuilt();
    }
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
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      ALevel *level = this->GetJGroup(j)->GetLevel(la);
      varname = "energy_" + std::to_string(energyIndex);
      params.Add(varname, level->GetE(), 0.1 * level->GetE(), ParameterRecord::LevelEnergy(j, la, energyIndex));
      bool isUnbound = false;
      for (int ir = 1; ir <= this->NumPairs(); ir++) {
        PPair *pair = this->GetPair(ir);
        if (pair->GetNucleus() != this->GetJGroup(j)->GetNucleus()) continue;   // thresholds of this nucleus only
        if (pair->GetPType() == 0 &&
            level->GetE() > (pair->GetSepE() + pair->GetExE())) isUnbound = true;
      }
      // Jakub's Fix: it can create problems with the fit with pyazr and brick
      if (!isUnbound) p.Fix(varname);
      if (level->EnergyFixed() && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);

      // Parameter settings will be applied by ParameterLimitsManager during fit
      int widthIndex = 1;
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) {
        varname = "width_" + std::to_string(energyIndex) + "_" + std::to_string(widthIndex);
        widthIndex++;
        params.Add(varname, level->GetGamma(ch), 0.1 * level->GetGamma(ch), ParameterRecord::Width(j, la, ch, energyIndex));
        // A width given as 0 is not fitted, in every basis: the input decides, not the transformed amplitude, which
        // in the formal basis is 0 or rounding by the order of the sums (external_capture.md section 9).  Without a
        // transformation the input is the parameter.
        if (this->IsTransformedIn() ? level->ZeroInputWidth(ch) : level->GetGamma(ch) == 0.0) p.Fix(varname);
        if (level->ChannelFixed(ch) && !p.Parameter(p.Index(varname)).IsFixed()) p.Fix(varname);
        // Apply Wigner Limit bounds if flag is enabled.  GetWignerLimit()
        // returns gamma_W^2 (MeV); the fit parameter at this point is the
        // reduced-width AMPLITUDE gamma (MeV^1/2), so the bound must be
        // sqrt(gamma_W^2), not gamma_W^2 itself.
        if (config && (config->paramMask & Config::USE_WIGNER_LIMITS)) {
          AChannel *channel = this->GetJGroup(j)->GetChannel(ch);
          // GetWignerLimit returns gamma^2_W = 3 hbar^2/(2 mu a^2) in MeV -- the
          // limit on the reduced width *squared*.  The fit parameter is the
          // reduced width *amplitude* in MeV^(1/2), so the bound is its square
          // root; using gamma^2_W directly bounds a MeV^(1/2) quantity by a MeV
          // number and lets theta^2 = gamma^2/gamma^2_W run past 1 by a factor
          // of gamma^2_W itself.
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
  // Isotopic abundances, one per compound nucleus and sample (sample by
  // sample), right after the level block (so EData::FillMnParams, which
  // takes its offsets from the vector size, needs no change); only when the
  // project carries an <abundances> block.  A fixed abundance is still a
  // parameter, so the layout does not depend on which ones vary.
  if (hasAbundanceBlock_) {
    for (unsigned int s = 0; s < samples_.size(); s++) {
      for (unsigned int k = 0; k < nuclei_.size(); k++) {
        const SampleAbundance &e = samples_[s].entries[k];
        std::string name = this->AbundanceParamName((int)s + 1, (int)k + 1);
        params.Add(name, e.abundance, e.sigma > 0.0 ? e.sigma : 0.1 * e.abundance + 1e-3,
                   ParameterRecord::Abundance((int)s + 1, samples_[s].name, (int)k + 1, nuclei_[k].Z, nuclei_[k].A));
        p.SetLowerLimit(name, 0.0);
        if (!e.vary) p.Fix(name);
      }
    }
  }
}

/*!
 * Fills the CNuc object from the Minuit parameter array.
 */

void CNuc::FillCompoundFromParams(const vector_r &p) {
  int i = 0;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      ALevel *level = this->GetJGroup(j)->GetLevel(la);
      level->SetFitE(p[i]);
      i++;
      std::vector<double> g;
      for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++) {
        level->SetFitGamma(ch, p[i]);
        level->SetECGamma(ch, p[i]);
        g.push_back(p[i]);
        i++;
      }
      level->SetSqrtNFFactor(SqrtNFFactor(*this, this->GetJGroup(j), level, 1.0, g));
    }
  }
  if (hasAbundanceBlock_)
    for (unsigned int s = 0; s < samples_.size(); s++)
      for (unsigned int k = 0; k < nuclei_.size(); k++) samples_[s].entries[k].abundance = p[i++];
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
      for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
        ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
        if (report && theLevel->IsInRMatrix() && theLevel->GetParkNorm() <= 0.0) WarnParkNormNotPositive(*this, configure, j, la);
      }
  }
  if (!(configure.paramMask & Config::USE_BRUNE_FORMALISM)) {
    // Formal -> Brune by the fixed-point iteration of physics/ brune_from_formal (level_parameters.md section 4), on
    // the fit values.  A level not in the R-matrix keeps its values.
    const phys::FormalToBruneOptions options;  // the engine's 1e-6 MeV and 1000 iterations
    for (int j = 1; j <= this->NumJGroups(); j++)
      for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
        ALevel *theLevel = this->GetJGroup(j)->GetLevel(la);
        if (theLevel->IsInRMatrix() && this->GetJGroup(j)->IsInRMatrix()) continue;
        theLevel->SetTransformE(theLevel->GetFitE());
        theLevel->SetTransformIterations(0);
        for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++)
          theLevel->SetTransformGamma(ch, theLevel->GetFitGamma(ch));
      }
    const CNucSystem structure = BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure);
    const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
    for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
      const int j = structure.map.groupOfBlock[b];
      JGroup *theJGroup = this->GetJGroup(j);
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
      for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
        this->GetJGroup(j)->GetLevel(la)->SetTransformIterations(0);
        this->GetJGroup(j)->GetLevel(la)->SetTransformE(this->GetJGroup(j)->GetLevel(la)->GetFitE());
        for (int ch = 1; ch <= this->GetJGroup(j)->NumChannels(); ch++)
          this->GetJGroup(j)->GetLevel(la)->SetTransformGamma(ch, this->GetJGroup(j)->GetLevel(la)->GetFitGamma(ch));
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
  const CNucSystem all = BuildSystemFromCNuc(*this, configure, CNucSystemContent::AllGroups);
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, all.map, configure);
  for (std::size_t b = 0; b < all.system.num_blocks(); b++) {
    JGroup *theJGroup = this->GetJGroup(all.map.groupOfBlock[b]);
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
    JGroup *theJGroup = this->GetJGroup(j);
    for (int la = 1; la <= theJGroup->NumLevels(); la++) {
      ALevel *theLevel = theJGroup->GetLevel(la);
      params.push_back(theLevel->GetTransformE());
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) params.push_back(theLevel->GetObservedValue(ch));
    }
  }
  // Abundances are the same in both bases; they follow the level block as in FillMnParams.
  if (hasAbundanceBlock_)
    for (unsigned int s = 0; s < samples_.size(); s++)
      for (unsigned int k = 0; k < nuclei_.size(); k++) params.push_back(samples_[s].entries[k].abundance);

  return params;
}

namespace {

/*!
 * Park's numerical shift-function derivatives by their exact arguments (the
 * pair's parameters, l and the level energy, compared bit for bit, and the
 * settings that change the shift function at those arguments: GSL or own
 * Coulomb functions, a nuclear potential, the Coulomb table, relativistic
 * kinematics -- pyazr sessions with different settings share one process), shared by
 * all compound nuclei of the process: the working copies of a fit (evaluation,
 * gradient, residual Jacobian) need the same values at the same parameters.
 * An entry is what the computation returns, so a hit changes nothing.
 * Cleared by CNuc::InvalidateShiftFunctions (a new nuclear potential) and
 * when it grows past kShiftDerivMemoMax.
 */
struct ShiftDerivKey {
  int z1, z2, l, two_j;  ///< two_j: PPair::SpinOrbitTwoJ
  bool gsl;
  std::vector<double> hybrid;  ///< azure2::physics::PotentialKey of the pair's potential
  bool table, relativistic;
  double redmass, mass, radius, threshold, energy;
  bool operator<(const ShiftDerivKey &o) const {
    return std::tie(z1, z2, l, two_j, gsl, hybrid, table, relativistic, redmass, mass, radius, threshold, energy) <
        std::tie(o.z1, o.z2, o.l, o.two_j, o.gsl, o.hybrid, o.table, o.relativistic, o.redmass, o.mass, o.radius, o.threshold,
                 o.energy);
  }
};
struct ShiftDerivValue {
  bool has1 = false, has2 = false;
  double d1 = 0.0, d2 = 0.0;
};
std::mutex g_shiftDerivMutex;
std::map<ShiftDerivKey, ShiftDerivValue> g_shiftDerivMemo;
const std::size_t kShiftDerivMemoMax = 1 << 16;

}  // namespace

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
 * them through an exact-argument memo (g_shiftDerivMemo).
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
    if (this->GetJGroup(j)->IsInRMatrix()) {
      JGroup *theJGroup = this->GetJGroup(j);
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
            PPair *thePair = this->GetPair(theChannel->GetPairNum());
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
      const ShiftDerivKey key{d.pair->GetZ(1), d.pair->GetZ(2), d.l, d.two_j, gslCoulomb, azure2::physics::PotentialKey(theCoulombFunction.HybridPotential()),
                              theCoulombFunction.UsesTable(), Kinematics::Relativistic(), d.pair->GetRedMass(),
                              d.pair->GetM(1) + d.pair->GetM(2), d.pair->GetChRad(), threshold, d.energy};
      bool need1 = d.first, need2 = d.second;
      {
        std::lock_guard<std::mutex> lock(g_shiftDerivMutex);
        auto it = g_shiftDerivMemo.find(key);
        if (it != g_shiftDerivMemo.end()) {
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
      std::lock_guard<std::mutex> lock(g_shiftDerivMutex);
      if (g_shiftDerivMemo.size() >= kShiftDerivMemoMax) g_shiftDerivMemo.clear();
      ShiftDerivValue &v = g_shiftDerivMemo[key];
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
        if (theChannel->GetRadType() == 'P' && this->GetPair(theChannel->GetPairNum())->GetPType() == 0)
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
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++) {
      energyIndex[{j, la}] = i;
      i += 1 + this->GetJGroup(j)->NumChannels();
    }
  // physics/ park_from_brune, brune_from_park on each R-matrix block (level_parameters.md Eq. 2.3); the errors
  // scale with the amplitudes, by J^1/2 or J^-1/2.
  const CNucSystem structure = BuildSystemFromCNuc(*this, configure, CNucSystemContent::Structure);
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(*this, structure.map, configure);
  for (std::size_t b = 0; b < structure.system.num_blocks(); b++) {
    const int j = structure.map.groupOfBlock[b];
    const std::vector<int> &levels = structure.map.levelOfBlockLevel[b];
    const int nc = this->GetJGroup(j)->NumChannels();
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
 * Forgets the cached shift functions of every level.
 */

void CNuc::InvalidateShiftFunctions() {
  {
    std::lock_guard<std::mutex> lock(g_shiftDerivMutex);
    g_shiftDerivMemo.clear();
  }
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      this->GetJGroup(j)->GetLevel(la)->InvalidateShiftFunctions();
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
    JGroup *theJGroup = this->GetJGroup(j);
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

PPair *CNuc::GetPair(int pairNum) {
  PPair *b = &pairs_[pairNum - 1];
  return b;
}

const PPair *CNuc::GetPair(int pairNum) const {
  return &pairs_[pairNum - 1];
}

PPair *CNuc::GetPairFromKey(int key) {
  const int num = GetPairNumFromKey(key);
  return num > 0 ? GetPair(num) : nullptr;
}

/*!
 * Returns a pointer to the \f$ J^\pi \f$ group specified by a position in the JGroup vector.
 */

JGroup *CNuc::GetJGroup(int jGroupNum) {
  JGroup *b = &jgroups_[jGroupNum - 1];
  return b;
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
