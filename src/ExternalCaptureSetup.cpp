// The external-capture amplitudes of the data points: their computation at every point and sub-point, or their
// reading from the intEC file of an earlier run, the count and the signature that decide whether that file belongs
// to this calculation, and the progress bar of the amplitude loops.
#include <atomic>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "CompoundSystem.h"
#include "Config.h"
#include "CoulombTable.h"
#include "ECTable.h"
#include "ExternalCaptureExit.h"
#include "WorkingCopy.h"
#include "GSLException.h"
#include "azure2/support/FloatingPointPrecision.h"

const azure2::ExitRow &ExternalCaptureExit(const CompoundSystem &compound, int entrancePair, int exitPair, int segmentKey) {
  const azure2::EntrancePathways &P = compound.Pathways().Entrance(entrancePair);
  const int exit = P.ExitOf(exitPair);
  if (!exit)
    throw std::runtime_error("data segment #" + std::to_string(segmentKey) + " asks for external capture from pair " +
                             std::to_string(compound.GetPair(entrancePair)->GetPairKey()) + " to pair " +
                             std::to_string(compound.GetPair(exitPair)->GetPairKey()) +
                             ", but no reaction pathway leads from the one to the other");
  return P.exits[exit - 1];
}

namespace {

/*!
 * 64-bit FNV-1a over a stream of values, each written as text.  Doubles are
 * printed at fixed precision (%.12e) rather than hashed bit for bit, so the
 * signature survives a value that went through a text file and back, while
 * still moving for any change that matters to the integrals.
 */
class ECSignatureHash {
 public:
  void Add(const char *text) {
    for (const char *c = text; *c; c++) {
      hash_ ^= (unsigned char)*c;
      hash_ *= 1099511628211ULL;
    }
    // A separator, so that ("1","23") and ("12","3") hash differently.
    hash_ ^= (unsigned char)';';
    hash_ *= 1099511628211ULL;
  }
  void Add(double value) {
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "%.12e", value);
    Add(buffer);
  }
  void Add(int value) {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%d", value);
    Add(buffer);
  }
  std::string Hex() const {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%016" PRIx64, hash_);
    return buffer;
  }

 private:
  uint64_t hash_ = 14695981039346656037ULL;
};

/*!
 * Walks the external-capture amplitudes of one set of segments in exactly the
 * order CalculateECAmplitudes and InitializeLeaves write them, and
 * returns how many there are.  With a hash, it also feeds in what each block of
 * amplitudes depends on, and collects the pairs involved so their properties
 * can be added afterwards.
 */
/// A leaf seen as the walk sees a segment: its pairs, its key and its points.
struct LeafView {
  const Segment *segment;
  int term;
  int GetEntranceKey() const { return segment->TermEntranceKey(term); }
  int GetExitKey() const { return segment->TermExitKey(term); }
  int GetSegmentKey() const { return segment->TermKey(term); }
  UnitRange<Unit> GetPoints(WorkingCopy &data) const { return segment->TermPoints(data, term); }
};

/// The leaves of `data` in the order InitializeLeaves writes their amplitudes.
std::vector<LeafView> LeafViews(WorkingCopy &data) {
  std::vector<LeafView> views;
  for (const WorkingCopy::LeafRef &leaf : data.Data().LeafRefs()) views.push_back({&data.Data().GetSegments()[leaf.first], leaf.second});
  return views;
}

template <class Segments>
long long WalkECAmplitudes(WorkingCopy &data, Segments &segments, const CompoundSystem *compound, ECSignatureHash *hash,
                           std::set<int> *pairs) {
  long long count = 0;
  for (auto segment = segments.begin(); segment != segments.end(); segment++) {
    int aa = compound->GetPairNumFromKey(segment->GetEntranceKey());
    if (compound->GetPair(aa)->GetPType() == 20) continue;
    if (!compound->GetPair(aa)->IsEntrance()) continue;
    const azure2::EntrancePathways &P = compound->Pathways().Entrance(aa);
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      for (int la = 1; la <= compound->GetJGroup(j)->NumLevels(); la++) {
        if (!compound->GetJGroup(j)->GetLevel(la)->IsECLevel()) continue;
        const ALevel *ecLevel = compound->GetJGroup(j)->GetLevel(la);
        int ir = compound->GetPairNumFromKey(segment->GetExitKey());
        if (ecLevel->GetECPairNum() != ir) continue;
        const azure2::ExitRow &exit = ExternalCaptureExit(*compound, aa, ir, segment->GetSegmentKey());
        if (hash) {
          // Which block this is, and the final state it captures to.  The
          // level energy enters every integral (Unit::CalculateECAmplitudes);
          // the ANC does not -- it multiplies the amplitude later -- so it is
          // deliberately left out and can be varied without invalidating.
          hash->Add("segment");
          hash->Add(segment->GetSegmentKey());
          hash->Add(segment->GetEntranceKey());
          hash->Add(segment->GetExitKey());
          hash->Add(compound->GetJGroup(j)->GetJ());
          hash->Add(ecLevel->GetECEnergy());  // the energy the integrals use (Brune; cnuc_split.md section 8)
          pairs->insert(aa);
          // The pathways in the table's order, which is the order of the trees the signature was defined on: a file
          // written before the table keeps its signature.
          for (int k = 1; k <= exit.num_groups; k++) {
            const azure2::SpinGroupRow &group = P.Group(exit, k);
            hash->Add(group.s);
            for (int ecm = 1; ecm <= group.num_external; ecm++) {
              const azure2::ExternalRow &pathway = P.External(group, ecm);
              const AChannel *finalChannel = compound->GetJGroup(j)->GetChannel(pathway.final_channel);
              hash->Add(pathway.l);
              hash->Add(pathway.mult);
              hash->Add((int)pathway.rad_type);
              hash->Add(pathway.J);
              hash->Add(finalChannel->GetL());
              hash->Add(finalChannel->GetS());
              hash->Add(compound->GetPair(finalChannel->GetPairNum())->GetPairKey());
              pairs->insert(finalChannel->GetPairNum());
              // A channel's own radius (<channelRadii>) only where it is not its pair's, which the pairs' part of the
              // signature holds: a file written before the records keeps its signature.
              {
                const azure2::InternalRow &through = P.internal[pathway.through_row];
                const AChannel *entranceChannel = compound->GetJGroup(through.group)->GetChannel(through.ch);
                const double finalRadius = finalChannel->GetRadius(), entranceRadius = entranceChannel->GetRadius();
                if (finalRadius != compound->GetPair(finalChannel->GetPairNum())->GetPairRadius() ||
                    entranceRadius != compound->GetPair(entranceChannel->GetPairNum())->GetPairRadius()) {
                  hash->Add("channel radii");
                  hash->Add(finalRadius);
                  hash->Add(entranceRadius);
                }
              }
              // An unbound final state is a narrow resonance (ECIntegral::IntegrateUnbound): marked only then, so
              // every other file keeps its signature.
              {
                const PPair *fp = compound->GetPair(finalChannel->GetPairNum());
                if (ecLevel->GetECEnergy() > fp->GetSepE() + fp->GetExE()) {
                  hash->Add("unbound final state v1");
                  // with the capture from the potential's interior (ECIntegral::Interior) since 2026-11-21
                  const auto &potential = fp->GetHybridPotential(finalChannel->GetL());
                  if (potential && potential->interior_boundary) hash->Add("unbound final state, interior v1");
                }
              }
              hash->Add(pathway.channel_capture ? 1 : 0);
              if (pathway.channel_capture) {
                const azure2::InternalRow &chanCap = P.internal[pathway.cc_row];
                const AChannel *initial = compound->GetJGroup(chanCap.group)->GetChannel(chanCap.chp);
                hash->Add(initial->GetL());
                hash->Add(initial->GetS());
              }
            }
          }
        }
        // Every point carries its own amplitudes, mapped or not (see
        // CalculateECAmplitudes), so the intEC file does not depend on mapping.
        for (auto point = segment->GetPoints(data).begin(); point != segment->GetPoints(data).end(); point++) {
          for (int k = 1; k <= exit.num_groups; k++)
            count += (long long)P.Group(exit, k).num_external * (1 + (long long)point->NumSubPoints());
          if (hash) {
            // The energies the amplitudes are evaluated at.  The angle does
            // not enter an EC amplitude, so it is not part of the signature.
            // The lab energy is included as well as the c.m. one because a
            // leaf's c.m. energy is converted again, with its own
            // entrance pair, after this point (InitializeLeaves).
            hash->Add("point");
            hash->Add(point->GetLabEnergy());
            hash->Add(point->GetCMEnergy());
            hash->Add(point->NumSubPoints());
            for (auto subPoint = point->GetSubPoints(data).begin(); subPoint != point->GetSubPoints(data).end(); subPoint++)
              hash->Add(subPoint->GetCMEnergy());
          }
        }
      }
    }
  }
  return count;
}

}  // namespace

/*!
 * Counts the external-capture amplitudes this model expects to read from, or
 * write to, an intEC file: those of the data segments, which
 * CalculateECAmplitudes writes, followed by those of the leaves,
 * which InitializeLeaves appends.  The walk mirrors both loops
 * exactly; if either changes, WalkECAmplitudes must change with it.
 *
 * A count only catches a file built for a different number of points.  The
 * signature (ECSignature) is what catches the same number of points at
 * different energies.
 */

long long WorkingCopy::CountECAmplitudes(const CompoundSystem *compound, const Config &configure) {
  std::vector<LeafView> leaves = LeafViews(*this);
  return WalkECAmplitudes(*this, Data().GetSegments(), compound, nullptr, nullptr) + WalkECAmplitudes(*this, leaves, compound, nullptr, nullptr);
}

/*!
 * The signature of the external-capture integrals this calculation needs.
 *
 * An intEC file records amplitudes at the sub-point energies of one grid and
 * carries no record of which grid.  Energy straggling, a target thickness or
 * the adaptive-grid settings move those energies without changing how many
 * there are; a channel radius, the Coulomb-function routine or the hybrid
 * potential change the integrand at the same energies.  Either way the count
 * check passes and the stale amplitudes are reused, silently.  The signature
 * covers all of it: every energy an amplitude is evaluated at, the structure
 * of each capture block, the level energy of each final state, the entrance
 * and final pairs (masses, charges, spins, separation and excitation
 * energies, channel radius, hybrid potential), and the options that select
 * how the Coulomb functions and the integrals are computed.
 *
 * Energy shifts: the amplitudes are computed once, here, at the energies the
 * points have at initialization, and are not recomputed when a fit moves a
 * segment's energy shift (Segment::UpdatePointEnergiesWithShift refreshes the
 * penetrabilities and phases, not the EC amplitudes).  The file is therefore
 * the initial state only, and the signature is taken over the energies at the
 * time the file is written or read, which is the same stage of the same
 * calculation in both cases.
 */

std::string WorkingCopy::ECSignature(const CompoundSystem *compound, const Config &configure) {
  ECSignatureHash hash;
  // Bump the version if what goes into the signature changes, so files written
  // under the old definition are recomputed rather than wrongly matched.
  // v2 (2026-10-23): the long-wavelength switch is part of the signature, and
  // the M1 amplitudes changed (sqrt(16/3) prefactor, no E1 charge function), so
  // every v1 file is recomputed once.
  // v3 (2026-10-26): the M1 spin-term sign for a cluster whose partner has
  // half-integer spin (ECIntegral.cpp), so every v2 file is recomputed once.
  // v4 (2026-10-26): the M1 orbital term is dropped on channel-spin-flip pathways (ECIntegral.cpp), so every
  // v3 file with such pathways is recomputed once.
  // v5 (2026-10-30): the file holds every digit (17 significant, was 6), so every v4 file is recomputed once.
  // v6 (2026-11-09): the Whittaker function of the final state is physics/coulomb's whittaker, not GSL's U
  // (coulomb.md section 1.4); the integrals move at ~1e-10 relative, so every v5 file is recomputed once.
  // v7 (2026-11-16): the Clebsch-Gordan and Racah coefficients of the amplitudes are those of physics/, not GSL's
  // 3j and 6j (AngCoeff.cpp); the integrals move at ~1e-14 relative, so every v6 file is recomputed once and a
  // result does not depend on whether a file from before was found.
  // v8 (2026-11-16): F of a neutral entrance pair (eta = 0) is the Riccati-Bessel function, no longer an upward
  // recurrence that was wrong for l >= 2 at rho << 1 (docs/upgrade/coulomb.md section 10); the integrals of a neutral
  // pair move (by ~1e-15 in the test projects), so every v7 file is recomputed once.
  // v9 (2026-11-18): the M1 amplitudes take the proton mass of CODATA 2018 (physics/ kProtonMass_MeV, 938.27208816 MeV,
  // 1.007276466617 u; was 1.00727638 u, 8.6e-8 lower; the owner's answer of 2026-11-18), so every v8 file is
  // recomputed once.
  hash.Add("AZURE2 EC signature v9");
  // The place of the option use_gsl_coul, retired on 2026-11-16 (docs/upgrade/coulomb.md section 9.1): a run with the
  // program's own Coulomb functions hashed 0 here, and they are the only ones now; the constant 0 is kept.
  hash.Add(0);
  hash.Add((configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX) ? 1 : 0);
  // The convection current changes exact-form electric terms only; long-wavelength signatures stay as in v2.
  if (configure.ecConvectionCurrent && !(configure.paramMask & Config::USE_LONGWAVELENGTH_APPROX))
    hash.Add("convection current");
  hash.Add((configure.paramMask & Config::CALCULATE_REACTION_RATE) ? 1 : 0);
  hash.Add(configure.useHybridMethod ? 1 : 0);
  // Relativistic kinematics (relativistic=1 in <config>) change eta, rho and the flux factor of every amplitude.
  // Added only when on, so every non-relativistic file keeps its signature and is not recomputed; a relativistic
  // run no longer reads back a file written without the switch (2026-10-27, checklist t3-25).
  if (configure.relativistic) hash.Add("relativistic");
  // The Coulomb waves of the integrals come from the Coulomb table when one is used and can be read (CoulFunc), else
  // from the exact functions.  The table is named by what it holds (CoulombTable::ContentIdentity: header and a hash
  // of its values), not by its file, so a copy of the same table elsewhere reuses this file, and a table with other
  // values, or none, does not (tests/ec_table_signature).
  const CoulombTable *table = configure.useCoulombTable ? CoulombTable::Shared(configure.coulombTablePath) : nullptr;
  hash.Add(table ? ("coulomb table " + table->ContentIdentity()).c_str() : "exact Coulomb functions");

  std::set<int> pairs;
  hash.Add("data segments");
  WalkECAmplitudes(*this, Data().GetSegments(), compound, &hash, &pairs);
  hash.Add("component segments");
  std::vector<LeafView> leaves = LeafViews(*this);
  WalkECAmplitudes(*this, leaves, compound, &hash, &pairs);

  for (std::set<int>::const_iterator it = pairs.begin(); it != pairs.end(); ++it) {
    const PPair *pair = compound->GetPair(*it);
    hash.Add("pair");
    hash.Add(pair->GetPairKey());
    hash.Add(pair->GetPType());
    for (int i = 1; i <= 2; i++) {
      hash.Add(pair->GetZ(i));
      hash.Add(pair->GetM(i));
      hash.Add(pair->GetJ(i));
      hash.Add(pair->GetPi(i));
      hash.Add(pair->GetG(i));
    }
    hash.Add(pair->GetSepE());
    hash.Add(pair->GetExE());
    // the split between the two particles (column 32) only when it is set, so an older cache keeps its signature
    if (pair->GetLightExE() != 0.) {
      hash.Add("light excitation");
      hash.Add(pair->GetLightExE());
    }
    hash.Add(pair->GetPairRadius());
    // The hybrid potential bends the Coulomb functions of this pair, and so
    // both the entrance-channel phases and the final-state integrand.
    const std::optional<azure2::physics::HybridPotential> &potential = pair->GetHybridPotential();
    hash.Add(potential ? 1 : 0);
    if (potential)
      for (double v : azure2::physics::PotentialKey(potential)) hash.Add(v);
    // potentials of single l (records with an l list): the rest of the pair's key, nothing for a pair without them
    const std::vector<double> all = pair->HybridPotentialKey(), base = azure2::physics::PotentialKey(potential);
    for (std::size_t i = base.size(); i < all.size(); i++) hash.Add(all[i]);
  }
  return hash.Hex();
}

std::string WorkingCopy::ECSignaturePath(const std::string &integralsFile) {
  return integralsFile + ".sig";
}

namespace {

// The 25-cell bar of the EC amplitude loops (CalculateECAmplitudes, InitializeLeaves) after `done` of
// `total` points: a cell per 4 %, and the percentage of the cells filled.
std::string ECProgressBar(int done, int total) {
  std::string progress = " [";
  double percent = 0.;
  for (int j = 1; j <= 25; j++) {
    if (done >= percent * total && percent < 1.) {
      percent += 0.04;
      progress += '*';
    } else
      progress += ' ';
  }
  std::ostringstream text;
  text << progress << "] " << percent * 100 << '%';
  return text.str();
}

}  // namespace

// Redraws the bar from inside the OpenMP loop of those two loops.  Every thread finishes points and wants to
// redraw, so the counter is atomic and the redraw runs under one named critical section (ec_progress, which the
// loops' error messages share): one thread at a time writes a whole bar, never interleaved with another's.  Up
// to dev 0119771 the threads wrote their bars piece by piece without a lock ("Segment #\r\tSegment #1111"),
// tests/progress_bars.
void WorkingCopy::ECProgressRedraw(const Config &configure, const char *segmentKeyOut, const std::atomic<int> &pointIndex,
                             int numPoints, time_t &startTime) {
#pragma omp critical(ec_progress)
  {
    if (difftime(time(NULL), startTime) > 0.25) {
      startTime = time(NULL);
      if (!(configure.paramMask & Config::USE_API)) {
        configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut << std::setw(0)
                            << ECProgressBar(pointIndex.load(), numPoints);
        configure.outStream.flush();
      }
    }
  }
}

/*!
 * If external capture amplitudes are to be calculated, Unit::CalculateECAmplitudes is
 * called for each point with a corresponding external capture component in the WorkingCopy object.
 * Otherwise, the amplitudes are read from the specified file.
 */

int WorkingCopy::CalculateECAmplitudes(const CompoundSystem *compound, const Config &configure) {
  std::ifstream in;
  std::ofstream out;
  std::string outputfile;
  if (configure.paramMask & Config::CALCULATE_WITH_DATA)
    outputfile = configure.outputdir + "intEC.dat";
  else
    outputfile = configure.outputdir + "intEC.extrap";

  // An intEC file records amplitudes at the sub-point energies of the grid it
  // was built for.  Reusing one built for a different grid, or for different
  // inputs at the same energies, returns amplitudes that belong to another
  // calculation -- silently, and with a result that looks like physics.  Two
  // checks guard the reuse: the number of amplitudes in the file, and the
  // signature recorded beside it (ECSignature).  If either disagrees, say so
  // and recompute rather than proceed.
  // The walk of the signature is the first to resolve the exit of every data segment and leaf with external
  // capture; an entrance without a pathway to the photon pair is refused here, before any amplitude is computed.
  try {
    ecSignature_ = ECSignature(compound, configure);
  } catch (const std::runtime_error &e) {
    configure.outStream << "ERROR: " << e.what() << "." << std::endl;
    return -1;
  }
  ecOutputFile_ = outputfile;
  bool usePrevious = (configure.paramMask & Config::USE_PREVIOUS_INTEGRALS) != 0;
  if (usePrevious) {
    long long expected = CountECAmplitudes(compound, configure);
    std::ifstream check(configure.OpenPath(configure.integralsfile).c_str());
    if (!check) {
      configure.outStream << "Could not open external capture file '"
                          << configure.integralsfile
                          << "'; the integrals will be recalculated." << std::endl;
      usePrevious = false;
    } else {
      long long found = 0;
      complex dummy(0.0, 0.0);
      while (check >> dummy) ++found;
      check.close();
      if (found != expected) {
        configure.outStream << "WARNING: '" << configure.integralsfile
                            << "' holds " << found << " external capture amplitudes but this"
                            << " calculation needs " << expected << "." << std::endl
                            << "         The file belongs to a different set of data segments or"
                            << " integration points." << std::endl
                            << "         Recalculating the integrals." << std::endl;
        usePrevious = false;
      }
    }
  }
  if (usePrevious) {
    std::string sigFile = ECSignaturePath(configure.integralsfile);
    std::ifstream sigIn(configure.OpenPath(sigFile).c_str());
    std::string recorded, line;
    while (sigIn && std::getline(sigIn, line)) {
      std::istringstream fields(line);
      std::string tag, version, value;
      if (fields >> tag >> version >> value && tag == "AZURE2-EC-SIGNATURE") recorded = version + " " + value;
    }
    if (recorded.empty()) {
      // A file from before signatures existed (or one copied without its
      // sidecar): there is no record of what it was computed for, and the
      // cases that matter are exactly the ones the count cannot see.
      configure.outStream << "WARNING: '" << configure.integralsfile
                          << "' has no signature file ('" << sigFile << "')," << std::endl
                          << "         so there is no record of the energies and inputs its"
                          << " integrals were computed for." << std::endl
                          << "         Recalculating the integrals." << std::endl;
      usePrevious = false;
    } else if (recorded != "1 " + ecSignature_) {
      configure.outStream << "WARNING: '" << configure.integralsfile
                          << "' was computed for a different integration grid or different inputs" << std::endl
                          << "         (signature " << recorded.substr(recorded.find(' ') + 1)
                          << ", this calculation needs " << ecSignature_ << ")." << std::endl
                          << "         The number of amplitudes matches, but their energies or integrands do not:"
                          << std::endl
                          << "         straggling, a target thickness, the integration-grid settings, a channel"
                          << std::endl
                          << "         radius, a level energy, the Coulomb-function routine or the hybrid potential"
                          << std::endl
                          << "         may have changed.  Recalculating the integrals." << std::endl;
      usePrevious = false;
    }
  }
  ecUsePrevious_ = usePrevious;
  ECTable::ResetCounters(configure);

  if (usePrevious)
    in.open(configure.OpenPath(configure.integralsfile).c_str());
  else {
    // The signature is written only once the file is complete, after the
    // leaves are appended (InitializeLeaves).  Remove
    // the old one first, so that a run interrupted part way leaves a file
    // with no signature -- which is recomputed -- rather than one whose old
    // signature still vouches for it.
    std::remove(configure.OpenPath(ECSignaturePath(outputfile)).c_str());
    out.open(configure.OpenPath(outputfile).c_str());
    if (!out) configure.outStream << "Could not write to EC Amplitude File." << std::endl;
    // Every digit (17 significant): the default 6 left a run that reads the file back ~1e-6 off the run that
    // wrote it (tests/13N chi2 by 0.62, 2026-10-30).
    out.precision(17);
  }
  int sumSegmentI = 0;
  int numSumSegments = 0;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (segment->IsTotalCapture()) {
      numSumSegments = segment->IsTotalCapture();
      sumSegmentI = 1;  // Start at 1 since total capture is now handled within the segment
    } else {
      sumSegmentI = 0;
      numSumSegments = 0;
    }
    char segmentKeyOut[256];
    if (segment->IsTotalCapture())
      snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d (Total: %d)", segment->GetSegmentKey(), numSumSegments);
    else
      snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d", segment->GetSegmentKey());
    int aa = compound->GetPairNumFromKey(segment->GetEntranceKey());
    if (compound->GetPair(aa)->GetPType() == 20) continue;
    if (compound->GetPair(aa)->IsEntrance()) {
      const azure2::EntrancePathways &P = compound->Pathways().Entrance(aa);
      for (int j = 1; j <= compound->NumJGroups(); j++) {
        for (int la = 1; la <= compound->GetJGroup(j)->NumLevels(); la++) {
          if (compound->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
            const ALevel *ecLevel = compound->GetJGroup(j)->GetLevel(la);
            int ir = compound->GetPairNumFromKey(segment->GetExitKey());
            if (ecLevel->GetECPairNum() == ir) {
              if (!usePrevious) {
                if (!(configure.paramMask & Config::USE_API)) {
                  configure.outStream << "\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [                         ] 0%";
                  configure.outStream.flush();
                }
                int numPoints = segment->NumPoints();
                std::atomic<int> pointIndex(0);
                time_t startTime = time(NULL);
                bool localStop = false;
#pragma omp parallel for shared(configure, localStop, pointIndex, startTime)
                for (int i = 1; i <= numPoints; i++) {
                  const azure2::support::ExtendedPrecisionScope fpPrecision;
                  if (configure.StopRequested() || localStop) continue;
                  Unit *point = segment->GetPoint(*this, i);
                  // Mapped points compute their own amplitudes too: the
                  // amplitudes are a one-off setup cost, and the intEC file
                  // then has the same layout with and without point mapping.
                  {
                    try {
                      point->CalculateECAmplitudes(this, compound, configure, channelTable_.get());
                    } catch (GSLException e) {
#pragma omp critical(ec_progress)
                      {
                        configure.outStream << e.what() << std::endl;
                        localStop = true;
                      }
                    }
                  }
                  ++pointIndex;
                  ECProgressRedraw(configure, segmentKeyOut, pointIndex, numPoints, startTime);
                }
                if (configure.StopRequested() || localStop) {
                  if (out.is_open()) out.close();
                  if (in.is_open()) in.close();
                  return -1;
                }
                if (!(configure.paramMask & Config::USE_API))
                  configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [*************************] 100%" << std::endl;
              }
              const azure2::ExitRow &exit = ExternalCaptureExit(*compound, aa, ir, segment->GetSegmentKey());
              for (UnitIterator point = segment->GetPoints(*this).begin();
                   point < segment->GetPoints(*this).end(); point++) {
                {
                  for (int k = 1; k <= exit.num_groups; k++) {
                    for (int ecm = 1; ecm <= P.Group(exit, k).num_external; ecm++) {
                      if (!usePrevious) {
                        if (out.is_open()) out << point->GetECAmplitude(k, ecm) << std::endl;
                        for (Site *subPoint = point->GetSubPoints(*this).begin();
                             subPoint < point->GetSubPoints(*this).end(); subPoint++)
                          if (out.is_open()) out << subPoint->GetECAmplitude(k, ecm) << std::endl;
                      } else {
                        complex ecAmplitude(0.0, 0.0);
                        in >> ecAmplitude;
                        point->AddECAmplitude(k, ecm, ecAmplitude, configure);
                        // Every slot of the point is in the file; the slot's writer holds it (SlotGroup.h).
                        for (int j = 1; j <= point->NumSubPoints(); j++) {
                          ecAmplitude = complex(0.0, 0.0);
                          in >> ecAmplitude;
                          if (WritesSlot(*point, j)) point->GetSubPoint(*this, j)->AddECAmplitude(k, ecm, ecAmplitude, configure);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (sumSegmentI == numSumSegments) {
      sumSegmentI = 0;
      numSumSegments = 0;
    }
  }
  if (out.is_open()) {
    out.flush();
    out.close();
  }
  if (!usePrevious && configure.useECTable && !(configure.paramMask & Config::USE_API) &&
      ECTable::Lookups(configure) + ECTable::DirectCalls(configure) > 0)
    configure.outStream << "\t" << ECTable::Report(configure) << std::endl;
  // Remember where the regular-segment integrals end so the component-segment
  // integrals (appended to the same file) can be resumed from here on read.
  if (in.is_open()) {
    ecReadPos_ = in.tellg();
    in.close();
  }
  return 0;
}

int WorkingCopy::LeafECAmplitudes(const CompoundSystem *compound, const Config &configure, bool usePrevious, std::ifstream &in,
                                  std::ofstream &out) {
  for (const std::pair<int, int> &leaf : Draft().leafRefs_) {
    int aa = compound->GetPairNumFromKey(Draft().segments_[leaf.first].TermEntranceKey(leaf.second));
    if (compound->GetPair(aa)->GetPType() != 20 && compound->GetPair(aa)->IsEntrance()) {
      const azure2::EntrancePathways &P = compound->Pathways().Entrance(aa);
      for (int j = 1; j <= compound->NumJGroups(); j++) {
        for (int la = 1; la <= compound->GetJGroup(j)->NumLevels(); la++) {
          if (compound->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
            const ALevel *ecLevel = compound->GetJGroup(j)->GetLevel(la);
            int ir = compound->GetPairNumFromKey(Draft().segments_[leaf.first].TermExitKey(leaf.second));
            if (ecLevel->GetECPairNum() == ir) {
              if (!usePrevious) {
                char segmentKeyOut[256];
                snprintf(segmentKeyOut, sizeof(segmentKeyOut), "%d", Draft().segments_[leaf.first].TermKey(leaf.second));
                if (!(configure.paramMask & Config::USE_API)) {
                  configure.outStream << "\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [                         ] 0%";
                  configure.outStream.flush();
                }
                int numPoints = (int)Draft().segments_[leaf.first].TermPoints(*this, leaf.second).size();
                std::atomic<int> pointIndex(0);
                time_t startTime = time(NULL);
                bool localStop = false;
#pragma omp parallel for shared(configure, localStop, pointIndex, startTime)
                for (int i = 1; i <= numPoints; i++) {
                  const azure2::support::ExtendedPrecisionScope fpPrecision;
                  if (configure.StopRequested() || localStop) continue;
                  Unit *point = Draft().segments_[leaf.first].GetTermPoint(*this, leaf.second, i);
                  try {
                    point->CalculateECAmplitudes(this, compound, configure, channelTable_.get());
                  } catch (GSLException e) {
#pragma omp critical(ec_progress)
                    {
                      configure.outStream << e.what() << std::endl;
                      localStop = true;
                    }
                  } catch (...) {
#pragma omp critical(ec_progress)
                    {
                      configure.outStream << "Warning: EC amplitude calculation failed for component segment point" << std::endl;
                    }
                  }
                  ++pointIndex;
                  ECProgressRedraw(configure, segmentKeyOut, pointIndex, numPoints, startTime);
                }
                if (configure.StopRequested() || localStop) {
                  if (out.is_open()) out.close();
                  if (in.is_open()) in.close();
                  return -1;
                }
                if (!(configure.paramMask & Config::USE_API))
                  configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut
                                      << std::setw(0) << " [*************************] 100%" << std::endl;
              }
              // Write the freshly calculated amplitudes to the file, or read the
              // previously saved ones back in.  The iteration order matches the
              // regular-segment loop in CalculateECAmplitudes exactly.
              const azure2::ExitRow &exit = ExternalCaptureExit(
                  *compound, aa, ir, Draft().segments_[leaf.first].TermKey(leaf.second));
              for (UnitIterator point = Draft().segments_[leaf.first].TermPoints(*this, leaf.second).begin();
                   point < Draft().segments_[leaf.first].TermPoints(*this, leaf.second).end(); point++) {
                for (int k = 1; k <= exit.num_groups; k++) {
                  for (int ecm = 1; ecm <= P.Group(exit, k).num_external; ecm++) {
                    if (!usePrevious) {
                      if (out.is_open()) out << point->GetECAmplitude(k, ecm) << std::endl;
                      for (Site *subPoint = point->GetSubPoints(*this).begin();
                           subPoint < point->GetSubPoints(*this).end(); subPoint++)
                        if (out.is_open()) out << subPoint->GetECAmplitude(k, ecm) << std::endl;
                    } else {
                      complex ecAmplitude(0.0, 0.0);
                      in >> ecAmplitude;
                      point->AddECAmplitude(k, ecm, ecAmplitude, configure);
                      for (int j = 1; j <= point->NumSubPoints(); j++) {
                        ecAmplitude = complex(0.0, 0.0);
                        in >> ecAmplitude;
                        if (WritesSlot(*point, j)) point->GetSubPoint(*this, j)->AddECAmplitude(k, ecm, ecAmplitude, configure);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return 0;
}
