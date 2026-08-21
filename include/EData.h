#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"
#include <deque>
#include <ios>

class CNuc;
struct BandData;
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

/// An AZURE data object

/*!
 * The EData object is the top level data object in AZURE.  It is the container object for a vector of ESegment objects.
 */

class EData {
 public:
  EData();
  /// Number of data segments.
  int NumSegments() const;
  /// Build the segments from the segment and data files. Returns -1 on failure.
  int Fill(const Config &, CNuc *);
  /// Build points from <segmentsTest> grids instead of data files, for a run with no data. Returns -1 on failure.
  int MakePoints(const Config &, CNuc *);
  /// Fit iterations used so far.
  int Iterations() const;
  /// Number of target-effect definitions attached to the data.
  int NumTargetEffects() const;
  /// Index at which the normalization parameters start in the Minuit vector.
  int GetNormParamOffset() const;
  /// Index at which the energy-shift parameters start in the Minuit vector.
  int GetEnergyShiftParamOffset() const;
  /// Read the target-effects input file and build the TargetEffect objects.
  int ReadTargetEffectsFile(const Config &, CNuc *);
  /// Is this a fit? AZURECalc clones the compound nucleus and data per thread only when it is.
  bool IsFit() const;
  /// Is this a Minos error-analysis call? Suppresses the transformation and file output.
  bool IsErrorAnalysis() const;
  /// Does a segment with this key exist?
  bool IsSegmentKey(int);
  void SetFit(bool);
  void SetErrorAnalysis(bool);
  /// Count one more fit iteration.
  void Iterate();
  /// Reset the iteration counter to zero.
  void ResetIterations();
  /// Initialize every point: EPoint::Initialize over the whole data set.
  int Initialize(CNuc *, const Config &);
  /// Append a segment.
  void AddSegment(ESegment);
  /// Print the data as read, or as generated for a run without data.
  void PrintData(const Config &);
  /// EPoint::CalcLegendreP for every point.
  void CalcLegendreP(int, CNuc *);
  void PrintLegendreP(const Config &);
  /// EPoint::CalcEDependentValues for every point -- penetrabilities, shifts and Coulomb phases at each energy.
  int CalcEDependentValues(CNuc *, const Config &);
  void PrintEDependentValues(const Config &, CNuc *);
  /// EPoint::CalcCoulombAmplitude for every point.
  void CalcCoulombAmplitude(CNuc *);
  void PrintCoulombAmplitude(const Config &, CNuc *);
  /// Write AZUREOut_*, chiSquared.out and the rest of the run's output files.
  void WriteOutputFiles(const Config &, bool = false, const BandData * = nullptr);
  /// External-capture amplitudes for every point that has an EC component.
  int CalculateECAmplitudes(CNuc *, const Config &);
  /// How many external-capture amplitudes this model expects in an intEC file.
  /// Mirrors CalculateECAmplitudes exactly, so a mismatch means the cached file belongs to a different model.
  long long CountECAmplitudes(CNuc *, const Config &);
  /// Set up the entrance/exit pairs of the component segments in the compound nucleus.
  int InitializeComponentSegments(CNuc *, const Config &);
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey);
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey, double fixedAngle);
  /// Map points at equal energies onto one another, so a shared energy is calculated once.
  void MapData();
  /// Append a target-effect definition.
  void AddTargetEffect(TargetEffect);
  void SetNormParamOffset(int);
  void SetEnergyShiftParamOffset(int);
  /// Seed the Minuit parameter array with the normalizations and energy shifts.
  void FillMnParams(ROOT::Minuit2::MnUserParameters &);
  /// Write the normalizations back from a Minuit parameter vector.
  void FillNormsFromParams(const vector_r &);
  /// Write the energy shifts back from a Minuit parameter vector.
  void FillEnergyShiftsFromParams(const vector_r &, EData *data = nullptr, CNuc *theCNuc = nullptr, const Config *configure = nullptr);
  /// Drop the last segment.
  void DeleteLastSegment();
  /// Segment \p i, 1-based.
  ESegment *GetSegment(int);
  /// The segment with this key, or null if there is none.
  ESegment *GetSegmentFromKey(int);
  EData *Clone() const;
  TargetEffect *GetTargetEffect(int);
  EDataIterator begin();
  EDataIterator end();
  std::vector<ESegment> &GetSegments();

 private:
  std::vector<TargetEffect> targetEffects_;
  std::vector<ESegment> segments_;
  std::deque<ESegment> componentSegments_;  // Separate storage for component segments (deque avoids pointer invalidation)
  int iterations_;
  int normParamOffset_;
  int energyShiftParamOffset_;
  bool isFit_;
  bool isErrorAnalysis_;
  std::streampos ecReadPos_;  // File offset where component-segment EC integrals begin in the intEC file
};

#endif
