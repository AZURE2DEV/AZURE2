#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"
#include <deque>

class CNuc;
struct BandData;
namespace ROOT {
  namespace Minuit2 {
    class MnUserParameters;
  }
}

///An AZURE data object

/*!
 * The EData object is the top level data object in AZURE.  It is the container object for a vector of ESegment objects.
 */

class EData {
 public:
  EData();
  int NumSegments() const;
  int Fill(const Config&,CNuc*);
  int MakePoints(const Config&,CNuc*);
  int Iterations() const;
  int NumTargetEffects() const;
  int GetNormParamOffset() const;
  int GetEnergyShiftParamOffset() const;
  int ReadTargetEffectsFile(const Config&,CNuc*);
  bool IsFit() const;
  bool IsErrorAnalysis() const;
  bool IsSegmentKey(int);
  void SetFit(bool);
  void SetErrorAnalysis(bool);
  void Iterate();
  void ResetIterations();
  int Initialize(CNuc*,const Config&);
  void AddSegment(ESegment);
  void PrintData(const Config&);
  void CalcLegendreP(int,CNuc*);
  void PrintLegendreP(const Config&);
  int CalcEDependentValues(CNuc*,const Config&);
  void PrintEDependentValues(const Config&,CNuc*);
  void CalcCoulombAmplitude(CNuc*);
  void PrintCoulombAmplitude(const Config&,CNuc*); 
  void WriteOutputFiles(const Config&,bool=false,const BandData* =nullptr);
  int CalculateECAmplitudes(CNuc*,const Config&);
  int InitializeComponentSegments(CNuc*,const Config&);
  ESegment* CreateComponentSegment(const ESegment& baseSegment, int entranceKey, int exitKey);
  ESegment* CreateComponentSegment(const ESegment& baseSegment, int entranceKey, int exitKey, double fixedAngle);
  void MapData();
  void AddTargetEffect(TargetEffect);
  void SetNormParamOffset(int);
  void SetEnergyShiftParamOffset(int);
  void FillMnParams(ROOT::Minuit2::MnUserParameters&);
  void FillNormsFromParams(const vector_r &);
  void FillEnergyShiftsFromParams(const vector_r &, EData *data = nullptr, CNuc* theCNuc = nullptr, const Config* configure = nullptr);
  void DeleteLastSegment();
  ESegment *GetSegment(int);
  ESegment *GetSegmentFromKey(int);
  EData *Clone() const;
  TargetEffect *GetTargetEffect(int);
  EDataIterator begin();
  EDataIterator end();
  std::vector<ESegment>& GetSegments();
  
 private:
  std::vector<TargetEffect> targetEffects_;
  std::vector<ESegment> segments_;
  std::deque<ESegment> componentSegments_;  // Separate storage for component segments (deque avoids pointer invalidation)
  int iterations_;
  int normParamOffset_;
  int energyShiftParamOffset_;
  bool isFit_;
  bool isErrorAnalysis_;
};

#endif
