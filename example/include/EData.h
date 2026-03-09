#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"

class CNuc;
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
  void CalcLegendreP(int);
  void PrintLegendreP(const Config&);
  int CalcEDependentValues(CNuc*,const Config&);
  void PrintEDependentValues(const Config&,CNuc*);
  void CalcCoulombAmplitude(CNuc*);
  void PrintCoulombAmplitude(const Config&,CNuc*); 
  void WriteOutputFiles(const Config&,bool=false);
  int CalculateECAmplitudes(CNuc*,const Config&);
  void MapData();
  void AddTargetEffect(TargetEffect);
  void SetNormParamOffset(int);
  void FillMnParams(ROOT::Minuit2::MnUserParameters&);
  void FillNormsFromParams(const vector_r &);
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
  int iterations_;
  int normParamOffset_;
  bool isFit_;
  bool isErrorAnalysis_;
};

#endif
