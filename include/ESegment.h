#ifndef ESEGMENT_H
#define ESEGMENT_H

#include "EPoint.h"
#include <vector>
#include <mutex>
#include <memory>

class EData;
class ExtrapLine;
class SegLine;
class CNuc;
class Config;

enum OperationType {
    SUM = 0,
    RATIO = 1
};

///An AZURE data segment

/*!
 * An AZURE data segment is specified by an entrance and exit particle pair key, as well as a range of 
 * energy and angle value and a data file name. The segment also contains flags specifing they type of data point it contains. The ESegment object is the container object for a vector of EData objects.
 */

class ESegment {
 public:
  ESegment(SegLine); 
  ESegment(ExtrapLine); 
  bool IsInSegment(EPoint);
  bool IsDifferential() const;
  bool IsCMDifferential() const;
  bool IsPhase() const;
  bool IsTargetEffect() const;
  bool IsVaryNorm() const;
  bool IsAngularDist() const;
  bool IsUPOS() const;
  int IsTotalCapture() const;
  int NumPoints() const;
  int GetEntranceKey() const;
  int GetExitKey() const;
  int Fill(CNuc*,EData*,const Config&);
  int GetL() const;
  int GetTargetEffectNum() const;
  int GetSegmentKey() const;
  int GetMaxAngDistOrder() const;
  int GetSecondaryDecayL() const;
  double GetMinEnergy() const;
  double GetMaxEnergy() const;
  double GetMinAngle() const;
  double GetMaxAngle() const;
  double GetSegmentChiSquared() const;
  double GetEStep() const;
  double GetAStep() const;
  double GetJ() const;
  double GetNorm() const;
  double GetNominalNorm() const;
  double GetNormError() const;
  double GetIc() const;
  double GetDelta() const;
  double GetEnergyShift() const;
  double GetLastEnergyShift() const;
  double GetNominalEnergyShift() const;
  double GetEnergyShiftError() const;
  bool IsVaryEnergyShift() const;
  bool IsAdvanced() const;
  int GetLegacyOperationType() const;
  std::string GetComponentsList() const;
  std::string GetDataFile() const;
  void AddPoint(EPoint);
  void SetSegmentChiSquared(double);
  void SetTargetEffectNum(int);
  void SetSegmentKey(int);
  void SetNorm(double);
  void SetEnergyShift(double);
  void SetLastEnergyShift(double);
  void UpdatePointEnergiesWithShift(CNuc* theCNuc = NULL, const Config* configure = NULL);
  void SetExitKey(int);
  void SetEntranceKey(int);
  void SetIsTotalCapture(int);
  void SetVaryNorm(bool);
  void SetMinAngle(double);
  void SetMaxAngle(double);
  EPoint *GetPoint(int);
  std::vector<EPoint>& GetPoints();
  
  // Advanced segment composition methods
  void AddComponent(int entranceKey, int exitKey);
  void AddComponentSegment(ESegment* componentSegment);
  void SetOperationType(OperationType operation);
  OperationType GetOperationType() const;
  bool HasComponents() const;
  const std::vector<ESegment*>& GetComponentSegments() const;
  void ClearComponents();
  
  // Calculate theoretical cross section including components
  double CalculateTheoreticalCrossSection(int pointIndex, CNuc* cnuc, const Config& configure, EData* edata);
 private:
  bool isdifferential_;
  bool iscmdifferential_;
  bool isphase_;
  bool isTargetEffect_;
  bool varyNorm_;
  bool isAngDist_;
  bool isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
  int isTotalCapture_;
  int entrancekey_;
  int exitkey_;
  int l_;
  int targetEffectNum_;
  int segmentKey_;
  int maxAngDistOrder_;
  double min_e_;
  double max_e_;
  double min_a_;
  double max_a_;
  double e_step_;
  double a_step_;
  double segment_chi_squared_;
  double j_;
  double dataNorm_;
  double dataNormNominal_;
  double dataNormError_;
  double energyShift_;
  double lastEnergyShift_; // To track last applied energy shift
  double energyShiftNominal_;
  double energyShiftError_;
  bool varyEnergyShift_;
  bool isAdvanced_;
  int operationType_;
  std::string componentsList_;
  std::string datafile_;
  std::vector<EPoint> points_;
  
  // Advanced segment composition
  std::vector<ESegment*> componentSegments_;
  OperationType segmentOperationType_;
  mutable std::shared_ptr<std::mutex> componentCalculationMutex_;
};

#endif
