#include "CNuc.h"
#include "Config.h"
#include "DataLine.h"
#include "EData.h"
#include "ESegment.h"
#include "ExtrapLine.h"
#include "SegLine.h"

/*!
 * This constructor is used if the segment contains actual experimental data.  The segment
 * is created with reference to an entry in the segments data file.
 */

ESegment::ESegment(SegLine segLine) {
  entrancekey_=segLine.entranceKey();
  exitkey_=segLine.exitKey();
  min_e_=segLine.minE();
  max_e_=segLine.maxE();
  min_a_=segLine.minA(); 
  max_a_=segLine.maxA();
  e_step_=0.0;
  a_step_=0.0;
  segment_chi_squared_=0.0;
  if(segLine.isDiff()==1 || segLine.isDiff()==4) isdifferential_=true;
  else isdifferential_=false;
  if(segLine.isDiff()==4) iscmdifferential_=true;
  else iscmdifferential_=false;
  if(segLine.isDiff()==2) {
    isphase_=true;
    j_=segLine.phaseJ();
    l_=segLine.phaseL();
  } else {
    isphase_=false;
    j_=0.0;
    l_=0;
  }
  isTotalCapture_ =  (segLine.isDiff()==3) ? 1 : 0;
  isAngDist_=false;
  maxAngDistOrder_=0;
  datafile_=segLine.dataFile();
  dataNorm_= dataNormNominal_ = segLine.dataNorm();
  dataNormError_=segLine.dataNormError();
  energyShift_=energyShiftNominal_=segLine.energyShift();
  lastEnergyShift_=0.0;
  energyShiftError_=segLine.energyShiftError();
  if(segLine.varyEnergyShift()==1) varyEnergyShift_=true;
  else varyEnergyShift_=false;
  if(segLine.varyNorm()==1) varyNorm_=true;
  else varyNorm_=false;
  isAdvanced_=segLine.isAdvanced() == 1;
  operationType_=segLine.operationType();
  componentsList_=segLine.componentsList();
  targetEffectNum_=0;
  isTargetEffect_=false;
  
  // Initialize advanced segment composition
  segmentOperationType_ = SUM;  // Default to SUM operation
  componentCalculationMutex_ = std::make_shared<std::mutex>();
}

/*!
 * This constructor is used if the segment contains no actual data, and the points must be created
 * (such as in an extrapolation).  The segment is created with reference to an entry in the segments
 * extrapolation file.
 */

ESegment::ESegment(ExtrapLine extrapLine) {
  entrancekey_=extrapLine.entranceKey();
  exitkey_=extrapLine.exitKey();
  min_e_=extrapLine.minE();
  max_e_=extrapLine.maxE();
  min_a_=extrapLine.minA(); 
  max_a_=extrapLine.maxA();
  e_step_=extrapLine.eStep();
  a_step_=extrapLine.aStep();
  segment_chi_squared_=0.0;
  if(extrapLine.isDiff()==1 || extrapLine.isDiff()==5) isdifferential_=true;
  else isdifferential_=false;
  if(extrapLine.isDiff()==5) iscmdifferential_=true;
  else iscmdifferential_=false;
  if(extrapLine.isDiff()==2) {
    isphase_=true;
    j_=extrapLine.phaseJ();
    l_=extrapLine.phaseL();
  } else {
    isphase_=false;
    j_=0.0;
    l_=0;
  }
  if(extrapLine.isDiff()==3) {
    isAngDist_=true;
    maxAngDistOrder_=extrapLine.maxAngDistOrder();
  } else {
    isAngDist_=false;
    maxAngDistOrder_=0;
  }
  isTotalCapture_ =  (extrapLine.isDiff()==4) ? 1 : 0;
  datafile_="";
  dataNorm_= dataNormNominal_ = 1.;
  dataNormError_=0.;
  energyShift_=energyShiftNominal_=0.0;
  energyShiftError_=0.0;
  lastEnergyShift_=0.0;
  varyEnergyShift_=false;
  varyNorm_=false;
  isAdvanced_=false;
  operationType_=0;
  componentsList_="";
  targetEffectNum_=0;
  isTargetEffect_=false;
  
  // Initialize advanced segment composition
  segmentOperationType_ = SUM;  // Default to SUM operation
  componentCalculationMutex_ = std::make_shared<std::mutex>();
}

/*!
 * Returns true if a point is with the specified angle and energy ranges of a segment, otherwise
 * returns false.
 */

bool ESegment::IsInSegment(EPoint point) {
  bool b=false;
  if(point.GetLabEnergy()>=this->GetMinEnergy()&&
     point.GetLabEnergy()<=this->GetMaxEnergy()) {
    if(this->IsDifferential()) {
      if(point.GetLabAngle()>=this->GetMinAngle()&&
	 point.GetLabAngle()<=this->GetMaxAngle()) b=true;
    } else b=true;
  }
  return b;
}

/*!
 * Returns true if the segment is differential cross section, otherwise returns false.
 */

bool ESegment::IsDifferential() const {
  return isdifferential_;
}

/*!
 * Returns true if the segment is differential C.M. cross section, otherwise returns false.
 */

bool ESegment::IsCMDifferential() const {
  return iscmdifferential_;
}

/*!
 * Returns true if the segment is phase shift, otherwise returns false.
 */

bool ESegment::IsPhase() const {
  return isphase_;
}

/*!
 * Returns the number of total capture segments to be summed.  Should be zero if the segment is
 * not total capture, otherwise the parameter should be the number of segments in the sum 
 * (inclusive of the current segment).
 */

int ESegment::IsTotalCapture() const {
  return isTotalCapture_;
}

/*!
 * Returns true if the segment is angular distribution.
 */

bool ESegment::IsAngularDist() const {
  return isAngDist_;
}

/*!
 * Returns true if the segment has a corresponding TargetEffect object, otherwise
 * returns false.
 */

bool ESegment::IsTargetEffect() const {
  return isTargetEffect_;
}

/*!
 * Returns true if the normalization parameter for the segment is to be fit, otherwise
 * returns false.
 */

bool ESegment::IsVaryNorm() const {
  return varyNorm_;
}

/*!
 * Returns the number of data point objects in the segment.
 */

int ESegment::NumPoints() const {
  return points_.size();
}

/*!
 * Returns the entrance particle pair key of the segment.
 */

int ESegment::GetEntranceKey() const {
  return entrancekey_;
}

/*!
 * Returns the exit particle pair key of the segment.
 */

int ESegment::GetExitKey() const {
  return exitkey_;
}

/*!
 * Fills the segment with points from the specified data file according to the
 * maximum and minimum energy and angle ranges.  If the data 
 * in the lab frame, conversions are performed to the center of mass frame.
 */

int ESegment::Fill(CNuc *theCNuc, EData *theData, const Config& configure) {
  std::string infile=this->GetDataFile();
  std::ifstream in(infile.c_str());
  if(!in) return -1;
  while(!in.eof()) {
    DataLine line(in);
    if(!in.eof()) {
      EPoint NewEPoint(line,this);
      if(this->IsInSegment(NewEPoint)) {
	this->AddPoint(NewEPoint);
	PPair *entrancePair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
	PPair *exitPair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
	this->GetPoint(this->NumPoints())->SetParentData(theData);
	if(entrancePair->GetPType()==20) {
           this->GetPoint(this->NumPoints())->ConvertDecayEnergy(exitPair);
        } else if (this->IsCMDifferential()) {
//          this->GetPoint(this->NumPoints())->ConvertExcitationEnergy(entrancePair);
          this->GetPoint(this->NumPoints())->ConvertLabEnergy(entrancePair); 
        } else if(!this->IsCMDifferential()){
          this->GetPoint(this->NumPoints())->ConvertLabEnergy(entrancePair);
        }
	if(exitPair->GetPType()==0&&this->IsDifferential()&&!this->IsPhase()&&!this->IsCMDifferential()) {
	  if(this->GetEntranceKey()==this->GetExitKey()) {
	    this->GetPoint(this->NumPoints())->ConvertLabAngle(entrancePair);
	  } else {
	    this->GetPoint(this->NumPoints())->ConvertLabAngle(entrancePair,exitPair,configure);
	  }
	  this->GetPoint(this->NumPoints())->ConvertCrossSection(entrancePair,exitPair);
	}
//        if(exitPair->GetPType()==0&&!this->IsDifferential()&&!this->IsPhase()&&this->IsCMDifferential()) {
//          this->GetPoint(this->NumPoints())->ConvertCMAngle(entrancePair,exitPair,configure);
//        }
	if(exitPair->GetPType()==10&&this->IsDifferential()&&!this->IsPhase()&&!this->IsCMDifferential()) {
	  this->GetPoint(this->NumPoints())->ConvertLabAngleGammas(entrancePair);
	  this->GetPoint(this->NumPoints())->ConvertCrossSectionGammas(entrancePair);
	}
      }
    }
  }
  in.close();
  return 0;
}

/*!
 * Returns the orbital angular momentum value for the segment.  Applies only if the segment is phase shift.
 */

int ESegment::GetL() const {
  return l_;
}

/*!
 * Returns the position of the corresponding TargetEffect object in the vector
 *  of the parent EData object.
 */

int ESegment::GetTargetEffectNum() const {
  return targetEffectNum_;
}

/*! 
 * Returns the segment key for the current ESegment object.  The segment key is the order of 
 * the segment specified in the input file, INCLUDING non-active segments.
 */

int ESegment::GetSegmentKey() const {
  return segmentKey_;
}

/*!
 * Returns the maximum polynomial order if segment is angular distribution.
 */

int ESegment::GetMaxAngDistOrder() const {
  return maxAngDistOrder_;
}

/*!
 * Returns the minimum energy of the segment (lab frame).
 */

double ESegment::GetMinEnergy() const {
  return min_e_;
}

/*!
 * Returns the maximum energy of the segment (lab frame).
 */

double ESegment::GetMaxEnergy() const {
  return max_e_;
}

/*!
 * Returns the minimum angle of the segment (lab frame).
 */

double ESegment::GetMinAngle() const {
  return min_a_;
}

/*!
 * Returns the maximum angle of the segment (lab frame).
 */

double ESegment::GetMaxAngle() const {
  return max_a_;
}

/*!
 * Returns the chi-squared value of the segment after the fitting process.
 */

double ESegment::GetSegmentChiSquared() const {
  return segment_chi_squared_;
}

/*!
 * Returns the chi-squared value due to the normalization for each segment after the fitting process.
 */

//double ESegment::GetNormChiSquared() const {
//  return norm_chi_squared_;
//}

/*!
 * Returns the energy step to take when creating points in a segment.  Only applies for extrapolation segments.
 */

double ESegment::GetEStep() const {
  return e_step_;
}

/*!
 * Returns the angle step to take when creating points in a segment.  Only applies for extrapolation segments.
 */

double ESegment::GetAStep() const {
  return a_step_;
}

/*!
 * Returns the total spin of the segment.  Only applies if the segment is phase shift.
 */

double ESegment::GetJ() const {
  return j_;
}

/*!
 * Returns the normalization parameter for the data segment.
 */

double ESegment::GetNorm() const {
  return dataNorm_;
}

/*!
 * Returns the nominal normalization parameter for the data segment.
 */

double ESegment::GetNominalNorm() const {
  return dataNormNominal_;
}

/*!
 * Returns the normalization error for the data segment.
 */

double ESegment::GetNormError() const {
  return dataNormError_;
}

/*!
 * Returns the name of the data file from which to read.
 */

std::string ESegment::GetDataFile() const {
  return datafile_;
}

/*!
 * Returns the energy shift value for the data segment.
 */

double ESegment::GetEnergyShift() const {
  return energyShift_;
}

/*!
 * Returns the last applied energy shift value for the data segment.
 */
double ESegment::GetLastEnergyShift() const {
  return lastEnergyShift_;
}

/*!
 * Returns the nominal energy shift value for the data segment.
 */

double ESegment::GetNominalEnergyShift() const {
  return energyShiftNominal_;
}

/*!
 * Returns the energy shift error for the data segment.
 */

double ESegment::GetEnergyShiftError() const {
  return energyShiftError_;
}

/*!
 * Returns true if the energy shift is varied in the fit.
 */

bool ESegment::IsVaryEnergyShift() const {
  return varyEnergyShift_;
}

/*!
 * Returns true if this is an advanced segment (sum/ratio).
 */

bool ESegment::IsAdvanced() const {
  return isAdvanced_;
}

/*!
 * Returns the operation type (0 for sum, 1 for ratio).
 */

int ESegment::GetLegacyOperationType() const {
  return operationType_;
}

/*!
 * Returns the semicolon-separated list of components.
 */

std::string ESegment::GetComponentsList() const {
  return componentsList_;
}

/*!
 * Sets the energy shift value for the data segment.
 */

void ESegment::SetEnergyShift(double energyShift) {
  energyShift_=energyShift;
}

/*!
 * Sets the last applied energy shift value for the data segment.
 */
void ESegment::SetLastEnergyShift(double lastEnergyShift) {
  lastEnergyShift_=lastEnergyShift;
}

/*!
 * Updates the energies of all points in this segment based on the current energy shift.
 */

void ESegment::UpdatePointEnergiesWithShift(CNuc* theCNuc, const Config* configure) {
  // No lock needed since each thread works on independent EData clones
  // after fixing EData::Clone to properly remap component segment pointers
  
  for(int i = 0; i < NumPoints(); i++) {
    EPoint* point = GetPoint(i+1);
    if(point && point->GetOriginalEnergy() > 0) {
      double originalEnergy = point->GetOriginalEnergy();
      double shiftedEnergy = originalEnergy + energyShift_;
      
      // Don't allow energies below 0.01 MeV (AZURE2 may crash)
      if(shiftedEnergy < 0.01) {
        // FIX: we do not want to change the energy so abruptly, so leave the last shifted energy
        //shiftedEnergy = originalEnergy;
        continue;
      }
      
      // Set the shifted energy
      point->SetLabEnergy(shiftedEnergy);

      // Recalculate energy dependent values
      if(theCNuc && configure) {
        PPair* entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
        PPair* exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));

        // Apply the same conversion functions as in ESegment::Fill
        if(entrancePair->GetPType()==20) {
          point->ConvertDecayEnergy(exitPair);
        } else if (this->IsCMDifferential()) {
          point->ConvertLabEnergy(entrancePair); 
        } else if(!this->IsCMDifferential()){
          point->ConvertLabEnergy(entrancePair);
        }

        if(exitPair->GetPType()==0 && this->IsDifferential() && !this->IsPhase() && !this->IsCMDifferential()) {
          if(this->GetEntranceKey()==this->GetExitKey()) {
            point->ConvertLabAngle(entrancePair);
          } else {
            point->ConvertLabAngle(entrancePair,exitPair,*configure);
          }
          point->ConvertCrossSection(entrancePair,exitPair);
        }

        if(exitPair->GetPType()==10 && this->IsDifferential() && !this->IsPhase() && !this->IsCMDifferential()) {
          point->ConvertLabAngleGammas(entrancePair);
          point->ConvertCrossSectionGammas(entrancePair);
        }

        if( this->IsDifferential() || this->IsCMDifferential()) {
          point->ClearLegendrePolynomials();
          point->CalcLegendreP(configure->maxLOrder, theCNuc, NULL);
        }
        point->RecalcEDependentValues(theCNuc,*configure);

      }
    }
  }
}

/*!
 * Adds a data point to the segment.
 */

void ESegment::AddPoint(EPoint point) {
  points_.push_back(point);
}

/*!
 * Sets the chi squared value for the segment during the fitting process.  
 */

void ESegment::SetSegmentChiSquared(double chiSquared) {
  segment_chi_squared_=chiSquared;
}

/*!
 * Sets the chi squared value for the segment normalization during the fitting process.  
 */

//void ESegment::SetNormChiSquared(double chiSquared) {
//  norm_chi_squared_=chiSquared;
//}

/*!
 * Sets the position of the corresponding TargetEffect object in the vector
 * of the parent EData object.
 */

void ESegment::SetTargetEffectNum(int targetEffectNum) {
  targetEffectNum_=targetEffectNum;
  isTargetEffect_=true;
}

/*! 
 * Sets the segment key for the current ESegment object.
 */

void ESegment::SetSegmentKey (int segmentKey) {
  segmentKey_=segmentKey;
}

/*! 
 * Sets the normalization parameter for the segment.
 */ 

void ESegment::SetNorm(double norm) {
  dataNorm_=norm;
}

/*!
 * Sets the exit pair key to the given value.
 */

void ESegment::SetExitKey(int key) {
  exitkey_=key;
  for(int i=1;i<=this->NumPoints();i++) 
    this->GetPoint(i)->SetExitKey(key);
}

/*!
 * Sets the entrance key for the segment and all its points
 */
void ESegment::SetEntranceKey(int key) {
  entrancekey_=key;
  for(int i=1;i<=this->NumPoints();i++) 
    this->GetPoint(i)->SetEntranceKey(key);
}

/*!
 * Sets the number of total capture segments to be summed.  Should be zero if the segment is
 * not total capture, otherwise the parameter should be the number of segments in the sum 
 * (inclusive of the current segment).
 */

void ESegment::SetIsTotalCapture(int num) {
  isTotalCapture_=num;
}

/*!
 * Set the flag determining if the normalization is varied.
 */

void ESegment::SetVaryNorm(bool varyNorm) {
  varyNorm_=varyNorm;

}
/*!
 * Returns a pointer to the data point object specified by a position in the EPoint vector.
 */

EPoint *ESegment::GetPoint(int pointNum) {
  EPoint *b=&points_[pointNum-1];
  return b;
}

/*!
 * Returns a reference to the vector of EPoint objects.
 */

std::vector<EPoint>& ESegment::GetPoints() {
  return points_;
}

// Advanced segment composition methods

/*!
 * Adds a component segment defined by entrance and exit keys
 */
void ESegment::AddComponent(int entranceKey, int exitKey) {
  // This method is now deprecated - use AddComponentSegment instead
  // Components should be full ESegment copies, not just entrance/exit keys
}

/*!
 * Adds a component segment (full ESegment copy)
 */
void ESegment::AddComponentSegment(ESegment* componentSegment) {
  if (componentSegment) {
    componentSegments_.push_back(componentSegment);
  }
}

/*!
 * Sets the operation type for combining components (SUM or RATIO)
 */
void ESegment::SetOperationType(OperationType operation) {
  segmentOperationType_ = operation;
}

/*!
 * Gets the operation type for combining components
 */
OperationType ESegment::GetOperationType() const {
  return segmentOperationType_;
}

/*!
 * Returns true if this segment has additional components
 */
bool ESegment::HasComponents() const {
  return !componentSegments_.empty();
}

/*!
 * Returns a reference to the component segments vector
 */
const std::vector<ESegment*>& ESegment::GetComponentSegments() const {
  return componentSegments_;
}

/*!
 * Clears all components from the segment
 */
void ESegment::ClearComponents() {
  // Clean up component segments (they are owned by EData, so just clear pointers)
  componentSegments_.clear();
}

/*!
 * Calculates the theoretical cross section for this segment, including 
 * combination with any additional components according to the operation type
 */
double ESegment::CalculateTheoreticalCrossSection(int pointIndex, CNuc* cnuc, const Config& configure, EData* edata) {
  if (pointIndex < 0 || pointIndex >= NumPoints()) {
    return 0.0;
  }
  
  // Calculate the base segment's theoretical cross section
  EPoint* basePoint = GetPoint(pointIndex + 1);
  if (!basePoint) {
    return 0.0;
  }
  
  // Calculate base point if not already calculated
  if (!basePoint->IsMapped()) {
    try {
      basePoint->Calculate(cnuc, configure);
    } catch (...) {
      return 0.0;
    }
  }
  
  double baseValue = basePoint->GetFitCrossSection();
  
  // If no components, return base value
  if (!HasComponents()) {
    return baseValue;
  }
  
  
  // Calculate component values and combine according to operation type
  if (segmentOperationType_ == SUM) {
    double sum = baseValue;
    for (const auto& componentSegment : componentSegments_) {
      if (componentSegment && pointIndex < componentSegment->NumPoints()) {
        EPoint* componentPoint = componentSegment->GetPoint(pointIndex + 1);
        if (componentPoint) {
          // No lock needed since each thread works on independent EData clones
          if (!componentPoint->IsMapped()) {
            try {
              componentPoint->Calculate(cnuc, configure);
            } catch (...) {
              continue; // Skip this component if calculation fails
            }
          }
          sum += componentPoint->GetFitCrossSection();
        }
      }
    }
    return sum;
  } else if (segmentOperationType_ == RATIO) {
    if (componentSegments_.empty()) {
      return baseValue;
    }
    
    // For ratio: base / first_component
    ESegment* firstComponent = componentSegments_[0];
    if (firstComponent && pointIndex < firstComponent->NumPoints()) {
      EPoint* denominatorPoint = firstComponent->GetPoint(pointIndex + 1);
      if (denominatorPoint) {
        
        // No lock needed since each thread works on independent EData clones
        if (!denominatorPoint->IsMapped()) {
          try {
            denominatorPoint->Calculate(cnuc, configure);
          } catch (...) {
            return 0.0;
          }
        }
        
        double denominatorValue = denominatorPoint->GetFitCrossSection();
        if (denominatorValue == 0.0) {
          return 0.0;
        }
        
        return baseValue / denominatorValue;
      }
    }
    return 0.0;
  }
  
  return baseValue;
}
