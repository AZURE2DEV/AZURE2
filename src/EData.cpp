#include "AZUREOutput.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ECAmplitudeCache.h"
#include "ExtrapLine.h"
#include "SegLine.h"
#include "AdaptiveIntegrationGrid.h"
#include "Minuit2/MnUserParameters.h"
#include "GSLException.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <time.h>
#include <unordered_map>

/*!
 * The EData object has a private attribute containing the number of iterations needed to find the best fit parameters.
 * At creation, this attribute is set to 0.
 */

EData::EData() {
  iterations_=0;
  normParamOffset_=0;
  energyShiftParamOffset_=0;
  isFit_=true;
  isErrorAnalysis_=false;
}

/*!
 * Returns the number of segment objects in the ESegment vector.
 */

int EData::NumSegments() const {
  return segments_.size();
}

/*!
 * This function fills the data object with segments from the segment data files.
 * After a segment is created, the ESegment::Fill method is called for that segment.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int EData::Fill(const Config& configure, CNuc *theCNuc) {
  std::ifstream in(configure.configfile.c_str());
  if(!in) return -1;
  std::string line="";
  while(line!="<segmentsData>"&&!in.eof()) getline(in,line);
  if(line!="<segmentsData>") return -1;
  line="";
  int numTotalSegments=0;
  while(!in.eof()&&line!="</segmentsData>") {
    getline(in,line);
    bool empty=true;
    for(unsigned int i=0;i<line.size();++i) 
      if(line[i]!=' '&&line[i]!='\t') {
	empty=false;
	break;
      }
    if(empty==true) continue;
    if(!in.eof()&&line!="</segmentsData>") {
      std::istringstream stm;
      stm.str(line);
      SegLine segment(stm);
      if(stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      numTotalSegments++;
      if(segment.isActive()==1) {
	ESegment NewSegment(segment);
	
	if(theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
	  theCNuc->GetPair(theCNuc->GetPairNumFromKey(NewSegment.GetEntranceKey()))->SetEntrance();
	  bool isValidTotal=false;
	  if(NewSegment.GetExitKey()==-1) {
	    for(int i = 1;i<=theCNuc->NumPairs();i++) {
	      if(theCNuc->GetPair(i)->GetPType()==10) {
		isValidTotal = true;
		break;
	      }
	    }
	  }
	  if(isValidTotal||theCNuc->IsPairKey(NewSegment.GetExitKey())) {
	    NewSegment.SetSegmentKey(numTotalSegments);
	    this->AddSegment(NewSegment);
	    
	    // Fill the segment with data first
	    if(this->GetSegment(this->NumSegments())->Fill(theCNuc,this,configure)==-1) {
              configure.outStream << "WARNING: Could Not Fill Segment #" << this->NumSegments() 
			<< " from file." << std::endl;
	      this->DeleteLastSegment();
	    } else if(this->GetSegment(this->NumSegments())->NumPoints()==0) {
	      configure.outStream << "WARNING: Segment #" << numTotalSegments
				  << " is empty and will not be used." << std::endl;
	      this->DeleteLastSegment();
	    } else {
	      // Handle advanced segments - parse components and add them to the segment
	      // This must be done AFTER the segment is filled with data points
	      if (NewSegment.IsAdvanced()) {
	        // Set the operation type from the segment data
	        OperationType operation = (NewSegment.GetLegacyOperationType() == 0) ? SUM : RATIO;
	        ESegment* filledSegment = this->GetSegment(this->NumSegments());
	        filledSegment->SetOperationType(operation);
	        
	        // Parse components and create full segment copies
	        std::string componentsStr = NewSegment.GetComponentsList();
	        if (!componentsStr.empty()) {
	          // Components are stored as "Entrance: X, Exit: Y;Entrance: A, Exit: B;..." or with optional "Angle: Z"
	          std::istringstream stream(componentsStr);
	          std::string component;
	          while (std::getline(stream, component, ';')) {
	            if (!component.empty()) {
	              // Parse "Entrance: X, Exit: Y" format (with optional ", Angle: Z")
	              size_t entrancePos = component.find("Entrance: ");
	              size_t exitPos = component.find("Exit: ");
	              size_t anglePos = component.find("Angle: ");

	              if (entrancePos != std::string::npos && exitPos != std::string::npos) {
	                entrancePos += 10; // Length of "Entrance: "
	                size_t commaPos = component.find(", Exit: ");
	                if (commaPos != std::string::npos) {
	                  int entranceKey = std::stoi(component.substr(entrancePos, commaPos - entrancePos));
	                  exitPos += 6; // Length of "Exit: "

	                  // Check if there's an angle specification
	                  ESegment* componentSegment = nullptr;
	                  if (anglePos != std::string::npos) {
	                    // Parse the exit key up to ", Angle:"
	                    size_t angleCommaPos = component.find(", Angle: ");
	                    int exitKey = std::stoi(component.substr(exitPos, angleCommaPos - exitPos));

	                    // Parse the angle value
	                    anglePos += 7; // Length of "Angle: "
	                    double fixedAngle = std::stod(component.substr(anglePos));

	                    // Create component segment with fixed angle
	                    componentSegment = this->CreateComponentSegment(*filledSegment, entranceKey, exitKey, fixedAngle);
	                  } else {
	                    // No angle specified, use standard method
	                    int exitKey = std::stoi(component.substr(exitPos));
	                    componentSegment = this->CreateComponentSegment(*filledSegment, entranceKey, exitKey);
	                  }

	                  if (componentSegment) {
	                    filledSegment->AddComponentSegment(componentSegment);
	                  }
	                }
	              }
	            }
	          }
	        }
	      }
	      
	      int thisSegmentNum = this->NumSegments();
	      if(this->GetSegment(thisSegmentNum)->IsTotalCapture()) {
	        ESegment* baseSegment = this->GetSegment(thisSegmentNum);
	        int numCapturePairs=0;
	        
	        // Set operation type to SUM for total capture
	        baseSegment->SetOperationType(SUM);
	        
	        for(int i = 1;i<=theCNuc->NumPairs();i++) {
	          if(theCNuc->GetPair(i)->GetPType()==10) {
	            numCapturePairs++;
	            if(numCapturePairs==1) {
	              // For the first capture pair, set the base segment's exit key
	              baseSegment->SetExitKey(theCNuc->GetPair(i)->GetPairKey());
	            } else {
	              // For additional capture pairs, create component segments instead of separate segments
	              ESegment* componentSegment = this->CreateComponentSegment(*baseSegment, 
	                                                                       baseSegment->GetEntranceKey(), 
	                                                                       theCNuc->GetPair(i)->GetPairKey());
	              if (componentSegment) {
	                baseSegment->AddComponentSegment(componentSegment);
	              }
	            }
	          }
	        }
	        baseSegment->SetIsTotalCapture(numCapturePairs);
	      } 
	    }
	  } else {
	    if(NewSegment.GetExitKey()==-1) {
	      configure.outStream << "WARNING: Total capture specified but no capture pair exists." 
				  << std::endl;
	    } else {
	      configure.outStream << "WARNING: Pair key " << NewSegment.GetExitKey() 
				  << " not in compound nucleus." << std::endl;
	    }
	  }
	} else configure.outStream << "WARNING: Pair key " << NewSegment.GetEntranceKey() 
			 << " not in compound nucleus." << std::endl;
      }
    }
  }

  if(line!="</segmentsData>") return -1;

  in.close();

  if(this->NumSegments()>0) {
    if(this->ReadTargetEffectsFile(configure,theCNuc)==-1) return -1;
    this->MapData();
  }

  return 0; 
}

/*!
 * If the AZURE calculation is not data driven, this function is called in place
 * of the EData::Fill function to create points at specified energies and angles.
 * Returns -1 if the input files could not be read, otherwise returns 0.
 */

int EData::MakePoints(const Config& configure, CNuc *theCNuc) {
  std::ifstream in(configure.configfile.c_str());
  if(!in) return -1;
  std::string line = "";
  while(line!="<segmentsTest>"&&!in.eof()) getline(in,line);
  if(line!="<segmentsTest>") return -1;
  line="";
  int numTotalSegments=0;
  while(!in.eof()&&line!="</segmentsTest>") {
    getline(in,line);
    bool empty=true;
    for(unsigned int i=0;i<line.size();++i) 
      if(line[i]!=' '&&line[i]!='\t') {
	empty=false;
	break;
      }
    if(empty==true) continue;
    if(!in.eof()&&line!="</segmentsTest>") {
      std::istringstream stm;
      stm.str(line);
      ExtrapLine segment(stm);
      if(stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      numTotalSegments++;
      if(segment.isActive()==1) {
	ESegment NewSegment(segment);
	if(theCNuc->IsPairKey(NewSegment.GetEntranceKey())) {
	  bool isValidTotal=false;
	  if(NewSegment.GetExitKey()==-1) {
	    for(int i = 1;i<=theCNuc->NumPairs();i++) {
	      if(theCNuc->GetPair(i)->GetPType()==10) {
		isValidTotal = true;
		break;
	      }
	    }
	  }
	  if(isValidTotal||theCNuc->IsPairKey(NewSegment.GetExitKey())) {
	    NewSegment.SetSegmentKey(numTotalSegments);
	    this->AddSegment(NewSegment);
	    ESegment *theSegment=this->GetSegment(this->NumSegments());

	    theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()))->SetEntrance();
	    PPair *entrancePair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetEntranceKey()));
	    PPair *exitPair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(theSegment->GetExitKey()));
	    double aStep=theSegment->GetAStep();
	    double eStep=theSegment->GetEStep();
	    for(double angle=theSegment->GetMinAngle();
		angle<=theSegment->GetMaxAngle();angle+=aStep) {
	      for(double energy=theSegment->GetMinEnergy();
		  energy<=theSegment->GetMaxEnergy();energy+=eStep) {
		EPoint NewPoint(angle,energy,theSegment);
		theSegment->AddPoint(NewPoint);
		EPoint *thePoint=theSegment->GetPoint(theSegment->NumPoints());
		thePoint->SetParentData(this);
		if(entrancePair->GetPType()==20) thePoint->ConvertDecayEnergy(exitPair);
		else if(!theSegment->IsCMDifferential()) thePoint->ConvertLabEnergy(entrancePair);
    else if(theSegment->IsCMDifferential()) thePoint->ConvertLabEnergy(entrancePair);
		if(exitPair->GetPType()==0&&theSegment->IsDifferential()&&
		   !theSegment->IsPhase()&&!theSegment->IsAngularDist()&&!theSegment->IsCMDifferential()) {
		  if(theSegment->GetEntranceKey()==theSegment->GetExitKey()) {
		    thePoint->ConvertLabAngle(entrancePair);
		  } else {
		    thePoint->ConvertLabAngle(entrancePair,exitPair,configure);
		  }
		  thePoint->ConvertCrossSection(entrancePair,exitPair);
		}
		if(exitPair->GetPType()==10&&theSegment->IsDifferential()&&
		  !theSegment->IsPhase()&&!theSegment->IsAngularDist()&&!theSegment->IsCMDifferential()){
		    thePoint->ConvertLabAngleGammas(entrancePair);
		    thePoint->ConvertCrossSectionGammas(entrancePair);
		}
		if(eStep==0.0) break;
	      }
	      if(aStep==0.0) break;
	    }

	    // Handle advanced segments (sum/ratio of components) - AFTER points are added
	    if(NewSegment.IsAdvanced()) {
	      int operation = NewSegment.GetOperationType();
	      theSegment->SetOperationType((OperationType)operation);

	      // Parse components and create full segment copies
	      std::string componentsStr = NewSegment.GetComponentsList();
	      if (!componentsStr.empty()) {
	        // Components are stored as "Entrance: X, Exit: Y;Entrance: A, Exit: B;..." or with optional "Angle: Z"
	        std::istringstream stream(componentsStr);
	        std::string component;
	        int componentCount = 0;
	        while (std::getline(stream, component, ';')) {
	          if (!component.empty()) {
	            // Parse "Entrance: X, Exit: Y" format (with optional ", Angle: Z")
	            size_t entrancePos = component.find("Entrance: ");
	            size_t exitPos = component.find("Exit: ");
	            size_t anglePos = component.find("Angle: ");

	            if (entrancePos != std::string::npos && exitPos != std::string::npos) {
	              entrancePos += 10; // Length of "Entrance: "
	              size_t commaPos = component.find(", Exit: ");
	              if (commaPos != std::string::npos) {
	                int entranceKey = std::stoi(component.substr(entrancePos, commaPos - entrancePos));
	                exitPos += 6; // Length of "Exit: "

	                // Check if there's an angle specification
	                ESegment* componentSegment = nullptr;
	                if (anglePos != std::string::npos) {
	                  // Parse the exit key up to ", Angle:"
	                  size_t angleCommaPos = component.find(", Angle: ");
	                  int exitKey = std::stoi(component.substr(exitPos, angleCommaPos - exitPos));

	                  // Parse the angle value
	                  anglePos += 7; // Length of "Angle: "
	                  double fixedAngle = std::stod(component.substr(anglePos));

	                  // Create component segment with fixed angle
	                  componentSegment = this->CreateComponentSegment(*theSegment, entranceKey, exitKey, fixedAngle);
	                } else {
	                  // No angle specified, use standard method
	                  int exitKey = std::stoi(component.substr(exitPos));
	                  componentSegment = this->CreateComponentSegment(*theSegment, entranceKey, exitKey);
	                }

	                if (componentSegment) {
	                  theSegment->AddComponentSegment(componentSegment);
	                  componentCount++;
	                }
	              }
	            }
	          }
	        }
	      } else {
	      }
	    }
	    if(theSegment->NumPoints()==0) {
	      configure.outStream << "WARNING: Extrapolation segment #" << numTotalSegments 
				  << " is empty and will not be used." << std::endl;
	      this->DeleteLastSegment();
	    } else {
	      int thisSegmentNum = this->NumSegments();
	      if(this->GetSegment(thisSegmentNum)->IsTotalCapture()) {
	        ESegment* baseSegment = this->GetSegment(thisSegmentNum);
	        int numCapturePairs=0;
	        
	        // Set operation type to SUM for total capture
	        baseSegment->SetOperationType(SUM);
	        
	        for(int i = 1;i<=theCNuc->NumPairs();i++) {
	          if(theCNuc->GetPair(i)->GetPType()==10) {
	            numCapturePairs++;
	            if(numCapturePairs==1) {
	              // For the first capture pair, set the base segment's exit key
	              baseSegment->SetExitKey(theCNuc->GetPair(i)->GetPairKey());
	            } else {
	              // For additional capture pairs, create component segments instead of separate segments
	              ESegment* componentSegment = this->CreateComponentSegment(*baseSegment, 
	                                                                       baseSegment->GetEntranceKey(), 
	                                                                       theCNuc->GetPair(i)->GetPairKey());
	              if (componentSegment) {
	                baseSegment->AddComponentSegment(componentSegment);
	              }
	            }
	          }
	        }
	        baseSegment->SetIsTotalCapture(numCapturePairs);
	      } 
	    }
	  } else {
	    if(NewSegment.GetExitKey()==-1) {
	      configure.outStream << "WARNING: Total capture specified but no capture pair exists." 
				  << std::endl;
	    } else {
	      configure.outStream << "WARNING: Pair key " << NewSegment.GetExitKey() 
				  << " not in compound nucleus." << std::endl;
	    }
	  }
	} else configure.outStream << "WARNING: Pair key " << NewSegment.GetEntranceKey() 
				   << " not in compound nucleus." << std::endl;
      }
    }
  }  

  if(line!="</segmentsTest>") return -1;
  
  in.close();

  if(this->NumSegments()>0) {
    if(this->ReadTargetEffectsFile(configure,theCNuc)==-1) return -1;
    this->MapData();
  }

  return 0; 
}

/*!
 * Returns the number of fit iterations needed to minimize the parameters to the data.
 */

int EData::Iterations() const {
  return iterations_;
}

/*!
 * Returns the number of TargetEffect objects contained in the present object.
 */

int EData::NumTargetEffects() const {
  return targetEffects_.size();
}

/*!
 * Returns the offset of the normalization paramters in the Minuit parameter vector.
 */
int EData::GetNormParamOffset() const {
  return normParamOffset_;
}

/*!
 * Returns the offset of the energy shift parameters in the Minuit parameter vector.
 */
int EData::GetEnergyShiftParamOffset() const {
  return energyShiftParamOffset_;
}

/*!
 * Reads the target effects input file and creates the TargetEffect objects 
 * to be applied to the data.
 */

int EData::ReadTargetEffectsFile(const Config& configure, CNuc *compound) {
  std::ifstream in(configure.configfile.c_str());
  if(!in) return -1;
  std::string line="";
  while(line!="<targetInt>"&&!in.eof()) getline(in,line);
  if(line!="<targetInt>") return -1;
  line="";
  while(line!="</targetInt>"&&!in.eof()) {
    getline(in,line);
    bool empty=true;
    for(unsigned int i=0;i<line.size();++i) 
      if(line[i]!=' '&&line[i]!='\t') {
	empty=false;
	break;
      }
    if(empty==true) continue;
    if(line!="</targetInt>"&&!in.eof()){
      std::istringstream stm;
      stm.str(line);
      TargetEffect targetEffect(stm,configure);
      if(stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)) return -1;
      if(targetEffect.IsActive()) {
	this->AddTargetEffect(targetEffect);
	TargetEffect *thisTargetEffect=this->GetTargetEffect(this->NumTargetEffects());
	std::vector<int> segmentsList = thisTargetEffect->GetSegmentsList();
	for(int i = 1;i<=segmentsList.size();i++) { 
	  if(this->IsSegmentKey(segmentsList[i-1])) {
	    ESegment* segment = this->GetSegmentFromKey(segmentsList[i-1]);
	    if(segment) {
	      segment->SetTargetEffectNum(this->NumTargetEffects());
	      // If segment has components, set the target effect to these as well
	      if(segment->HasComponents()) {
	        for(auto component : segment->GetComponentSegments()) {
	          component->SetTargetEffectNum(this->NumTargetEffects());
	        }
	      }
	    }
	  }
	}
      }
    }
  }
  if(line!="</targetInt>") return -1;
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(segment->GetEntranceKey()));
    PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(segment->GetExitKey()));
    double cmConversion;
    if(entrancePair->GetPType()==20) cmConversion = (exitPair->GetM(1)+exitPair->GetM(2))/exitPair->GetM(2);
    else cmConversion = entrancePair->GetM(2)/(entrancePair->GetM(1)+entrancePair->GetM(2));
    
    if(segment->IsTargetEffect()) {
      TargetEffect *targetEffect = this->GetTargetEffect(segment->GetTargetEffectNum());
      double sigma = targetEffect->GetSigma();
      targetEffect->SetSigma(cmConversion*sigma);
      
      for(EPointIterator point=segment->GetPoints().begin();point<segment->GetPoints().end();point++) {
	      point->SetTargetEffectNum(segment->GetTargetEffectNum());
	      
        if(targetEffect->IsTargetIntegration()||targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
	        double forwardDepth=0.0;
	        double backwardDepth=0.0;
	  
          if(targetEffect->IsTargetIntegration()) {
            double totalM=entrancePair->GetM(1)+entrancePair->GetM(2);
            double targetThickness = cmConversion*targetEffect->TargetThickness(point->GetLabEnergy(),configure);
            point->SetTargetThickness(targetThickness);
            if(targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
              if(targetEffect->IsConvCoefficients()) {
                backwardDepth=targetThickness+targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure)*5.0;
                forwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure)*5.0;
              }
              else {
                backwardDepth=targetThickness+targetEffect->convolutionRange*targetEffect->GetSigma()*5.0;
                forwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma()*5.0;
              }
              // Add straggling range extension at the deepest layer
              if(targetEffect->IsStraggling()) {
                double targetThickness_keV = targetThickness * 1000.0;  // Convert MeV to keV
                double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                double stragglingSigma = stragglingSigma_keV / 1000.0;  // Convert back to MeV
                double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                backwardDepth += stragglingRange;
              }
            }
            else {
              backwardDepth=targetThickness;
              forwardDepth=0.0;
              // Add straggling range extension for target integration only
              if(targetEffect->IsStraggling()) {
                double targetThickness_keV = targetThickness * 1000.0;
                double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                double stragglingSigma = stragglingSigma_keV / 1000.0;
                double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                backwardDepth += stragglingRange;
              }
            }
	        } 
    
          else if(targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
            if(targetEffect->IsConvCoefficients()){
              backwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure);
              forwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure);
            }
            else {
              backwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma();
              forwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma();
            }
          }

          // Generate adaptive integration grid
          double startEnergy = point->GetCMEnergy() + forwardDepth;
          double endEnergy ;

          // Check if target thickness is zero (for convolution-only cases)
          if(targetEffect->IsTargetIntegration()) {
            double targetThickness = point->GetTargetThickness();
            if(targetThickness < 1.0e-10) {
              // Target thickness is effectively zero - integrate down to 0.001 MeV
              endEnergy = 0.001;
              // Set density to 1e24 to prevent division issues
              targetEffect->SetDensity(1.0e24);
            } else {
              // Normal case - use backward depth
              endEnergy = point->GetCMEnergy() - backwardDepth;
              // Safety check: if backwardDepth > energy, set endEnergy to minimum
              if(endEnergy < 0.001) {
                endEnergy = 0.001;
              }
            }
          } else {
            // Convolution or ConvCoefficients - use backward depth
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if(endEnergy < 0.001) {
              endEnergy = 0.001;
            }
          }

          std::vector<double> energyGrid;
          int numPoints = targetEffect->NumSubPoints();
          if(configure.useAdaptiveGrid) {
            AdaptiveIntegrationGrid::GridConfig gridConfig;
            gridConfig.maxPoints = numPoints;
            gridConfig.entranceKey = segment->GetEntranceKey();
            gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
            gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
            gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
            AdaptiveIntegrationGrid gridGenerator(gridConfig);
            energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
          } else {
            double step = (startEnergy - endEnergy) / numPoints;
            for(int i = 0; i <= numPoints; i++)
              energyGrid.push_back(startEnergy - i * step);
          }

          for(size_t i=0; i<energyGrid.size(); i++) {
            double subEnergy = energyGrid[i];
            EPoint subPoint(point->GetCMAngle(),subEnergy,&*segment);
            if(targetEffect->IsTargetIntegration()) {
              double stoppingPower=cmConversion*targetEffect->GetStoppingPowerEq()->Evaluate(configure,subEnergy/cmConversion);
              subPoint.SetStoppingPower(stoppingPower);
            }
            point->AddSubPoint(subPoint);
          }

	      }
      }
    }
    if( segment->HasComponents() ) {
      for(auto component : segment->GetComponentSegments()) {
        if(component->IsTargetEffect()) {
          TargetEffect *targetEffect = this->GetTargetEffect(component->GetTargetEffectNum());
          double sigma = targetEffect->GetSigma();
          targetEffect->SetSigma(cmConversion*sigma);
          
          for(EPointIterator point=component->GetPoints().begin();point<component->GetPoints().end();point++) {
            point->SetTargetEffectNum(component->GetTargetEffectNum());
            
            if(targetEffect->IsTargetIntegration()||targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
              double forwardDepth=0.0;
              double backwardDepth=0.0;
        
              if(targetEffect->IsTargetIntegration()) {
                double totalM=entrancePair->GetM(1)+entrancePair->GetM(2);
                double targetThickness = cmConversion*targetEffect->TargetThickness(point->GetLabEnergy(),configure);
                point->SetTargetThickness(targetThickness);
                if(targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
                  if(targetEffect->IsConvCoefficients()) {
                    backwardDepth=targetThickness+targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure)*5.0;
                    forwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure)*5.0;
                  }
                  else {
                    backwardDepth=targetThickness+targetEffect->convolutionRange*targetEffect->GetSigma()*5.0;
                    forwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma()*5.0;
                  }
                  // Add straggling range extension at the deepest layer (component segments)
                  if(targetEffect->IsStraggling()) {
                    double targetThickness_keV = targetThickness * 1000.0;  // Convert MeV to keV
                    double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                    double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                    double stragglingSigma = stragglingSigma_keV / 1000.0;  // Convert back to MeV
                    double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                    backwardDepth += stragglingRange;
                  }
                }
                else {
                  backwardDepth=targetThickness;
                  forwardDepth=0.0;
                  // Add straggling range extension for target integration only (component segments)
                  if(targetEffect->IsStraggling()) {
                    double targetThickness_keV = targetThickness * 1000.0;
                    double stragglingCoeff = targetEffect->GetStragglingCoefficient();
                    double stragglingSigma_keV = stragglingCoeff * std::sqrt(targetThickness_keV);
                    double stragglingSigma = stragglingSigma_keV / 1000.0;
                    double stragglingRange = targetEffect->convolutionRange * stragglingSigma;
                    backwardDepth += stragglingRange;
                  }
                }
              }

              else if(targetEffect->IsConvolution()||targetEffect->IsConvCoefficients()) {
                if(targetEffect->IsConvCoefficients()){
                  backwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure);
                  forwardDepth=targetEffect->convolutionRange*targetEffect->CalculateSigma(point->GetLabEnergy(),configure);
                }
                else {
                  backwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma();
                  forwardDepth=targetEffect->convolutionRange*targetEffect->GetSigma();
                }
              }

              // Generate adaptive integration grid for component segment
              double startEnergy = point->GetCMEnergy() + forwardDepth;
              double endEnergy;

              // Check if target thickness is zero (for convolution-only cases)
              if(targetEffect->IsTargetIntegration()) {
                double targetThickness = point->GetTargetThickness();
                if(targetThickness < 1.0e-10) {
                  // Target thickness is effectively zero - integrate down to 0.001 MeV
                  endEnergy = 0.001;
                  // Set density to 1e24 to prevent division issues
                  targetEffect->SetDensity(1.0e24);
                } else {
                  // Normal case - use backward depth
                  endEnergy = point->GetCMEnergy() - backwardDepth;
                  // Safety check: if backwardDepth > energy, set endEnergy to minimum
                  if(endEnergy < 0.001) {
                    endEnergy = 0.001;
                  }
                }
              } else {
                // Convolution or ConvCoefficients - use backward depth
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if(endEnergy < 0.001) {
                  endEnergy = 0.001;
                }
              }

              std::vector<double> energyGrid;
              int numPoints = targetEffect->NumSubPoints();
              if(configure.useAdaptiveGrid) {
                AdaptiveIntegrationGrid::GridConfig gridConfig;
                gridConfig.maxPoints = numPoints;
                gridConfig.entranceKey = segment->GetEntranceKey();
                gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
                gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
                gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
                AdaptiveIntegrationGrid gridGenerator(gridConfig);
                energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
              } else {
                double step = (startEnergy - endEnergy) / numPoints;
                for(int i = 0; i <= numPoints; i++)
                  energyGrid.push_back(startEnergy - i * step);
              }

              for(size_t i=0; i<energyGrid.size(); i++) {
                double subEnergy = energyGrid[i];
                EPoint subPoint(point->GetCMAngle(),subEnergy,&*component);
                if(targetEffect->IsTargetIntegration()) {
                  double stoppingPower=cmConversion*targetEffect->GetStoppingPowerEq()->Evaluate(configure,subEnergy/cmConversion);
                  subPoint.SetStoppingPower(stoppingPower);
                }
                point->AddSubPoint(subPoint);
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

/*!
 * Returns true if the data is to be fit, otherwise returns false.  Used in the AZURECalc function class
 * to determine if a clone of the CNuc and EData objects should be made for thread safety.
 */

bool EData::IsFit() const {
  return isFit_;
}

/*!
 * Returns true if the call to function is for error analysis via Minos, otherwise returns false.  Used in the AZURECalc function
 * class to suppress transformation and file output during error analysis.
 */

bool EData::IsErrorAnalysis() const {
  return isErrorAnalysis_;
}

/*!
 * Sets the boolean indicating if the data is to be fit by AZURECalc function class. Used in AZUREMain function
 * class before calls to Minuit and AZURECalc.
 */

/*!
 * Returns true if the specified segment key exists corresponds to a segment in the ESegment vector, 
 * otherwise returns false.
 */

bool EData::IsSegmentKey(int segmentKey) {
  bool isKey=false;
  for (ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    if(segment->GetSegmentKey()==segmentKey) {
      isKey=true;
      break;
    }
  }
  return isKey;
}

/*!
 * Sets an internal variable specifying if the data is to be fit by Minuit.  Needed to determine cloning behavior in AZURECalc for 
 * thread safety.
 */

void EData::SetFit(bool fit) {
  isFit_=fit;
}

/*!
 * Sets the boolean indicating if the call to the function is for error analysis via Minos.
 */

void EData::SetErrorAnalysis(bool errorAnalysis) {
  isErrorAnalysis_=errorAnalysis;
}

/*!
 * This function updates the number of fit iterations per iteration during the fitting process.
 */

void EData::Iterate(){
  iterations_++;
}

/*!
 * This function sets the number of iterations to zero.
 */

void EData::ResetIterations(){
  iterations_=0;
}

/*!
 * This function is identical in role to the EPoint::Initialize function, except that it initializes
 * and entire EData object instead of a single EPoint object.
 */

int EData::Initialize(CNuc *compound,const Config &configure) {
  //Calculate channel lo-matrix and channel penetrability for each channel at each local energy
  configure.outStream << "Calculating Lo-Matrix, Phases, and Penetrabilities..." << std::endl;
  if(this->CalcEDependentValues(compound,configure)==-1) return -1;
  if((configure.fileCheckMask|configure.screenCheckMask) & Config::CHECK_ENERGY_DEP) 
    this->PrintEDependentValues(configure,compound); 
  //Calculate legendre polynomials for each data point
  configure.outStream << "Calculating Legendre Polynomials..." << std::endl;
  this->CalcLegendreP(configure.maxLOrder,compound);
  if((configure.fileCheckMask|configure.screenCheckMask) & Config::CHECK_LEGENDRE) 
    this->PrintLegendreP(configure);

  //Calculate Coulomb Amplitudes
  configure.outStream << "Calculating Coulomb Amplitudes..." << std::endl;
  this->CalcCoulombAmplitude(compound);
  if((configure.fileCheckMask|configure.screenCheckMask) & Config::CHECK_COUL_AMPLITUDES) {
    this->PrintCoulombAmplitude(configure,compound);
  }

  //Calculate new ec amplitudes
  if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
    configure.outStream << "Calculating External Capture Amplitudes..." << std::endl;
    if(this->CalculateECAmplitudes(compound,configure)==-1) return -1;
  }
  
  //Initialize component segments - ensure their entrance/exit pairs are properly set up
  configure.outStream << "Initializing Component Segments..." << std::endl;
  if(this->InitializeComponentSegments(compound,configure)==-1) return -1;
  
  return 0;
}

/*!
 * Adds a segment to the ESegment vector.
 */

void EData::AddSegment(ESegment segment) {
  segments_.push_back(segment);
}

/*!
 * Prints the data point after the object is filled or points are created.
 */

void EData::PrintData(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if(configure.fileCheckMask & Config::CHECK_DATA) {
    std::string outfile=configure.checkdir+"data.chk";
    fbuffer.open(outfile.c_str(),std::ios::out);
    sbuffer = &fbuffer;
  } else if(configure.screenCheckMask & Config::CHECK_DATA) sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if(((configure.fileCheckMask & Config::CHECK_DATA)&&fbuffer.is_open())||
     (configure.screenCheckMask & Config::CHECK_DATA)) {
    out << std::endl
	<< "************************************" << std::endl
	<< "*            Segments              *" << std::endl
	<< "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
	<< std::setw(17) << "Segment Key #"
	<< std::setw(17) << "Entrance Key #"
	<< std::setw(13) << "Exit Key #"
	<< std::setw(12) << "Min Energy"
	<< std::setw(12) << "Max Energy"
	<< std::setw(11) << "Min Angle"
	<< std::setw(11) << "Max Angle"
	<< std::setw(25) << "Data File"
	<< std::endl; 
    for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
      out << std::setw(11) << segment-GetSegments().begin()+1
	  << std::setw(17) << segment->GetSegmentKey()
	  << std::setw(17) << segment->GetEntranceKey()
	  << std::setw(13) << segment->GetExitKey()
	  << std::setw(12) << segment->GetMinEnergy()
	  << std::setw(12) << segment->GetMaxEnergy()
	  << std::setw(11) << segment->GetMinAngle()
	  << std::setw(11) << segment->GetMaxAngle()
	  << std::setw(25) << segment->GetDataFile()
	  << std::endl;
    }  
    out << std::endl
	<< "************************************" << std::endl
	<< "*               Data               *" << std::endl
	<< "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
	<< std::setw(14) << "Data Point #"
	<< std::setw(15) << "Lab Energy"
	<< std::setw(15) << "CM Energy"
	<< std::setw(15) << "Angle"
	<< std::setw(20) << "Cross Section"
	<< std::setw(22) << "Cross Section Error"
        << std::setw(12) << "Map Point"
	<< std::setw(18) << "# of Subpoints"
	<< std::setw(18) << "Low Sub Energy"
	<< std::setw(18) << "High Sub Energy"
	<< std::endl;
    for(EDataIterator data=begin();data!=end();data++) {
      out << std::setw(11) << data.segment()-GetSegments().begin()+1
	  << std::setw(14) << data.point()-(data.segment()->GetPoints()).begin()+1
	  << std::setw(15) << data.point()->GetLabEnergy() 
	  << std::setw(15) << data.point()->GetCMEnergy() 
	  << std::setw(15) << data.point()->GetCMAngle()
	  << std::setw(20) << data.point()->GetCMCrossSection()
	  << std::setw(22) << data.point()->GetCMCrossSectionError();
      if(data.point()->IsMapped()){
	EnergyMap map=data.point()->GetMap();
	char tempMap[25];
	sprintf(tempMap,"(%d,%d)",map.segment,map.point);
	out << std::setw(12) <<  tempMap << std::endl;
      } else
	out << std::setw(12) << "Not Mapped"
	    << std::setw(18) << data.point()->NumSubPoints();
      if(data.point()->IsTargetEffect()&&
	 (data.point()->GetParentData()->GetTargetEffect(data.point()->GetTargetEffectNum())->IsConvolution()||
	  data.point()->GetParentData()->GetTargetEffect(data.point()->GetTargetEffectNum())->IsTargetIntegration())) {
	out << std::setw(18) << data.point()->GetSubPoint(data.point()->NumSubPoints())->GetCMEnergy() 
	    << std::setw(18) << data.point()->GetSubPoint(1)->GetCMEnergy();
      }
      out << std::endl;
      if(data.point()==data.segment()->GetPoints().end()-1) out << std::endl;
    }
  } else configure.outStream << "Could not write data check file." << std::endl;
  out.flush();
  if(fbuffer.is_open()) fbuffer.close();
}

/*!
 * Calls EPoint::CalcLegendreP for each point in the entire EData object.
 */

void EData::CalcLegendreP(int maxL,CNuc *theCNuc) { 
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) { 
    TargetEffect *effect = (segment->IsTargetEffect() && 
			    this->GetTargetEffect(segment->GetTargetEffectNum())->IsQCoefficients()) ?
      this->GetTargetEffect(segment->GetTargetEffectNum()) : NULL;
#pragma omp parallel for 
    for(int i=1;i<=segment->NumPoints();i++) {
      EPoint* point=segment->GetPoint(i);
      point->CalcLegendreP(maxL, theCNuc, effect);    
    }
  }
}

/*!
 * Prints the Legendre polynomials for each point in the EData object.
 */ 

void EData::PrintLegendreP(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if(configure.fileCheckMask & Config::CHECK_LEGENDRE) {
    std::string outfile=configure.checkdir+"legendre.chk";
    fbuffer.open(outfile.c_str(),std::ios::out);
    sbuffer = &fbuffer;
  } else if(configure.screenCheckMask & Config::CHECK_LEGENDRE) sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if(((configure.fileCheckMask & Config::CHECK_LEGENDRE)&&fbuffer.is_open())||
     (configure.screenCheckMask & Config::CHECK_LEGENDRE)) {
    out << std::endl
    << "************************************" << std::endl
    << "*       Legendre Polynomials       *" << std::endl
    << "************************************" << std::endl;
    out << std::setw(10) << "Segment #"
    << std::setw(10) << "Point #"
    << std::setw(15) << "CM Energy"
    << std::setw(15) << "Angle"
    << std::setw(5)  << "L"
    << std::setw(15) << "Leg. Poly." << std::endl;
    for(EDataIterator data=begin();data!=end();data++) {
      for(int lOrder=0;lOrder<=data.point()->GetMaxLOrder();lOrder++) {
	out << std::setw(10) << data.segment()-GetSegments().begin()+1
	    << std::setw(10) << data.point()-(data.segment()->GetPoints()).begin()+1
	    << std::setw(15) << data.point()->GetCMEnergy()
	    << std::setw(15) << data.point()->GetCMAngle()
	    << std::setw(5)  << lOrder
	    << std::setw(15) << data.point()->GetLegendreP(lOrder) << std::endl;	
      }
    }
  } else configure.outStream << "Could not write legendre polynomials check file." << std::endl;
  out.flush();
  if(fbuffer.is_open()) fbuffer.close();    
}

/*!
 * Calls EPoint::CalcEDependentValues for each point in the entire EData object.
 */

int EData::CalcEDependentValues(CNuc *theCNuc,const Config& configure) {
  for(ESegmentIterator segment=GetSegments().begin(); segment<GetSegments().end(); segment++) {
    bool localStop = false;
#pragma omp parallel for shared(localStop,configure)
    for(int i=1;i<=segment->NumPoints();i++) {
      if(configure.stopFlag||localStop) continue;
      EPoint *point = segment->GetPoint(i);
      if(!(point->IsMapped())){
	try {
	  point->CalcEDependentValues(theCNuc,configure);
	} catch(GSLException e) {
#pragma omp critical 
	  { 
	    std::cout<<point->GetLabEnergy()<<"\t"<<point->GetLabAngle()<<std::endl;
	    configure.outStream << e.what() << std::endl;
	    localStop = true;
	  }
	}
      }
    }
    if(configure.stopFlag||localStop) return -1;
  }
  return 0;
}

/*!
 * Prints the values calculated by EPoint::CalcEDependentValues for each point in the entire EData object.
 */

void EData::PrintEDependentValues(const Config &configure,CNuc *theCNuc) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if(configure.fileCheckMask & Config::CHECK_ENERGY_DEP) {
    std::string outfile=configure.checkdir+"lomatrixandpene.chk";
    fbuffer.open(outfile.c_str(),std::ios::out);
    sbuffer = &fbuffer;
  } else if(configure.screenCheckMask & Config::CHECK_ENERGY_DEP) sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if(((configure.fileCheckMask & Config::CHECK_ENERGY_DEP)&&fbuffer.is_open())||
     (configure.screenCheckMask & Config::CHECK_ENERGY_DEP)) {
    out << std::endl
    << "************************************" << std::endl
    << "*  Lo Matrix and Penetrabilities   *" << std::endl
    << "************************************" << std::endl;
    out << std::setw(10) << "Seg_#"
    << std::setw(10) << "Point_#"
    << std::setw(5)  << "j" 
    << std::setw(5)  << "ch" 
    << std::setw(5)  << "l"
    << std::setw(15) << "E_chan"
    << std::setw(15) << "sqrt_pene" 
    << std::setw(25) << "Lo" 
    << std::setw(25) << "expHSP" << std::endl;
    for(EDataIterator data=begin();data!=end();data++) {
      double inEnergy=data.point()->GetCMEnergy()
        +theCNuc->GetPair(theCNuc->GetPairNumFromKey(data.segment()->GetEntranceKey()))->GetSepE();
      for(int j=1;j<=theCNuc->NumJGroups();j++) {
	if(theCNuc->GetJGroup(j)->IsInRMatrix()) {
	  JGroup *theJGroup=theCNuc->GetJGroup(j);
	  for(int ch=1;ch<=theJGroup->NumChannels();ch++) {
	    AChannel *theChannel=theJGroup->GetChannel(ch);
	    PPair *thePair=theCNuc->GetPair(theChannel->GetPairNum());
	    int lValue=theChannel->GetL();
	    double localEnergy=inEnergy-thePair->GetSepE()-thePair->GetExE();
            complex loElement = data.point()->GetLoElement(j,ch);
            complex expHSP = data.point()->GetExpHardSpherePhase(j,ch);
	    out << std::setw(10) << data.segment()-GetSegments().begin()+1
		<< std::setw(10) << data.point()-(data.segment()->GetPoints()).begin()+1
		<< std::setw(5) << j 
		<< std::setw(5) << ch  
		<< std::setw(5) << lValue 
		<< std::setw(15) << localEnergy 
		<< std::setw(15) << data.point()->GetSqrtPenetrability(j,ch)
		<< std::setw(25) << loElement.real() << " " << loElement.imag() 
                << std::setw(25) << expHSP.real() << " " << expHSP.imag() << std::endl;
	  }           
	}
      }
    }   
  } else configure.outStream << "Could not write lo-matrix and penetrabilities check file." << std::endl;
  out.flush();
  if(fbuffer.is_open()) fbuffer.close();    
}

/*!
 * Calls EPoint::CalcCoulombAmplitude for each point in the entire EData object.
 */

void EData::CalcCoulombAmplitude(CNuc *theCNuc) {
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
#pragma omp parallel for
    for(int i=1;i<=segment->NumPoints();i++) {
      EPoint* point = segment->GetPoint(i);
      point->CalcCoulombAmplitude(theCNuc);
    }
  }
}

/*!
 * Prints the values calculated by EPoint::CalcCoulombAmplitude for each point in the entire EData object.
 */

void EData::PrintCoulombAmplitude(const Config &configure,CNuc *theCNuc) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if(configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES) {
    std::string outfile=configure.checkdir+"coulombamplitudes.chk";
    fbuffer.open(outfile.c_str(),std::ios::out);
    sbuffer = &fbuffer;
  } else if(configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES) sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if(((configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES)&&fbuffer.is_open())||
     (configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES)) {
    out << std::endl
	<< "************************************" << std::endl
	<< "*        Coulomb Amplitudes        *" << std::endl
	<< "************************************" << std::endl;
    out << std::setw(10) << "segment #"
	<< std::setw(10) << "point #"
	<< std::setw(10) << "aa"
	<< std::setw(15) << "cmenergy"
	<< std::setw(15) << "angle"
	<< std::setw(25) << "coulomb amplitude"
	<< std::endl;
    for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
      if(segment->GetEntranceKey()==segment->GetExitKey()) {
	for(EPointIterator point=segment->GetPoints().begin();point<segment->GetPoints().end();point++) {
	  out << std::setw(10) << segment-GetSegments().begin()+1
	      << std::setw(10) << point-segment->GetPoints().begin()+1
	      << std::setw(10) << theCNuc->GetPairNumFromKey(segment->GetEntranceKey())
	      << std::setw(15) << point->GetCMEnergy()
	      << std::setw(15) << point->GetCMAngle()
	      << std::setw(25) << point->GetCoulombAmplitude()
	      << std::endl;
	}
      }
    }
  } else configure.outStream << "Could not write coulomb amplitudes check file." << std::endl;
  out.flush();
  if(fbuffer.is_open()) fbuffer.close();
}

/*!
 * Writes the output files for the calculation.  The output files are all in center of mass frame, and contain
 * columns for energy, angle, calculated cross section, calculated s-factor, experimental cross section and error
 * and experimental s-factor and error.
 */

void EData::WriteOutputFiles(const Config &configure, bool isFit) {
  AZUREOutput output(configure.outputdir);
  std::ofstream chiOut;
  if(!isFit&&(configure.paramMask & Config::CALCULATE_WITH_DATA)) {
    std::string chiOutFile = configure.outputdir+"chiSquared.out";
    chiOut.open(chiOutFile.c_str());
  }
  chiOut << "Segment#," << " Chi-Squared, " << " N, " << " Norm, " << " Norm-Chi-Squared " << std::endl; 
  
  if(!(configure.paramMask & Config::CALCULATE_WITH_DATA)) output.SetExtrap();
  bool isVaryNorm=false;
  double totalChiSquared=0.;
  double totalNormChiSquared=0.;
  double totalN=0.;
  //use kinflag.dat to control output option
  std::ifstream kinoption;
  kinoption.open("kinflag.dat");
  int kinflag=0; // 0 cm output, 1 lab output
  if(kinoption){
    kinoption >> kinflag;
    kinoption.close();
  }
  if(kinflag!=0) configure.outStream<<"Using alternate output format..."<<std::endl;
  ESegmentIterator firstSumIterator = GetSegments().end();
  for(ESegmentIterator segment=GetSegments().begin();
      segment<GetSegments().end();segment++) {
    if(segment->IsTotalCapture()) {
      firstSumIterator=segment;
      // With component segments, no need to skip - total capture is handled within the segment
    } 
    if(segment->IsVaryNorm()) isVaryNorm=true;
    int aa=segment->GetEntranceKey();
    int ir=segment->GetExitKey();
    std::filebuf* buf;
    if(firstSumIterator!=GetSegments().end()) buf=output(aa,-1);
    else {
      if(segment->IsAngularDist()&&
	 !(configure.paramMask & Config::CALCULATE_WITH_DATA)) buf=output(aa,ir,true);
      else buf=output(aa,ir);
    }
    std::ostream out(buf);	
    ESegmentIterator thisSegment = segment;
    if(firstSumIterator!=GetSegments().end()) thisSegment = firstSumIterator;
     
    if(kinflag==0){
      for(EPointIterator point=segment->GetPoints().begin();point<segment->GetPoints().end();point++) {
        out.precision(10);
        if(segment->IsAngularDist()) {
	  out << std::setw(18) << std::scientific << point->GetCMEnergy();
	  for(int i = 0;i<point->GetNumAngularDists();i++) out << std::setw(18) << point->GetAngularDist(i);
	  out << std::endl;
        } else {
	  double fitCrossSection=point->GetFitCrossSection();
	  if(firstSumIterator!=GetSegments().end()) {
	    int pointIndex=point-segment->GetPoints().begin()+1;
	    for(ESegmentIterator it=firstSumIterator;it<segment;it++)
	      fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
	  }
	  out << std::setw(18) << std::scientific << point->GetCMEnergy()
	      << std::setw(18) << std::scientific << point->GetExcitationEnergy()
	      << std::setw(18) << std::scientific << point->GetCMAngle()
	      << std::setw(18) << std::scientific << fitCrossSection
	      << std::setw(18) << std::scientific << fitCrossSection*point->GetSFactorConversion();
	  if(!output.IsExtrap()) {
	    double dataNorm=thisSegment->GetNorm();
	    out << std::setw(18) << std::scientific << point->GetCMCrossSection()*dataNorm
	        << std::setw(18) << std::scientific << point->GetCMCrossSectionError()*dataNorm
	        << std::setw(18) << std::scientific << point->GetCMCrossSection()*dataNorm*point->GetSFactorConversion()
	        << std::setw(18) << std::scientific << point->GetCMCrossSectionError()*dataNorm*point->GetSFactorConversion()
	        << std::endl;
	  } else out << std::endl;
        }
      }
    }
    if(kinflag==1){
      for(EPointIterator point=segment->GetPoints().begin();point<segment->GetPoints().end();point++) {
        out.precision(10);
        if(segment->IsAngularDist()) {
	  out << std::setw(18) << std::scientific << point->GetCMEnergy();
	  for(int i = 0;i<point->GetNumAngularDists();i++) out << std::setw(18) << point->GetAngularDist(i);
	  out << std::endl;
        } else {
	  double fitCrossSection=point->GetFitCrossSection()/point->GetCrossSectionKinFactor();
	  if(firstSumIterator!=GetSegments().end()) {
	    int pointIndex=point-segment->GetPoints().begin()+1;
	    for(ESegmentIterator it=firstSumIterator;it<segment;it++)
	      fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
	  }
	  out << std::setw(18) << std::scientific << point->GetLabEnergy()
	      << std::setw(18) << std::scientific << point->GetExcitationEnergy()
	      << std::setw(18) << std::scientific << point->GetLabAngle()
	      << std::setw(18) << std::scientific << fitCrossSection
	      << std::setw(18) << std::scientific << fitCrossSection*point->GetSFactorConversion();
	  if(!output.IsExtrap()) {
	    double dataNorm=thisSegment->GetNorm();
	    out << std::setw(18) << std::scientific << point->GetLabCrossSection()*dataNorm
	        << std::setw(18) << std::scientific << point->GetLabCrossSectionError()*dataNorm
	        << std::setw(18) << std::scientific << point->GetLabCrossSection()*dataNorm*point->GetSFactorConversion()
	        << std::setw(18) << std::scientific << point->GetLabCrossSectionError()*dataNorm*point->GetSFactorConversion()
	        << std::endl;
	  } else out << std::endl;
        }
      }
    }
    if(!isFit&&(configure.paramMask & Config::CALCULATE_WITH_DATA)) {
      totalChiSquared+=(thisSegment->GetSegmentChiSquared());
      totalNormChiSquared+=1;
      totalN+=thisSegment->NumPoints();
      chiOut << thisSegment->GetSegmentKey()
             << "," 
	     << thisSegment->GetSegmentChiSquared()
             << ","
             << thisSegment->NumPoints()
             << ","
             << thisSegment->GetNorm()
             << ","
             //<< thisSegment->GetNormChiSquared()
	     << std::endl;
    }
    out<<std::endl<<std::endl;out.flush();
    firstSumIterator=GetSegments().end();
  }
  
  
  if(!isFit&&(configure.paramMask & Config::CALCULATE_WITH_DATA)) {
    chiOut << "Total-Chi-Squared: " 
	      << totalChiSquared
              << " Total-Norm-Chi-Squared: "
              << totalNormChiSquared 
              << " Total-N: "
              << totalN
              << std::endl << std::endl;
    chiOut.flush();chiOut.close();
  }   
  if(isVaryNorm) {
    std::string outputfile=configure.outputdir+"normalizations.out";
    std::ofstream out(outputfile.c_str());
    if(out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
	if(segment->IsVaryNorm()) out << segment->GetSegmentKey() << ","
                                      << segment->GetDataFile() << ","
				      << segment->GetMinAngle() << ","
				      << segment->GetMaxAngle() << ",";
                                      if(segment->IsCMDifferential()==true){
                                         out << "CM,";
                                      } else out << "Lab,";
                                  out << segment->GetMinEnergy() << ","
                                      << segment->GetMaxEnergy() << ",";
                                      if(segment->IsCMDifferential()==true){
                                         out << "CM,";
                                      } else out << "Lab,";
				  out << segment->GetNorm() << ","
                                      << segment->GetEnergyShift() << std::endl;
      }
      out.flush();
      out.close();
    } else configure.outStream << "Could not write normalization file." << std::endl;
  }
  
  // Check if any energy shifts are varied
  bool isVaryEnergyShift=false;
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    if(segment->IsVaryEnergyShift()) {
      isVaryEnergyShift=true;
      break;
    }
  }
  
  // Write energy shifts file
  if(isVaryEnergyShift) {
    std::string outputfile=configure.outputdir+"shifts.out";
    std::ofstream out(outputfile.c_str());
    if(out) {
      out.precision(4);
      out << std::scientific;
      out << "segment_key_#, file_name, low_angle_bound, high_angle_bound, aframe, low_energy_bound, high_energy_bound, eframe, norm, shift" << std::endl;
      for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
	if(segment->IsVaryEnergyShift()) out << segment->GetSegmentKey() << ","
                                             << segment->GetDataFile() << ","
			             << segment->GetMinAngle() << ","
			             << segment->GetMaxAngle() << ",";
                                             if(segment->IsCMDifferential()==true){
                                                out << "CM,";
                                             } else out << "Lab,";
                                         out << segment->GetMinEnergy() << ","
                                             << segment->GetMaxEnergy() << ",";
                                             if(segment->IsCMDifferential()==true){
                                                out << "CM,";
                                             } else out << "Lab,";
		 out << segment->GetNorm() << ","
                                             << segment->GetEnergyShift() << std::endl;
      }
      out.flush();
      out.close();
    } else configure.outStream << "Could not write energy shifts file." << std::endl;
  }
}

/*!
 * If external capture amplitudes are to be calculated, EPoint::CalculateECAmplitudes is
 * called for each point with a corresponding external capture component in the EData object.
 * Otherwise, the amplitudes are read from the specified file.
 */

int EData::CalculateECAmplitudes(CNuc *theCNuc,const Config& configure) {
  std::ifstream in;
  std::ofstream out;
  std::string outputfile;
  if(configure.paramMask & Config::CALCULATE_WITH_DATA) outputfile=configure.outputdir+"intEC.dat";
  else outputfile=configure.outputdir+"intEC.extrap";
  if(configure.paramMask & Config::USE_PREVIOUS_INTEGRALS) in.open(configure.integralsfile.c_str());
  else {
    out.open(outputfile.c_str());
    if(!out) configure.outStream << "Could not write to EC Amplitude File." << std::endl;
  }
  int sumSegmentI=0;
  int numSumSegments=0;
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    if(segment->IsTotalCapture()) {
      numSumSegments = segment->IsTotalCapture();
      sumSegmentI=1; // Start at 1 since total capture is now handled within the segment
    } else {
      sumSegmentI=0;
      numSumSegments=0;
    }
    char segmentKeyOut[256];
    if(segment->IsTotalCapture()) sprintf(segmentKeyOut,"%d (Total: %d)",segment->GetSegmentKey(),numSumSegments);
    else sprintf(segmentKeyOut,"%d",segment->GetSegmentKey());
    int aa=theCNuc->GetPairNumFromKey(segment->GetEntranceKey());
    if(theCNuc->GetPair(aa)->GetPType()==20) continue;
    if(theCNuc->GetPair(aa)->IsEntrance()) {
      PPair *entrancePair=theCNuc->GetPair(aa);
      for(int j=1;j<=theCNuc->NumJGroups();j++) {
	for(int la=1;la<=theCNuc->GetJGroup(j)->NumLevels();la++) {
	  if(theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
	    ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
	    int ir=theCNuc->GetPairNumFromKey(segment->GetExitKey());
	    if(ecLevel->GetECPairNum()==ir) {
	      if(!(configure.paramMask & Config::USE_PREVIOUS_INTEGRALS)) {
		configure.outStream << "\tSegment #" << std::setw(12) << segmentKeyOut 
		          << std::setw(0) << " [                         ] 0%";configure.outStream.flush();
		int numPoints=segment->NumPoints();
		int pointIndex=0;
		time_t startTime = time(NULL);
		bool localStop = false;
#pragma omp parallel for shared(configure,localStop) 
		for(int i=1;i<=numPoints;i++) {
		  if(configure.stopFlag||localStop) continue;
		  EPoint *point = segment->GetPoint(i);
		  if(!(point->IsMapped())) {
		    try {
		      point->CalculateECAmplitudes(theCNuc,configure);
		    } catch(GSLException e) {
#pragma omp critical
		      {
			configure.outStream << e.what() << std::endl;
			localStop=true;
		      }
		    }
		  }
		  ++pointIndex;
		  if(difftime(time(NULL),startTime)>0.25) {
		    startTime=time(NULL);
		    std::string progress=" [";
		    double percent=0.;
		    for(int j = 1;j<=25;j++) {
		      if(pointIndex>=percent*numPoints&&percent<1.) {
			percent+=0.04;
			progress+='*';
		      } else progress+=' ';
		    } progress+="] ";
		    configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut 
					<< std::setw(0) << progress << percent*100 << '%';configure.outStream.flush();
		  }
		}
		if(configure.stopFlag||localStop) {
		  if(out.is_open()) out.close();
		  if(in.is_open()) in.close();
		  return -1;
		}
		configure.outStream << "\r\tSegment #" << std::setw(12) << segmentKeyOut 
				    << std::setw(0) << " [*************************] 100%" << std::endl;
	      }
	      for(EPointIterator point=segment->GetPoints().begin();
		  point<segment->GetPoints().end();point++) {
		if(!(point->IsMapped())) {
		  for(int k=1;k<=entrancePair->GetDecay(ir)->NumKGroups();k++) {
		    for(int ecm=1;ecm<=entrancePair->GetDecay(ir)->GetKGroup(k)->NumECMGroups();ecm++) {
		      if(!(configure.paramMask & Config::USE_PREVIOUS_INTEGRALS)) {
			if(out.is_open()) out << point->GetECAmplitude(k,ecm) << std::endl;
			for(EPointIterator subPoint=point->GetSubPoints().begin();
			    subPoint<point->GetSubPoints().end();subPoint++)
			  if(out.is_open()) out << subPoint->GetECAmplitude(k,ecm) << std::endl;
		      } else {
			complex ecAmplitude(0.0,0.0);
			in >> ecAmplitude;
			point->AddECAmplitude(k,ecm,ecAmplitude);
			for(EPointIterator subPoint=point->GetSubPoints().begin();
			    subPoint<point->GetSubPoints().end();subPoint++) {
			  ecAmplitude=complex(0.0,0.0);
			  in >> ecAmplitude;
			  subPoint->AddECAmplitude(k,ecm,ecAmplitude);
			}
		      }
		      for(EPointMapIterator mappedPoint=point->GetMappedPoints().begin();
			  mappedPoint<point->GetMappedPoints().begin();mappedPoint++) {
			(*mappedPoint)->AddECAmplitude(k,ecm,point->GetECAmplitude(k,ecm));
			for(int i=1;i<=point->NumSubPoints();i++) {
			  (*mappedPoint)->GetSubPoint(i)->
			    AddECAmplitude(k,ecm,point->GetSubPoint(i)->GetECAmplitude(k,ecm));
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
    if(sumSegmentI==numSumSegments) {
      sumSegmentI=0;
      numSumSegments=0;
    }
  }
  if(out.is_open()) {
    out.flush();
    out.close();
  }
  if(in.is_open()) in.close();
  return 0;
}

/*!
 * Initialize component segments by ensuring their entrance/exit pairs are properly set up
 * in the compound nucleus for calculations.
 */

int EData::InitializeComponentSegments(CNuc *theCNuc, const Config& configure) {
  // Initialize all component segments with the same process as regular segments
  configure.outStream << "Initializing Component Segments (" << componentSegments_.size() << " components)..." << std::endl;
  
  // First, ensure entrance/exit pairs are properly set up in the compound nucleus
  for(auto& componentSegment : componentSegments_) {
    int entranceKey = componentSegment.GetEntranceKey();
    int exitKey = componentSegment.GetExitKey();
    
    // Ensure entrance pair is set as entrance if it exists
    if(theCNuc->IsPairKey(entranceKey)) {
      int pairNum = theCNuc->GetPairNumFromKey(entranceKey);
      theCNuc->GetPair(pairNum)->SetEntrance();
    }
    
    // Check if exit pair is valid (either -1 for total capture or existing pair)
    bool isValidExit = false;
    if(exitKey == -1) {
      // Check for total capture validity
      for(int i = 1; i <= theCNuc->NumPairs(); i++) {
        if(theCNuc->GetPair(i)->GetPType() == 10) {
          isValidExit = true;
          break;
        }
      }
    } else {
      isValidExit = theCNuc->IsPairKey(exitKey);
    }
    
    if(!isValidExit) {
      configure.outStream << "Warning: Invalid exit key " << exitKey 
                         << " for component segment" << std::endl;
    }
  }
  
  // Initialize component segments using the same process as regular segments
  // Apply pair-dependent conversions for component segments
  for(auto& componentSegment : componentSegments_) {
    int entranceKey = componentSegment.GetEntranceKey();
    int exitKey = componentSegment.GetExitKey();

    PPair *entrancePair = nullptr;
    PPair *exitPair = nullptr;

    if(theCNuc->IsPairKey(entranceKey)) {
      entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(entranceKey));
    }
    if(theCNuc->IsPairKey(exitKey)) {
      exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(exitKey));
    }

    // Apply the same conversions that are done for test segment points
    for(int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint *point = componentSegment.GetPoint(i);

      if(entrancePair && exitPair) {
        if(entrancePair->GetPType()==20) {
          point->ConvertDecayEnergy(exitPair);
        } else if(!componentSegment.IsCMDifferential()) {
          point->ConvertLabEnergy(entrancePair);
        } else if(componentSegment.IsCMDifferential()) {
          point->ConvertLabEnergy(entrancePair);
        }

        if(exitPair->GetPType()==0 && componentSegment.IsDifferential() &&
           !componentSegment.IsPhase() && !componentSegment.IsAngularDist() && !componentSegment.IsCMDifferential()) {
          if(entranceKey == exitKey) {
            point->ConvertLabAngle(entrancePair);
          } else {
            point->ConvertLabAngle(entrancePair, exitPair, configure);
          }
          point->ConvertCrossSection(entrancePair, exitPair);
        }

        if(exitPair->GetPType()==10 && componentSegment.IsDifferential() &&
           !componentSegment.IsPhase() && !componentSegment.IsAngularDist() && !componentSegment.IsCMDifferential()) {
          point->ConvertLabAngleGammas(entrancePair);
          point->ConvertCrossSectionGammas(entrancePair);
        }
      }
    }
  }

  // CalcEDependentValues for component segments
  for(auto& componentSegment : componentSegments_) {
    bool localStop = false;
    for(int i = 1; i <= componentSegment.NumPoints(); i++) {
      if(configure.stopFlag || localStop) continue;
      EPoint *point = componentSegment.GetPoint(i);
      if(!(point->IsMapped())) {
        try {
          point->CalcEDependentValues(theCNuc, configure);
        } catch(GSLException e) {
          configure.outStream << "Component segment calculation error: " << e.what() << std::endl;
          localStop = true;
        }
      }
    }
    if(configure.stopFlag || localStop) return -1;
  }
  
  // CalcLegendreP for component segments
  for(auto& componentSegment : componentSegments_) {
    // Get TargetEffect with Q-coefficients if applicable
    TargetEffect* effect = NULL;
    if(componentSegment.IsTargetEffect()) {
      TargetEffect* te = this->GetTargetEffect(componentSegment.GetTargetEffectNum());
      if(te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    for(int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint* point = componentSegment.GetPoint(i);
      point->CalcLegendreP(configure.maxLOrder, theCNuc, effect);
    }
  }
  
  // CalcCoulombAmplitude for component segments
  for(auto& componentSegment : componentSegments_) {
    for(int i = 1; i <= componentSegment.NumPoints(); i++) {
      EPoint* point = componentSegment.GetPoint(i);
      point->CalcCoulombAmplitude(theCNuc);
    }
  }
  
  // Calculate EC amplitudes for component segments if external capture is enabled
  if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
    configure.outStream << "Calculating EC Amplitudes for Component Segments..." << std::endl;
    for(auto& componentSegment : componentSegments_) {
      int aa = theCNuc->GetPairNumFromKey(componentSegment.GetEntranceKey());
      if(theCNuc->GetPair(aa)->GetPType() != 20 && theCNuc->GetPair(aa)->IsEntrance()) {
        PPair *entrancePair = theCNuc->GetPair(aa);
        for(int j = 1; j <= theCNuc->NumJGroups(); j++) {
          for(int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
            if(theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
              ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
              int ir = theCNuc->GetPairNumFromKey(componentSegment.GetExitKey());
              if(ecLevel->GetECPairNum() == ir) {
                for(int i = 1; i <= componentSegment.NumPoints(); i++) {
                  EPoint *point = componentSegment.GetPoint(i);
                  if(!(point->IsMapped())) {
                    try {
                      point->CalculateECAmplitudes(theCNuc, configure);
                    } catch(...) {
                      configure.outStream << "Warning: EC amplitude calculation failed for component segment point" << std::endl;
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

/*!
 * Creates a component segment as a full copy of the base segment with different entrance/exit keys
 */
ESegment* EData::CreateComponentSegment(const ESegment& baseSegment, int entranceKey, int exitKey) {
  // Create a copy of the base segment with different entrance/exit keys
  // Store component segments separately to avoid iterator invalidation

  // Add a copy of the base segment to the component segments vector
  this->componentSegments_.push_back(baseSegment);
  ESegment* componentSegment = &componentSegments_.back();

  // Set the new entrance/exit keys for the component segment
  componentSegment->SetEntranceKey(entranceKey);
  componentSegment->SetExitKey(exitKey);

  // Set a unique segment key to avoid cache conflicts
  // Use negative keys for component segments to distinguish from regular segments
  static int componentSegmentKeyCounter = -1000;
  componentSegment->SetSegmentKey(componentSegmentKeyCounter--);

  // Clear any existing components to avoid circular references
  componentSegment->ClearComponents();

  // The component segment now has the same data points as the base segment
  // but with different entrance/exit keys, so it will calculate different
  // theoretical cross sections when initialized

  return componentSegment;
}

/*!
 * Creates a component segment with a fixed angle for ratio denominators
 */
ESegment* EData::CreateComponentSegment(const ESegment& baseSegment, int entranceKey, int exitKey, double fixedAngle) {
  // Create a copy of the base segment with different entrance/exit keys
  this->componentSegments_.push_back(baseSegment);
  ESegment* componentSegment = &componentSegments_.back();

  // Set the new entrance/exit keys for the component segment
  componentSegment->SetEntranceKey(entranceKey);
  componentSegment->SetExitKey(exitKey);

  // Set the fixed angle by constraining min and max to the same value
  componentSegment->SetMinAngle(fixedAngle);
  componentSegment->SetMaxAngle(fixedAngle);

  // CRITICAL: Override all point angles to the fixed angle
  // This ensures calculations use the fixed angle, not the original data angles
  std::vector<EPoint>& points = componentSegment->GetPoints();
  for(auto& point : points) {
    point.SetCMAngle(fixedAngle);
    point.SetLabAngle(fixedAngle); // Also set lab angle to be consistent
  }

  // Set a unique segment key to avoid cache conflicts
  static int componentSegmentKeyCounter = -1000;
  componentSegment->SetSegmentKey(componentSegmentKeyCounter--);

  // Clear any existing components to avoid circular references
  componentSegment->ClearComponents();

  return componentSegment;
}

/*!
 * This function determined what points should be mapped to another to reduce
 * redundant calculations at like energies.
 */

void EData::MapData() {
  // Disable all point mapping - just clear existing mappings
  for(ESegmentIterator segment=GetSegments().begin(); segment<GetSegments().end(); segment++) {
    for(EPointIterator point=segment->GetPoints().begin(); point<segment->GetPoints().end(); point++) {
      point->ClearMapping();
      point->ClearLocalMappedPoints();
    }
  }
}

/*!
 * Adds a TargetEffect object to the vector contained within the present object.
 */

void EData::AddTargetEffect(TargetEffect targetEffect) {
  targetEffects_.push_back(targetEffect);
}

/*!
 * Sets the normalization parameter offset in the parameter vector.
 */  

void EData::SetNormParamOffset(int offset) {
  normParamOffset_=offset;
}

/*!
 * Sets the energy shift parameter offset in the parameter vector.
 */  

void EData::SetEnergyShiftParamOffset(int offset) {
  energyShiftParamOffset_=offset;
}

/*!
 * Fills the Minuit parameter array from initial values in the EData object.
 */

void EData::FillMnParams(ROOT::Minuit2::MnUserParameters &p) {
  SetNormParamOffset(p.Params().size());
  char varname[50];
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    if(segment->IsVaryNorm()) {
      sprintf(varname,"segment_%d_norm",segment->GetSegmentKey());
      p.Add(varname,segment->GetNorm(),segment->GetNorm()*0.05);
      p.SetLowerLimit(varname,0.0);
      
      // Parameter settings will be applied by ParameterLimitsManager during fit
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
  
  // Add energy shift parameters FOR ALL SEGMENTS (like normalization)
  SetEnergyShiftParamOffset(p.Params().size());
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    // Always add energy shift parameter, regardless of IsVaryEnergyShift()
    sprintf(varname,"segment_%d_energy_shift",segment->GetSegmentKey());
    double stepSize = (segment->GetEnergyShiftError() > 0.0) ? segment->GetEnergyShiftError() * 0.1 : 0.001;
    p.Add(varname,segment->GetEnergyShift(),stepSize);

    if(!segment->IsVaryEnergyShift()) {
      p.Fix(varname); // Fix parameter if not varying
    }

  }
}


/*!
 * Deletes the last segment from the segment vector.
 */

void EData::DeleteLastSegment() {
  segments_.pop_back();
}

/*!
 * Fills the Normalizations from the Minuit parameter array.
 */

void EData::FillNormsFromParams(const vector_r &p) {
  int i=GetNormParamOffset();
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    if(segment->IsVaryNorm()) {
      segment->SetNorm(p[i]); 
      i++;
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
}

/*!
 * Fills the Energy Shifts from the Minuit parameter array.
 */

void EData::FillEnergyShiftsFromParams(const vector_r &p, EData *data, CNuc* theCNuc, const Config* configure) {
  int i=GetEnergyShiftParamOffset();
  int k = 0;
  bool anyEnergyShifted = false;
  
  if(data) {
    for(ESegmentIterator segment=data->GetSegments().begin();segment<data->GetSegments().end();segment++) {
      k++;

      // Always apply energy shift since parameters are always present for all segments
      if(segment->IsVaryEnergyShift() || p[i] != 0.0) {

        // Check if energy is the same, if so, continue
        if(segment->GetLastEnergyShift() == p[i]) {
          i++;
          continue;
        }

        segment->SetEnergyShift(p[i]);
        segment->SetLastEnergyShift(p[i]);
        segment->UpdatePointEnergiesWithShift(theCNuc, configure);
        anyEnergyShifted = true;
        
        // Apply the same energy shift to ALL component segments of this master segment
        if(segment->HasComponents()) {
          for(ESegment* componentSegment : segment->GetComponentSegments()) {
            if(componentSegment) {
              componentSegment->SetEnergyShift(segment->GetEnergyShift());
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
            }
          }
        }
      }
      
      // Apply the same energy shift to all component segments in total capture segments
      if(segment->IsTotalCapture() && segment->HasComponents()) {
        if(segment->IsVaryEnergyShift() || p[i] != 0.0) {
          // Apply energy shift to all component segments
          const std::vector<ESegment*>& componentSegments = segment->GetComponentSegments();
          for(ESegment* componentSegment : componentSegments) {
            if(componentSegment) {

              // Check if energy is the same, if so, continue
              if(componentSegment->GetLastEnergyShift() == p[i]) {
                i++;
                continue;
              }

              componentSegment->SetEnergyShift(p[i]);
              componentSegment->SetLastEnergyShift(p[i]);
              componentSegment->UpdatePointEnergiesWithShift(theCNuc, configure);
              anyEnergyShifted = true;
            }
          }
        }
      }
      
      // Always increment i since energy shift parameters are now present for ALL segments
      i++;
    }
    
    // Critical fix: Rebuild energy-based mapping system after energy shifts
    if(anyEnergyShifted) {
      data->MapData();
    }
    
    return;
  }
  for(ESegmentIterator segment=GetSegments().begin();segment<GetSegments().end();segment++) {
    // Always apply energy shift since parameters are always present for all segments
    if(segment->IsVaryEnergyShift() || p[i] != 0.0) {
      segment->SetEnergyShift(p[i]); 
      segment->UpdatePointEnergiesWithShift();
      anyEnergyShifted = true;
      
      // CRITICAL FIX: Apply the same energy shift to ALL component segments of this master segment
      if(segment->HasComponents()) {
        for(ESegment* componentSegment : segment->GetComponentSegments()) {
          if(componentSegment) {
            componentSegment->SetEnergyShift(segment->GetEnergyShift());
            componentSegment->UpdatePointEnergiesWithShift();
          }
        }
      }
    }
    
    // CRITICAL FIX: Apply the same energy shift to all component segments in total capture segments
    if(segment->IsTotalCapture() && segment->HasComponents()) {
      if(segment->IsVaryEnergyShift() || p[i] != 0.0) {
        // Apply energy shift to all component segments
        const std::vector<ESegment*>& componentSegments = segment->GetComponentSegments();
        for(ESegment* componentSegment : componentSegments) {
          if(componentSegment) {
            componentSegment->SetEnergyShift(p[i]);
            componentSegment->UpdatePointEnergiesWithShift();
            anyEnergyShifted = true;
          }
        }
      }
    }
    
    // Always increment i since energy shift parameters are now present for ALL segments
    i++;
  }
  
  // Critical fix: Rebuild energy-based mapping system after energy shifts
  if(anyEnergyShifted) {
    this->MapData();
  }
}

/*!
 * Returns a pointer to a segment specified by a position in the ESegment vector.
 */

ESegment *EData::GetSegment(int segmentNum) {
  ESegment *b=&segments_[segmentNum-1];
  return b;
}

/*!
 * Returns a pointer to a segment based on the segment key, as opposed to a position in the ESegment vector.
 */

ESegment *EData::GetSegmentFromKey(int segmentKey) {
  int segmentNumber=1;
  while(segmentNumber <= this->NumSegments()) {
    if(segmentKey==this->GetSegment(segmentNumber)->GetSegmentKey())
      break;
    else segmentNumber++;
  }
  if(segmentNumber<= this->NumSegments()) return this->GetSegment(segmentNumber);
  else return NULL;
}

/*!
 * Creates a new copy of the EData object in memory and returns a pointer to the new object.
 * Used in AZURECalc function class for thread safety.
 */

EData *EData::Clone() const {
  EData *dataCopy = new EData();

  // Explicitly copy all member variables
  dataCopy->iterations_ = this->iterations_;
  dataCopy->normParamOffset_ = this->normParamOffset_;
  dataCopy->energyShiftParamOffset_ = this->energyShiftParamOffset_;
  dataCopy->isFit_ = this->isFit_;
  dataCopy->isErrorAnalysis_ = this->isErrorAnalysis_;
  dataCopy->targetEffects_ = this->targetEffects_;
  dataCopy->segments_ = this->segments_;
  dataCopy->componentSegments_ = this->componentSegments_;

  // Build a mapping from component segment keys to cloned component segments
  std::unordered_map<int, ESegment*> clonedComponentMap;
  for (size_t i = 0; i < dataCopy->componentSegments_.size(); i++) {
    int segmentKey = dataCopy->componentSegments_[i].GetSegmentKey();
    clonedComponentMap[segmentKey] = &(dataCopy->componentSegments_[i]);
  }
  
  // Update component segment pointers in cloned regular segments
  for (size_t i = 0; i < dataCopy->segments_.size(); i++) {
    auto& segment = dataCopy->segments_[i];
    const auto& originalSegment = this->segments_[i];

    if (originalSegment.HasComponents()) {
      // Clear the old component pointers (they point to original EData)
      segment.ClearComponents();

      // Preserve the operation type
      segment.SetOperationType(originalSegment.GetOperationType());

      // Add the remapped component segments by finding them in the cloned componentSegments_
      for (const auto& originalComponent : originalSegment.GetComponentSegments()) {
        // Find this component in the original componentSegments_ to get its position
        int componentIndex = -1;
        for (size_t j = 0; j < this->componentSegments_.size(); j++) {
          if (&this->componentSegments_[j] == originalComponent) {
            componentIndex = j;
            break;
          }
        }

        // Add the corresponding cloned component
        if (componentIndex >= 0 && componentIndex < (int)dataCopy->componentSegments_.size()) {
          segment.AddComponentSegment(&dataCopy->componentSegments_[componentIndex]);
        }
      }
    }
  }
  
  for(EDataIterator data=dataCopy->begin();data!=dataCopy->end();data++) {
    data.point()->SetParentData(dataCopy);
    data.point()->ClearLocalMappedPoints();
    data.point()->ClearECAmplitudes();
  }
  // Also fix parentData for component segment points
  for (size_t i = 0; i < dataCopy->componentSegments_.size(); i++) {
    for (int j = 1; j <= dataCopy->componentSegments_[i].NumPoints(); j++) {
      dataCopy->componentSegments_[i].GetPoint(j)->SetParentData(dataCopy);
      dataCopy->componentSegments_[i].GetPoint(j)->ClearLocalMappedPoints();
      dataCopy->componentSegments_[i].GetPoint(j)->ClearECAmplitudes();
    }
  }
  for(EDataIterator data=dataCopy->begin();data!=dataCopy->end();data++) {
    if(data.point()->IsMapped()) {
      EnergyMap pointMap = data.point()->GetMap();
      dataCopy->GetSegment(pointMap.segment)->GetPoint(pointMap.point)->AddLocalMappedPoint(&*data.point());
    }
  }

  return dataCopy;
}

/*!
 * Returns a pointer to the specified TargetEffect object.
 */

TargetEffect *EData::GetTargetEffect(int effectNumber) {
  TargetEffect *temp;
  if(effectNumber<=targetEffects_.size()) 
    temp=&targetEffects_[effectNumber-1];
  else return temp=NULL;
  return temp;
}

/*!
 * Returns an EDataIterator referring to the first data point in the set.
 */

EDataIterator EData::begin() {
  return EDataIterator(&segments_);
}

/*!
 * Returns an EDataIterator referring to one object past the last data point in the set.
 */

EDataIterator EData::end() {
  EDataIterator it(&segments_);
  return it.SetEnd();
}

/*!
 * Returns a reference to the vector of ESegment objects.
 */

std::vector<ESegment>& EData::GetSegments() {
  return segments_;
}

