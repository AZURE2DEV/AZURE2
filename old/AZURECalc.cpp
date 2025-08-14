#include "AZUREOutput.h"
#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include <iostream>
#include <iomanip>

double AZURECalc::operator()(const vector_r& p) const {
  int thisIteration=data()->Iterations();
  data()->Iterate();
  bool isFit=data()->IsFit();

  CNuc * localCompound = NULL;
  EData *localData = NULL;
  if(isFit) {
    localCompound = compound()->Clone();
    localData = data()->Clone();
	//std::cout << "AzureCalc, Cloned" << std::endl;
  } else {
    localCompound = compound();
    localData = data();
  }

  //Fill Compound Nucleus From Minuit Parameters
  
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
//  std::cout << "Filled local Compound and Data in AZURECalc" << std::endl;
  
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  //loop over segments and points
  double chiSquared=0.0;
  double segmentChiSquared=0.0;

  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  for(EDataIterator data=localData->begin();data!=localData->end();data++) {
    if(data.segment()->GetPoints().begin()==data.point()) {
      segmentChiSquared=0.0;
      if(data.segment()->IsTotalCapture()) {
	firstSumIterator=data.segment();
	lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
      } 
    }
    if(!data.point()->IsMapped()) data.point()->Calculate(localCompound,configure());
    if(firstSumIterator!=localData->GetSegments().end()&&
       data.segment()!=lastSumIterator) continue;
    if((data.segment()->GetSegmentKey()>=1&&data.segment()->GetSegmentKey()<=184)&&
       data.segment()->GetSegmentKey()%2!=0)  continue;
    if(data.segment()->GetSegmentKey()==185) continue;
    if(data.segment()->GetSegmentKey()==187) continue;
    if(data.segment()->GetSegmentKey()==189) continue;
    double fitCrossSection=data.point()->GetFitCrossSection();
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
	fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
      thisSegment = firstSumIterator;
    }
    double dataNorm=thisSegment->GetNorm();
    double CrossSection=data.point()->GetCMCrossSection()*dataNorm;
    double CrossSectionError=data.point()->GetCMCrossSectionError()*dataNorm;
    if((data.segment()->GetSegmentKey()>=1&&data.segment()->GetSegmentKey()<=184)&&
       data.segment()->GetSegmentKey()%2==0) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      EPoint* previousPoint = (data.segment()-1)->GetPoint(pointIndex);
      fitCrossSection/=previousPoint->GetFitCrossSection();
      CrossSection/=previousPoint->GetCMCrossSection();
      CrossSectionError=CrossSection*pow(pow(data.point()->GetCMCrossSectionError()/data.point()->GetCMCrossSection(),2.0)+
					 pow(previousPoint->GetCMCrossSectionError()/previousPoint->GetCMCrossSection(),2.0),0.5);
    }
    if((data.segment()->GetSegmentKey()>=185&&data.segment()->GetSegmentKey()<=190)&&
       data.segment()->GetSegmentKey()%2==0) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      EPoint* previousPoint = (data.segment()-1)->GetPoint(pointIndex);
      fitCrossSection+=previousPoint->GetFitCrossSection();
    } 
    if(data.segment()->GetSegmentKey()>=191&&data.segment()->GetSegmentKey()<=204) {
      fitCrossSection=data.point()->GetFitE1CrossSection();
    }
    if(data.segment()->GetSegmentKey()>=205&&data.segment()->GetSegmentKey()<=215) {
      fitCrossSection=data.point()->GetFitE2CrossSection();
    }    
    double chi=(fitCrossSection-CrossSection)/CrossSectionError;
    double pointChiSquared=pow(chi,2.0);
    //the log value is an alternate method used in Sivia that is suppose to decrease the fit dependance
    //on outlyer data points, which is no doubt an issue. The -1.0 is put in because the Sivia function
    //is to be maximized instead of minimized
    //I'm not sure that multiplying by -1 is what I want. Maybe I should instead take the inverse or 
    //actually change the call to MINUIT to a maximize flag
    if(pointChiSquared > 1.0e-16){  
      pointChiSquared = -1.0*log((1.0-exp(-1.0*pointChiSquared/2.0))/pointChiSquared);
    }

    if(!isFit) data.point()->SetPointChiSquared(pointChiSquared);
    
    segmentChiSquared+=pointChiSquared;
    
    if(data.segment()->GetPoints().end()-1==data.point()) {
    
      if(!isFit) thisSegment->SetSegmentChiSquared(segmentChiSquared);
//      segmentChiSquared=segmentChiSquared/thisSegment->NumPoints();        
      
      
      if(data.segment()==lastSumIterator) {
	firstSumIterator=localData->GetSegments().end();
	lastSumIterator=localData->GetSegments().end();
      }
      double dataNormNominal=thisSegment->GetNominalNorm(); //initial value of normalization 
      double dataNormError=dataNormNominal/100.*thisSegment->GetNormError();
      //I also used this alternate form of the goodness of fit parameter for the data set normalizations
      if(dataNormError!=0.) {
          double normChiSquared = pow((dataNorm-dataNormNominal)/dataNormError,2.0);
          if(!isFit) thisSegment->SetNormChiSquared(normChiSquared);
//          if(normChiSquared > 1.0e-18){	
//	    segmentChiSquared += -1.0*log((1.0-exp(-1.0*normChiSquared/2.0))/normChiSquared);
//	  } else segmentChiSquared += normChiSquared;
          segmentChiSquared += normChiSquared; 
      }
      
      chiSquared+=segmentChiSquared;
    }
    
  }
//  std::cout << "Finished chi2 calculation" << std::endl;

  //get parameter uncertainties initally defined by the user and compare with parameters
  //varied by MINUIT

//  std::cout << "Filled local " << std::endl;
  for (size_t i=0;i<localCompound->GetUserDefinedUncertainties().size();i++){
  
    double dp = localCompound->GetUserDefinedUncertainties().at(i);

    if(dp != 0.0){
      double paramChiSquared = pow((p.at(i)-localCompound->GetNominalParamValues().at(i))/dp,2.0);
      if(paramChiSquared > 1.0e-15){
        paramChiSquared = -1.0*log((1.0-exp(-1.0*paramChiSquared/2.0))/paramChiSquared);     
      }
      chiSquared += paramChiSquared;
    }  
  }
  
//  std::cout << "Finished chi2 comparison of parameter uncertainties" << std::endl;
  
//write output for parameter chi2  
/*  if(!isFit){
    AZUREOutput output(configure().outputdir);
    std::ofstream paramChiOut;
    std::string paramChiOutFile = configure().outputdir+"parameterChiSquared.out";
    paramChiOut.open(paramChiOutFile.c_str());
    paramChiOut << "Fitted Parameter"
                << "\t"
                << "Initial Parameter"
                << "\t"
                << "User Uncertainty"
                << "\t"
	      << "chi2"
                << std::endl;
    for (size_t i=0;i<localCompound->GetUserDefinedUncertainties().size();i++){            
      double dp = fabs(localCompound->GetUserDefinedUncertainties().at(i));
      if(dp != 0.0){
        double paramChiSquared = pow((p.at(i)-localCompound->GetNominalParamValues().at(i))/dp,2.0);
        paramChiOut << p.at(i) << "\t" << localCompound->GetNominalParamValues().at(i) << "\t" << dp << "\t" << paramChiSquared << std::endl;
      }
    }
    paramChiOut.flush();paramChiOut.close();
  }
  
  std::cout << "Finished writting output for parameterChiSquared.out" << std::endl;
*/  
  if(!localData->IsErrorAnalysis()&&thisIteration!=0) {
    if(thisIteration%100==0) configure().outStream
			       << "\r\tIteration: " << std::setw(6) << thisIteration
			       << " Chi-Squared: " << chiSquared << std::endl;  configure().outStream.flush();

    if(thisIteration%1000==0) {
      localData->WriteOutputFiles(configure(),isFit);
      localCompound->TransformOut(configure());
      localCompound->PrintTransformParams(configure());
    }
  }
  
  if(isFit) {
    delete localCompound;
    delete localData;
  }
//  std::cout << "Finished AZURECalc section" << std::cout;
  
  if(configure().stopFlag&&isFit) return 0.;
  else return chiSquared;
}
