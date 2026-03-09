#include "AZURECalc.h"
#include "AZUREMain.h"
#include "AZUREParams.h"
#include "Config.h"
#include "ReactionRate.h"
#include "Minuit2/MnPrint.h"
#include "GSLException.h"
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>

int AZUREMain::operator()(){
  //Fill compound nucleus from nucfile
  configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if(compound()->Fill(configure())==-1) {
    configure().outStream << "Could not fill compound nucleus from file." 
	      << std::endl;
    return -1;
  } else if(compound()->NumPairs()==0 || compound()->NumJGroups()==0) {
    configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl; 
    return -1;
  } 
  if((configure().screenCheckMask|configure().fileCheckMask) & 
     Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    //Fill the data object from the segments and data file
    //  Compound object is passed to the function for pair key verification and
    //  center of mass conversions, s-factor conversions, etc.
    configure().outStream << "Filling Data Structures..." << std::endl;
    if(configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if(data()->Fill(configure(),compound())==-1) {
	configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	configure().outStream << "There is no data provided." << std::endl;
	return -1;
      }
    } else {
      if(data()->MakePoints(configure(),compound())==-1) {
	configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	configure().outStream << "Extrapolation segments produce no data." << std::endl;
	return -1;
      }
    } 
    if((configure().fileCheckMask|configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());
  } else {
    if(!compound()->IsPairKey(configure().rateParams.entrancePair)||!compound()->IsPairKey(configure().rateParams.exitPair)) {
      configure().outStream << "Reaction rate pairs do not exist in compound nucleus." << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  //Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << std::endl
			  << "Calculation was aborted." << std::endl;
    return -1;
  }

  //Create new parameters for minuit, fill them from compound nucleus object and data file.
  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams());
  data()->FillMnParams(params.GetMinuitParams());
  if(!(configure().paramMask & Config::USE_PREVIOUS_PARAMETERS)) {
    configure().outStream << "Creating New param.par File..." << std::endl;
    params.WriteUserParameters(configure(),false);
  } else {
    configure().outStream << "Reading User Parameter File..." << std::endl;
    params.ReadUserParameters(configure());
  }

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    //Initialize data object
    if(data()->Initialize(compound(),configure())==-1) {
      configure().outStream << std::endl
			  << "Calculation was aborted." << std::endl;
      return -1;
    }
 
    //Declare a new instance of FCNBase
    AZURECalc theFunc(data(),compound(),configure());
    theFunc.SetErrorDef(1.0);
    
    if(configure().paramMask & Config::PERFORM_FIT) {
      //Call Minuit for function minimization, write minimized parameters to params
      if(configure().paramMask & Config::USE_AMATRIX) configure().outStream << "Performing A-Matrix Fit..." << std::endl; 
      else configure().outStream << "Performing R-Matrix Fit..." << std::endl;
      data()->SetFit(true);
      ROOT::Minuit2::MnMigrad migrad(theFunc,params.GetMinuitParams());
      ROOT::Minuit2::FunctionMinimum min=migrad(50000);
      if(configure().paramMask & Config::PERFORM_ERROR_ANALYSIS) {
	configure().outStream << std::endl 
		  << "Performing parameter error analysis with Up=" <<  configure().chiVariance << "." << std::endl;
	data()->SetErrorAnalysis(true);
	theFunc.SetErrorDef(configure().chiVariance);
	ROOT::Minuit2::MnMinos minos(theFunc,min);
	std::vector<std::pair<double,double> > errors;
	for(int i = 0; i<params.GetMinuitParams().Params().size(); i++) { 
	  configure().outStream << "\tParameter " << i+1 << "..." << std::endl;
	  if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
	    std::pair< double, double > error=minos(i);
	    errors.push_back(error);
	  } else errors.push_back(std::pair< double, double > (0.,0.));
	}

        // New output of the covariance matrix
        char filename[256];
        sprintf(filename,"%scovariance_matrix.out",configure().outputdir.c_str());
        std::ofstream out;
        out.open(filename);
        if(out) {

          // Write header information
          params.GetMinuitParams()=min.UserParameters();
          out << "Parameter List" << std::endl << std::endl;
          for(int i = 0; i<params.GetMinuitParams().Params().size(); i++) { 
            out << std::setw(20) << params.GetMinuitParams().GetName(i)
                << " (" << std::setw(3) << i << ")"	 
                << std::scientific << std::setw(20) <<  params.GetMinuitParams().Value(i);
            if(params.GetMinuitParams().Parameter(i).IsFixed()) {
	      out << " Fixed" << std::endl;
            }
	    else {
              out << " Fitted by Minuit" << std::endl;    
            }
	  }

	  std::cout << min.UserCovariance();

	  std::vector <double> CovarianceData;
          CovarianceData=min.UserCovariance().Data();
          int parameterTable[100];
          int size=0;

          // Header for Covariance Matrix
          out << std::endl <<  "Covariance Matrix" << std::endl << std::endl;
          out << std::setw(5) << " ";
          for(int i = 0; i<params.GetMinuitParams().Params().size(); i++) {
            if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
              out << std::setw(15) << i ;
              parameterTable[size]=i;
              size=size+1;
            }
          }
          out << std::endl;

	  std::cout << "covariance length " << CovarianceData.size() << std::endl;

          // Covariance matrix
          for(int i = 0; i<size; i++) {
            out << std::setw(5) << parameterTable[i];
            for(int j = 0; j<size; j++) {
              int k;
              if(i<=j) { k=((j+1)*j)/2+i; }
              if(i>j)  { k=((i+1)*i)/2+j; }
	      std::cout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k << std::endl;
              out << std::setw(15) << CovarianceData[k];
            }
            out << std::endl;
	  }

          // Header for Correlation Matrix
          out << std::endl <<  "Correlation Matrix" << std::endl << std::endl;
          out << std::setw(5) << " ";
          for(int i = 0; i<params.GetMinuitParams().Params().size(); i++) {
            if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
              out << std::setw(15) << i ;
            }
          }
          out << std::endl;

          // Correlation matrix
          for(int i = 0; i<size; i++) {
            out << std::setw(5) << parameterTable[i];
            for(int j = 0; j<size; j++) {
              int k;
              if(i<=j) { k=((j+1)*j)/2+i; }
              if(i>j)  { k=((i+1)*i)/2+j; }
              int jdiag=((j+2)*(j+1))/2-1;
              int idiag=((i+2)*(i+1))/2-1;
	      std::cout << k << "  " << CovarianceData[k] << std::endl;
              out << std::setw(15) << std::fixed << CovarianceData[k]/sqrt(fabs(CovarianceData[jdiag]*CovarianceData[idiag]));
            }
            out << std::endl;
	  }

          out.flush();
          out.close();
        } else std::cout << "Could not save " << configure().outputdir.c_str() 
                         << "covariance_matrix.out file." << std::endl;   

        // ECS this line added to set the azure variables to the minimised fit parameters
        params.GetMinuitParams()=min.UserParameters();
	params.WriteParameterErrors(errors,configure());

      }
      params.GetMinuitParams()=min.UserParameters();
      params.WriteUserParameters(configure(),true);
    } else {
      if(configure().paramMask & Config::USE_AMATRIX) configure().outStream << "Performing A-Matrix Calculation..." << std::endl; 
      else configure().outStream << "Performing R-Matrix Calculation..." << std::endl; 
    }
    data()->SetFit(false);
    data()->SetErrorAnalysis(false);
    double chiSquared=theFunc(params.GetMinuitParams().Params());
    if(configure().paramMask & Config::CALCULATE_WITH_DATA) {
      configure().outStream << std::endl << std::endl;
      for(ESegmentIterator segment=data()->GetSegments().begin();
	  segment<data()->GetSegments().end();segment++) {
	configure().outStream << "Segment #"
		  << segment->GetSegmentKey() 
		  << " Chi-Squared/N: "
		  << segment->GetSegmentChiSquared()/segment->NumPoints()
		  << std::endl;
	if(segment->IsTotalCapture()) segment+=segment->IsTotalCapture()-1;
      }
      configure().outStream << "Total Chi-Squared: " 
		<< chiSquared << std::endl << std::endl;
    }

    //Write Output Files
    configure().outStream << "Writing output files..." << std::endl;
    data()->WriteOutputFiles(configure());
  } else {
    //Calculate Reaction Rate
    // This uses the adaptive integration routines of GSL.  As the energy stepsize is 
    // unknown until integration, AZURE is called not in segment control mode but as 
    // an energy dependent function.  As every data point must be reinitialized (new energy 
    // dependent terms calculated) the routine is slow.  This should be a tradeoff 
    // for good accuracy.
    configure().outStream << "Performing reaction rate calculation..." << std::endl;
    ReactionRate reactionRate(compound(),params.GetMinuitParams().Params(),configure(),
			      configure().rateParams.entrancePair,configure().rateParams.exitPair);
    if(configure().paramMask & Config::USE_BRUNE_FORMALISM) 
      compound()->CalcShiftFunctions(configure());
    if(configure().rateParams.useFile)
      reactionRate.CalculateFileRates();
    else 
      reactionRate.CalculateRates();
    reactionRate.WriteRates();
  }

  configure().outStream << "Performing final parameter transformation..." << std::endl;
  try {
    compound()->TransformOut(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << "Problem with output transformation.  Aborting." 
			  << std::endl;
    return -1;
  }
  compound()->PrintTransformParams(configure());

  return 0;
}
