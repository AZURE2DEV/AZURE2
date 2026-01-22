#include "AZURECalc.h"
#include "AZUREMain.h"
#include "AZUREParams.h"
#include "Config.h"
#include "ReactionRate.h"
#include "ParameterLimitsManager.h"
#include "ECAmplitudeCache.h"
#include "Minuit2/MnPrint.h"
#include "GSLException.h"
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <fstream>
#include <sstream>
#include <map>

#ifdef USE_NLOPT
#include <nlopt.h>

// Structure to pass data to the nlopt objective function
struct NLoptData {
  const AZURECalc* theFunc;
  std::vector<double>* full_params;        // Full parameter vector (includes fixed params)
  std::vector<unsigned>* variable_indices; // Indices of variable parameters
};

// NLopt objective function wrapper
// This wrapper receives only the VARIABLE parameters from nlopt,
// reconstructs the full parameter vector with fixed parameters,
// and passes it to AZURECalc
double nlopt_objective(unsigned n, const double* x, double* grad, void* f_data) {
  NLoptData* data = static_cast<NLoptData*>(f_data);

  // Reconstruct full parameter vector by inserting variable params at correct positions
  for(unsigned j = 0; j < n; j++) {
    unsigned i = (*data->variable_indices)[j];
    (*data->full_params)[i] = x[j];
  }

  // Calculate chi-squared with full parameter vector
  double chi_squared = (*data->theFunc)(*data->full_params);

  // Gradient computation not implemented (using derivative-free algorithms)
  if (grad) {
    // NLopt requires gradient to be set to NULL or computed
    // For derivative-free algorithms, this won't be called
  }

  return chi_squared;
}
#endif

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
      
    // Cache will be populated dynamically as integrals are calculated
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
  compound()->FillMnParams(params.GetMinuitParams(), &configure());
  data()->FillMnParams(params.GetMinuitParams());
  if(!(configure().paramMask & Config::USE_PREVIOUS_PARAMETERS)) {
    configure().outStream << "Creating New param.par File..." << std::endl;
    params.WriteUserParameters(configure(),false);
  } else {
    configure().outStream << "Reading User Parameter File..." << std::endl;
    params.ReadUserParameters(configure());
  }

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    // Initialize shared EC amplitude cache for interpolation
    if(configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
      InitializeECAmplitudeCache();
      // Clear cache to ensure fresh start for this calculation
      if(g_ecAmplitudeCache) {
        g_ecAmplitudeCache->Clear();
      }
      configure().outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
    }
    
    //Initialize data object
    if(data()->Initialize(compound(),configure())==-1) {
      configure().outStream << std::endl
			  << "Calculation was aborted." << std::endl;
      return -1;
    }

    // Create ParameterLimitsManager and apply settings
    configure().outStream << "Applying parameter settings with Parameter Manager..." << std::endl;
    ParameterLimitsManager limitsManager(&configure(), compound(), data(), &params);
    limitsManager.ReadParameterSettings();
    limitsManager.ApplyAllParameterSettings(params.GetMinuitParams());
  
    //Declare a new instance of FCNBase
    AZURECalc theFunc(data(),compound(),configure(),&limitsManager);
    //theFunc.InitializePools();
    theFunc.SetErrorDef(1.0);
    
    if(configure().paramMask & Config::PERFORM_FIT) {
      //Call minimizer for function minimization, write minimized parameters to params
      if(configure().paramMask & Config::USE_AMATRIX) configure().outStream << "Performing A-Matrix Fit..." << std::endl;
      else configure().outStream << "Performing R-Matrix Fit..." << std::endl;
      data()->SetFit(true);

#ifdef USE_NLOPT
      if(configure().paramMask & Config::USE_NLOPT_MINIMIZER) {
        // Use NLopt for minimization
        configure().outStream << "Using NLopt minimizer..." << std::endl;

        // Build full parameter vector and identify variable parameters
        unsigned n_total = params.GetMinuitParams().Params().size();
        std::vector<double> full_params(n_total);
        std::vector<unsigned> variable_indices;

        // Extract all parameters and identify which are variable
        for(unsigned i = 0; i < n_total; i++) {
          full_params[i] = params.GetMinuitParams().Value(i);
          if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
            variable_indices.push_back(i);
          }
        }

        // Only pass variable parameters to nlopt
        unsigned n_variable = variable_indices.size();
        std::vector<double> x(n_variable);
        std::vector<double> lb(n_variable);
        std::vector<double> ub(n_variable);

        // Extract initial values and bounds for variable parameters only
        for(unsigned j = 0; j < n_variable; j++) {
          unsigned i = variable_indices[j];
          x[j] = params.GetMinuitParams().Value(i);

          if(params.GetMinuitParams().Parameter(i).HasLowerLimit())
            lb[j] = params.GetMinuitParams().Parameter(i).LowerLimit();
          else
            lb[j] = -1e10; // Large negative number for unbounded

          if(params.GetMinuitParams().Parameter(i).HasUpperLimit())
            ub[j] = params.GetMinuitParams().Parameter(i).UpperLimit();
          else
            ub[j] = 1e10; // Large positive number for unbounded
        }

        configure().outStream << "Optimizing " << n_variable << " parameters (out of "
                              << n_total << " total)..." << std::endl;

        // Map algorithm index to nlopt algorithm enum
        nlopt_algorithm algorithms[] = {
          NLOPT_LN_SBPLX,       // 0: SBPLX (Subplex)
          NLOPT_LN_COBYLA,      // 1: COBYLA
          NLOPT_LN_BOBYQA,      // 2: BOBYQA
          NLOPT_LN_NEWUOA,      // 3: NEWUOA
          NLOPT_LN_PRAXIS,      // 4: PRAXIS
          NLOPT_LN_NELDERMEAD   // 5: Nelder-Mead
        };
        const char* algorithm_names[] = {
          "SBPLX", "COBYLA", "BOBYQA", "NEWUOA", "PRAXIS", "Nelder-Mead"
        };

        int alg_index = configure().nloptAlgorithm;

        configure().outStream << "Using NLopt " << algorithm_names[alg_index]
                              << " algorithm..." << std::endl;

        // Create NLopt optimizer with selected algorithm
        nlopt_opt opt = nlopt_create(algorithms[alg_index], n_variable);

        // Set up objective function with wrapper data
        NLoptData nlopt_data;
        nlopt_data.theFunc = &theFunc;
        nlopt_data.full_params = &full_params;
        nlopt_data.variable_indices = &variable_indices;
        nlopt_set_min_objective(opt, nlopt_objective, &nlopt_data);

        // Set bounds
        nlopt_set_lower_bounds(opt, lb.data());
        nlopt_set_upper_bounds(opt, ub.data());

        // Set stopping criteria
        nlopt_set_maxeval(opt, 50000);
        nlopt_set_ftol_rel(opt, 1e-8);
        nlopt_set_xtol_rel(opt, 1e-8);

        // Perform optimization
        double minf;
        nlopt_result result = nlopt_optimize(opt, x.data(), &minf);

        if(result < 0) {
          configure().outStream << std::endl;
          configure().outStream << "NLopt optimization failed: "
                                << nlopt_result_to_string(result) << std::endl;
        } else {
          configure().outStream << std::endl;
          configure().outStream << "NLopt optimization succeeded: "
                                << nlopt_result_to_string(result) << std::endl;
          configure().outStream << "Final chi-squared: " << minf << std::endl;
        }

        // Update only the variable parameters with optimized values
        // Fixed parameters remain unchanged in full_params
        for(unsigned j = 0; j < n_variable; j++) {
          unsigned i = variable_indices[j];
          params.GetMinuitParams().SetValue(i, x[j]);
        }

        nlopt_destroy(opt);

        // Note: Error analysis with NLopt is not directly available
        // Skip Minos error analysis for NLopt
        if(configure().paramMask & Config::PERFORM_ERROR_ANALYSIS) {
          configure().outStream << std::endl
                    << "Error analysis is not available with NLopt minimizer." << std::endl;
          configure().outStream << "Please use Minuit2 for error analysis." << std::endl;
        }

        params.WriteUserParameters(configure(),true);
      } else {
#endif
        // Use Minuit2 for minimization (default)
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
          std::vector<int> parameterTable;
          parameterTable.reserve(params.GetMinuitParams().Params().size());

          // Header for Covariance Matrix
          out << std::endl <<  "Covariance Matrix" << std::endl << std::endl;
          out << std::setw(5) << " ";
          for(int i = 0; i<params.GetMinuitParams().Params().size(); i++) {
            if(!params.GetMinuitParams().Parameter(i).IsFixed()) {
              out << std::setw(15) << i ;
              parameterTable.push_back(i);
            }
          }
          out << std::endl;
          int size = parameterTable.size();

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
#ifdef USE_NLOPT
      }
#endif
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
	// With component segments, no need to skip - total capture is handled within the segment
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

  // Cleanup shared EC amplitude cache
  if(configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
    if(g_ecAmplitudeCache) {
      //g_ecAmplitudeCache->PrintStats();
    }
    CleanupECAmplitudeCache();
  }

  return 0;
}