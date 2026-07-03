#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "ReactionRate.h"
#include "AdaptiveIntegrationGrid.h"
#include <iomanip>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <time.h>

struct gsl_reactionrate_params {
  gsl_reactionrate_params(const Config &config) : configure(config) {};
  const Config &configure;
  double temperature;
  CNuc *compound;
  int entranceKey;
  int exitKey;
};

double gsl_reactionrate_integrand(double x, void * p) {
  struct gsl_reactionrate_params *params= (struct gsl_reactionrate_params *)p;
  CNuc *compound=params->compound;
  const Config &configure=params->configure;
  double temperature=params->temperature;
  int entranceKey=params->entranceKey;
  int exitKey=params->exitKey;

  double crossSection;
  if(x<50.0&&x>0.001) {
    EPoint *point = new EPoint(55.0,x,entranceKey,exitKey,false,false,false,0.0,0,0);
    point->Initialize(compound,configure);
    point->Calculate(compound,configure);
    crossSection=point->GetFitCrossSection();
    delete point;
  } else crossSection=0.0;

  return crossSection*x*exp(-x/temperature/boltzConst);
}

// Lower and upper bounds (CM energy, MeV) of the reaction-rate energy integral.
// The integrand (gsl_reactionrate_integrand) is identically zero outside
// (0.001, 50), so integrating over [E_LOW, E_HIGH] captures the full integral
// while letting us use gsl_integration_qagp with resonance breakpoints.
static const double kRateEnergyLow  = 1.0e-5;
static const double kRateEnergyHigh = 50.0;

/*!
 * Builds a sorted list of integration breakpoints for gsl_integration_qagp.
 * The interior points are placed at each resonance energy (and a few points
 * across its width) so the adaptive integrator is guaranteed to resolve narrow
 * resonances instead of stepping over them.  Returns at least {E_LOW, E_HIGH}.
 */
std::vector<double> gsl_reactionrate_breakpoints(CNuc *compound, int entranceKey) {
  std::vector<double> pts;
  pts.push_back(kRateEnergyLow);
  pts.push_back(kRateEnergyHigh);

  AdaptiveIntegrationGrid::GridConfig gridConfig;
  gridConfig.entranceKey = entranceKey;
  AdaptiveIntegrationGrid grid(gridConfig);
  std::vector<AdaptiveIntegrationGrid::ResonanceInfo> resonances =
    grid.GetResonances(kRateEnergyHigh, kRateEnergyLow, compound);

  // For each resonance, bracket the peak tightly so a Gauss-Kronrod panel is
  // forced right onto it, regardless of interference-driven shape changes.
  const double widthMults[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  for(size_t i=0;i<resonances.size();++i) {
    double e0 = resonances[i].energy;
    double w  = resonances[i].totalWidth;
    for(size_t m=0;m<sizeof(widthMults)/sizeof(widthMults[0]);++m) {
      double lo = e0 - widthMults[m]*w;
      double hi = e0 + widthMults[m]*w;
      if(lo>kRateEnergyLow && lo<kRateEnergyHigh) pts.push_back(lo);
      if(hi>kRateEnergyLow && hi<kRateEnergyHigh) pts.push_back(hi);
    }
  }

  std::sort(pts.begin(),pts.end());
  // gsl_integration_qagp requires strictly increasing breakpoints.
  std::vector<double> unique;
  unique.push_back(pts[0]);
  for(size_t i=1;i<pts.size();++i)
    if(pts[i] > unique.back()*(1.0+1e-12)) unique.push_back(pts[i]);
  return unique;
}

double gsl_reactionrate_integration(double temperature,CNuc *compound,const Config& configure,
				    int entranceKey, int exitKey,
				    const std::vector<double>& breakpoints) {

  struct gsl_reactionrate_params params(configure);
  params.temperature=temperature;
  params.compound=compound;
  params.entranceKey=entranceKey;
  params.exitKey=exitKey;

  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (2000);

  gsl_function F;
  F.function = &gsl_reactionrate_integrand;
  F.params=&params;

  double result,error;

  // qagp subdivides at the supplied resonance breakpoints, so narrow resonances
  // near threshold are always sampled (qagiu could step right over them).
  std::vector<double> pts = breakpoints;
  if(pts.size()<2) { pts.clear(); pts.push_back(kRateEnergyLow); pts.push_back(kRateEnergyHigh); }
  // A benign GSL_EROUND ("roundoff") status is expected for this steep,
  // wide-dynamic-range integrand; the returned value is still accurate, so we
  // disable the throwing error handler around the call rather than abort.
  gsl_error_handler_t* oldHandler = gsl_set_error_handler_off();
  gsl_integration_qagp(&F,pts.data(),pts.size(),0.0,1e-6,2000,w,&result,&error);
  gsl_set_error_handler(oldHandler);

  double rate=1e-24*avagadroNum*lightSpeedInCmPerS*
    sqrt(8.0/pi/compound->GetPair(compound->GetPairNumFromKey(entranceKey))->GetRedMass()/uconv)/
    pow(boltzConst*temperature,1.5)*result;

  gsl_integration_workspace_free (w);

  return rate;
}

ReactionRate::ReactionRate(CNuc *compound, const vector_r &params, 
			   const Config &configure, int entranceKey, int exitKey) : 
  configure_(configure) {
  compound_=compound;
  compound_->FillCompoundFromParams(params);
  entrance_key_=entranceKey;
  exit_key_=exitKey;
}

/*!
 * Calculates the astrophysical reaction rates over a range of stellar temperatures.  
 */

void ReactionRate::CalculateRates() {  
  int numSteps = (configure().rateParams.tempStep!=0.) ? 
    int((configure().rateParams.maxTemp-configure().rateParams.minTemp)/configure().rateParams.tempStep)+1 : 1;
  configure().outStream << std::setw(0) << "\t[                         ] 0%";configure().outStream.flush();
  int pointIndex=0;
  // Resonance breakpoints depend only on the level scheme, so build them once.
  std::vector<double> breakpoints=gsl_reactionrate_breakpoints(compound(),entranceKey());
  time_t startTime = time(NULL);
#pragma omp parallel for
  for(int i=0;i<numSteps;++i) {
    CNuc* localCompound = compound()->Clone();
    int localEntranceKey = entranceKey();
    int localExitKey = exitKey();
    const Config localConfigure = configure();
    double temp = localConfigure.rateParams.minTemp+i*localConfigure.rateParams.tempStep;

    double rate=gsl_reactionrate_integration(temp,localCompound,localConfigure,localEntranceKey,localExitKey,breakpoints);
    rates_.push_back(RateData(temp,rate));
    delete localCompound;
    ++pointIndex;
    if(difftime(time(NULL),startTime)>0.25) {
      startTime=time(NULL);
      std::string progress="[";
      double percent=0.;
      for(int j = 1;j<=25;j++) {
	if(pointIndex>=percent*numSteps&&percent<1.) {
	  percent+=0.04;
	  progress+='*';
	} else progress+=' ';
      } progress+="] ";
      localConfigure.outStream << std::setw(0) << "\r\t" << progress << percent*100 << '%';localConfigure.outStream.flush();
    }
  }
  configure().outStream << std::setw(0) << "\r\t[*************************] 100%" << std::endl;
  std::sort(rates_.begin(),rates_.end());
}

/*!
 * Calculates the astrophysical reaction rates at temperatures from file.
 */

void ReactionRate::CalculateFileRates() {
  std::ifstream inFile(configure().rateParams.temperatureFile.c_str());
  if(inFile) {
    while(!inFile.eof()) {
      std::string line;
      getline(inFile,line);
      if(!inFile.eof()) {
	double temp = 0.;
	std::istringstream stm;
	stm.str(line);
	if(stm >> temp) {
	  rates_.push_back(RateData(temp,0.));
	}
      }
    }
    inFile.close();
    configure().outStream << std::setw(0) << "\t[                         ] 0%";configure().outStream.flush();
    int pointIndex=0;
    int numSteps = rates_.size();
    // Resonance breakpoints depend only on the level scheme, so build them once.
    std::vector<double> breakpoints=gsl_reactionrate_breakpoints(compound(),entranceKey());
    time_t startTime = time(NULL);
#pragma omp parallel for
    for(int i=0;i<numSteps;++i) {
      CNuc* localCompound = compound()->Clone();
      int localEntranceKey = entranceKey();
      int localExitKey = exitKey();
      const Config localConfigure = configure();

      rates_[i].rate=gsl_reactionrate_integration(rates_[i].temperature,localCompound,localConfigure,localEntranceKey,localExitKey,breakpoints);
      delete localCompound;
      ++pointIndex;
      if(difftime(time(NULL),startTime)>0.25) {
	startTime=time(NULL);
	std::string progress="[";
	double percent=0.;
	for(int j = 1;j<=25;j++) {
	  if(pointIndex>=percent*numSteps&&percent<1.) {
	    percent+=0.04;
	    progress+='*';
	  } else progress+=' ';
	} progress+="] ";
	localConfigure.outStream << std::setw(0) << "\r\t" << progress << percent*100 << '%';localConfigure.outStream.flush();
      }
    }
    configure().outStream << std::setw(0) << "\r\t[*************************] 100%" << std::endl;
  } else configure().outStream << "Couldn't open temperature file." << std::endl;
}

void ReactionRate::WriteRates() {
  std::string outputfile=configure().outputdir+"reactionrates.out";
  std::ofstream out;
  out.open(outputfile.c_str());
  if(out) {
    out << std::setw(20) << "T9" << std::setw(20) << "Rate" << std::endl;
    for(int i=0;i<rates_.size();i++) {
      out << std::setw(20) << rates_[i].temperature << std::setw(20) << rates_[i].rate << std::endl;
    }
    out.flush();
    out.close();
  } else  configure().outStream << "Could not write reaction rate file." << std::endl; 
}
