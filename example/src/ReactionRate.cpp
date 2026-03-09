#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "ReactionRate.h"
#include <iomanip>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include <algorithm>
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

double gsl_reactionrate_integration(double temperature,CNuc *compound,const Config& configure,
				    int entranceKey, int exitKey) {

  struct gsl_reactionrate_params params(configure);
  params.temperature=temperature;
  params.compound=compound;
  params.entranceKey=entranceKey;
  params.exitKey=exitKey;
  
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &gsl_reactionrate_integrand;
  F.params=&params;

  double result,error;
  
  gsl_integration_qagiu(&F,0.00001,0.0,1e-4,1000,w,&result,&error);

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
  time_t startTime = time(NULL);
#pragma omp parallel for
  for(int i=0;i<numSteps;++i) {
    CNuc* localCompound = compound()->Clone();
    int localEntranceKey = entranceKey();
    int localExitKey = exitKey();
    const Config localConfigure = configure(); 
    double temp = localConfigure.rateParams.minTemp+i*localConfigure.rateParams.tempStep;

    double rate=gsl_reactionrate_integration(temp,localCompound,localConfigure,localEntranceKey,localExitKey);
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
    time_t startTime = time(NULL);
#pragma omp parallel for
    for(int i=0;i<numSteps;++i) {
      CNuc* localCompound = compound()->Clone();
      int localEntranceKey = entranceKey();
      int localExitKey = exitKey();
      const Config localConfigure = configure(); 
      
      rates_[i].rate=gsl_reactionrate_integration(rates_[i].temperature,localCompound,localConfigure,localEntranceKey,localExitKey);
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
