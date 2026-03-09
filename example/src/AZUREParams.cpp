#include "AZUREParams.h"
#include "Config.h"

/*!
 * This function returns the MnUserParameters object used by Minuit to store 
 * the fit parameters.
 */

ROOT::Minuit2::MnUserParameters &AZUREParams::GetMinuitParams() {
  return params_;
}

/*!
 * This function reads the user specified parameters from a given file.  
 * These parameters are formal R-matrix parameters, and overwrite any
 * initial parameters determined from the nuclear input file.
 */

void AZUREParams::ReadUserParameters(const Config& configure) {
  std::vector<std::string> names;
  vector_r values;
  vector_r errors;
  std::vector<bool> fixed;
  std::string tempname,tempfixed,tempfixed_nows;
  double tempvalue,temperror;

  std::ifstream in;
  in.open(configure.paramfile.c_str());
  if(in) {
    while(!in.eof()) {
      in >> tempname >> tempvalue >> temperror; getline(in,tempfixed);
      if(!in.eof()) {
	names.push_back(tempname);
	values.push_back(tempvalue);
	errors.push_back(temperror);
	tempfixed_nows.clear();
	for(int i=0;i<tempfixed.length();i++) 
	  if(tempfixed[i]!=' '&&tempfixed[i]!='\t') 
	    tempfixed_nows.push_back(tempfixed[i]);
	if(tempfixed_nows=="fixed") fixed.push_back(true);
	else fixed.push_back(false);
      }
    }
    in.close();
  } else configure.outStream << "Could not read user parameter file." << std::endl;
  
  for(int i=0;i<GetMinuitParams().Params().size();i++) {
    for(int ii=0;ii<names.size();ii++) {
      if(GetMinuitParams().GetName(i)==names[ii]) {
	if(GetMinuitParams().Value(i)==0.0&&values[ii]!=0.0&&
	   GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Release(i);
	if(GetMinuitParams().Value(i)!=0.0&&values[ii]==0.0&&
	   !GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Fix(i);
	GetMinuitParams().SetValue(i,values[ii]);
	GetMinuitParams().SetError(i,errors[ii]);
	if(fixed[ii]&&!GetMinuitParams().Parameter(i).IsFixed()) GetMinuitParams().Fix(i);
      }
    }
  }
}

/*!
 * This function writes the formal R-matrix parameters to a file.
 */

void AZUREParams::WriteUserParameters(const Config& configure, bool fitParameters) {
  char filename[256];
  if(fitParameters) sprintf(filename,"%sparam.sav",configure.outputdir.c_str());
  else sprintf(filename,"%sparam.par",configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if(out) {
    out.precision(7);
    for(int i=0;i<GetMinuitParams().Params().size();i++) {
      out << std::setw(20) << GetMinuitParams().GetName(i) 
	  << std::scientific << std::setw(20) <<  GetMinuitParams().Value(i) 
	  << std::scientific << std::setw(20) <<  GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else configure.outStream << "Could not save param.par file." << std::endl;
}

/*!
 * This function writes the parameter errors to a file if Minos has been invoked. 
 */

void AZUREParams::WriteParameterErrors(const std::vector<std::pair<double,double> > &errors,const Config& configure) {
  char filename[256];
  sprintf(filename,"%sparam.errors",configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if(out) {
    out.precision(7);
    for(int i=0;i<GetMinuitParams().Params().size();i++) {
      out << std::setw(20) << GetMinuitParams().GetName(i) 
	  << std::scientific << std::setw(20) <<  GetMinuitParams().Value(i) 
	  << std::scientific << std::setw(20) <<  fabs(errors[i].first) 
	  << std::scientific << std::setw(20) <<  fabs(errors[i].second) 
	  << std::endl;
    }
    out.flush();
    out.close();
  } else configure.outStream << "Could not save param.errors file." << std::endl;
}
