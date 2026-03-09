#include "AZUREOutput.h"

/*!
 * The AZUREOutput object is created with reference to an output directory.
 */

AZUREOutput::AZUREOutput(std::string outputdir) {
  outputdir_=outputdir;
  is_extrap_=false;
}

/*! 
 * On destruction of the AZUREOutput instance, each AZUREFBuffer object is also destroyed.
 */

AZUREOutput::~AZUREOutput() {
  for(int i=0;i<azurefbuffers_.size();i++) delete azurefbuffers_[i];
}

/*!
 * Returns true if the output is an extrapolation.
 */

bool AZUREOutput::IsExtrap() const {
  return is_extrap_;
}

/*!
 * The parenthesis operator is defined so that the instance of AZUREOutput can be called as a function.  
 * The instance is called using a reference to an entrance and exit pair key combination.  The function 
 * tests if there is a pointer to an AZUREFBuffer object in a vector.  If such a pointer exists,
 * the pointer to the actual file buffer contained in the corresponding AZUREFBuffer object is returned.
 * Otherwise, a new AZUREFBuffer object is created with the entrance and exit key, a pointer to that object
 * is stored in a vector, and the pointer to the actual new file buffer is returned.
 */

std::filebuf *AZUREOutput::operator()(int entranceKey, int exitKey, bool isAngDist) {
  int c=this->IsAZUREFBuffer(entranceKey,exitKey,isAngDist);
  if(!c) {
    AZUREFBuffer *d = new AZUREFBuffer(entranceKey,exitKey,this->GetOutputDir(),this->IsExtrap(),isAngDist);
    this->AddAZUREFBuffer(d);
    c=this->IsAZUREFBuffer(entranceKey,exitKey,isAngDist);
  }
  return this->GetAZUREFBuffer(c)->GetFBuffer();
}

/*!
 * Returns the number of pointers to AZUREFBuffer objects stored in the vector.
 */

int AZUREOutput::NumAZUREFBuffers() const {
  return azurefbuffers_.size();
}

/*!
 * Tests if a pointer to an AZUREFBuffer object corresponding to the specified entrance
 * and exit keys exists in the vector.  If such a pointer exists, the position in the vector
 * is returned.  Otherwise, the function returns 0.
 */

int AZUREOutput::IsAZUREFBuffer(int entranceKey, int exitKey, bool isAngDist) {
  bool c=false;
  int d=0;
  while(!c&&d<this->NumAZUREFBuffers()) {
    if(this->GetAZUREFBuffer(d+1)->IsAngDist()==isAngDist&&
       this->GetAZUREFBuffer(d+1)->GetEntranceKey()==entranceKey&&
       this->GetAZUREFBuffer(d+1)->GetExitKey()==exitKey) c=true;
    d++;
  }
  if(c) return d;
  else return 0;
}

/*!
 * Returns the output directory for the AZUREOutput object.
 */

std::string AZUREOutput::GetOutputDir() const {
  return outputdir_;
};

/*!
 * Adds a pointer to an AZUREFBuffer object to the vector.
 */

void AZUREOutput::AddAZUREFBuffer(AZUREFBuffer *azureFBuffer) {
  azurefbuffers_.push_back(azureFBuffer);
}

/*!
 * Sets the extrapolation flag to true.
 */

void AZUREOutput::SetExtrap() {
  is_extrap_=true;
}

/*!
 * Returns a pointer to the AZUREFBuffer object specified by a position in the vector.
 */

AZUREFBuffer *AZUREOutput::GetAZUREFBuffer(int fBufferNum) {
  AZUREFBuffer *b=azurefbuffers_[fBufferNum-1];
  return b;
}


