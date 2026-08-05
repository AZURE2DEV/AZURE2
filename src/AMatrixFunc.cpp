#include "AMatrixFunc.h"
#include "PolarizationFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "AZUREGrad.h"
#include "Decay.h"
#include "KGroup.h"
#include "KLGroup.h"
#include "MGroup.h"
#include "Interference.h"
#include "JGroup.h"
#include "ALevel.h"
#include "AChannel.h"
#include "PPair.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * The AMatrixFunc object is created with reference to a CNuc object.
 */

AMatrixFunc::AMatrixFunc(CNuc* compound, const Config &configure) :
  configure_(&configure), compound_(compound),
  cached_max_levels_(0), cached_max_channels_(0) {}

/*!
 * Returns the A-Matrix of a J-Group, materializing it from the stored
 * factorization on first use.  Returns an empty matrix if the J-Group has no
 * usable factorization.
 */

const matrix_c &AMatrixFunc::GetAMatrix(int jGroupNum) const {
  static const matrix_c emptyMatrix;
  if(jGroupNum<1||jGroupNum>(int)solvers_.size()) return emptyMatrix;
  if(!solver_valid_[jGroupNum-1]) return emptyMatrix;
  return solvers_[jGroupNum-1].Inverse();
}

/*!
 * Returns an A-Matrix element specified by positions in the JGroup and ALevel vectors.
 * The level numbers are the original ones; they are mapped onto the compacted
 * indices of the levels actually included in the R-Matrix.
 */

complex AMatrixFunc::GetAMatrixElement(int jGroupNum, int lambdaNum, int muNum) const {
  if(jGroupNum<1||jGroupNum>(int)solvers_.size()) return complex(0.0,0.0);
  if(!solver_valid_[jGroupNum-1]) return complex(0.0,0.0);
  const std::vector<int> &map = level_active_index_[jGroupNum-1];
  if(lambdaNum<1||lambdaNum>=(int)map.size()) return complex(0.0,0.0);
  if(muNum<1||muNum>=(int)map.size()) return complex(0.0,0.0);
  const int row=map[lambdaNum];
  const int col=map[muNum];
  if(row==0||col==0) return complex(0.0,0.0);
  return solvers_[jGroupNum-1].Inverse()[row-1][col-1];
}

/*!
 * Returns a pointer to an entire A-Matrix specified by a position in the JGroup vector.
 */

matrix_c *AMatrixFunc::GetJSpecAInvMatrix(int jGroupNum) {
  matrix_c *b=&a_inv_matrices_[jGroupNum-1];
  return b;
}

/*!
 * Clears all matrices associated with the AMatrixFunc object.
 */

void AMatrixFunc::ClearMatrices() {
  tmatrix_.clear();
  ec_tmatrix_.clear();
  chan_cap_cache_.clear();

  const int numJGroups = compound()->NumJGroups();

  // Grow the per-JGroup storage rather than rebuilding it, so an instance reused
  // across energy points reallocates nothing after its first point.
  if ((int)a_inv_matrices_.size() != numJGroups) {
    a_inv_matrices_.resize(numJGroups);
    solvers_.resize(numJGroups);
    solver_valid_.resize(numJGroups);
    level_active_index_.resize(numJGroups);
    gamma_vectors_.resize(numJGroups);
    gamma_valid_.resize(numJGroups);
    channel_solves_.resize(numJGroups);
    channel_solve_valid_.resize(numJGroups);
  }

  for (int j = 1; j <= numJGroups; ++j) {
    solver_valid_[j-1] = 0;
    if (!compound()->GetJGroup(j)->IsInRMatrix()) {
      // keep empty placeholders to preserve external indexing
      level_active_index_[j-1].clear();
      a_inv_matrices_[j-1].clear();
      gamma_valid_[j-1].clear();
      channel_solve_valid_[j-1].clear();
      continue;
    }

    JGroup* jg = compound()->GetJGroup(j);
    const int numLevels = jg->NumLevels();
    const int numChannels = jg->NumChannels();

    // Build original->active map (1-based for readability in existing code)
    std::vector<int>& map = level_active_index_[j-1];
    map.assign(numLevels + 1, 0);
    int act = 0;
    for (int la = 1; la <= numLevels; ++la) {
      if (jg->GetLevel(la)->IsInRMatrix()) map[la] = ++act;
    }

    // Pre-size to dense (act x act) with zeros; no push_backs later
    matrix_c& M = a_inv_matrices_[j-1];
    M.resize(act);
    for (int r = 0; r < act; ++r) M[r].assign(act, complex(0.0, 0.0));

    gamma_vectors_[j-1].resize(numChannels);
    gamma_valid_[j-1].assign(numChannels, 0);
    channel_solves_[j-1].resize(numChannels);
    channel_solve_valid_[j-1].assign(numChannels, 0);
  }
}


/*!
 * This function creates the inverted A-Matrix from the parameters in the CNuc object.
 */

void AMatrixFunc::FillMatrices (EPoint *point) {
  double inEnergy;
  if(compound()->
     GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->
     GetPType()==20)
    inEnergy=point->GetCMEnergy()+
      compound()->
      GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->
      GetSepE()+
      compound()->
      GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->
      GetExE();
  else inEnergy = point->GetCMEnergy()+
	 compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetSepE()+
	 compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetExE();
  for(int j=1;j<=compound()->NumJGroups();j++) {
    if(compound()->GetJGroup(j)->IsInRMatrix()) {
      // Cache JGroup pointer to avoid repeated calls
      JGroup *jGroup = compound()->GetJGroup(j);
      int numLevels = jGroup->NumLevels();
      int numChannels = jGroup->NumChannels();
      
      // Ensure pre-allocated buffers are large enough - resize only when needed
      if(numLevels > cached_max_levels_ || numChannels > cached_max_channels_) {
        cached_max_levels_ = std::max(numLevels, cached_max_levels_);
        cached_max_channels_ = std::max(numChannels, cached_max_channels_);
        
        levelGammas_.resize(cached_max_levels_ + 1);
        levelEnergies_.resize(cached_max_levels_ + 1);
        shiftFunctions_.resize(cached_max_levels_ + 1);
        
        for(int la = 0; la <= cached_max_levels_; la++) {
          levelGammas_[la].resize(cached_max_channels_ + 1);
          shiftFunctions_[la].resize(cached_max_channels_ + 1);
        }
      }
      
      // Pre-cache gamma values for all levels and channels using pre-allocated buffers
      for(int la=1; la<=numLevels; la++) {
        if(jGroup->GetLevel(la)->IsInRMatrix()) {
          ALevel *level = jGroup->GetLevel(la);
          levelEnergies_[la] = level->GetFitE();
          
          for(int ch=1; ch<=numChannels; ch++) {
            levelGammas_[la][ch] = level->GetFitGamma(ch);
            shiftFunctions_[la][ch] = level->GetShiftFunction(ch);
          }
        }
      }
      
      for(int la=1;la<=numLevels;la++) {
	if(jGroup->GetLevel(la)->IsInRMatrix()) {
	  ALevel *level=jGroup->GetLevel(la);
	  for(int lap=1;lap<=numLevels;lap++) {
	    if(jGroup->GetLevel(lap)->IsInRMatrix()) {
	      ALevel *levelp=jGroup->GetLevel(lap);
	      complex sum(0.0,0.0);
	      for(int ch=1;ch<=numChannels;ch++) {
		double gammaCh=levelGammas_[la][ch];
		double gammaChp=levelGammas_[lap][ch];
		
		// Early termination for effectively zero gamma values
		if(fabs(gammaCh) < 1.0e-12 || fabs(gammaChp) < 1.0e-12) continue;
		
		complex loElement=point->GetLoElement(j,ch);
		sum+=gammaCh*gammaChp*loElement;
		
		// Cache channel pointer and radiation type to avoid repeated calls
		AChannel *channel = jGroup->GetChannel(ch);
		char radType = channel->GetRadType();
		
		if((radType == 'M' || radType == 'E') && 
		   la==lap &&
		   (configure().paramMask & Config::USE_RMC_FORMALISM)) 
		  sum+=complex(0.0,1.0)*gammaCh*gammaChp;
		if((configure().paramMask & Config::USE_BRUNE_FORMALISM) && radType=='P') {
		  sum+=gammaCh*gammaChp*channel->GetBoundaryCondition();
		  if(la==lap) sum-=gammaCh*gammaChp*shiftFunctions_[la][ch];
		  else sum-=gammaCh*gammaChp*
			 (shiftFunctions_[la][ch]*(inEnergy-levelEnergies_[lap])-shiftFunctions_[lap][ch]*(inEnergy-levelEnergies_[la]))/
			 (levelEnergies_[la]-levelEnergies_[lap]);				
		}
	      }
	      if(la==lap) {
		double resenergy=levelEnergies_[la];
		this->AddAInvMatrixElement(j,la,lap,resenergy-inEnergy-sum);
	      } else this->AddAInvMatrixElement(j,la,lap,-sum);
	    }
	  }
	}     
      }
    }
  }
}

/*!
 * This function factors the inverse A-Matrix of each J-Group.  The A-Matrix
 * itself is not formed: CalculateTMatrix needs only bilinear forms of it, which
 * the factorization supplies with one triangular solve per channel, and the
 * callers that do need individual elements get them from GetAMatrixElement,
 * which materializes the inverse on demand.
 */

void AMatrixFunc::InvertMatrices() {
  for(int j=1;j<=compound()->NumJGroups();j++) {
    if(compound()->GetJGroup(j)->IsInRMatrix()) {
      matrix_c *theAInvMatrix = this->GetJSpecAInvMatrix(j);
      // Add validation to catch corrupted matrices before processing
      if(theAInvMatrix->empty()) {
        continue; // Skip empty matrices
      }
      solver_valid_[j-1] = solvers_[j-1].Decompose(*theAInvMatrix) ? 1 : 0;
      // Anything derived from the previous factorization is now stale.
      gamma_valid_[j-1].assign(gamma_valid_[j-1].size(),0);
      channel_solve_valid_[j-1].assign(channel_solve_valid_[j-1].size(),0);
    }
  }
  chan_cap_cache_.clear();
}

/*!
 * Returns the reduced widths of one channel over the active levels of a J-Group.
 * Widths below 1e-12 are zeroed so that the bilinear form reproduces, term for
 * term, the level pairs the explicit double loop used to skip.
 */

const vector_r &AMatrixFunc::GetGammaVector(int jGroupNum, int chNum) {
  vector_r &gammas = gamma_vectors_[jGroupNum-1][chNum-1];
  if(!gamma_valid_[jGroupNum-1][chNum-1]) {
    JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum-1];
    gammas.assign(solvers_[jGroupNum-1].size(),0.0);
    for(int la=1;la<(int)map.size();la++) {
      const int row=map[la];
      if(row==0) continue;
      double gamma=jGroup->GetLevel(la)->GetFitGamma(chNum);
      if(fabs(gamma)<1.0e-12) continue;
      gammas[row-1]=gamma;
    }
    gamma_valid_[jGroupNum-1][chNum-1]=1;
  }
  return gammas;
}

/*!
 * Returns A times the reduced-width vector of one channel, cached for the
 * duration of the energy point.  Every pathway sharing an exit channel reuses
 * the same solve, which is what replaces the old per-pathway double loop.
 */

const vector_c &AMatrixFunc::GetChannelSolve(int jGroupNum, int chNum) {
  vector_c &solve = channel_solves_[jGroupNum-1][chNum-1];
  if(!channel_solve_valid_[jGroupNum-1][chNum-1]) {
    const vector_r &gammas = this->GetGammaVector(jGroupNum,chNum);
    const int n = (int)gammas.size();
    solve.resize(n);
    for(int i=0;i<n;i++) solve[i]=gammas[i];
    if(n>0) solvers_[jGroupNum-1].Solve(&solve[0]);
    channel_solve_valid_[jGroupNum-1][chNum-1]=1;
  }
  return solve;
}

/*!
 * The level matrix is complex symmetric, so A is too and the bilinear form may
 * contract either index with the cached solve.
 */

complex AMatrixFunc::GetUBilinear(int jGroupNum, int chNum, int chpNum) {
  if(jGroupNum<1||jGroupNum>(int)solvers_.size()) return complex(0.0,0.0);
  if(!solver_valid_[jGroupNum-1]) return complex(0.0,0.0);
  // Fill the solve cache first: it populates the gamma cache for chpNum, so
  // taking the chNum reference afterwards keeps it clear that nothing can move
  // out from under it.
  const vector_c &solve = this->GetChannelSolve(jGroupNum,chpNum);
  const vector_r &gammas = this->GetGammaVector(jGroupNum,chNum);
  complex sum(0.0,0.0);
  for(size_t i=0;i<solve.size();i++) sum+=gammas[i]*solve[i];
  return sum;
}

/*!
 * Channel capture excludes levels with a negligible width in the internal
 * channel from both indices, which is a mask on the level set rather than the
 * per-channel cutoff, so this path builds its vectors fresh.  Results are
 * memoized because the same pathway recurs across k-groups.
 */

complex AMatrixFunc::GetChannelCaptureBilinear(int jGroupNum, int chNum, int chpNum,
                                               int maskChannel) {
  if(jGroupNum<1||jGroupNum>(int)solvers_.size()) return complex(0.0,0.0);
  for(size_t i=0;i<chan_cap_cache_.size();i++) {
    const ChanCapEntry &entry = chan_cap_cache_[i];
    if(entry.jGroupNum==jGroupNum&&entry.chNum==chNum&&
       entry.chpNum==chpNum&&entry.maskChannel==maskChannel) return entry.value;
  }
  complex value(0.0,0.0);
  if(solver_valid_[jGroupNum-1]) {
    JGroup *jGroup = compound()->GetJGroup(jGroupNum);
    const std::vector<int> &map = level_active_index_[jGroupNum-1];
    const int n = solvers_[jGroupNum-1].size();
    chan_cap_entrance_.assign(n,complex(0.0,0.0));
    chan_cap_exit_.assign(n,complex(0.0,0.0));
    for(int la=1;la<(int)map.size();la++) {
      const int row=map[la];
      if(row==0) continue;
      ALevel *level=jGroup->GetLevel(la);
      if(maskChannel&&fabs(level->GetFitGamma(maskChannel))<1.0e-8) continue;
      chan_cap_entrance_[row-1]=level->GetFitGamma(chNum);
      chan_cap_exit_[row-1]=level->GetFitGamma(chpNum);
    }
    value=solvers_[jGroupNum-1].Bilinear(chan_cap_entrance_,chan_cap_exit_);
  }
  ChanCapEntry entry;
  entry.jGroupNum=jGroupNum;
  entry.chNum=chNum;
  entry.chpNum=chpNum;
  entry.maskChannel=maskChannel;
  entry.value=value;
  chan_cap_cache_.push_back(entry);
  return value;
}

/*!
 * This function calculates the T-Matrix for each reaction pathway based on the A-Matrix.
 */

void AMatrixFunc::CalculateTMatrix(EPoint *point) {
  // Cache frequently accessed values to avoid repeated function calls
  int entranceKey = point->GetEntranceKey();
  int exitKey = point->GetExitKey();
  int aa = compound()->GetPairNumFromKey(entranceKey);
  int exitPairNum = compound()->GetPairNumFromKey(exitKey);
  
  int irEnd;
  int irStart;
  bool isRMC=false;
  if((configure().paramMask & Config::USE_RMC_FORMALISM) && 
     compound()->GetPair(exitPairNum)->GetPType()==10) {
    irStart=1;
    irEnd=compound()->GetPair(aa)->NumDecays();
    isRMC=true;
  } else {
    irStart=0;
    int numDecays = compound()->GetPair(aa)->NumDecays();
    while(irStart<numDecays) {
      irStart++;
      if(compound()->GetPair(aa)->GetDecay(irStart)->GetPairNum()==exitPairNum) break;
    }
    irEnd=irStart;
  }
  for(int ir=irStart;ir<=irEnd;ir++) {
    Decay *theDecay=compound()->GetPair(aa)->GetDecay(ir);
    for(int k=1;k<=theDecay->NumKGroups();k++) {
      for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
	MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	// Cache frequently accessed values to avoid repeated function calls
	int jNum = theMGroup->GetJNum();
	int chNum = theMGroup->GetChNum();
	int chpNum = theMGroup->GetChpNum();
	
	JGroup *theJGroup=compound()->GetJGroup(jNum);
	AChannel *entranceChannel=theJGroup->GetChannel(chNum);
	AChannel *exitChannel=theJGroup->GetChannel(chpNum);
	
	// Cache phase calculations 
	complex coulombPhaseEn = point->GetExpCoulombPhase(jNum,chNum);
	complex hardSpherePhaseEn = point->GetExpHardSpherePhase(jNum,chNum);
	complex coulombPhaseEx = point->GetExpCoulombPhase(jNum,chpNum);
	complex hardSpherePhaseEx = point->GetExpHardSpherePhase(jNum,chpNum);
	complex sqrtPenEn = point->GetSqrtPenetrability(jNum,chNum);
	complex sqrtPenEx = point->GetSqrtPenetrability(jNum,chpNum);
	
	complex uphase = coulombPhaseEn * hardSpherePhaseEn * coulombPhaseEx * hardSpherePhaseEx;
	complex umatrix = 2.0*complex(0.0,1.0) * sqrtPenEn * sqrtPenEx *
	  this->GetUBilinear(jNum,chNum,chpNum);

	complex tphase = coulombPhaseEn * coulombPhaseEn;
	complex tmatrix;
	if(isRMC) this->AddTMatrixElement(k,m,complex(0.0,-1.0)*umatrix,ir);
	else {
	  if(chNum==chpNum) {
	    tmatrix=tphase-uphase*(1.0+umatrix);
	  } else tmatrix=-uphase*umatrix;
	  this->AddTMatrixElement(k,m,tmatrix);
	}
      }
      for(int m=1;m<=theDecay->GetKGroup(k)->NumECMGroups();m++) {
	ECMGroup *theECMGroup=theDecay->GetKGroup(k)->GetECMGroup(m);
	ALevel *finalLevel=compound()->GetJGroup(theECMGroup->GetJGroupNum())
	  ->GetLevel(theECMGroup->GetLevelNum());
	double ecNormParam=finalLevel->GetFitGamma(theECMGroup->GetFinalChannel())*
	  finalLevel->GetSqrtNFFactor()*finalLevel->GetECConversionFactor(theECMGroup->GetFinalChannel());
	// Use energy-shift aware EC amplitude calculation if energy shifts are active
	complex ecAmplitude;
	if(configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
	  // Use the new method that accounts for energy shifts through interpolation
	  ecAmplitude = point->GetECAmplitudeWithShift(k, m, compound(), configure());
	} else {
	  // Use cached amplitude for backward compatibility
	  ecAmplitude = point->GetECAmplitude(k, m);
	}
	complex tmatrix = ecNormParam * ecAmplitude;
	if(theECMGroup->IsChannelCapture()) {
	  int internalChannel=theECMGroup->GetIntChannelNum();
	  MGroup *chanMGroup=compound()->GetPair(aa)->GetDecay(theECMGroup->GetChanCapDecay())
	    ->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
	  int maskChannel = (configure().paramMask & Config::IGNORE_ZERO_WIDTHS) ?
	    internalChannel : 0;
	  complex umatrix=2.0*complex(0.0,1.0)*
	    point->GetSqrtPenetrability(chanMGroup->GetJNum(),chanMGroup->GetChNum())*
	    this->GetChannelCaptureBilinear(chanMGroup->GetJNum(),chanMGroup->GetChNum(),
	                                    chanMGroup->GetChpNum(),maskChannel);
	  tmatrix=tmatrix*umatrix;
	}
	this->AddECTMatrixElement(k,m,tmatrix);
      }
    }
  }
}

/*!
 * This function adds an inverse A-Matrix element specified by positions in the JGroup and ALevel vectors.
 */

void AMatrixFunc::AddAInvMatrixElement(int jGroupNum, int lambdaNum, int muNum, complex aMatrixElement) {
  // Basic bounds
  if (jGroupNum < 1 || jGroupNum > (int)a_inv_matrices_.size()) return;

  // If this J-group wasn’t included, skip
  if (level_active_index_.empty() || level_active_index_[jGroupNum-1].empty()) return;

  // Translate original level indices to compact active indices
  const auto& map = level_active_index_[jGroupNum-1];
  if (lambdaNum < 1 || lambdaNum >= (int)map.size()) return;
  if (muNum     < 1 || muNum     >= (int)map.size()) return;

  const int row = map[lambdaNum];
  const int col = map[muNum];
  if (row == 0 || col == 0) return; // either level is inactive; logic says: ignore

  // Validate value
  if (!std::isfinite(aMatrixElement.real()) || !std::isfinite(aMatrixElement.imag())) return;

  // Direct assignment into pre-sized dense matrix
  a_inv_matrices_[jGroupNum-1][row-1][col-1] = aMatrixElement;
}

/*!
 * Reverse-mode adjoint of one energy point's cross section (angle-integrated or
 * differential), accumulating dlnL/dE and dlnL/dgamma into `accum`.
 *
 * Adjoint convention (consistent throughout, see PLAN.md "Gotchas"):
 *   For a complex intermediate z and the real scalar lnL we store the cotangent
 *     zbar := dlnL/dRe(z) + i dlnL/dIm(z).
 *   - A linear map z = c * a (c complex constant) back-propagates as
 *         abar += conj(c) * zbar.
 *   - The contribution of z to a *real* parameter p is
 *         dlnL/dp += Re( conj(zbar) * dz/dp ).
 *   - The holomorphic inverse A = M^{-1} back-propagates as
 *         Mbar = -A^H * Abar * A^H.
 *
 * Forward chain reproduced here (per point, standard A-matrix):
 *   T(k,m) = [tphase] - uphase ( [1] + umatrix )         ([..] only if chNum==chpNum)
 *   umatrix = sum_{l,l'} 2i sqrtP_en sqrtP_ex gamma_en,l gamma_ex,l' A_j[l][l']
 *   A_j = M_j^{-1},   M_j[r][c] = (E_r - inE) delta_rc - sum_ch gamma_r,ch gamma_c,ch L_ch
 * with the cross-section -> T step depending on whether the point is
 * angle-integrated (quadratic |T_temp|^2 sum) or differential (Legendre
 * interference sum + elastic Coulomb interference).  Returns false (touching
 * nothing) for configurations outside the supported standard path.
 */
bool AMatrixFunc::PointAdjoint(EPoint* point, double fitBar, GradAccum& accum,
                              const vector_matrix_r* shiftDeriv, int xsComponent) {
  const complex I(0.0, 1.0);

  // --- Common gating. ---
  if(configure().paramMask & Config::USE_RMC_FORMALISM) return false;
  if(point->IsAngularDist()) return false;        // angular-dist coefficients
  // E1/E2 component selection only applies to the angle-integrated capture XS.
  if(xsComponent != 0 && (point->IsDifferential() || point->IsPhase())) return false;

  // Brune formalism: the boundary/shift terms change the M-construction kernel
  // (handled in the structural gamma step via cached shift values) and add
  // explicit E_lambda dependence through S(E_lambda) (handled in the structural
  // energy step using shiftDeriv = dS/dE).
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  int aa = compound()->GetPairNumFromKey(point->GetEntranceKey());
  int exitPairNum = compound()->GetPairNumFromKey(point->GetExitKey());
  // Beta-delayed particle emission is supported below; its observable is a plain
  // incoherent sum over pathways, so only the T cotangents differ from the
  // cross-section case. E1/E2 component selection has no meaning for it.
  const bool isBetaDelayed = (compound()->GetPair(aa)->GetPType() == 20);
  if(isBetaDelayed && xsComponent != 0) return false;

  // Active-level compaction unsupported here: require every level of every
  // in-R-matrix JGroup to itself be in the R matrix (so original idx == active).
  for(int j = 1; j <= compound()->NumJGroups(); j++) {
    JGroup* jg = compound()->GetJGroup(j);
    if(!jg->IsInRMatrix()) continue;
    for(int la = 1; la <= jg->NumLevels(); la++)
      if(!jg->GetLevel(la)->IsInRMatrix()) return false;
  }

  // Locate the decay matching the exit pair (mirror CalculateCrossSection).
  int ir = 0;
  while(ir < compound()->GetPair(aa)->NumDecays()) {
    ir++;
    if(compound()->GetPair(aa)->GetDecay(ir)->GetPairNum() == exitPairNum) break;
  }
  Decay* theDecay = compound()->GetPair(aa)->GetDecay(ir);

  // External capture (both direct and channel-capture) is supported.
  bool hasEC = false;
  for(int k = 1; k <= theDecay->NumKGroups(); k++)
    if(theDecay->GetKGroup(k)->NumECMGroups() > 0) { hasEC = true; break; }

  // Beam energy in the compound system (mirrors FillMatrices); needed for the
  // Brune off-diagonal boundary/shift kernel.
  double inEnergy;
  if(compound()->GetPair(aa)->GetPType()==20) {
    int ep = compound()->GetPairNumFromKey(point->GetExitKey());
    inEnergy = point->GetCMEnergy() + compound()->GetPair(ep)->GetSepE()
             + compound()->GetPair(ep)->GetExE();
  } else {
    inEnergy = point->GetCMEnergy() + compound()->GetPair(aa)->GetSepE()
             + compound()->GetPair(aa)->GetExE();
  }

  const double geom = point->GetGeometricalFactor();
  const double i1i2 = compound()->GetPair(aa)->GetI1I2Factor();
  const bool isDiff = point->IsDifferential();
  const bool isPhase = point->IsPhase();
  if(isPhase && aa != ir) return false;  // phase output is for elastic only

  // === Cross-section -> T:  build the T-matrix cotangents. ===
  // tBar  : cotangents of the internal T(k,m).
  // ecBar : cotangents of the external-capture T(k,ecm).
  int nK = theDecay->NumKGroups();
  std::vector<vector_c> tBar(nK), ecBar(nK);
  for(int k = 1; k <= nK; k++) {
    tBar[k-1].assign(theDecay->GetKGroup(k)->NumMGroups(), complex(0.0,0.0));
    ecBar[k-1].assign(theDecay->GetKGroup(k)->NumECMGroups(), complex(0.0,0.0));
  }

  // Keep only the contributions feeding the requested cross-section component
  // (0 = full, 1 = E1, 2 = E2); E1/E2 select electric channels of that L.
  auto includeInternal = [&](MGroup* mg) -> bool {
    if(xsComponent == 0) return true;
    AChannel* ex = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChpNum());
    return ex->GetRadType() == 'E' && ex->GetL() == xsComponent;
  };
  auto includeEC = [&](ECMGroup* ecg) -> bool {
    if(xsComponent == 0) return true;
    return ecg->GetRadType() == 'E' && ecg->GetMult() == xsComponent;
  };

  if(isBetaDelayed) {
    // ---- Beta-delayed particle emission. GenMatrixFunc forms
    //        sum = sum_{k,m} 25 |T(k,m)|^2 ,   model = Re(sum)/100
    //      with no geometrical or spin factor and no external-capture term --
    //      an incoherent sum over pathways, with no interference between them.
    //      Hence d(model)/dT*(k,m) = (25/100) T(k,m), and the factor of two is
    //      the cotangent convention used throughout this function. ----
    for(int k = 1; k <= nK; k++)
      for(int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++)
        tBar[k-1][m-1] += fitBar * (2.0 * 25.0 / 100.0) * this->GetTMatrixElement(k, m);
  } else if(point->IsAnalyzingPower()) {
    // ---- Vector analyzing power. A_y is built from the channel-spin amplitude
    //      matrix, which is linear in T(k,m) with the coefficients AddPathway
    //      applies, so the whole derivative is a contraction against those same
    //      coefficients. The Coulomb amplitude carries no parameter dependence
    //      and drops out. ----
    if(point->NumSubPoints() > 0) return false;   // see the note below
    if(compound()->GetPair(aa)->GetPType() != 0 ||
       compound()->GetPair(exitPairNum)->GetPType() != 0) return false;

    Polarization::AmplitudeMatrix M(compound(), point, aa, exitPairNum);
    for(int k = 1; k <= nK; k++) {
      for(int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
        M.AddPathway(mg->GetJNum(), mg->GetChNum(), mg->GetChpNum(),
                     this->GetTMatrixElement(k, m));
      }
    }
    if(aa == exitPairNum) M.AddCoulomb(point->GetCoulombAmplitude());
    if(M.size() == 0) return false;

    const std::vector<complex> bar = M.AnalyzingPowerBar();
    for(int k = 1; k <= nK; k++) {
      for(int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
        if(!includeInternal(mg)) continue;
        tBar[k-1][m-1] += fitBar * M.PathwayAdjoint(mg->GetJNum(), mg->GetChNum(),
                                                    mg->GetChpNum(), bar);
      }
    }
  } else if(isPhase) {
    // ---- Phase shift: model = (90/pi) * arg(U) [+ const for identical pairs],
    //      U = sum_{matching k,m} (expCP^2 - T(k,m)) / expCP^2, matched by
    //      (J, l, elastic channel) = (segment J, segment L). ----
    double segJ = point->GetJ();
    int segL = point->GetL();
    struct PM { int k; int m; complex invExpCP2; };
    std::vector<PM> matches;
    complex U(0.0,0.0);
    for(int k = 1; k <= nK; k++) {
      for(int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
        MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
        JGroup* jg = compound()->GetJGroup(mg->GetJNum());
        AChannel* en = jg->GetChannel(mg->GetChNum());
        AChannel* ex = jg->GetChannel(mg->GetChpNum());
        if(jg->GetJ() == segJ && en->GetL() == segL && en == ex) {
          complex expCP = point->GetExpCoulombPhase(mg->GetJNum(), mg->GetChNum());
          complex expCP2 = expCP * expCP;
          U += (expCP2 - this->GetTMatrixElement(k, m)) / expCP2;
          matches.push_back(PM{k, m, 1.0 / expCP2});
        }
      }
    }
    double absU2 = std::norm(U);
    if(absU2 > 0.0) {
      // phase = (90/pi) arg(U) -> Ubar = fitBar (90/pi)(-Im U + i Re U)/|U|^2.
      complex Ubar = fitBar * (90.0 / pi) *
                     complex(-std::imag(U), std::real(U)) / absU2;
      // U = sum [1 - T/expCP^2]  ->  dU/dT = -1/expCP^2.
      for(const PM& pm : matches)
        tBar[pm.k-1][pm.m-1] += conj(-pm.invExpCP2) * Ubar;
    }
  } else if(!isDiff) {
    // ---- Angle-integrated: model = (1/100) sum_k sum_temp w_temp |T_temp|^2,
    //      w_temp = geom (2J+1) I1I2,  T_temp = sum_m T(k,m) sharing (J,l,l'). ----
    for(int k = 1; k <= nK; k++) {
      // Mirror the angle-integrated UPOS over-counting guard.
      if(aa != ir && compound()->GetPair(exitPairNum)->GetPType() == 0) {
        if(theDecay->GetKGroup(k)->GetSp() != theDecay->GetKGroup(k)->GetSp2()) continue;
      }
      int numM = theDecay->GetKGroup(k)->NumMGroups();
      int numEC = theDecay->GetKGroup(k)->NumECMGroups();
      struct Group { double j; int l; int lp; complex tsum; };
      std::vector<Group> groups;
      std::vector<int> mGroupIdx(numM + 1, -1);
      std::vector<int> ecGroupIdx(numEC + 1, -1);
      auto findGroup = [&](double jv, int lv, int lpv) -> int {
        for(int g = 0; g < (int)groups.size(); g++)
          if(groups[g].j == jv && groups[g].l == lv && groups[g].lp == lpv) return g;
        groups.push_back(Group{jv, lv, lpv, complex(0.0,0.0)});
        return (int)groups.size() - 1;
      };
      // Internal pathways T(k,m), grouped by (J,l,l').
      for(int m = 1; m <= numM; m++) {
        MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
        if(!includeInternal(mg)) continue;
        JGroup* jg = compound()->GetJGroup(mg->GetJNum());
        int gi = findGroup(jg->GetJ(),
                           jg->GetChannel(mg->GetChNum())->GetL(),
                           jg->GetChannel(mg->GetChpNum())->GetL());
        groups[gi].tsum += this->GetTMatrixElement(k, m);
        mGroupIdx[m] = gi;
      }
      // External-capture pathways T(k,ecm), grouped by (J, L, mult) and merged
      // with the matching internal group (forward: same temp T-matrix list).
      for(int m = 1; m <= numEC; m++) {
        ECMGroup* ecg = theDecay->GetKGroup(k)->GetECMGroup(m);
        if(!includeEC(ecg)) continue;
        int gi = findGroup(ecg->GetJ(), ecg->GetL(), ecg->GetMult());
        groups[gi].tsum += this->GetECTMatrixElement(k, m);
        ecGroupIdx[m] = gi;
      }
      std::vector<complex> tBarGroup(groups.size());
      for(int g = 0; g < (int)groups.size(); g++) {
        double w = geom * (2.0 * groups[g].j + 1.0) * i1i2;
        tBarGroup[g] = fitBar * (2.0 / 100.0) * w * groups[g].tsum;
      }
      for(int m = 1; m <= numM; m++)  if(mGroupIdx[m]  >= 0) tBar[k-1][m-1]  = tBarGroup[mGroupIdx[m]];
      for(int m = 1; m <= numEC; m++) if(ecGroupIdx[m] >= 0) ecBar[k-1][m-1] = tBarGroup[ecGroupIdx[m]];
    }
  } else {
    // ---- Differential: model = (Re(CT)+Re(RT)+Re(IT))/100. ----
    //   RT = (geom I1I2 rtFactor/pi) Re( sum_{kL,inter} weight T1 conj(T2) P_L )
    //   IT = (geom itFactor/sqrt(pi)) Re( i sum_{k,m: en==ex} statSpin coulombAmp conj(T) P_l )  (elastic)
    //   CT = const (no T dependence).
    // Supports the standard case (entrance==exit pair, or gamma exit) with
    // external-capture interferences (RR/ER/RE/EE), and the inelastic-particle
    // UPOS / normal-angular-distribution case.
    PPair* exitPairDiff = compound()->GetPair(exitPairNum);
    double rtFactor = 1.0, itFactor = 1.0;
    if(aa == ir && compound()->GetPair(aa)->IsIdentical()) { rtFactor = 4.0; itFactor = 2.0; }
    double cR = geom * i1i2 * rtFactor / (pi * 100.0);

    // Fetch a T value (internal 'R' or external-capture 'E') and route a
    // cotangent contribution to the right accumulator.
    auto Tval = [&](int K, int idx, char rad) -> complex {
      return (rad == 'E') ? this->GetECTMatrixElement(K, idx)
                          : this->GetTMatrixElement(K, idx);
    };
    auto Tadd = [&](int K, int idx, char rad, complex v) {
      if(rad == 'E') ecBar[K-1][idx-1] += v; else tBar[K-1][idx-1] += v;
    };
    // model += W Re(T1 conj T2)  ->  T1bar += fitBar W T2 ;  T2bar += fitBar W T1.
    auto addQuad = [&](int K, int m1, char r1, int m2, char r2, double W) {
      complex T1 = Tval(K, m1, r1), T2 = Tval(K, m2, r2);
      Tadd(K, m1, r1, fitBar * W * T2);
      Tadd(K, m2, r2, fitBar * W * T1);
    };

    const bool inelasticParticle = (aa != ir && exitPairDiff->GetPType() == 0);
    for(int kL = 1; kL <= theDecay->NumKLGroups(); kL++) {
      int K = theDecay->GetKLGroup(kL)->GetK();
      int lOrder = theDecay->GetKLGroup(kL)->GetLOrder();
      double legP = point->GetLegendreP(lOrder);
      KGroup* kg = theDecay->GetKGroup(K);
      for(int inter = 1; inter <= theDecay->GetKLGroup(kL)->NumInterferences(); inter++) {
        Interference* itf = theDecay->GetKLGroup(kL)->GetInterference(inter);
        std::string type = itf->GetInterferenceType();  // "RR","ER","RE","EE"
        if(type.size() != 2) return false;
        char rad1 = type[0], rad2 = type[1];            // 'R' internal, 'E' external
        int m1 = itf->GetM1(), m2 = itf->GetM2();

        if(inelasticParticle) {
          // Inelastic particle exit: pick the normal or UPOS sub-case.
          double sp1 = kg->GetSp(), sp2 = kg->GetSp2();
          int lp1 = compound()->GetJGroup(kg->GetMGroup(m1)->GetJNum())
                      ->GetChannel(kg->GetMGroup(m1)->GetChpNum())->GetL();
          int lp2 = compound()->GetJGroup(kg->GetMGroup(m2)->GetJNum())
                      ->GetChannel(kg->GetMGroup(m2)->GetChpNum())->GetL();
          if(sp1 == sp2 && !point->IsUPOS()) {
            addQuad(K, m1, rad1, m2, rad2, cR * itf->GetZ1Z2() * legP);
          }
          if(lp1 == lp2 && point->IsUPOS() && type == "RR") {
            // R_L coefficient (constant w.r.t. parameters), even lOrder only.
            double finalL = (double)point->GetSecondaryDecayL();
            double Ic = point->GetIc();
            double j2f = compound()->GetPair(ir)->GetJ(2);
            double delta = point->GetDelta();
            double R_L = 0.0;
            if((int)lOrder % 2 == 0) {
              R_L = this->GetRk(j2f, finalL, finalL, Ic, lOrder);
              if(delta != 0.0) {
                double finalLp = finalL + 1.0;
                double R_LLp  = this->GetRk(j2f, finalL, finalLp, Ic, lOrder);
                double R_LpLp = this->GetRk(j2f, finalLp, finalLp, Ic, lOrder);
                if(R_LpLp != 0.0)
                  R_L = (R_L + 2.0*delta*R_LLp + delta*delta*R_LpLp) / (1.0 + delta*delta);
              }
            }
            double W = cR * itf->GetZ1Z2_UPOS() * std::pow(2.0*lOrder+1.0, 0.5) / 4.0 * R_L * legP;
            addQuad(K, m1, 'R', m2, 'R', W);
          }
        } else {
          // Standard case (entrance==exit pair, or gamma exit).
          addQuad(K, m1, rad1, m2, rad2, cR * itf->GetZ1Z2() * legP);
        }
      }
    }

    if(aa == ir) {  // elastic Coulomb-nuclear interference
      complex coulombAmp = point->GetCoulombAmplitude();
      complex cI = I * geom * itFactor / (std::sqrt(pi) * 100.0);
      for(int k = 1; k <= nK; k++) {
        for(int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
          MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
          AChannel* en = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChNum());
          AChannel* ex = compound()->GetJGroup(mg->GetJNum())->GetChannel(mg->GetChpNum());
          if(en == ex) {
            int l = en->GetL();
            complex B = cI * mg->GetStatSpinFactor() * coulombAmp * point->GetLegendreP(l);
            // model += Re(B conj T)  ->  Tbar += fitBar B.
            tBar[k-1][m-1] += fitBar * B;
          }
        }
      }
    }
  }

  // Per-JGroup A-matrix cotangent accumulators (active == original indexing).
  std::vector<matrix_c> aBar(compound()->NumJGroups());
  std::vector<bool> aBarUsed(compound()->NumJGroups(), false);
  const complex I2 = 2.0 * I;
  auto ensureABar = [&](int jNum, int numLevels) {
    if(!aBarUsed[jNum-1]) {
      aBar[jNum-1].assign(numLevels, vector_c(numLevels, complex(0.0,0.0)));
      aBarUsed[jNum-1] = true;
    }
  };

  // === External-capture T -> parameters. ===
  // Forward:  T_ec(k,ecm) = gamma_fc * sqrtNF * conv * ecAmplitude  [* umatrix_cc],
  //   gamma_fc = final level's reduced width for the final channel (explicit),
  //   sqrtNF   = final level's 1/sqrt(NF) factor (depends on its gammas),
  //   conv     = EC conversion factor (constant),
  //   ecAmplitude treated constant w.r.t. E/gamma (energy shifts handled by FD).
  // For channel capture an extra factor umatrix_cc = sum 2i sqrtP gamma gamma A
  // over an internal pathway multiplies T_ec; it back-propagates to that
  // pathway's A-matrix and gammas just like the standard umatrix.
  if(hasEC) {
    const bool useEC = (configure().paramMask & Config::USE_EXTERNAL_CAPTURE);
    const bool ignoreZero = (configure().paramMask & Config::IGNORE_ZERO_WIDTHS);
    for(int k = 1; k <= nK; k++) {
      int numEC = theDecay->GetKGroup(k)->NumECMGroups();
      for(int m = 1; m <= numEC; m++) {
        complex eb = ecBar[k-1][m-1];
        if(eb == complex(0.0,0.0)) continue;
        ECMGroup* ecg = theDecay->GetKGroup(k)->GetECMGroup(m);
        int jf = ecg->GetJGroupNum();
        int lf = ecg->GetLevelNum();
        int fc = ecg->GetFinalChannel();
        ALevel* finalLevel = compound()->GetJGroup(jf)->GetLevel(lf);
        JGroup* jgf = compound()->GetJGroup(jf);
        double sqrtNF = finalLevel->GetSqrtNFFactor();
        double conv = finalLevel->GetECConversionFactor(fc);
        complex ecAmp = useEC ? point->GetECAmplitudeWithShift(k, m, compound(), configure())
                              : point->GetECAmplitude(k, m);

        // base = gamma_fc * sqrtNF * conv * ecAmp  (the part without umatrix_cc).
        double gamma_fc = finalLevel->GetFitGamma(fc);
        complex base = gamma_fc * sqrtNF * conv * ecAmp;

        // "effective amplitude" carrying the umatrix_cc factor for the
        // gamma_fc / sqrtNF derivatives (= ecAmp for direct capture).
        complex effAmp = ecAmp;

        if(ecg->IsChannelCapture()) {
          int internalChannel = ecg->GetIntChannelNum();
          MGroup* ccMg = compound()->GetPair(aa)->GetDecay(ecg->GetChanCapDecay())
                           ->GetKGroup(ecg->GetChanCapKGroup())->GetMGroup(ecg->GetChanCapMGroup());
          int jcc = ccMg->GetJNum();
          int ccChN = ccMg->GetChNum();
          int ccChpN = ccMg->GetChpNum();
          JGroup* jgcc = compound()->GetJGroup(jcc);
          int ccLevels = jgcc->NumLevels();
          complex sqrtPcc = point->GetSqrtPenetrability(jcc, ccChN);  // single sqrtP

          // Forward umatrix_cc (mirrors the IGNORE_ZERO_WIDTHS skips).
          complex umatrix_cc(0.0,0.0);
          for(int la = 1; la <= ccLevels; la++) {
            if(!jgcc->GetLevel(la)->IsInRMatrix()) continue;
            if(internalChannel && ignoreZero &&
               std::fabs(jgcc->GetLevel(la)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
            for(int lap = 1; lap <= ccLevels; lap++) {
              if(!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if(internalChannel && ignoreZero &&
                 std::fabs(jgcc->GetLevel(lap)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
              umatrix_cc += I2 * sqrtPcc *
                jgcc->GetLevel(la)->GetFitGamma(ccChN) *
                jgcc->GetLevel(lap)->GetFitGamma(ccChpN) *
                this->GetAMatrixElement(jcc, la, lap);
            }
          }
          effAmp = ecAmp * umatrix_cc;

          // umatrix_cc cotangent:  T_ec = base * umatrix_cc -> ubar = conj(base)*eb.
          complex uBarcc = conj(base) * eb;
          ensureABar(jcc, ccLevels);
          for(int la = 1; la <= ccLevels; la++) {
            if(!jgcc->GetLevel(la)->IsInRMatrix()) continue;
            if(internalChannel && ignoreZero &&
               std::fabs(jgcc->GetLevel(la)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
            double gEn = jgcc->GetLevel(la)->GetFitGamma(ccChN);
            for(int lap = 1; lap <= ccLevels; lap++) {
              if(!jgcc->GetLevel(lap)->IsInRMatrix()) continue;
              if(internalChannel && ignoreZero &&
                 std::fabs(jgcc->GetLevel(lap)->GetFitGamma(internalChannel)) < 1.0e-8) continue;
              double gEx = jgcc->GetLevel(lap)->GetFitGamma(ccChpN);
              complex Aval = this->GetAMatrixElement(jcc, la, lap);
              complex w = I2 * sqrtPcc * gEn * gEx;          // weight of A[la][lap]
              aBar[jcc-1][la-1][lap-1] += conj(w) * uBarcc;
              accum.AddGamma(jcc, la,  ccChN,  std::real(conj(uBarcc) * (I2 * sqrtPcc * gEx * Aval)));
              accum.AddGamma(jcc, lap, ccChpN, std::real(conj(uBarcc) * (I2 * sqrtPcc * gEn * Aval)));
            }
          }
        }

        // Explicit gamma_fc:  d T_ec / d gamma_fc = sqrtNF * conv * effAmp.
        accum.AddGamma(jf, lf, fc, std::real(conj(eb) * (sqrtNF * conv * effAmp)));

        // sqrtNF couples to every NF channel d of the final level:
        //   nFSum = 1 + sum_d A_d gamma_d^2,  A_d = 2 chRad redMass uconv/hbarc^2 NF_d,
        //   sqrtNF = nFSum^{-1/2},  d sqrtNF / d gamma_d = -A_d gamma_d sqrtNF^3.
        complex pre = gamma_fc * conv * effAmp;  // d T_ec / d sqrtNF
        int numNF = finalLevel->NumNFIntegrals();
        double sqrtNF3 = sqrtNF * sqrtNF * sqrtNF;
        for(int d = 1; d <= numNF; d++) {
          PPair* pd = compound()->GetPair(jgf->GetChannel(d)->GetPairNum());
          double Ad = 2.0 * pd->GetChRad() * pd->GetRedMass() * uconv / (hbarc * hbarc)
                    * finalLevel->GetNFIntegral(d);
          double gd = finalLevel->GetFitGamma(d);
          double dSqrtNF = -Ad * gd * sqrtNF3;
          accum.AddGamma(jf, lf, d, std::real(conj(eb) * (pre * dSqrtNF)));
        }
      }
    }
  }

  // === T -> U -> A (+ explicit gamma), back-propagating tBar over all pathways. ===
  for(int k = 1; k <= nK; k++) {
    int numM = theDecay->GetKGroup(k)->NumMGroups();
    for(int m = 1; m <= numM; m++) {
      complex tb = tBar[k-1][m-1];
      if(tb == complex(0.0,0.0)) continue;
      MGroup* mg = theDecay->GetKGroup(k)->GetMGroup(m);
      int jNum = mg->GetJNum();
      int chNum = mg->GetChNum();
      int chpNum = mg->GetChpNum();
      JGroup* jg = compound()->GetJGroup(jNum);

      complex uphase = point->GetExpCoulombPhase(jNum, chNum) *
                       point->GetExpHardSpherePhase(jNum, chNum) *
                       point->GetExpCoulombPhase(jNum, chpNum) *
                       point->GetExpHardSpherePhase(jNum, chpNum);
      complex sqrtPenEn = point->GetSqrtPenetrability(jNum, chNum);
      complex sqrtPenEx = point->GetSqrtPenetrability(jNum, chpNum);

      // d T / d umatrix = -uphase  ->  ubar = conj(-uphase) * tBar.
      complex uBar = conj(-uphase) * tb;

      int numLevels = jg->NumLevels();
      if(!aBarUsed[jNum-1]) {
        aBar[jNum-1].assign(numLevels, vector_c(numLevels, complex(0.0,0.0)));
        aBarUsed[jNum-1] = true;
      }
      complex pen2 = 2.0 * I * sqrtPenEn * sqrtPenEx;  // common factor of w_{l,l'}

      for(int la = 1; la <= numLevels; la++) {
        if(!jg->GetLevel(la)->IsInRMatrix()) continue;
        double gammaEn = jg->GetLevel(la)->GetFitGamma(chNum);
        for(int lap = 1; lap <= numLevels; lap++) {
          if(!jg->GetLevel(lap)->IsInRMatrix()) continue;
          double gammaEx = jg->GetLevel(lap)->GetFitGamma(chpNum);
          complex Aval = this->GetAMatrixElement(jNum, la, lap);
          complex w = pen2 * gammaEn * gammaEx;          // weight of A[la][lap]

          // A cotangent:  Abar += conj(w) * ubar.
          aBar[jNum-1][la-1][lap-1] += conj(w) * uBar;

          // Explicit gamma paths (umatrix depends on gamma_en and gamma_ex too):
          //   d umatrix / d gamma_en,la += pen2 * gamma_ex * A
          //   d umatrix / d gamma_ex,lap += pen2 * gamma_en * A
          accum.AddGamma(jNum, la,  chNum,  std::real(conj(uBar) * (pen2 * gammaEx * Aval)));
          accum.AddGamma(jNum, lap, chpNum, std::real(conj(uBar) * (pen2 * gammaEn * Aval)));
        }
      }
    }
  }

  // === A -> M lemma + structural E / gamma, once per JGroup. ===
  for(int jNum = 1; jNum <= compound()->NumJGroups(); jNum++) {
    if(!aBarUsed[jNum-1]) continue;
    JGroup* jg = compound()->GetJGroup(jNum);
    const matrix_c& A = this->GetAMatrix(jNum);
    const matrix_c& Ab = aBar[jNum-1];
    int n = (int)A.size();
    if(n == 0) continue;

    // Conjugate transpose A^H.
    matrix_c Ah(n, vector_c(n, complex(0.0,0.0)));
    for(int r = 0; r < n; r++)
      for(int c = 0; c < n; c++)
        Ah[r][c] = conj(A[c][r]);

    // Mbar = -A^H * Abar * A^H.
    matrix_c tmp(n, vector_c(n, complex(0.0,0.0)));  // tmp = Abar * A^H
    for(int r = 0; r < n; r++)
      for(int c = 0; c < n; c++) {
        complex s(0.0,0.0);
        for(int t = 0; t < n; t++) s += Ab[r][t] * Ah[t][c];
        tmp[r][c] = s;
      }
    matrix_c Mbar(n, vector_c(n, complex(0.0,0.0)));
    for(int r = 0; r < n; r++)
      for(int c = 0; c < n; c++) {
        complex s(0.0,0.0);
        for(int t = 0; t < n; t++) s += Ah[r][t] * tmp[t][c];
        Mbar[r][c] = -s;
      }

    // Structural energy gradient.  The (E_r - inE) delta term gives
    //   d M[a][a]/d E_a = 1  ->  dlnL/dE_a += Re(Mbar[a][a]).
    for(int la = 1; la <= n; la++)
      accum.AddE(jNum, la, std::real(Mbar[la-1][la-1]));

    // Brune adds explicit E_lambda dependence through the boundary/shift kernel
    //   M[r][c] = ... - sum_{ch P} g_r,ch g_c,ch K_ch[r][c],
    //   K_ch[r][c] = L + B - D_ch[r][c],  D diagonal = S_r,  off-diag = Q_{rc}.
    // so dlnL/dE_nu += sum_{r,c} Re(conj(Mbar[r][c]) * sum_ch g_r g_c dD/dE_nu).
    // dD/dE only couples to E_r and E_c; using S'(E) = shiftDeriv.
    if(brune && shiftDeriv) {
      const std::vector<std::vector<double>>& sd = (*shiftDeriv)[jNum-1];
      int numChannels = jg->NumChannels();
      for(int ch = 1; ch <= numChannels; ch++) {
        if(jg->GetChannel(ch)->GetRadType() != 'P') continue;
        for(int r = 1; r <= n; r++) {
          double gr  = jg->GetLevel(r)->GetFitGamma(ch);
          if(gr == 0.0) continue;
          double Sr  = jg->GetLevel(r)->GetShiftFunction(ch);
          double dSr = sd[r-1][ch-1];
          double Er  = jg->GetLevel(r)->GetFitE();
          for(int c = 1; c <= n; c++) {
            double gc = jg->GetLevel(c)->GetFitGamma(ch);
            double g = gr * gc;
            if(g == 0.0) continue;
            if(r == c) {
              // D = S_r ;  dD/dE_r = S'_r.
              accum.AddE(jNum, r, std::real(Mbar[r-1][r-1]) * g * dSr);
            } else {
              double Sc  = jg->GetLevel(c)->GetShiftFunction(ch);
              double dSc = sd[c-1][ch-1];
              double Ec  = jg->GetLevel(c)->GetFitE();
              double den = Er - Ec;
              double Q   = (Sr * (inEnergy - Ec) - Sc * (inEnergy - Er)) / den;
              double dQdEr = (dSr * (inEnergy - Ec) + Sc - Q) / den;
              double dQdEc = (-Sr - dSc * (inEnergy - Er) + Q) / den;
              double reM = std::real(Mbar[r-1][c-1]);
              accum.AddE(jNum, r, reM * g * dQdEr);
              accum.AddE(jNum, c, reM * g * dQdEc);
            }
          }
        }
      }
    }

    // Structural gamma gradient.  M[r][c] = (E_r-inE) delta - sum_ch g_r,ch g_c,ch K_ch[r][c],
    // so  dlnL/d gamma_{a,d} += -sum_c Re( conj(Mbar[a][c]+Mbar[c][a]) * K_d[a][c] ) * gamma_{c,d}.
    // The kernel K_d[r][c] is the channel-d contribution to the level matrix:
    //   non-Brune (or non-'P'):  K_d = L_d  (the Lo-matrix element, indep. of r,c)
    //   Brune 'P' channel:       K_d[r][c] = L_d + B_d - D_d[r][c],
    //     D_d[r][c] = (r==c) ? S_r,d
    //                        : (S_r,d (inE-E_c) - S_c,d (inE-E_r)) / (E_r - E_c)
    // exactly mirroring FillMatrices.
    int numChannels = jg->NumChannels();
    for(int d = 1; d <= numChannels; d++) {
      complex Ld = point->GetLoElement(jNum, d);
      AChannel* channel = jg->GetChannel(d);
      bool bruneP = brune && (channel->GetRadType() == 'P');
      double Bd = bruneP ? channel->GetBoundaryCondition() : 0.0;
      for(int a = 1; a <= n; a++) {
        double g = 0.0;
        for(int c = 1; c <= n; c++) {
          complex Kd = Ld;
          if(bruneP) {
            double Sa = jg->GetLevel(a)->GetShiftFunction(d);
            double Sc = jg->GetLevel(c)->GetShiftFunction(d);
            double Dd;
            if(a == c) Dd = Sa;
            else {
              double Ea = jg->GetLevel(a)->GetFitE();
              double Ec = jg->GetLevel(c)->GetFitE();
              Dd = (Sa * (inEnergy - Ec) - Sc * (inEnergy - Ea)) / (Ea - Ec);
            }
            Kd = Ld + Bd - Dd;
          }
          complex sym = Mbar[a-1][c-1] + Mbar[c-1][a-1];
          double gammaCd = jg->GetLevel(c)->GetFitGamma(d);
          g += -std::real(conj(sym) * Kd) * gammaCd;
        }
        accum.AddGamma(jNum, a, d, g);
      }
    }
  }

  return true;
}
