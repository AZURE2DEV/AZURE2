#ifndef GRADFREEMASK_H
#define GRADFREEMASK_H

#include "AMatrixFunc.h"
#include "AZUREGrad.h"

class CompoundSystem;

/// The level parameters the formal basis's adjoint must resolve, per J-group: those free in the packed vector, the
/// photon amplitudes of a gamma_decay level with a free input, and the derived slots of isospin levels with a free
/// source (GradFreeMask.cpp).  For the residual Jacobian and the model gradients of AZUREGrad.cpp.
AMatrixFunc::FreeMask BuildFreeMask(const CompoundSystem *compound, const ParamIndexMap &pmap);

#endif
