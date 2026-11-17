#ifndef INFERENCETYPES_H
#define INFERENCETYPES_H

/*!
 * The engine's names of the inference/ library's types (docs/upgrade/inference.md): what a fit parameter is
 * (ParameterKind.h), where the level block and every other parameter sit in the vector (LevelParameterLayout.h,
 * ParameterIndexMap.h), the posterior's support (PosteriorSupport.h), one Gaussian prior (GaussianPrior.h), the list of
 * the diagonal priors (DiagonalPriors.h) and the prior covariance (PriorCovariance.h).  The rule of the objective's
 * term list (TermList.h) is named by Objective.h, which instantiates it.  The library is below the engine and
 * includes nothing of it; the engine and the binding spell its types without the namespace, as ChannelSpace.h does
 * for physics/'s constants.
 */
#include "azure2/inference/DiagonalPriors.h"
#include "azure2/inference/GaussianPrior.h"
#include "azure2/inference/LevelParameterLayout.h"
#include "azure2/inference/ParameterIndexMap.h"
#include "azure2/inference/ParameterKind.h"
#include "azure2/inference/PosteriorSupport.h"
#include "azure2/inference/PriorCovariance.h"

using azure2::inference::AddPriorsChi2;
using azure2::inference::AddPriorsCurvature;
using azure2::inference::AddPriorsGradient;
using azure2::inference::DiagonalPrior;
using azure2::inference::DiagonalPriors;
using azure2::inference::GaussianPriorChi2;
using azure2::inference::GaussianPriorCurvature;
using azure2::inference::GaussianPriorGradient;
using azure2::inference::IsLevelSlot;
using azure2::inference::IsRMatrixParameter;
using azure2::inference::KindRow;
using azure2::inference::LevelParameterLayout;
using azure2::inference::ParamDesc;
using azure2::inference::ParamIndexMap;
using azure2::inference::ParamKind;
using azure2::inference::ParamKindRow;
using azure2::inference::ParamKindTable;
using azure2::inference::ParameterRecord;
using azure2::inference::PosteriorSupport;
using azure2::inference::PriorCovariance;
using azure2::inference::PriorsChi2;
using azure2::inference::PriorSource;

#endif  // INFERENCETYPES_H
