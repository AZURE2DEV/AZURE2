#ifndef INFERENCETYPES_H
#define INFERENCETYPES_H

/*!
 * The engine's names of the inference/ library's types (docs/upgrade/inference.md): what a fit parameter is
 * (ParameterKind.h), where the level block and every other parameter sit in the vector (LevelParameterLayout.h,
 * ParameterIndexMap.h) and the posterior's support (PosteriorSupport.h).  The library is below the engine and
 * includes nothing of it; the engine and the binding spell its types without the namespace, as ChannelSpace.h does
 * for physics/'s constants.
 */
#include "azure2/inference/LevelParameterLayout.h"
#include "azure2/inference/ParameterIndexMap.h"
#include "azure2/inference/ParameterKind.h"
#include "azure2/inference/PosteriorSupport.h"

using azure2::inference::IsLevelSlot;
using azure2::inference::KindRow;
using azure2::inference::LevelParameterLayout;
using azure2::inference::ParamDesc;
using azure2::inference::ParamIndexMap;
using azure2::inference::ParamKind;
using azure2::inference::ParamKindRow;
using azure2::inference::ParamKindTable;
using azure2::inference::ParameterRecord;
using azure2::inference::PosteriorSupport;

#endif  // INFERENCETYPES_H
