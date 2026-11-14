#ifndef COMPOUNDSYSTEMPATHWAYS_H
#define COMPOUNDSYSTEMPATHWAYS_H

#include <vector>

#include "PathwayTable.h"

class CompoundSystem;

/// The external-capture final states of a CompoundSystem after Fill, in J-group and level order: what physics/
/// reaction_pathways needs of each (2J, parity, photon partition, multipoles, channels) and its J-group and level.
std::vector<azure2::PathwayFinalState> ECFinalStates(const CompoundSystem &compound);

#endif
