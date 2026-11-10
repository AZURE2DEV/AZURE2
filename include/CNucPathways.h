#ifndef CNUCPATHWAYS_H
#define CNUCPATHWAYS_H

#include <vector>

#include "PathwayTable.h"

class CNuc;

/// The external-capture final states of a CNuc after Fill, in J-group and level order: what physics/
/// reaction_pathways needs of each (2J, parity, photon partition, multipoles, channels) and its J-group and level.
std::vector<azure2::PathwayFinalState> ECFinalStates(CNuc &compound);

#endif
