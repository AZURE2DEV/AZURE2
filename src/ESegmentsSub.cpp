#include "ESegmentsSub.h"
#include "EData.h"
#include "ESegment.h"
#include "CNuc.h"
#include "Config.h"
#include <iostream>

ESegmentsSub::ESegmentsSub(int entranceKey, int exitKey)
    : entranceKey_(entranceKey), exitKey_(exitKey) {
}

ESegmentsSub::~ESegmentsSub() {
}

double ESegmentsSub::CalculateTheoretical(int pointIndex, CNuc* cnuc, const Config& configure, EData* edata) const {
    if (!edata || !cnuc) {
        return 0.0;
    }
    
    // Component segments calculate theoretical cross section for different entrance/exit pairs
    // using the same energy/angle as the main segment. The calculation should be performed
    // by the main R-Matrix calculation engine with the component's entrance/exit keys.
    
    // Find any segment to get energy/angle reference point - we need the energy and angle
    // but will calculate using our own entrance/exit keys
    if (edata->NumSegments() == 0) {
        return 0.0;
    }
    
    ESegment* firstSegment = edata->GetSegment(1);
    if (!firstSegment || pointIndex >= firstSegment->NumPoints()) {
        return 0.0;
    }
    
    EPoint* referencePoint = firstSegment->GetPoint(pointIndex + 1);
    if (!referencePoint) {
        return 0.0;
    }
    
    // Create a temporary point with our entrance/exit keys and the reference energy/angle
    // Use the detailed constructor: EPoint(energy, angle, entranceKey, exitKey, differential, phase, angDist, j, l, maxLOrder)
    double energy = referencePoint->GetCMEnergy();
    double angle = referencePoint->GetCMAngle();
    bool differential = referencePoint->IsDifferential();
    bool phase = referencePoint->IsPhase();
    bool angDist = referencePoint->IsAngularDist();
    double j = referencePoint->GetJ();
    int l = referencePoint->GetL();
    int maxLOrder = referencePoint->GetMaxLOrder();
    
    EPoint tempPoint(energy, angle, entranceKey_, exitKey_, differential, phase, angDist, j, l, maxLOrder);
    
    try {
        // Calculate the theoretical cross section for this component
        tempPoint.Calculate(cnuc, configure);
        return tempPoint.GetFitCrossSection();
    } catch (...) {
        return 0.0;
    }
}