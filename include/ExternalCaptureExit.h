#ifndef EXTERNALCAPTUREEXIT_H
#define EXTERNALCAPTUREEXIT_H

class CNuc;

/*!
 * The position (1-based) of the exit to photon pair `exitPair` among the pathways from pair `entrancePair` (both
 * CNuc pair numbers, 1-based), for the external-capture amplitudes of data segment `segmentKey`.  The exits of an
 * entrance skip every pair it has no pathway to (a beta pair, a pair of another compound nucleus), so a pair number
 * is not a position.  Throws std::runtime_error naming the data segment and the two pair keys when the entrance has
 * no pathway to that pair.  ExternalCaptureSetup.cpp.
 */
int ExternalCaptureExit(const CNuc &compound, int entrancePair, int exitPair, int segmentKey);

#endif
