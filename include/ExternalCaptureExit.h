#ifndef EXTERNALCAPTUREEXIT_H
#define EXTERNALCAPTUREEXIT_H

class CNuc;
namespace azure2 {
struct ExitRow;
}

/*!
 * The exit to photon pair `exitPair` among the pathways from pair `entrancePair` (both CNuc pair numbers, 1-based) in
 * the pathway table, for the external-capture amplitudes of data segment `segmentKey`.  The exits of an entrance skip
 * every pair it has no pathway to (a beta pair, a pair of another compound nucleus), so it is found by the pair's
 * number (EntrancePathways::ExitOf), never taken at that number as a position.  Throws std::runtime_error naming
 * the data segment and the two pair keys when the entrance has no pathway to that pair.  ExternalCaptureSetup.cpp.
 */
const azure2::ExitRow &ExternalCaptureExit(const CNuc &compound, int entrancePair, int exitPair, int segmentKey);

#endif
