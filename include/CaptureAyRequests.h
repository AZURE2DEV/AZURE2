#ifndef CAPTUREAYREQUESTS_H
#define CAPTUREAYREQUESTS_H

class CompoundSystem;
class WorkingCopy;
struct CompoundSetUpRequests;

/*!
 * The capture analyzing-power terms the data of a run need (cnuc_split.md section 11.5, decision 2): every term of an
 * analyzing-power data set or extrapolation -- its own and each leaf (component, sample nucleus) -- whose exit pair
 * is a photon pair, asked for in `requests` (CompoundSetUpRequests::AddCaptureAy).  Called by the data's set-up after
 * the data are filled (WorkingCopy::Fill or MakePoints); the compound's last part builds the terms into the pathway
 * table once.
 */
void RequestCaptureAyTerms(const CompoundSystem &compound, const WorkingCopy &data, CompoundSetUpRequests &requests);

#endif
