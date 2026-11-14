#ifndef CAPTUREAYREQUESTS_H
#define CAPTUREAYREQUESTS_H

class CompoundSystem;
class WorkingCopy;
struct CompoundSetUpRequests;

/*!
 * The capture analyzing-power terms the data of a run need (cnuc_split.md section 11.5, decision 2): every term of an
 * analyzing-power data set or extrapolation -- its own and each leaf (component, sample nucleus) -- whose exit pair
 * is a photon pair, requested on `compound` with CompoundSystem::RequestCaptureAy.  Called after the data are filled (WorkingCopy::Fill
 * or MakePoints) and before CompoundSystem::Initialize, which builds the terms into the pathway table once.
 */
void RequestCaptureAyTerms(CompoundSystem &compound, WorkingCopy &data);
/// The same terms asked for in `requests` (CompoundSystemSetUp.h), from the data's set-up.
void RequestCaptureAyTerms(const CompoundSystem &compound, const WorkingCopy &data, CompoundSetUpRequests &requests);

#endif
