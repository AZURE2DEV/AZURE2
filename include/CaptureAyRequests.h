#ifndef CAPTUREAYREQUESTS_H
#define CAPTUREAYREQUESTS_H

class CNuc;
class EData;

/*!
 * The capture analyzing-power terms the data of a run need (cnuc_split.md section 11.5, decision 2): every term of an
 * analyzing-power data set or extrapolation -- its own and each leaf (component, sample nucleus) -- whose exit pair
 * is a photon pair, requested on `compound` with CNuc::RequestCaptureAy.  Called after the data are filled (EData::Fill
 * or MakePoints) and before CNuc::Initialize, which builds the terms into the pathway table once.
 */
void RequestCaptureAyTerms(CNuc &compound, EData &data);

#endif
