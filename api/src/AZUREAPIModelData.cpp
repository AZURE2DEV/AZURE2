// AZUREAPI: the model values of the <modelData> data (AZUREAPI.h ModelDataRWA; ModelData.h).
#include "AZUREAPI.h"

#include "CNuc.h"
#include "Config.h"
#include "ModelData.h"

std::vector<AZUREAPI::ModelDatumValue> AZUREAPI::ModelDataRWA(const vector_r &rwaParams) const {
  vector_r full = all_rwa_;
  int k = 0;
  for (int i = 0; i < (int)all_rwa_.size(); ++i)
    if (!fixed_[i]) full[i] = rwaParams[k++];
  CNuc *lc = compound();
  lc->FillCompoundFromParams(full);
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());
  std::vector<ModelDatumValue> out;
  for (const ModelData::Value &v : ModelData::SectionValues(lc, configure(), lc->GetModelData())) {
    const azure2::experiment::ModelDatum &d = *v.datum;
    ModelDatumValue r;
    r.name = d.name;
    r.quantity = d.quantityName;
    r.unit = d.unit;
    r.source = d.source;
    r.model = v.model / d.unitFactor;
    r.value = d.value / d.unitFactor;
    r.error = d.error / d.unitFactor;
    r.residual = v.residual;
    r.hasImag = d.hasImag;
    if (d.hasImag) {
      r.modelImag = v.modelImag / d.unitFactor;
      r.imag = d.imag / d.unitFactor;
      r.imagError = d.imagError / d.unitFactor;
      r.residualImag = v.residualImag;
    }
    r.chi2 = v.chi2;
    out.push_back(r);
  }
  return out;
}
