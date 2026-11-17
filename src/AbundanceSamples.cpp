// The samples of the <abundances> block and the abundance block of the parameter vectors (AbundanceSamples.h).
#include "AbundanceSamples.h"

#include <cmath>

#include "AZUREParams.h"
#include "Minuit2/MnUserParameters.h"

/*!
 * The abundances of the samples (AbundancesSpec): one row per compound nucleus and sample, or a sample's sum-prior
 * width.  A nucleus without a row keeps 1.0, fixed.  A row for a (Z, A) the pairs do not form is reported and
 * ignored.  Without the block no abundance becomes a parameter, so every project written before the block existed
 * keeps its parameter layout, param.par and <parameterSettings> indices.
 */
int AbundanceSamples::Set(std::ostream &messages, const azure2::experiment::AbundanceSpecs &specs,
                          const std::vector<std::pair<int, int> > &nucleiZA) {
  present_ = specs.abundances.present;
  numNuclei_ = (int)nucleiZA.size();
  samples_.clear();
  if (!present_) return 0;
  // A sample is created the first time a row sets one of its values; the unnamed sample is "".
  auto sampleOf = [&](const std::string &name) -> AbundanceSample & {
    for (unsigned int s = 0; s < samples_.size(); s++)
      if (samples_[s].name == name) return samples_[s];
    AbundanceSample sample;
    sample.name = name;
    sample.entries.assign(nucleiZA.size(), SampleAbundance());
    for (unsigned int i = 0; i < nucleiZA.size(); i++) {
      sample.entries[i].z = nucleiZA[i].first;
      sample.entries[i].a = nucleiZA[i].second;
    }
    samples_.push_back(sample);
    return samples_.back();
  };
  for (const azure2::experiment::AbundanceRowSpec &row : specs.abundances.rows) {
    if (row.isSum) {
      sampleOf(row.sample).sumSigma = row.sumSigma;
      sampleOf(row.sample).sumCentre = row.sumCentre.Or(1.0);
      continue;
    }
    bool found = false;
    for (unsigned int i = 0; i < nucleiZA.size(); i++) {
      if (nucleiZA[i].first != row.z || nucleiZA[i].second != row.a) continue;
      SampleAbundance &e = sampleOf(row.sample).entries[i];
      e.start = row.value;
      e.nominal = row.centre.Or(row.value);  // the prior's centre as written (parameter_kinds.md section 4)
      e.vary = row.vary;
      e.sigma = row.sigma;
      e.isotope = row.isotope;
      found = true;
    }
    if (!found)
      messages << "WARNING: <abundances> names Z=" << row.z << " A=" << row.a
                          << ", which no pair of the project forms; ignored." << std::endl;
  }
  if (specs.failed) {
    if (!specs.error.empty()) messages << specs.error << std::endl;
    return -1;
  }
  // An empty block still means "abundances are parameters": one unnamed sample at the defaults.
  if (samples_.empty()) sampleOf("");
  return 0;
}

int AbundanceSamples::Find(const std::string &name) const {
  for (unsigned int s = 0; s < samples_.size(); s++)
    if (samples_[s].name == name) return (int)s + 1;
  return 0;
}

int AbundanceSamples::Default() const {
  int unnamed = this->Find("");
  if (unnamed > 0) return unnamed;
  return samples_.size() == 1 ? 1 : 0;
}

void AbundanceSamples::Of(int ordinal, int &sample, int &nucleus) const {
  sample = (ordinal - 1) / numNuclei_ + 1;
  nucleus = (ordinal - 1) % numNuclei_ + 1;
}

std::string AbundanceSamples::ParamName(int sample, int nucleus) const {
  const SampleAbundance &e = samples_[sample - 1].entries[nucleus - 1];
  return azure2::experiment::AbundanceParameterName(samples_[sample - 1].name, e.isotope, e.z, e.a);
}

double AbundanceSamples::SumSigma(int sample, double globalSumSigma) const {
  const double own = samples_[sample - 1].sumSigma;
  return own >= 0.0 ? own : globalSumSigma;
}

std::vector<double> AbundanceSamples::StartValues() const {
  std::vector<double> values;
  if (!present_) return values;
  for (const AbundanceSample &s : samples_)
    for (const SampleAbundance &e : s.entries) values.push_back(e.start);
  return values;
}

double AbundanceSamples::Value(const std::vector<double> &values, int sample, int nucleus) const {
  if (sample < 1 || sample > (int)samples_.size()) return 1.0;
  return values[Ordinal(sample, nucleus) - 1];
}

void AbundanceSamples::AddParameters(const std::vector<double> &values, AZUREParams &params) const {
  if (!present_) return;
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  for (unsigned int s = 0; s < samples_.size(); s++) {
    for (int k = 0; k < numNuclei_; k++) {
      const SampleAbundance &e = samples_[s].entries[k];
      const double w = values[Ordinal((int)s + 1, k + 1) - 1];
      std::string name = this->ParamName((int)s + 1, k + 1);
      params.Add(name, w, e.sigma > 0.0 ? e.sigma : 0.1 * w + 1e-3,
                 ParameterRecord::Abundance((int)s + 1, samples_[s].name, k + 1, e.z, e.a));
      p.SetLowerLimit(name, 0.0);
      if (!e.vary) p.Fix(name);
    }
  }
}

void AbundanceSamples::AppendTo(const vector_r &fit, vector_r &observed) const {
  const std::size_t offset = observed.size();
  for (int k = 0; k < NumParams(); k++) observed.push_back(fit[offset + (std::size_t)k]);
}

void AbundanceSamples::ReadFrom(const vector_r &vector, std::size_t offset, std::vector<double> &values) const {
  values.resize((std::size_t)NumParams());
  for (int k = 0; k < NumParams(); k++) values[(std::size_t)k] = vector[offset + (std::size_t)k];
}
