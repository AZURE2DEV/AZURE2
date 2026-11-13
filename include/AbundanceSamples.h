#ifndef ABUNDANCESAMPLES_H
#define ABUNDANCESAMPLES_H

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "Constants.h"
#include "azure2/experiment/SetupSpecs.h"

class AZUREParams;

/// The isotopic abundance of one compound nucleus's target in one sample, as the project states it.
struct SampleAbundance {
  double start = 1.0;       ///< the value the project file gives (the fit's start value)
  double nominal = 1.0;     ///< the centre of the prior as written, else the value read (parameter_kinds.md section 4)
  double sigma = 0.0;       ///< prior width; 0 = no prior
  bool vary = false;
  bool priorReplaced = false;  ///< the prior covariance replaces this prior (docs/upgrade/objective.md D3)
  std::string isotope;  ///< format 3: the isotope as written, which names the parameter ("" = the old name)
  int z = 0, a = 0;     ///< the compound nucleus
};

/// A sample: one set of isotopic (or stoichiometric) abundances, one entry per compound nucleus.  The <abundances>
/// block holds one unnamed sample (the lines without a "<name>:" prefix) and any number of named ones; a segment
/// chooses its sample with "natural 1" (the unnamed, or the only, sample) or "natural <name>".  Each sample has its
/// own sum-to-one penalty width: the block's "<name>: sum <s>" line, else the config's abundanceSumSigma.
struct AbundanceSample {
  std::string name;                      ///< "" for the unnamed sample
  std::vector<SampleAbundance> entries;  ///< index nucleus-1
  double sumSigma = -1.0;                ///< < 0: use the config's abundanceSumSigma
  double sumCentre = 1.0;                ///< the centre of the sum-to-one prior: as written, else 1
};

/*!
 * The samples of the project's <abundances> block (the targets' isotopic abundances), as set up: their names, the
 * file's values, the priors and the sum-to-one widths.  The current abundances are a vector of values held by the
 * evaluation, one per (sample, nucleus) in the order of the abundance parameters (Ordinal).
 *
 * The one place of the abundance block in the parameter vectors: one abundance per (sample, compound nucleus), sample
 * by sample, right after the level block, in the fit's vector and in the observed vector alike (an abundance is the
 * same number in both bases).  AddParameters writes it into the fit's vector, AppendTo into an observed vector,
 * ReadFrom reads it back from either.  A new per-sample quantity of the parameter vectors is added here.
 */
class AbundanceSamples {
 public:
  /// The samples from their specs, for compound nuclei with these (Z, A) (index nucleus-1); -1 (with the reader's
  /// message, on \p messages) when the reader stopped.  A row for a (Z, A) no nucleus has is reported and ignored.
  int Set(std::ostream &messages, const azure2::experiment::AbundanceSpecs &specs,
          const std::vector<std::pair<int, int> > &nucleiZA);
  /// Did the project carry an <abundances> block?  Only then are the abundances parameters (one per nucleus and
  /// sample, after the level block), so a project without the block keeps its parameter layout.
  bool Present() const { return present_; }
  /// Number of samples (0 without the block).
  int NumSamples() const { return (int)samples_.size(); }
  /// Sample \p s, 1-based.
  const AbundanceSample &Sample(int s) const { return samples_[s - 1]; }
  /// The prior covariance replaces the prior of (sample, nucleus), both 1-based (ReplaceDiagonalPriors).
  void ReplacePrior(int sample, int nucleus) { samples_[sample - 1].entries[nucleus - 1].priorReplaced = true; }
  /// 1-based index of the sample with this name, 0 if there is none.
  int Find(const std::string &name) const;
  /// The sample "natural 1" means: the unnamed one if present, else the only one; 0 if ambiguous or none.
  int Default() const;
  /// Number of abundance parameters: samples x nuclei with the block, else 0.
  int NumParams() const { return present_ ? (int)samples_.size() * numNuclei_ : 0; }
  /// 1-based ordinal of the abundance parameter of (sample, nucleus), both 1-based: samples are blocks of nuclei.
  int Ordinal(int sample, int nucleus) const { return (sample - 1) * numNuclei_ + nucleus; }
  /// The (sample, nucleus) of a 1-based abundance ordinal.
  void Of(int ordinal, int &sample, int &nucleus) const;
  /// The parameter's name: "<sample>.isotope.<isotope>" (format 3), else "abundance_[<sample>_]<Z>_<A>".
  std::string ParamName(int sample, int nucleus) const;
  /// The sum-to-one width of a sample: its own "sum" line, else \p globalSumSigma.
  double SumSigma(int sample, double globalSumSigma) const;
  /// The file's values, in parameter order: the current values before any parameter vector is read.
  std::vector<double> StartValues() const;
  /// The current abundance of (sample, nucleus) in \p values (1 when there is no such sample).
  double Value(const std::vector<double> &values, int sample, int nucleus) const;
  /// The sum of the current abundances of one sample.
  double Sum(const std::vector<double> &values, int sample) const;
  /// The abundance penalty of chi-squared, summed over samples: every entry with a prior contributes
  /// ((w - w0)/sigma)^2, and each sample whose width s is > 0 contributes ((sum w) - c)^2 / s^2, c its centre.
  double Penalty(const std::vector<double> &values, double globalSumSigma) const;
  /// The sum-to-one terms alone (the sampler: the per-nucleus priors are its normalised log-priors).
  double SumPenalty(const std::vector<double> &values, double globalSumSigma) const;

  // -- the abundance block of the parameter vectors ----------------------------------------------------------------
  /// Writes the block into the fit's vector: one parameter per (sample, nucleus) at \p values, with its record, a
  /// lower limit 0 and fixed unless varied (a fixed abundance is still a parameter, so the layout does not depend on
  /// which ones vary).  Nothing without the block.
  void AddParameters(const std::vector<double> &values, AZUREParams &params) const;
  /// Writes the block into an observed vector that holds the level block: the abundances of the fit's vector
  /// \p fit it was made from, at the same positions.
  void AppendTo(const vector_r &fit, vector_r &observed) const;
  /// Reads the block of a vector (fit's or observed) whose level block has \p offset entries into \p values.
  void ReadFrom(const vector_r &vector, std::size_t offset, std::vector<double> &values) const;

 private:
  bool present_ = false;
  int numNuclei_ = 0;
  std::vector<AbundanceSample> samples_;
};

#endif
