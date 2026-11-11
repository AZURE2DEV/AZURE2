// The data covariance (segment line "covariance <file>"): reading the covariance files, one group per file with
// its Cholesky factor, the chi2 correction r^T M^-1 r - sum r_i^2 with the per-point seeds, the whitened residual
// rows (and Jacobian columns) of the least-squares interface, and the summary file covariance.out.
#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Config.h"
#include "DataCovariance.h"
#include "EData.h"
#include "ModelPass.h"

std::shared_ptr<DataCovariance::File> EData::CovarianceFileFor(const std::string &path, const Config &configure) {
  const std::string key = DataCovariance::CanonicalPath(path);
  auto it = covFiles_.find(key);
  if (it != covFiles_.end()) return it->second;
  std::string err;
  std::shared_ptr<DataCovariance::File> f = DataCovariance::Read(path, err, configure.cacheInputFiles);
  if (!f) {
    configure.outStream << "ERROR: Covariance file: " << err << "." << std::endl;
    covarianceError_ = true;
    return nullptr;
  }
  configure.outStream << "Covariance file " << path << ": " << f->N << " rows over " << f->files.size()
                      << " data file(s)" << (f->rel.empty() ? std::string() : ", " + std::to_string(f->rel.size()) +
                                                                                   " relvector term(s)")
                      << "." << std::endl;
  covFiles_[key] = f;
  return f;
}

/*!
 * One group per covariance file: its segments' points in the file's row order, the
 * absolute part of their submatrix, the relvector fractions on their rows, and the
 * Cholesky factor of M = D^-1 C D^-1.  Without relvector terms the factor is final here;
 * with them it is checked here with the data in place of the model and frozen at the
 * first evaluation (DataCovariance::Group::EnsureFactored).
 */
int EData::ResolveCovarianceGroups(const Config &configure) {
  covGroups_.clear();
  std::map<std::string, int> groupOf;
  struct Item { int row, seg, pid; };
  std::vector<std::vector<Item>> items;
  std::vector<std::shared_ptr<DataCovariance::File>> files;
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (seg->GetCovarianceFile().empty()) continue;
    const std::string key = DataCovariance::CanonicalPath(seg->GetCovarianceFile());
    auto fit = covFiles_.find(key);
    if (fit == covFiles_.end()) return -1;  // Fill reported it
    auto git = groupOf.find(key);
    int g;
    if (git == groupOf.end()) {
      g = (int)items.size();
      groupOf[key] = g;
      items.emplace_back();
      files.push_back(fit->second);
      auto grp = std::make_shared<DataCovariance::Group>();
      grp->path = seg->GetCovarianceFile();
      covGroups_.push_back(grp);
    } else {
      g = git->second;
    }
    seg->SetCovarianceGroup(g);
    covGroups_[g]->numSegments++;
    const std::string canonicalData = DataCovariance::CanonicalPath(seg->GetDataFile());
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      const int row = fit->second->RowOf(canonicalData, seg->FileRecord(pid));
      if (row < 0) return -1;  // cannot happen: Fill checked the file and the record count
      items[g].push_back({row, s, pid});
    }
  }
  if (covGroups_.empty()) {
    covFiles_.clear();
    return 0;
  }
  const std::string modelFile = configure.outputdir + "covarianceModel.out";
  bool anyRelative = false;
  for (size_t g = 0; g < covGroups_.size(); g++) {
    DataCovariance::Group &grp = *covGroups_[g];
    const DataCovariance::File &f = *files[g];
    std::vector<Item> &it = items[g];
    std::stable_sort(it.begin(), it.end(), [](const Item &a, const Item &b) { return a.row < b.row; });
    for (size_t k = 1; k < it.size(); k++)
      if (it[k].row == it[k - 1].row) {
        configure.outStream << "ERROR: Row " << it[k].row + 1 << " of covariance file " << grp.path
                            << " is fitted twice, in segments #" << this->GetSegment(it[k - 1].seg)->GetSegmentKey()
                            << " and #" << this->GetSegment(it[k].seg)->GetSegmentKey()
                            << " (overlapping energy or angle windows on one data file)." << std::endl;
        return -1;
      }
    const int n = (int)it.size();
    grp.n = n;
    grp.N = f.N;
    grp.members.resize(n);
    grp.rows.resize(n);
    grp.sigma.resize(n);
    grp.data.resize(n);
    grp.label.resize(n);
    for (int a = 0; a < n; a++) {
      grp.members[a] = {it[a].seg, it[a].pid};
      grp.rows[a] = it[a].row;
      ESegment *seg = this->GetSegment(it[a].seg);
      EPoint *pt = seg->GetPoint(it[a].pid + 1);
      grp.sigma[a] = seg->LabError(it[a].pid);
      grp.data[a] = seg->LabValue(it[a].pid);
      std::ostringstream lab;
      lab << seg->GetSegmentKey() << " " << pt->GetLabEnergy() << " " << pt->GetLabAngle();
      grp.label[a] = lab.str();
    }
    grp.modelFile = modelFile;
    // B (diagonal and blocks) factored by connected components; it must be positive
    // definite by itself on the fitted rows.
    const int bad = grp.SetUp(f);
    if (bad >= 0) {
      configure.outStream << "ERROR: Covariance file " << grp.path << ": its diagonal and block terms are not positive"
                          << " definite on the fitted rows (at file row " << grp.rows[bad] + 1
                          << "). Every fitted row needs a variance of its own from a diagonal or block term;"
                          << " vector and relvector terms alone are not enough." << std::endl;
      return -1;
    }
    if (grp.HasRelative())
      anyRelative = true;  // the low-rank part waits for the model of the first evaluation
    else
      grp.Factor(grp.data);
    configure.outStream << "Covariance " << grp.path << ": " << n << " of " << f.N << " rows fitted, in "
                        << grp.numSegments << " segment(s); " << grp.NumComponents()
                        << " diagonal/block component(s) (largest " << grp.LargestComponent() << " rows), "
                        << grp.absolute.size() << " vector term(s)"
                        << (grp.HasRelative() ? ", " + std::to_string(grp.rel.size()) +
                                                    " relvector term(s) on the model frozen at the first evaluation"
                                              : std::string())
                        << "." << std::endl;
    if (f.numSigmaMismatch > 0)
      configure.outStream << "WARNING: " << f.numSigmaMismatch << " point(s) of " << grp.path
                          << " have an error column more than 1 % from sqrt(C_ii) (worst "
                          << 100.0 * f.worstSigmaMismatch << " %; first: " << f.firstSigmaMismatch
                          << "); sqrt(C_ii) is used." << std::endl;
  }
  if (anyRelative) std::ofstream(modelFile.c_str(), std::ios::trunc);
  covFiles_.clear();
  return 0;
}

double EData::CovarianceChi2Correction(azure2::Results &results, std::vector<double> *segChis) {
  double correction = 0.0;
  std::vector<double> r, z;
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    r.assign(g.n, 0.0);
    double diag = 0.0;
    for (int a = 0; a < g.n; a++) {
      ESegment *seg = this->GetSegment(g.members[a].first);
      EPoint *pt = seg->GetPoint(g.members[a].second + 1);
      r[a] = azure2::PointResidual(*seg, g.members[a].second, results, seg->GetNorm()).Pull();
      diag += r[a] * r[a];
    }
    g.EnsureFactored(r);
    z = r;
    g.Whiten(z.data());
    double q = 0.0;
    for (int a = 0; a < g.n; a++) q += z[a] * z[a];
    g.WhitenTranspose(z.data());
    correction += q - diag;
    for (int a = 0; a < g.n; a++)
      this->GetSegment(g.members[a].first)->GetPoint(g.members[a].second + 1)->SetCovarianceSeed(results, z[a]);
    // Per-segment chi2: the exact partition sum_i r_i (M^-1 r)_i of the group's chi2.
    std::map<int, double> part;
    for (int a = 0; a < g.n; a++) part[g.members[a].first] += r[a] * z[a];
    for (const auto &ps : part) {
      results.objective.segmentChi2[ps.first] = ps.second;
      if (segChis && ps.first < (int)segChis->size()) (*segChis)[ps.first] = ps.second;
    }
  }
  return correction;
}

void EData::WhitenResidualRows(std::vector<double> &residuals, std::vector<double> *jacobian, int nCols) {
  if (covGroups_.empty()) return;
  std::vector<int> rowStart(this->NumSegments() + 2, 0);
  int row = 0;
  for (int s = 1; s <= this->NumSegments(); s++) {
    rowStart[s] = row;
    row += this->GetSegment(s)->NumPoints();
  }
  if (row > (int)residuals.size()) return;  // not the data-row layout
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    std::vector<int> rows(g.n);
    for (int a = 0; a < g.n; a++) rows[a] = rowStart[g.members[a].first] + g.members[a].second;
    std::vector<double> y(g.n);
    for (int a = 0; a < g.n; a++) y[a] = residuals[rows[a]];
    g.EnsureFactored(y);  // the first evaluation of a run with relvector terms freezes the model here
    g.Whiten(y.data());
    for (int a = 0; a < g.n; a++) residuals[rows[a]] = y[a];
    if (!jacobian || nCols <= 0) continue;
    std::vector<double> &J = *jacobian;
#pragma omp parallel for schedule(dynamic, 4)
    for (int c = 0; c < nCols; c++) {
      std::vector<double> col(g.n);
      bool any = false;
      for (int a = 0; a < g.n; a++) any |= (col[a] = J[(size_t)rows[a] * nCols + c]) != 0.0;
      if (!any) continue;
      g.Whiten(col.data());
      for (int a = 0; a < g.n; a++) J[(size_t)rows[a] * nCols + c] = col[a];
    }
  }
}

/*!
 * covariance.out: per covariance group, the fitted rows, the chi2 r^T M^-1 r of the
 * current fit values and, for comparison, the sum of r_i^2 (the chi2 the same points
 * would have without their correlations).
 */
void EData::WriteCovarianceSummary(const Config &configure, const azure2::Results &results) {
  std::ofstream out((configure.outputdir + "covariance.out").c_str());
  if (!out) return;
  out << "# file, fitted rows, rows, segments, chi2 (r^T M^-1 r), sum of r_i^2 (correlations ignored), relvector terms"
      << std::endl;
  out.precision(12);
  out << std::scientific;
  std::vector<double> r;
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) {
    DataCovariance::Group &g = *gp;
    r.assign(g.n, 0.0);
    double diag = 0.0;
    for (int a = 0; a < g.n; a++) {
      ESegment *seg = this->GetSegment(g.members[a].first);
      EPoint *pt = seg->GetPoint(g.members[a].second + 1);
      r[a] = azure2::PointResidual(*seg, g.members[a].second, results, seg->GetNorm()).Pull();
      diag += r[a] * r[a];
    }
    g.EnsureFactored(r);
    g.Whiten(r.data());
    double q = 0.0;
    for (int a = 0; a < g.n; a++) q += r[a] * r[a];
    out << g.path << ", " << g.n << ", " << g.N << ", " << g.numSegments << ", " << q << ", " << diag << ", "
        << g.rel.size() << std::endl;
  }
}

bool EData::HasRelativeCovariance() const {
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_)
    if (gp->HasRelative()) return true;
  return false;
}

void EData::RefreezeCovariance() {
  for (const std::shared_ptr<DataCovariance::Group> &gp : covGroups_) gp->Refreeze();
}

void EData::PrepareCovarianceSeeds(CNuc *theCNuc, const Config &configure, azure2::Results &results) {
  if (covGroups_.empty()) return;
  std::vector<std::pair<ESegment *, int>> work;
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *seg = this->GetSegment(s);
    if (seg->GetCovarianceGroup() < 0) continue;
    for (int pid = 0; pid < seg->NumPoints(); pid++) work.push_back({seg, pid});
  }
  // Each point on its own, without point mapping (as the gradient's own point loops): a
  // mapped point would otherwise read a host value that may not be calculated yet.
#pragma omp parallel for schedule(dynamic) if (!omp_in_parallel())
  for (int k = 0; k < (int)work.size(); k++) {
    ScopedPointMappingOff noPointMapping;
    azure2::EvaluateModelPoint(*this, *work[k].first, work[k].second, theCNuc, configure, results, azure2::LeafValues::Calculate);
  }
  this->CovarianceChi2Correction(results);
}
