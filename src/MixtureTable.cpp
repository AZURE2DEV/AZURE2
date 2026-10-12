#include "MixtureTable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "PathToken.h"

std::shared_ptr<const MixtureTable> MixtureTable::Load(const std::string &path, std::string &error) {
  std::ifstream in(path.c_str());
  if (!in) {
    error = "cannot open \"" + path + "\"";
    return std::shared_ptr<const MixtureTable>();
  }
  std::shared_ptr<MixtureTable> table(new MixtureTable());
  table->file_ = path;
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    lineNo++;
    size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream row(line);
    double e, s;
    if (!(row >> e)) continue;  // blank or comment-only line
    if (!(row >> s) || !std::isfinite(e) || !std::isfinite(s)) {
      std::ostringstream msg;
      msg << "\"" << path << "\" line " << lineNo << ": expected  E_lab[MeV]  sigma[b]";
      error = msg.str();
      return std::shared_ptr<const MixtureTable>();
    }
    if (!table->e_.empty() && e <= table->e_.back()) {
      std::ostringstream msg;
      msg << "\"" << path << "\" line " << lineNo << ": energies must increase strictly (" << e
          << " after " << table->e_.back() << ")";
      error = msg.str();
      return std::shared_ptr<const MixtureTable>();
    }
    table->e_.push_back(e);
    table->sigma_.push_back(s);
  }
  if (table->e_.size() < 2) {
    error = "\"" + path + "\" holds fewer than two rows";
    return std::shared_ptr<const MixtureTable>();
  }
  return table;
}

double MixtureTable::Interpolate(double eLab) const {
  if (e_.empty()) return 0.0;
  if (eLab <= e_.front()) return sigma_.front();
  if (eLab >= e_.back()) return sigma_.back();
  std::vector<double>::const_iterator hi = std::upper_bound(e_.begin(), e_.end(), eLab);
  size_t j = hi - e_.begin();  // e_[j-1] < eLab <= e_[j]
  double frac = (eLab - e_[j - 1]) / (e_[j] - e_[j - 1]);
  return sigma_[j - 1] + frac * (sigma_[j] - sigma_[j - 1]);
}

bool MixtureTable::Covers(double lo, double hi) const {
  return !e_.empty() && lo >= e_.front() && hi <= e_.back();
}

bool MixtureSpec::Read(std::istream &stream) {
  MixtureSpec spec;
  int n = 0;
  if (!(stream >> spec.selfWeight >> n) || n < 0 || !std::isfinite(spec.selfWeight)) return false;
  for (int i = 0; i < n; i++) {
    std::string file;
    double weight;
    if (!ReadPathToken(stream, file) || !(stream >> weight) || !std::isfinite(weight)) return false;
    spec.tables.push_back(std::make_pair(file, weight));
  }
  spec.has = true;
  *this = spec;
  return true;
}
