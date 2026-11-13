// The reader of a data covariance file (CovarianceFile.h): the one implementation of its grammar.
#include "azure2/experiment/CovarianceFile.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace azure2 {
namespace experiment {

std::string CanonicalPath(const std::string &path) {
  std::error_code ec;
  std::filesystem::path p = std::filesystem::weakly_canonical(std::filesystem::path(path), ec);
  if (ec) return std::filesystem::path(path).lexically_normal().string();
  return p.string();
}

int CovarianceFile::FileIndex(const std::string &canonicalFile) const {
  for (size_t k = 0; k < canonical.size(); k++)
    if (canonical[k] == canonicalFile) return (int)k;
  return -1;
}

int CovarianceFile::RowOf(const std::string &canonicalFile, int record) const {
  const int k = FileIndex(canonicalFile);
  if (k < 0 || record < 1 || record > counts[k]) return -1;
  return offsets[k] + record - 1;
}

namespace {

// Whitespace-separated tokens with '#' comments; a token may be a "quoted path".
class Tokens {
 public:
  explicit Tokens(std::istream &in) : in_(in) {}
  bool Next(std::string &tok) {
    while (true) {
      while (pos_ < line_.size() && std::isspace((unsigned char)line_[pos_])) pos_++;
      if (pos_ < line_.size() && line_[pos_] != '#') break;
      if (!std::getline(in_, line_)) return false;
      lineNo_++;
      pos_ = 0;
    }
    if (line_[pos_] == '"') {
      size_t end = line_.find('"', pos_ + 1);
      if (end == std::string::npos) end = line_.size();
      tok = line_.substr(pos_ + 1, end - pos_ - 1);
      pos_ = std::min(line_.size(), end + 1);
      return true;
    }
    size_t end = pos_;
    while (end < line_.size() && !std::isspace((unsigned char)line_[end]) && line_[end] != '#') end++;
    tok = line_.substr(pos_, end - pos_);
    pos_ = end;
    return true;
  }
  bool Number(double &x) {
    std::string tok;
    if (!Next(tok)) return false;
    char *e = nullptr;
    x = std::strtod(tok.c_str(), &e);
    if (e == tok.c_str() || *e != '\0' || !std::isfinite(x)) {
      bad_ = tok;
      return false;
    }
    return true;
  }
  int Line() const { return lineNo_; }
  const std::string &Bad() const { return bad_; }

 private:
  std::istream &in_;
  std::string line_;
  size_t pos_ = 0;
  int lineNo_ = 0;
  std::string bad_;
};

// Data rows of a data file, counted as Segment::Fill counts them (DataLine: neither blank nor '#').
int CountRecords(std::istream &in) {
  std::string raw;
  int n = 0;
  while (std::getline(in, raw)) {
    // CR-only line ends are split as DataLine does.
    std::stringstream parts(raw);
    std::string piece;
    while (std::getline(parts, piece, '\r')) {
      std::string::size_type k = piece.find_first_not_of(" \t");
      if (k == std::string::npos || piece[k] == '#') continue;
      n++;
    }
  }
  return n;
}

}  // namespace

std::shared_ptr<CovarianceFile> ReadCovarianceFile(const std::string &path, const CovarianceFileAccess &access,
                                                   std::string &err) {
  auto resolve = [&](const std::string &name) { return access.resolve ? access.resolve(name) : name; };
  auto open = [&](const std::string &name) -> std::unique_ptr<std::istream> {
    if (access.open) return access.open(resolve(name));
    return std::unique_ptr<std::istream>(new std::ifstream(resolve(name)));
  };
  std::unique_ptr<std::istream> inStream = open(path);
  if (!inStream || !*inStream) {
    err = "cannot open the covariance file " + path;
    return nullptr;
  }
  std::istream &in = *inStream;
  auto f = std::make_shared<CovarianceFile>();
  f->path = path;
  Tokens t(in);
  std::string tok;
  auto where = [&]() { return path + " (line " + std::to_string(t.Line()) + ")"; };
  bool haveTerm = false;
  auto needNumbers = [&](const std::string &what, size_t count, std::vector<double> &v) -> bool {
    v.resize(count);
    for (size_t k = 0; k < count; k++)
      if (!t.Number(v[k])) {
        err = what + " in " + where() + ": expected " + std::to_string(count) + " numbers, " +
              (t.Bad().empty() ? std::string("the file ended") : "found \"" + t.Bad() + "\"") + " after " +
              std::to_string(k);
        return false;
      }
    return true;
  };
  auto startMatrix = [&]() -> bool {
    if (f->N == 0) {
      err = where() + ": the matrix terms must follow the \"file <data file> <rows>\" lines";
      return false;
    }
    if (f->diag.empty()) f->diag.assign(f->N, 0.0);
    haveTerm = true;
    return true;
  };
  // Symmetry of an input square block a (m x m): |a_ij - a_ji| <= 1e-6 sqrt(|a_ii a_jj|).
  auto checkSymmetric = [&](const std::vector<double> &a, size_t off, int m, const std::string &what) -> bool {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < i; j++) {
        const double x = a[off + (size_t)i * m + j], y = a[off + (size_t)j * m + i];
        const double scale = std::max({std::sqrt(std::fabs(a[off + (size_t)i * m + i] * a[off + (size_t)j * m + j])),
                                       std::fabs(x), std::fabs(y)});
        if (std::fabs(x - y) > 1.0e-6 * scale) {
          std::ostringstream msg;
          msg.precision(10);
          msg << what << " in " << where() << " is not symmetric: element (" << i + 1 << "," << j + 1 << ") " << x
              << " vs (" << j + 1 << "," << i + 1 << ") " << y;
          err = msg.str();
          return false;
        }
      }
    return true;
  };
  std::vector<double> v;
  while (t.Next(tok)) {
    if (tok == "file") {
      if (haveTerm) {
        err = where() + ": every \"file\" line must come before the matrix terms";
        return nullptr;
      }
      std::string name;
      double count = 0;
      if (!t.Next(name) || !t.Number(count) || count < 1 || count != std::floor(count)) {
        err = where() + ": \"file\" needs a data file name and its number of data rows";
        return nullptr;
      }
      const std::string canon = CanonicalPath(resolve(name));
      if (f->FileIndex(canon) >= 0) {
        err = where() + ": data file " + name + " is listed twice";
        return nullptr;
      }
      std::unique_ptr<std::istream> data = open(name);
      const int actual = data && *data ? CountRecords(*data) : -1;
      if (actual < 0) {
        err = where() + ": cannot open data file " + name + " (paths are relative to the project's directory)";
        return nullptr;
      }
      if (actual != (int)count) {
        err = where() + ": data file " + name + " has " + std::to_string(actual) + " data rows, the header says " +
              std::to_string((int)count);
        return nullptr;
      }
      f->files.push_back(name);
      f->canonical.push_back(canon);
      f->counts.push_back((int)count);
      f->offsets.push_back(f->N);
      f->N += (int)count;
    } else if (tok == "size") {
      double n = 0;
      if (!t.Number(n) || (int)n != f->N) {
        err = where() + ": \"size\" " + std::to_string((long long)n) + " does not match the " + std::to_string(f->N) +
              " rows of the file lines";
        return nullptr;
      }
    } else if (tok == "block") {
      if (!startMatrix()) return nullptr;
      double m = 0;
      if (!t.Number(m) || m < 1 || m != std::floor(m) || m > f->N) {
        err = where() + ": \"block\" needs its size n (1.." + std::to_string(f->N) + ")";
        return nullptr;
      }
      const int n = (int)m;
      std::vector<double> idx;
      if (!needNumbers("block row indices", n, idx)) return nullptr;
      std::vector<int> rows(n);
      std::vector<char> seen(f->N, 0);
      for (int k = 0; k < n; k++) {
        if (idx[k] != std::floor(idx[k]) || idx[k] < 1 || idx[k] > f->N) {
          err = where() + ": block row index " + std::to_string(idx[k]) + " is not a row 1.." + std::to_string(f->N);
          return nullptr;
        }
        rows[k] = (int)idx[k] - 1;
        if (seen[rows[k]]++) {
          err = where() + ": block row " + std::to_string(rows[k] + 1) + " appears twice";
          return nullptr;
        }
      }
      if (!needNumbers("block values", (size_t)n * n, v)) return nullptr;
      if (!checkSymmetric(v, 0, n, "block")) return nullptr;
      CovarianceBlock b;
      b.rows = rows;
      b.values.resize((size_t)n * n);
      for (int a = 0; a < n; a++)
        for (int c = 0; c < n; c++) b.values[(size_t)a * n + c] = 0.5 * (v[(size_t)a * n + c] + v[(size_t)c * n + a]);
      f->blocks.push_back(std::move(b));
    } else if (tok == "diagonal") {
      if (!startMatrix()) return nullptr;
      if (!needNumbers("diagonal", f->N, v)) return nullptr;
      for (int a = 0; a < f->N; a++) f->diag[a] += v[a];
    } else if (tok == "vector") {
      if (!startMatrix()) return nullptr;
      if (!needNumbers("vector", f->N, v)) return nullptr;
      f->vectors.push_back(v);
    } else if (tok == "relvector") {
      if (!startMatrix()) return nullptr;
      if (!needNumbers("relvector", f->N, v)) return nullptr;
      f->rel.push_back(v);
    } else if (tok == "matrix") {
      if (!startMatrix()) return nullptr;
      if (!needNumbers("matrix", (size_t)f->N * f->N, v)) return nullptr;
      if (!checkSymmetric(v, 0, f->N, "matrix")) return nullptr;
      CovarianceBlock b;  // a dense block on every row
      b.rows.resize(f->N);
      for (int a = 0; a < f->N; a++) b.rows[a] = a;
      b.values.resize((size_t)f->N * f->N);
      for (int a = 0; a < f->N; a++)
        for (int c = 0; c < f->N; c++)
          b.values[(size_t)a * f->N + c] = 0.5 * (v[(size_t)a * f->N + c] + v[(size_t)c * f->N + a]);
      f->blocks.push_back(std::move(b));
    } else {
      err = where() + ": unknown term \"" + tok + "\" (file, size, block, diagonal, vector, relvector, matrix)";
      return nullptr;
    }
  }
  if (f->N == 0) {
    err = path + ": no \"file <data file> <rows>\" lines";
    return nullptr;
  }
  if (!haveTerm) {
    err = path + ": no matrix terms (block, diagonal, vector, relvector, matrix)";
    return nullptr;
  }
  f->diagB = f->diag;
  for (const CovarianceBlock &b : f->blocks)
    for (size_t a = 0; a < b.rows.size(); a++) f->diagB[b.rows[a]] += b.values[a * b.rows.size() + a];
  return f;
}

double CovarianceFile::Variance(int row, double model) const {
  double c = diagB[row];
  for (const std::vector<double> &v : vectors) c += v[row] * v[row];
  for (const std::vector<double> &fr : rel) c += fr[row] * fr[row] * model * model;
  return c;
}

}  // namespace experiment
}  // namespace azure2
