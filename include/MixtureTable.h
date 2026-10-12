#ifndef MIXTURETABLE_H
#define MIXTURETABLE_H

#include <istream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/// A tabulated cross section sigma(E_lab) that a "mixture" segment adds to its model.

/*!
 * Natural-abundance and compound samples give a measured cross section that is the
 * abundance-weighted sum over the isotopes present, each with its own compound
 * nucleus.  One AZURE2 project describes one compound nucleus, so the other
 * isotopes enter a segment as fixed tables, read from a file with two columns,
 *
 *   E_lab [MeV]   sigma [b]      (# comments and blank lines ignored)
 *
 * in strictly increasing energy.  The energy is the projectile's lab energy,
 * the same axis the data file uses; the cross section is in the frame and units
 * AZURE2 compares the segment in (barns for angle-integrated and total cross
 * sections, c.m. b/sr for a differential one).  Values are interpolated
 * linearly; outside the table they are clamped to the end values, and whether
 * that is allowed is the caller's decision (see ESegment::ComputeMixtureExternal).
 */
class MixtureTable {
 public:
  /*!
   * Reads the table; on failure returns a null pointer and describes the
   * problem in \p error.
   */
  static std::shared_ptr<const MixtureTable> Load(const std::string &path, std::string &error);
  /// Linear interpolation at a lab energy, clamped to the end values outside the table.
  double Interpolate(double eLab) const;
  /// Does the table span [lo, hi]?
  bool Covers(double lo, double hi) const;
  double MinEnergy() const { return e_.empty() ? 0.0 : e_.front(); }
  double MaxEnergy() const { return e_.empty() ? 0.0 : e_.back(); }
  const std::string &File() const { return file_; }
  size_t Size() const { return e_.size(); }

 private:
  std::string file_;
  std::vector<double> e_;
  std::vector<double> sigma_;
};

/*!
 * The "mixture" keyword block of a segment line, as read from the file:
 *
 *   mixture <w_self> <n> {"<file>" <w_i>} x n
 *
 * The segment's own R-matrix value is multiplied by w_self (the abundance of
 * this project's target isotope, times its stoichiometry), and each table is
 * added with its weight after being folded with the segment's resolution
 * function.  Weights are plain numbers, not fit parameters.
 */
struct MixtureSpec {
  bool has = false;
  double selfWeight = 1.0;
  std::vector<std::pair<std::string, double> > tables;
  /*!
   * Reads the block following the keyword.  Returns false (and leaves the
   * spec unset) when the block is malformed.
   */
  bool Read(std::istream &stream);
};

#endif
