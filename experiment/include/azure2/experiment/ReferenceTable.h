#ifndef AZURE2_EXPERIMENT_REFERENCETABLE_H
#define AZURE2_EXPERIMENT_REFERENCETABLE_H

#include <iosfwd>
#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * A tabulated reference cross section sigma_ref(E) (quantity ratio_to_reference, field 'reference <file> [unit <u>]';
 * docs/upgrade/data_kinds.md section 4.9): the data are sigma(E) / sigma_ref(E).  The file has two columns, the
 * energy in the data set's energy role, frame and unit (an extrapolation: its grid's frame and unit) and the cross
 * section in the field's unit (b by default); '#' starts a comment.  Build reads it once and keeps it in the spec
 * with the energies taken to the entrance pair's lab frame in MeV and the values to b, so that the engine looks it
 * up at a point's (or sub-point's) lab energy whatever frame the file was written in.
 */
struct ReferenceTable {
  std::string file;                    ///< as written in the record (for --describe and messages)
  std::vector<double> energyLab_MeV;   ///< strictly increasing
  std::vector<double> sigma_b;         ///< > 0
  bool Empty() const { return energyLab_MeV.empty(); }
};

/// The rows of a reference file: (energy, value) as written, or an error ("" on success): a line with fewer than two
/// numbers, energies that do not increase strictly, a value <= 0 or not finite, fewer than two rows.
std::string ReadReferenceRows(std::istream &in, std::vector<double> &energy, std::vector<double> &value);

/// sigma_ref (b) at a lab energy (MeV): linear in ln sigma against E between the two rows around it (the end rows
/// included); NaN outside [first, last] energy (the set-up refuses data and grids there; a fold's sub-point or a
/// shifted energy outside makes the point's model NaN).
double ReferenceAt(const ReferenceTable &table, double eLab_MeV);

}  // namespace experiment
}  // namespace azure2

#endif
