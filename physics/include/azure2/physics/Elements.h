#ifndef AZURE2_PHYSICS_ELEMENTS_H
#define AZURE2_PHYSICS_ELEMENTS_H

#include <string>

namespace azure2 {
namespace physics {

/*!
 * Element symbols by charge number Z, 0 (the neutron, "n") to 118 (Og), as in the AME2020 mass table
 * (pyazr/data/mass_1.mas20, which the test physics_elements compares symbol by symbol).  The one table of the C++
 * code: format-3 reaction labels (`15N(p,a)12C`, experiment/) are read with it (docs/upgrade/project_format3.md
 * section 7.3).
 */
constexpr int kMaxElementZ = 118;

/// The symbol of charge number z ("n" for 0, "H", "He", ...), or "" outside 0..kMaxElementZ.
std::string ElementSymbol(int z);

/// The charge number of a symbol as written in a nuclide label (case-sensitive: "He", not "HE"); -1 if unknown.
/// "n" gives 0.
int ElementZ(const std::string &symbol);

}  // namespace physics
}  // namespace azure2

#endif
