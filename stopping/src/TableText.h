#ifndef AZURE2_STOPPING_TABLETEXT_H
#define AZURE2_STOPPING_TABLETEXT_H

#include <string>

namespace azure2 {
namespace stopping {

/// The text of stopping/data/srim/SCOEF03.DAT and SNUC03.dat as compiled in by CMake (TableText.cpp.in), carriage
/// returns removed; empty when the build did not find the files (the owner may choose not to ship them).
std::string EmbeddedScoefText();
std::string EmbeddedSnucText();

}  // namespace stopping
}  // namespace azure2

#endif
