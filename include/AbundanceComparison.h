#ifndef ABUNDANCECOMPARISON_H
#define ABUNDANCECOMPARISON_H
// The comparison of the abundances and THM coherent values held by the compound and by the working copy while both
// are kept (removed with the compound's copies).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

inline void CompareCopies(const std::vector<double> &compound, const std::vector<double> &data, const char *what,
                          const char *where) {
  if (compound.size() == data.size() &&
      (compound.empty() || std::memcmp(compound.data(), data.data(), compound.size() * sizeof(double)) == 0))
    return;
  std::printf("%s COMPARISON: the compound's and the working copy's differ at %s (%zu, %zu values)\n", what, where,
              compound.size(), data.size());
  for (std::size_t k = 0; k < compound.size() || k < data.size(); k++)
    std::printf("  %zu: %.17g %.17g\n", k, k < compound.size() ? compound[k] : -1.0, k < data.size() ? data[k] : -1.0);
  std::fflush(stdout);
  std::fprintf(stderr, "%s COMPARISON failed at %s\n", what, where);
  std::abort();
}

#endif
