// Element symbols (Elements.h) against the AME2020 mass table pyazr/data/mass_1.mas20, the file the masses come from
// (format a1,i3,i5,i5,i5,1x,a3,...: Z in columns 10-14, the symbol in 22-24): every Z of the file has our symbol,
// every Z of ours is in the file, and ElementZ inverts ElementSymbol.
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>

#include "azure2/physics/Elements.h"

using namespace azure2::physics;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: physics_elements_test mass_1.mas20\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::printf("cannot read %s\n", argv[1]);
    return 2;
  }
  std::map<int, std::set<std::string>> fileSymbols;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.size() < 123) continue;
    int z = 0;
    try {
      std::size_t used = 0;
      z = std::stoi(line.substr(9, 5), &used);
      std::stoi(line.substr(14, 5));
    } catch (...) {
      continue;
    }
    std::string symbol;
    for (char c : line.substr(20, 3))
      if (c != ' ') symbol += c;
    if (symbol.empty()) continue;
    fileSymbols[z].insert(symbol);
  }
  int failures = 0, checks = 0;
  auto check = [&](const std::string &what, bool ok) {
    checks++;
    if (!ok) {
      std::printf("  FAIL  %s\n", what.c_str());
      failures++;
    }
  };
  check("the file has Z = 0 .. 118", fileSymbols.size() == static_cast<std::size_t>(kMaxElementZ + 1) &&
                                         fileSymbols.begin()->first == 0 && fileSymbols.rbegin()->first == kMaxElementZ);
  for (const auto &kv : fileSymbols) {
    check("Z = " + std::to_string(kv.first) + ": one symbol in the file", kv.second.size() == 1);
    check("Z = " + std::to_string(kv.first) + ": " + ElementSymbol(kv.first) + " = " + *kv.second.begin(),
          ElementSymbol(kv.first) == *kv.second.begin());
  }
  for (int z = 0; z <= kMaxElementZ; z++) check("ElementZ(ElementSymbol(" + std::to_string(z) + "))", ElementZ(ElementSymbol(z)) == z);
  check("outside the table", ElementSymbol(-1).empty() && ElementSymbol(kMaxElementZ + 1).empty() && ElementZ("Xx") == -1 &&
                                 ElementZ("HE") == -1 && ElementZ("") == -1);
  std::printf("%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
