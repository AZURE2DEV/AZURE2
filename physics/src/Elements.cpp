// Element symbols by charge number (Elements.h).
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace physics {

namespace {

// IUPAC symbols, Z = 0 (the neutron, as AME2020 writes it) to 118.
const char *const kSymbols[kMaxElementZ + 1] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
    "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
    "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

}  // namespace

std::string ElementSymbol(int z) { return (z < 0 || z > kMaxElementZ) ? std::string() : std::string(kSymbols[z]); }

int ElementZ(const std::string &symbol) {
  for (int z = 0; z <= kMaxElementZ; z++)
    if (symbol == kSymbols[z]) return z;
  return -1;
}

}  // namespace physics
}  // namespace azure2
