#ifndef ELECTRONBINDING_H
#define ELECTRONBINDING_H

/*!
 * Total electron binding energy of the neutral atom, Z = 1..92, in eV: the sum of all its ionization
 * energies, NIST Atomic Spectra Database (fetched 2026-10-21).
 */
inline double TotalElectronBindingEv(int z) {
  static const double kBinding[92] = {
    13.598, 79.005, 203.486, 399.149, 670.984, 1030.108,
    1486.058, 2043.843, 2715.891, 3511.697, 4419.969, 5451.072,
    6604.952, 7888.522, 9305.790, 10859.598, 12556.443, 14400.856,
    16382.422, 18510.103, 20787.893, 23221.324, 25819.344, 28582.028,
    31514.095, 34619.321, 37900.131, 41354.904, 45002.586, 48840.900,
    52861.386, 57077.608, 61470.775, 66068.036, 70863.003, 75865.753,
    81048.431, 86440.143, 92029.404, 97825.640, 103830.374, 110065.684,
    116487.808, 123150.963, 130033.151, 137141.044, 144480.603, 152046.287,
    159851.438, 167892.692, 176144.825, 184652.775, 193401.698, 202399.388,
    211628.360, 221101.222, 230824.919, 240804.538, 251060.196, 261584.448,
    272377.390, 283454.942, 294803.500, 306422.726, 318334.267, 330558.316,
    343075.932, 355895.024, 369029.959, 382457.206, 396209.334, 410249.555,
    424569.672, 438987.004, 454144.434, 469448.058, 485001.507, 500896.259,
    517177.249, 533832.424, 550807.349, 568143.519, 585860.985, 603926.518,
    622296.318, 641057.438, 660204.973, 679733.026, 699585.066, 719895.575,
    740488.990, 761522.794,
  };
  return (z >= 1 && z <= 92) ? kBinding[z - 1] : -1.0;
}

/*!
 * Adiabatic limit of the electron-screening potential for a projectile Z1 on a neutral target Z2 (eV):
 * the binding the united atom gains, U_ad = B(Z1 + Z2) - B(Z1) - B(Z2).  934.7 eV for alpha + 13C.
 * -1 for a neutral particle or Z1 + Z2 > 92.
 */
inline double AdiabaticScreeningEv(int z1, int z2) {
  if (z1 <= 0 || z2 <= 0 || z1 + z2 > 92) return -1.0;
  return TotalElectronBindingEv(z1 + z2) - TotalElectronBindingEv(z1) - TotalElectronBindingEv(z2);
}

#endif

