#ifndef GAMMA_DECAY_FILL_H
#define GAMMA_DECAY_FILL_H

#include <complex>
#include <vector>

/// What a fill leaves in the level state for one gamma_decay record (LevelState::GammaDecayFillOf; GammaDecayMap.h):
/// its numbers and, per slot (photon channel), the parts of the engine's amplitude
/// gamma_int = sign(g) (g^2 - (Im gamma^ext)^2)^1/2 - Re gamma^ext.
struct GammaDecayFill {
  std::vector<double> q;                         ///< the record's numbers, slot order (Gamma_gamma in eV)
  std::vector<double> total;                     ///< g_k, the amplitude of the record's observed width (MeV^1/2)
  std::vector<std::complex<double>> external;    ///< gamma^ext_k (MeV^1/2); 0 without external part
  std::vector<unsigned char> floor;              ///< 1: g_k^2 < (Im gamma^ext_k)^2, gamma_int = -Re gamma^ext_k
  /// d gamma^ext_k / d gamma_lambda c per channel of the level's group, and d gamma^ext_k / d gamma~_f d per channel
  /// of the final state's group (N_f's dependence included); empty without external part.
  std::vector<std::vector<std::complex<double>>> dLevel, dFinal;
};

#endif
