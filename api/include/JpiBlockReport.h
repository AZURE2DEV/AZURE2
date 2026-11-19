#ifndef JPI_BLOCK_REPORT_H
#define JPI_BLOCK_REPORT_H

#include <complex>
#include <vector>

/*!
 * One J^pi block of a session's level scheme on a grid of excitation energies (AZUREAPI::JpiBlockAllRWA;
 * docs/upgrade/physics_tab.md section 4): the channel functions of its channels, the level matrix A, the collision
 * and T matrices of its open particle channels and det A^-1, at each E_x.  Plain data, as the engine's arrays give them:
 * pyazr's Session.jpi_block and the front ends' facade (azure2::front::PhysicsSession::Block) read it.
 */
struct JpiBlockReport {
  int jgroup = 0;  ///< 1-based, as parameters.out
  double J = 0.0;
  int pi = 1;
  bool brune = false, park = false;  ///< the parametrization of the levels
  /// The levels in the R-matrix (ALevel::IsInRMatrix), 1-based in the group, and their energies E_lambda (MeV, E_x
  /// of the compound nucleus, in the session's parametrization: the formal pole energy of a standard level).
  std::vector<int> levels;
  std::vector<double> levelEnergy;
  struct Channel {
    int channel = 0;      ///< 1-based in the group
    int pairKey = 0;
    char radType = 'P';   ///< 'P' particle, 'E' or 'M' photon (and its multipolarity l), 'B' beta
    int l = 0;
    double s = 0.0;
    double radius = 0.0;     ///< a_c (fm); 0 for a photon or beta channel
    double threshold = 0.0;  ///< E_x of the pair's threshold (separation + excitation energy, MeV)
    double boundary = 0.0;   ///< B_c (0 where the level matrix does not read it, as the engine)
    bool identical = false;  ///< the pair's particles are identical (the phase shift's range)
  };
  std::vector<Channel> channels;
  std::vector<double> ex;  ///< the grid, E_x (MeV)
  /// Per energy e and channel c, [e * nc + c]: E_c = E_x - threshold (MeV, c.m. of the channel's pair), P_c, S_c (the
  /// shift function at a_c; the closed-channel value below threshold), L_c - B_c = S_c - B_c + i P_c as the level
  /// matrix reads it, exp(i omega_c) (omega_c = sigma_l - sigma_0) and exp(-i phi_c).  A photon or beta channel: E_c
  /// as said, every other entry 0 (its penetrability is not a channel function of this table).
  std::vector<double> ec, P, S;
  std::vector<std::complex<double>> lo, expCoulomb, expHardSphere;
  /// Per energy, [e * n * n + l * n + m]: A_lm (MeV^-1).
  std::vector<std::complex<double>> A;
  /// Per energy, [e * nc * nc + c * nc + c']: U_cc' and T_cc' of the pathway c -> c' (physics/CollisionMatrix.h),
  /// for open particle channels c and c' (E_c > 0); NaN where either is closed, a photon or a beta channel.
  std::vector<std::complex<double>> U, T;
  /// Per energy and channel, [e * nc + c]: the nuclear phase shift delta_c (deg) and the inelasticity eta_c = |U_cc|
  /// of an open particle channel, as the engine's phase_shift data kind forms them (GenMatrixFunc: U = (e^{2 i omega} -
  /// T_cc)/e^{2 i omega}, delta = arg U / 2 in (-90, 90], [0, 180) for identical particles; physics/Observables.h);
  /// NaN for any other channel.
  std::vector<double> delta, eta;
  /// Per energy: det A^-1 (MeV^n) and whether A^-1 was singular there (A, U, T then zero).
  std::vector<std::complex<double>> detInverse;
  std::vector<char> singular;
};

#endif
