// The level matrix and the collision matrix of one J^pi block at one energy for pyazr
// (pyazr.physics.collision_matrix; docs/upgrade/physics_functions.md, stage 4): the arrays of
// physics/CollisionMatrix.h as numpy arrays.  No session, no project and no physics in this file: it checks the
// shapes, forms L_c - B_c and the two phase factors from the caller's channel functions, and calls
// azure2::physics::collision_matrix.
//
// The levels are formal (standard R-matrix) parameters, Brune's (2002) or Park's (2021): the two alternative forms
// read the shift function at each level energy (and its derivative where two levels coincide), Park's the norm J_l.
// A channel is a particle channel or a photon channel (L_c - B_c = 0: it does not enter the level matrix, but for
// the Reich-Moore term of use_rmc).  Every choice is an argument of fill_level_matrix; this file names them.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/CollisionMatrix.h"
#include "azure2/physics/LevelPairShift.h"

namespace py = pybind11;
namespace phys = azure2::physics;

namespace {

using RealArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using ComplexArray = py::array_t<std::complex<double>>;

/// A per-channel array: one finite value per channel, or an error naming the argument.
const double *ChannelValues(const RealArray &a, py::ssize_t nc, const char *name) {
  if (a.ndim() != 1 || a.shape(0) != nc)
    throw std::invalid_argument(std::string("collision_matrix: ") + name + " must hold one value per channel (" +
                                std::to_string(nc) + ", the second dimension of reduced_widths)");
  for (py::ssize_t c = 0; c < nc; c++)
    if (!std::isfinite(a.data()[c]))
      throw std::invalid_argument(std::string("collision_matrix: ") + name + " of channel " + std::to_string(c) +
                                  " is not finite");
  return a.data();
}

ComplexArray Square(const std::vector<std::complex<double>> &values, py::ssize_t n) {
  ComplexArray out({n, n});
  for (py::ssize_t i = 0; i < n * n; i++) out.mutable_data()[i] = values[i];
  return out;
}

/// A per-level, per-channel array (levels x channels) of finite values, or an error naming the argument.
const double *LevelChannelValues(const RealArray &a, py::ssize_t n, py::ssize_t nc, const char *name) {
  if (a.ndim() != 2 || a.shape(0) != n || a.shape(1) != nc)
    throw std::invalid_argument(std::string("collision_matrix: ") + name + " must have the shape (levels, channels) = (" +
                                std::to_string(n) + ", " + std::to_string(nc) + "), as reduced_widths");
  for (py::ssize_t i = 0; i < n * nc; i++)
    if (!std::isfinite(a.data()[i]))
      throw std::invalid_argument(std::string("collision_matrix: ") + name + " of level " + std::to_string(i / nc) +
                                  ", channel " + std::to_string(i % nc) + " is not finite");
  return a.data();
}

/// What a photon channel must not carry: its L_c - B_c is zero and it has no phases (the engine's channel tables).
void RequireZeroOfPhoton(double value, const char *name, py::ssize_t c) {
  if (value != 0.0)
    throw std::invalid_argument(std::string("collision_matrix: ") + name + " of channel " + std::to_string(c) +
                                " must be 0: it is a photon channel (photon_channels), which has no shift "
                                "function, boundary value or phase");
}

py::dict CollisionMatrixOf(const RealArray &level_energies, const RealArray &reduced_widths, double energy,
                           const RealArray &penetrability, const RealArray &shift, const py::object &boundary,
                           const RealArray &hard_sphere_phase, const RealArray &coulomb_phase,
                           const std::string &parametrization, const py::object &shift_at_levels,
                           const py::object &shift_derivative_at_levels_per_MeV, const py::object &park_norms,
                           const py::object &photon_channels, bool use_rmc) {
  const bool formal = parametrization == "formal", park = parametrization == "park";
  const bool brune = parametrization == "brune" || park;  // the engine's two flags: Park is set with Brune
  if (!formal && !brune)
    throw std::invalid_argument("collision_matrix: parametrization '" + parametrization +
                                "' is not known; allowed: 'formal', 'brune', 'park'");
  if (level_energies.ndim() != 1)
    throw std::invalid_argument("collision_matrix: level_energies must be one-dimensional (one energy per level)");
  const py::ssize_t n = level_energies.shape(0);
  if (reduced_widths.ndim() != 2 || reduced_widths.shape(0) != n)
    throw std::invalid_argument("collision_matrix: reduced_widths must have the shape (levels, channels), with " +
                                std::to_string(n) + " levels as level_energies");
  const py::ssize_t nc = reduced_widths.shape(1);
  if (!std::isfinite(energy)) throw std::invalid_argument("collision_matrix: energy is not finite");
  for (py::ssize_t l = 0; l < n; l++) {
    if (!std::isfinite(level_energies.data()[l]))
      throw std::invalid_argument("collision_matrix: the energy of level " + std::to_string(l) + " is not finite");
    for (py::ssize_t c = 0; c < nc; c++)
      if (!std::isfinite(reduced_widths.data()[l * nc + c]))
        throw std::invalid_argument("collision_matrix: the reduced width of level " + std::to_string(l) +
                                    ", channel " + std::to_string(c) + " is not finite");
  }
  const double *P = ChannelValues(penetrability, nc, "penetrability");
  const double *S = ChannelValues(shift, nc, "shift");
  const double *phi = ChannelValues(hard_sphere_phase, nc, "hard_sphere_phase");
  const double *omega = ChannelValues(coulomb_phase, nc, "coulomb_phase");

  // The channel kinds.  The level matrix asks only "particle or photon" ('E' stands for either multipole type).
  std::vector<char> radType(nc, 'P');
  bool anyPhoton = false;
  if (!photon_channels.is_none()) {
    const auto photon = py::cast<py::array_t<bool, py::array::c_style | py::array::forcecast>>(photon_channels);
    if (photon.ndim() != 1 || photon.shape(0) != nc)
      throw std::invalid_argument("collision_matrix: photon_channels must hold one True or False per channel (" +
                                  std::to_string(nc) + ", the second dimension of reduced_widths)");
    for (py::ssize_t c = 0; c < nc; c++)
      if (photon.data()[c]) {
        radType[c] = 'E';
        anyPhoton = true;
      }
  }
  if (use_rmc && !formal)
    throw std::invalid_argument("collision_matrix: use_rmc needs parametrization 'formal' (the engine's rule: the "
                                "Reich-Moore term is defined with formal parameters), not '" + parametrization + "'");
  if (use_rmc && !anyPhoton)
    throw std::invalid_argument("collision_matrix: use_rmc adds a term for each photon channel, and photon_channels "
                                "names none");

  // B_c: the constant of the formal parameters.  Brune's and Park's level matrices do not contain it (the level
  // energies take its place), so it is refused there instead of being read and dropped.
  RealArray boundaryValues;
  std::vector<double> noBoundary(nc, 0.0);
  const double *B = noBoundary.data();
  if (formal) {
    if (boundary.is_none())
      throw std::invalid_argument("collision_matrix: boundary (B_c, one value per channel) is required for "
                                  "parametrization 'formal'");
    boundaryValues = py::cast<RealArray>(boundary);
    B = ChannelValues(boundaryValues, nc, "boundary");
  } else if (!boundary.is_none()) {
    throw std::invalid_argument("collision_matrix: boundary is not an argument of parametrization '" +
                                parametrization + "': its level matrix has S_c(E_l) (shift_at_levels) where the "
                                "formal one has B_c");
  }

  // S_c(E_l) and dS_c/dE at E_l: Brune and Park.  The derivative is read only for two levels closer than
  // kCoincidentLevelGap_MeV (LevelPairShift.h: the chord between them becomes its Taylor form), so it is asked
  // for only then; the second derivative enters no element of the matrix (the gradients read it) and is zero here.
  RealArray levelShifts, levelSlopes;
  std::vector<double> zeros;
  const double *SL = nullptr, *dSL = nullptr, *J = nullptr;
  RealArray norms;
  if (brune) {
    if (shift_at_levels.is_none())
      throw std::invalid_argument("collision_matrix: shift_at_levels (S_c at each level energy, levels x channels) "
                                  "is required for parametrization '" + parametrization + "'");
    levelShifts = py::cast<RealArray>(shift_at_levels);
    SL = LevelChannelValues(levelShifts, n, nc, "shift_at_levels");
    zeros.assign(static_cast<std::size_t>(n * nc), 0.0);
    dSL = zeros.data();
    if (!shift_derivative_at_levels_per_MeV.is_none()) {
      levelSlopes = py::cast<RealArray>(shift_derivative_at_levels_per_MeV);
      dSL = LevelChannelValues(levelSlopes, n, nc, "shift_derivative_at_levels_per_MeV");
    } else {
      for (py::ssize_t l = 0; l < n; l++)
        for (py::ssize_t m = l + 1; m < n; m++)
          if (phys::levels_coincide(level_energies.data()[l], level_energies.data()[m]))
            throw std::invalid_argument(
                "collision_matrix: levels " + std::to_string(l) + " and " + std::to_string(m) +
                " are closer than " + std::to_string(phys::kCoincidentLevelGap_MeV * 1e6) +
                " eV: parametrization '" + parametrization + "' then reads shift_derivative_at_levels_per_MeV "
                "(dS_c/dE at each level energy, levels x channels), which was not given");
    }
    if (park) {
      if (park_norms.is_none())
        throw std::invalid_argument("collision_matrix: park_norms (J_l = 1 - sum_c gamma_lc^2 dS_c/dE at E_l over "
                                    "the particle channels, one per level) is required for parametrization 'park'");
      norms = py::cast<RealArray>(park_norms);
      if (norms.ndim() != 1 || norms.shape(0) != n)
        throw std::invalid_argument("collision_matrix: park_norms must hold one value per level (" +
                                    std::to_string(n) + ")");
      for (py::ssize_t l = 0; l < n; l++)
        if (!std::isfinite(norms.data()[l]))
          throw std::invalid_argument("collision_matrix: park_norms of level " + std::to_string(l) + " is not finite");
      J = norms.data();
    }
  } else {
    for (const auto &given : {std::make_pair(&shift_at_levels, "shift_at_levels"),
                              std::make_pair(&shift_derivative_at_levels_per_MeV, "shift_derivative_at_levels_per_MeV"),
                              std::make_pair(&park_norms, "park_norms")})
      if (!given.first->is_none())
        throw std::invalid_argument(std::string("collision_matrix: ") + given.second +
                                    " is not an argument of parametrization 'formal' (allowed with 'brune' or "
                                    "'park'; park_norms with 'park')");
    zeros.assign(static_cast<std::size_t>(n * nc), 0.0);
    SL = zeros.data();  // fill_level_matrix forms the pointer of each level's row; no element is read
  }
  if (brune && !park && !park_norms.is_none())
    throw std::invalid_argument("collision_matrix: park_norms is not an argument of parametrization 'brune' "
                                "(allowed with 'park')");

  std::vector<std::complex<double>> lo(nc), coulomb(nc), hardSphere(nc);
  std::vector<double> sqrtP(nc);
  for (py::ssize_t c = 0; c < nc; c++) {
    if (P[c] < 0.0)
      throw std::invalid_argument("collision_matrix: the penetrability of channel " + std::to_string(c) +
                                  " is negative (0 for a closed channel)");
    if (radType[c] != 'P') {
      RequireZeroOfPhoton(S[c], "shift", c);
      RequireZeroOfPhoton(B[c], "boundary", c);
      RequireZeroOfPhoton(phi[c], "hard_sphere_phase", c);
      RequireZeroOfPhoton(omega[c], "coulomb_phase", c);
      if (brune)
        for (py::ssize_t l = 0; l < n; l++) RequireZeroOfPhoton(SL[l * nc + c], "shift_at_levels", c);
    }
    lo[c] = radType[c] == 'P' ? std::complex<double>(S[c] - B[c], P[c]) : std::complex<double>(0.0, 0.0);
    coulomb[c] = std::exp(std::complex<double>(0.0, omega[c]));
    hardSphere[c] = std::exp(std::complex<double>(0.0, -phi[c]));
    sqrtP[c] = std::sqrt(P[c]);
  }
  phys::LevelMatrixTables tables;
  tables.n = static_cast<int>(n);
  tables.nc = static_cast<int>(nc);
  tables.E = level_energies.data();
  tables.gamma = reduced_widths.data();
  tables.shift = SL;
  tables.shift_d1 = dSL ? dSL : zeros.data();
  tables.shift_d2 = zeros.data();
  tables.park_norm = J;
  tables.lo = lo.data();
  tables.boundary = B;
  tables.rad_type = radType.data();
  tables.brune = brune;
  tables.park = park;
  tables.rmc = use_rmc;
  phys::ChannelTablesAtE channels;
  channels.coulomb_phase = coulomb.data();
  channels.hard_sphere_phase = hardSphere.data();
  channels.sqrt_penetrability = sqrtP.data();
  const phys::CollisionMatrix result = phys::collision_matrix(tables, energy, channels);
  py::dict out;
  out["A"] = Square(result.A, n);
  out["U"] = Square(result.U, nc);
  out["T"] = Square(result.T, nc);
  out["singular"] = result.singular;
  return out;
}

}  // namespace

void BindCollisionMatrix(py::module_ &m) {
  m.def("collision_matrix", &CollisionMatrixOf, py::arg("level_energies"), py::arg("reduced_widths"),
        py::arg("energy"), py::arg("penetrability"), py::arg("shift"), py::arg("boundary"),
        py::arg("hard_sphere_phase"), py::arg("coulomb_phase"), py::arg("parametrization") = "formal",
        py::arg("shift_at_levels") = py::none(), py::arg("shift_derivative_at_levels_per_MeV") = py::none(),
        py::arg("park_norms") = py::none(), py::arg("photon_channels") = py::none(), py::arg("use_rmc") = false,
        "The level matrix A (levels x levels) and the collision and T-matrices U, T (channels x channels) of one "
        "J^pi block at one energy, from physics/CollisionMatrix.h: a dict of A, U, T (complex arrays) and singular "
        "(True: the inverse level matrix has no inverse and A, U, T are zero).  level_energies (MeV), "
        "reduced_widths (MeV^1/2, levels x channels), energy (MeV, the scale of the level energies), and per "
        "channel at that energy: penetrability P_c (0: closed), shift S_c, boundary B_c (formal; None otherwise), "
        "hard_sphere_phase phi_c (rad), coulomb_phase omega_c (rad).  parametrization 'formal', 'brune' or 'park': "
        "the last two read shift_at_levels (S_c(E_l), levels x channels), shift_derivative_at_levels_per_MeV "
        "(dS_c/dE at E_l, MeV^-1; needed for two levels closer than 1 eV) and, Park, park_norms (J_l).  "
        "photon_channels: True per photon channel (no shift, boundary or phase; its penetrability is the caller's "
        "factor); use_rmc: the Reich-Moore term i gamma^2 of each photon channel on the diagonal (formal).  "
        "pyazr.physics.collision_matrix is the documented call.");
}
