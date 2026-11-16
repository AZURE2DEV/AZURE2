// The level matrix and the collision matrix of one J^pi block at one energy for pyazr
// (pyazr.physics.collision_matrix; docs/upgrade/physics_functions.md, stage 4): the arrays of
// physics/CollisionMatrix.h as numpy arrays.  No session, no project and no physics in this file: it checks the
// shapes, forms L_c - B_c and the two phase factors from the caller's channel functions, and calls
// azure2::physics::collision_matrix.
//
// The levels are formal (standard R-matrix) parameters and every channel a particle channel; Brune's or Park's
// parameters and photon channels are arguments of fill_level_matrix that this door does not offer yet.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/CollisionMatrix.h"

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

py::dict CollisionMatrixOf(const RealArray &level_energies, const RealArray &reduced_widths, double energy,
                           const RealArray &penetrability, const RealArray &shift, const RealArray &boundary,
                           const RealArray &hard_sphere_phase, const RealArray &coulomb_phase) {
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
  const double *B = ChannelValues(boundary, nc, "boundary");
  const double *phi = ChannelValues(hard_sphere_phase, nc, "hard_sphere_phase");
  const double *omega = ChannelValues(coulomb_phase, nc, "coulomb_phase");
  std::vector<std::complex<double>> lo(nc), coulomb(nc), hardSphere(nc);
  std::vector<double> sqrtP(nc);
  const std::vector<char> radType(nc, 'P');
  for (py::ssize_t c = 0; c < nc; c++) {
    if (P[c] < 0.0)
      throw std::invalid_argument("collision_matrix: the penetrability of channel " + std::to_string(c) +
                                  " is negative (0 for a closed channel)");
    lo[c] = std::complex<double>(S[c] - B[c], P[c]);
    coulomb[c] = std::exp(std::complex<double>(0.0, omega[c]));
    hardSphere[c] = std::exp(std::complex<double>(0.0, -phi[c]));
    sqrtP[c] = std::sqrt(P[c]);
  }
  phys::LevelMatrixTables tables;
  tables.n = static_cast<int>(n);
  tables.nc = static_cast<int>(nc);
  tables.E = level_energies.data();
  tables.gamma = reduced_widths.data();
  tables.lo = lo.data();
  tables.boundary = B;
  tables.rad_type = radType.data();
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
        py::arg("hard_sphere_phase"), py::arg("coulomb_phase"),
        "The level matrix A (levels x levels) and the collision and T-matrices U, T (channels x channels) of one "
        "J^pi block of formal levels and particle channels at one energy, from physics/CollisionMatrix.h: a dict "
        "of A, U, T (complex arrays) and singular (True: the inverse level matrix has no inverse and A, U, T are "
        "zero).  level_energies (MeV), reduced_widths (MeV^1/2, levels x channels), energy (MeV, the scale of the "
        "level energies), and per channel at that energy: penetrability P_c (0: closed), shift S_c, boundary B_c, "
        "hard_sphere_phase phi_c (rad), coulomb_phase omega_c (rad).  pyazr.physics.collision_matrix is the "
        "documented call.");
}
