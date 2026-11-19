// AZUREAPI: one J^pi block of the session's levels on a grid of excitation energies (JpiBlockReport.h;
// docs/upgrade/physics_tab.md section 4): pyazr's Session.jpi_block, the facade's PhysicsSession::Block.
//
// Nothing here is a formula of its own.  The levels are the session's LevelState at the vector p, filled as
// ObservedLevelRowsAllRWA fills it (CompoundSystem::FillCompoundFromParams, then CalcShiftFunctions); the arrays of the
// block are gathered as AMatrixFunc::FillLevelMatrix gathers them (the active levels, E, gamma, the shifts at the
// levels, Park's norms, B_c of a coupled particle channel of Brune and Park); the channel functions at E_x are the
// stage's channel rows, ChannelTable::Compute (P, S, exp(-i phi), exp(i(sigma_l - sigma_0))), with the zero-energy
// shift exactly at a threshold (ScatteringLength::ZeroEnergyShift), as Site::CalcEDependentValues makes them; A, U and
// T are physics::collision_matrix (the level route, which the engine's channel route agrees with to round-off), the
// phase shift and inelasticity the engine's phase_shift data kind (physics/Observables.h), and det A^-1
// physics::level_matrix_determinant of fill_level_matrix's matrix.
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>

#include "AChannel.h"
#include "ALevel.h"
#include "AZUREAPI.h"
#include "ChannelTable.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "JGroup.h"
#include "JpiBlockReport.h"
#include "LevelState.h"
#include "PPair.h"
#include "ScatteringLength.h"
#include "azure2/physics/CollisionMatrix.h"
#include "azure2/physics/LevelMatrixDeterminant.h"
#include "azure2/physics/Observables.h"

namespace {

/// The pair's type in PPair::GetPType: 0 particle, 10 photon, 20 beta.
constexpr int kParticlePair = 0;

}  // namespace

JpiBlockReport AZUREAPI::JpiBlockAllRWA(const vector_r &p, int jgroup, const vector_r &excitation) {
  if (p.size() != all_rwa_.size())
    throw std::invalid_argument("J^pi block: expected every parameter (" + std::to_string(all_rwa_.size()) +
                                " values, fixed ones included), got " + std::to_string(p.size()));
  const CompoundSystem *lc = compound();
  if (jgroup < 1 || jgroup > lc->NumJGroups())
    throw std::invalid_argument("J^pi block: no J group " + std::to_string(jgroup) + " (the project has " +
                                std::to_string(lc->NumJGroups()) + ", numbered from 1 as parameters.out)");
  const JGroup *group = lc->GetJGroup(jgroup);
  if (!group->IsInRMatrix())
    throw std::invalid_argument("J^pi block: J group " + std::to_string(jgroup) + " has no level in the R-matrix");
  lc->FillCompoundFromParams(p, state());
  lc->CalcShiftFunctions(configure(), state());
  const LevelState &st = state();
  const bool park = (configure().paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  JpiBlockReport r;
  r.jgroup = jgroup;
  r.J = group->GetJ();
  r.pi = group->GetPi();
  r.brune = brune;
  r.park = park;
  for (int la = 1; la <= group->NumLevels(); la++)
    if (group->GetLevel(la)->IsInRMatrix()) {
      r.levels.push_back(la);
      r.levelEnergy.push_back(st.E(jgroup, la));
    }
  const int n = (int)r.levels.size(), nc = group->NumChannels();

  // The level parameters, as AMatrixFunc::FillLevelMatrix gathers them.
  const size_t cells = (size_t)n * nc;
  std::vector<double> E(n), gamma(cells), shift(cells), d1(cells, 0.0), d2(cells, 0.0), parkNorm(n, 1.0);
  std::vector<double> boundary(nc, 0.0);
  std::vector<char> radType(nc);
  for (int i = 0; i < n; i++) {
    const int la = r.levels[i];
    E[i] = st.E(jgroup, la);
    for (int ch = 1; ch <= nc; ch++) {
      gamma[(size_t)i * nc + ch - 1] = st.Gamma(jgroup, la, ch);
      shift[(size_t)i * nc + ch - 1] = st.Shift(jgroup, la, ch);
      if (park || brune) {
        d1[(size_t)i * nc + ch - 1] = st.ShiftDerivative(jgroup, la, ch);
        d2[(size_t)i * nc + ch - 1] = st.ShiftSecondDerivative(jgroup, la, ch);
      }
    }
    if (park) parkNorm[i] = st.ParkNorm(jgroup, la);
  }
  for (int ch = 1; ch <= nc; ch++) {
    const AChannel *channel = group->GetChannel(ch);
    const PPair *pair = lc->GetPair(channel->GetPairNum());
    radType[ch - 1] = channel->GetRadType();
    JpiBlockReport::Channel c;
    c.channel = ch;
    c.pairKey = pair->GetPairKey();
    c.radType = channel->GetRadType();
    c.l = channel->GetL();
    c.s = channel->GetS();
    c.threshold = pair->GetSepE() + pair->GetExE();
    c.identical = pair->IsIdentical();
    if (pair->GetPType() == kParticlePair) c.radius = channel->GetRadius();
    if ((park || brune) && radType[ch - 1] == 'P')
      for (int i = 0; i < n; i++)
        if (std::fabs(gamma[(size_t)i * nc + ch - 1]) >= azure2::physics::kWidthCutoff) {
          boundary[ch - 1] = channel->GetBoundaryCondition();
          break;
        }
    c.boundary = boundary[ch - 1];
    r.channels.push_back(c);
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::complex<double> cnan(nan, nan);
  const size_t ne = excitation.size();
  r.ex.assign(excitation.begin(), excitation.end());
  r.ec.assign(ne * nc, 0.0);
  r.P.assign(ne * nc, 0.0);
  r.S.assign(ne * nc, 0.0);
  r.lo.assign(ne * nc, 0.0);
  r.expCoulomb.assign(ne * nc, 0.0);
  r.expHardSphere.assign(ne * nc, 0.0);
  r.A.assign(ne * n * n, 0.0);
  r.U.assign(ne * nc * nc, cnan);
  r.T.assign(ne * nc * nc, cnan);
  r.delta.assign(ne * nc, nan);
  r.eta.assign(ne * nc, nan);
  r.detInverse.assign(ne, cnan);
  r.singular.assign(ne, 0);
  std::vector<std::complex<double>> lo(nc), expCP(nc), expHSP(nc);
  std::vector<double> sqrtP(nc);
  for (size_t e = 0; e < ne; e++) {
    const double ex = excitation[e];
    std::vector<bool> particle(nc, false);
    for (int ch = 1; ch <= nc; ch++) {
      const AChannel *channel = group->GetChannel(ch);
      const PPair *pair = lc->GetPair(channel->GetPairNum());
      const double ec = ex - pair->GetSepE() - pair->GetExE();
      const size_t k = e * nc + ch - 1;
      r.ec[k] = ec;
      lo[ch - 1] = 0.0, expCP[ch - 1] = 1.0, expHSP[ch - 1] = 1.0, sqrtP[ch - 1] = 0.0;
      if (pair->GetPType() != kParticlePair) continue;  // a photon or beta channel: lo = 0, as Site makes it
      particle[ch - 1] = ec > 0.0;  // U and T are given between open particle channels
      const int l = channel->GetL();
      if (std::fabs(ec) <= ScatteringLength::kThresholdTolerance) {
        r.S[k] = ScatteringLength::ZeroEnergyShift(pair, l, channel->GetRadius());
        lo[ch - 1] = std::complex<double>(r.S[k] - channel->GetBoundaryCondition(), 0.0);
      } else {
        const ChannelTable::Row row = ChannelTable::Compute(pair, l, pair->SpinOrbitTwoJ(group->GetJ()),
                                                            channel->GetRadius(), ec, configure());
        r.P[k] = row.P, r.S[k] = row.S;
        lo[ch - 1] = std::complex<double>(row.S - channel->GetBoundaryCondition(), row.P);
        sqrtP[ch - 1] = std::sqrt(row.P);
        expCP[ch - 1] = row.expCP;
        expHSP[ch - 1] = row.expHSP;
      }
      r.lo[k] = lo[ch - 1];
      r.expCoulomb[k] = expCP[ch - 1];
      r.expHardSphere[k] = expHSP[ch - 1];
    }
    azure2::physics::LevelMatrixTables tables;
    tables.n = n;
    tables.nc = nc;
    tables.E = E.data();
    tables.gamma = gamma.data();
    tables.shift = shift.data();
    tables.shift_d1 = d1.data();
    tables.shift_d2 = d2.data();
    tables.park_norm = parkNorm.data();
    tables.lo = lo.data();
    tables.boundary = boundary.data();
    tables.rad_type = radType.data();
    tables.brune = brune;
    tables.park = park;
    tables.rmc = (configure().paramMask & Config::USE_RMC_FORMALISM);
    azure2::physics::LevelMatrixSolver::matrix_c inverse;
    azure2::physics::fill_level_matrix(tables, ex, inverse);
    r.detInverse[e] = azure2::physics::level_matrix_determinant(inverse);
    azure2::physics::ChannelTablesAtE functions;
    functions.coulomb_phase = expCP.data();
    functions.hard_sphere_phase = expHSP.data();
    functions.sqrt_penetrability = sqrtP.data();
    const azure2::physics::CollisionMatrix m = azure2::physics::collision_matrix(tables, ex, functions);
    r.singular[e] = m.singular ? 1 : 0;
    for (size_t i = 0; i < m.A.size(); i++) r.A[e * n * n + i] = m.A[i];
    for (int c = 0; c < nc; c++)
      for (int cp = 0; cp < nc; cp++) {
        if (!particle[c] || !particle[cp]) continue;
        r.U[e * nc * nc + (size_t)c * nc + cp] = m.u(c, cp);
        r.T[e * nc * nc + (size_t)c * nc + cp] = m.t(c, cp);
      }
    for (int c = 0; c < nc; c++) {
      if (!particle[c]) continue;
      // the engine's phase_shift data kind (GenMatrixFunc::CalculateBareCrossSection)
      const std::complex<double> u = azure2::physics::collision_from_t_matrix(expCP[c], m.t(c, c));
      r.delta[e * nc + c] = azure2::physics::phase_shift_degrees(u, r.channels[c].identical);
      r.eta[e * nc + c] = std::abs(u);
    }
  }
  return r;
}
