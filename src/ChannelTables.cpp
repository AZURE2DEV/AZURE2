// The energy-dependent tables of a point: its channel rows (L^0 elements, square roots of the penetrabilities,
// Coulomb and hard-sphere phase factors, THM vertex pieces), their recalculation after an energy change, their
// sharing between points at the same energy, the Coulomb amplitude, and the getters the matrix functions read.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "AChannel.h"
#include "CompoundSystem.h"
#include "ChannelTable.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "IntegratedFermiFunc.h"
#include "JGroup.h"
#include "Kinematics.h"
#include "azure2/physics/AstrophysicalFactor.h"
#include "PPair.h"
#include "ScatteringLength.h"
#include "ShftFunc.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"

/// Returns the geometrical cross section factor \f$ \frac{\pi}{k^2} \f$.

double Site::GetGeometricalFactor() const {
  return values_.geofactor;
}

/*!
 * Returns the multiplicative conversion factor from cross section to s-factor.
 */

double Site::GetSFactorConversion() const {
  return values_.sfactorconv;
}

/*!
 * Returns the square root of the penetrability for a channel specified by
 * positions in the JGroup and subsequent AChannel vectors.
 */

double Site::GetSqrtPenetrability(int jGroupNum, int channelNum) const {
  return values_.edep->penetrabilities[values_.edep->At(jGroupNum, channelNum)];
}

/*!
 * Returns the THM entrance transfer form factor M_l = (b-1) j_l - rho dj_l/drho
 * (+ the external Coulomb term, when stored) for the channel specified by positions in the JGroup and AChannel vectors,
 * assembled with the given boundary b (the per-level shift function under the
 * Brune formalism, or the channel boundary constant otherwise). Returns 0 if
 * not stored (a channel of another pair than the entrance pair, or a non-THM point); the vertex reads the entrance
 * pair's channels only (THMMatrixFunc).
 */

complex Site::GetThmFormFactor(int jGroupNum, int channelNum, complex boundary) const {
  if (!values_.thm_vertex) return 0.0;
  const ThmVertexRows &t = *values_.thm_vertex;
  const int k = t.Find(jGroupNum, channelNum);
  if (k < 0) return 0.0;
  return (boundary - 1.0) * t.jl[k] - t.rhodjl[k] + t.coul[k];
}

double Site::GetThmVertexEps() const { return values_.thm_vertex ? values_.thm_vertex->eps : 0.0; }

double Site::GetThmVertexEps(int node) const {
  const ThmPsTable *t = values_.thm_ps.get();
  return (t && node >= 0 && node < (int)t->eps.size()) ? t->eps[node] : 0.0;
}

complex Site::GetThmFormFactor(int jGroupNum, int channelNum, complex boundary, int node) const {
  const ThmPsTable *t = values_.thm_ps.get();
  if (!t || node < 0 || node >= t->nodes) return 0.0;
  if (jGroupNum < 1 || jGroupNum >= (int)t->offset.size() || channelNum < 1) return 0.0;
  const int flat = t->offset[jGroupNum - 1] + channelNum - 1;
  if (flat >= t->offset[jGroupNum]) return 0.0;
  const int k = t->slot[flat];
  if (k < 0) return 0.0;
  const size_t i = (size_t)k * t->nodes + node;
  const complex coul = t->coul.empty() ? complex(0.0, 0.0) : t->coul[i];
  return (boundary - 1.0) * t->jl[i] - t->rhodjl[i] + coul;
}

/*!
 * Returns the \f$ L_o \f$ diagonal matrix element for a channel specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex Site::GetLoElement(int jGroupNum, int channelNum) const {
  return values_.edep->lo_elements[values_.edep->At(jGroupNum, channelNum)];
}

/*!
 * Returns the factor \f$ \exp(\omega_c) \f$ where \f$ \omega_c \f$ is the Coulomb
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex Site::GetExpCoulombPhase(int jGroupNum, int channelNum) const {
  return values_.edep->coulombphase[values_.edep->At(jGroupNum, channelNum)];
}

/*!
 * Returns the factor \f$ \exp(\delta_c) \f$ where \f$ \delta_c \f$ is the hard sphere
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex Site::GetExpHardSpherePhase(int jGroupNum, int channelNum) const {
  return values_.edep->hardspherephase[values_.edep->At(jGroupNum, channelNum)];
}

/*!
 * Returns the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

complex Site::GetCoulombAmplitude() const {
  return values_.coulombamplitude;
}

complex Site::GetCoulombAmplitude(double s) const {
  // Every pair but an identical one with spin returns the stored amplitude
  // untouched, so the spin-0 and distinguishable results cannot move.
  if (!values_.spinDependentCoulomb) return values_.coulombamplitude;
  const int twoS = (int)std::lround(2.0 * s);
  // Channel spins of two identical particles are integers.
  const double sign = ((twoS / 2) % 2 == 0) ? 1.0 : -1.0;
  return values_.coulombdirect + sign * values_.coulombexchange;
}

/*!
 * Sets the geometrical cross section factor \f$ \frac{\pi}{k^2} \f$.
 */

void Site::SetGeometricalFactor(double geoFactor) {
  values_.geofactor = geoFactor;
}

/*!
 * Sets the multiplicative conversion from cross section to s-factor.
 */

void Site::SetSFactorConversion(double conversion) {
  values_.sfactorconv = conversion;
}

namespace {
/* Below E = -B the half-off-shell momentum is imaginary and the form factor
   (hence the HOES cross section) is left identically zero.  Physically the QF
   relative energy satisfies E + B > 0, so points here should only be far
   Gaussian-tail sub-points.  Reported once per session, so that a wrong
   binding energy cannot silently null a whole segment. */
void WarnThmBelowB(const Config &configure, double energy, double bindingE) {
  if (configure.thm.belowBWarned->exchange(true)) return;
  configure.outStream << "WARNING: THM point at E_cm = " << energy << " MeV lies below E = -B (B = " << bindingE
                      << " MeV); the transfer form factor and HOES cross section are set to zero there. Check the "
                      << "binding energy if this is unexpected. (warning shown once)" << std::endl;
}
}  // namespace

/*!
 * This function calculates several energy dependent quantities simultaniously.
 * This includes the geometrical cross section, the s-factor conversion, the \f$ L_o \f$ matrix
 * elements, the square root of the penetrability, and the exponentials of the Coulomb phase shifts
 * and hard sphere phase shfits.
 */

double Site::CompoundEnergy(const CompoundSystem *compound) {
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()));
  const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(this->GetExitKey()));
  if (entrancePair->GetPType() == 20) return this->GetCMEnergy() + exitPair->GetSepE() + exitPair->GetExE();
  return this->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
}

void Site::ChannelKeys(WorkingCopy *data, const CompoundSystem *compound, std::vector<ChannelTable::Key> &keys) {
  // The walk of CalcEDependentValues: the blocks of the entrance pair's compound nucleus, their particle channels,
  // E_c by the same expression; the exact threshold of a scattering length is not a row.
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()));
  const double inEnergy = this->CompoundEnergy(compound);
  const int nucleus = entrancePair->GetNucleus();
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *theJGroup = compound->GetJGroup(j);
    if (!theJGroup->IsInRMatrix() || theJGroup->GetNucleus() != nucleus) continue;
    for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
      const AChannel *theChannel = theJGroup->GetChannel(ch);
      const PPair *thePair = compound->GetPair(theChannel->GetPairNum());
      if (thePair->GetPType() != 0) continue;
      double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
      if (fabs(localEnergy) <= ScatteringLength::kThresholdTolerance) continue;
      keys.push_back({theChannel->GetPairNum(), theChannel->GetL(), thePair->SpinOrbitTwoJ(theJGroup->GetJ()), localEnergy});
    }
  }
}

std::tuple<int, std::uint64_t, int> Site::ChannelRowsKey(const CompoundSystem *compound) {
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()));
  const int nucleus = entrancePair->GetNucleus();
  const double inEnergy = this->CompoundEnergy(compound);
  std::uint64_t bits;
  std::memcpy(&bits, &inEnergy, sizeof bits);
  // The walk of CalcEDependentValues: a beta-decay channel's penetrability reads the exit pair's charge.
  int betaCharge = -1;
  for (int j = 1; j <= compound->NumJGroups() && betaCharge < 0; j++) {
    const JGroup *group = compound->GetJGroup(j);
    if (!group->IsInRMatrix() || group->GetNucleus() != nucleus) continue;
    for (int ch = 1; ch <= group->NumChannels(); ch++)
      if (compound->GetPair(group->GetChannel(ch)->GetPairNum())->GetPType() == 20) {
        const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(this->GetExitKey()));
        betaCharge = exitPair->GetZ(1) + exitPair->GetZ(2);
        break;
      }
  }
  return std::make_tuple(nucleus, bits, betaCharge);
}

void Site::CalcEDependentValues(WorkingCopy *data, const CompoundSystem *compound, const Config &configure,
                                const ChannelRowSource *rows, const Site *rowsOf) {
  AZUREProfile::Timer profileTimer(configure.profile, AZUREProfile::kEDependent);
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()));
  const PPair *exitPair = compound->GetPair(compound->GetPairNumFromKey(this->GetExitKey()));

  double inEnergy;
  double geofactor;
  double sfactorconv;
  if (compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()))->GetPType() == 20) {
    inEnergy = this->CompoundEnergy(compound);
    geofactor = 1.;
    sfactorconv = 1.;
  } else {
    inEnergy = this->CompoundEnergy(compound);
    geofactor = Kinematics::Channel(entrancePair).PiOverK2(this->GetCMEnergy());
    // The output's S factor column: physics::sfactor_conversion, and E sigma for a neutral pair (eta = 0).
    const azure2::physics::ChannelKinematics kin = Kinematics::Channel(entrancePair).Physics();
    sfactorconv = (kin.z1 * kin.z2 == 0) ? this->GetCMEnergy()
                                         : azure2::physics::sfactor_conversion(kin, this->GetCMEnergy());
  }
  // THM (HOES) cross sections are in arbitrary units and extend below the
  // entrance threshold, where the Gamow conversion is undefined (NaN).
  if (this->IsTHM()) sfactorconv = 1.;
  this->SetGeometricalFactor(geofactor);
  this->SetSFactorConversion(sfactorconv);
  // Electron screening: f = exp(pi eta U_e / E), so d ln f / d U_e = pi eta / E (charged particle pairs only).
  values_.screen_k = 0.0;
  if (entrancePair->GetPType() == 0 && entrancePair->GetZ(1) * entrancePair->GetZ(2) != 0 && this->GetCMEnergy() > 0.0) {
    const double eta = Kinematics::Channel(entrancePair).Eta(this->GetCMEnergy());
    values_.screen_k = pi * eta / this->GetCMEnergy();
  }

  // Every table below is rebuilt from scratch: a second call (a parent
  // point's CalcEDependentValues reaches its sub-points again, e.g. after an
  // energy shift) used to append a second copy to the sub-points' tables.
  // The channel rows are emptied in place (ClearEDependentRows), keeping their
  // storage, as RecalcEDependentValues does, or taken whole from rowsOf; the
  // THM tables are made anew.
  if (rowsOf)
    values_.edep = rowsOf->values_.edep;
  else
    this->ClearEDependentRows();
  values_.thm_ps.reset();
  values_.thm_vertex.reset();
  std::shared_ptr<ThmVertexRows> vertexRows;
  if (this->IsTHM()) vertexRows = std::make_shared<ThmVertexRows>();
  std::shared_ptr<ThmPsTable> psTable;
  bool psCoul = false;
  std::vector<ThmSpectatorWindow::Node> psNodes;
  // The Coulomb terms C_l at the window nodes depend on l, the energy and the
  // node only, not on the channel: one set per l for this point (several
  // channels of the entrance pair share an l; ThmCoulombTerm's small memo
  // cannot hold them once a window has more than 32 nodes).
  std::map<int, std::vector<complex>> psCoulByL;
  if (this->IsTHM() && settings_->thmWindow) {
    // The accepted directions at this point's energy (or, for a folding
    // sub-point beyond the reach of the window, at the nearest data point).
    settings_->thmWindow->Nodes(this->GetCMEnergy(), psNodes);
    psTable = std::make_shared<ThmPsTable>();
    psTable->nodes = (int)psNodes.size();
    psTable->offset.assign(compound->NumJGroups() + 1, 0);
    for (const ThmSpectatorWindow::Node &k : psNodes) psTable->weight.push_back(k.weight);
  }
  // Only the groups of the entrance pair's compound nucleus take part in this
  // point: a group of another nucleus in the file has thresholds on a
  // different energy axis, and its rows stay empty (every consumer skips it).
  const int nucleus = entrancePair->GetNucleus();
  std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>> ownRows;
  // ownRows: the distinct rows this point used, from the tables or computed; rows->used receives them once.
  auto channelRow = [&](int pairNum, const PPair *pair, int l, int two_j, double e) -> ChannelTable::Row {
    const ChannelTable::Key key{pairNum, l, two_j, e};
    for (const auto &kr : ownRows)
      if (kr.first == key) return kr.second;
    const ChannelTable::Row *r = rows ? rows->Find(key) : nullptr;
    ownRows.emplace_back(key, r ? *r : ChannelTable::Compute(pair, l, two_j, e, configure));
    return ownRows.back().second;
  };
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    if (psTable) psTable->offset[j - 1] = (int)psTable->slot.size();
    if (compound->GetJGroup(j)->IsInRMatrix() && compound->GetJGroup(j)->GetNucleus() == nucleus) {
      const JGroup *theJGroup = compound->GetJGroup(j);
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        const AChannel *theChannel = theJGroup->GetChannel(ch);
        const PPair *thePair = compound->GetPair(theChannel->GetPairNum());
        int lValue = theChannel->GetL();
        double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
        if (rowsOf) {
          // the table is rowsOf's (the same ChannelRowsKey): nothing to compute
        } else if (thePair->GetPType() == 0) {
          if (fabs(localEnergy) <= ScatteringLength::kThresholdTolerance) {
            // Exactly at the channel's threshold (the entrance channel of a scattering length):
            // the E -> 0 limit, P = 0 and S_l(0) from the zero-energy Coulomb solutions, in
            // place of the Whittaker log-derivative at zero binding, which is ill-conditioned.
            complex loElement(ScatteringLength::ZeroEnergyShift(thePair, lValue) - theChannel->GetBoundaryCondition(), 0.0);
            this->AddLoElement(j, ch, loElement);
            this->AddSqrtPenetrability(j, ch, 0.0);
            this->AddExpCoulombPhase(j, ch, 1.0);
            this->AddExpHardSpherePhase(j, ch, 1.0);
          } else {
            // P, S, exp(-i phi), exp(i(sigma_l - sigma_0)) of (pair, l, E_c): from the stage's channel table, or
            // computed once per (pair, l, E_c) of this point (docs/upgrade/channel_table.md).
            const ChannelTable::Row row = channelRow(theChannel->GetPairNum(), thePair, lValue, thePair->SpinOrbitTwoJ(theJGroup->GetJ()), localEnergy);
            complex loElement(row.S - theChannel->GetBoundaryCondition(), row.P);
            this->AddLoElement(j, ch, loElement);
            this->AddSqrtPenetrability(j, ch, sqrt(row.P));
            this->AddExpCoulombPhase(j, ch, row.expCP);
            this->AddExpHardSpherePhase(j, ch, row.expHSP);
          }
        } else if (thePair->GetPType() == 10) {
          complex loElement = complex(0.0, 0.0);
          this->AddLoElement(j, ch, loElement);
          double sqrtPene = (configure.paramMask & Config::USE_RMC_FORMALISM) ? 1. : pow(localEnergy / hbarc, (double)lValue + 0.5);
          this->AddSqrtPenetrability(j, ch, sqrtPene);
          this->AddExpCoulombPhase(j, ch, 1.0);
          this->AddExpHardSpherePhase(j, ch, 1.0);
        } else if (thePair->GetPType() == 20) {
          complex loElement = complex(0.0, 0.0);
          this->AddLoElement(j, ch, loElement);
          IntegratedFermiFunc fermiFunc(thePair->GetZ(1));
          double endPointE = thePair->GetSepE() - inEnergy;
          double sqrtPene = (1. + endPointE / electronMassOfFermiFunction <= 1.) ? 0. : sqrt(fermiFunc(1. + endPointE / electronMassOfFermiFunction, exitPair->GetZ(1) + exitPair->GetZ(2), thePair->GetChRad()));
          this->AddSqrtPenetrability(j, ch, sqrtPene);
          this->AddExpCoulombPhase(j, ch, 1.0);
          this->AddExpHardSpherePhase(j, ch, 1.0);
        }
        // Boundary-independent pieces of the THM entrance transfer form factor
        // M_l (mrmpy _form_factor). Stored for every channel (0 for non-entrance
        // channels) so they stay index-aligned with the penetrabilities. M_l
        // itself is assembled at the entrance vertex (THMMatrixFunc) with the
        // per-level boundary, which floats with the fit under Brune. The interior
        // LoElement above is unchanged -- the form factor replaces only the
        // entrance *vertex*, not the interior.
        double thmJl = 0.0;
        double thmRhoDjl = 0.0;
        complex thmCoul(0.0, 0.0);
        const ThmSpectatorWindow *window = this->IsTHM() ? settings_->thmWindow : nullptr;
        if (window) {
          // Spectator-momentum window (ps=..., ThmLineshape.h): the pieces at
          // every node p_k of this energy, with T_k = p_k^2/2 mu_sx added to
          // E + B (a window and spectatorEnergy for the same pair are refused
          // at startup).  Only the entrance-pair channels are stored
          // (ThmPsTable: every other channel reads 0); the single-node pieces
          // stay 0.
          if (thePair == entrancePair && thePair->GetPType() == 0) {
            psTable->slot.push_back(psTable->nSlots++);
            std::vector<complex> *coulOfL = nullptr;
            bool newL = false;
            if (configure.thm.coulombIntegral && thePair->GetZ(1) * thePair->GetZ(2) != 0) {
              auto at = psCoulByL.emplace(lValue, std::vector<complex>(psNodes.size(), complex(0.0, 0.0)));
              coulOfL = &at.first->second;
              newL = at.second;
            }
            const bool firstSlot = psTable->eps.empty();
            for (size_t k = 0; k < psNodes.size(); k++) {
              double jl = 0.0, rhoDjl = 0.0;
              complex coul(0.0, 0.0);
              double muMeV = thePair->GetRedMass() * uconv;
              double bindingE = thePair->GetBindingEnergy() + psNodes[k].es;
              if (firstSlot) psTable->eps.push_back(std::max(localEnergy + bindingE, 0.0));
              if (!(localEnergy + bindingE > 0.0)) {
                WarnThmBelowB(configure, localEnergy, bindingE);
              } else {
                ThmBesselParts(lValue, muMeV, localEnergy, bindingE, thePair->GetChRad(), jl, rhoDjl);
                if (coulOfL) {
                  if (newL)
                    (*coulOfL)[k] = ThmCoulombTerm(thePair, lValue, localEnergy,
                                                   ThmRho(muMeV, localEnergy, bindingE, 1.0));
                  coul = (*coulOfL)[k];
                }
              }
              psTable->jl.push_back(jl);
              psTable->rhodjl.push_back(rhoDjl);
              psTable->coul.push_back(coul);
              if (coul != complex(0.0, 0.0)) psCoul = true;
            }
          } else {
            psTable->slot.push_back(-1);
          }
        } else if (this->IsTHM() && thePair == entrancePair && thePair->GetPType() == 0) {
          double muMeV = thePair->GetRedMass() * uconv;
          // The spectator's kinetic energy raises the half-off-shell momentum
          // above its quasi-free value (Typel & Baur eq. 11); 0 by default.
          double bindingE = thePair->GetBindingEnergy() +
                            configure.thm.SpectatorEnergy(this->GetEntranceKey());
          if (localEnergy + bindingE > 0.0) {
            if (vertexRows) vertexRows->eps = localEnergy + bindingE;
            ThmBesselParts(lValue, muMeV, localEnergy, bindingE, thePair->GetChRad(),
                           thmJl, thmRhoDjl);
            if (configure.thm.coulombIntegral && thePair->GetZ(1) * thePair->GetZ(2) != 0)
              thmCoul = ThmCoulombTerm(thePair, lValue, localEnergy,
                                       ThmRho(muMeV, localEnergy, bindingE, 1.0));
          } else {
            WarnThmBelowB(configure, localEnergy, bindingE);
          }
        }
        // Only a THM point reads them, and only for the channels of its entrance pair (THMMatrixFunc; for any other
        // GetThmFormFactor returns 0), so only those are kept: one entry per entrance channel, not one per channel
        // of the compound nucleus.
        if (vertexRows && thePair == entrancePair) {
          vertexRows->channel.push_back(ch);
          vertexRows->jl.push_back(thmJl);
          vertexRows->rhodjl.push_back(thmRhoDjl);
          vertexRows->coul.push_back(thmCoul);
        }
      }
    }
    if (vertexRows) vertexRows->groupStart.push_back((int)vertexRows->channel.size());
  }
  if (vertexRows) {
    vertexRows->channel.shrink_to_fit();
    vertexRows->jl.shrink_to_fit();
    vertexRows->rhodjl.shrink_to_fit();
    vertexRows->coul.shrink_to_fit();
    values_.thm_vertex = std::move(vertexRows);
  }
  if (psTable) psTable->offset.back() = (int)psTable->slot.size();
  if (psTable && !psTable->slot.empty()) {  // no channel at all: no nodes, as before
    if (!psCoul) psTable->coul.clear();
    psTable->slot.shrink_to_fit();
    psTable->jl.shrink_to_fit();
    psTable->rhodjl.shrink_to_fit();
    psTable->coul.shrink_to_fit();
    psTable->weight.shrink_to_fit();
    psTable->eps.shrink_to_fit();
    values_.thm_ps = psTable;
  }
  // The rows grew channel by channel (push_back): give back the spare capacity,
  // about a quarter of the tables' memory (natural neon: 2.3 GB of 9.6 GB).  A
  // point is filled here once per energy (setup, grid re-placement); the shift
  // step refills its rows in place every evaluation and keeps their capacity.
  if (rows && rows->used) rows->used->insert(rows->used->end(), ownRows.begin(), ownRows.end());
  if (values_.edep && values_.edep.SoleOwner()) {
    EDependentTables &t = MutableEDep();
    t.ShrinkToFit();
  }
  values_.eDependentEnergy = this->GetCMEnergy();
  values_.eDependentValid = true;
}

/*!
 * Empties the per-JGroup rows of the energy-dependent arrays, keeping their storage.
 */
void Site::ClearEDependentRows() {
  if (!values_.edep) return;
  if (!values_.edep.SoleOwner()) {
    // Shared with another point: start a table of our own with the same rows,
    // empty, rather than copying values that are about to be replaced.
    auto own = SharedTable<EDependentTables>::Make();
    EDependentTables &t = own.Mutable();
    t.rowStart.reserve(values_.edep->rowStart.size());
    t.lo_elements.reserve(values_.edep->lo_elements.size());
    t.penetrabilities.reserve(values_.edep->penetrabilities.size());
    t.coulombphase.reserve(values_.edep->coulombphase.size());
    t.hardspherephase.reserve(values_.edep->hardspherephase.size());
    values_.edep = std::move(own);  // drops the shared table (release; SharedTable.h)
    return;
  }
  MutableEDep().Clear();  // our own: emptied in place, keeping its storage
}

/*!
 * Recalculates energy-dependent values using the current (possibly shifted) energy.
 * This is needed when energy shifts are applied after initialization.
 */
void Site::RecalcEDependentValues(WorkingCopy *data, const CompoundSystem *compound, const Config &configure, const ChannelRowSource *rows) {
  // Nothing to do while the energy is the one the values were computed at.
  // Every sub-point of a target-effect point used to pay for this on every
  // evaluation whenever external capture was in use (Coulomb and Whittaker
  // functions for each channel, serial): 40 s per evaluation for two
  // convolved 16N beta-delayed alpha spectra on the 12C+alpha model.
  if (EDependentCurrent()) return;
  // Clear existing energy-dependent values first.  The per-JGroup rows are
  // emptied in place rather than destroyed, so CalcEDependentValues refills the
  // same storage instead of freeing and re-allocating every row of every point
  // on every energy-shift step.  That churn, done on the pooled data sets by
  // whichever OpenMP thread happens to evaluate them, fragmented the malloc
  // arenas and made resident memory creep up for the whole fit.  The Add*
  // methods index the rows by J-group, so the values stored are unchanged.
  this->ClearEDependentRows();
  // Recalculate with current energy.
  this->CalcEDependentValues(data, compound, configure, rows);

  // Also recalculate Coulomb amplitude with current energy
  this->CalcCoulombAmplitude(data, compound);
}

/*!
 * This point's energy-dependent tables, writable.  The tables are shared
 * between copies of the point (WorkingCopy clones, sub-points that kept a node); if
 * another point still holds them, this point gets its own copy first, so a
 * write never changes what another point -- another thread's clone -- reads.
 * SharedTable::SoleOwner orders the other holders' reads before the write.
 */

Unit::EDependentTables &Site::MutableEDep() {
  if (!values_.edep)
    values_.edep = SharedTable<EDependentTables>::Make();
  else if (!values_.edep.SoleOwner())
    values_.edep = SharedTable<EDependentTables>::Make(*values_.edep);
  return values_.edep.Mutable();
}

/*!
 * This point's Legendre polynomials, writable (copy-on-write, as MutableEDep).
 */

vector_r &Site::MutableLegendre() {
  if (!values_.legendreP)
    values_.legendreP = SharedTable<vector_r>::Make();
  else if (!values_.legendreP.SoleOwner())
    values_.legendreP = SharedTable<vector_r>::Make(*values_.legendreP);
  return values_.legendreP.Mutable();
}

/*!
 * The point's energy-dependent values (CalcEDependentValues, CalcLegendreP) as
 * a share: the tables by pointer, the scalars by value.
 */

Unit::EDependentShare Site::ShareEDependentValues() const {
  EDependentShare share;
  share.energy = this->GetCMEnergy();
  share.angle = this->GetCMAngle();
  share.entranceKey = this->GetEntranceKey();
  share.exitKey = this->GetExitKey();
  share.geofactor = values_.geofactor;
  share.sfactorconv = values_.sfactorconv;
  share.screen_k = values_.screen_k;
  if (values_.eDependentValid && values_.eDependentEnergy == this->GetCMEnergy()) {  // else empty: AdoptEDependentValues refuses it
    share.tables = values_.edep;
    share.legendre = values_.legendreP;
    share.thmPs = values_.thm_ps;
    share.thmVertex = values_.thm_vertex;
  }
  return share;
}

/*!
 * Takes the values of a share made by a point at exactly the same c.m. energy,
 * c.m. angle and pairs -- what CalcEDependentValues and CalcLegendreP would have
 * computed here.  Returns false, and changes nothing, if any of them differ or
 * the share is empty.
 */

bool Site::AdoptEDependentValues(const EDependentShare &share) {
  if (!share.tables || !share.legendre || share.energy != this->GetCMEnergy() || share.angle != this->GetCMAngle() ||
      share.entranceKey != this->GetEntranceKey() || share.exitKey != this->GetExitKey())
    return false;
  values_.geofactor = share.geofactor;
  values_.sfactorconv = share.sfactorconv;
  values_.screen_k = share.screen_k;
  values_.edep = share.tables;
  values_.legendreP = share.legendre;
  values_.thm_ps = share.thmPs;
  values_.thm_vertex = share.thmVertex;
  values_.eDependentEnergy = this->GetCMEnergy();
  values_.eDependentValid = true;
  return true;
}

/*!
 * Adds an \f$ L_o \f$ matrix element with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Site::AddLoElement(int jGroupNum, int channelNum, complex loElement) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.lo_elements, t.Slot(jGroupNum, channelNum), loElement);
}

/*!
 * Adds a square root of penetrability  with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Site::AddSqrtPenetrability(int jGroupNum, int channelNum, double sqrtPene) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.penetrabilities, t.Slot(jGroupNum, channelNum), sqrtPene);
}

/*!
 * Adds an exponential of the Coulomb phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Site::AddExpCoulombPhase(int jGroupNum, int channelNum, complex expShift) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.coulombphase, t.Slot(jGroupNum, channelNum), expShift);
}

/*!
 * Adds an exponential of the hard sphere phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Site::AddExpHardSpherePhase(int jGroupNum, int channelNum, complex expShift) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.hardspherephase, t.Slot(jGroupNum, channelNum), expShift);
}

/*!
 * Calculates the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void Site::CoulombAmplitudesAt(const PPair *entrancePair, double energy, double angle, complex &direct, complex &exchange,
                               bool &spinDependent) {
  int z1 = entrancePair->GetZ(1);
  int z2 = entrancePair->GetZ(2);
  double redmass = entrancePair->GetRedMass();
  double eta = Kinematics::Channel(entrancePair).rel ? Kinematics::Channel(entrancePair).Eta(energy)
                                                     : sqrt(uconv / 2.) * fstruc * z1 * z2 * sqrt(redmass / energy);
  double cal = (1.0 / (2.0 * sqrt(pi))) * eta * (1.0 / pow(sin(angle * pi / 360.0), 2.));
  double cex = 2.0 * eta * log(sin(angle * pi / 360.0));
  complex calpha(cal * cos(cex), -cal * sin(cex));
  // Identical-particle Mott amplitude: the exchange term f_C(pi - theta).
  // Exchanging two identical particles of spin j multiplies a channel-spin
  // state |s v> by (-1)^(2j-s) and the wave function must pick up
  // (-1)^(2j), so the exchange term enters channel spin s with (-1)^s --
  // not with one boson/fermion sign for all s. Only for j = 0 (s = 0 alone)
  // do the two coincide, and there the symmetrized amplitude is returned and
  // reproduces the spin-0 Mott cross section. For j != 0 the
  // direct and exchange terms are kept apart and combined per channel spin
  // by GetCoulombAmplitude(s).
  spinDependent = false;
  exchange = complex(0., 0.);
  if (entrancePair->IsIdentical()) {
    double cal_p = (1.0 / (2.0 * sqrt(pi))) * eta * (1.0 / pow(cos(angle * pi / 360.0), 2.));
    double cex_p = 2.0 * eta * log(cos(angle * pi / 360.0));
    complex calpha_p(cal_p * cos(cex_p), -cal_p * sin(cex_p));
    if (entrancePair->HasSpinDependentExchange()) {
      spinDependent = true;
      exchange = calpha_p;
    } else {
      calpha = calpha + double(entrancePair->GetIdenticalSign()) * calpha_p;
    }
  }
  direct = calpha;
}

void Site::CalcCoulombAmplitude(WorkingCopy *data, const CompoundSystem *compound) {
  // A slot without an angle has no Coulomb amplitude: its readers hold theirs (RowSet.h).
  if (is_sub_point_ && AngleFree()) return;
  if (this->GetEntranceKey() == this->GetExitKey()) {
    const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(this->GetEntranceKey()));
    complex direct, exchange;
    CoulombAmplitudesAt(entrancePair, this->GetCMEnergy(), this->GetCMAngle(), direct, exchange, values_.spinDependentCoulomb);
    if (values_.spinDependentCoulomb) values_.coulombdirect = direct, values_.coulombexchange = exchange;
    this->SetCoulombAmplitude(direct);
  } else {
    values_.spinDependentCoulomb = false;
    this->SetCoulombAmplitude(complex(0., 0.));
  }
}

/*!
 * Sets the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void Site::SetCoulombAmplitude(complex amplitude) {
  values_.coulombamplitude = amplitude;
}
