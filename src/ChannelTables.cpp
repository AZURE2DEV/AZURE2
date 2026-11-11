// The energy-dependent tables of a point: its channel rows (L^0 elements, square roots of the penetrabilities,
// Coulomb and hard-sphere phase factors, THM vertex pieces), their recalculation after an energy change, their
// sharing between points at the same energy, the Coulomb amplitude, and the getters the matrix functions read.
#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "AChannel.h"
#include "CNuc.h"
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

double Unit::GetGeometricalFactor() const {
  return geofactor_;
}

/*!
 * Returns the multiplicative conversion factor from cross section to s-factor.
 */

double Unit::GetSFactorConversion() const {
  return sfactorconv_;
}

/*!
 * Returns the square root of the penetrability for a channel specified by
 * positions in the JGroup and subsequent AChannel vectors.
 */

double Unit::GetSqrtPenetrability(int jGroupNum, int channelNum) const {
  return edep_->penetrabilities[edep_->At(jGroupNum, channelNum)];
}

/*!
 * Returns the THM entrance transfer form factor M_l = (b-1) j_l - rho dj_l/drho
 * (+ the external Coulomb term, when stored) for the channel specified by positions in the JGroup and AChannel vectors,
 * assembled with the given boundary b (the per-level shift function under the
 * Brune formalism, or the channel boundary constant otherwise). Returns 0 if
 * not stored (e.g. non-entrance channels, or non-THM points).
 */

complex Unit::GetThmFormFactor(int jGroupNum, int channelNum, complex boundary) const {
  if (!edep_) return 0.0;
  const EDependentTables &t = *edep_;
  if (jGroupNum - 1 >= (int)t.thm_jl.size()) return 0.0;
  if (channelNum - 1 >= (int)t.thm_jl[jGroupNum - 1].size()) return 0.0;
  return (boundary - 1.0) * t.thm_jl[jGroupNum - 1][channelNum - 1] - t.thm_rhodjl[jGroupNum - 1][channelNum - 1]
         + t.thm_coul[jGroupNum - 1][channelNum - 1];
}

complex Unit::GetThmFormFactor(int jGroupNum, int channelNum, complex boundary, int node) const {
  const ThmPsTable *t = thm_ps_.get();
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

complex Unit::GetLoElement(int jGroupNum, int channelNum) const {
  return edep_->lo_elements[edep_->At(jGroupNum, channelNum)];
}

/*!
 * Returns the factor \f$ \exp(\omega_c) \f$ where \f$ \omega_c \f$ is the Coulomb
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex Unit::GetExpCoulombPhase(int jGroupNum, int channelNum) const {
  return edep_->coulombphase[edep_->At(jGroupNum, channelNum)];
}

/*!
 * Returns the factor \f$ \exp(\delta_c) \f$ where \f$ \delta_c \f$ is the hard sphere
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex Unit::GetExpHardSpherePhase(int jGroupNum, int channelNum) const {
  return edep_->hardspherephase[edep_->At(jGroupNum, channelNum)];
}

/*!
 * Returns the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

complex Unit::GetCoulombAmplitude() const {
  return coulombamplitude_;
}

complex Unit::GetCoulombAmplitude(double s) const {
  // Every pair but an identical one with spin returns the stored amplitude
  // untouched, so the spin-0 and distinguishable results cannot move.
  if (!spinDependentCoulomb_) return coulombamplitude_;
  const int twoS = (int)std::lround(2.0 * s);
  // Channel spins of two identical particles are integers.
  const double sign = ((twoS / 2) % 2 == 0) ? 1.0 : -1.0;
  return coulombdirect_ + sign * coulombexchange_;
}

/*!
 * Sets the geometrical cross section factor \f$ \frac{\pi}{k^2} \f$.
 */

void Unit::SetGeometricalFactor(double geoFactor) {
  geofactor_ = geoFactor;
}

/*!
 * Sets the multiplicative conversion from cross section to s-factor.
 */

void Unit::SetSFactorConversion(double conversion) {
  sfactorconv_ = conversion;
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

double Unit::CompoundEnergy(CNuc *theCNuc) {
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
  if (entrancePair->GetPType() == 20) return this->GetCMEnergy() + exitPair->GetSepE() + exitPair->GetExE();
  return this->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
}

void Unit::ChannelKeys(WorkingCopy *data, CNuc *theCNuc, std::vector<ChannelTable::Key> &keys, bool includeSubPoints) {
  // The walk of CalcEDependentValues: the blocks of the entrance pair's compound nucleus, their particle channels,
  // E_c by the same expression; the exact threshold of a scattering length is not a row.
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  const double inEnergy = this->CompoundEnergy(theCNuc);
  const int nucleus = entrancePair->GetNucleus();
  for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
    JGroup *theJGroup = theCNuc->GetJGroup(j);
    if (!theJGroup->IsInRMatrix() || theJGroup->GetNucleus() != nucleus) continue;
    for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
      AChannel *theChannel = theJGroup->GetChannel(ch);
      PPair *thePair = theCNuc->GetPair(theChannel->GetPairNum());
      if (thePair->GetPType() != 0) continue;
      double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
      if (fabs(localEnergy) <= ScatteringLength::kThresholdTolerance) continue;
      keys.push_back({theChannel->GetPairNum(), theChannel->GetL(), thePair->SpinOrbitTwoJ(theJGroup->GetJ()), localEnergy});
    }
  }
  if (includeSubPoints)
    for (int i = 1; i <= this->NumSubPoints(); i++) this->GetSubPoint(*data, i)->ChannelKeys(data, theCNuc, keys, true);
}

void Unit::CalcEDependentValues(WorkingCopy *data, CNuc *theCNuc, const Config &configure, bool includeSubPoints,
                                  WorkingCopy *shareWithMappedOf, const ChannelRowSource *rows) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kEDependent);
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));

  double inEnergy;
  double geofactor;
  double sfactorconv;
  if (theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()))->GetPType() == 20) {
    inEnergy = this->CompoundEnergy(theCNuc);
    geofactor = 1.;
    sfactorconv = 1.;
  } else {
    inEnergy = this->CompoundEnergy(theCNuc);
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
  screen_k_ = 0.0;
  if (entrancePair->GetPType() == 0 && entrancePair->GetZ(1) * entrancePair->GetZ(2) != 0 && this->GetCMEnergy() > 0.0) {
    const double eta = Kinematics::Channel(entrancePair).Eta(this->GetCMEnergy());
    screen_k_ = pi * eta / this->GetCMEnergy();
  }

  // Every table below is rebuilt from scratch: a second call (a parent
  // point's CalcEDependentValues reaches its sub-points again, e.g. after an
  // energy shift) used to append a second copy to the sub-points' tables.
  // The tables, THM vertex tables included, are emptied in place
  // (ClearEDependentRows), keeping their storage, as RecalcEDependentValues does.
  this->ClearEDependentRows();
  thm_ps_.reset();
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
    psTable->offset.assign(theCNuc->NumJGroups() + 1, 0);
    for (const ThmSpectatorWindow::Node &k : psNodes) psTable->weight.push_back(k.weight);
  }
  // Only the groups of the entrance pair's compound nucleus take part in this
  // point: a group of another nucleus in the file has thresholds on a
  // different energy axis, and its rows stay empty (every consumer skips it).
  const int nucleus = entrancePair->GetNucleus();
  std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>> ownRows;
  // ownRows: the distinct rows this point used, from the tables or computed; rows->used receives them once.
  auto channelRow = [&](int pairNum, PPair *pair, int l, int two_j, double e) -> ChannelTable::Row {
    const ChannelTable::Key key{pairNum, l, two_j, e};
    for (const auto &kr : ownRows)
      if (kr.first == key) return kr.second;
    const ChannelTable::Row *r = rows ? rows->Find(key) : nullptr;
    ownRows.emplace_back(key, r ? *r : ChannelTable::Compute(pair, l, two_j, e, configure));
    return ownRows.back().second;
  };
  for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
    if (psTable) psTable->offset[j - 1] = (int)psTable->slot.size();
    if (theCNuc->GetJGroup(j)->IsInRMatrix() && theCNuc->GetJGroup(j)->GetNucleus() == nucleus) {
      JGroup *theJGroup = theCNuc->GetJGroup(j);
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        AChannel *theChannel = theJGroup->GetChannel(ch);
        PPair *thePair = theCNuc->GetPair(theChannel->GetPairNum());
        int lValue = theChannel->GetL();
        double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
        if (thePair->GetPType() == 0) {
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
            for (size_t k = 0; k < psNodes.size(); k++) {
              double jl = 0.0, rhoDjl = 0.0;
              complex coul(0.0, 0.0);
              double muMeV = thePair->GetRedMass() * uconv;
              double bindingE = thePair->GetBindingEnergy() + psNodes[k].es;
              if (!(localEnergy + bindingE > 0.0)) {
                WarnThmBelowB(configure, localEnergy, bindingE);
              } else {
                ThmBesselParts(lValue, muMeV, localEnergy, bindingE, thePair->GetChRad(), jl, rhoDjl);
                if (coulOfL) {
                  if (newL)
                    (*coulOfL)[k] = ThmCoulombTerm(thePair, lValue, localEnergy,
                                                   ThmRho(muMeV, localEnergy, bindingE, 1.0),
                                                   !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
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
            ThmBesselParts(lValue, muMeV, localEnergy, bindingE, thePair->GetChRad(),
                           thmJl, thmRhoDjl);
            if (configure.thm.coulombIntegral && thePair->GetZ(1) * thePair->GetZ(2) != 0)
              thmCoul = ThmCoulombTerm(thePair, lValue, localEnergy,
                                       ThmRho(muMeV, localEnergy, bindingE, 1.0),
                                       !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
          } else {
            WarnThmBelowB(configure, localEnergy, bindingE);
          }
        }
        // Only a THM point reads them (GetThmFormFactor returns 0 for a
        // missing entry): a classic point storing three zeros per channel cost
        // 0.6 GB with dev's hard-sphere groups (includeAllJPi=1,
        // tests/convolution_unify, 430k sub-points).
        if (this->IsTHM()) this->AddThmFormFactor(j, ch, thmJl, thmRhoDjl, thmCoul);
      }
    }
  }
  if (psTable) psTable->offset.back() = (int)psTable->slot.size();
  if (psTable && !psTable->slot.empty()) {  // no channel at all: no nodes, as before
    if (!psCoul) psTable->coul.clear();
    psTable->slot.shrink_to_fit();
    psTable->jl.shrink_to_fit();
    psTable->rhodjl.shrink_to_fit();
    psTable->coul.shrink_to_fit();
    psTable->weight.shrink_to_fit();
    thm_ps_ = psTable;
  }
  // The rows grew channel by channel (push_back): give back the spare capacity,
  // about a quarter of the tables' memory (natural neon: 2.3 GB of 9.6 GB).  A
  // point is filled here once per energy (setup, grid re-placement); the shift
  // step refills its rows in place every evaluation and keeps their capacity.
  if (rows && rows->used) rows->used->insert(rows->used->end(), ownRows.begin(), ownRows.end());
  if (edep_ && edep_.SoleOwner()) {
    EDependentTables &t = MutableEDep();
    t.ShrinkToFit();  // the THM vertex tables too
  }
  if (includeSubPoints) {
    for (int i = 1; i <= this->NumSubPoints(); i++) {
      this->GetSubPoint(*data, i)->CalcEDependentValues(data, theCNuc, configure, true, nullptr, rows);
    }
  }
  // The points mapped onto this one share its tables -- at setup and on a grid
  // re-placement, where nothing else writes them meanwhile.  Not on an energy
  // shift (RecalcEDependentValues, shareWithMapped false): there every mapped
  // point is re-energised by its own segment, a point mapped from another
  // segment possibly on another thread at the same moment
  // (WorkingCopy::FillEnergyShiftsFromParams), and two threads assigning the same
  // shared_ptr freed the tables under the other's copy-on-write
  // (heap-use-after-free, AddressSanitizer 2026-11-06; macOS CI aborts).
  const azure2::MappedPoints mappedPts = this->PlanMapped(shareWithMappedOf);
  for (int i = 1; i <= mappedPts.size(); i++) {
    Unit *mappedPoint = mappedPts[i - 1];
    if (mappedPoint->OwnsEDependentValues()) continue;
    mappedPoint->geofactor_ = geofactor_;
    mappedPoint->sfactorconv_ = sfactorconv_;
    mappedPoint->screen_k_ = screen_k_;
    mappedPoint->edep_ = edep_;  // shared; either point copies before it writes
    mappedPoint->thm_ps_ = thm_ps_;
    if (!includeSubPoints) continue;
    for (int ii = 1; ii <= this->NumSubPoints(); ii++) {
      Unit *subMappedPoint = mappedPoint->GetSubPoint(*data, ii);
      subMappedPoint->geofactor_ = this->GetSubPoint(*data, ii)->geofactor_;
      subMappedPoint->sfactorconv_ = this->GetSubPoint(*data, ii)->sfactorconv_;
      subMappedPoint->screen_k_ = this->GetSubPoint(*data, ii)->screen_k_;
      subMappedPoint->edep_ = this->GetSubPoint(*data, ii)->edep_;
      subMappedPoint->thm_ps_ = this->GetSubPoint(*data, ii)->thm_ps_;
    }
  }
  eDependentEnergy_ = this->GetCMEnergy();
  eDependentValid_ = true;
}

/*!
 * Empties the per-JGroup rows of the energy-dependent arrays, keeping their storage.
 */
void Unit::ClearEDependentRows() {
  if (!edep_) return;
  if (!edep_.SoleOwner()) {
    // Shared with another point: start a table of our own with the same rows,
    // empty, rather than copying values that are about to be replaced.
    auto own = SharedTable<EDependentTables>::Make();
    EDependentTables &t = own.Mutable();
    t.rowStart.reserve(edep_->rowStart.size());
    t.lo_elements.reserve(edep_->lo_elements.size());
    t.penetrabilities.reserve(edep_->penetrabilities.size());
    t.coulombphase.reserve(edep_->coulombphase.size());
    t.hardspherephase.reserve(edep_->hardspherephase.size());
    edep_ = std::move(own);  // drops the shared table (release; SharedTable.h)
    return;
  }
  MutableEDep().Clear();  // our own: emptied in place, keeping its storage
}

/*!
 * Recalculates energy-dependent values using the current (possibly shifted) energy.
 * This is needed when energy shifts are applied after initialization.
 */
void Unit::RecalcEDependentValues(WorkingCopy *data, CNuc *theCNuc, const Config &configure, bool includeSubPoints,
                                    const ChannelRowSource *rows) {
  // Nothing to do while the energy is the one the values were computed at.
  // Every sub-point of a target-effect point used to pay for this on every
  // evaluation whenever external capture was in use (Coulomb and Whittaker
  // functions for each channel, serial): 40 s per evaluation for two
  // convolved 16N beta-delayed alpha spectra on the 12C+alpha model.
  if (eDependentValid_ && this->GetCMEnergy() == eDependentEnergy_) return;
  // Clear existing energy-dependent values first.  The per-JGroup rows are
  // emptied in place rather than destroyed, so CalcEDependentValues refills the
  // same storage instead of freeing and re-allocating every row of every point
  // on every energy-shift step.  That churn, done on the pooled data sets by
  // whichever OpenMP thread happens to evaluate them, fragmented the malloc
  // arenas and made resident memory creep up for the whole fit.  The Add*
  // methods index the rows by J-group, so the values stored are unchanged.
  this->ClearEDependentRows();
  // CalcEDependentValues also refills the sub-points (appending to their rows),
  // so empty theirs too: otherwise every call doubled their rows, and the
  // doubled storage would now be kept.  The callers re-energise and recalculate
  // each sub-point afterwards in any case.
  //
  // Segment::UpdatePointEnergiesWithShift passes includeSubPoints = false: it
  // moves each sub-point to its shifted energy and recalculates it right
  // after, so recomputing the sub-points here, still at their OLD energies,
  // only doubled the cost of every shift step (all Coulomb functions of every
  // sub-point, twice, serially).
  if (includeSubPoints)
    for (int i = 1; i <= this->NumSubPoints(); i++) this->GetSubPoint(*data, i)->ClearEDependentRows();

  // Recalculate with current energy.  The points mapped onto this one are not
  // touched: each is re-energised itself, from its host's values
  // (CopyEDependentValuesFrom) or anew (Segment::UpdatePointEnergyWithShift).
  this->CalcEDependentValues(data, theCNuc, configure, includeSubPoints, nullptr, rows);

  // Also recalculate Coulomb amplitude with current energy
  this->CalcCoulombAmplitude(data, theCNuc, includeSubPoints);
}

/*!
 * This point's energy-dependent tables, writable.  The tables are shared
 * between copies of the point (WorkingCopy clones, the points mapped onto it); if
 * another point still holds them, this point gets its own copy first, so a
 * write never changes what another point -- another thread's clone -- reads.
 * SharedTable::SoleOwner orders the other holders' reads before the write.
 */

Unit::EDependentTables &Unit::MutableEDep() {
  if (!edep_)
    edep_ = SharedTable<EDependentTables>::Make();
  else if (!edep_.SoleOwner())
    edep_ = SharedTable<EDependentTables>::Make(*edep_);
  return edep_.Mutable();
}

/*!
 * This point's Legendre polynomials, writable (copy-on-write, as MutableEDep).
 */

vector_r &Unit::MutableLegendre() {
  if (!legendreP_)
    legendreP_ = SharedTable<vector_r>::Make();
  else if (!legendreP_.SoleOwner())
    legendreP_ = SharedTable<vector_r>::Make(*legendreP_);
  return legendreP_.Mutable();
}

/*!
 * The point's energy-dependent values (CalcEDependentValues, CalcLegendreP) as
 * a share: the tables by pointer, the scalars by value.
 */

Unit::EDependentShare Unit::ShareEDependentValues() const {
  EDependentShare share;
  share.energy = this->GetCMEnergy();
  share.angle = this->GetCMAngle();
  share.entranceKey = this->GetEntranceKey();
  share.exitKey = this->GetExitKey();
  share.geofactor = geofactor_;
  share.sfactorconv = sfactorconv_;
  share.screen_k = screen_k_;
  if (eDependentValid_ && eDependentEnergy_ == this->GetCMEnergy()) {  // else empty: AdoptEDependentValues refuses it
    share.tables = edep_;
    share.legendre = legendreP_;
    share.thmPs = thm_ps_;
  }
  return share;
}

/*!
 * Takes the values of a share made by a point at exactly the same c.m. energy,
 * c.m. angle and pairs -- what CalcEDependentValues and CalcLegendreP would have
 * computed here.  Returns false, and changes nothing, if any of them differ or
 * the share is empty.
 */

bool Unit::AdoptEDependentValues(const EDependentShare &share) {
  if (!share.tables || !share.legendre || share.energy != this->GetCMEnergy() || share.angle != this->GetCMAngle() ||
      share.entranceKey != this->GetEntranceKey() || share.exitKey != this->GetExitKey())
    return false;
  geofactor_ = share.geofactor;
  sfactorconv_ = share.sfactorconv;
  screen_k_ = share.screen_k;
  edep_ = share.tables;
  legendreP_ = share.legendre;
  thm_ps_ = share.thmPs;
  eDependentEnergy_ = this->GetCMEnergy();
  eDependentValid_ = true;
  return true;
}

/*!
 * Copies the energy-dependent values (geometrical factor, S-factor conversion,
 * L_o elements, penetrabilities, Coulomb and hard-sphere phases) from another
 * point at exactly the same c.m. energy and with the same entrance pair --
 * a host of point mapping (WorkingCopy::MapData) or its sub-point.  Those values do
 * not depend on the angle; the Coulomb amplitude does, and is not copied.
 * Returns false, and copies nothing, if the energies differ.
 */
bool Unit::CopyEDependentValuesFrom(const Unit &other) {
  if (other.GetCMEnergy() != this->GetCMEnergy() || other.GetEntranceKey() != this->GetEntranceKey() ||
      other.GetExitKey() != this->GetExitKey() || !other.eDependentValid_ ||
      other.eDependentEnergy_ != other.GetCMEnergy())
    return false;
  geofactor_ = other.geofactor_;
  sfactorconv_ = other.sfactorconv_;
  screen_k_ = other.screen_k_;
  edep_ = other.edep_;  // shared; either point copies before it writes
  thm_ps_ = other.thm_ps_;  // immutable, shared like the tables
  eDependentEnergy_ = other.eDependentEnergy_;
  eDependentValid_ = true;
  return true;
}

/*!
 * Adds an \f$ L_o \f$ matrix element with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Unit::AddLoElement(int jGroupNum, int channelNum, complex loElement) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.lo_elements, t.Slot(jGroupNum, channelNum), loElement);
}

/*!
 * Adds a square root of penetrability  with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Unit::AddSqrtPenetrability(int jGroupNum, int channelNum, double sqrtPene) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.penetrabilities, t.Slot(jGroupNum, channelNum), sqrtPene);
}

/*!
 * Adds the boundary-independent pieces j_l(rho) and rho dj_l/drho of the THM
 * entrance transfer form factor with reference to positions in the JGroup and
 * subsequent AChannel vectors. Filled in lockstep with the penetrabilities
 * (one entry per channel), so a non-entrance channel stores 0.
 */

void Unit::AddThmFormFactor(int jGroupNum, int channelNum, double jl, double rhoDjl, complex coul) {
  vector_r d;
  vector_c dc;
  EDependentTables &t = MutableEDep();
  while (jGroupNum > t.thm_jl.size()) t.thm_jl.push_back(d);
  while (jGroupNum > t.thm_rhodjl.size()) t.thm_rhodjl.push_back(d);
  while (jGroupNum > t.thm_coul.size()) t.thm_coul.push_back(dc);
  t.thm_jl[jGroupNum - 1].push_back(jl);
  t.thm_rhodjl[jGroupNum - 1].push_back(rhoDjl);
  t.thm_coul[jGroupNum - 1].push_back(coul);
  assert(channelNum == (int)t.thm_jl[jGroupNum - 1].size());
}

/*!
 * Adds an exponential of the Coulomb phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Unit::AddExpCoulombPhase(int jGroupNum, int channelNum, complex expShift) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.coulombphase, t.Slot(jGroupNum, channelNum), expShift);
}

/*!
 * Adds an exponential of the hard sphere phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void Unit::AddExpHardSpherePhase(int jGroupNum, int channelNum, complex expShift) {
  EDependentTables &t = MutableEDep();
  EDependentTables::Put(t.hardspherephase, t.Slot(jGroupNum, channelNum), expShift);
}

/*!
 * Calculates the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void Unit::CalcCoulombAmplitude(WorkingCopy *data, CNuc *theCNuc, bool includeSubPoints) {
  if (this->GetEntranceKey() == this->GetExitKey()) {
    PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    int z1 = entrancePair->GetZ(1);
    int z2 = entrancePair->GetZ(2);
    double redmass = entrancePair->GetRedMass();
    double energy = this->GetCMEnergy();
    double angle = this->GetCMAngle();
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
    // do the two coincide, and there the symmetrized amplitude is stored as
    // before and reproduces the spin-0 Mott cross section. For j != 0 the
    // direct and exchange terms are kept apart and combined per channel spin
    // by GetCoulombAmplitude(s).
    spinDependentCoulomb_ = false;
    if (entrancePair->IsIdentical()) {
      double cal_p = (1.0 / (2.0 * sqrt(pi))) * eta * (1.0 / pow(cos(angle * pi / 360.0), 2.));
      double cex_p = 2.0 * eta * log(cos(angle * pi / 360.0));
      complex calpha_p(cal_p * cos(cex_p), -cal_p * sin(cex_p));
      if (entrancePair->HasSpinDependentExchange()) {
        spinDependentCoulomb_ = true;
        coulombdirect_ = calpha;
        coulombexchange_ = calpha_p;
      } else {
        calpha = calpha + double(entrancePair->GetIdenticalSign()) * calpha_p;
      }
    }
    this->SetCoulombAmplitude(calpha);
  } else {
    spinDependentCoulomb_ = false;
    this->SetCoulombAmplitude(complex(0., 0.));
  }
  if (!includeSubPoints) return;
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(*data, i)->CalcCoulombAmplitude(data, theCNuc);
  }
}

/*!
 * Sets the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void Unit::SetCoulombAmplitude(complex amplitude) {
  coulombamplitude_ = amplitude;
}
