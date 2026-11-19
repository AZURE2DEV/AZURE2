#ifndef AZURE2_PHYSICS_CONSTANTS_H
#define AZURE2_PHYSICS_CONSTANTS_H

namespace azure2 {
namespace physics {

/// Physical constants, CODATA 2018 (the engine's include/Constants.h takes its values from here).
constexpr double kHbarC_MeV_fm = 197.3269804;    ///< hbar c (MeV fm), exact in the 2019 SI
constexpr double kAtomicMassUnit_MeV = 931.49410242;  ///< MeV per u
constexpr double kFineStructure = 1.00 / 137.035999084;  ///< alpha
constexpr double kElectronMass_MeV = 0.51099895;  ///< m_e c^2 (MeV); also the beta-decay Fermi function's argument
constexpr double kProtonMass_MeV = 938.27208816;   ///< m_p c^2 (MeV); the M1 capture amplitude's nuclear magneton
constexpr double kElementaryCharge_C = 1.602176634e-19;  ///< e (C), exact in the 2019 SI; per-charge yields (format 3)

}  // namespace physics
}  // namespace azure2

#endif
