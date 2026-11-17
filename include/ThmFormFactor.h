#ifndef AZURE2_THMFORMFACTOR_H
#define AZURE2_THMFORMFACTOR_H

#include <string>
#include <vector>

#include "azure2/physics/transfer/FormFactor.h"

/*!
 * The full form factor of a THM experiment (`vertexModel=ff`; docs/source/theory/thm_implementation.rst, "Full form
 * factor"): per level lambda of every J group with a channel of the entrance pair, the volume term Delta_lambda(p)
 * that the transfer vertex adds to the plane-wave surface term M_l,
 *
 *   vertex_lambda = gamma_lambda,c [ M_l(rho; B_lambda) + Delta_lambda(p) ],
 *   Delta_lambda(p) = (e_lambda - eps)/(hbar^2/2mu) Int_0^a u_lambda r j_l(p r) dr / u_lambda(a),   eps = hbar^2 p^2/2mu,
 *
 * u_lambda the Woods-Saxon level function of node count n_r with a u'(a)/u(a) = B_lambda = S_c(E_lambda) at
 * e_lambda = E_lambda - threshold (azure2/physics/transfer/FormFactor.h).
 *
 * Built once at set-up (ThmNorm.cpp, ThmBuildFormFactor) and constant: vertexModel=ff needs the energies of these
 * levels fixed (the review of M. La Cognata: resonance energies at their literature values), so the level functions
 * do not change with the fit.  Per level the overlap g(p) is a Chebyshev series on [0, p_max] (p_max from the highest
 * E + B + T_s of the experiment's points); outside it, and for a level whose energy or boundary differs from the
 * table's (not reached while the energies are fixed), Delta is computed directly from the level function.
 */
class ThmFormFactorTable {
 public:
  struct Slot {
    int j = 0, level = 0, channel = 0;  ///< 1-based
    int l = 0, nodes = 0;
    bool wildermuth = false;  ///< the node count from the Wildermuth rule (else ffNodes=)
    double Ex = 0.0;          ///< E_lambda (MeV, excitation energy of the compound nucleus) the level was solved at
    double B = 0.0;           ///< its boundary value S_c(E_lambda)
    azure2::physics::transfer::form_factor::LevelFunction f;
    azure2::physics::transfer::form_factor::Chebyshev g;  ///< Int u r j_l(p r) dr / u(a) on [0, pMax]
  };

  std::string experiment;
  double r0 = 1.25, a0 = 0.65;  ///< fm
  bool volume = true;           ///< ffVolume=off: Delta = 0
  double pMax = 0.0;            ///< fm^-1: the range of the Chebyshev series
  double threshold = 0.0;       ///< MeV: the entrance pair's threshold (S + E_x of its heavier nucleus)
  azure2::physics::transfer::form_factor::Well well;
  std::vector<Slot> slots;

  /// The slot of (J group j, level la, channel ch), or null.
  const Slot *Find(int j, int la, int ch) const;
  /*!
   * Delta of a slot at eps = E + B_xs + T_s (MeV; 0 for eps <= 0 or ffVolume=off).  `Ex` and `B` are the level's
   * energy and boundary at the evaluation; when they are not the table's, the level function is solved for them here
   * (`*resolved` is then set).
   */
  double Delta(const Slot &s, double eps, double Ex, double B, bool *resolved = nullptr) const;
  /// R = 1 + Delta/M_l of a slot at eps (for the set-up report and the tests); 1 where M_l = 0.
  double Ratio(const Slot &s, double eps) const;
};

class Config;
class CompoundSystem;
class WorkingCopy;
struct ThmExperiment;
struct ThmGroup;

/*!
 * vertexModel=ff of experiment `x` (ThmNorm.cpp, WorkingCopy::BuildThmGroups): checks the options (vertex=perlevel,
 * the Brune basis, one entrance particle pair, fixed energies of the levels the table covers), takes each level's node
 * count from ffNodes= or the Wildermuth rule, solves the level functions at the energies of <levels>, writes the table
 * to the set-up output and gives it to the group and its segments.  False (with the reason written to the output) when
 * it cannot.
 */
bool ThmBuildFormFactor(WorkingCopy &data, const Config &configure, const CompoundSystem *compound,
                        const ThmExperiment &x, ThmGroup &group, const std::string &where);

#endif
