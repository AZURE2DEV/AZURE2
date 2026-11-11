#ifndef AZURE2_TERMSETTINGS_H
#define AZURE2_TERMSETTINGS_H

struct ThmWeightTable;
class ThmDistortion;
class ThmDwVertex;
struct ThmAngleWindow;
struct ThmCoherentBackground;
struct ThmLineshape;
struct ThmSpectatorWindow;

/*!
 * What every calculation unit of one term of a data set computes (docs/upgrade/data_structure.md section 1.4): the
 * observable and its options, the pairs, the THM attachments.  One record per term (a segment's own term, or one of
 * its leaves), built by the segment at set-up and shared, immutable, by the units of the term and their sub-points and
 * by every copy of them (std::shared_ptr<const TermSettings>); a unit holds no copy of these settings.
 */
struct TermSettings {
  int entranceKey = 0, exitKey = 0;  ///< <pairs> keys of the term
  bool differential = false;         ///< dsigma/dOmega (lab or c.m.)
  bool phase = false;                ///< phase shift
  bool thm = false;                  ///< THM (half-off-shell)
  bool angularDist = false;          ///< Legendre coefficients (an angular distribution or coefficient data)
  bool analyzingPower = false;       ///< A_y
  bool total = false;                ///< the total cross section of the entrance pair (optical theorem)
  bool upos = false;                 ///< unobserved primary, observed secondary
  bool coefficientAbsolute = false;  ///< angular-coefficient data are A_k (b/sr), not a_k
  int maxAngDistOrder = 0;
  double j = 0.0;                    ///< total spin of a phase shift
  int l = 0;                         ///< orbital angular momentum of a phase shift
  double phaseS = -1.0;              ///< channel spin of a phase shift (< 0: none given)
  int secondaryDecayL = 0;           ///< UPOS
  double ic = 0.0, delta = 0.0;      ///< UPOS: final spin, mixing ratio
  int crossSectionComponent = 0;     ///< 0 full, 1 E1, 2 E2
  const ThmWeightTable *thmWeight = nullptr;  ///< THM attachments, owned by the segment's shared pointers
  const ThmLineshape *thmLineshape = nullptr;
  const ThmSpectatorWindow *thmWindow = nullptr;
  const ThmDistortion *thmDistortion = nullptr;
  const ThmDwVertex *thmDw = nullptr;
  const ThmAngleWindow *thmAngle = nullptr;
  const ThmCoherentBackground *thmCoherent = nullptr;
};

#endif
