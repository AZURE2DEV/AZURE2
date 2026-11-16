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
  double lineWindow = 0.0;           ///< > 0: capture into a broad unbound final state, E_f in E_R +- this (keV; LineWindow.h)
  const ThmWeightTable *thmWeight = nullptr;  ///< THM attachments, owned by the segment's shared pointers
  const ThmLineshape *thmLineshape = nullptr;
  const ThmSpectatorWindow *thmWindow = nullptr;
  const ThmDistortion *thmDistortion = nullptr;
  const ThmDwVertex *thmDw = nullptr;
  const ThmAngleWindow *thmAngle = nullptr;
  const ThmCoherentBackground *thmCoherent = nullptr;
};

/// Whether two weight tables of THM data (thm_weight) give the same w(E): the same table, or equal rows (ConfigThm.cpp).
bool SameThmWeight(const ThmWeightTable *a, const ThmWeightTable *b);

/// Whether the sites of two terms compute the same observable at an energy (the sharing key of the mesh, SlotGroup.h):
/// every setting the energy pass reads (GenMatrixFunc, the external capture, THM), all but those of
/// kReaderOnlyTermSettings below, which only the reader applies (Site::ValueFor).
///
/// THM (a site holds sigma_HOES at E_n without an angle, SiteObservable.h): every attachment enters the value at the
/// node, so every one is compared.  The attachments of an experiment (line shape, spectator window, distortion, DW
/// vertex, angular window, coherent background) are built once per <thm> experiment and held by each of its data sets
/// (ThmNorm.cpp): the same object is the same experiment, which is what "equal settings" means for them (a window
/// holds its experiment's data energies, a background its parameters).  The weight table is read per data set
/// (thm_weight), so two data sets of one experiment that name equal tables hold two objects: it is compared by value
/// (SameThmWeight), and such data sets share their sites.
inline bool SameObservable(const TermSettings &a, const TermSettings &b) {
  return a.entranceKey == b.entranceKey && a.exitKey == b.exitKey && a.differential == b.differential &&
         a.phase == b.phase && a.thm == b.thm && a.angularDist == b.angularDist && a.analyzingPower == b.analyzingPower &&
         a.total == b.total && a.upos == b.upos && a.j == b.j && a.l == b.l && a.phaseS == b.phaseS &&
         a.secondaryDecayL == b.secondaryDecayL && a.ic == b.ic && a.delta == b.delta &&
         a.crossSectionComponent == b.crossSectionComponent && a.lineWindow == b.lineWindow &&
         SameThmWeight(a.thmWeight, b.thmWeight) &&
         a.thmLineshape == b.thmLineshape && a.thmWindow == b.thmWindow && a.thmDistortion == b.thmDistortion &&
         a.thmDw == b.thmDw && a.thmAngle == b.thmAngle && a.thmCoherent == b.thmCoherent;
}

/// The fields of TermSettings that only the reader applies, and that SameObservable therefore does not compare.  Every
/// other field must be compared there: a field in neither place would let two different observables share their sites
/// silently (tests/consistency/repo_consistency_test.py, check 37, reads both lists).
inline constexpr const char *kReaderOnlyTermSettings[] = {"coefficientAbsolute", "maxAngDistOrder"};

#endif
