#ifndef AZURE2_THMGROUP_H
#define AZURE2_THMGROUP_H

#include <memory>
#include <string>
#include <vector>

struct ThmLineshape;
struct ThmSpectatorWindow;
struct ThmAngleWindow;
struct ThmCoherentBackground;
class ThmDistortion;
class ThmDwVertex;

/*!
 * A THM experiment (<thm> experiment[<name>] ..., ThmExperiment.h) as the
 * data see it: the segments that share one profiled norm and background.
 */
struct ThmGroup {
  std::string name;
  int terms = 0;              ///< background terms, 0..3
  std::vector<int> segments;  ///< 1-based indices into the segments, ascending
  /// One segment and no background: exactly the per-segment profile of
  /// ESegment::ProfileNormChiSquared, which then handles it.
  bool trivial = false;
  /// Coulomb line shape (lineshape=on), or null; shared with its segments.
  std::shared_ptr<const ThmLineshape> lineshape;
  /// Spectator-momentum window (ps=hulthen|gauss|table), or null; shared
  /// with its segments.
  std::shared_ptr<const ThmSpectatorWindow> window;
  /// Distortion factor R(E) (distortion=...), or null; shared with its segments.
  std::shared_ptr<const ThmDistortion> distortion;
  /// Distorted-wave entrance vertex (vertexModel=dw), or null; shared with
  /// its segments.  It replaces R(E): distortion is then null.
  std::shared_ptr<const ThmDwVertex> dwVertex;
  /// Angular window of the exit pair (theta=thmin-thmax), or null; shared
  /// with its segments.
  std::shared_ptr<const ThmAngleWindow> angle;
  /// Entrance pair key of its segments (0 if they differ).
  int pairKey = 0;
  /// Coherent background (cbackground=), or null; shared with its segments.
  std::shared_ptr<const ThmCoherentBackground> coherent;
};
/*!
 * The fit parameters of the THM coherent backgrounds (cbackground=), the
 * last block of the parameter vector after the energy shifts (none without
 * the key): Re c0, Im c0 [, Re c1, Im c1] of each combination, in the order
 * of the experiments and their terms.  `value` is the current value.
 */
struct ThmCoherentParam {
  std::string name;
  std::string experiment;
  double value = 0.0;
  bool fixed = false;
  int jGroup = -1, entrance = -1, exit = -1;  ///< the combination's J group and channels in it (1-based)
  int term = -1;                              ///< 0 Re c0, 1 Im c0, 2 Re c1, 3 Im c1
};

#endif
