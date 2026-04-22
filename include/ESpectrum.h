#ifndef ESPECTRUM_H
#define ESPECTRUM_H

#include "EPoint.h"
#include <vector>

class CNuc;
class Config;
class ESegment;

/*!
 * A pre-computed cross-section spectrum on a resonance-aware energy grid.
 *
 * One ESpectrum covers a unique combination of reaction channel and observable
 * type (entrance/exit pair, angle, data-type flags).  Multiple ESegments that
 * share the same combination but span different energy ranges are merged into
 * one ESpectrum whose grid covers the union of those ranges.
 *
 * Workflow:
 *   1. EData::BuildSpectra   – creates grids (called once from Initialize)
 *   2. EData::InitializeSpectra – pre-computes Lo-matrix, Legendre P, Coulomb
 *   3. EData::CalculateSpectra  – computes cross-sections (call each fit iter)
 *   4. ESpectrum::Interpolate   – fast O(log N) linear interpolation
 */
class ESpectrum {
public:

    /*!
     * Uniquely identifies a spectrum: reaction channel + observable type.
     * Angle is rounded to 3 decimal places so that nominally identical angles
     * from different segments map to the same ESpectrum.
     * For angular-distribution segments, angle is set to 0 because the full
     * Legendre expansion, not a single angle, characterises the observable.
     */
    struct Key {
        int    entranceKey;
        int    exitKey;
        double angle;           ///< CM angle in degrees (0 for ang-dist segments)
        bool   isDifferential;
        bool   isCMDifferential;
        bool   isPhase;
        bool   isAngDist;
        int    isTotalCapture;
        int    lValue;          ///< partial-wave L (meaningful for phase shifts)
        double jValue;          ///< partial-wave J (meaningful for phase shifts)
        int    maxAngDistOrder;

        bool operator==(const Key& o) const;
        bool operator<(const Key& o) const;

        /// Build a Key from segment metadata; cmAngle is the CM angle of the
        /// first data point in the segment (ignored for ang-dist segments).
        static Key FromSegment(const ESegment& seg, double cmAngle = 0.0);
    };

    explicit ESpectrum(const Key& key);

    /// Build the resonance-aware energy grid and populate the EPoint vector.
    /// Called once during EData::BuildSpectra.
    void BuildGrid(double minCMEnergy, double maxCMEnergy,
                   CNuc* compound, const Config& configure);

    /// Pre-compute Lo-matrix elements, Legendre polynomials, Coulomb amplitudes
    /// and (optionally) EC amplitudes for every grid point.
    /// Called once per run from EData::InitializeSpectra.
    void InitializePoints(CNuc* compound, const Config& configure);

    /// Compute cross-sections at all grid points using the current parameter set.
    /// Must be called after CNuc::TransformOut each fit iteration.
    void Calculate(CNuc* compound, const Config& configure);

    /// Linear interpolation of the fit cross-section at an arbitrary CM energy.
    /// Clamps to the boundary value when energy is outside the grid range.
    double Interpolate(double cmEnergy) const;

    const Key& GetKey() const { return key_; }
    int NumPoints() const     { return static_cast<int>(points_.size()); }
    EPoint*       GetPoint(int i);       ///< 1-based
    const EPoint* GetPoint(int i) const; ///< 1-based

private:
    Key               key_;
    std::vector<EPoint> points_; ///< sorted ascending by CM energy
};

#endif
