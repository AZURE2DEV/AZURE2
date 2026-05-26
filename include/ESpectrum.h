#ifndef ESPECTRUM_H
#define ESPECTRUM_H

#include "Convolution.h"
#include "EPoint.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

class CNuc;
class Config;
class ESegment;

/*!
 * A pre-computed cross-section spectrum on a uniform CM energy grid.
 *
 * Each ESpectrum caches a vector of EPoints over its energy range. The
 * energy-dependent quantities that depend only on the channel structure
 * (Lo matrix, Legendre polynomials, Coulomb amplitudes) are computed
 * once during BuildGrid. Each fit iteration only needs to call
 * RecalculateGrid which performs the cheap R-matrix step per point.
 *
 * Evaluation at an arbitrary CM energy is done by linear interpolation
 * (CalculateIntrinsic). When a convolution is active, Evaluate applies it
 * on top of the interpolated intrinsic.
 */
class ESpectrum {
public:

    struct Key {
        int    entranceKey;
        int    exitKey;
        double angle;
        bool   isDifferential;
        bool   isCMDifferential;
        bool   isPhase;
        bool   isAngDist;
        int    isTotalCapture;
        int    lValue;
        double jValue;
        int    maxAngDistOrder;

        bool operator==(const Key& o) const;
        bool operator<(const Key& o) const;

        static Key FromSegment(const ESegment& seg, double cmAngle = 0.0);
    };

    explicit ESpectrum(const Key& key);

    // Custom copy semantics: grid is copied, convolution is reset (it owns
    // FFTW plans which are not safe to share / copy).
    ESpectrum(const ESpectrum& other);
    ESpectrum& operator=(const ESpectrum& other);
    ESpectrum(ESpectrum&&) noexcept = default;
    ESpectrum& operator=(ESpectrum&&) noexcept = default;

    const Key& GetKey() const { return key_; }

    void Setup(
        CNuc* compound,
        const Config& configure,
        Convolution::Range intrinsicEvaluationRange,
        std::vector<std::function<double(double)>> kernelFunctions);

    /// Allocate the grid and pre-compute everything that depends only on
    /// channel structure (NOT on the current fit parameters).
    void BuildGrid(CNuc* compound, const Config& configure, std::size_t n);

    /// Recompute the fit cross section at every grid point. Cheap; only
    /// performs the R-matrix step that depends on level parameters.
    void RecalculateGrid(CNuc* compound, const Config& configure);

    std::size_t NumGridPoints() const { return gridPoints_.size(); }
    Convolution::Range GridRange() const { return gridRange_; }

    std::size_t DetermineConvolutionGridSize(Convolution::Range intrinsicEvaluationRange) const;

    void SetMaximumConvolutionStepSize(double stepSize);
    double MaximumConvolutionStepSize() const;

    /// Interpolated intrinsic cross section at an arbitrary CM energy.
    double CalculateIntrinsic(double cmEnergy) const;

    void EnableConvolution(std::size_t n,
                           Convolution::Range inputRange,
                           std::vector<std::function<double(double)>> kernelFunctions);

    void DisableConvolution();

    bool HasConvolution() const;

    void MarkConvolutionDirty();

    void UpdateConvolution();

    double Evaluate(double cmEnergy) const;

private:

    Key key_;

    std::vector<EPoint> gridPoints_;
    Convolution::Range gridRange_{0.0, 0.0};

    std::unique_ptr<Convolution> convolution_;
    bool convolutionDirty_{true};

    // Hardcoded for now; with cell-averaged FFT sampling a small value is fine.
    double maximumConvolutionStepSize_{0.001};

};

#endif
