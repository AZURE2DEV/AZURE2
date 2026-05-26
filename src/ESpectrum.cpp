#include "ESpectrum.h"
#include "ESegment.h"
#include "CNuc.h"
#include "Config.h"
#include <cmath>

// ─── Key ─────────────────────────────────────────────────────────────────────

bool ESpectrum::Key::operator==(const Key& o) const {
    return entranceKey     == o.entranceKey
        && exitKey         == o.exitKey
        && angle           == o.angle
        && isDifferential  == o.isDifferential
        && isCMDifferential== o.isCMDifferential
        && isPhase         == o.isPhase
        && isAngDist       == o.isAngDist
        && isTotalCapture  == o.isTotalCapture
        && lValue          == o.lValue
        && jValue          == o.jValue
        && maxAngDistOrder == o.maxAngDistOrder;
}

bool ESpectrum::Key::operator<(const Key& o) const {
    if (entranceKey      != o.entranceKey)      return entranceKey      < o.entranceKey;
    if (exitKey          != o.exitKey)          return exitKey          < o.exitKey;
    if (angle            != o.angle)            return angle            < o.angle;
    if (isDifferential   != o.isDifferential)   return isDifferential   < o.isDifferential;
    if (isCMDifferential != o.isCMDifferential) return isCMDifferential < o.isCMDifferential;
    if (isPhase          != o.isPhase)          return isPhase          < o.isPhase;
    if (isAngDist        != o.isAngDist)        return isAngDist        < o.isAngDist;
    if (isTotalCapture   != o.isTotalCapture)   return isTotalCapture   < o.isTotalCapture;
    if (lValue           != o.lValue)           return lValue           < o.lValue;
    if (jValue           != o.jValue)           return jValue           < o.jValue;
    return maxAngDistOrder < o.maxAngDistOrder;
}

ESpectrum::Key ESpectrum::Key::FromSegment(const ESegment& seg, double cmAngle) {
    Key k;
    k.entranceKey     = seg.GetEntranceKey();
    k.exitKey         = seg.GetExitKey();
    k.isDifferential  = seg.IsDifferential();
    k.isCMDifferential= seg.IsCMDifferential();
    k.isPhase         = seg.IsPhase();
    k.isAngDist       = seg.IsAngularDist();
    k.isTotalCapture  = seg.IsTotalCapture();
    k.lValue          = seg.GetL();
    k.jValue          = seg.GetJ();
    k.maxAngDistOrder = seg.GetMaxAngDistOrder();
    k.angle = seg.IsAngularDist() ? 0.0
                                  : std::round(cmAngle * 1000.0) / 1000.0;
    return k;
}

// ─── ESpectrum ───────────────────────────────────────────────────────────────

ESpectrum::ESpectrum(const Key& key) : key_(key) {}

ESpectrum::ESpectrum(const ESpectrum& other)
    : key_(other.key_),
      gridPoints_(other.gridPoints_),
      gridRange_(other.gridRange_),
      convolution_(nullptr),
      convolutionDirty_(true),
      maximumConvolutionStepSize_(other.maximumConvolutionStepSize_)
{}

ESpectrum& ESpectrum::operator=(const ESpectrum& other) {
    if (this != &other) {
        key_                        = other.key_;
        gridPoints_                 = other.gridPoints_;
        gridRange_                  = other.gridRange_;
        convolution_.reset();
        convolutionDirty_           = true;
        maximumConvolutionStepSize_ = other.maximumConvolutionStepSize_;
    }
    return *this;
}

//------------------------------------------------
void ESpectrum::Setup(
    CNuc* compound,
    const Config& configure,
    Convolution::Range intrinsicEvaluationRange,
    std::vector<std::function<double(double)>> kernelFunctions)
{
    const std::size_t n = DetermineConvolutionGridSize(intrinsicEvaluationRange);

    gridRange_ = intrinsicEvaluationRange;
    BuildGrid(compound, configure, n);

    if (kernelFunctions.empty()) {
        DisableConvolution();
        return;
    }

    EnableConvolution(
        n,
        intrinsicEvaluationRange,
        std::move(kernelFunctions)
    );
}

//------------------------------------------------
std::size_t ESpectrum::DetermineConvolutionGridSize(
    Convolution::Range intrinsicEvaluationRange) const
{
    const double range =
        intrinsicEvaluationRange.second - intrinsicEvaluationRange.first;

    if (range <= 0.0) {
        return 2;
    }

    std::size_t n = 2;

    while (range / static_cast<double>(n - 1) > maximumConvolutionStepSize_) {
        n *= 2;
    }

    return n;
}

//------------------------------------------------
void ESpectrum::SetMaximumConvolutionStepSize(double stepSize) {
    maximumConvolutionStepSize_ = stepSize;
}

double ESpectrum::MaximumConvolutionStepSize() const {
    return maximumConvolutionStepSize_;
}

//------------------------------------------------
void ESpectrum::BuildGrid(CNuc* compound, const Config& configure, std::size_t n)
{
    gridPoints_.clear();
    if (!compound || n < 2) return;

    const double eMin = gridRange_.first;
    const double eMax = gridRange_.second;
    if (!(eMax > eMin)) return;

    gridPoints_.reserve(n);

    const double dx = (eMax - eMin) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        const double e = eMin + dx * static_cast<double>(i);
        gridPoints_.emplace_back(key_.angle, e,
                                 key_.entranceKey, key_.exitKey,
                                 key_.isDifferential, key_.isPhase,
                                 key_.isAngDist, key_.jValue,
                                 key_.lValue, key_.maxAngDistOrder);
    }

    for (auto& point : gridPoints_) {
        point.CalcEDependentValues(compound, configure);
        point.CalcLegendreP(Config::maxLOrder, compound, nullptr);
        point.CalcCoulombAmplitude(compound);
        if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
            point.CalculateECAmplitudes(compound, configure);
        }
    }
}

//------------------------------------------------
void ESpectrum::RecalculateGrid(CNuc* compound, const Config& configure)
{
    if (!compound) return;
    for (auto& point : gridPoints_) {
        point.Calculate(compound, configure);
    }
}

//------------------------------------------------
double ESpectrum::CalculateIntrinsic(double cmEnergy) const
{
    if (gridPoints_.empty()) return 0.0;
    if (gridPoints_.size() == 1) return gridPoints_.front().GetFitCrossSection();

    const double eMin = gridRange_.first;
    const double eMax = gridRange_.second;
    const std::size_t n = gridPoints_.size();

    if (cmEnergy <= eMin) return gridPoints_.front().GetFitCrossSection();
    if (cmEnergy >= eMax) return gridPoints_.back().GetFitCrossSection();

    const double dx = (eMax - eMin) / static_cast<double>(n - 1);
    const double fIndex = (cmEnergy - eMin) / dx;

    std::size_t i = static_cast<std::size_t>(fIndex);
    if (i >= n - 1) i = n - 2;

    const double t = fIndex - static_cast<double>(i);
    const double a = gridPoints_[i].GetFitCrossSection();
    const double b = gridPoints_[i + 1].GetFitCrossSection();

    return a + t * (b - a);
}

//------------------------------------------------
void ESpectrum::EnableConvolution(
    std::size_t n,
    Convolution::Range inputRange,
    std::vector<std::function<double(double)>> kernelFunctions)
{
    convolution_ = std::make_unique<Convolution>(n, inputRange);

    std::vector<Convolution::Factor> factors;
    factors.reserve(kernelFunctions.size() + 1);

    factors.push_back(
        Convolution::Factor::MakeIntrinsic(
            [this](double e) {
                return this->CalculateIntrinsic(e);
            }
        )
    );

    for (auto& kernelFunction : kernelFunctions) {
        factors.push_back(
            Convolution::Factor::MakeKernel(std::move(kernelFunction))
        );
    }

    convolution_->SetFactors(factors);

    MarkConvolutionDirty();
}

//------------------------------------------------
void ESpectrum::DisableConvolution()
{
    convolution_.reset();
    convolutionDirty_ = false;
}

//------------------------------------------------
bool ESpectrum::HasConvolution() const
{
    return static_cast<bool>(convolution_);
}

//------------------------------------------------
void ESpectrum::MarkConvolutionDirty()
{
    convolutionDirty_ = true;
}

//------------------------------------------------
void ESpectrum::UpdateConvolution()
{
    if (!convolution_) {
        return;
    }

    if (!convolutionDirty_) {
        return;
    }

    convolution_->MarkFactorDirty(0);

    if (convolution_->UpdateFactorSamples()) {
        convolution_->Rebuild();
        convolution_->NormaliseToIntrinsicArea();
    }

    convolutionDirty_ = false;
}

//------------------------------------------------
double ESpectrum::Evaluate(double cmEnergy) const
{
    if (convolution_) {
        return convolution_->Evaluate(cmEnergy);
    }

    return CalculateIntrinsic(cmEnergy);
}
