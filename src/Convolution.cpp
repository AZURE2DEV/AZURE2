#include "Convolution.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

Convolution::Factor Convolution::Factor::MakeIntrinsic(
    std::function<double(double)> fn
)
{
    Factor factor;
    factor.function = std::move(fn);
    factor.dirty = true;
    factor.scaleByDx = false;
    return factor;
}

Convolution::Factor Convolution::Factor::MakeKernel(
    std::function<double(double)> fn
)
{
    Factor factor;
    factor.function = std::move(fn);
    factor.dirty = true;
    factor.scaleByDx = true;
    return factor;
}

Convolution::Convolution(std::size_t n, Range inputRange)
    : n_(n), inputRange_(inputRange)
{
    if (n_ <= 1) {
        throw std::logic_error("Convolution: grid size must be larger than 1.");
    }

    if (inputRange_.second <= inputRange_.first) {
        throw std::logic_error("Convolution: invalid input range.");
    }

    if (n_ > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::logic_error("Convolution: grid size exceeds FFTW int limit.");
    }

    dx_ = (inputRange_.second - inputRange_.first)
        / static_cast<double>(n_ - 1);

    observed_.assign(n_, 0.0);
    combinedSpec_.assign(n_ / 2 + 1, Complex{0.0, 0.0});

    inversePlan_ = fftw_plan_dft_c2r_1d(
        static_cast<int>(n_),
        reinterpret_cast<fftw_complex*>(combinedSpec_.data()),
        observed_.data(),
        FFTW_MEASURE
    );

    if (!inversePlan_) {
        throw std::runtime_error("Convolution: failed to create inverse FFTW plan.");
    }

    maskRangeFirst_ = inputRange_.first;
    maskRangeSecond_ = inputRange_.second;
}

Convolution::~Convolution()
{
    for (auto& plan : forwardPlans_) {
        if (plan) {
            fftw_destroy_plan(plan);
        }
    }

    if (inversePlan_) {
        fftw_destroy_plan(inversePlan_);
    }
}

void Convolution::SetConvolutionInputRange(Range inputRange)
{
    if (inputRange.second <= inputRange.first) {
        throw std::logic_error("Convolution::SetConvolutionInputRange: invalid range.");
    }

    inputRange_ = inputRange;

    dx_ = (inputRange_.second - inputRange_.first)
        / static_cast<double>(n_ - 1);

    cellMaskDirty_ = true;
    maskRangeFirst_ = inputRange_.first;
    maskRangeSecond_ = inputRange_.second;

    MarkAllFactorsDirty();

    for (std::size_t m = 0; m < shiftBins_.size(); ++m) {
        if (m == 0) {
            shiftBins_[m] = 0;
        } else if (0.0 >= inputRange_.first && 0.0 <= inputRange_.second) {
            const double bins = (0.0 - inputRange_.first) / dx_;
            const auto shift = static_cast<std::ptrdiff_t>(std::llround(bins));
            const auto nSigned = static_cast<std::ptrdiff_t>(n_);

            shiftBins_[m] = (shift % nSigned + nSigned) % nSigned;
        } else {
            shiftBins_[m] = 0;
        }
    }
}

void Convolution::SetFactors(const std::vector<Factor>& factors)
{
    if (factors.empty()) {
        throw std::logic_error(
            "Convolution::SetFactors: at least one factor is required."
        );
    }

    factors_ = factors;

    for (const auto& factor : factors_) {
        if (!factor.function) {
            throw std::logic_error(
                "Convolution::SetFactors: factor has empty function."
            );
        }
    }

    const std::size_t factorCount = factors_.size();

    samples_.assign(factorCount, std::vector<double>(n_, 0.0));
    spectra_.assign(
        factorCount,
        std::vector<Complex>(n_ / 2 + 1, Complex{0.0, 0.0})
    );
    shiftBins_.assign(factorCount, 0);

    for (auto& plan : forwardPlans_) {
        if (plan) {
            fftw_destroy_plan(plan);
        }
    }

    forwardPlans_.clear();
    forwardPlans_.reserve(factorCount);

    for (std::size_t m = 0; m < factorCount; ++m) {
        fftw_plan plan = fftw_plan_dft_r2c_1d(
            static_cast<int>(n_),
            samples_[m].data(),
            reinterpret_cast<fftw_complex*>(spectra_[m].data()),
            FFTW_MEASURE
        );

        if (!plan) {
            throw std::runtime_error(
                "Convolution::SetFactors: failed to create forward FFTW plan."
            );
        }

        forwardPlans_.push_back(plan);

        if (m == 0) {
            shiftBins_[m] = 0;
            factors_[m].scaleByDx = false;
        } else if (0.0 >= inputRange_.first && 0.0 <= inputRange_.second) {
            const double bins = (0.0 - inputRange_.first) / dx_;
            const auto shift = static_cast<std::ptrdiff_t>(std::llround(bins));
            const auto nSigned = static_cast<std::ptrdiff_t>(n_);

            shiftBins_[m] = (shift % nSigned + nSigned) % nSigned;
        } else {
            shiftBins_[m] = 0;
        }

        factors_[m].dirty = true;
    }
}

void Convolution::SetCellAveragingIntervals(
    std::vector<CellAveragingInterval> intervals
)
{
    cellAveragingIntervals_ = std::move(intervals);
    cellMaskDirty_ = true;
    MarkAllFactorsDirty();
}

void Convolution::MarkFactorDirty(std::size_t index)
{
    if (index >= factors_.size()) {
        throw std::out_of_range("Convolution::MarkFactorDirty: index out of range.");
    }

    factors_[index].dirty = true;
}

void Convolution::MarkAllFactorsDirty()
{
    for (auto& factor : factors_) {
        factor.dirty = true;
    }
}

void Convolution::RecalculateCellMaskIfNeeded()
{
    if (!cellMaskDirty_
        && maskRangeFirst_ == inputRange_.first
        && maskRangeSecond_ == inputRange_.second
        && cellAveragingMask_.size() == n_) {
        return;
    }

    cellAveragingMask_.assign(n_, 0);

    for (std::size_t i = 0; i < n_; ++i) {
        const double x = inputRange_.first + static_cast<double>(i) * dx_;

        int largestNSubSamples = 0;

        for (const auto& interval : cellAveragingIntervals_) {
            const double xMin = std::get<0>(interval);
            const double xMax = std::get<1>(interval);
            const int nSub = std::get<2>(interval);

            if (x >= xMin && x < xMax) {
                largestNSubSamples = std::max(largestNSubSamples, nSub);
            }
        }

        if (largestNSubSamples <= 1) {
            largestNSubSamples = 0;
        }

        cellAveragingMask_[i] = largestNSubSamples;
    }

    maskRangeFirst_ = inputRange_.first;
    maskRangeSecond_ = inputRange_.second;
    cellMaskDirty_ = false;
}

bool Convolution::UpdateFactorSamples()
{
    if (factors_.empty()) {
        return false;
    }

    RecalculateCellMaskIfNeeded();

    std::vector<std::size_t> updated;
    updated.reserve(factors_.size());

    for (std::size_t m = 0; m < factors_.size(); ++m) {
        if (factors_[m].dirty) {
            updated.push_back(m);
        }
    }

    if (updated.empty()) {
        return false;
    }

    const int* maskPtr = cellAveragingMask_.empty()
        ? nullptr
        : cellAveragingMask_.data();

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n_; ++i) {
        const double x = inputRange_.first + static_cast<double>(i) * dx_;
        const int nSub = maskPtr ? maskPtr[i] : 0;

        for (std::size_t ui = 0; ui < updated.size(); ++ui) {
            const std::size_t m = updated[ui];

            const auto& factor = factors_[m];
            auto& output = samples_[m];

            double value = 0.0;

            if (nSub >= 2) {
                const double subDx = dx_ / static_cast<double>(nSub);
                const double x0 = x - dx_ / 2.0;

                double sum = 0.0;

                for (int k = 0; k < nSub; ++k) {
                    const double xx =
                        x0 + (static_cast<double>(k) + 0.5) * subDx;

                    sum += factor.function(xx);
                }

                value = sum / static_cast<double>(nSub);
            } else {
                value = factor.function(x);
            }

            if (factor.scaleByDx) {
                value *= dx_;
            }

            output[i] = value;
        }
    }

    for (const std::size_t m : updated) {
        fftw_execute(forwardPlans_[m]);
        factors_[m].dirty = false;
    }

    return true;
}

void Convolution::Rebuild()
{
    if (factors_.empty()) {
        return;
    }

    const std::size_t kCount = n_ / 2 + 1;
    const double pi = std::acos(-1.0);
    const double twoPiOverN = 2.0 * pi / static_cast<double>(n_);

    combinedSpec_.assign(kCount, Complex{1.0, 0.0});

    #pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < kCount; ++k) {
        Complex product{1.0, 0.0};

        for (std::size_t m = 0; m < spectra_.size(); ++m) {
            const auto& spectrum = spectra_[m];
            const std::ptrdiff_t shift = shiftBins_[m];

            const double theta =
                twoPiOverN
                * static_cast<double>(shift)
                * static_cast<double>(k);

            const Complex phase(std::cos(theta), std::sin(theta));

            product *= spectrum[k] * phase;
        }

        combinedSpec_[k] = product;
    }

    fftw_execute(inversePlan_);

    const double invN = 1.0 / static_cast<double>(n_);

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < observed_.size(); ++i) {
        observed_[i] *= invN;
    }
}

double Convolution::Evaluate(double x) const
{
    if (observed_.empty()) {
        return 0.0;
    }

    if (x <= inputRange_.first) {
        return observed_.front();
    }

    if (x >= inputRange_.second) {
        return observed_.back();
    }

    const double t = (x - inputRange_.first) / dx_;
    const std::size_t i = static_cast<std::size_t>(t);
    const double frac = t - static_cast<double>(i);

    const double y1 = observed_[i];
    const double y2 = observed_[(i + 1 < observed_.size()) ? i + 1 : i];

    return y1 + frac * (y2 - y1);
}

void Convolution::NormaliseToIntrinsicArea(std::size_t intrinsicIndex)
{
    if (intrinsicIndex >= samples_.size()) {
        throw std::out_of_range(
            "Convolution::NormaliseToIntrinsicArea: intrinsic index out of range."
        );
    }

    const double intrinsicArea = Trapezoid(samples_[intrinsicIndex]);
    const double observedArea = Trapezoid(observed_);

    if (std::abs(observedArea) > 0.0) {
        const double scale = intrinsicArea / observedArea;

        for (double& value : observed_) {
            value *= scale;
        }
    }
}

const std::vector<double>& Convolution::ObservedSamples() const
{
    return observed_;
}

const std::vector<std::vector<double>>& Convolution::RawSamples() const
{
    return samples_;
}

double Convolution::Dx() const
{
    return dx_;
}

Convolution::Range Convolution::InputRange() const
{
    return inputRange_;
}

std::size_t Convolution::Size() const
{
    return n_;
}

double Convolution::Trapezoid(const std::vector<double>& y) const
{
    double sum = 0.0;

    for (std::size_t i = 1; i < y.size(); ++i) {
        sum += 0.5 * (y[i - 1] + y[i]);
    }

    return sum * dx_;
}