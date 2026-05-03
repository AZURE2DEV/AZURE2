#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <complex>
#include <cstddef>
#include <fftw3.h>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

class Convolution
{
public:
    using Complex = std::complex<double>;
    using Range = std::pair<double, double>;
    using CellAveragingInterval = std::tuple<double, double, int>;

    struct Factor
    {
        std::function<double(double)> function;
        bool dirty = true;
        bool scaleByDx = true;

        static Factor MakeIntrinsic(std::function<double(double)> fn);
        static Factor MakeKernel(std::function<double(double)> fn);
    };

    Convolution(std::size_t n, Range inputRange);
    ~Convolution();

    Convolution(const Convolution&) = delete;
    Convolution& operator=(const Convolution&) = delete;

    void SetFactors(const std::vector<Factor>& factors);
    void SetConvolutionInputRange(Range inputRange);

    void SetCellAveragingIntervals(std::vector<CellAveragingInterval> intervals);

    void MarkFactorDirty(std::size_t index);
    void MarkAllFactorsDirty();

    bool UpdateFactorSamples();
    void Rebuild();

    double Evaluate(double x) const;
    void NormaliseToIntrinsicArea(std::size_t intrinsicIndex = 0);

    const std::vector<double>& ObservedSamples() const;
    const std::vector<std::vector<double>>& RawSamples() const;

    double Dx() const;
    Range InputRange() const;
    std::size_t Size() const;

private:
    void RecalculateCellMaskIfNeeded();
    double Trapezoid(const std::vector<double>& y) const;

    std::size_t n_{0};
    Range inputRange_{0.0, 0.0};
    double dx_{0.0};

    std::vector<Factor> factors_;

    std::vector<std::vector<double>> samples_;
    std::vector<std::vector<Complex>> spectra_;
    std::vector<double> observed_;
    std::vector<Complex> combinedSpec_;

    std::vector<fftw_plan> forwardPlans_;
    fftw_plan inversePlan_{nullptr};

    std::vector<std::ptrdiff_t> shiftBins_;

    std::vector<int> cellAveragingMask_;
    std::vector<CellAveragingInterval> cellAveragingIntervals_;
    bool cellMaskDirty_{true};
    double maskRangeFirst_{0.0};
    double maskRangeSecond_{0.0};
};

#endif