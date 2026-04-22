#include "ESpectrum.h"
#include "ESegment.h"
#include "CNuc.h"
#include "Config.h"
#include "AdaptiveIntegrationGrid.h"
#include <algorithm>
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
    // Angular-distribution observables are characterised by the full Legendre
    // expansion, not a single angle; collapse all to angle = 0.
    k.angle = seg.IsAngularDist() ? 0.0
                                  : std::round(cmAngle * 1000.0) / 1000.0;
    return k;
}

// ─── ESpectrum ───────────────────────────────────────────────────────────────

ESpectrum::ESpectrum(const Key& key) : key_(key) {}

void ESpectrum::BuildGrid(double minCMEnergy, double maxCMEnergy,
                          CNuc* compound, const Config& configure) {
    points_.clear();

    // 2 % margin so resonances sitting exactly at a segment boundary are
    // still covered by the fine grid on both sides.
    const double range  = maxCMEnergy - minCMEnergy;
    const double margin = range > 0.0 ? 0.02 * range : 0.001;
    const double lo     = std::max(minCMEnergy - margin, 1.0e-4);
    const double hi     = maxCMEnergy + margin;

    AdaptiveIntegrationGrid::GridConfig cfg;
    cfg.entranceKey           = key_.entranceKey;
    cfg.maxPoints             = 2000;
    cfg.baseEnergyStep        = range > 0.0 ? range / 100.0 : 0.01;
    cfg.resonanceWidthMultiplier = 5.0;
    cfg.pointsPerWidth        = 20.0;

    AdaptiveIntegrationGrid grid(cfg);
    // GenerateGrid returns energies in descending order.
    std::vector<double> energies = grid.GenerateGrid(hi, lo, compound);

    std::sort(energies.begin(), energies.end());
    energies.erase(std::unique(energies.begin(), energies.end()), energies.end());

    points_.reserve(energies.size());
    for (double e : energies) {
        points_.emplace_back(key_.angle, e,
                             key_.entranceKey, key_.exitKey,
                             key_.isDifferential, key_.isPhase,
                             key_.isAngDist, key_.jValue,
                             key_.lValue, key_.maxAngDistOrder);
    }
}

void ESpectrum::InitializePoints(CNuc* compound, const Config& configure) {
    for (auto& pt : points_) {
        pt.CalcEDependentValues(compound, configure);
        pt.CalcLegendreP(Config::maxLOrder, compound, nullptr);
        pt.CalcCoulombAmplitude(compound);
        if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE)
            pt.CalculateECAmplitudes(compound, configure);
    }
}

void ESpectrum::Calculate(CNuc* compound, const Config& configure) {
    for (auto& pt : points_)
        pt.Calculate(compound, configure);
}

double ESpectrum::Interpolate(double cmEnergy) const {
    if (points_.empty()) return 0.0;

    // Binary search: first point whose energy >= cmEnergy.
    auto it = std::lower_bound(points_.begin(), points_.end(), cmEnergy,
                               [](const EPoint& p, double e) {
                                   return p.GetCMEnergy() < e;
                               });

    if (it == points_.end())   return points_.back().GetFitCrossSection();
    if (it == points_.begin()) return points_.front().GetFitCrossSection();

    const EPoint& hi = *it;
    const EPoint& lo = *std::prev(it);

    const double eLo   = lo.GetCMEnergy();
    const double eHi   = hi.GetCMEnergy();
    const double span  = eHi - eLo;
    if (span < 1.0e-15) return lo.GetFitCrossSection();

    const double t = (cmEnergy - eLo) / span;
    return lo.GetFitCrossSection() + t * (hi.GetFitCrossSection() - lo.GetFitCrossSection());
}

EPoint* ESpectrum::GetPoint(int i) {
    return &points_[i - 1];
}

const EPoint* ESpectrum::GetPoint(int i) const {
    return &points_[i - 1];
}
