// The kinds of observable a site's record holds, their record layouts and how a reader contracts them
// (SiteObservable.h).
#include "SiteObservable.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "PPair.h"
#include "Site.h"
#include "TermSettings.h"

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// kLegendre: sum_L B_L P_L(cos theta) with the reader's row; for an angular-coefficient reader its order's
/// coefficient, a_k = B_k / (f A_0) or A_k = B_k / f (f the reader's EntranceSymmetry).
double ContractLegendre(const double *rec, const Site &reader, int slot) {
  const double *r = rec + azure2::Results::AngularOffset();
  if (std::isnan(r[0])) return kNaN;  // poisoned, never written
  const int n = (int)r[0];
  const double *b = r + 1;
  if (reader.IsCoefficient()) {
    const int k = reader.GetCoefficientOrder();
    const double ak = ((k >= 0 && k < n && k <= reader.GetMaxAngDistOrder()) ? b[k] : 0.0) / reader.EntranceSymmetry();
    return reader.GetCoefficientAbsolute() ? ak : ak / rec[azure2::Results::kA0];
  }
  const double *p = reader.Row(slot);
  double sum = 0.0;
  for (int l = 0; l < n; l++)
    if (b[l] != 0.0) sum += b[l] * p[l];
  return sum;
}

/// kElastic: the Coulomb term with the reader's amplitudes at the site's energy, the nuclear B_L and the
/// Coulomb-nuclear interference with the reader's row (SiteObservable.h).
double ContractElastic(const double *rec, const Site &reader, int slot) {
  using namespace elastic_record;
  const double *r = rec + azure2::Results::AngularOffset();
  if (std::isnan(r[0]) || std::isnan(r[kNumL])) return kNaN;
  complex fd, fx;
  reader.CoulombAt(slot, fd, fx);
  const int n = (int)r[0], nl = (int)r[kNumL];
  // the exchange terms only for identical particles with spin (f_x = 0, c_x = 0 and X_l = 0 otherwise)
  const bool exchange = r[kCross] != 0.0;
  double value = r[kCoulomb] * std::norm(fd);
  if (exchange) value += r[kCoulomb] * std::norm(fx) + r[kCross] * std::real(fd * std::conj(fx));
  const double *p = reader.Row(slot);
  for (int l = 0; l < n; l++)
    if (r[1 + l] != 0.0) value += r[1 + l] * p[l];
  // Re(f D) = Re f Re D - Im f Im D
  double it = 0.0;
  for (int l = 0; l < nl; l++) {
    double t = fd.real() * r[kDirect + 2 * l] - fd.imag() * r[kDirect + 2 * l + 1];
    if (exchange) t += fx.real() * r[kExchange + 2 * l] - fx.imag() * r[kExchange + 2 * l + 1];
    it += t * p[l];
  }
  return value + it;
}

}  // namespace

const SiteObservableRow kSiteObservables[3] = {
    {SiteObservable::kValue, "value at the site's angle", 0, nullptr},
    {SiteObservable::kLegendre, "Legendre coefficients B_L", Config::maxLOrder + 1, ContractLegendre},
    {SiteObservable::kElastic, "elastic amplitudes", elastic_record::kSlots, ContractElastic},
};

SiteObservable ObservableOf(const TermSettings &t) {
  if (t.phase || t.total || t.analyzingPower || t.thm) return SiteObservable::kValue;
  if (t.angularDist || (t.differential && t.entranceKey != t.exitKey)) return SiteObservable::kLegendre;
  if (t.differential) return SiteObservable::kElastic;
  return SiteObservable::kValue;
}

double LegendreSymmetryFactor(const PPair &entrance, bool elastic) {
  if (!entrance.IsIdentical()) return 1.0;
  return elastic ? 4.0 : entrance.GetEntranceSymmetryFactor();
}

SiteObservable Site::Observable() const { return ObservableOf(*settings_); }

int Site::ResultAngularCapacity() const { return kSiteObservables[(int)Observable()].angularSlots; }

double Site::ValueFor(const Site &reader, int slot, const azure2::Results &results) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  return row.contract ? row.contract(Rec(results), reader, slot) : GetFitCrossSection(results);
}

double Site::ContractRecord(SiteObservable kind, const double *record, const Site &reader, int slot) {
  return kSiteObservables[(int)kind].contract(record, reader, slot);
}

void Site::CoulombAt(int slot, complex &direct, complex &exchange) const {
  if (slot > 0) {
    if (!node_coulomb_)
      throw std::logic_error("Site::CoulombAt: a reader of elastic slots without their Coulomb amplitudes (segment " +
                             std::to_string(segment_key_) + "): Unit::CalcNodeRows was not called");
    direct = (*node_coulomb_)[2 * (std::size_t)(slot - 1)];
    exchange = (*node_coulomb_)[2 * (std::size_t)(slot - 1) + 1];
    return;
  }
  direct = coulombamplitude_;
  exchange = spinDependentCoulomb_ ? coulombexchange_ : complex(0.0, 0.0);
}

void Site::SetLegendreCoefficients(azure2::Results &results, const vector_r &coefficients) const {
  if ((int)coefficients.size() > std::min<int>(num_ang_slots_, Config::maxLOrder + 1))
    throw std::logic_error("Site::SetLegendreCoefficients: " + std::to_string(coefficients.size()) + " coefficients, " +
                           std::to_string(num_ang_slots_) + " result slots (segment " + std::to_string(segment_key_) +
                           ")");
  WriteLegendre(Rec(results) + azure2::Results::AngularOffset(), coefficients);
}

void Site::WriteLegendre(double *r, const vector_r &coefficients) {
  r[0] = (double)coefficients.size();
  std::copy(coefficients.begin(), coefficients.end(), r + 1);
}

void Site::SetElasticCoefficients(azure2::Results &results, const ElasticCoefficients &e) const {
  using namespace elastic_record;
  if (num_ang_slots_ != kSlots || (int)e.nuclear.size() > kOrders || (int)e.direct.size() > kOrders ||
      e.exchange.size() > e.direct.size())
    throw std::logic_error("Site::SetElasticCoefficients: " + std::to_string(e.direct.size()) + " orders, " +
                           std::to_string(num_ang_slots_) + " result slots (segment " + std::to_string(segment_key_) +
                           ")");
  WriteElastic(Rec(results) + azure2::Results::AngularOffset(), e);
}

void Site::WriteElastic(double *r, const ElasticCoefficients &e) {
  using namespace elastic_record;
  WriteLegendre(r, e.nuclear);
  r[kCoulomb] = e.coulomb;
  r[kCross] = e.cross;
  r[kNumL] = (double)e.direct.size();
  for (std::size_t l = 0; l < e.direct.size(); l++) {
    r[kDirect + 2 * l] = e.direct[l].real(), r[kDirect + 2 * l + 1] = e.direct[l].imag();
    const complex x = l < e.exchange.size() ? e.exchange[l] : complex(0.0, 0.0);
    r[kExchange + 2 * l] = x.real(), r[kExchange + 2 * l + 1] = x.imag();
  }
}

void Site::ScaleObservable(azure2::Results &results, double f) const {
  if (!num_ang_slots_) return;
  double *r = Rec(results) + azure2::Results::AngularOffset();
  const int n = (int)r[0];
  for (int l = 1; l <= n; l++) r[l] *= f;
  if (Observable() != SiteObservable::kElastic) return;
  using namespace elastic_record;
  r[kCoulomb] *= f;
  r[kCross] *= f;
  for (int l = 0; l < (int)r[kNumL]; l++)
    for (int i = 0; i < 2; i++) r[kDirect + 2 * l + i] *= f, r[kExchange + 2 * l + i] *= f;
}

void Site::PoisonObservable(azure2::Results &results) const {
  if (!num_ang_slots_) return;
  double *r = Rec(results) + azure2::Results::AngularOffset();
  std::fill(r + 1, r + 1 + num_ang_slots_, kNaN);
  r[0] = (double)std::min<int>(num_ang_slots_, Config::maxLOrder + 1);
  if (Observable() == SiteObservable::kElastic) r[elastic_record::kNumL] = elastic_record::kOrders;
}

double Site::CoefficientOfOrder(const Site &reader, const azure2::Results &results, int k) const {
  // As far as the reader's distribution is expanded (max_ang_dist_order); B_k / f is a_k A_0.
  const int n = NumLegendreCoefficients(results);
  const double bk = (k >= 0 && k < n && k <= reader.GetMaxAngDistOrder()) ? LegendreCoefficients(results)[k] : 0.0;
  return bk / reader.EntranceSymmetry();
}
