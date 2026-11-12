// The kinds of observable a site's record holds, their record layouts and how a reader contracts them
// (SiteObservable.h).
#include "SiteObservable.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "AngCoeff.h"
#include "CNuc.h"
#include "PathwayTable.h"
#include "PPair.h"
#include "PolarizationFunc.h"
#include "Site.h"
#include "TermSettings.h"

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// kLegendre: sum_L B_L P_L(cos theta) with the reader's row.
double ContractLegendre(const double *r, const AngularRow &row, double *) {
  if (std::isnan(r[0])) return kNaN;  // poisoned, never written
  const int n = (int)r[0];
  const double *b = r + 1;
  const double *p = row.p;
  double sum = 0.0;
  for (int l = 0; l < n; l++)
    if (b[l] != 0.0) sum += b[l] * p[l];
  return sum;
}

/// kElastic: the Coulomb term with the reader's amplitudes at the site's energy, the nuclear B_L and the
/// Coulomb-nuclear interference with the reader's row (SiteObservable.h).
double ContractElastic(const double *r, const AngularRow &row, double *) {
  using namespace elastic_record;
  if (std::isnan(r[kNumB]) || std::isnan(r[kNumL])) return kNaN;
  const complex fd = row.fd, fx = row.fx;
  const int n = (int)r[kNumB], nl = (int)r[kNumL];
  const double *b = r + kFirst, *d = b + n, *x = d + 2 * nl;
  // the exchange terms only for identical particles with spin (f_x = 0, c_x = 0 and X_l = 0 otherwise)
  const bool exchange = r[kCross] != 0.0;
  double value = r[kCoulomb] * std::norm(fd);
  if (exchange) value += r[kCoulomb] * std::norm(fx) + r[kCross] * std::real(fd * std::conj(fx));
  const double *p = row.p;
  for (int l = 0; l < n; l++)
    if (b[l] != 0.0) value += b[l] * p[l];
  // Re(f D) = Re f Re D - Im f Im D
  double it = 0.0;
  for (int l = 0; l < nl; l++) {
    double t = fd.real() * d[2 * l] - fd.imag() * d[2 * l + 1];
    if (exchange) t += fx.real() * x[2 * l] - fx.imag() * x[2 * l + 1];
    it += t * p[l];
  }
  return value + it;
}

/// A_y of particle coefficients (ay_record, kParticle): the amplitude matrix from C and the row's Y_l'^mu, the
/// Coulomb amplitude on its diagonal for elastic scattering, then Polarization::AmplitudeMatrix::AnalyzingPowerAy.
double ParticleAnalyzingPower(const double *r, const AngularRow &row) {
  using namespace ay_record;
  Polarization::AmplitudeMatrix m(0.5 * r[kSpins], 0.5 * r[kSpins + 1], 0.5 * r[kSpins + 2], 0.5 * r[kSpins + 3]);
  const int na = (int)r[kNumAmplitudes], nl = (int)r[kNumLp];
  if ((std::size_t)na != m.size()) return kNaN;
  const int muMax = Polarization::MaxProjection(0.5 * r[kSpins], 0.5 * r[kSpins + 1], 0.5 * r[kSpins + 2],
                                                0.5 * r[kSpins + 3]);
  vector_c own;
  const complex *y = row.y;
  int yMu = row.muMax;
  if (!y || yMu < muMax) {  // a reader without its table (its own site): Y at its angle
    Polarization::HarmonicTable(row.theta, nl > 0 ? nl - 1 : 0, muMax, own);
    y = own.data(), yMu = muMax;
  }
  const double *c = r + kCoefficients;
  for (int a = 0; a < na; a++) {
    const Polarization::Amplitude &amp = m.AmplitudeAt(a);
    const int mu = (int)std::lround(amp.v - amp.vp);
    complex value(0.0, 0.0);
    for (int lp = std::abs(mu); lp < nl; lp++) {
      const double *ca = c + 2 * ((std::size_t)a * nl + lp);
      value += complex(ca[0], ca[1]) * Polarization::AmplitudeMatrix::HarmonicOf(y, yMu, lp, mu);
    }
    m.SetAmplitude(a, value);
  }
  if (r[kElastic] != 0.0) m.AddCoulomb(row.fd, row.fx);
  return m.AnalyzingPowerAy();
}

/// A_y of capture coefficients (ay_record, kCapture): sum_k n_k P_k^1 / sum_k d_k P_k over the orders of the
/// reader's row, 0 where the denominator vanishes.
double CaptureAnalyzingPower(const double *r, const AngularRow &row) {
  using namespace ay_record;
  const int nk = std::min((int)r[kNumOrders], row.length);
  const double *d = r + kOrders, *n = d + (int)r[kNumOrders];
  vector_r own;
  const double *p1 = row.p1;
  if (!p1) {
    own.assign((std::size_t)std::max(nk, 1), 0.0);
    for (int k = 1; k < nk; k++) own[k] = AngCoeff::LegendreP1(k, std::cos(row.theta));
    p1 = own.data();
  }
  double num = 0.0, den = 0.0;
  for (int k = 0; k < nk; k++) {
    if (d[k] != 0.0) den += d[k] * row.p[k];
    if (n[k] != 0.0) num += n[k] * p1[k];
  }
  if (std::fabs(den) < 1.e-300) return 0.0;
  return num / den;
}

/// kAnalyzingPower: the cross section from its record (kElastic's or kLegendre's), A_y into `second`.
double ContractAnalyzingPower(const double *r, const AngularRow &row, double *second) {
  using namespace ay_record;
  if (std::isnan(r[kWeight]) || std::isnan(r[kForm])) {
    if (second) *second = kNaN;
    return kNaN;
  }
  const double *w = r + (int)r[kWeight];
  const double sigma = r[kElastic] != 0.0 ? ContractElastic(w, row, nullptr) : ContractLegendre(w, row, nullptr);
  if (second)
    *second = (int)r[kForm] == kCapture ? CaptureAnalyzingPower(r, row) : ParticleAnalyzingPower(r, row);
  return sigma;
}

int LegendreSlots(const RecordCapacity &c) { return c.orders; }
int ElasticSlots(const RecordCapacity &c) { return elastic_record::Slots(c.orders); }
/// The larger of the two forms, each with the cross section's record (its count included) after it.
int AnalyzingPowerSlots(const RecordCapacity &c) {
  using namespace ay_record;
  const int weight = 1 + std::max(c.orders, elastic_record::Slots(c.orders));
  const int particle = kCoefficients + 2 * c.amplitudes * c.orders;
  const int capture = kOrders + 2 * c.orders;
  return std::max(particle, capture) + weight - 1;
}

void ScaleLegendre(double *r, double f) {
  const int n = (int)r[0];
  for (int l = 1; l <= n; l++) r[l] *= f;
}

void ScaleElastic(double *r, double f) {
  using namespace elastic_record;
  const int end = kFirst + (int)r[kNumB] + 4 * (int)r[kNumL];
  r[kCoulomb] *= f;
  r[kCross] *= f;
  for (int i = kFirst; i < end; i++) r[i] *= f;
}

/// The cross section that weights A_y (an A_y sub-point's: GenMatrixFunc::CalculateCrossSection); A_y is a ratio.
void ScaleAnalyzingPower(double *r, double f) {
  using namespace ay_record;
  double *w = r + (int)r[kWeight];
  if (r[kElastic] != 0.0)
    ScaleElastic(w, f);
  else
    ScaleLegendre(w, f);
}

/// NaN coefficients under a count that reads them all.
void PoisonLegendre(double *r, int slots) {
  std::fill(r + 1, r + 1 + slots, kNaN);
  r[0] = (double)std::min<int>(slots, Config::maxLOrder + 1);
}

/// NaN counts: every reader is NaN (ContractElastic, ContractAnalyzingPower).
void PoisonAll(double *r, int slots) { std::fill(r, r + 1 + slots, kNaN); }

/// A reader of one Legendre coefficient (kLegendre): a_k = B_k / (f A_0) or A_k = B_k / f (f the reader's
/// EntranceSymmetry), as far as its distribution is expanded.
double CoefficientValue(const double *rec, const Site &reader) {
  const double *r = rec + azure2::Results::AngularOffset();
  if (std::isnan(r[0])) return kNaN;
  const int n = (int)r[0];
  const int k = reader.GetCoefficientOrder();
  const double ak = ((k >= 0 && k < n && k <= reader.GetMaxAngDistOrder()) ? r[1 + k] : 0.0) / reader.EntranceSymmetry();
  return reader.GetCoefficientAbsolute() ? ak : ak / rec[azure2::Results::kA0];
}

/// The reader's angular row at its slot `slot`: from its row set, or its own (slot 0, or no set).
AngularRow RowOf(const Site &reader, int slot) {
  AngularRow row;
  const azure2::RowSet *set = slot > 0 ? reader.RowSetOf() : nullptr;
  if (set) return set->RowAt(reader.RowSetIndex(slot));
  row.p = reader.Row(slot);
  row.length = reader.GetMaxLOrder() + 1;
  reader.CoulombAt(slot, row.fd, row.fx);
  row.theta = reader.GetCMAngle() * pi / 180.0;
  return row;
}

const double *Angular(const double *rec) { return rec + azure2::Results::AngularOffset(); }

}  // namespace

const SiteObservableRow kSiteObservables[4] = {
    {SiteObservable::kValue, "value at the site's angle", nullptr, nullptr, 1, nullptr, nullptr},
    {SiteObservable::kLegendre, "Legendre coefficients B_L", LegendreSlots, ContractLegendre, 1, ScaleLegendre,
     PoisonLegendre},
    {SiteObservable::kElastic, "elastic amplitudes", ElasticSlots, ContractElastic, 1, ScaleElastic, PoisonAll},
    {SiteObservable::kAnalyzingPower, "analysing power and its cross section", AnalyzingPowerSlots,
     ContractAnalyzingPower, 2, ScaleAnalyzingPower, PoisonAll},
};

RecordCapacity RecordCapacityOf(const CNuc &compound) {
  RecordCapacity c;
  c.orders = 1;
  const azure2::PathwayTable &pathways = compound.Pathways();
  for (int pair = 1; pair <= pathways.NumPairs(); pair++)
    for (const azure2::ExitRow &exit : pathways.Entrance(pair).exits) c.orders = std::max(c.orders, exit.num_orders);
  c.orders = std::min(c.orders, Config::maxLOrder + 1);
  int states = 0;
  for (int pair = 1; pair <= compound.NumPairs(); pair++) {
    const PPair *p = compound.GetPair(pair);
    if (p->GetPType() == 0)
      states = std::max(states, (int)std::lround((2.0 * p->GetJ(1) + 1.0) * (2.0 * p->GetJ(2) + 1.0)));
  }
  c.amplitudes = states * states;
  return c;
}

SiteObservable ObservableOf(const TermSettings &t) {
  if (t.phase || t.total || t.thm) return SiteObservable::kValue;
  if (t.analyzingPower) return SiteObservable::kAnalyzingPower;
  if (t.angularDist || (t.differential && t.entranceKey != t.exitKey)) return SiteObservable::kLegendre;
  if (t.differential) return SiteObservable::kElastic;
  return SiteObservable::kValue;
}

double LegendreSymmetryFactor(const PPair &entrance, bool elastic) {
  if (!entrance.IsIdentical()) return 1.0;
  return elastic ? 4.0 : entrance.GetEntranceSymmetryFactor();
}

SiteObservable Site::Observable() const { return ObservableOf(*settings_); }

int Site::ResultAngularCapacity(const RecordCapacity &capacity) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  return row.angularSlots ? row.angularSlots(capacity) : 0;
}

double Site::ValueFor(const Site &reader, int slot, const azure2::Results &results) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  if (!row.contract) return GetFitCrossSection(results);
  if (reader.IsCoefficient()) return CoefficientValue(Rec(results), reader);
  return row.contract(Angular(Rec(results)), RowOf(reader, slot), nullptr);
}

double Site::AnalyzingPowerFor(const Site &reader, int slot, const azure2::Results &results) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  if (row.values < 2) return GetAnalyzingPower(results);
  double ay = 0.0;
  row.contract(Angular(Rec(results)), RowOf(reader, slot), &ay);
  return ay;
}

void Site::ValuesWith(const azure2::RowSet &set, int k, const azure2::Results &results, double *values) const {
  const SiteObservableRow &row = kSiteObservables[(int)set.kind];
  values[0] = row.contract(Angular(Rec(results)), set.RowAt(k), row.values > 1 ? values + 1 : nullptr) /
              set.JacobianAt(k);
}

void Site::CoulombAt(int slot, complex &direct, complex &exchange) const {
  if (slot > 0 && rowset_ && !rowset_->coulomb.empty()) {
    rowset_->CoulombAt(rowset_at_ + slot - 1, direct, exchange);
    return;
  }
  if (slot > 0 && Observable() == SiteObservable::kElastic) {
      throw std::logic_error("Site::CoulombAt: a reader of elastic slots without their Coulomb amplitudes (segment " +
                             std::to_string(segment_key_) + "): it has no row set (WorkingCopy::BuildRowSets)");
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
  if (e.exchange.size() > e.direct.size() ||
      elastic_record::kFirst - 1 + (int)e.nuclear.size() + 4 * (int)e.direct.size() > num_ang_slots_)
    throw std::logic_error("Site::SetElasticCoefficients: " + std::to_string(e.nuclear.size()) + " B_L and " +
                           std::to_string(e.direct.size()) + " D_l, " + std::to_string(num_ang_slots_) +
                           " result slots (segment " + std::to_string(segment_key_) + ")");
  WriteElastic(Rec(results) + azure2::Results::AngularOffset(), e);
}

void Site::WriteElastic(double *r, const ElasticCoefficients &e) {
  using namespace elastic_record;
  const std::size_t n = e.nuclear.size(), m = e.direct.size();
  r[kNumB] = (double)n;
  r[kNumL] = (double)m;
  r[kCoulomb] = e.coulomb;
  r[kCross] = e.cross;
  double *b = r + kFirst, *d = b + n, *x = d + 2 * m;
  std::copy(e.nuclear.begin(), e.nuclear.end(), b);
  for (std::size_t l = 0; l < m; l++) {
    d[2 * l] = e.direct[l].real(), d[2 * l + 1] = e.direct[l].imag();
    const complex xl = l < e.exchange.size() ? e.exchange[l] : complex(0.0, 0.0);
    x[2 * l] = xl.real(), x[2 * l + 1] = xl.imag();
  }
}

void Site::SetAnalyzingPowerCoefficients(azure2::Results &results, const AnalyzingPowerCoefficients &a) const {
  using namespace ay_record;
  double *r = Rec(results) + azure2::Results::AngularOffset();
  const int body = a.capture ? kOrders + 2 * (int)a.denominator.size()
                             : kCoefficients + 2 * (int)a.coefficients.size();
  const int weight = a.elastic ? elastic_record::kFirst + (int)a.weightElastic.nuclear.size() +
                                     4 * (int)a.weightElastic.direct.size()
                               : 1 + (int)a.weightLegendre.size();
  if (body + weight > 1 + num_ang_slots_ || a.weightElastic.exchange.size() > a.weightElastic.direct.size())
    throw std::logic_error("Site::SetAnalyzingPowerCoefficients: " + std::to_string(body + weight) + " slots, " +
                           std::to_string(1 + num_ang_slots_) + " in the record (segment " +
                           std::to_string(segment_key_) + ")");
  r[kWeight] = body;
  r[kForm] = a.capture ? kCapture : kParticle;
  r[kElastic] = a.elastic ? 1.0 : 0.0;
  if (a.capture) {
    const std::size_t k = a.denominator.size();
    r[kNumOrders] = (double)k;
    std::copy(a.denominator.begin(), a.denominator.end(), r + kOrders);
    std::copy(a.numerator.begin(), a.numerator.end(), r + kOrders + k);
  } else {
    for (int i = 0; i < 4; i++) r[kSpins + i] = 2.0 * a.spins[i];
    r[kNumAmplitudes] = a.amplitudes;
    r[kNumLp] = a.orders;
    for (std::size_t i = 0; i < a.coefficients.size(); i++)
      r[kCoefficients + 2 * i] = a.coefficients[i].real(), r[kCoefficients + 2 * i + 1] = a.coefficients[i].imag();
  }
  if (a.elastic)
    WriteElastic(r + body, a.weightElastic);
  else
    WriteLegendre(r + body, a.weightLegendre);
}

void Site::SetIsotropicObservable(azure2::Results &results, double value) const {
  if (Observable() == SiteObservable::kElastic) {
    ElasticCoefficients e;
    e.nuclear = {value};
    SetElasticCoefficients(results, e);
  } else {
    SetLegendreCoefficients(results, {value});
  }
}

void Site::ScaleObservable(azure2::Results &results, double f) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  if (num_ang_slots_ && row.scale) row.scale(Rec(results) + azure2::Results::AngularOffset(), f);
}

void Site::PoisonObservable(azure2::Results &results) const {
  const SiteObservableRow &row = kSiteObservables[(int)Observable()];
  if (num_ang_slots_ && row.poison) row.poison(Rec(results) + azure2::Results::AngularOffset(), num_ang_slots_);
}

double Site::CoefficientOfOrder(const Site &reader, const azure2::Results &results, int k) const {
  // As far as the reader's distribution is expanded (max_ang_dist_order); B_k / f is a_k A_0.
  const int n = NumLegendreCoefficients(results);
  const double bk = (k >= 0 && k < n && k <= reader.GetMaxAngDistOrder()) ? LegendreCoefficients(results)[k] : 0.0;
  return bk / reader.EntranceSymmetry();
}
