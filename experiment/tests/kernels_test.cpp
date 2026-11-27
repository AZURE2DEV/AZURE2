// The table of apparatus kernels (azure2/experiment/Kernels.h) on nodes laid out by hand: every row is complete, the
// weights the gradient reads (the kernel on LinearForm) reproduce the forward value, and each kernel against a closed
// form (a constant survives every normalised kernel; a thin target of constant stopping averages a linear cross
// section to its value at the mid-depth; a thick target integrates a constant to 1e-24 sigma dE/epsilon).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void Check(const std::string &what, bool ok) {
  std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

void Close(const std::string &what, double got, double want, double tolerance) {
  const double scale = std::max(1.0, std::fabs(want));
  const bool ok = std::fabs(got - want) <= tolerance * scale;
  std::printf("%s %s: %.17g (want %.17g)\n", ok ? "  ok  " : "  FAIL", what.c_str(), got, want);
  if (!ok) failures++;
}

/// Nodes from hi down to lo, n intervals, with a stopping per node.
struct Grid {
  std::vector<double> energy, stopping;
  KernelNodes nodes;
  Grid(double hi, double lo, int n, double eps) {
    for (int i = 0; i <= n; i++) energy.push_back(hi - (hi - lo) * i / n);
    stopping.assign(energy.size(), eps);
    nodes.n = static_cast<int>(energy.size());
    nodes.energy = energy.data();
    nodes.stopping = stopping.data();
  }
};

double ConstantSigma(const void *owner, const void *, double, double cmFactor) {
  return cmFactor * *static_cast<const double *>(owner);
}

/// The fold of `values` and the dot product of its LinearForm weights with them.
void FoldBoth(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const Grid &g,
              const std::vector<double> &values, double &value, double &fromWeights) {
  EffectKernels ek;
  ek.resolution = &k;
  value = Fold(ek, Broadening::kApply, p, point, g.nodes, values.data());
  std::vector<LinearForm> unit;
  for (int j = 1; j <= g.nodes.n; j++) unit.push_back(LinearForm::Unit(j));
  const LinearForm form = FoldForm(ek, Broadening::kApply, p, point, g.nodes, unit.data());
  fromWeights = 0.0;
  for (const auto &t : form.Terms()) fromWeights += t.second * values[t.first - 1];
}


/// The value of a form on `values` (node j at values[j - 1]); a node listed twice counts twice.
double Evaluate(const LinearForm &form, const std::vector<double> &values) {
  double sum = 0.0;
  for (const auto &t : form.Terms()) sum += t.second * values[t.first - 1];
  return sum;
}

/// LinearForm against plain arithmetic: the same expression built on doubles and on forms, in the shapes the
/// kernels use (interpolations of two nodes, scaled, summed over a window that moves up, down and apart).
void LinearFormArithmetic() {
  std::printf("the linear form\n");
  const int n = 400;
  std::vector<double> values(n);
  for (int j = 0; j < n; j++) values[j] = std::sin(0.37 * j) + 0.01 * j;
  auto node = [&](int j) { return values[j - 1]; };
  double want = 0.0;
  LinearForm got = 0.0;
  // windows of interpolated values: upward, downward (the low end grows), disjoint from the sum so far, overlapping
  const int starts[] = {200, 150, 320, 10, 190};
  for (int w = 0; w < 5; w++) {
    double inner = 0.0;
    LinearForm innerForm = 0.0;
    for (int i = 0; i < 40; i++) {
      const int j = (w == 1) ? starts[w] + 40 - i : starts[w] + i;
      const double f = 0.1 + 0.02 * i, weight = 1.0 / (1.0 + i);
      inner += (node(j) * (1.0 - f) + node(j + 1) * f) * weight;
      innerForm += (LinearForm::Unit(j) * (1.0 - f) + LinearForm::Unit(j + 1) * f) * weight;
    }
    inner /= 3.0, innerForm /= 3.0;
    want += inner * (0.5 + w), got += innerForm * (0.5 + w);
  }
  Close("sums of scaled interpolations over moving windows", Evaluate(got, values), want, 1e-13);
  LinearForm two = LinearForm::Unit(7) * 2.0 + LinearForm::Unit(9) * -3.0;
  Check("two terms stay two terms", two.Terms().size() == 2);
  Close("a small form added to a large one, and a large one to a small one",
        Evaluate(got + two, values) + Evaluate(two + got, values), 2.0 * (want + 2.0 * node(7) - 3.0 * node(9)), 1e-13);
  LinearForm same = LinearForm::Unit(5) + LinearForm::Unit(5) + LinearForm::Unit(5);
  Close("one node three times", Evaluate(same, values), 3.0 * node(5), 1e-15);
  Check("a sum of forms is valid; a nonzero constant is not", got.Valid() && !(got + LinearForm(1.0)).Valid());
  Check("zero is the empty form", LinearForm(0.0).Terms().empty() && LinearForm(0.0).Valid());
}

}  // namespace

/// The window search of GridQuadrature against a walk over every interval: the bisection may only leave out
/// intervals the walk's own test skips (nodes descending, some repeated; windows inside, outside, empty, NaN).
void WindowSearch() {
  std::printf("the window of GridQuadrature by bisection\n");
  unsigned state = 12345u;
  auto next = [&state]() {
    state = state * 1664525u + 1013904223u;
    return (state >> 8) / 16777216.0;
  };
  const double nan = std::numeric_limits<double>::quiet_NaN();
  bool same = true, bracket = true;
  int windows = 0;
  for (int n : {0, 1, 2, 3, 7, 64, 501}) {
    std::vector<double> energy;
    double e = 2.0;
    for (int i = 0; i < n; i++) {
      energy.push_back(e);
      if (next() > 0.15) e -= 0.001 + 0.01 * next();  // else: a repeated node
    }
    KernelNodes nodes;
    nodes.n = n;
    nodes.energy = energy.data();
    for (int trial = 0; trial < 400; trial++) {
      double lo = e - 0.05 + (2.05 - e + 0.1) * next(), hi = lo + (trial % 3 == 0 ? 0.004 : 1.0) * next();
      if (trial % 50 == 1) hi = lo;
      if (trial % 50 == 2) std::swap(lo, hi);
      if (trial % 50 == 3) lo = nan;
      if (trial % 50 == 4) hi = nan;
      if (trial % 50 == 5 && n > 2) lo = energy[n / 2], hi = energy[n / 3];  // on nodes
      // the interval that brackets an energy: the first the walk finds
      const double at = trial % 50 == 6 && n > 0 ? energy[(std::size_t)(next() * n) % n] : trial % 50 == 7 ? nan : lo;
      int walkedTo = 0;
      for (int j = 1; j < n && !walkedTo; j++)
        if (at <= nodes.Energy(j) && at >= nodes.Energy(j + 1)) walkedTo = j;
      if (BracketingInterval(nodes, at) != walkedTo) bracket = false;
      const int first = FirstIntervalBelow(nodes, hi), end = EndIntervalAbove(nodes, lo);
      for (int j = 1; j < n; j++) {
        const bool walked = !(nodes.Energy(j + 1) >= hi || nodes.Energy(j) <= lo);
        if (walked && !(j >= first && j < end)) same = false;
      }
      if (first < 1 || end > std::max(n, 1)) same = false;
      windows++;
    }
  }
  Check("every interval the walk integrates is inside the bisected range (" + std::to_string(windows) + " windows)", same);
  Check("the interval that brackets an energy is the first a walk finds (on nodes, repeated nodes, outside, NaN)", bracket);
  // and the quadrature itself: a window of a few intervals of many
  Grid g(2.0, 1.0, 4000, 1.0);
  std::vector<double> sigma;
  for (double x : g.energy) sigma.push_back(1.0 + 2.0 * x);  // linear: its interpolant is exact
  auto one = [](double) { return 1.0; };
  double num = 0.0, den = 0.0;
  GridQuadrature<double>(g.nodes, sigma.data(), 1.49990, 1.50035, false, one, one, num, den);
  Close("a narrow window on 4000 intervals: its width", den, 0.00045, 1e-13);
  Close("... and the integral of 1 + 2 E", num, 0.00045 + (1.50035 * 1.50035 - 1.49990 * 1.49990), 1e-13);
}

int main() {
  LinearFormArithmetic();
  std::printf("the table\n");
  std::set<std::string> names;
  for (const Kernel &k : Kernels()) {
    Check(k.name + ": a unique name", names.insert(k.name).second);
    Check(k.name + ": a format-3 keyword", !k.keyword.empty());
    Check(k.name + ": a support", k.support != nullptr);
    if (k.stage == Kernel::Stage::kResolution)
      Check(k.name + ": fold and weights", k.fold && k.weights && !k.broaden);
    else
      Check(k.name + ": the broadening on both value types, d/dT and its weights",
            k.broaden && k.broadenForm && k.broadenDT && k.broadeningWeights && !k.fold);
    Check(k.name + ": KernelNamed finds it", KernelNamed(k.name) == &k);
    if (k.tableForm) {
      const std::string field = k.keyword.substr(0, k.keyword.find(' '));
      const std::string form = k.keyword.substr(k.keyword.find(' ') + 1);
      Check(k.name + ": KernelOfForm finds it", KernelOfForm(field.substr(0, field.find('.')),
                                                             field.substr(field.find('.') + 1), form) == &k);
      for (const KernelParameter &p : k.parameters)
        Check(k.name + "." + p.name + ": a dimension and a unit", !p.dimension.empty() && !p.unit.empty());
    }
  }
  Check("no kernel named 'none'", KernelNamed("none") == nullptr);

  const double sigmaLab = 0.004;  // MeV
  KernelParameters p;
  p.sigma.at = ConstantSigma;
  p.sigma.owner = &sigmaLab;
  KernelPoint point;
  point.energy = 1.0;
  point.labEnergy = 1.2;
  point.nominalLabEnergy = 1.2;
  point.cmFactor = 1.0 / 1.2;

  std::printf("Gaussian beam spread\n");
  {
    const Kernel &k = *KernelNamed("gaussian");
    p.gaussian = true;
    Grid g(1.0 + 5.0 * sigmaLab, 1.0 - 5.0 * sigmaLab, 40, 1.0);
    double v, w;
    FoldBoth(k, p, point, g, std::vector<double>(g.energy.size(), 3.0), v, w);
    Close("a constant survives", v, 3.0, 1e-14);
    Close("the weights give the same value", w, v, 1e-14);
    std::vector<double> linear;
    for (double e : g.energy) linear.push_back(2.0 + 5.0 * (e - 1.0));
    FoldBoth(k, p, point, g, linear, v, w);
    Close("a linear cross section about the centre is its value there", v, 2.0, 1e-13);
    Close("the weights give the same value", w, v, 1e-14);
  }

  std::printf("thin target, constant stopping\n");
  {
    const Kernel &k = *KernelNamed("target");
    p.gaussian = false;
    const double eps = 1e-20;      // c.m. MeV cm^2 (10 eV per 1e15 atoms/cm^2)
    const double density = 1e18;   // atoms/cm^2
    point.density = density;
    point.targetThickness = eps * density;  // 0.01 MeV
    Grid g(1.0, 1.0 - point.targetThickness, 20, eps);
    std::vector<double> linear;
    for (double e : g.energy) linear.push_back(2.0 + 5.0 * (e - 1.0));
    double v, w;
    FoldBoth(k, p, point, g, linear, v, w);
    Close("the depth average of a linear cross section is its value at mid-depth", v, 2.0 - 5.0 * 0.005, 1e-12);
    Close("the weights give the same value", w, v, 1e-14);
    p.thickTarget = true;
    FoldBoth(k, p, point, g, std::vector<double>(g.energy.size(), 3.0), v, w);
    Close("thick: 1e-24 sigma dE / epsilon", v, 1e-24 * 3.0 * point.targetThickness / eps, 1e-12);
    p.thickTarget = false;
  }

  std::printf("rectangular resolution\n");
  {
    // window 20 keV lab about E0 = 1 MeV c.m.; f = 1/1.2, so 1/60 MeV c.m.; nodes every 1/6000 MeV over 3 windows
    const Kernel &k = *KernelNamed("rectangular");
    double width = 0.02;  // MeV, lab
    p.values = &width;
    p.numValues = 1;
    const double f = point.cmFactor, w = f * width;
    Grid g(1.0 + 1.5 * w, 1.0 - 1.5 * w, 300, 1.0);
    auto at = [&](double (*s)(double)) {
      std::vector<double> v;
      for (double e : g.energy) v.push_back(s(e));
      return v;
    };
    double v, wv;
    FoldBoth(k, p, point, g, at([](double e) { return 2.0 + 5.0 * (e - 1.0); }), v, wv);
    Close("sigma linear in E: the average is sigma(E0)", v, 2.0, 1e-13);
    Close("the weights give the same value", wv, v, 1e-14);
    // sigma'' = 2 a: sigma(E0) + sigma'' w^2/24, to the linear interpolation between nodes (h^2 sigma''/12 per
    // interval averaged: h = w/100, a relative 1e-4 of the curvature term)
    const std::vector<double> quad = at([](double e) { return 2.0 + 3.0 * (e - 1.0) + 400.0 * (e - 1.0) * (e - 1.0); });
    FoldBoth(k, p, point, g, quad, v, wv);
    const double h = (g.energy[0] - g.energy[1]);
    Close("sigma quadratic: sigma(E0) + sigma'' w^2/24 (+ sigma'' h^2/12 of the interpolation)", v,
          2.0 + 800.0 * w * w / 24.0 + 800.0 * h * h / 12.0, 1e-12);
    // a step from 0 to 1 between two nodes just below E0 + w/4: the node values ramp across one interval
    Grid s(1.0 + 1.5 * w, 1.0 - 1.5 * w, 300, 1.0);
    const double eStep = s.energy[125];  // a node 0.25 w above E0 (125 intervals of w/100 below the top, 1.5 w)
    std::vector<double> step;
    for (double e : s.energy) step.push_back(e >= eStep ? 1.0 : 0.0);
    FoldBoth(k, p, point, s, step, v, wv);
    // the window [E0 - w/2, E0 + w/2] holds (E0 + w/2 - eStep) of value 1 and half an interval h of the ramp
    Close("a step: the fraction of the window above it", v, (1.0 + 0.5 * w - eStep + 0.5 * h) / w, 1e-12);
    // d/dw against central differences, on the quadratic (window edges between nodes)
    width = 0.0201;
    const double d = k.derivative(p, point, g.nodes, quad.data(), 0);
    double plus, minus;
    width = 0.0201 * (1.0 + 1e-6);
    plus = k.fold(p, point, g.nodes, quad.data());
    width = 0.0201 * (1.0 - 1e-6);
    minus = k.fold(p, point, g.nodes, quad.data());
    Close("dY/dw against central differences", d, (plus - minus) / (2.0 * 0.0201e-6), 1e-6);
    SupportPoint sp;
    sp.energy = 1.0;
    sp.cmFactor = f;
    double hi, lo;
    EffectKernels ek;
    ek.resolution = &k;
    width = 0.02;
    EffectSupport(ek, p, sp, hi, lo);
    Check("support: the window [E0 - f w/2, E0 + f w/2]", hi == 1.0 + 0.5 * w && lo == 1.0 - 0.5 * w);
    // a fitted width: the grid spans the layout value, the fold the present one
    double layout = 0.03;
    p.layoutValues = &layout;
    EffectSupport(ek, p, sp, hi, lo);
    Check("support of a fitted width: the window at its layout value", hi == 1.0 + 0.5 * f * layout && lo == 1.0 - 0.5 * f * layout);
    p.layoutValues = nullptr;
    p.values = nullptr;
    p.numValues = 0;
  }

  std::printf("triangular resolution\n");
  {
    // base 20 keV lab about E0 = 1 MeV c.m.; nodes every w/100 c.m. over 3 bases, E0 a node
    const Kernel &k = *KernelNamed("triangular");
    double width = 0.02;  // MeV, lab
    p.values = &width;
    p.numValues = 1;
    const double f = point.cmFactor, w = f * width, h = 0.5 * w;
    Grid g(1.0 + 1.5 * w, 1.0 - 1.5 * w, 300, 1.0);
    auto at = [&](double (*s)(double)) {
      std::vector<double> v;
      for (double e : g.energy) v.push_back(s(e));
      return v;
    };
    double v, wv;
    FoldBoth(k, p, point, g, at([](double e) { return 2.0 + 5.0 * (e - 1.0); }), v, wv);
    Close("sigma linear in E: the average is sigma(E0)", v, 2.0, 1e-13);
    Close("the weights give the same value", wv, v, 1e-14);
    // sigma = 2 + 3 x + 400 x^2: sigma(E0) + 400 <x^2>, <x^2> = h^2/6 for the triangle, plus the linear interpolation's
    // h_n^2/12 per interval of the curvature (h_n the node spacing)
    const std::vector<double> quad = at([](double e) { return 2.0 + 3.0 * (e - 1.0) + 400.0 * (e - 1.0) * (e - 1.0); });
    FoldBoth(k, p, point, g, quad, v, wv);
    const double hn = g.energy[0] - g.energy[1];
    Close("sigma quadratic: sigma(E0) + sigma'' h^2/12 (+ sigma'' h_n^2/12 of the interpolation)", v,
          2.0 + 800.0 * h * h / 12.0 + 800.0 * hn * hn / 12.0, 1e-12);
    // d/dw against central differences, on the quadratic (edges between nodes)
    width = 0.0201;
    const double d = k.derivative(p, point, g.nodes, quad.data(), 0);
    width = 0.0201 * (1.0 + 1e-6);
    const double plus = k.fold(p, point, g.nodes, quad.data());
    width = 0.0201 * (1.0 - 1e-6);
    const double minus = k.fold(p, point, g.nodes, quad.data());
    Close("dY/dw against central differences", d, (plus - minus) / (2.0 * 0.0201e-6), 1e-6);
    // and against the closed form d/dw [400 h^2/6] = 400 h f/6 (the interpolation term does not depend on w)
    Close("dY/dw: the closed form 400 h f / 6", d, 400.0 * (0.5 * f * 0.0201) * f / 6.0, 1e-6);
    width = 0.02;
    SupportPoint sp;
    sp.energy = 1.0;
    sp.cmFactor = f;
    double hi, lo;
    EffectKernels ek;
    ek.resolution = &k;
    EffectSupport(ek, p, sp, hi, lo);
    Check("support: the base [E0 - f w/2, E0 + f w/2]", hi == 1.0 + h && lo == 1.0 - h);
    p.values = nullptr;
    p.numValues = 0;
  }

  std::printf("support\n");
  {
    SupportPoint s;
    s.energy = 1.0;
    s.cmFactor = 1.0;
    s.minimumEnergy = 0.001;
    double hi, lo;
    EffectKernels ek;
    ek.resolution = KernelNamed("gaussian");
    p.gaussian = true;
    EffectSupport(ek, p, s, hi, lo);
    Close("Gaussian: E0 + 5 sigma", hi, 1.0 + 5.0 * sigmaLab, 1e-15);
    Close("Gaussian: E0 - 5 sigma", lo, 1.0 - 5.0 * sigmaLab, 1e-15);
    ek.resolution = nullptr;
    EffectSupport(ek, p, s, hi, lo);
    Check("no kernel: the point itself", hi == 1.0 && lo == 1.0);
  }

  // The rule of include/ShiftedFold.h (the owner's answer of 2026-11-18): a point measured at E_n with the shift s
  // folds as the same point without a shift at E_n + s, its reconstructed-energy window moved by s.  With the
  // detailed-balance weight the norm is the weight at the point's own energy, the true one; at the nominal energy (the
  // kernel before) the shifted point came out 7.7e-4 off its twin here (D ~ E / E_gamma^2, s / E = 0.25 %).
  std::printf("beam profile: a shifted point is its twin at the shifted energy\n");
  {
    const Kernel &k = *KernelNamed("beam_profile");
    const std::vector<BeamProfileComponent> profile = {{1.21, 0.02, -1.5, 1.0}};  // lab MeV
    KernelParameters q;
    q.profile = &profile;
    q.tpcSigma = 0.003;
    q.detailedBalance = true;
    const double f = 1.0 / 1.2, s = 0.003;  // the c.m./lab factor, the shift (lab MeV)
    KernelPoint shifted;
    shifted.cmFactor = f;
    shifted.nominalLabEnergy = 1.2;
    shifted.labEnergy = 1.2 + s;
    shifted.energy = f * (1.2 + s);
    shifted.binLow = f * 1.196;  // the window as stored: unshifted, c.m.
    shifted.binHigh = f * 1.204;
    shifted.photoQ = 1.9;          // Q_gamma (MeV)
    shifted.photoMass = 13.0 * 931.494;
    KernelPoint twin = shifted;
    twin.nominalLabEnergy = twin.labEnergy;
    twin.binLow = f * (1.196 + s);
    twin.binHigh = f * (1.204 + s);
    Grid g(f * 1.32, f * 1.08, 480, 1.0);
    for (const bool constant : {true, false}) {
      std::vector<double> sigma;
      for (double e : g.energy) sigma.push_back(constant ? 3.0 : 1.0 + 40.0 * (e - 1.0) * (e - 1.0));
      double vs, ws, vt, wt;
      FoldBoth(k, q, shifted, g, sigma, vs, ws);
      FoldBoth(k, q, twin, g, sigma, vt, wt);
      const std::string what = constant ? "a constant" : "a parabola";
      Close(what + ": the shifted point = its twin (relative)", vs / vt, 1.0, 1e-12);
      Close(what + ": the weights give the same value", ws, vs, 1e-13);
    }
  }

  WindowSearch();

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
