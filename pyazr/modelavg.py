"""Model averaging over a set of fitted model variants -- the bookkeeping only.

A THM (or any R-matrix) result depends on choices the data cannot make: the
channel radius, the entrance vertex (plane or distorted waves), the optical
potentials, the spectator-momentum window.  Each choice is a *variant*: a
complete model, fitted on its own.  This module turns a list of fitted
variants into one number per quantity, with its statistical and its model
spread kept apart.  It does no fitting: a variant is a record the caller
fills in after its own fit (``scripts/thm_model_average.py`` is one driver).

For variant ``i`` with weight ``w_i`` (normalised over the variants that have
the quantity), value ``m_i`` and statistical variance ``s_i^2``::

    mean            m   = sum_i w_i m_i
    stat. variance  S^2 = sum_i w_i s_i^2          (weighted mean of the variances)
    model variance  M^2 = sum_i w_i (m_i - m)^2    (weighted variance of the means)
    total variance  T^2 = S^2 + M^2

Weights, ``w_i ∝ prior_i exp(-IC_i / 2)`` (normalised; only differences of
``IC`` matter):

``"aic"`` (default)  ``IC = chi2/s + 2k``              (Akaike)
``"bic"``            ``IC = chi2/s + k ln N``           (Schwarz)
``"chi2"``           ``IC = chi2/s``
``"flat"``           every variant 1 (times its prior)

``s`` is the optional rescale factor (``rescale=``): ``None`` or ``1`` leaves
chi2 alone; ``"best"`` divides every chi2 by ``max(1, chi2/nu)`` of the
variant with the lowest unscaled criterion; a number is used as given.  When
the data errors are underestimated, chi2 differences are inflated by the same
factor and Akaike weights are over-confident (one variant takes all the
weight); ``rescale="best"`` is the usual remedy.  The weights compare models
of the *same* data: variants with different ``N`` trigger a warning.

Parameters are matched across variants by name.  :func:`parameter_label`
gives an R-matrix parameter a name that the ``.azr`` alone determines
(``E[2-#1]`` for the energy of the first 2- level, ``G[2-#1;p1;L1;S1]`` for
its width in the channel of file pair key 1, orbital L = 1, channel spin
S = 1), so it is the same in every variant of one level scheme and
:func:`write_averaged_azr` can find it in a template ``.azr``.

Pure Python: no numpy at import or at run time, so this module and its test
run where numpy is missing.
"""

import csv
import importlib.util
import json
import math
import os
import warnings
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = ["Variant", "Averaged", "ModelAverage", "model_average",
           "parameter_label", "level_label", "width_label",
           "azr_parameter_slots", "write_averaged_azr", "WEIGHT_METHODS"]

WEIGHT_METHODS = ("aic", "bic", "chi2", "flat")


# ---------------------------------------------------------------------------
#  Names
# ---------------------------------------------------------------------------

def _jpi(J, parity):
    J = float(J)
    j = str(int(J)) if J.is_integer() else f"{int(round(2 * J))}/2"
    return f"{j}{'+' if int(parity) > 0 else '-'}"


def _num(v):
    """1.0 -> "1", 0.5 -> "0.5" (channel spins, as short text)."""
    v = float(v)
    return str(int(v)) if v.is_integer() else f"{v:g}"


def level_label(J, parity, n):
    """The name of level ``n`` (1-based, in file order among the *active*
    levels of that J^pi -- the engine's level number within its J group)."""
    return f"{_jpi(J, parity)}#{int(n)}"


def width_label(J, parity, n, pair_key, L, S):
    """The name of a width: level, file pair key, L and channel spin S."""
    return f"G[{level_label(J, parity, n)};p{int(pair_key)};L{int(L)};S{_num(S)}]"


def parameter_label(param, pairs=None):
    """The model-averaging name of a :class:`pyazr.parameters.Parameter`.

    ``E[<J^pi>#<n>]`` for a level energy, ``G[<J^pi>#<n>;p<key>;L<l>;S<s>]``
    for a width, ``norm[<segment key>]`` / ``shift[<segment key>]`` for a
    segment's normalization / energy shift, the engine's own name otherwise
    (``cbkg_*``).  ``<n>`` is the engine's level number within the J^pi group
    (``Parameter.level``), ``<key>`` the pair key the ``.azr`` writes: pass
    ``pairs=session.pairs`` so the engine's pair number is translated (the
    two differ when ``<levels>`` does not introduce the pairs in key order;
    see :meth:`pyazr.azrfile.AzrModel.apply_fit`).
    """
    kind = getattr(param, "kind", None)
    if kind == "energy":
        return f"E[{level_label(param.J, param.parity, param.level)}]"
    if kind == "width":
        key = param.pair
        if pairs is not None:
            key = {p.number: p.key for p in pairs}.get(param.pair, param.pair)
        return width_label(param.J, param.parity, param.level, key,
                           param.L, param.S or 0.0)
    if kind in ("norm", "shift"):
        return f"{kind}[{param.segment_key}]"
    return str(param.name)


# ---------------------------------------------------------------------------
#  Variants and the average
# ---------------------------------------------------------------------------

def _split_value(v):
    """``x`` or ``(x, sigma)`` -> (float, float or None)."""
    if isinstance(v, (tuple, list)):
        if len(v) != 2:
            raise ValueError(f"expected a value or (value, sigma), got {v!r}")
        return float(v[0]), (None if v[1] is None else float(v[1]))
    return float(v), None


class Variant:
    """One fitted model variant.

    ``label``      a short unique name
    ``settings``   dict of what defines the variant (radius, vertexModel, ...)
    ``chi2``       the minimised objective (data chi2 plus any penalty rows)
    ``npoints``    N, the number of data points (residual rows)
    ``nfree``      k, the number of free parameters
    ``values``     ``{name: best-fit value}``
    ``errors``     optional ``{name: 1-sigma statistical error}``
    ``covariance`` optional ``(names, matrix)``: a statistical covariance
                   over ``names`` (list of lists); its diagonal fills
                   ``errors`` for names that have none
    ``derived``    optional ``{name: value or (value, sigma)}``, quantities
                   computed from the fit (a strength, S(E) at an energy)
    ``prior``      prior weight (default 1), e.g. < 1 for an extreme radius
    """

    def __init__(self, label, settings=None, chi2=0.0, npoints=0, nfree=0,
                 values=None, errors=None, covariance=None, derived=None,
                 prior=1.0):
        self.label = str(label)
        self.settings = dict(settings or {})
        self.chi2 = float(chi2)
        self.npoints = int(npoints)
        self.nfree = int(nfree)
        self.values = {str(k): float(v) for k, v in (values or {}).items()}
        self.errors = {str(k): float(v) for k, v in (errors or {}).items()
                       if v is not None}
        self.covariance = None
        if covariance is not None:
            names, mat = covariance
            names = [str(n) for n in names]
            mat = [[float(c) for c in row] for row in mat]
            if len(mat) != len(names) or any(len(r) != len(names) for r in mat):
                raise ValueError(f"variant {label}: covariance is not "
                                 f"{len(names)} x {len(names)}")
            self.covariance = (names, mat)
            for i, n in enumerate(names):
                if n not in self.errors and mat[i][i] >= 0:
                    self.errors[n] = math.sqrt(mat[i][i])
        self.derived = {}
        self.derived_errors = {}
        for k, v in (derived or {}).items():
            val, sig = _split_value(v)
            self.derived[str(k)] = val
            if sig is not None:
                self.derived_errors[str(k)] = sig
        self.prior = float(prior)
        if not self.prior >= 0:
            raise ValueError(f"variant {label}: prior must be >= 0")
        if not math.isfinite(self.chi2):
            raise ValueError(f"variant {label}: chi2 is not finite")

    @property
    def dof(self):
        return self.npoints - self.nfree

    def to_dict(self):
        return dict(label=self.label, settings=self.settings, chi2=self.chi2,
                    npoints=self.npoints, nfree=self.nfree, values=self.values,
                    errors=self.errors, derived=self.derived,
                    derived_errors=self.derived_errors, prior=self.prior,
                    covariance=(None if self.covariance is None else
                                {"names": self.covariance[0],
                                 "matrix": self.covariance[1]}))

    @classmethod
    def from_dict(cls, d):
        cov = d.get("covariance")
        if cov is not None:
            cov = (cov["names"], cov["matrix"])
        derived = dict(d.get("derived") or {})
        for k, s in (d.get("derived_errors") or {}).items():
            if k in derived:
                derived[k] = (derived[k], s)
        return cls(d["label"], d.get("settings"), d["chi2"], d["npoints"],
                   d["nfree"], d.get("values"), d.get("errors"), cov, derived,
                   d.get("prior", 1.0))

    def __repr__(self):
        return (f"Variant({self.label!r}, chi2={self.chi2:.4g}, N={self.npoints}, "
                f"k={self.nfree}, {len(self.values)} values, "
                f"{len(self.derived)} derived)")


class Averaged:
    """The average of one quantity.

    ``mean``; ``stat``, ``model``, ``total`` are standard deviations (the
    square roots of S^2, M^2 and S^2 + M^2 of the module docstring).  ``stat``
    (and so ``total``) is ``None`` when a variant that carries the quantity
    has no statistical error for it.  ``coverage`` is the summed weight of the
    variants that carry the quantity (1 when all do; the average is over
    those, with their weights renormalised); ``values`` lists
    ``(label, weight, value, sigma)`` of the contributing variants.
    """

    __slots__ = ("name", "mean", "stat", "model", "total", "coverage",
                 "values")

    def __init__(self, name, mean, stat, model, coverage, values):
        self.name = name
        self.mean = mean
        self.stat = stat
        self.model = model
        self.total = (None if stat is None else math.hypot(stat, model))
        self.coverage = coverage
        self.values = values

    @property
    def min(self):
        return min(v[2] for v in self.values)

    @property
    def max(self):
        return max(v[2] for v in self.values)

    def to_dict(self):
        return dict(name=self.name, mean=self.mean, stat=self.stat,
                    model=self.model, total=self.total, coverage=self.coverage,
                    min=self.min, max=self.max,
                    values=[dict(label=l, weight=w, value=v, sigma=s)
                            for l, w, v, s in self.values])

    def __repr__(self):
        s = "?" if self.stat is None else f"{self.stat:.4g}"
        t = "?" if self.total is None else f"{self.total:.4g}"
        return (f"Averaged({self.name}: {self.mean:.6g} stat {s} model "
                f"{self.model:.4g} total {t})")


def _average(name, items):
    """items: [(label, weight, value, sigma or None)] with weight > 0 or not."""
    wsum = sum(w for _, w, _, _ in items)
    if wsum <= 0:
        return None
    mean = sum(w * v for _, w, v, _ in items) / wsum
    model2 = sum(w * (v - mean) ** 2 for _, w, v, _ in items) / wsum
    if any(s is None for _, _, _, s in items):
        stat = None
    else:
        stat = math.sqrt(sum(w * s * s for _, w, _, s in items) / wsum)
    return Averaged(name, mean, stat, math.sqrt(max(model2, 0.0)), wsum,
                    [(l, w / wsum, v, s) for l, w, v, s in items])


class ModelAverage:
    """The result of :func:`model_average`.

    ``variants``    the input :class:`Variant` list
    ``method``      the weighting, ``rescale`` the factor applied to chi2
    ``weights``     one dict per variant, in input order: label, chi2, N, k,
                    chi2/nu, ic, delta (IC - min IC), prior, weight
    ``parameters``  ``{name: Averaged}`` over every name in ``values``
    ``derived``     ``{name: Averaged}`` over every derived quantity
    """

    def __init__(self, variants, method, rescale, weights, parameters, derived):
        self.variants = variants
        self.method = method
        self.rescale = rescale
        self.weights = weights
        self.parameters = parameters
        self.derived = derived

    def weight(self, label):
        for w in self.weights:
            if w["label"] == label:
                return w["weight"]
        raise KeyError(label)

    @property
    def values(self):
        """``{name: mean}`` of the averaged parameters (for the writer)."""
        return {k: a.mean for k, a in self.parameters.items()}

    def covariance(self, names=None):
        """``(names, stat, model, total)`` averaged covariance matrices over
        ``names`` (default: the parameters every variant has).  The stat part
        is ``sum_i w_i C_i`` (a variant without a covariance contributes the
        diagonal of its errors; one without errors raises), the model part
        ``sum_i w_i (m_i - m)(m_i - m)^T``."""
        if names is None:
            names = sorted(set.intersection(
                *[set(v.values) for v in self.variants])) if self.variants else []
        names = list(names)
        n = len(names)
        mean = [self.parameters[k].mean for k in names]
        stat = [[0.0] * n for _ in range(n)]
        model = [[0.0] * n for _ in range(n)]
        wsum = 0.0
        for v, wrow in zip(self.variants, self.weights):
            w = wrow["weight"]
            if w <= 0:
                continue
            if any(k not in v.values for k in names):
                raise ValueError(f"variant {v.label} lacks some of {names}")
            wsum += w
            d = [v.values[k] - m for k, m in zip(names, mean)]
            c = None
            if v.covariance is not None:
                idx = {k: i for i, k in enumerate(v.covariance[0])}
                if all(k in idx for k in names):
                    c = [[v.covariance[1][idx[a]][idx[b]] for b in names]
                         for a in names]
            if c is None:
                if any(k not in v.errors for k in names):
                    raise ValueError(f"variant {v.label} has no statistical "
                                     f"error for some of {names}")
                c = [[(v.errors[a] ** 2 if i == j else 0.0)
                      for j, _ in enumerate(names)] for i, a in enumerate(names)]
            for i in range(n):
                for j in range(n):
                    stat[i][j] += w * c[i][j]
                    model[i][j] += w * d[i] * d[j]
        if wsum > 0:
            stat = [[x / wsum for x in r] for r in stat]
            model = [[x / wsum for x in r] for r in model]
        total = [[stat[i][j] + model[i][j] for j in range(n)] for i in range(n)]
        return names, stat, model, total

    # -- output ---------------------------------------------------------------

    def weights_table(self):
        """The weights as text, one row per variant."""
        rows = [f"{'variant':<28} {'chi2':>11} {'N':>5} {'k':>4} {'chi2/nu':>8} "
                f"{'dIC':>9} {'prior':>6} {'weight':>8}"]
        for w in self.weights:
            rn = "-" if w["chi2_nu"] is None else f"{w['chi2_nu']:.3f}"
            rows.append(f"{w['label']:<28} {w['chi2']:11.3f} {w['npoints']:5d} "
                        f"{w['nfree']:4d} {rn:>8} {w['delta']:9.3f} "
                        f"{w['prior']:6.3g} {w['weight']:8.4f}")
        return "\n".join(rows)

    def table(self, names=None, derived=True):
        """Averaged quantities as text: mean, stat, model, total, min, max."""
        def fmt(x):
            return "-" if x is None else f"{x:.5g}"
        items = []
        if names is None:
            items = list(self.parameters.values())
        else:
            items = [self.parameters[n] for n in names if n in self.parameters]
        if derived:
            items += list(self.derived.values())
        width = max([len(a.name) for a in items] + [8])
        rows = [f"{'quantity':<{width}} {'mean':>12} {'stat':>11} {'model':>11} "
                f"{'total':>11} {'min':>12} {'max':>12} {'cover':>6}"]
        for a in items:
            rows.append(f"{a.name:<{width}} {fmt(a.mean):>12} {fmt(a.stat):>11} "
                        f"{fmt(a.model):>11} {fmt(a.total):>11} {fmt(a.min):>12} "
                        f"{fmt(a.max):>12} {a.coverage:6.3f}")
        return "\n".join(rows)

    def summary(self):
        r = self.rescale
        head = (f"Model average over {len(self.variants)} variants, weights "
                f"{self.method}" + ("" if r in (None, 1.0) else
                                   f", chi2 divided by {r:.4g}"))
        return "\n\n".join([head, self.weights_table(), self.table()])

    def to_dict(self):
        return dict(method=self.method, rescale=self.rescale,
                    weights=self.weights,
                    parameters={k: a.to_dict() for k, a in self.parameters.items()},
                    derived={k: a.to_dict() for k, a in self.derived.items()},
                    variants=[v.to_dict() for v in self.variants])

    def write_json(self, path):
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=1)
        return path

    def write_csv(self, path):
        """One row per averaged quantity: kind, name, mean, stat, model,
        total, min, max, coverage."""
        with open(path, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["kind", "name", "mean", "stat", "model", "total",
                         "min", "max", "coverage"])
            for kind, group in (("parameter", self.parameters),
                                ("derived", self.derived)):
                for a in group.values():
                    wr.writerow([kind, a.name, repr(a.mean),
                                 "" if a.stat is None else repr(a.stat),
                                 repr(a.model),
                                 "" if a.total is None else repr(a.total),
                                 repr(a.min), repr(a.max), repr(a.coverage)])
        return path


def model_average(variants: Sequence[Variant], weights="aic", rescale=None,
                  priors: Optional[Dict[str, float]] = None) -> ModelAverage:
    """Average a list of fitted :class:`Variant` records.

    ``weights``: ``"aic"`` (default), ``"bic"``, ``"chi2"`` or ``"flat"``
    (see the module docstring).  ``rescale``: ``None``, ``"best"`` or a
    number dividing every chi2 before the weights are formed.  ``priors``:
    optional ``{label: prior}`` overriding the variants' own.  Returns a
    :class:`ModelAverage`.
    """
    variants = list(variants)
    if not variants:
        raise ValueError("model_average needs at least one variant.")
    labels = [v.label for v in variants]
    if len(set(labels)) != len(labels):
        raise ValueError("variant labels must be unique.")
    method = str(weights).lower()
    if method not in WEIGHT_METHODS:
        raise ValueError(f"weights must be one of {WEIGHT_METHODS}, got {weights!r}")
    priors = dict(priors or {})
    unknown = set(priors) - set(labels)
    if unknown:
        raise ValueError(f"priors for unknown variants: {sorted(unknown)}")
    if len({v.npoints for v in variants}) > 1 and method != "flat":
        warnings.warn("the variants do not have the same number of points; "
                      "information-criterion weights compare models of the "
                      "same data only.", stacklevel=2)

    def ic(v, s):
        c = v.chi2 / s
        if method == "aic":
            return c + 2.0 * v.nfree
        if method == "bic":
            return c + v.nfree * math.log(max(v.npoints, 1))
        if method == "chi2":
            return c
        return 0.0

    pri = [priors.get(v.label, v.prior) for v in variants]
    if rescale is None:
        s = 1.0
    elif isinstance(rescale, str):
        if rescale.lower() != "best":
            raise ValueError(f"rescale must be None, 'best' or a number, got {rescale!r}")
        # among the variants that take part (a prior-0 variant is excluded)
        best = min((v for v, p in zip(variants, pri) if p > 0),
                   key=lambda v: ic(v, 1.0), default=None)
        s = (max(1.0, best.chi2 / best.dof) if best is not None and best.dof > 0 else 1.0)
    else:
        s = float(rescale)
        if not s > 0:
            raise ValueError("rescale must be > 0")

    ics = [ic(v, s) for v in variants]
    live = [c for c, p in zip(ics, pri) if p > 0]
    if not live:
        raise ValueError("every variant has prior weight 0.")
    icmin = min(live)
    raw = [(p * math.exp(-0.5 * (c - icmin)) if p > 0 else 0.0)
           for c, p in zip(ics, pri)]
    tot = sum(raw)
    w = [r / tot for r in raw]
    table = [dict(label=v.label, chi2=v.chi2, npoints=v.npoints,
                  nfree=v.nfree,
                  chi2_nu=(v.chi2 / v.dof if v.dof > 0 else None),
                  ic=c, delta=c - icmin, prior=p, weight=wi,
                  settings=v.settings)
             for v, c, p, wi in zip(variants, ics, pri, w)]

    def collect(get_values, get_errors):
        names = []
        for v in variants:
            for k in get_values(v):
                if k not in names:
                    names.append(k)
        out = {}
        for k in names:
            # a prior-0 variant is excluded: not in the mean, the
            # stat (a missing error), or the min/max range
            items = [(v.label, wi, get_values(v)[k], get_errors(v).get(k))
                     for v, wi, p in zip(variants, w, pri) if p > 0 and k in get_values(v)]
            a = _average(k, items)
            if a is not None:
                out[k] = a
        return out

    params = collect(lambda v: v.values, lambda v: v.errors)
    derived = collect(lambda v: v.derived, lambda v: v.derived_errors)
    return ModelAverage(variants, method, (None if s == 1.0 else s), table,
                        params, derived)


# ---------------------------------------------------------------------------
#  Writing the averaged parameters into a template .azr
# ---------------------------------------------------------------------------

def _azr_model_class():
    """AzrModel, also when this file is loaded on its own (no package)."""
    if __package__:
        from .azrfile import AzrModel
        return AzrModel
    spec = importlib.util.spec_from_file_location(
        "azrfile", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "azrfile.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.AzrModel


def azr_parameter_slots(model):
    """``{name: setter}`` for every level energy and width of an AzrModel.

    Names as :func:`parameter_label` gives them; levels numbered as the engine
    numbers them (active levels only, file order within a J^pi; see
    :meth:`pyazr.azrfile.AzrModel.engine_level_keys`).  A setter takes the
    *physical* value (level energy in MeV; partial width in eV, or ANC in
    fm^-1/2 for a closed channel), as a ``<levels>`` line holds it.
    """
    slots = {}
    for (_, n), lv in model.engine_level_keys().items():
        slots[f"E[{level_label(lv.J, lv.parity, n)}]"] = lv.set_energy
        for c in lv.channels:
            name = width_label(lv.J, lv.parity, n, c.pair, c.L, c.S)
            slots[name] = (lambda v, c=c: setattr(c, "gamma", float(v)))
    return slots


def write_averaged_azr(template, path, values, output_dir=None):
    """Write the averaged parameter values into a copy of ``template``.

    ``template``: a ``.azr`` path or an AzrModel (not modified).  ``values``:
    a :class:`ModelAverage` (its means are written) or ``{name: value}``.
    Level energies and widths are matched by name (:func:`parameter_label`)
    and set through AzrModel's own setters, so the result opens in the GUI
    and in the engine like any other project; everything else in the file
    (radii, the ``<thm>`` block, data) is the template's.  A name the
    template has no slot for (a norm, a cbkg amplitude, a level the template
    lacks) is skipped with a warning.  ``output_dir`` optionally repoints the
    ``<config>`` output directory.

    The result is a representative model, not a fit: its chi2 is not the
    variants' and need not be close to the best one when the variants
    disagree.  Returns ``(path, written, skipped)``.
    """
    AzrModel = _azr_model_class()
    if isinstance(template, (str, os.PathLike)):
        model = AzrModel.from_file(os.fspath(template))
    else:
        # a private copy, so the caller's model stays untouched
        model = _copy_model(AzrModel, template)
    if isinstance(values, ModelAverage):
        values = values.values
    slots = azr_parameter_slots(model)
    written, skipped = [], []
    for name, v in values.items():
        setter = slots.get(name)
        if setter is None:
            skipped.append(name)
            continue
        setter(v)
        written.append(name)
    if skipped:
        warnings.warn(f"write_averaged_azr: {len(skipped)} name(s) not in the "
                      f"template, skipped: {', '.join(skipped[:10])}"
                      + (" ..." if len(skipped) > 10 else ""), stacklevel=2)
    if output_dir is not None:
        model.set_output_dir(output_dir)
    model.write(os.fspath(path))
    return os.fspath(path), written, skipped


def _copy_model(AzrModel, template):
    """A deep copy of an AzrModel through its own text (no shared levels)."""
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix=".azr")
    os.close(fd)
    try:
        with open(tmp, "w") as fh:
            fh.write(template.to_text())
        model = AzrModel.from_file(tmp)
    finally:
        os.remove(tmp)
    model.source = template.source
    return model
