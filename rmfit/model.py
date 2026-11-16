"""Model variants written as new ``.azr`` files.

The engine reads its model from a file, so every variant -- a shard holding a
subset of the data segments, an energy window, a level switched off -- is a
new ``.azr`` written next to the campaign and opened in a fresh session.  The
original is never modified.

Shards are the workhorse: chi-squared is exactly additive over data segments
and the engine's cost is per data point (and does not parallelize with OpenMP),
so splitting the segments across N single-thread sessions divides the cost of
one objective or Jacobian by ~N.  Big segments are cut into energy pieces so
the shards balance; pieces of one segment share its normalization, which the
:mod:`rmfit.engine` assembler sums back into one column.
"""
import os
import re
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import AZURE2_ROOT  # noqa: F401  (puts pyazr on sys.path)
from .azr3 import SegmentTable as SegmentSet, Segment

DIFFERENTIAL = {"differential", "differential-cm", "analyzing-power"}


# -- text helpers -------------------------------------------------------------

def _block(text: str, tag: str) -> Tuple[int, int]:
    """Character offsets of a block's body: after ``<tag>\\n`` up to ``</tag>``."""
    start = text.find(f"<{tag}>")
    end = text.find(f"</{tag}>")
    if start == -1 or end == -1:
        raise ValueError(f"no <{tag}> block")
    body_start = text.index("\n", start) + 1
    return body_start, end


def block_lines(text: str, tag: str) -> List[str]:
    a, b = _block(text, tag)
    return [l for l in text[a:b].splitlines() if l.strip()]


def replace_block(text: str, tag: str, lines: Sequence[str]) -> str:
    a, b = _block(text, tag)
    body = "".join(l + "\n" for l in lines)
    return text[:a] + body + text[b:]


def set_output_dir_text(text: str, path: str) -> str:
    path = os.path.abspath(path)
    if not path.endswith("/"):
        path += "/"
    os.makedirs(path, exist_ok=True)
    out = []
    done = False
    for line in text.splitlines(keepends=True):
        if not done and line.rstrip().endswith("#Full Path to Output Directory"):
            comment = line[line.index("#"):]
            line = f"{path:<100}{comment}"
            done = True
        out.append(line)
    if not done:
        raise ValueError("no output-directory line in <config>")
    return "".join(out)


# -- data-point counting -------------------------------------------------------

def read_points(path: str) -> np.ndarray:
    """(energy, angle, value, error) of every data row -- tolerant of comments
    and of files with only the first two columns (the rest come back NaN)."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#!":
                continue
            t = s.split()
            try:
                e, a = float(t[0]), float(t[1])
            except (ValueError, IndexError):
                continue
            v = float(t[2]) if len(t) > 2 and _isfloat(t[2]) else np.nan
            u = float(t[3]) if len(t) > 3 and _isfloat(t[3]) else np.nan
            rows.append((e, a, v, u))
    return np.array(rows, float).reshape(-1, 4)


def _isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def rows_in_segment(seg: Segment, base_dir: str) -> np.ndarray:
    """The data rows AZURE2 loads for this segment, sorted by energy
    (columns energy, angle, value, error)."""
    pts = read_points(os.path.join(base_dir, seg.data_file))
    if pts.size == 0:
        return pts.reshape(0, 4)
    e, a = pts[:, 0], pts[:, 1]
    sel = (e >= seg.energy_min) & (e <= seg.energy_max)
    if seg.observable in DIFFERENTIAL:
        tol = 1e-6
        sel &= (a >= seg.angle_min - tol) & (a <= seg.angle_max + tol)
    rows = pts[sel]
    return rows[np.argsort(rows[:, 0], kind="stable")]


def points_in_segment(seg: Segment, base_dir: str) -> np.ndarray:
    """Energies of the data rows AZURE2 loads for this segment (sorted)."""
    return rows_in_segment(seg, base_dir)[:, 0]


# -- target effects ----------------------------------------------------------------

@dataclass
class TargetEffectSpec:
    """One <targetInt> entry (TargetEffect::TargetEffect(std::string))."""
    enabled: bool
    keys: List[int]
    num_points: int
    is_convolution: bool
    sigma: float
    is_target_integration: bool
    conv_coeffs: Optional[List[float]]     # energy-dependent sigma equation coefficients
    conv_equation: str
    line: str
    sigma_lab: bool = False                # `sigmalab` keyword: expression lab in / lab out (AZURE2 dad58dc, 2026-10-11)
    # optional restriction to lab-energy windows (TargetEffect ranges_, transitionWidth_, autoTolerance_; 2026-10-15):
    # points outside every window are ordinary points, a positive blend width ramps each edge with a smoothstep
    ranges: List[Tuple[float, float]] = field(default_factory=list)
    blend_width: float = 0.0
    auto_tolerance: float = 0.0

    def blend_weight(self, e_lab):
        """TargetEffect::BlendWeight: 1 with no ranges; else the max over windows of a hard 0/1 (blend 0) or of the
        product of two smoothstep ramps of width blend_width centred on the window edges."""
        x = np.atleast_1d(np.asarray(e_lab, float))
        if not self.ranges:
            return np.ones_like(x)
        out = np.zeros_like(x)
        for lo, hi in self.ranges:
            if self.blend_width <= 0:
                w = ((x >= lo) & (x <= hi)).astype(float)
            else:
                up = np.clip((x - (lo - 0.5 * self.blend_width)) / self.blend_width, 0.0, 1.0)
                dn = np.clip(((hi + 0.5 * self.blend_width) - x) / self.blend_width, 0.0, 1.0)
                w = up * up * (3 - 2 * up) * dn * dn * (3 - 2 * dn)
            out = np.maximum(out, w)
        return out

    @property
    def pure_convolution(self) -> bool:
        return (self.is_convolution or self.conv_coeffs is not None) and not self.is_target_integration

    def sigma_eq(self, energy):
        """sigma from the convolution expression (AZURE2 equation syntax: a0..aN, x, ^, sqrt/exp/log...);
        falls back to the fixed sigma.  Until 2026-10-12 only sqrt(a2*x^3+a1*x^2+a0*x) was accepted."""
        x = np.asarray(energy, float)
        if self.conv_coeffs is None:
            return np.full_like(x, self.sigma)
        expr = self.conv_equation.replace(" ", "").replace("^", "**")
        ns = {f"a{i}": float(v) for i, v in enumerate(self.conv_coeffs)}
        ns.update(x=x, sqrt=np.sqrt, exp=np.exp, log=np.log, ln=np.log, log10=np.log10, abs=np.abs,
                  sin=np.sin, cos=np.cos, tan=np.tan, pow=np.power, pi=np.pi)
        val = eval(compile(expr, "<targetInt sigma>", "eval"), {"__builtins__": {}}, ns)
        return np.broadcast_to(np.asarray(val, float), x.shape).astype(float)


def _split_quoted(line: str) -> List[str]:
    out, i, n = [], 0, len(line)
    while i < n:
        if line[i].isspace():
            i += 1
        elif line[i] == '"':
            j = line.index('"', i + 1)
            out.append(line[i:j + 1])
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            out.append(line[i:j])
            i = j
    return out


def parse_seg_list(text: str) -> List[int]:
    """The segment keys of a ``<targetInt>`` key field, read as AZURE2 reads it (TargetEffect::GetSegmentsList):
    comma-separated keys and inclusive ranges, ``"41,47-50,54-56"`` -> [41, 47, 48, 49, 50, 54, 55, 56].

    2026-11-07 (the evaluator approved): rmfit used to read only plain comma lists.  A range made ``int("47-50")``
    raise, ``target_effects`` swallowed the error and ``_TINT_RE`` did not match, so every entry with a range
    was left out of every shard and its segments were evaluated UNINTEGRATED in sharded runs (37
    such entries in one model)."""
    out: List[int] = []
    for item in text.strip().strip('"').replace(" ", "").replace("\t", "").split(","):
        if not item:
            continue
        if "-" in item:
            lo, hi = item.split("-", 1)
            a, b = int(lo), int(hi)
            if b < a:
                raise ValueError(f"<targetInt> segment range {item!r} runs backwards")
            out.extend(range(a, b + 1))
        else:
            out.append(int(item))
    return out


def parse_target_effect(line: str) -> TargetEffectSpec:
    t = _split_quoted(line)
    enabled = t[0] == "1"
    keys = parse_seg_list(t[1])
    num_points = int(float(t[2]))
    is_conv = int(float(t[3])) == 1
    sigma = float(t[4])
    is_ti = int(float(t[5])) == 1
    i = 8                                   # density t[6], stopping eq t[7]
    n_par = int(float(t[i])); i += 1 + n_par
    coeffs, eq = None, ""
    ranges, blend, tol = [], 0.0, 0.0
    if i < len(t):
        n_q = int(float(t[i + 1])); i += 2 + n_q     # isQ, numQ, q...
        if i < len(t):
            is_cc = int(float(t[i])) == 1
            eq = t[i + 1].strip('"')
            n_cc = int(float(t[i + 2]))
            vals = [float(v) for v in t[i + 3:i + 3 + n_cc]]
            if is_cc:
                coeffs = vals
            # TargetEffect.cpp: up to four numbers (straggling flag, coefficient, resonanceWidthMultiplier,
            # pointsPerWidth), then optionally a quoted "lo-hi,lo-hi" window list, a blend width and an
            # auto-application tolerance; keyword blocks (sigmalab, udr, beamprofile) follow
            j = i + 3 + n_cc
            k = 0
            while j < len(t) and k < 4 and _isnum(t[j]):
                j += 1; k += 1
            if j < len(t) and t[j].startswith('"'):
                for item in t[j].strip('"').split(","):
                    dash = item.find("-", 1)
                    if dash < 0:
                        continue
                    lo, hi = float(item[:dash]), float(item[dash + 1:])
                    if hi > lo:
                        ranges.append((lo, hi))
                j += 1
                if j < len(t) and _isnum(t[j]):
                    blend = max(0.0, float(t[j])); j += 1
                    if j < len(t) and _isnum(t[j]):
                        tol = max(0.0, float(t[j]))
    return TargetEffectSpec(enabled, keys, num_points, is_conv, sigma, is_ti, coeffs, eq, line,
                            sigma_lab=any(tok.lower() == "sigmalab" for tok in t),
                            ranges=ranges, blend_width=blend, auto_tolerance=tol)


def target_effects(text: str) -> List[TargetEffectSpec]:
    out = []
    for line in block_lines(text, "targetInt"):
        try:
            out.append(parse_target_effect(line))
        except (ValueError, IndexError) as e:
            # an ENABLED entry rmfit cannot read would drop target integration from its segments without a
            # word (the 2026-11-07 range bug): say so instead of skipping quietly
            if line.split()[:1] == ["1"]:
                print(f"rmfit WARNING: <targetInt> entry not understood ({e}); its segments get no target "
                      f"effect in rmfit: {line.strip()[:120]}", flush=True)
            continue
    return out


def cm_conversion(text: str, entrance_key: int) -> float:
    """``E_cm / E_lab`` for a particle entrance pair, from the <levels> lines
    (EData.cpp: m2 / (m1 + m2))."""
    for line in block_lines(text, "levels"):
        t = line.split()
        if len(t) >= 20 and int(float(t[5])) == entrance_key:
            m1, m2 = float(t[17]), float(t[18])
            return m2 / (m1 + m2)
    raise ValueError(f"no <levels> line uses pair {entrance_key}")


# -- shards --------------------------------------------------------------------

@dataclass
class ConvSpec:
    """A pure-convolution segment evaluated on a fine grid instead of through
    AZURE2's per-point sub-point quadrature (which costs ~40x a normal point).

    The shard evaluates the *unconvolved* cross section on ``grid_cm`` (a
    synthetic data file with value 0 and error 1, norm fixed at 1, so the
    returned residuals *are* the fit values); the assembler folds them with
    the same kernel AZURE2 uses (since dad58dc, 2026-10-11: one Gaussian in
    c.m. energy per data point, sigma = TargetEffect::SigmaCM at the point,
    truncated at +-3 sigma and normalized by its own integral, for a constant
    sigma and for an expression alike; see build_conv_spec), and forms the
    real residuals against the segment's data.  Exact up to quadrature (AZURE2
    itself uses an adaptive grid + 2-point Gauss-Legendre with linear
    interpolation); on a large model the difference is ~1e-5 of the
    total objective, and every verified number comes from a single session
    anyway.
    """
    e_lab: np.ndarray            # data energies (lab)
    value: np.ndarray
    error: np.ndarray
    grid_cm: np.ndarray          # fine grid, c.m. energy, ascending
    cm_conv: float
    W: np.ndarray                # (n_data, n_grid) kernel weights

    @property
    def n_data(self):
        return int(self.e_lab.size)

    @property
    def n_grid(self):
        return int(self.grid_cm.size)


def build_conv_spec(seg: Segment, te: TargetEffectSpec, rows: np.ndarray, cm_conv: float,
                    points_per_sigma: int = 10, range_sigma: float = 3.0) -> ConvSpec:
    e_lab, value, error = rows[:, 0], rows[:, 2], rows[:, 3]
    e_cm = e_lab * cm_conv
    # AZURE2 after dad58dc (2026-10-11, "one Gaussian beam spread"): TargetEffect::SigmaCM gives ONE width per data
    # point, evaluated at the point -- a constant is lab and converted once (x cm_conv); an expression with
    # `sigmalab` is evaluated at the point's LAB energy and converted once; a legacy expression (no keyword) is
    # evaluated at the point's c.m. energy as a c.m. width.  Every kernel is normalised by its own integral over
    # +-3 sigma.  (Before 2026-10-12 rmfit reproduced the pre-dad58dc engine: constant-sigma kernels unnormalised,
    # i.e. 0.27 % low, and expressions evaluated per sub-point -- a sharded/single offset on large models.)
    if te.conv_coeffs is not None and te.sigma_lab:
        s_pt = te.sigma_eq(e_lab) * cm_conv
    elif te.conv_coeffs is not None:
        s_pt = te.sigma_eq(e_cm)
    else:
        s_pt = np.full(e_cm.size, te.sigma * cm_conv)
    s_pt = np.asarray(s_pt, float)
    if te.ranges:
        return _windowed_conv_spec(te, e_lab, value, error, e_cm, s_pt, cm_conv, points_per_sigma, range_sigma)
    if not np.all(s_pt > 0):
        raise ValueError(f"<targetInt> {te.keys}: non-positive Gaussian width at some points")
    lo = float(np.min(e_cm - range_sigma * s_pt))
    hi = float(np.max(e_cm + range_sigma * s_pt))
    h = float(np.min(s_pt)) / points_per_sigma
    grid = np.arange(lo, hi + h, h)
    dg = np.gradient(grid)                       # trapezoid-like weights
    W = np.zeros((e_cm.size, grid.size))
    for i, (ec, sp) in enumerate(zip(e_cm, s_pt)):
        inside = np.abs(grid - ec) <= range_sigma * sp
        w = np.where(inside, np.exp(-0.5 * ((grid - ec) / sp) ** 2) * dg, 0.0)
        W[i] = w / np.sum(w)
    return ConvSpec(e_lab, value, error, grid, cm_conv, W)


def _windowed_conv_spec(te, e_lab, value, error, e_cm, s_pt, cm_conv, points_per_sigma, range_sigma) -> ConvSpec:
    """A <targetInt> restricted to lab-energy windows (2026-10-15; AZURE2 TargetEffect ranges + EPoint blend):
    model = b * (Gaussian-folded) + (1 - b) * (bare at the point), b = te.blend_weight(E_lab), so a point outside every
    window is the bare cross section at its own energy.
    Grid: every data point's own c.m. energy, plus one UNIFORM lattice per cluster of windowed points (step = the
    cluster's smallest sigma / points_per_sigma, spanning the cluster's +-range_sigma sigma), instead of one global
    min(sigma)/10 grid over the whole segment.  Kernel: per windowed point a symmetric stencil e + k sigma/points_per_sigma,
    |k| <= points_per_sigma * range_sigma, Gaussian weights normalised to 1, each node linearly interpolated onto the grid
    -- exact for a cross section linear in E (mean preserved), width sigma (truncated at range_sigma)."""
    b = te.blend_weight(e_lab)
    conv = np.nonzero(b > 0)[0]
    if conv.size and np.any(~(s_pt[conv] > 0)):
        raise ValueError(f"<targetInt> {te.keys}: non-positive Gaussian width at a windowed point")
    parts = [e_cm]
    order = conv[np.argsort(e_cm[conv])]
    cl = []
    for i in order:                                  # clusters of overlapping +-range_sigma sigma intervals
        lo, hi = e_cm[i] - range_sigma * s_pt[i], e_cm[i] + range_sigma * s_pt[i]
        if cl and lo <= cl[-1][1]:
            cl[-1][1] = max(cl[-1][1], hi); cl[-1][2] = min(cl[-1][2], s_pt[i])
        else:
            cl.append([lo, hi, s_pt[i]])
    for lo, hi, smin in cl:
        h = smin / points_per_sigma
        n = int(np.ceil((hi - lo) / h)) + 1
        parts.append(lo + h * np.arange(n + 1))
    grid = np.unique(np.concatenate(parts))
    W = np.zeros((e_cm.size, grid.size))
    own = np.searchsorted(grid, e_cm)
    if not np.all(grid[own] == e_cm):
        raise RuntimeError("windowed conv grid lost a data energy")
    K = int(round(points_per_sigma * range_sigma))
    k = np.arange(-K, K + 1)
    gk = np.exp(-0.5 * (k / points_per_sigma) ** 2); gk /= gk.sum()
    for i in conv:
        x = e_cm[i] + k * (s_pt[i] / points_per_sigma)
        j = np.clip(np.searchsorted(grid, x, side="right") - 1, 0, grid.size - 2)
        t = (x - grid[j]) / (grid[j + 1] - grid[j])
        np.add.at(W[i], j, b[i] * gk * (1.0 - t))
        np.add.at(W[i], j + 1, b[i] * gk * t)
    W[np.arange(e_cm.size), own] += 1.0 - b
    return ConvSpec(e_lab, value, error, grid, cm_conv, W)


@dataclass
class Piece:
    base_key: int
    tokens: List[str]           # the <segmentsData> tokens, energy window edited
    npts: int                   # rows this piece contributes to the assembled residual
    cost: float
    cls: str
    e_lo: Optional[float] = None
    e_hi: Optional[float] = None
    kind: str = "data"          # data | conv
    conv: Optional[ConvSpec] = None
    data_file: Optional[str] = None   # synthetic data file (conv pieces)


@dataclass
class ShardSpec:
    name: str
    azr: str
    output_dir: str
    pieces: List[Piece] = field(default_factory=list)
    predicted: Optional[float] = None   # Jacobian-balanced plans: the predicted resid+Jacobian time (s)

    @property
    def seg_map(self) -> Dict[int, int]:
        """shard segment key (1-based, file order) -> base segment key"""
        return {i + 1: p.base_key for i, p in enumerate(self.pieces)}

    @property
    def npts(self):
        return sum(p.npts for p in self.pieces)

    @property
    def cost(self):
        return sum(p.cost for p in self.pieces)


def free_shift(seg: Segment) -> bool:
    """Any free energy-shift term: the constant, the sqrt(E) coefficient or a shift-function parameter.
    Such a segment is finite-differenced in the Jacobian, so it is never split and gets a shift-solo shard."""
    return bool(seg.vary_shift or getattr(seg, "vary_shift_sqrt", False)
                or any(v for _, _, v in getattr(seg, "shift_fn_params", ())))


def segment_class(seg: Segment) -> str:
    return f"{seg.entrance_key}->{seg.exit_key}:{seg.observable}"


#: Relative cost per data point by segment class, calibrated on a large
#: model (rmfit.model.calibrate_weights): fixed-angle differential segments
#: cost more per point than angle-integrated ones.
DEFAULT_WEIGHTS = {"2->2:differential": 8.5, "2->1:differential": 2.3}


_TINT_RE = re.compile(r'^(\s*\d+\s+)"([\d,\s\-]*)"(.*)$', re.S)     # ranges ("47-50") read since 2026-11-07


def target_int_entries(text: str) -> List[Tuple[bool, List[int], str]]:
    """The <targetInt> entries as ``(enabled, [segment keys], line)``.

    An entry's key field is a quoted list of segment keys and inclusive
    ranges (``"8,9"``, ``"41,47-50"``; :func:`parse_seg_list`), matched by
    number against *both* the data and the test segments
    (``EData::ReadTargetEffectsFile``).  :func:`retarget_int_line` writes
    plain comma lists, which AZURE2 reads the same way.
    """
    out = []
    for line in block_lines(text, "targetInt"):
        m = _TINT_RE.match(line)
        if not m:
            if line.split()[:1] == ["1"]:
                print(f"rmfit WARNING: <targetInt> key field not understood; entry left out of the shards: "
                      f"{line.strip()[:120]}", flush=True)
            continue
        enabled = m.group(1).split()[0] == "1"
        keys = parse_seg_list(m.group(2))
        out.append((enabled, keys, line))
    return out


def target_int_keys(text: str) -> Dict[str, str]:
    """``{segment key: targetInt line}`` for the enabled <targetInt> entries."""
    out = {}
    for enabled, keys, line in target_int_entries(text):
        if enabled:
            for k in keys:
                out[str(k)] = line
    return out


def retarget_int_line(line: str, new_keys: Sequence[int]) -> str:
    m = _TINT_RE.match(line)
    return f'{m.group(1)}"{",".join(str(k) for k in new_keys)}"{m.group(3)}'


def plan_pieces(base_azr: str, weights: Optional[Dict[str, float]] = None,
                target_cost: Optional[float] = None,
                no_split: Sequence[int] = (), conv_grid: bool = True,
                split_composite: Optional[bool] = None,
                seg_factors: Optional[Dict[int, float]] = None, default_factor: float = 1.0) -> List[Piece]:
    """One :class:`Piece` per active segment, big ones cut into energy pieces.

    A segment is never split if it carries a free energy shift, an enabled
    target-effect entry, or is a ratio composite -- or if its key is in
    ``no_split``.  A composite "sum" segment is split like a plain one
    (``split_composite``, default on; RMFIT_SPLIT_SUM=0 turns it off) and
    costs one calculation per component (2026-10-10: a merged model's
    neutron total cross section, a sum over 8 exits with ~13k points, sat whole in one
    shard and set the sharded objective at 25 s; split by component cost into
    15 pieces it is 3 s, the objective bit-identical).  A segment whose target effect is a pure
    convolution becomes a ``conv`` piece (see :class:`ConvSpec`) unless
    ``conv_grid`` is False.
    """
    text = open(base_azr).read()
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    weights = dict(DEFAULT_WEIGHTS, **(weights or {}))
    segs = SegmentSet.from_file(base_azr)
    raw = block_lines(text, "segmentsData")
    if len(raw) != len(segs):
        raise RuntimeError("segment count mismatch between parser and text")
    if split_composite is None:
        split_composite = os.environ.get("RMFIT_SPLIT_SUM", "1") != "0"
    tint = target_int_keys(text)
    effects = {k: te for te in target_effects(text) if te.enabled for k in te.keys}
    # 2026-10-22: a conv piece evaluates the bare cross section on an UNSHIFTED grid (its line gets shift 0) and
    # folds it at the unshifted data energies, so any energy shift -- a fixed nonzero constant, a sqrtshift, a
    # shiftfn, a shift group's -- was silently dropped; such segments stay plain pieces (the engine shifts them)
    from .covariance import split_keyword_tail
    shift_kw = {"sqrtshift", "shiftfn", "shiftgroup", "shiftfngroup"}
    pieces: List[Piece] = []
    for seg, line in zip(segs, raw):
        if not seg.active:
            continue
        rows = rows_in_segment(seg, base_dir)
        energies = rows[:, 0]
        n = int(energies.size)
        cls = segment_class(seg)
        w = weights.get(cls, 1.0)
        if seg_factors is not None:            # measured cost correction (Jacobian balancing): seconds per cost unit
            w *= float(seg_factors.get(int(seg.key), default_factor))
        if split_composite and seg.operation == "sum" and getattr(seg, "components", None):
            w *= max(1, len(seg.components))       # one calculation per component
        cost = n * w
        tokens = line.split()
        te = effects.get(seg.key)
        # a composite "sum" segment (the archive's total cross sections are
        # 1->1 plus 1->2 components) keeps its component list on the
        # synthetic grid line, so it may be a conv piece too
        # an auto-application tolerance (engine-side per-point decision) has no grid counterpart: such a segment
        # stays a plain piece and the shard's own engine applies its <targetInt> line
        shifted = (seg.vary_shift or seg.energy_shift != 0
                   or any(kw in shift_kw for kw, _ in split_keyword_tail(tokens)[1]))
        if (te is not None and te.pure_convolution and not shifted and not te.auto_tolerance
                and seg.operation in (None, "sum") and n > 0 and conv_grid):
            cmc = cm_conversion(text, seg.entrance_key)
            spec = build_conv_spec(seg, te, rows, cmc,
                                   points_per_sigma=float(os.environ.get("RMFIT_CONV_PPS", "10")))   # 2026-10-15: grid density, default unchanged
            pieces.append(Piece(seg.key, tokens, n, spec.n_grid * w * 1.2, cls,
                                seg.energy_min, seg.energy_max, kind="conv", conv=spec))
            continue
        splittable = (seg.key not in no_split and not free_shift(seg)
                      and str(seg.key) not in tint
                      and (seg.operation is None or (split_composite and seg.operation == "sum"))
                      and n >= 40)
        if target_cost is None or cost <= target_cost or not splittable:
            pieces.append(Piece(seg.key, tokens, n, cost, cls, seg.energy_min, seg.energy_max))
            continue
        k = int(np.ceil(cost / target_cost))
        # cut between distinct energies so no point straddles two pieces
        uniq = np.unique(energies)
        k = min(k, uniq.size)
        cuts = [uniq[int(round(j * uniq.size / k))] for j in range(1, k)]
        edges = [seg.energy_min]
        for c in cuts:
            prev = uniq[uniq < c]
            if prev.size == 0:
                continue
            edges.append(0.5 * (prev[-1] + c))
        edges.append(seg.energy_max)
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi <= lo:
                continue
            t = list(tokens)
            t[3] = repr(float(lo))
            t[4] = repr(float(hi))
            m = int(np.sum((energies >= lo) & (energies <= hi)))
            if m == 0:
                continue
            pieces.append(Piece(seg.key, t, m, m * w, cls, lo, hi))
    return pieces


def pack(pieces: List[Piece], nshards: int) -> List[List[Piece]]:
    """Longest-processing-time greedy bin packing by cost."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    load = np.zeros(nshards)
    for p in sorted(pieces, key=lambda p: -p.cost):
        j = int(np.argmin(load))
        bins[j].append(p)
        load[j] += p.cost
    # keep file order inside a shard (stable row order for bookkeeping)
    return [sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0))
            for b in bins]


def pack_units(units: List[List[Piece]], nshards: int) -> List[List[Piece]]:
    """:func:`pack` for units that must share a shard (a shift group's pieces): LPT by unit cost."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    load = np.zeros(nshards)
    for u in sorted(units, key=lambda u: -sum(p.cost for p in u)):
        j = int(np.argmin(load))
        bins[j].extend(u)
        load[j] += sum(p.cost for p in u)
    return [sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0))
            for b in bins]



# -- data directories ----------------------------------------------------------

def data_dirs(azr_path: str) -> List[str]:
    """Directories, relative to the model, that hold the files its segment lines read.

    Every ``<segmentsData>`` / ``<segmentsTest>`` line (active or not) names its file relative to the ``.azr``;
    the first path component of each such name that is an existing directory beside the model is returned, with
    ``data`` always first.  Before 2026-11-07 every rmfit working directory (shards, single-session variants,
    windows, level-scan and mixture members) linked ``data/`` only, so a model reading ``akdata/`` or ``thermal/``
    lost those segments there (4 of 75 in one fit; pyazr then refused the segment count)."""
    base_dir = os.path.dirname(os.path.abspath(azr_path))
    text = open(azr_path, encoding="latin-1").read()
    out = ["data"]
    for blk in ("segmentsData", "segmentsTest"):
        m = re.search(rf"<{blk}>(.*?)</{blk}>", text, re.S)
        if not m:
            continue
        for line in m.group(1).splitlines():
            for tok in line.split():
                tok = tok.strip('"')
                if "/" not in tok or os.path.isabs(tok):
                    continue
                first = os.path.normpath(tok).split(os.sep)[0]
                if first in ("", ".", "..") or first in out:
                    continue
                if os.path.isdir(os.path.join(base_dir, first)):
                    out.append(first)
    return out


def link_data_dirs(azr_path: str, workdir: str, src_dir: Optional[str] = None) -> List[str]:
    """Link every directory of :func:`data_dirs` from ``src_dir`` (default: the model's own directory) into
    ``workdir``, so a copy of the model written there finds all its files.  Existing entries are left alone.
    Returns the names linked or already present."""
    src_dir = os.path.abspath(src_dir or os.path.dirname(os.path.abspath(azr_path)))
    os.makedirs(workdir, exist_ok=True)
    done = []
    for sub in data_dirs(azr_path):
        src, dst = os.path.join(src_dir, sub), os.path.join(workdir, sub)
        if os.path.isdir(src) and not os.path.lexists(dst):
            os.symlink(src, dst)
        if os.path.lexists(dst):
            done.append(sub)
    return done


def write_shard(base_azr: str, name: str, pieces: List[Piece], workdir: str) -> ShardSpec:
    """Write ``<workdir>/<name>.azr`` holding only ``pieces`` as its data."""
    from .covariance import NormGroups, norm_field, strip_keywords
    text = open(base_azr).read()
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    groups = NormGroups.from_azr(base_azr)
    os.makedirs(workdir, exist_ok=True)
    # data files are referenced relative to the .azr: link every directory they live in (data/ and any other)
    link_data_dirs(base_azr, workdir)
    lines = []
    for i, p in enumerate(pieces, start=1):
        if p.kind == "conv":
            # synthetic grid segment: unconvolved cross section at the grid
            # energies, data 0 / error 1 / norm fixed 1 -> residual == fit
            c = p.conv
            fname = f"conv_seg{p.base_key}_{name}.dat"
            # the angle column must pass the segment's angle filter: a fixed-angle (differential) segment keeps only
            # rows at its angle, so the grid rows carry it (fixed-angle conv pieces once came out empty with 0.0);
            # an angle range (angle-integrated data, 0-180) keeps 0.0 as before
            a_lo, a_hi = float(p.tokens[5]), float(p.tokens[6])
            ang = a_lo if a_lo == a_hi else 0.0
            with open(os.path.join(workdir, fname), "w") as fh:
                for e in c.grid_cm / c.cm_conv:
                    fh.write(f"{e:.10f}  {ang!r}  0.0  1.0\n" if ang else f"{e:.10f}  0.0  0.0  1.0\n")
            t = list(p.tokens)
            e_lab = c.grid_cm / c.cm_conv
            t[3] = repr(float(e_lab[0] * (1 - 1e-9)))
            t[4] = repr(float(e_lab[-1] * (1 + 1e-9)))
            # norm 1, fixed, no error; no shift; data file
            isdiff = int(float(t[7]))
            j = 8 + (2 if isdiff == 2 else 0)
            t[j], t[j + 1], t[j + 2] = "1", "0", "0"
            k = j + 3
            while k < len(t) and _isnum(t[k]):
                k += 1
            t = t[:j + 3] + ["0", "0", "0", fname] + t[k + 1:]
            p.data_file = os.path.join(workdir, fname)
            lines.append("  ".join(strip_keywords(t, ("covariance", "normgroup"))))
        else:
            # 2026-10-22: no `covariance` (the evaluator whitens the assembled rows) and no `normgroup`
            # (a group's members may sit in different shards): a follower keeps its own norm, fixed,
            # and the evaluator rescales its rows to the leader's (rmfit.covariance)
            t = strip_keywords(p.tokens, ("covariance", "normgroup"))
            if p.base_key in groups.leader:
                t[norm_field(t) + 1] = "0"
            lines.append("  ".join(t))
    text = replace_block(text, "segmentsData", lines)
    text = replace_block(text, "segmentsTest", [])
    # every enabled target-effect entry that touches a data piece of this
    # shard, re-keyed to the shard's own segment numbers (pieces of one split
    # segment all inherit its entry); conv pieces carry none
    new_tint = []
    for enabled, keys, line in target_int_entries(text):
        if not enabled:
            continue
        mine = [i for i, p in enumerate(pieces, start=1) if p.base_key in keys and p.kind == "data"]
        if mine:
            new_tint.append(retarget_int_line(line, mine))
    text = replace_block(text, "targetInt", new_tint)
    # measured scattering lengths are a property of the whole model, not of a shard: the evaluator
    # carries them as residual rows once (engine.ScatLenRows, 2026-11-07); a shard that kept the block
    # would add the engine's term once per shard on the scalar path
    if "<scatteringLengths>" in text and "</scatteringLengths>" in text:
        text = replace_block(text, "scatteringLengths", [])
    outdir = os.path.join(workdir, f"output_{name}")
    text = set_output_dir_text(text, outdir)
    for f in ("intEC.dat", "intEC.extrap"):
        fp = os.path.join(outdir, f)
        if os.path.exists(fp):
            os.remove(fp)
    path = os.path.join(workdir, f"{name}.azr")
    with open(path, "w") as fh:
        fh.write(text)
    return ShardSpec(name=name, azr=path, output_dir=outdir, pieces=list(pieces))


def _isnum(s):
    return _isfloat(s)


def shift_columns(seg: Segment) -> int:
    """Number of finite-differenced (shift-type) Jacobian columns a segment carries: the constant shift, the
    sqrt(E) coefficient and each free shift-function parameter."""
    return (int(bool(seg.vary_shift)) + int(bool(getattr(seg, "vary_shift_sqrt", False)))
            + sum(1 for _, _, v in getattr(seg, "shift_fn_params", ()) if v))


def shard_shift_count(base_azr: str, pieces: Sequence[Piece], segs=None) -> int:
    """Shift-type columns a shard holding ``pieces`` finite-differences (each segment counted once)."""
    segs = segs if segs is not None else SegmentSet.from_file(base_azr)
    return sum(shift_columns(segs.by_key(k)) for k in sorted({int(p.base_key) for p in pieces}))


def pack_jac(units: List[List[Piece]], nshards: int, alpha: float, beta: float,
             kcount) -> Tuple[List[List[Piece]], List[float]]:
    """Greedy packing for the Jacobian (2026-10-26, the evaluator: "go ahead with the Jacobian-cost shard balancing").

    A shard's resid+Jacobian time is modelled as ``R (alpha + beta k)``: R its data cost (s, one residual
    evaluation), k its finite-differenced shift columns (each re-evaluates the WHOLE shard, ~beta R).  Units with
    shifts are placed first (largest first), then the rest; each unit goes where the resulting time is smallest, so
    the shifts spread out and the shards holding them stay light.  Returns the bins (file order) and their predicted
    times."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    R = np.zeros(nshards)
    K = np.zeros(nshards)
    cost = lambda u: float(sum(p.cost for p in u))
    order = sorted(units, key=lambda u: (kcount(u) == 0, -cost(u)))
    for u in order:
        r, k = cost(u), kcount(u)
        t = (R + r) * (alpha + beta * (K + k))
        j = int(np.argmin(t))
        bins[j].extend(u)
        R[j] += r
        K[j] += k
    loads = list(R * (alpha + beta * K))
    return ([sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0)) for b in bins], loads)


def build_shards(base_azr: str, nshards: int, workdir: str,
                 weights: Optional[Dict[str, float]] = None,
                 no_split: Sequence[int] = (), conv_grid: bool = True,
                 shift_solo: Optional[bool] = None, balance: Optional[dict] = None) -> List[ShardSpec]:
    """Plan, pack and write ``nshards`` shard files of ``base_azr``.

    ``balance`` (2026-10-26): a measured Jacobian cost model ``{"alpha", "beta", "factors": {segment key: s per
    cost unit}, "default_factor"}`` (Sharded measures it, ``engine.Sharded._rebalance``); the pieces are then
    costed in seconds and packed by :func:`pack_jac` (shift_solo is not used).

    ``shift_solo`` (default on; RMFIT_SHIFT_SOLO=0 turns it off): every
    segment with a free energy shift gets a shard of its own, the rest is
    packed into the others.  A shift is finite-differenced in the Jacobian
    and each step re-evaluates its whole shard (2026-10-10, a large merged
    model, 16 shards: shards holding the 6 shifted segments took 28-92 s of
    a 90 s Jacobian, the others 10-17 s; solo, the Jacobian is ~30 s).
    Used only while the shifted segments are fewer than half the shards."""
    if balance is not None:
        return _build_shards_balanced(base_azr, nshards, workdir, weights, no_split, conv_grid, balance)
    if shift_solo is None:
        shift_solo = os.environ.get("RMFIT_SHIFT_SOLO", "1") != "0"
    pieces = plan_pieces(base_azr, weights=weights, conv_grid=conv_grid)
    total = sum(p.cost for p in pieces)
    target = total / nshards
    from .covariance import shift_groups, shift_fn_groups
    sg = shift_groups(base_azr).units()
    # 2026-10-24: a shift-function group's members use the function of its first active member, which declares it;
    # a split member would declare it twice in one shard (an engine error) or not at all, so members are never split
    fg = shift_fn_groups(base_azr).units()
    both = sorted(set(sg) & set(fg))
    if both:
        raise ValueError(f"segments {both} carry both shiftgroup and shiftfngroup (the engine refuses that)")
    pieces = plan_pieces(base_azr, weights=weights, target_cost=target * 0.5,
                         no_split=tuple(no_split) + tuple(fg), conv_grid=conv_grid)
    unit_of = {k: ("g", l) for k, l in sg.items()}
    unit_of.update({k: ("f", l) for k, l in fg.items()})
    if unit_of:
        # 2026-10-22: a shift group's members follow the constant shift of its first member (the engine, within
        # one file), so all their pieces go to ONE shard, which keeps the keyword: the group is exact there, and a
        # finite-difference step of the group's shift moves every member.  A group with a free shift is one solo unit.
        # 2026-10-24: the same for a shift-function group (shiftfngroup); bins keep file order, so the declaring
        # member stays the first of its group in the shard file.
        units: Dict[object, List[Piece]] = {}
        for p in pieces:
            units.setdefault(unit_of.get(p.base_key, ("p", id(p))), []).append(p)
        units_l = list(units.values())
        solo_u = []
        if shift_solo:
            shift_keys = {s.key for s in SegmentSet.from_file(base_azr) if s.active and free_shift(s)}
            solo_u = [u for u in units_l if any(p.base_key in shift_keys for p in u)]
        if solo_u and len(solo_u) <= nshards // 2:
            rest = [u for u in units_l if all(u is not v for v in solo_u)]
            bins = [sorted(u, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0)) for u in solo_u] + \
                pack_units(rest, nshards - len(solo_u))
        else:
            bins = pack_units(units_l, nshards)
    else:
        solo = []
        if shift_solo:
            shift_keys = {s.key for s in SegmentSet.from_file(base_azr) if s.active and free_shift(s)}
            solo = [p for p in pieces if p.base_key in shift_keys]
        if solo and len(solo) <= nshards // 2:
            rest = [p for p in pieces if p.base_key not in shift_keys]
            bins = [[p] for p in solo] + pack(rest, nshards - len(solo))
        else:
            bins = pack(pieces, nshards)
    shards = []
    for j, b in enumerate(bins):
        if not b:
            continue
        shards.append(write_shard(base_azr, f"shard{j:02d}", b, workdir))
    return shards


def plan_balanced(base_azr: str, nshards: int, weights, no_split, conv_grid, balance: dict):
    """The bins and predicted times of a Jacobian-balanced plan (no files written)."""
    from .covariance import shift_groups, shift_fn_groups
    f = {int(k): float(v) for k, v in (balance.get("factors") or {}).items()}
    dflt = float(balance.get("default_factor", 1.0))
    alpha, beta = float(balance["alpha"]), float(balance["beta"])
    sg = shift_groups(base_azr).units()
    fg = shift_fn_groups(base_azr).units()
    p0 = plan_pieces(base_azr, weights=weights, conv_grid=conv_grid, seg_factors=f, default_factor=dflt)
    total = sum(p.cost for p in p0)
    frac = float(os.environ.get("RMFIT_BALANCE_PIECE", "0.35"))      # split target as a fraction of a shard's share
    pieces = plan_pieces(base_azr, weights=weights, target_cost=frac * total / nshards,
                         no_split=tuple(no_split) + tuple(fg), conv_grid=conv_grid, seg_factors=f,
                         default_factor=dflt)
    segs = SegmentSet.from_file(base_azr)
    shifted = {s.key for s in segs if s.active and free_shift(s)}
    unit_of = {k: ("g", l) for k, l in sg.items()}
    unit_of.update({k: ("f", l) for k, l in fg.items()})
    units: Dict[object, List[Piece]] = {}
    for p in pieces:
        key = unit_of.get(p.base_key) or (("s", p.base_key) if p.base_key in shifted else ("p", id(p)))
        units.setdefault(key, []).append(p)
    kc = lambda u: shard_shift_count(base_azr, u, segs)
    return pack_jac(list(units.values()), nshards, alpha, beta, kc)


def _build_shards_balanced(base_azr, nshards, workdir, weights, no_split, conv_grid, balance):
    bins, loads = plan_balanced(base_azr, nshards, weights, no_split, conv_grid, balance)
    shards = []
    for j, (b, t) in enumerate(zip(bins, loads)):
        if not b:
            continue
        spec = write_shard(base_azr, f"shard{j:02d}", b, workdir)
        spec.predicted = float(t)
        shards.append(spec)
    return shards


def calibrate_weights(shards: List[ShardSpec], times: Sequence[float]) -> Dict[str, float]:
    """Per-class cost per point from measured shard evaluation times
    (non-negative least squares), normalized so the cheapest class is 1."""
    from scipy.optimize import nnls
    classes = sorted({p.cls for s in shards for p in s.pieces})
    A = np.zeros((len(shards), len(classes)))
    for i, s in enumerate(shards):
        for p in s.pieces:
            A[i, classes.index(p.cls)] += p.npts
    c, _ = nnls(A, np.asarray(times, float))
    c = np.where(c > 0, c, np.nan)
    base = np.nanmin(c) if np.isfinite(c).any() else 1.0
    return {cls: (float(v / base) if np.isfinite(v) else 1.0) for cls, v in zip(classes, c)}


# -- other variants ------------------------------------------------------------

def single_variant(base_azr: str, name: str, workdir: str) -> str:
    """A copy of the base model with its own output directory (for an
    authoritative single session that must not touch the campaign's output/)."""
    text = open(base_azr).read()
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    os.makedirs(workdir, exist_ok=True)
    link_data_dirs(base_azr, workdir)
    # a covariance file named relative to the model but outside data/ is linked too
    from .covariance import line_keyword, segment_lines
    for t in segment_lines(base_azr):
        c = line_keyword(t, "covariance")
        if c and not os.path.isabs(c) and not os.path.normpath(c).startswith(("data" + os.sep, "..")):
            src, lnk = os.path.join(base_dir, c), os.path.join(workdir, c)
            if os.path.exists(src) and not os.path.lexists(lnk):
                os.makedirs(os.path.dirname(lnk) or workdir, exist_ok=True)
                os.symlink(src, lnk)
    outdir = os.path.join(workdir, f"output_{name}")
    text = set_output_dir_text(text, outdir)
    for f in ("intEC.dat", "intEC.extrap"):
        fp = os.path.join(outdir, f)
        if os.path.exists(fp):
            os.remove(fp)
    path = os.path.join(workdir, f"{name}.azr")
    with open(path, "w") as fh:
        fh.write(text)
    return path


def base_segments(base_azr: str) -> SegmentSet:
    return SegmentSet.from_file(base_azr)
