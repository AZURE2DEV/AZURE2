"""``rmfit review`` -- the post-fit "GUI eye" as checks (DESIGN.md Phase 1.3, 2026-11-07).

What the evaluator catches by looking at the AZURE2 GUI, as detectors on the per-point table of a finished fit:

- per-point summary in the archive's wording (fraction within 1 and 2 sigma, median |r|, how many points carry
  half of the chi-squared and where they sit) and the top chi-squared carriers;
- normalizations that moved AWAY from the data/fit mismatch of their own segment (e.g. a norm of 1.32 at
  median fit/data 0.64), fixed norms the data would move, and an optimal scale that trends with angle or energy
  inside one experiment (a shape problem dressed as a scale);
- residual runs: consecutive points on one side of the fit (bumpiness, a missing level, a wrong width);
- average-versus-structure: residuals that vanish when blurred over 50 keV while the raw ones are large (the fit
  follows the mean of structure the data resolve), and MeV-wide levels spanning several visible unfitted peaks;
- levels narrower than the data spacing that covers them (eV-wide levels living on the
  target-integration grid), levels supported by one data set or by none;
- optionally, on a campaign, the objective at a converged target-integration grid (``--grid-check``).

Two sources for the table: ``table_from_output(azr)`` reads ``output/AZUREOut_*.out`` + ``chiSquared.out`` +
``parameters.out`` (engine-free, any archived fit; the data columns there are n*d and n*e, the fit column bare --
azure2-eval skill, "Output files"); ``table_from_evaluator(ev, x)`` uses a single session (a campaign's incumbent).
Findings are candidates for the evaluator, never decisions (CLAUDE.md "Judging fit quality").
"""
import math
import os
import re
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .spec import Finding
from .backend import PYAZR3

# observables measured at an angle (pyazr.datasets._OBSERVABLE names).  "differential-cm" (isDiff 4) was missing until
# 2026-11-15: every c.m. differential set lost its angles, so the angle checks (norm-angle-trend, detector curves) skipped
# them.
ANGLE_OBS = ("differential", "differential-cm", "analyzing-power", "phase-shift")
# fit/data ratios mean nothing there (the observable names of pyazr 3 / pyazr 2)
RATIO_OBS_SKIP = ("analyzing_power", "phase_shift") if PYAZR3 else ("analyzing-power", "phase-shift")


def _at_angle(seg) -> bool:
    """True for a data set measured at an angle: pyazr 3's ``differential`` flag of the data-set table, pyazr 2's
    observable name in :data:`ANGLE_OBS`."""
    return bool(seg.differential) if PYAZR3 else seg.observable in ANGLE_OBS
BLUR_SCALES = (0.02, 0.05, 0.1, 0.2)                       # MeV in Ex
UNIT_MEV = {"meV": 1e-9, "eV": 1e-6, "keV": 1e-3, "MeV": 1.0}


# -- data structures ------------------------------------------------------------------------------------------

@dataclass
class Level:
    jpi: str
    energy: float                      # Ex, MeV
    gamma_mev: Optional[float] = None  # total width (sum of |partial| over particle + gamma channels), MeV
    channels: list = field(default_factory=list)   # (pair, L, S, width_mev or None)

    @property
    def label(self):
        g = "" if self.gamma_mev is None else f" G={_fmt_width(self.gamma_mev)}"
        return f"{self.jpi} {self.energy:.4f}{g}"


def _fmt_width(g_mev):
    if g_mev is None or not np.isfinite(g_mev):
        return "?"
    if g_mev >= 1.0:
        return f"{g_mev:.3g} MeV"
    if g_mev >= 1e-3:
        return f"{1e3 * g_mev:.3g} keV"
    if g_mev >= 1e-6:
        return f"{1e6 * g_mev:.3g} eV"
    return f"{1e9 * g_mev:.3g} meV"


@dataclass
class PointTable:
    key: np.ndarray        # segment key per point
    name: np.ndarray       # data set name (file stem) per point
    ecm: np.ndarray
    ex: np.ndarray
    ang: np.ndarray        # c.m. angle (NaN if unknown)
    fit: np.ndarray        # bare calculation c
    data: np.ndarray       # RAW data d (normalization removed)
    err: np.ndarray        # RAW error e
    norm: np.ndarray       # fitted normalization n per point
    notes: list = field(default_factory=list)

    @property
    def n(self):
        return int(self.key.size)

    @property
    def resid(self):
        """AZURE2's residual (c - n d) / (n e); chi2 = resid**2 (plain, correlations ignored)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (self.fit - self.norm * self.data) / (self.norm * self.err)
        return np.where(np.isfinite(r), r, 0.0)

    @property
    def chi2(self):
        return self.resid ** 2

    def keys(self):
        return [int(k) for k in np.unique(self.key)]

    def curves(self, min_points: int = 1) -> List[Tuple[int, float, np.ndarray]]:
        """(cached per table) ``(key, angle_label, indices sorted by Ex)`` per detector of each segment (one fixed lab detector drifts
        several degrees in c.m. over a long run -- cluster consecutive angles closer than 2 deg, as
        an archived campaign does)."""
        if getattr(self, "_curves", None) is None:
            out = []
            for k in self.keys():
                idx = np.flatnonzero(self.key == k)
                labs = detector_labels(self.ang[idx])
                for lab in np.unique(labs):
                    s = idx[labs == lab]
                    s = s[np.argsort(self.ex[s], kind="stable")]
                    out.append((k, float(lab), s))
            self._curves = out
        return [c for c in self._curves if c[2].size >= min_points]


def detector_labels(ang, gap: float = 2.0) -> np.ndarray:
    a = np.nan_to_num(np.asarray(ang, float), nan=-1.0)
    if a.size == 0:
        return a
    order = np.argsort(a)
    s = a[order]
    cid = np.concatenate([[0], np.cumsum(np.diff(s) > gap)])
    lab = np.empty_like(a)
    for c in np.unique(cid):
        lab[order[cid == c]] = round(float(s[cid == c].mean()), 1)
    return lab


def dataset_name(path: str) -> str:
    return os.path.splitext(os.path.basename(str(path)))[0]


def family_name(stem: str) -> str:
    """``meyer_84`` / ``meyer_114`` -> ``meyer``; ``SETB25_a1_30deg`` -> ``SETB25_a1`` (trailing numeric tokens
    with an optional unit are dropped)."""
    return re.sub(r"([_-]\d+(\.\d+)?(deg|MeV|keV|mev|kev)?)+$", "", stem) or stem


# -- loaders --------------------------------------------------------------------------------------------------

def read_chisq_out(path: str) -> Dict[int, dict]:
    """``{segment key: dict(chi2, n, norm, norm_chi2)}`` from ``chiSquared.out``."""
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        t = [s.strip() for s in line.split(",")]
        if len(t) >= 4 and t[0].isdigit():
            try:
                out[int(t[0])] = dict(chi2=float(t[1]), n=int(float(t[2])), norm=float(t[3]),
                                      norm_chi2=float(t[4]) if len(t) > 4 else 0.0)
            except ValueError:
                pass
    return out


def _blocks(path: str) -> List[np.ndarray]:
    blocks, cur = [], []
    for line in open(path):
        s = line.strip()
        if not s:
            if cur:
                blocks.append(np.array(cur, float)); cur = []
            continue
        try:
            cur.append([float(v) for v in s.split()])
        except ValueError:
            continue
    if cur:
        blocks.append(np.array(cur, float))
    return blocks


def output_dir_of(azr_path: str) -> str:
    """The output directory the model names (pyazr 2: ``<fileDir>`` of the .azr; pyazr 3: ``output_dir`` of
    ``<config>``; possibly relative to the model) when it holds a run's files, else ``output/`` next to it."""
    base = os.path.dirname(os.path.abspath(azr_path))
    if PYAZR3:
        try:
            import pyazr
            p3 = str(pyazr.Project.open(os.path.abspath(azr_path)).config().get("output_dir", "") or "")
        except Exception:
            p3 = ""
        m = None
    else:
        text = open(azr_path).read()
        m = re.search(r"<fileDir>\s*\n?(.*?)\n?\s*</fileDir>", text, re.S)
    if m or (PYAZR3 and p3):
        p = p3 if PYAZR3 else (m.group(1).strip().splitlines()[0].strip() if m.group(1).strip() else "")
        if p:
            p = os.path.expanduser(p)
            p = p if os.path.isabs(p) else os.path.join(base, p)
            if os.path.isdir(p) and any(f.startswith("AZUREOut_") for f in os.listdir(p)):
                return p
    return os.path.join(base, "output")


def table_from_output(azr_path: str, output_dir: Optional[str] = None) -> PointTable:
    """The point table of a finished CLI run (mode 1 or 2): ``AZUREOut_aa=<in>_R=<out>.out`` blocks in
    ``<segmentsData>`` order among the active segments of each pair; norms from ``chiSquared.out`` (the file's
    data columns are n*d and n*e; they are divided out here).  Raises if no output file is found."""
    from .model import base_segments, rows_in_segment
    azr_path = os.path.abspath(azr_path)
    base = os.path.dirname(azr_path)
    odir = output_dir or output_dir_of(azr_path)
    segs = base_segments(azr_path)
    chisq = read_chisq_out(os.path.join(odir, "chiSquared.out"))
    cols = defaultdict(list)
    notes = []
    groups: Dict[Tuple[int, int], list] = defaultdict(list)
    for s in segs:
        if s.active:
            groups[(int(s.entrance_key), int(s.exit_key))].append(s)
    nfiles = 0
    for (ain, rout), group in groups.items():
        tail = "TOTAL_CAPTURE" if rout < 0 else f"R={rout}"
        path = os.path.join(odir, f"AZUREOut_aa={ain}_{tail}.out")
        if not os.path.exists(path):
            notes.append(f"missing {os.path.basename(path)} for segments {[s.key for s in group]}")
            continue
        nfiles += 1
        blocks = _blocks(path)
        expected = []
        for s in group:
            n = chisq.get(s.key, {}).get("n")
            if n is None:
                rows = rows_in_segment(s, base)
                n = int(s.rows_kept if rows is None else rows.shape[0])      # None: pyazr 3, rows not placed
            expected.append(n)
        if len(blocks) != len(group) or any(b.shape[0] != n for b, n in zip(blocks, expected)):
            allrows = np.concatenate(blocks) if blocks else np.zeros((0, 9))
            if allrows.shape[0] == sum(expected):
                notes.append(f"{os.path.basename(path)}: {len(blocks)} blocks for {len(group)} segments; split by the "
                             f"chiSquared.out / data-file point counts instead")
                blocks, at = [], 0
                for n in expected:
                    blocks.append(allrows[at:at + n]); at += n
            else:
                notes.append(f"{os.path.basename(path)}: {allrows.shape[0]} rows vs {sum(expected)} expected for "
                             f"segments {[s.key for s in group]}; blocks assigned in order, counts may be off")
                blocks = (blocks + [np.zeros((0, 9))] * len(group))[:len(group)]
        for s, b in zip(group, blocks):
            if b.size == 0:
                continue
            n = float(chisq.get(s.key, {}).get("norm", s.norm) or 1.0)
            m = b.shape[0]
            cols["key"].append(np.full(m, s.key, int))
            cols["name"].append(np.array([dataset_name(s.data_file)] * m, dtype=object))
            cols["ecm"].append(b[:, 0]); cols["ex"].append(b[:, 1])
            cols["ang"].append(b[:, 2] if _at_angle(s) else np.full(m, np.nan))
            cols["fit"].append(b[:, 3])
            cols["data"].append(b[:, 5] / n); cols["err"].append(np.abs(b[:, 6]) / n)
            cols["norm"].append(np.full(m, n))
    if nfiles == 0:
        raise FileNotFoundError(f"no AZUREOut_*.out for the active segments of {azr_path} in {odir}")
    if not cols:
        raise FileNotFoundError(f"AZUREOut_*.out files in {odir} hold no rows for the active segments")
    tab = PointTable(**{k: np.concatenate(v) for k, v in cols.items()}, notes=notes)
    tab.notes.insert(0, f"source: {odir} ({nfiles} AZUREOut files, {tab.n} points, "
                        f"{'norms from chiSquared.out' if chisq else 'norms from the .azr (no chiSquared.out)'})")
    return tab


def table_from_evaluator(ev, x) -> PointTable:
    """The point table of a single session at ``x`` (``ev.data()`` raw data, ``ev.cross(x)`` the bare calculation,
    norms from the vector or the segment)."""
    d = ev.data()
    c = ev.cross(np.asarray(x, float))
    x = np.asarray(x, float)
    norm_of = {k.seg: float(x[j]) for j, k in enumerate(ev.keys) if k.kind == "norm"}
    cols = defaultdict(list)
    angles = d.get("angles") or []
    for i, (key, (a, b)) in enumerate(ev.row_ranges):
        if i >= len(d["cross"]):
            break
        seg = ev.segments.by_key(key)
        ex = np.asarray(d["excitation"][i], float); dat = np.asarray(d["cross"][i], float)
        err = np.abs(np.asarray(d["cross_err"][i], float)); fit = np.asarray(c[i], float)
        m = ex.size
        n = norm_of.get(key, float(getattr(seg, "norm", 1.0) or 1.0))
        cols["key"].append(np.full(m, key, int))
        cols["name"].append(np.array([dataset_name(seg.data_file) if seg is not None else str(key)] * m, dtype=object))
        cols["ecm"].append(np.asarray(d["energies"][i], float) if i < len(d.get("energies", [])) else np.full(m, np.nan))
        cols["ex"].append(ex)
        ang = np.asarray(angles[i], float) if i < len(angles) and len(angles[i]) == m else np.full(m, np.nan)
        if seg is not None and not _at_angle(seg):
            ang = np.full(m, np.nan)
        cols["ang"].append(ang)
        cols["fit"].append(fit); cols["data"].append(dat); cols["err"].append(err); cols["norm"].append(np.full(m, n))
    tab = PointTable(**{k: np.concatenate(v) for k, v in cols.items()})
    tab.notes.append(f"source: single session on {getattr(ev, 'azr', '?')} ({tab.n} points)")
    return tab


_LEVEL_RE = re.compile(r"J\s*=\s*([\d.]+)([+-])\s+E_level\s*=\s*([-\d.]+)\s*MeV")
_CHAN_RE = re.compile(r"R\s*=\s*(\d+)\s+[lL]\s*=\s*(\S+)\s+s\s*=\s*([\d.]+)\s+(G|C)\s*=\s*([-\d.eE+]+)\s*(\S+)")


def levels_from_parameters_out(path: str) -> List[Level]:
    """Levels with physical widths from ``output/parameters.out`` (one ``J = .. E_level = ..`` block per level;
    ``G = <value> <meV|eV|keV|MeV>`` per channel, ``C = <ANC>`` for a bound channel)."""
    levels: List[Level] = []
    for line in open(path):
        m = _LEVEL_RE.search(line)
        if m:
            j = float(m.group(1)); jpi = (f"{int(2*j)}/2" if abs(2 * j - round(2 * j)) < 1e-6 and int(round(2 * j)) % 2 else f"{int(round(j))}") + m.group(2)
            levels.append(Level(jpi=jpi, energy=float(m.group(3)), gamma_mev=0.0))
            continue
        m = _CHAN_RE.search(line)
        if m and levels:
            L = m.group(2); kind = m.group(4); val = float(m.group(5)); unit = m.group(6)
            w = abs(val) * UNIT_MEV[unit] if kind == "G" and unit in UNIT_MEV else None
            try:
                Lval = int(L)
            except ValueError:
                Lval = L
            levels[-1].channels.append((int(m.group(1)), Lval, float(m.group(3)), w))
            if w is not None:
                levels[-1].gamma_mev += w
    return levels


def levels_from_azr(azr_path: str) -> List[Level]:
    """Levels of the .azr with total widths from its ``gamma`` fields, which hold PHYSICAL values (eV for an open
    particle channel or a photon channel, an ANC for a closed one -- ``AzrModel.apply_fit`` docstring).  Open means
    the level lies above the pair's threshold Ex (``entrance_sep_energy - sep_energy + excitation``).  These are
    the values the file was last saved with, not necessarily the fit's: prefer ``parameters.out`` when it exists."""
    from .backend import AzrModel
    from .spec import pair_lines
    m = AzrModel.from_file(azr_path)
    thr = {}
    for key, c in pair_lines(m).items():
        try:
            if getattr(c, "is_photon", False) or (not PYAZR3 and c.ptype == 10):
                thr[key] = None            # photon: always "open"
            elif PYAZR3:
                thr[key] = float(c.threshold)              # Project.pairs(): separation energy + excitation
            else:
                thr[key] = float(c.entrance_sep_energy - c.sep_energy + c.excitation)
        except Exception:
            thr[key] = float("nan")
    if PYAZR3:
        levels = m.active_levels
    else:
        try:
            levels = m.active_levels()
        except Exception:
            levels = m.levels
    out = []
    for lv in levels:
        j = float(lv.J); pi = "+" if int(lv.parity) > 0 else "-"
        jpi = (f"{int(2*j)}/2" if int(round(2 * j)) % 2 else f"{int(round(j))}") + pi
        e = float(lv.energy); g = 0.0; chans = []
        for c in lv.channels:
            t = thr.get(c.pair, float("nan"))
            is_open = (t is None) or (np.isfinite(t) and e > t)
            w = abs(float(c.gamma)) * 1e-6 if is_open else None
            chans.append((int(c.pair), int(c.L), float(c.S), w))
            if w is not None:
                g += w
        out.append(Level(jpi=jpi, energy=e, gamma_mev=g, channels=chans))
    return out


def levels_from_evaluator(ev, x) -> List[Level]:
    """Levels of a single session at ``x``: energies from the vector / the level table, total widths from the
    free width keys (``local.level_total_widths``, eV)."""
    from .local import level_total_widths
    from .keys import level_group
    x = np.asarray(x, float)
    en = {lab: float(e) for lab, e in getattr(ev, "level_energy", {}).items()}
    for j, k in enumerate(ev.keys):
        if k.kind == "energy" and k.level in en:
            en[k.level] = float(x[j])
    gam = {}
    try:
        phys = np.asarray(ev.transform(x), float)
        try:
            gw = ev.wigner_widths(x)           # Gamma_W > 0 only for an open particle channel
        except Exception:
            gw = None
        for j, k in enumerate(ev.keys):
            if k.kind != "width" or not k.level or j >= phys.size or not np.isfinite(phys[j]):
                continue
            if gw is not None and gw.get(j, 0.0) <= 0:
                continue                       # closed channel: the value is an ANC (fm^-1/2), not a width
            if abs(phys[j]) > 1e8:             # 100 MeV: not a width whatever the engine says (seen for a level 3 keV below a threshold)
                continue
            gam[k.level] = gam.get(k.level, 0.0) + abs(float(phys[j]))
    except Exception:
        try:
            gam = level_total_widths(ev, x)
        except Exception:
            gam = {}
    out = []
    for lab, e in en.items():
        g = gam.get(lab)
        out.append(Level(jpi=level_group(lab).split("@")[0], energy=e, gamma_mev=None if g is None else float(g) * 1e-6))
    out.sort(key=lambda l: l.energy)
    return out


# -- helpers ----------------------------------------------------------------------------------------------------

def blur(e, y, s):
    """Gaussian kernel regression of ``y(e)`` at the points ``e`` (sorted), window +-4 sigma."""
    lo = np.searchsorted(e, e - 4 * s); hi = np.searchsorted(e, e + 4 * s, side="right")
    out = np.empty_like(y, dtype=float)
    for i in range(len(e)):
        w = np.exp(-0.5 * ((e[lo[i]:hi[i]] - e[i]) / s) ** 2)
        out[i] = w @ y[lo[i]:hi[i]] / w.sum()
    return out


def spearman(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return 0.0
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = math.sqrt(float(ra @ ra) * float(rb @ rb))
    return float(ra @ rb / den) if den > 0 else 0.0


def optimal_norm(c, d, e) -> Optional[float]:
    """The chi2-optimal normalization of one segment at fixed R-matrix parameters, AZURE2 convention
    (skill: u* = sum(c d / e^2) / sum(c^2 / e^2), n* = 1 / u*)."""
    ok = np.isfinite(c) & np.isfinite(d) & (e > 0)
    if ok.sum() < 2:
        return None
    den = float(np.sum(c[ok] ** 2 / e[ok] ** 2))
    if den <= 0:
        return None
    u = float(np.sum(c[ok] * d[ok] / e[ok] ** 2)) / den
    return 1.0 / u if u > 0 else None


def peaks(y, err, w: int = 3, nsig: float = 3.0) -> np.ndarray:
    """Indices of local maxima of ``y`` standing more than ``nsig`` errors above the lower of the two
    neighbourhood minima within ``w`` points (a plain prominence test; resonance-search/peakfind.py idea)."""
    y = np.asarray(y, float); err = np.asarray(err, float)
    out = []
    for i in range(1, y.size - 1):
        if not (y[i] > y[i - 1] and y[i] >= y[i + 1]):
            continue
        left = y[max(0, i - w):i].min(); right = y[i + 1:i + 1 + w].min()
        prom = y[i] - max(left, right)
        if prom > nsig * max(err[i], 1e-300):
            out.append(i)
    return np.array(out, int)


def _share_count(chi2_sorted_desc, share):
    tot = chi2_sorted_desc.sum()
    if tot <= 0:
        return 0
    cs = np.cumsum(chi2_sorted_desc)
    return int(np.searchsorted(cs, share * tot) + 1)


# -- the review ------------------------------------------------------------------------------------------------

# Width-tied closeness of two same-J^pi levels: a gap below this fraction of the NARROWER total width makes the pair a
# merge candidate for the evaluator (report only; nothing in rmfit enforces a spacing).  Origin (2026-11-08):
# a same-J^pi pair 1 keV apart with widths ~220 keV; the 10 % is the scale proposed then and is the evaluator's to
# change (rmfit review --close-fraction, spec budget "close_fraction"; 0 switches the check off).  Known-level check of
# the criterion: it flags that pair, and the evaluator's merge tests on the main line found the pair NOT redundant (the
# (a,n2g) data need both, an archived campaign) -- so the flag is a question, never a verdict.
CLOSE_GAP_FRACTION = 0.10


class Review:
    """Run the detectors on a point table; ``findings`` (``spec.Finding``) and markdown ``sections``."""

    def __init__(self, tab: PointTable, levels: Sequence[Level], segments=None, azr_text: str = "",
                 top: int = 20, close_fraction: Optional[float] = None):
        self.close_fraction = CLOSE_GAP_FRACTION if close_fraction is None else float(close_fraction)
        self.tab = tab
        self.levels = sorted(levels, key=lambda l: l.energy)
        self.segments = segments
        self.azr_text = azr_text
        self.top = top
        self.findings: List[Finding] = []
        self.sections: List[Tuple[str, List[str]]] = []
        ex = tab.ex[np.isfinite(tab.ex)]
        self.ex_lo = float(ex.min()) if ex.size else float("nan")
        self.ex_hi = float(ex.max()) if ex.size else float("nan")
        self.chi2_total = float(tab.chi2.sum())

    # ---- helpers
    def _seg(self, key):
        if self.segments is None:
            return None
        try:
            return self.segments.by_key(int(key))
        except Exception:
            for s in self.segments:
                if int(s.key) == int(key):
                    return s
        return None

    def _in_data(self, level: Level, margin: float = 0.5) -> bool:
        return (self.ex_lo - margin) <= level.energy <= (self.ex_hi + margin)

    def nearest_level(self, ex: float):
        if not self.levels:
            return None
        return min(self.levels, key=lambda l: abs(l.energy - ex))

    def add(self, severity, code, message, evidence=""):
        self.findings.append(Finding(severity, code, message, evidence))

    def run(self, figures_dir: Optional[str] = None, nfig: int = 6):
        self.per_point()
        self.top_points()
        self.norms()
        self.norm_trends()
        self.residual_runs()
        self.structure_scale()
        self.fit_zigzag()
        self.level_support()
        self.close_pairs()
        if figures_dir:
            try:
                self.figures(figures_dir, nfig)
            except Exception as err:          # a figure never stops the review
                self.sections.append(("Figures", [f"(figures failed: {err})"]))
        return self

    # ---- 1. per point ---------------------------------------------------------------------------------------
    def per_point(self):
        t = self.tab
        r = np.abs(t.resid); c2 = t.chi2
        rows = []
        for name in sorted(set(t.name), key=lambda nm: -float(c2[t.name == nm].sum()) / max(int((t.name == nm).sum()), 1)):
            sel = t.name == name
            rr = r[sel]; cc = c2[sel]
            order = np.argsort(-cc)
            k50 = _share_count(cc[order], 0.5); k80 = _share_count(cc[order], 0.8)
            idx = np.flatnonzero(sel)[order[:k50]]
            exs = t.ex[idx]; angs = t.ang[idx]
            where = f"Ex {exs.min():.3f}-{exs.max():.3f}" if exs.size else ""
            if np.isfinite(angs).any():
                a = angs[np.isfinite(angs)]
                where += f", {a.min():.0f}-{a.max():.0f} deg" if a.max() - a.min() > 1 else f", {a.mean():.0f} deg"
            rows.append(dict(name=name, n=int(sel.sum()), chi2=float(cc.sum()), per_pt=float(cc.sum()) / max(int(sel.sum()), 1),
                             f1=float((rr < 1).mean()), f2=float((rr < 2).mean()), med=float(np.median(rr)),
                             k50=k50, k80=k80, where=where, keys=sorted(set(int(k) for k in t.key[sel]))))
        self.per_point_rows = rows
        n = t.n
        L = [f"{n} points, chi2 {self.chi2_total:.1f} (plain, correlations ignored), chi2/N {self.chi2_total / max(n, 1):.2f}; "
             f"{100 * float((r < 1).mean()):.0f} % within 1 sigma, {100 * float((r < 2).mean()):.0f} % within 2 sigma, "
             f"median |r| {float(np.median(r)):.2f}; {_share_count(np.sort(c2)[::-1], 0.5)} points carry half of the chi2.", "",
             "| data set | N | chi2/N | chi2 | <1s | <2s | med |r| | pts for 50 % | pts for 80 % | where the 50 % sit | segments |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|"]
        for d in rows:
            segs = ",".join(str(k) for k in d["keys"][:6]) + ("..." if len(d["keys"]) > 6 else "")
            L.append(f"| {d['name']} | {d['n']} | {d['per_pt']:.2f} | {d['chi2']:.1f} | {100*d['f1']:.0f} % | {100*d['f2']:.0f} % | "
                     f"{d['med']:.2f} | {d['k50']} | {d['k80']} | {d['where']} | {segs} |")
        self.sections.append(("Per-point summary", L))
        for d in rows:
            if d["per_pt"] < 3 or d["n"] < 5:
                continue
            if d["k50"] <= max(3, 0.1 * d["n"]) and d["f2"] >= 0.5:
                self.add("LOW", "point-carried",
                         f"{d['name']}: most points described ({100*d['f2']:.0f} % within 2 sigma); {d['k50']} of {d['n']} points "
                         f"carry half of its chi2 ({d['where']})", f"chi2/N {d['per_pt']:.1f}; possible outliers or a sharp edge -- the evaluator's call")
            elif d["f2"] < 0.5:
                self.add("MEDIUM", "systematic-misfit",
                         f"{d['name']}: only {100*d['f2']:.0f} % of {d['n']} points within 2 sigma, median |r| {d['med']:.1f} -- a shape or "
                         f"scale problem across the set, not isolated points", f"chi2/N {d['per_pt']:.1f}; half of the chi2 in {d['k50']} points ({d['where']})")

    # ---- 2. top points
    def top_points(self):
        t = self.tab
        order = np.argsort(-t.chi2)[:self.top]
        share = float(t.chi2[order].sum()) / max(self.chi2_total, 1e-300)
        L = [f"the {len(order)} largest carry {100 * share:.1f} % of the chi2", "",
             "| chi2 | r | data set | seg | Ex | angle | data | fit | norm |", "|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
        for i in order:
            ang = "" if not np.isfinite(t.ang[i]) else f"{t.ang[i]:.1f}"
            L.append(f"| {t.chi2[i]:.1f} | {t.resid[i]:+.1f} | {t.name[i]} | {t.key[i]} | {t.ex[i]:.4f} | {ang} | "
                     f"{t.data[i]:.3e} +- {t.err[i]:.1e} | {t.fit[i]:.3e} | {t.norm[i]:.3f} |")
        self.sections.append((f"Top {len(order)} chi2 carriers", L))

    # ---- 3. normalizations --------------------------------------------------------------------------------------
    def norms(self):
        t = self.tab
        rows = []
        for k in t.keys():
            seg = self._seg(k)
            sel = t.key == k
            if seg is not None and seg.observable in RATIO_OBS_SKIP:
                continue
            c, d, e = t.fit[sel], t.data[sel], t.err[sel]
            n = float(np.median(t.norm[sel]))
            nstar = optimal_norm(c, d, e)
            pos = (d > 0) & np.isfinite(c)
            med = float(np.median(c[pos] / d[pos])) if pos.sum() >= 2 else float("nan")
            chi_n = float(t.chi2[sel].sum())
            chi_star = float(np.sum(((c - nstar * d) / (nstar * e)) ** 2)) if nstar else float("nan")
            vary = bool(getattr(seg, "vary_norm", False)) if seg is not None else None
            nerr = float(getattr(seg, "norm_error", 0.0) or 0.0) if seg is not None else 0.0
            if nerr > 1.0:
                nerr = nerr / 100.0        # AZURE2 reads percent; a legacy fraction > 1 would be absurd anyway
            rows.append(dict(key=k, name=str(t.name[sel][0]), n=int(sel.sum()), norm=n, nstar=nstar, med=med, chi_n=chi_n,
                             chi_star=chi_star, vary=vary, nerr=nerr, group=getattr(seg, "norm_group", None) if seg is not None else None,
                             cov=getattr(seg, "covariance", None) if seg is not None else None))
        self.norm_rows = rows
        L = ["n = fitted normalization (multiplies the data); n* = the value this segment alone would take at the fitted "
             "R-matrix parameters; median c/d = fit over raw data.  A free n on the other side of 1 from n* was pulled by "
             "something else: a norm group, a covariance file, the norm prior or another segment.", "",
             "| seg | data set | N | free | n | n* | median c/d | chi2 at n | chi2 at n* | sys % | group / cov |",
             "|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---|"]
        for d in rows:
            ns = "" if d["nstar"] is None else f"{d['nstar']:.3f}"
            cs = "" if not np.isfinite(d["chi_star"]) else f"{d['chi_star']:.1f}"
            gc = " ".join(s for s in (d["group"] or "", os.path.basename(d["cov"]) if d["cov"] else "") if s)
            L.append(f"| {d['key']} | {d['name']} | {d['n']} | {'' if d['vary'] is None else ('yes' if d['vary'] else 'fixed')} | "
                     f"{d['norm']:.3f} | {ns} | {d['med']:.3f} | {d['chi_n']:.1f} | {cs} | {100*d['nerr']:.1f} | {gc} |")
        self.sections.append(("Normalizations", L))
        fixed_wants: Dict[str, list] = {}
        for d in rows:
            ns = d["nstar"]
            if ns is None or not np.isfinite(ns):
                continue
            n = d["norm"]
            if d["vary"] and (n - 1) * (ns - 1) < 0 and abs(n - 1) > 0.05 and abs(ns - 1) > 0.05:
                why = ("norm group " + d["group"]) if d["group"] else ("covariance " + os.path.basename(d["cov"])) if d["cov"] else \
                    ("the norm prior" if d["nerr"] > 0 else "another segment sharing the parameter")
                self.add("MEDIUM", "norm-opposite",
                         f"{d['name']} (seg {d['key']}): fitted norm {n:.3f} but this segment alone wants n* = {ns:.3f} "
                         f"(median fit/data {d['med']:.2f}) -- the norm moved AWAY from the mismatch; likely pulled by {why}",
                         f"chi2 {d['chi_n']:.1f} at n vs {d['chi_star']:.1f} at n*")
            elif d["vary"] and abs(n / ns - 1) > 0.2 and d["chi_n"] - d["chi_star"] > 0.1 * max(d["chi_star"], 1.0):
                self.add("LOW", "norm-off-optimum",
                         f"{d['name']} (seg {d['key']}): free norm {n:.3f} vs its own optimum {ns:.3f}",
                         f"chi2 {d['chi_n']:.1f} -> {d['chi_star']:.1f} at n*; held by " + (("norm group " + d["group"]) if d["group"] else "its prior / shared rows"))
            elif d["vary"] is False and abs(ns - 1) > max(2 * d["nerr"], 0.10) and d["chi_n"] - d["chi_star"] > 0.3 * d["chi_n"]:
                fixed_wants.setdefault(d["name"], []).append(d)
        for name, ds in fixed_wants.items():
            ns = np.array([d["nstar"] for d in ds]); gain = sum(d["chi_n"] - d["chi_star"] for d in ds); tot = sum(d["chi_n"] for d in ds)
            nerr = max(d["nerr"] for d in ds)
            if ns.max() < 0.05:
                self.add("INFO", "fit-near-zero", f"{name}: the calculation is essentially zero where these data are not "
                         f"(fixed norm {ds[0]['norm']:.3f}; n* < 0.05 on {len(ds)} segment(s))",
                         f"chi2 {tot:.1f}; a channel the model does not feed (closed, or no level with that partial width) -- segments "
                         f"{','.join(str(d['key']) for d in ds[:8])}")
                continue
            self.add("INFO", "fixed-norm-wants",
                     f"{name}: fixed norm {ds[0]['norm']:.3f} on {len(ds)} segment(s); the data alone would take n* "
                     f"{ns.min():.3f}-{ns.max():.3f} (median {np.median(ns):.3f}" +
                     (f", beyond 2x the {100*nerr:.0f} % systematic)" if nerr > 0 else ", no systematic recorded)"),
                     f"chi2 {tot:.1f} -> {tot - gain:.1f} at n*; a scale offset the model absorbs elsewhere, or a shape problem if n* does not help "
                     f"(segments {','.join(str(d['key']) for d in ds[:8])}{'...' if len(ds) > 8 else ''})")

    # ---- 4. norm trends inside a family --------------------------------------------------------------------------
    def norm_trends(self):
        t = self.tab
        rows = getattr(self, "norm_rows", None) or []
        L = []
        by = {d["key"]: d for d in rows}
        fam: Dict[Tuple, list] = defaultdict(list)
        for k in t.keys():
            seg = self._seg(k)
            if k not in by or by[k]["nstar"] is None:
                continue
            sel = t.key == k
            ang = t.ang[sel]; ex = t.ex[sel]
            fam_key = (family_name(by[k]["name"]), int(seg.entrance_key) if seg is not None else 0, int(seg.exit_key) if seg is not None else 0)
            fam[fam_key].append(dict(key=k, nstar=by[k]["nstar"], ang=float(np.nanmean(ang)) if np.isfinite(ang).any() else float("nan"),
                                     ex=float(np.mean(ex)), name=by[k]["name"]))
        for (fname, ain, rout), members in fam.items():
            if len(members) < 3:
                continue
            ns = np.array([m["nstar"] for m in members]); rng = float(ns.max() / ns.min()) if ns.min() > 0 else float("inf")
            for axis, label in (("ang", "angle"), ("ex", "energy")):
                v = np.array([m[axis] for m in members])
                ok = np.isfinite(v)
                if ok.sum() < 3 or np.unique(np.round(v[ok], 1)).size < 3:
                    continue
                rho = spearman(v[ok], ns[ok])
                L.append(f"- {fname} (pair {ain}->{rout}, {ok.sum()} segments): n* {ns[ok].min():.3f}-{ns[ok].max():.3f} vs {label}, Spearman {rho:+.2f}")
                if abs(rho) >= 0.8 and rng > 1.10:
                    lo = members[int(np.argmin(np.where(ok, v, np.inf)))]; hi = members[int(np.argmax(np.where(ok, v, -np.inf)))]
                    self.add("MEDIUM", f"norm-{label}-trend",
                             f"{fname}: the optimal scale of its {ok.sum()} segments runs with {label} (Spearman {rho:+.2f}, "
                             f"n* {ns[ok].min():.3f}-{ns[ok].max():.3f}; {label} {min(lo[axis], hi[axis]):.1f}-{max(lo[axis], hi[axis]):.1f}) -- "
                             f"an {'angular-distribution' if label == 'angle' else 'energy-dependent'} shape problem dressed as normalization",
                             "one experiment has one scale; segments " + ",".join(str(m["key"]) for m in members))
        # inside one segment: optimal scale per detector
        for k in t.keys():
            seg = self._seg(k)
            if seg is not None and seg.observable in RATIO_OBS_SKIP:
                continue
            cur = [(lab, s) for kk, lab, s in t.curves(min_points=4) if kk == k and np.isfinite(t.ang[s]).all()]
            if len(cur) < 3:
                continue
            labs, ns = [], []
            for lab, s in cur:
                n_ = optimal_norm(t.fit[s], t.data[s], t.err[s])
                if n_ is not None and np.isfinite(n_) and n_ > 0:
                    labs.append(lab); ns.append(n_)
            if len(ns) < 3:
                continue
            labs = np.array(labs); ns = np.array(ns)
            rho = spearman(labs, ns); rng = float(ns.max() / ns.min())
            name = str(t.name[t.key == k][0])
            L.append(f"- {name} (seg {k}, {len(ns)} detectors): per-angle n* {ns.min():.3f}-{ns.max():.3f}, Spearman vs angle {rho:+.2f}")
            if abs(rho) >= 0.8 and rng > 1.15:
                self.add("MEDIUM", "norm-angle-trend",
                         f"{name} (seg {k}): the scale its detectors want runs from {ns.min():.3f} at {labs[np.argmin(ns)]:.0f} deg to "
                         f"{ns.max():.3f} at {labs[np.argmax(ns)]:.0f} deg (Spearman {rho:+.2f}) -- the angular distribution's shape is "
                         f"off, a single normalization cannot fix it", f"{len(ns)} detectors; one experiment has one scale")
        self.sections.append(("Normalization trends inside one experiment", L or ["- no family with three or more segments / detectors"]))

    # ---- 5. residual runs ---------------------------------------------------------------------------------------
    def residual_runs(self, z: float = 2.0, min_run: int = 4):
        t = self.tab
        r = t.resid; c2 = t.chi2
        runs = []
        for k, lab, s in t.curves(min_points=min_run):
            rr = r[s]; ex = t.ex[s]
            i = 0
            while i < s.size:
                if abs(rr[i]) > z:
                    sgn = 1 if rr[i] > 0 else -1; j = i
                    while j + 1 < s.size and abs(rr[j + 1]) > z and (rr[j + 1] > 0) == (sgn > 0):
                        j += 1
                    if j - i + 1 >= min_run:
                        idx = s[i:j + 1]
                        runs.append(dict(key=k, ang=lab, name=str(t.name[s[0]]), lo=float(ex[i]), hi=float(ex[j]), n=j - i + 1,
                                         sign=sgn, maxz=float(np.max(np.abs(rr[i:j + 1]))), chi2=float(c2[idx].sum())))
                    i = j + 1
                else:
                    i += 1
        runs.sort(key=lambda d: -d["chi2"])
        self.runs = runs
        L = [f"runs of >= {min_run} consecutive points (in Ex, per detector) more than {z:g} sigma on one side of the fit: "
             f"{len(runs)} found" + (f", together {100 * sum(d['chi2'] for d in runs) / max(self.chi2_total, 1e-300):.0f} % of the chi2" if runs else ""), "",
             "| chi2 | data set | seg | angle | Ex from | Ex to | pts | side | max |r| | nearest level |", "|---:|---|---:|---:|---:|---:|---:|---|---:|---|"]
        for d in runs[:self.top]:
            mid = 0.5 * (d["lo"] + d["hi"]); lv = self.nearest_level(mid)
            inside = lv is not None and d["lo"] - 0.5 * (d["hi"] - d["lo"]) <= lv.energy <= d["hi"] + 0.5 * (d["hi"] - d["lo"])
            d["level"] = lv.label if lv else ""; d["level_inside"] = inside
            ang = "" if d["ang"] < 0 else f"{d['ang']:.0f}"
            L.append(f"| {d['chi2']:.1f} | {d['name']} | {d['key']} | {ang} | {d['lo']:.4f} | {d['hi']:.4f} | {d['n']} | "
                     f"{'fit above' if d['sign'] > 0 else 'fit below'} | {d['maxz']:.1f} | {d['level'] + (' (inside)' if inside else ' (nearest)') if lv else ''} |")
        self.sections.append(("Residual runs", L))
        nflag = 0
        seen = set()
        for d in runs:
            if nflag >= 5:
                break
            if d["chi2"] >= 0.03 * self.chi2_total or d["n"] >= 8:
                lv = self.nearest_level(0.5 * (d["lo"] + d["hi"]))
                tag = (d["name"], lv.label if (lv is not None and d.get("level_inside")) else None)
                if tag in seen:
                    continue
                seen.add(tag)
                nflag += 1
                if d.get("level_inside"):
                    hint = f"level {lv.label} sits inside: its width, shape or interference is off"
                else:
                    hint = ("no level inside the run" + (f" (nearest {lv.label})" if lv else "") +
                            " -- a missing-level candidate or a background / data problem")
                ang = "" if d["ang"] < 0 else f", {d['ang']:.0f} deg"   # not nested: Python < 3.12
                self.add("MEDIUM", "residual-run",
                         f"{d['name']} (seg {d['key']}{ang}): {d['n']} consecutive points "
                         f"{'above' if d['sign'] < 0 else 'below'} the fit by more than {z:g} sigma at Ex {d['lo']:.3f}-{d['hi']:.3f} MeV; {hint}",
                         f"chi2 {d['chi2']:.1f} ({100 * d['chi2'] / max(self.chi2_total, 1e-300):.1f} % of the total), max |r| {d['maxz']:.1f}")

    # ---- 6. average versus structure -----------------------------------------------------------------------------
    def structure_scale(self):
        t = self.tab
        r = t.resid
        acc: Dict[str, Dict[float, list]] = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0, 0]))   # name -> D -> [sum r0^2, sum rD^2, n]
        data_only_peaks = []      # (ex, name, key, ang)
        for k, lab, s in t.curves(min_points=12):
            ex = t.ex[s]
            if np.unique(np.round(ex, 6)).size < 12:
                continue
            span = float(ex[-1] - ex[0])
            name = str(t.name[s[0]])
            rr = r[s]
            for D in BLUR_SCALES:
                if span < 6 * D:
                    continue
                m = (ex >= ex[0] + 1.5 * D) & (ex <= ex[-1] - 1.5 * D)
                if m.sum() < 3:
                    continue
                rD = blur(ex, rr, D)
                a = acc[name][D]
                a[0] += float(np.sum(rr[m] ** 2)); a[1] += float(np.sum(rD[m] ** 2)); a[2] += int(m.sum())
            # peaks the data show and the fit does not
            nd = t.norm[s] * t.data[s]; ne = t.norm[s] * t.err[s]
            pd = peaks(nd, ne); pf = peaks(t.fit[s], ne)
            if pd.size:
                spacing = float(np.median(np.diff(ex))) if ex.size > 1 else 0.0
                for i in pd:
                    tol = max(3 * spacing, 0.01)
                    if not pf.size or np.min(np.abs(ex[pf] - ex[i])) > tol:
                        # the fit must also miss the peak height: data above fit by > 2 sigma at the peak
                        if (nd[i] - t.fit[s][i]) > 2 * ne[i]:
                            data_only_peaks.append((float(ex[i]), name, k, lab))
        L = ["rms residual per data set, raw (r0) and after a Gaussian blur of the residuals over D in Ex (rD); rD/r0 well below 1 "
             "with a large r0 means the fit passes through the AVERAGE of structure the data resolve at scales below D", "",
             "| data set | " + " | ".join(f"D={int(1000*D)} keV: r0 / rD" for D in BLUR_SCALES) + " |", "|---|" + "---|" * len(BLUR_SCALES)]
        flagged = []
        for name in sorted(acc, key=lambda nm: -sum(v[2] for v in acc[nm].values())):
            cells = []
            for D in BLUR_SCALES:
                a = acc[name].get(D)
                if not a or a[2] == 0:
                    cells.append("-"); continue
                r0 = math.sqrt(a[0] / a[2]); rD = math.sqrt(a[1] / a[2])
                cells.append(f"{r0:.2f} / {rD:.2f}")
                if D == 0.05 and r0 > 2.0 and rD / max(r0, 1e-9) < 0.4:
                    flagged.append((name, r0, rD, a[2]))
            L.append(f"| {name} | " + " | ".join(cells) + " |")
        # cluster the data-only peaks across detectors (10 keV)
        clusters = []
        for ex, name, k, lab in sorted(data_only_peaks):
            if clusters and ex - clusters[-1]["ex_hi"] <= 0.010:
                c = clusters[-1]; c["ex_hi"] = ex; c["members"].append((name, k, lab)); c["ex"] = 0.5 * (c["ex_lo"] + ex)
            else:
                clusters.append(dict(ex=ex, ex_lo=ex, ex_hi=ex, members=[(name, k, lab)]))
        strong = [c for c in clusters if len(c["members"]) >= 2 or len(set(m[0] for m in c["members"])) >= 1 and len(c["members"]) >= 2]
        self.unfit_peaks = clusters
        if clusters:
            L += ["", f"peaks in the data (prominence > 3 sigma) with no peak in the fit nearby and the fit > 2 sigma low at the top: "
                      f"{len(clusters)} distinct energies" + (f"; seen in >= 2 detectors: {len(strong)}" if strong else ""),
                  "", "| Ex | detectors | data sets | nearest level |", "|---:|---:|---|---|"]
            for c in sorted(clusters, key=lambda c: -len(c["members"]))[:self.top]:
                lv = self.nearest_level(c["ex"])
                away = "" if lv is None else f" ({1e3 * abs(lv.energy - c['ex']):.0f} keV away)"   # not nested: Python < 3.12
                L.append(f"| {c['ex']:.4f} | {len(c['members'])} | {', '.join(sorted(set(m[0] for m in c['members'])))[:60]} | "
                         f"{lv.label if lv else ''}{away} |")
        self.sections.append(("Average versus structure", L))
        for name, r0, rD, n in flagged:
            self.add("MEDIUM", "average-not-structure",
                     f"{name}: residuals of rms {r0:.1f} sigma shrink to {rD:.1f} when blurred over 50 keV -- the fit follows the average, "
                     f"the data resolve structure the model does not have", f"{n} points; look at the excitation functions, not chi2/N")
        # MeV-wide levels spanning several unfitted peaks
        for lv in self.levels:
            if lv.gamma_mev is None or lv.gamma_mev < 0.3 or not self._in_data(lv, margin=0.0):
                continue
            inside = [c for c in strong if abs(c["ex"] - lv.energy) <= lv.gamma_mev]
            if len(inside) >= 3:
                self.add("HIGH" if lv.gamma_mev >= 1.0 else "MEDIUM", "broad-level-averages-peaks",
                         f"{lv.label}: a level {_fmt_width(lv.gamma_mev)} wide spans {len(inside)} visible unfitted peaks "
                         f"(Ex {', '.join(format(c['ex'], '.3f') for c in inside[:8])}{'...' if len(inside) > 8 else ''}) -- a broad pole fitting "
                         f"the average of narrow resonances", "each peak seen in >= 2 detectors with the fit > 2 sigma low at its top")

    # ---- 6b. zigzag in the calculated curve ("bumpiness": a saw-tooth over a narrow resonance) ------
    def fit_zigzag(self, min_steps: int = 4, min_rel: float = 0.005):
        t = self.tab
        found = []
        for k, lab, s in t.curves(min_points=min_steps + 2):
            ex = t.ex[s]; c = t.fit[s]
            if np.ptp(ex) < 1e-6 or not np.all(np.isfinite(c)) or np.any(c == 0):
                continue
            d = np.diff(c); rel = np.abs(d) / np.abs(c[:-1])
            sg = np.sign(d)
            i = 0
            while i < d.size - 1:
                if sg[i] != 0 and sg[i + 1] == -sg[i] and rel[i] > min_rel and rel[i + 1] > min_rel:
                    j = i + 1
                    while j + 1 < d.size and sg[j + 1] == -sg[j] and rel[j + 1] > min_rel:
                        j += 1
                    # bumpiness is a sub-resolution artefact: only runs over finely spaced points (<= 10 keV) count,
                    # coarse points alternating across real resonances are the data's sampling, not the curve's fault
                    if j - i + 1 >= min_steps and float(np.median(np.diff(ex[i:j + 2]))) <= 0.010:
                        amp = float(np.max(rel[i:j + 1]))
                        lo, hi = float(ex[i]), float(ex[j + 1])
                        narrow = [lv for lv in self.levels if lo <= lv.energy <= hi and lv.gamma_mev is not None
                                  and lv.gamma_mev < float(np.median(np.diff(ex)))]
                        found.append(dict(key=k, ang=lab, name=str(t.name[s[0]]), lo=lo, hi=hi, steps=j - i + 1, amp=amp, narrow=narrow))
                    i = j + 1
                else:
                    i += 1
        found.sort(key=lambda d: -d["amp"] * d["steps"])
        self.zigzags = found
        L = [f"alternating up/down steps of the CALCULATED curve between consecutive data energies spaced <= 10 keV (>= {min_steps} steps, "
             f"each > {100*min_rel:g} %): {len(found)} found" + ("" if not found else " -- a smooth model does not do this unless a level "
             "narrower than the point spacing sits there; otherwise it is the integration or convolution grid"), ""]
        if found:
            L += ["| data set | seg | angle | Ex from | Ex to | steps | max step | narrow levels inside |", "|---|---:|---:|---:|---:|---:|---:|---|"]
            for d in found[:self.top]:
                ang = "" if d["ang"] < 0 else f"{d['ang']:.0f}"
                L.append(f"| {d['name']} | {d['key']} | {ang} | {d['lo']:.4f} | {d['hi']:.4f} | {d['steps']} | {100*d['amp']:.1f} % | "
                         f"{'; '.join(lv.label for lv in d['narrow'][:3])} |")
        self.sections.append(("Zigzag in the calculated curve", L))
        seen = set()
        for d in found:
            if len(seen) >= 5:
                break
            if d["name"] in seen:
                continue
            seen.add(d["name"])
            hint = (f"levels narrower than the point spacing inside ({'; '.join(lv.label for lv in d['narrow'][:3])}): the curve may be real "
                    f"but is under-sampled by the data -- check the target-integration / convolution grid" if d["narrow"] else
                    "no narrow level inside: a numerical artefact of the integration or convolution grid, not physics")
            ang = "" if d["ang"] < 0 else f", {d['ang']:.0f} deg"   # not nested: Python < 3.12
            self.add("MEDIUM", "fit-zigzag",
                     f"{d['name']} (seg {d['key']}{ang}): the calculated curve alternates up and down "
                     f"over {d['steps']} consecutive steps at Ex {d['lo']:.4f}-{d['hi']:.4f} MeV (steps up to {100*d['amp']:.1f} %); {hint}",
                     "the GUI shows it as bumpiness on the fit curve")

    # ---- 7. level support ---------------------------------------------------------------------------------------
    def level_support(self):
        t = self.tab
        ti_keys = set()
        if PYAZR3:
            # the data sets the engine folds (their record names a target, a beam or a detector)
            ti_keys = {int(s.key) for s in (self.segments or []) if getattr(s, "active", True) and getattr(s, "folded", False)}
        elif self.azr_text:
            from .model import target_int_entries
            try:
                for active, keys, _ in target_int_entries(self.azr_text):
                    if active:
                        ti_keys.update(int(k) for k in keys)
            except Exception:
                pass
        L = ["| level | Ex | Gamma | data sets within +-max(3 Gamma, 5 keV) | points | local Ex spacing | target-int segments |",
             "|---|---:|---:|---|---:|---:|---|"]
        nflag = 0
        for lv in self.levels:
            if not self._in_data(lv, margin=0.0):
                continue
            g = lv.gamma_mev
            half = max(3 * g, 0.005) if g is not None else 0.005
            sel = np.abs(t.ex - lv.energy) <= half
            names = sorted(set(str(n) for n in t.name[sel]))
            keys = sorted(set(int(k) for k in t.key[sel]))
            spac = []
            for k in keys:
                for kk, lab, s in t.curves(min_points=3):
                    if kk != k:
                        continue
                    ex = t.ex[s]
                    near = np.abs(ex - lv.energy) <= max(10 * half, 0.05)
                    if near.sum() >= 2:
                        spac.append(float(np.median(np.diff(ex[near]))))
            sp = float(np.median(spac)) if spac else float("nan")
            ti = sorted(k for k in keys if k in ti_keys)
            L.append(f"| {lv.jpi} | {lv.energy:.4f} | {_fmt_width(g)} | {', '.join(names)[:70]} | {int(sel.sum())} | "
                     f"{'' if not np.isfinite(sp) else f'{1e3*sp:.1f} keV'} | {','.join(map(str, ti))} |")
            if g is None:
                continue
            if g > 0 and np.isfinite(sp) and g < 0.5 * sp and sel.sum() > 0:
                nflag += 1
                # a level under 2 % of the point spacing is a delta function on the data grid: it acts on the objective
                # only through the integration (levels of 13 eV and 0.8 eV against 1.5-2 keV spacing)
                sev = "MEDIUM" if g < 0.02 * sp else "LOW"
                self.add(sev, "level-below-resolution",
                         f"{lv.label}: {'under 2 % of' if sev == 'MEDIUM' else 'narrower than half'} the local data spacing ({1e3*sp:.1f} keV)" +
                         (f"; segments {','.join(map(str, ti))} are target-integrated, so the objective can live on the integration "
                          f"grid rather than on the data -- confirm with a converged grid (rmfit review --grid-check)" if ti else
                          "; the data cannot constrain its width -- check what the level is doing for the fit"),
                         f"{int(sel.sum())} points within +-{1e3*half:.1f} keV from {len(names)} data set(s)")
            elif sel.sum() == 0 and g < 1.0:
                self.add("LOW", "level-without-data", f"{lv.label}: no data point within +-{1e3*half:.1f} keV",
                         "a level the data do not see directly (or a background pole placed inside the data range)")
            elif len(names) == 1 and g < 0.3:
                self.add("LOW", "one-dataset-level", f"{lv.label}: supported by one data set only ({names[0]})",
                         f"{int(sel.sum())} points within +-{1e3*half:.1f} keV; provisional until a second set or the literature confirms it")
        self.sections.append(("Levels inside the data range", L))

    # ---- close same-J^pi pairs (gap tied to the width) ------------------------------------------------------
    def close_pairs(self):
        """Neighbouring levels of one J^pi inside the data range whose gap is below ``close_fraction`` x the narrower
        total width.  Report only: such a pair is a merge candidate for the evaluator (two levels that close are told
        apart only by a different pattern of partial widths), it is not moved, tied or merged here."""
        f = self.close_fraction
        if not f or f <= 0:
            return
        by = {}
        for lv in self.levels:
            if lv.gamma_mev is not None and lv.gamma_mev > 0 and np.isfinite(lv.gamma_mev) and self._in_data(lv, margin=0.0):
                by.setdefault(lv.jpi, []).append(lv)
        L = [f"| J^pi | Ex (1) | Gamma (1) | Ex (2) | Gamma (2) | gap | gap / narrower Gamma |", "|---|---:|---:|---:|---:|---:|---:|"]
        for jpi, ls in by.items():
            ls = sorted(ls, key=lambda l: l.energy)
            for a, b in zip(ls, ls[1:]):
                gap = b.energy - a.energy; g = min(a.gamma_mev, b.gamma_mev)
                if gap >= f * g:
                    continue
                gtxt = f"{1e3 * gap:.3g} keV"
                L.append(f"| {jpi} | {a.energy:.5f} | {_fmt_width(a.gamma_mev)} | {b.energy:.5f} | {_fmt_width(b.gamma_mev)} | {gtxt} | {gap / g:.3g} |")

                def chans(l):
                    return ", ".join(f"p{c[0]} L{c[1]} S{c[2]:g} {_fmt_width(c[3])}" for c in l.channels if c[3]) or "no open-channel widths listed"
                self.add("MEDIUM", "level-close-width",
                         f"two {jpi} levels at {a.energy:.5f} and {b.energy:.5f} MeV: gap {gtxt} is "
                         f"{100 * gap / g:.2g} % of the narrower width ({_fmt_width(g)}), below the {100 * f:g} % scale -- "
                         f"a merge candidate for the evaluator (compare the partial-width patterns; test with a merged refit AND "
                         f"a refit of the pair from the merged solution before deciding)",
                         f"({a.energy:.5f}: {chans(a)}) vs ({b.energy:.5f}: {chans(b)})")
        if len(L) > 2:
            self.sections.append((f"Close same-J^pi pairs (gap below {100 * f:g} % of the narrower width)", L))

    # ---- figures ------------------------------------------------------------------------------------------------
    def figures(self, out_dir: str, nfig: int = 6):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(out_dir, exist_ok=True)
        t = self.tab
        paths = []
        # residual map
        fig, ax = plt.subplots(figsize=(10, 3.2))
        r = t.resid
        ax.scatter(t.ex, r, s=4, c="0.6", label="|r| <= 3")
        big = np.abs(r) > 3
        ax.scatter(t.ex[big], r[big], s=8, c="C3", label="|r| > 3")
        for lv in self.levels:
            if self._in_data(lv, margin=0.0):
                ax.axvline(lv.energy, color="C0", lw=0.4, alpha=0.5)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("Ex (MeV)"); ax.set_ylabel("(fit - n d) / (n e)")
        ax.set_ylim(-max(8, np.percentile(np.abs(r), 99.5)), max(8, np.percentile(np.abs(r), 99.5))); ax.legend(loc="upper right", fontsize=7)
        ax.set_title("residuals of every point; vertical lines = levels", fontsize=9)
        p = os.path.join(out_dir, "review_residuals.png"); fig.tight_layout(); fig.savefig(p, dpi=110); plt.close(fig); paths.append(p)
        # the worst data sets, one figure each: excitation functions per detector (up to 8), data +- err and fit
        rows = getattr(self, "per_point_rows", [])[:nfig]
        for d in rows:
            # excitation functions per detector; a segment at one energy (an angular distribution) is one panel vs angle
            curves = []
            for k in d["keys"]:
                sel = np.flatnonzero(t.key == k)
                if np.ptp(t.ex[sel]) < 2e-3 and np.isfinite(t.ang[sel]).all():
                    curves.append((k, -1.0, sel[np.argsort(t.ang[sel])]))
                else:
                    curves += [(kk, lab, s) for kk, lab, s in t.curves(min_points=2) if kk == k]
            curves.sort(key=lambda c: -float(t.chi2[c[2]].sum()))
            curves = curves[:8]
            if not curves:
                continue
            ncol = 2 if len(curves) > 1 else 1; nrow = int(math.ceil(len(curves) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.6 * nrow), squeeze=False)
            for ax, (k, lab, s) in zip(axes.flat, curves):
                ex = t.ex[s]; nd = t.norm[s] * t.data[s]; ne = t.norm[s] * t.err[s]; c = t.fit[s]
                single_e = np.ptp(ex) < 2e-3 and np.isfinite(t.ang[s]).all()
                xx = t.ang[s] if single_e else ex
                o = np.argsort(xx); xx, nd, ne, c = xx[o], nd[o], ne[o], c[o]
                if single_e:
                    ax.text(0.02, 0.92, f"Ex {ex.mean():.4f} MeV", transform=ax.transAxes, fontsize=7)
                ax.errorbar(xx, nd, yerr=ne, fmt="o", ms=2.5, color="k", lw=0.6, label="data x n")
                ax.plot(xx, c, "-", color="C3", lw=1.0, label="fit")
                pos = nd[nd > 0]
                if pos.size and pos.max() / max(pos.min(), 1e-300) > 30 and (c > 0).all():
                    ax.set_yscale("log")
                if not single_e:
                    for lv in self.levels:
                        if ex.min() <= lv.energy <= ex.max():
                            ax.axvline(lv.energy, color="C0", lw=0.4, alpha=0.5)
                ax.set_xlabel("angle (deg)" if single_e else "Ex (MeV)", fontsize=8)
                ax.set_title(f"seg {k}" + ("" if lab < 0 else f"  {lab:.0f} deg") + f"  chi2 {float(t.chi2[s].sum()):.0f} / {s.size}", fontsize=8)
                ax.tick_params(labelsize=7)
            for ax in list(axes.flat)[len(curves):]:
                ax.axis("off")
            axes.flat[0].legend(fontsize=7)
            fig.suptitle(f"{d['name']}: chi2/N {d['per_pt']:.1f}, {100*d['f2']:.0f} % within 2 sigma", fontsize=9)
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", d["name"])
            p = os.path.join(out_dir, f"review_{safe}.png"); fig.tight_layout(); fig.savefig(p, dpi=110); plt.close(fig); paths.append(p)
        self.figure_paths = paths
        self.sections.append(("Figures", [f"- `{os.path.relpath(p, os.path.dirname(out_dir))}`" for p in paths]))

    # ---- report -------------------------------------------------------------------------------------------------
    def write(self, path: str, title: str = "Fit review (rmfit review)", header_lines: Sequence[str] = ()) -> dict:
        sev_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "INFO": 3}
        fs = sorted(self.findings, key=lambda f: sev_order.get(f.severity, 9))
        counts = {s: sum(1 for f in fs if f.severity == s) for s in sev_order}
        L = [f"# {title}", "", f"{time.strftime('%Y-%m-%d %H:%M')}  " + ", ".join(f"{k} {v}" for k, v in counts.items()), ""]
        L += list(header_lines)
        L += [f"- {n}" for n in self.tab.notes]
        L += ["", "Findings are candidates for the evaluator (CLAUDE.md: judge per point, leave exclusions and data treatment to the "
              "evaluator).  Plain chi2 throughout: data covariances and the robust loss are ignored here.", ""]
        L += [f.line() for f in fs] or ["- no findings"]
        for title_, lines in self.sections:
            L += ["", f"## {title_}", ""] + list(lines)
        with open(path, "w") as fh:
            fh.write("\n".join(L) + "\n")
        return counts


# -- entry points ---------------------------------------------------------------------------------------------

def review_output(azr_path: str, output_dir: Optional[str] = None, out_dir: Optional[str] = None,
                  figures: int = 6, top: int = 20, close_fraction: Optional[float] = None) -> dict:
    """Engine-free review of a finished CLI run: ``AZUREOut_*.out`` + ``chiSquared.out`` + ``parameters.out``.
    Writes ``<out_dir>/REVIEW.md`` and ``<out_dir>/figures/review_*.png`` (default ``<azr dir>/review/``)."""
    from .model import base_segments
    azr_path = os.path.abspath(azr_path)
    odir = output_dir or output_dir_of(azr_path)
    tab = table_from_output(azr_path, odir)
    pout = os.path.join(odir, "parameters.out")
    if os.path.exists(pout):
        levels = levels_from_parameters_out(pout); lsrc = pout
    else:
        levels = levels_from_azr(azr_path); lsrc = azr_path + " (energies only)"
    segs = base_segments(azr_path)
    out_dir = out_dir or os.path.join(os.path.dirname(azr_path), "review")
    os.makedirs(out_dir, exist_ok=True)
    rv = Review(tab, levels, segs, open(azr_path).read(), top=top, close_fraction=close_fraction).run(os.path.join(out_dir, "figures") if figures else None, figures)
    path = os.path.join(out_dir, "REVIEW.md")
    counts = rv.write(path, header_lines=[f"model `{azr_path}`", f"levels from `{lsrc}` ({len(levels)})", ""])
    return dict(path=path, counts=counts, findings=rv.findings, review=rv)


def review_campaign(campaign, grid_check: bool = False, figures: int = 6, top: int = 20, cid=None, path=None,
                    close_fraction: Optional[float] = None) -> dict:
    """Review of a campaign's incumbent through a fresh single session; writes ``<campaign>/REVIEW.md`` and
    ``figures/review_*.png``.  ``grid_check`` re-scores the incumbent on a converged target-integration grid
    (20/50, ``Campaign.set_grid_setting``) and flags a shift of the objective.
    ``cid`` reviews that ledger candidate instead (on its own structure's model), ``path`` writes the report elsewhere
    (2026-11-10: the drive reviews the init baseline once, to tell a HIGH flag the base model already carried from a new one)."""
    from .engine import SingleSession
    from .model import single_variant, base_segments
    from .keys import take
    from .diagnose import theta2_table
    c = campaign
    if cid is None:
        cid, x, obj, keys = c.incumbent()
        base_azr = c.base_azr
    else:
        x, keys = c.ledger.load(cid)
        base_azr = c.ledger.structure(c.ledger.candidate(cid)["structure_id"])["azr"]
        obj = float("nan")
    src = single_variant(base_azr, "review" if path is None else "review_" + os.path.splitext(os.path.basename(path))[0], c.vdir)
    with SingleSession(src, cwd=c.vdir, flags=c.flags, name="review", objective=c.objective) as ev:
        xb = take(x, keys, ev.keys, ev.x0)
        tab = table_from_evaluator(ev, xb)
        levels = levels_from_evaluator(ev, xb)
        sc = ev.score(xb)
        if obj != obj:
            obj = sc.objective
        header = [f"campaign `{c.dir}`, {'incumbent ' if path is None else ''}candidate {cid}, objective {obj:.3f}",
                  f"single-session score: objective {sc.objective:.3f} = chi2 {sc.chi2:.3f} + penalty {sc.penalty:.3f}", f"model `{base_azr}`", ""]
        th = theta2_table(ev, xb)[:10]
    rv = Review(tab, levels, base_segments(base_azr), open(base_azr).read(), top=top, close_fraction=close_fraction)
    rv.run(os.path.join(c.dir, "figures") if figures else None, figures)
    if th:
        rv.sections.append(("Largest theta^2 (free widths vs the Wigner limit)",
                            ["| key | theta2 | rwa |", "|---|---:|---:|"] + [f"| {d['key']} | {d['theta2']:.3f} | {d['rwa']:.4g} |" for d in th]))
    if grid_check:
        try:
            gsrc = single_variant(c.base_azr, "review_grid", c.vdir)
            nchg = c.set_grid_setting(gsrc, 20.0, 50.0)
            if nchg:
                with SingleSession(gsrc, cwd=c.vdir, flags=c.flags, name="review_grid", objective=c.objective) as ev2:
                    sc2 = ev2.score(take(x, keys, ev2.keys, ev2.x0))
                rel = (sc2.objective - sc.objective) / max(abs(sc.objective), 1e-300)
                rv.sections.append(("Target-integration grid check",
                                    [f"objective {sc.objective:.3f} at the campaign grid, {sc2.objective:.3f} at the converged grid 20/50 "
                                     f"({100*rel:+.2f} %, {nchg} targetInt lines changed)"]))
                if abs(rel) > 0.01:
                    rv.add("HIGH", "grid-exploited", f"the objective moves by {100*rel:+.1f} % on a converged target-integration grid "
                           f"({sc.objective:.1f} -> {sc2.objective:.1f}) -- narrow levels are living on the integration grid, not on the data "
                           f"(seen once: the objective doubled)", "re-fit at the converged grid before reporting any number")
                elif abs(rel) > 0.002:
                    rv.add("MEDIUM", "grid-sensitive", f"the objective moves by {100*rel:+.2f} % on a converged target-integration grid",
                           f"{sc.objective:.1f} -> {sc2.objective:.1f}; acceptable for screening, not for a reported fit")
            else:
                rv.sections.append(("Target-integration grid check", ["no active <targetInt> line: nothing to check"]))
        except Exception as err:
            rv.sections.append(("Target-integration grid check", [f"(failed: {err})"]))
    path = path or os.path.join(c.dir, "REVIEW.md")
    counts = rv.write(path, header_lines=header)
    c.log(f"review: {path} " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    return dict(path=path, counts=counts, findings=rv.findings, review=rv)


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="engine-free review of a finished AZURE2 run (AZUREOut_*.out): per-point summary, "
                                            "norm direction and trends, residual runs, average-vs-structure, level support")
    p.add_argument("azr")
    p.add_argument("--output", help="output directory (default: the .azr's <fileDir> or output/ next to it)")
    p.add_argument("--out", help="where REVIEW.md and figures/ go (default <azr dir>/review/)")
    p.add_argument("--figures", type=int, default=6, help="figures for the N worst data sets (0 = none)")
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--close-fraction", type=float, default=None,
                   help=f"flag same-J^pi pairs with a gap below this fraction of the narrower width (default {CLOSE_GAP_FRACTION}; 0 = off)")
    a = p.parse_args(argv)
    r = review_output(a.azr, a.output, a.out, figures=a.figures, top=a.top, close_fraction=a.close_fraction)
    for f in r["findings"]:
        print(f.line())
    print("wrote", r["path"], r["counts"])


if __name__ == "__main__":
    main()
