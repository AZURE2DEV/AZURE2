"""Energy-window variants and the energy ladder.

A window variant keeps only the data inside an excitation-energy range and
fixes every level outside it (plus a margin), so a local question -- does a
level belong at 10.55 MeV, which J^pi, which interference sign -- is answered
on a model with a few thousand points and a few dozen free parameters, in
seconds per evaluation, before the full model is asked.

The ladder (an archived campaign) is the same device used
for a from-scratch fit: Ex <= 7.5 / 8.5 / 9.5 / 10.3 / 11.2 MeV, each step
freeing the levels up to its cutoff cumulatively and warm-starting from the
previous step's parameters (moved between variants by key).
"""
import dataclasses
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import AZURE2_ROOT  # noqa: F401
from .backend import PYAZR3
from .model import rows_in_segment
if PYAZR3:
    from .azr3 import SegmentTable
else:
    from .model import block_lines, replace_block, set_output_dir_text
    from .model import cm_conversion as _cm_conversion_p2
    from pyazr import AzrModel
    from pyazr.datasets import SegmentSet


def _pair(project_or_path, pair_key: int):
    import pyazr
    project = (pyazr.Project.open(os.path.abspath(project_or_path)) if isinstance(project_or_path, str)
               else project_or_path)
    for q in project.pairs():
        if int(q.key) == int(pair_key):
            return q
    raise ValueError(f"no pair with key {pair_key}")


def cm_conversion(text, pair_key: int) -> float:
    """E_cm / E_lab of a pair with the light particle as the beam (pyazr 2: :func:`rmfit.model.cm_conversion` on the
    .azr text; pyazr 3: from ``Project.pairs()``)."""
    return _cm_conversion_p3(text, pair_key) if PYAZR3 else _cm_conversion_p2(text, pair_key)


def _cm_conversion_p3(project_or_path, pair_key: int) -> float:
    """E_cm / E_lab of a pair with the light particle as the beam: M_heavy / (M_light + M_heavy)."""
    q = _pair(project_or_path, pair_key)
    return float(q.heavy_mass) / (float(q.light_mass) + float(q.heavy_mass))


def pair_separation(text: str, pair_key: int) -> float:
    """pyazr 3 / pyazr 2 versions below (rmfit.backend): ``text`` is the .azr text (pyazr 2) or the project's path or Project (pyazr 3)."""
    return (pair_separation_p3(text, pair_key) if PYAZR3 else pair_separation_p2(text, pair_key))

def pair_separation_p2(text: str, pair_key: int) -> float:
    """Separation energy (MeV) of a particle pair from the <levels> lines."""
    for line in block_lines(text, "levels"):
        t = line.split()
        if len(t) >= 23 and int(float(t[5])) == pair_key:
            return float(t[22])
    raise ValueError(f"no <levels> line uses pair {pair_key}")

def pair_separation_p3(project_or_path, pair_key: int) -> float:
    """Separation energy (MeV) of a particle pair (``Project.pairs()``)."""
    return float(_pair(project_or_path, pair_key).separation_energy)


def excitation_bounds_to_lab(text: str, pair_key: int, ex_lo: float, ex_hi: float) -> Tuple[float, float]:
    """pyazr 3 / pyazr 2 versions below (rmfit.backend): ``text`` is the .azr text (pyazr 2) or the project's path or Project (pyazr 3)."""
    return (excitation_bounds_to_lab_p3(text, pair_key, ex_lo, ex_hi) if PYAZR3 else excitation_bounds_to_lab_p2(text, pair_key, ex_lo, ex_hi))

def excitation_bounds_to_lab_p2(text: str, pair_key: int, ex_lo: float, ex_hi: float) -> Tuple[float, float]:
    cm = cm_conversion(text, pair_key)
    sep = pair_separation(text, pair_key)
    return max(0.0, (ex_lo - sep) / cm), max(0.0, (ex_hi - sep) / cm)

def excitation_bounds_to_lab_p3(project_or_path, pair_key: int, ex_lo: float, ex_hi: float) -> Tuple[float, float]:
    cm = cm_conversion(project_or_path, pair_key)
    sep = pair_separation(project_or_path, pair_key)
    return max(0.0, (ex_lo - sep) / cm), max(0.0, (ex_hi - sep) / cm)


def window_variant(base_azr: str, path: str, ex_lo: float, ex_hi: float,
                   level_margin: float = 0.0, fix_outside: bool = True,
                   output_dir: Optional[str] = None,
                   free_levels_below: Optional[float] = None) -> Dict[str, object]:
    """pyazr 3 / pyazr 2 versions below (rmfit.backend): ``text`` is the .azr text (pyazr 2) or the project's path or Project (pyazr 3)."""
    return (window_variant_p3(base_azr, path, ex_lo, ex_hi, level_margin, fix_outside, output_dir, free_levels_below) if PYAZR3 else window_variant_p2(base_azr, path, ex_lo, ex_hi, level_margin, fix_outside, output_dir, free_levels_below))

def window_variant_p2(base_azr: str, path: str, ex_lo: float, ex_hi: float,
                   level_margin: float = 0.0, fix_outside: bool = True,
                   output_dir: Optional[str] = None,
                   free_levels_below: Optional[float] = None) -> Dict[str, object]:
    """Write ``path``: ``base_azr`` with every active segment capped to
    ``[ex_lo, ex_hi]`` (in excitation energy, converted per entrance pair),
    segments entirely outside deactivated, and -- if ``fix_outside`` -- every
    level outside ``[ex_lo - margin, ex_hi + margin]`` fixed (energy and all
    widths).  ``free_levels_below`` (ladder mode) instead fixes every level
    above that energy and leaves the rest as in the base file.

    Returns a summary dict (segments kept/dropped, levels fixed).
    """
    text = open(base_azr).read()
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    segs = SegmentSet.from_file(base_azr)
    raw = block_lines(text, "segmentsData")
    kept, dropped, lines = 0, 0, []
    for seg, line in zip(segs, raw):
        t = line.split()
        if not seg.active:
            lines.append(line)
            continue
        lo, hi = excitation_bounds_to_lab(text, seg.entrance_key, ex_lo, ex_hi)
        new_lo, new_hi = max(seg.energy_min, lo), min(seg.energy_max, hi)
        # a segment left without data points is dropped by the engine with a
        # warning, which would throw the segment bookkeeping off: deactivate it
        n_in = 0
        if new_hi > new_lo:
            capped = dataclasses.replace(seg, energy_min=new_lo, energy_max=new_hi)
            n_in = int(rows_in_segment(capped, base_dir).shape[0])
        if n_in == 0:
            t[0] = "0"
            dropped += 1
        else:
            t[3], t[4] = repr(float(new_lo)), repr(float(new_hi))
            kept += 1
        lines.append("  ".join(t))
    text = replace_block(text, "segmentsData", lines)
    text = replace_block(text, "segmentsTest", [])
    if output_dir is not None:
        text = set_output_dir_text(text, output_dir)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    mdl = AzrModel.from_file(tmp)
    fixed_levels = []
    for lv in mdl.levels:
        if not lv.active:
            continue
        if free_levels_below is not None:
            outside = lv.energy > free_levels_below
        else:
            outside = fix_outside and (lv.energy < ex_lo - level_margin or lv.energy > ex_hi + level_margin)
        if outside:
            lv.set_fixed(True)
            for c in lv.channels:
                c.channel_fixed = True
            fixed_levels.append(f"{lv.jpi}@{lv.energy:.4f}")
    mdl.write(path)
    os.remove(tmp)
    return dict(path=path, segments_kept=kept, segments_dropped=dropped, levels_fixed=fixed_levels)

def window_variant_p3(base_azr: str, path: str, ex_lo: float, ex_hi: float,
                   level_margin: float = 0.0, fix_outside: bool = True,
                   output_dir: Optional[str] = None,
                   free_levels_below: Optional[float] = None) -> Dict[str, object]:
    """Write ``path``: ``base_azr`` with every active data set capped to
    ``[ex_lo, ex_hi]`` (in excitation energy, converted per entrance pair),
    data sets entirely outside switched off (``active no``), the extrapolations removed, and -- if
    ``fix_outside`` -- every level outside ``[ex_lo - margin, ex_hi + margin]`` fixed (energy and all
    widths).  ``free_levels_below`` (ladder mode) instead fixes every level
    above that energy and leaves the rest as in the base file.

    A data set whose rows rmfit cannot place (:func:`rmfit.azr3.read_points`) is refused: its window cannot be
    converted.  Returns a summary dict (data sets kept/dropped, levels fixed).
    """
    import pyazr
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    project = pyazr.Project.open(os.path.abspath(base_azr))
    segs = SegmentTable.from_file(base_azr)
    kept, dropped = 0, 0
    for seg in segs:
        if not seg.active:
            continue
        lo, hi = excitation_bounds_to_lab(project, seg.entrance_key, ex_lo, ex_hi)
        new_lo, new_hi = max(seg.energy_min, lo), min(seg.energy_max, hi)
        # a data set left without points would throw the bookkeeping off: switch it off
        n_in = 0
        if new_hi > new_lo:
            capped = dataclasses.replace(seg, energy_min=new_lo, energy_max=new_hi)
            rows = rows_in_segment(capped, base_dir)
            if rows is None:
                raise NotImplementedError(f"window_variant: data set {seg.name}: its file's columns "
                                          f"({seg.get('columns')!r}) are not lab energies rmfit can window")
            n_in = int(rows.shape[0])
        if n_in == 0:
            project.set(seg.name, active="no")
            dropped += 1
        else:
            project.set(seg.name, energy_window=(repr(float(new_lo)), repr(float(new_hi)), "MeV"))
            kept += 1
    for name, rec in list(project.records().items()):
        if rec.record_kind == "extrapolation":
            project.remove(name)
    if output_dir is not None:
        project.set_config(output_dir=output_dir.rstrip("/") + "/")
    fixed_levels = []
    for lv in project.levels():
        if not lv.active:
            continue
        if free_levels_below is not None:
            outside = lv.energy > free_levels_below
        else:
            outside = fix_outside and (lv.energy < ex_lo - level_margin or lv.energy > ex_hi + level_margin)
        if outside:
            project.set_level(lv, fixed=True)
            for c in lv.channels:
                if not c.fixed:
                    project.set_channel(lv, c.pair, c.l, c.s, fixed=True)
            fixed_levels.append(f"{lv.jpi}@{lv.energy:.4f}")
    project.save(path)
    return dict(path=path, segments_kept=kept, segments_dropped=dropped, levels_fixed=fixed_levels)


LADDER = [7.5, 8.5, 9.5, 10.3, 11.2]


def ladder_variants(base_azr: str, workdir: str, cutoffs: Sequence[float] = LADDER,
                    ex_lo: float = 0.0) -> List[Dict[str, object]]:
    """One variant per ladder step: data up to the cutoff, levels above it fixed."""
    os.makedirs(workdir, exist_ok=True)
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    from .model import link_data_dirs
    link_data_dirs(base_azr, workdir)
    out = []
    for i, cut in enumerate(cutoffs):
        name = f"step{i+1}_Ex{cut:g}"
        path = os.path.join(workdir, name + ".azr")
        info = window_variant(base_azr, path, ex_lo, cut, output_dir=os.path.join(workdir, f"output_{name}"),
                              free_levels_below=cut)
        info["cutoff"] = cut
        out.append(info)
    return out
