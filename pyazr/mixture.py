"""Helpers for mixture segments (natural-abundance and compound targets).

A mixture segment's model is ``w_self * (this project's cross section) +
sum_i w_i * fold(T_i)`` where each ``T_i`` is a tabulated cross section of
another isotope (an ENDF/JEFF point-wise table, or another AZURE2 project's
extrapolation) read from a two-column file::

    # provenance
    E_lab[MeV]   sigma[b]

in the projectile's LAB energy -- the axis the data file uses -- and in the
frame AZURE2 compares the segment in (barns for angle-integrated and total
cross sections, c.m. b/sr for a differential one).  AZURE2 folds each table
with the segment's own resolution function (``ESegment::ComputeMixtureExternal``).

This module writes such tables.  :func:`extrap_to_lab_table` converts an
``AZUREOut_*.extrap`` file, whose first column is the c.m. energy, with the
entrance pair's masses (:func:`pair_masses` reads them from the other
project's ``.azr``).  :func:`write_table` writes any ``(E_lab, sigma)`` arrays.
The block on the segment line is edited with ``AzrModel.set_mixture`` and
``set_test_mixture``.
"""
import os
import time

import numpy as np


def pair_masses(model, pair_key):
    """``(M1, M2)`` in u of particle pair ``pair_key`` of an ``AzrModel``
    (the light and the heavy particle, as on the ``<levels>`` lines)."""
    for level in model.levels:
        for ch in level.channels:
            if ch.pair == int(pair_key):
                return ch.M1, ch.M2
    raise KeyError(f"no channel of pair {pair_key} in the model.")


def lab_from_cm(e_cm, m1, m2):
    """Projectile lab energy from the entrance-channel c.m. energy,
    non-relativistic, as ``EPoint::ConvertLabEnergy``: E_cm = E_lab*M2/(M1+M2)."""
    return np.asarray(e_cm, dtype=float) * (m1 + m2) / m2


def read_extrap(path, block=0):
    """``(E, sigma)`` of one blank-line-separated block of an
    ``AZUREOut_*.extrap`` (or ``.out``) file: column 1 and column 4."""
    blocks, cur = [], []
    with open(path) as f:
        for line in f:
            if not line.strip():
                if cur:
                    blocks.append(cur); cur = []
                continue
            if line.lstrip().startswith("#"):
                continue
            cur.append(line.split())
    if cur:
        blocks.append(cur)
    if block >= len(blocks):
        raise IndexError(f"{path} has {len(blocks)} block(s); no block {block}.")
    rows = np.array([[float(r[0]), float(r[3])] for r in blocks[block]])
    return rows[:, 0], rows[:, 1]


def write_table(e_lab, sigma, out_path, header=None):
    """Write a mixture table.  Energies must increase strictly (rows are
    sorted and exact duplicates dropped first)."""
    e = np.asarray(e_lab, dtype=float)
    s = np.asarray(sigma, dtype=float)
    if e.shape != s.shape or e.ndim != 1:
        raise ValueError("e_lab and sigma must be 1-d arrays of equal length.")
    order = np.argsort(e, kind="stable")
    e, s = e[order], s[order]
    keep = np.concatenate(([True], np.diff(e) > 0))
    e, s = e[keep], s[keep]
    if len(e) < 2:
        raise ValueError("a mixture table needs at least two rows.")
    with open(out_path, "w") as f:
        f.write(f"# mixture table written by pyazr {time.strftime('%Y-%m-%d %H:%M')}\n")
        if header:
            for line in str(header).splitlines():
                f.write("# " + line + "\n")
        f.write("# E_lab[MeV]  sigma[b]\n")
        for ei, si in zip(e, s):
            f.write(f"{ei:.10e} {si:.10e}\n")
    return out_path


def extrap_to_lab_table(extrap_path, m1, m2, out_path, kinflag=False, block=0,
                        header=None):
    """Convert an AZURE2 extrapolation into a mixture table.

    ``extrap_path`` is an ``AZUREOut_aa=..._R=....extrap`` file written with
    the ``<targetInt>`` lines INACTIVE (a bare cross section: the table is
    folded again by the segment that uses it).  Column 1 is the entrance
    c.m. energy unless the file was written with ``kinflag.dat = 1``
    (``kinflag=True``), in which case it already holds the lab energy.
    ``m1, m2`` are the entrance pair's masses in u (:func:`pair_masses`).
    The grid should resolve every resonance of the other isotope: AZURE2
    interpolates the table linearly between its rows.
    """
    e, s = read_extrap(extrap_path, block)
    e_lab = e if kinflag else lab_from_cm(e, m1, m2)
    prov = f"from {os.path.abspath(extrap_path)}"
    if not kinflag:
        prov += f"; E_lab = E_cm*(M1+M2)/M2 with M1 = {m1!r}, M2 = {m2!r}"
    if header:
        prov = str(header) + "\n" + prov
    return write_table(e_lab, s, out_path, header=prov)


# -- ENDF-6 point-wise cross sections ------------------------------------------

def _endf_float(tok):
    """ENDF-6 number: Fortran style, often without the 'E' (1.234567+5)."""
    tok = tok.strip()
    if not tok:
        return 0.0
    try:
        return float(tok)
    except ValueError:
        pass
    body = tok[1:] if tok[0] in "+-" else tok
    for i in range(len(body) - 1, 0, -1):
        if body[i] in "+-" and body[i - 1] not in "eE":
            mant, exp = tok[:len(tok) - len(body) + i], body[i:]
            return float(mant + "e" + exp)
    raise ValueError(f"not an ENDF number: {tok!r}")


def _endf_sections(path):
    """``{(MF, MT): [lines]}`` of an ENDF-6 tape (its first material)."""
    sections = {}
    with open(path) as f:
        for line in f:
            if len(line) < 75:
                continue
            try:
                mf, mt = int(line[70:72]), int(line[72:75])
            except ValueError:
                continue
            if mf == 0 or mt == 0:
                continue
            sections.setdefault((mf, mt), []).append(line)
    return sections


def endf_resonance_ranges(path):
    """The resonance ranges declared in MF2/MT151 of an ENDF-6 file, as a
    list of ``(E_low_eV, E_high_eV, LRU, LRF)``.  LRU 0 means only a
    scattering radius is given, i.e. the MF3 cross sections are already
    point-wise; LRU 1 or 2 means MF3 holds only a background and the
    resonance range has to be reconstructed by a processing code."""
    lines = _endf_sections(path).get((2, 151))
    if not lines:
        return []
    ranges = []
    nis = int(lines[0][44:55])            # [ZA, AWR, 0, 0, NIS, 0]
    k = 1
    for _ in range(nis):
        ner = int(lines[k][44:55]); k += 1   # [ZAI, ABN, 0, LFW, NER, 0]
        for _ in range(ner):
            el, eh = _endf_float(lines[k][0:11]), _endf_float(lines[k][11:22])
            lru, lrf = int(lines[k][22:33]), int(lines[k][33:44])
            ranges.append((el, eh, lru, lrf))
            k += 1
            # Skip the range's body: cannot be done without parsing each LRF, so
            # stop after the first range unless it is the scattering-radius-only case.
            if lru != 0:
                return ranges
            k += 1  # the SPI/AP line of an LRU = 0 range
    return ranges


def read_endf_mf3(path, mt=1):
    """``(E_lab_MeV, sigma_b, interpolation)`` of the MF3 section ``mt`` of an
    ENDF-6 file; ``interpolation`` is the list of ``(NBT, INT)`` pairs of its
    TAB1 record (``INT == 2`` is lin-lin, what a mixture table assumes)."""
    lines = _endf_sections(path).get((3, mt))
    if not lines:
        raise KeyError(f"{path}: no MF3/MT{mt} section.")
    head = lines[1]                       # [QM, QI, 0, LR, NR, NP]
    nr, npts = int(head[44:55]), int(head[55:66])
    k = 2
    interp = []
    pairs = []
    while len(pairs) < nr:
        row = lines[k]; k += 1
        for j in range(3):
            if len(pairs) < nr:
                pairs.append((int(row[22 * j:22 * j + 11]), int(row[22 * j + 11:22 * j + 22])))
    interp = pairs
    e, s = [], []
    while len(e) < npts:
        row = lines[k]; k += 1
        for j in range(3):
            if len(e) < npts:
                e.append(_endf_float(row[22 * j:22 * j + 11]))
                s.append(_endf_float(row[22 * j + 11:22 * j + 22]))
    return np.array(e) * 1e-6, np.array(s), interp


def endf_to_table(path, out_path, mt=1, e_min=None, e_max=None, header=None):
    """Write the point-wise MF3/MT``mt`` cross section of an ENDF-6 file as a
    mixture table (lab MeV, barns), optionally cut to ``[e_min, e_max]``
    (the table's end rows are kept just outside the cut).

    Refuses a file whose MF2 declares a resolved or unresolved resonance
    range (LRU 1 or 2): there the MF3 section holds only a smooth background,
    and the resonance cross section must first be reconstructed with a
    processing code (NJOY/PREPRO), or the evaluation taken from a service
    that serves reconstructed point-wise data.  Also refuses a non-linear
    interpolation scheme, which AZURE2's linear interpolation of the table
    would not reproduce.
    """
    ranges = endf_resonance_ranges(path)
    bad = [r for r in ranges if r[2] != 0]
    if bad:
        raise ValueError(f"{path}: MF2 declares a resonance range {bad[0][0]:g}-{bad[0][1]:g} eV with "
                         f"LRU={bad[0][2]} (LRF={bad[0][3]}); its MF3/MT{mt} is a background only -- "
                         "reconstruct the point-wise cross section first.")
    e, s, interp = read_endf_mf3(path, mt)
    if any(i != 2 for _, i in interp):
        raise ValueError(f"{path}: MF3/MT{mt} uses interpolation scheme(s) {sorted(set(i for _, i in interp))}; "
                         "only lin-lin (2) can be used as a mixture table as is.")
    keep = np.ones(len(e), dtype=bool)
    if e_min is not None:
        lo = np.searchsorted(e, e_min, side="right") - 1
        keep &= np.arange(len(e)) >= max(lo, 0)
    if e_max is not None:
        hi = np.searchsorted(e, e_max, side="left") + 1
        keep &= np.arange(len(e)) <= min(hi, len(e) - 1)
    prov = f"ENDF-6 file {os.path.abspath(path)}, MF3 MT{mt}, {int(keep.sum())} of {len(e)} points, lin-lin"
    if header:
        prov = str(header) + "\n" + prov
    return write_table(e[keep], s[keep], out_path, header=prov)


# -- Breit-Wigner resonance ranges (ENDF LRF 1/2): the total cross section ------

def _endf_lrf2_blocks(lines, k):
    """Parse one LRU=1, LRF=1/2 range starting at line index k (the SPI/AP
    line).  Returns ``(spi, ap, naps, blocks, k_next)`` with blocks a list of
    ``(awri, qx, l, lrx, resonances[ER, AJ, GT, GN, GG, GF])``."""
    spi, ap = _endf_float(lines[k][0:11]), _endf_float(lines[k][11:22])
    nls = int(lines[k][44:55]); k += 1
    blocks = []
    for _ in range(nls):
        awri, qx = _endf_float(lines[k][0:11]), _endf_float(lines[k][11:22])
        l, lrx, nrs = int(lines[k][22:33]), int(lines[k][33:44]), int(lines[k][55:66]); k += 1
        res = []
        for _ in range(nrs):
            res.append([_endf_float(lines[k][11 * j:11 * j + 11]) for j in range(6)]); k += 1
        blocks.append((awri, qx, l, lrx, np.array(res)))
    return spi, ap, blocks, k


def _penetrability_shift_phase(l, rho):
    """Hard-sphere penetrability P_l, shift S_l and phase phi_l (ENDF-102 App. D, Table D.1)."""
    r2 = rho * rho
    if l == 0:
        return rho, 0.0 * rho, rho
    if l == 1:
        return rho * r2 / (1 + r2), -1 / (1 + r2), rho - np.arctan(rho)
    if l == 2:
        d = 9 + 3 * r2 + r2 * r2
        return rho ** 5 / d, -(18 + 3 * r2) / d, rho - np.arctan(3 * rho / (3 - r2))
    if l == 3:
        d = 225 + 45 * r2 + 6 * r2 * r2 + r2 ** 3
        return rho ** 7 / d, -(675 + 90 * r2 + 6 * r2 * r2) / d, rho - np.arctan(rho * (15 - r2) / (15 - 6 * r2))
    if l == 4:
        d = 11025 + 1575 * r2 + 135 * r2 * r2 + 10 * r2 ** 3 + r2 ** 4
        return (rho ** 9 / d, -(44100 + 4725 * r2 + 270 * r2 * r2 + 10 * r2 ** 3) / d,
                rho - np.arctan(rho * (105 - 10 * r2) / (105 - 45 * r2 + r2 * r2)))
    raise ValueError(f"l = {l} not supported")


def breit_wigner_total(e_ev, spi, ap, blocks, naps=1):
    """Total cross section (b) of an ENDF SLBW/MLBW (LRF 1/2) resonance range
    at the lab energies ``e_ev`` (eV), ENDF-102 section D.1.1-D.1.2.  The
    total is additive over resonances in both formalisms (the inter-level
    interference of MLBW enters only the elastic), with the energy-dependent
    neutron width, the level shift and the hard-sphere phase of each l.
    ``ap`` (scattering radius, 1e-12 cm) gives the phase; the penetrability
    radius is ``ap`` too when NAPS=1, else 0.123 AWRI^(1/3) + 0.08.
    """
    e = np.asarray(e_ev, dtype=float)
    sigma = np.zeros_like(e)
    for awri, qx, l, lrx, res in blocks:
        kfac = 2.196771e-3 * awri / (awri + 1.0)         # k in 1e12/cm per sqrt(eV)
        k = kfac * np.sqrt(np.abs(e))
        a_pen = ap if naps == 1 else 0.123 * awri ** (1.0 / 3.0) + 0.08
        P, S, phi = _penetrability_shift_phase(l, k * a_pen)
        _, _, phi = _penetrability_shift_phase(l, k * ap)
        pot = 4 * np.pi / k ** 2 * (2 * l + 1) * np.sin(phi) ** 2 * 1e24 * 1e-24  # k^2 in 1e24/cm^2 -> barns
        sigma += pot
        for er, aj, gt, gn, gg, gf in res:
            g = (2 * aj + 1) / (2 * (2 * spi + 1))
            kr = kfac * np.sqrt(abs(er))
            Pr, Sr, _ = _penetrability_shift_phase(l, kr * a_pen)
            gn_e = gn * P / Pr
            gtot = gn_e + gg + gf + (gt - gn - gg - gf if lrx else 0.0)
            eshift = er + (Sr - S) * gn / (2 * Pr)
            d = e - eshift
            sigma += (np.pi / k ** 2) * g * gn_e * (gtot * np.cos(2 * phi) + 2 * d * np.sin(2 * phi)) / (d * d + gtot * gtot / 4)
    return sigma


def endf_reconstructed_total(path, grid_ev=None, points_per_width=40.0, max_points=400000):
    """``(E_lab_MeV, sigma_total_b, info)``: the total cross section of an
    ENDF-6 file with its Breit-Wigner resonance range (LRF 1/2) reconstructed
    at 0 K and added to the MF3/MT1 background, joined to the point-wise
    MF3/MT1 above the range.  ``grid_ev`` is an optional energy grid for the
    range; by default a union of a uniform grid and dense windows (+-15 total
    widths at ``points_per_width`` per width) around every resonance.  Other
    formalisms (Reich-Moore, R-matrix limited, unresolved) are refused.
    """
    sections = _endf_sections(path)
    lines = sections.get((2, 151))
    e_mf3, s_mf3, interp = read_endf_mf3(path, 1)
    if any(i != 2 for _, i in interp):
        raise ValueError(f"{path}: MF3/MT1 is not lin-lin")
    if not lines:
        return e_mf3, s_mf3, {"ranges": []}
    nis = int(lines[0][44:55]); k = 1
    info = {"ranges": []}
    e_all, s_all = [], []
    for _ in range(nis):
        ner = int(lines[k][44:55]); k += 1
        for _ in range(ner):
            el, eh = _endf_float(lines[k][0:11]), _endf_float(lines[k][11:22])
            lru, lrf, nro, naps = (int(lines[k][22:33]), int(lines[k][33:44]), int(lines[k][44:55]), int(lines[k][55:66])); k += 1
            if lru == 0:
                k += 1
                continue
            if lru == 2:
                # Unresolved range: MF3 holds the smooth (average) cross section there (TENDL files
                # set LSSF=1), and the average parameters only serve self-shielding.  It follows the
                # resolved range(s), so nothing after it is needed.
                info["ranges"].append({"E_low_eV": el, "E_high_eV": eh, "LRU": 2, "note": "unresolved: MF3 used"})
                break
            if lru != 1 or lrf not in (1, 2) or nro != 0:
                raise ValueError(f"{path}: resonance range {el:g}-{eh:g} eV with LRU={lru}, LRF={lrf}, NRO={nro} "
                                 "is not Breit-Wigner; reconstruct it with a processing code.")
            spi, ap, blocks, k = _endf_lrf2_blocks(lines, k)
            nres = sum(len(b[4]) for b in blocks)
            info["ranges"].append({"E_low_eV": el, "E_high_eV": eh, "LRF": lrf, "spin": spi, "AP": ap, "NAPS": naps,
                                   "l_blocks": [(b[2], len(b[4])) for b in blocks]})
            if grid_ev is None:
                g = [np.linspace(max(el, 1.0), eh, 4001)]
                for _, _, _, _, res in blocks:
                    for er, aj, gt, gn, gg, gf in res:
                        w = max(gt, 1.0)
                        n = int(min(30 * points_per_width, 4000))
                        g.append(np.linspace(er - 15 * w, er + 15 * w, n))
                grid = np.unique(np.concatenate(g))
                grid = grid[(grid >= max(el, 1.0)) & (grid <= eh)]
            else:
                grid = np.asarray(grid_ev, dtype=float)
                grid = grid[(grid >= el) & (grid <= eh)]
            if len(grid) > max_points:
                raise ValueError(f"{len(grid)} grid points; raise max_points or pass grid_ev")
            res_part = breit_wigner_total(grid, spi, ap, blocks, naps)
            bkg = np.interp(grid, e_mf3 * 1e6, s_mf3)
            e_all.append(grid); s_all.append(res_part + bkg)
            info["n_resonances"] = info.get("n_resonances", 0) + nres
    if not e_all:
        return e_mf3, s_mf3, info
    e_res = np.concatenate(e_all); s_res = np.concatenate(s_all)
    e_hi = e_res.max()
    above = e_mf3 * 1e6 > e_hi
    e_out = np.concatenate([e_res, e_mf3[above] * 1e6]) * 1e-6
    s_out = np.concatenate([s_res, s_mf3[above]])
    order = np.argsort(e_out)
    return e_out[order], s_out[order], info


def endf_total_to_table(path, out_path, e_min=None, e_max=None, header=None, **kw):
    """As :func:`endf_to_table` for MT=1, but with a Breit-Wigner resonance
    range reconstructed (:func:`endf_reconstructed_total`) instead of refused."""
    e, s, info = endf_reconstructed_total(path, **kw)
    keep = np.ones(len(e), dtype=bool)
    if e_min is not None:
        keep &= np.arange(len(e)) >= max(np.searchsorted(e, e_min, side="right") - 1, 0)
    if e_max is not None:
        keep &= np.arange(len(e)) <= min(np.searchsorted(e, e_max, side="left") + 1, len(e) - 1)
    prov = (f"ENDF-6 file {os.path.abspath(path)}, MF3 MT1 with the MF2 Breit-Wigner range(s) "
            f"{[(r['E_low_eV'], r['E_high_eV'], 'LRU', r.get('LRU', 1), 'LRF', r.get('LRF')) for r in info['ranges']]} reconstructed at 0 K "
            f"(pyazr.mixture.breit_wigner_total, {info.get('n_resonances', 0)} resonances), {int(keep.sum())} points")
    if header:
        prov = str(header) + "\n" + prov
    write_table(e[keep], s[keep], out_path, header=prov)
    return info
