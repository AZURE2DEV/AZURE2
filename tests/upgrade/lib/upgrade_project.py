"""Shared helpers for the upgrade test projects T1 and T3 (tests/upgrade/).

A project is described in Python (particle pairs, levels with their channels,
data segments on lab-energy grids, test segments), written as an ``.azr`` with
the true parameters, and its synthetic data are the engine's own observables at
those parameters with Gaussian noise from a stated seed.

The ``.azr`` written is project format 3 (docs/upgrade/project_format3.md): the
model is composed in the old grammar here (the <levels> columns, which format 3
keeps, and one old segment line per data set and extrapolation) and converted by
the engine's converter, ``AZURE2 --convert OLD NEW --to 3 --exact-grids``, so the
committed T1 and T3 are exactly the conversions the generators write.  Nothing here depends
on the machine: the AZURE2 binary comes from an argument, ``$AZURE2_BIN`` or
``$AZURE2_BINARY``, else the newest ``build*/src/AZURE2`` of the checkout.

Only the standard library and numpy are needed (no pyazr), so the generators
and the golden comparison run wherever the CLI runs.

The old lines are those of docs/source/reference/data_formats.rst (``<levels>``:
31 columns, ``io/include/azure2/io/LevelRows.h``; ``<segmentsData>``:
``convert/include/azure2/convert/SegLine.h``; ``<segmentsTest>``: the converter's
extrapolation lines).
"""

import glob
import math
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

# Atomic masses (u), AME2020 (Wang et al., Chin. Phys. C 45, 030003 (2021)).
MASS = {
    "n": 1.00866491595, "p": 1.00782503223, "e": 0.000548579909,
    "4He": 4.00260325413, "7Li": 7.0160034366, "7Be": 7.016928717,
    "8Li": 8.022486246, "8Be": 8.005305102, "12C": 12.0, "13N": 13.00573861,
}
U_MEV = 931.49410242  # MeV per u (CODATA 2018)


def separation(m1, m2, mc):
    """S = m1 + m2 - M (MeV) from atomic masses in u."""
    return (m1 + m2 - mc) * U_MEV


# ----------------------------------------------------------------------------- model


@dataclass
class Pair:
    """One particle pair (columns 13-21, 23, 27-31 of a level line)."""
    key: int
    j1: float
    pi1: int
    j2: float
    pi2: int
    m1: float
    m2: float
    z1: int
    z2: int
    S: float          # ground-state separation energy (MeV); Q value of a beta pair
    ptype: int = 0    # 0 particle, 10 gamma, 20 beta
    R: float = 0.0    # channel radius (fm)
    e2: float = 0.0   # heavy-particle excitation (MeV); final level of a gamma pair
    g1: float = 0.0
    g2: float = 0.0
    ecmask: int = 0   # external-capture multipolarities: E1 = 1, M1 = 2, E2 = 4


@dataclass
class Channel:
    pair: Pair
    s2: int           # 2 x channel spin
    l2: int           # 2 x orbital angular momentum (2 x multipolarity for a photon)
    width: float      # eV (open), fm^-1/2 ANC (closed), amplitude (beta)
    fixed: bool = True


@dataclass
class Level:
    J: float
    pi: int
    E: float          # excitation energy (MeV)
    channels: List[Channel]
    fixed: bool = True
    name: str = ""


@dataclass
class DataSegment:
    """A <segmentsData> line plus the points it is generated on."""
    name: str                     # data/<name>.dat
    entrance: int
    exit: int
    isdiff: int
    points: List[Tuple[float, float]]   # (E_lab MeV, angle deg) as in the data file
    noise: Tuple[str, float]      # ("rel", fraction) or ("abs", value)
    phase_jl: Optional[Tuple[float, int]] = None
    tail: str = ""                # keyword blocks after the file name
    keep: Optional[object] = None  # callable(truth_value) -> bool, drops points (e.g. phase near +-90 deg)
    amin: float = 0.0
    amax: float = 180.0

    def line(self):
        es = [p[0] for p in self.points]
        lo, hi = min(es) - 1e-6, max(es) + 1e-6
        f = [1, self.entrance, self.exit, f"{lo:.6f}", f"{hi:.6f}", self.amin, self.amax, self.isdiff]
        if self.isdiff == 2:
            f += [self.phase_jl[0], self.phase_jl[1]]
        f += [1, 0, 0, 0, 0, 0, f"data/{self.name}.dat", 0, 0]
        s = " ".join(str(x) for x in f)
        return s + (" " + self.tail if self.tail else "")


@dataclass
class TestSegment:
    """A <segmentsTest> line: E range (lab, MeV) and step, angle range and step, observable."""
    entrance: int
    exit: int
    emin: float
    emax: float
    estep: float
    amin: float = 0.0
    amax: float = 0.0
    astep: float = 0.0
    isdiff: int = 0
    extra: str = ""   # "J l" of a phase shift, the order of an angular distribution, keyword blocks

    def line(self):
        f = [1, self.entrance, self.exit, self.emin, self.emax, self.estep, self.amin, self.amax, self.astep,
             self.isdiff]
        s = " ".join(str(x) for x in f)
        return s + (" " + self.extra if self.extra else "")


@dataclass
class Project:
    name: str
    pairs: List[Pair]
    levels: List[Level]
    data: List[DataSegment]
    tests: List[TestSegment] = field(default_factory=list)
    potential: str = ("useHybridPotential=0\npotentialType=0\nV0=150\nR=3.6\na=0.6\n")
    config_extra: str = ""
    scattering_lengths: str = ""
    rate: Tuple[int, int] = (1, 2)    # <lastRun> entrance and exit pair of the reaction rate

    # -- free parameters, in the order AZURE2 numbers them (level energy, then its channels' widths)
    def free_parameters(self):
        out = []
        for i, lev in enumerate(self.levels, 1):
            if not lev.fixed:
                out.append((f"{lev.name} E (MeV)", ("E", i - 1, None)))
            for c, ch in enumerate(lev.channels):
                if not ch.fixed:
                    out.append((f"{lev.name} width pair {ch.pair.key} 2s={ch.s2} 2l={ch.l2}", ("W", i - 1, c)))
        return out

    def param_names(self):
        """[(param.par name, ref or None)] in the engine's order: J-pi groups by first appearance, levels in file
        order within a group, channels in the group's channel order; ref is set for the free parameters."""
        groups, order = {}, []
        for i, lev in enumerate(self.levels):
            k = (lev.J, lev.pi)
            if k not in groups:
                groups[k] = []
                order.append(k)
            groups[k].append(i)
        chans = {}
        for k in order:
            chans[k] = []
            for i in groups[k]:
                for ch in self.levels[i].channels:
                    if (ch.pair.key, ch.s2, ch.l2) not in chans[k]:
                        chans[k].append((ch.pair.key, ch.s2, ch.l2))
        out, n = [], 0
        for k in order:
            for i in groups[k]:
                n += 1
                lev = self.levels[i]
                out.append((f"energy_{n}", None if lev.fixed else ("E", i, None)))
                idx = {(c.pair.key, c.s2, c.l2): c_i for c_i, c in enumerate(lev.channels)}
                for j, key in enumerate(chans[k], 1):
                    c_i = idx.get(key)
                    free = c_i is not None and not lev.channels[c_i].fixed
                    out.append((f"width_{n}_{j}", ("W", i, c_i) if free else None))
        return out

    def get(self, ref):
        kind, i, c = ref
        return self.levels[i].E if kind == "E" else self.levels[i].channels[c].width

    def set(self, ref, value):
        kind, i, c = ref
        if kind == "E":
            self.levels[i].E = value
        else:
            self.levels[i].channels[c].width = value

    # -- writer
    def level_lines(self):
        """Every level gets a line for every channel of its J-pi group (zero and fixed where the model has none),
        so the engine has nothing to fill in and says nothing about it."""
        S_in = self.pairs[0].S
        group = {}
        for lev in self.levels:
            g = group.setdefault((lev.J, lev.pi), [])
            for ch in lev.channels:
                k = (ch.pair.key, ch.s2, ch.l2)
                if k not in [(c.pair.key, c.s2, c.l2) for c in g]:
                    g.append(ch)
        out = []
        for lid, lev in enumerate(self.levels, 1):
            have = {(c.pair.key, c.s2, c.l2): c for c in lev.channels}
            chans = [have.get((c.pair.key, c.s2, c.l2), Channel(c.pair, c.s2, c.l2, 0.0, True))
                     for c in group[(lev.J, lev.pi)]]   # in the J-group's channel order
            for ch in chans:
                p = ch.pair
                cols = [lev.J, lev.pi, _num(lev.E), int(lev.fixed), 1, p.key, ch.s2, ch.l2, lid, 1, int(ch.fixed),
                        _num(ch.width), p.j1, p.pi1, p.j2, p.pi2, _num(p.e2), _num(p.m1), _num(p.m2), p.z1, p.z2,
                        _num(S_in), _num(p.S), 0, 0, "0.0", p.ptype, _num(p.R), _num(p.g1), _num(p.g2), p.ecmask]
                out.append(" ".join(str(x) for x in cols))
            out.append("")
        return "\n".join(out)

    def old_text(self):
        """The project in the old grammar, the converter's input (write)."""
        cfg = ["true #Perform A-Matrix Calculation", "output/ #Full Path to Output Directory",
               "checks/ #Full Path to Checks Directory"] + ["none"] * 8
        txt = "<config>\n" + "\n".join(cfg) + "\n" + self.config_extra + "</config>\n"
        txt += "<potential>\n" + self.potential + "</potential>\n"
        txt += "<levels>\n" + self.level_lines() + "\n</levels>\n"
        txt += "<segmentsData>\n" + "".join(d.line() + "\n" for d in self.data) + "</segmentsData>\n"
        txt += "<segmentsTest>\n" + "".join(t.line() + "\n" for t in self.tests) + "</segmentsTest>\n"
        if self.scattering_lengths:
            txt += "<scatteringLengths>\n" + self.scattering_lengths + "</scatteringLengths>\n"
        txt += "<targetInt>\n</targetInt>\n"
        txt += ("<lastRun>\n9997\n\"\"\n\"\"\n%d %d\n1 \"temperatures.dat\"\n-1. -1. -1.\n</lastRun>\n"
                % self.rate)
        return txt

    def write(self, directory, binary, azr_name=None, data=None):
        """Write <directory>/<name>.azr (format 3) and, if given, first the data files {segment name: rows}.

        The old text goes to a temporary file in directory, which the converter (binary) reads with the data files
        beside it and removes nothing else; the conversion is <name>.azr."""
        os.makedirs(os.path.join(directory, "data"), exist_ok=True)
        if data is not None:
            for seg in self.data:
                write_data(os.path.join(directory, "data", seg.name + ".dat"), data[seg.name], seg)
        name = (azr_name or self.name) + ".azr"
        old = f".{name}.old"
        with open(os.path.join(directory, old), "w") as f:
            f.write(self.old_text())
        try:
            if os.path.exists(os.path.join(directory, name)):
                os.remove(os.path.join(directory, name))
            p = subprocess.run([binary, "--convert", old, name, "--to", "3", "--exact-grids"], cwd=directory,
                               capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1", LC_ALL="C"))
            if p.returncode != 0:
                raise RuntimeError(f"AZURE2 --convert of {name} failed:\n" + (p.stdout + p.stderr)[-3000:])
        finally:
            os.remove(os.path.join(directory, old))
        return os.path.join(directory, name)


def _num(x):
    if isinstance(x, (int, np.integer)):
        return str(x)
    return repr(float(x)) if abs(x) < 1e15 else "%.10e" % x


def write_data(path, rows, seg, header=""):
    with open(path, "w") as f:
        if header:
            f.write(header)
        for e, a, v, dv in rows:
            f.write(f"{e:.6f} {a:.4f} {v:.10e} {dv:.6e}\n")


# ----------------------------------------------------------------------------- engine


def project_copy(name, work):
    """Copy tests/upgrade/<name> (the .azr files, data/, temperatures.dat) into work, with output/ and checks/;
    returns the path of <name>.azr there.  Tests run on such copies and never write into the source tree."""
    src = os.path.join(REPO, "tests", "upgrade", name)
    for f in os.listdir(src):
        if f.endswith(".azr") or f == "temperatures.dat":
            shutil.copy(os.path.join(src, f), work)
    shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    os.makedirs(os.path.join(work, "checks"), exist_ok=True)
    return os.path.join(work, name + ".azr")


def find_binary(arg=None):
    """The AZURE2 CLI: the argument, $AZURE2_BIN, $AZURE2_BINARY, else the newest build*/src/AZURE2."""
    for cand in (arg, os.environ.get("AZURE2_BIN"), os.environ.get("AZURE2_BINARY")):
        if cand:
            if not os.path.isfile(cand):
                raise SystemExit(f"AZURE2 binary {cand!r} does not exist")
            return os.path.abspath(cand)
    found = sorted(glob.glob(os.path.join(REPO, "build*", "src", "AZURE2")), key=os.path.getmtime)
    if not found:
        raise SystemExit("no AZURE2 binary: pass it, or set AZURE2_BIN (or AZURE2_BINARY)")
    return found[-1]


def run(binary, workdir, azr, command, flags=(), threads=2, timeout=1800):
    """Run the CLI in workdir by `command` (a list: the command and its options); returns the log text."""
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    os.makedirs(os.path.join(workdir, "checks"), exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), LC_ALL="C")
    p = subprocess.run([binary, *command, "--no-gui", *flags, os.path.basename(azr)], cwd=workdir,
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, env=env, timeout=timeout)
    log = p.stdout + p.stderr
    with open(os.path.join(workdir, "run.log"), "w") as f:
        f.write(log)
    if p.returncode != 0:
        raise RuntimeError(f"AZURE2 exited with {p.returncode} in {workdir}:\n" + log[-3000:])
    return log


CALC = ["calculate"]           # calculate segments from data
FIT = ["fit"]                  # fit segments from data (Minuit), no band
EXTRAP = ["extrapolate"]       # calculate segments without data, no band


def rate_command(entrance, exit_key):
    """The rate of the pairs entrance -> exit_key at the temperatures of temperatures.dat."""
    return ["rate", f"--rate-entrance={entrance}", f"--rate-exit={exit_key}", "--temperature-file=temperatures.dat"]


def read_blocks(path):
    """An AZUREOut file: list of segments, each an array of rows (blank-line separated)."""
    blocks, cur = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                if cur:
                    blocks.append(np.array(cur, float))
                    cur = []
                continue
            try:
                cur.append([float(x) for x in s.split()])
            except ValueError:
                continue
    if cur:
        blocks.append(np.array(cur, float))
    return blocks


def segment_outputs(outdir, segments):
    """The AZUREOut rows of each data segment, in segment order (one file per entrance/exit pair)."""
    files = {}
    out = []
    for seg in segments:
        key = (seg.entrance, seg.exit)
        if key not in files:
            files[key] = read_blocks(os.path.join(outdir, f"AZUREOut_aa={seg.entrance}_R={seg.exit}.out"))
        out.append(files[key].pop(0))
    return out


def read_chi2(path):
    seg = {}
    total = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("Total-Chi-Squared:"):
                total = float(s.split()[1])
            elif s and s[0].isdigit():
                p = s.split(",")
                seg[int(p[0])] = (float(p[1]), int(p[2]))
    return total, seg


# ----------------------------------------------------------------------------- synthetic data


def truth_values(binary, project, threads=2):
    """The engine's observable at every data point, in the data file's frame and units.

    The data files are first written with value 1 and uncertainty 1; the engine converts them to the c.m. frame,
    so its data column (6) is the lab-to-c.m. factor of each point and the model in the data frame is column 4
    divided by it.  Points are checked to come back in file order (c.m. energy proportional to lab energy).
    """
    work = tempfile.mkdtemp(prefix=f"{project.name}_truth_")
    try:
        ones = {s.name: [(e, a, 1.0, 1.0) for e, a in s.points] for s in project.data}
        azr = project.write(work, binary, data=ones)
        with open(os.path.join(work, "temperatures.dat"), "w") as f:
            f.write("0.1\n")
        run(binary, work, azr, CALC, threads=threads)
        rows = segment_outputs(os.path.join(work, "output"), project.data)
        vals = {}
        for seg, r in zip(project.data, rows):
            if len(r) != len(seg.points):
                raise RuntimeError(f"{seg.name}: {len(r)} output rows for {len(seg.points)} points")
            elab = np.array([p[0] for p in seg.points])
            ratio = r[:, 0] / elab
            if np.ptp(ratio) > 1e-3 * abs(ratio.mean()):
                raise RuntimeError(f"{seg.name}: output rows are not in file order")
            vals[seg.name] = r[:, 3] / r[:, 5]
        return vals
    finally:
        shutil.rmtree(work, ignore_errors=True)


def synthesize(project, truth, seed):
    """Gaussian noise on the truth: ("rel", f) -> sigma = f |truth|; ("abs", s) -> sigma = s."""
    rng = np.random.default_rng(seed)
    data = {}
    for seg in project.data:
        t = truth[seg.name]
        kind, size = seg.noise
        sig = size * np.abs(t) if kind == "rel" else np.full_like(t, size)
        z = rng.normal(0.0, 1.0, t.size)
        data[seg.name] = [(e, a, v + s * zz, s) for (e, a), v, s, zz in zip(seg.points, t, sig, z)]
    return data


def grid(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [round(lo + i * step, 6) for i in range(n + 1)]


def read_param_file(path):
    """param.par / param.sav: {name: (value, error)}; the error of a fitted parameter is Minuit's (param.sav)."""
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3 and not p[0].startswith("#"):
                out[p[0]] = (float(p[1]), float(p[2]))
    return out
