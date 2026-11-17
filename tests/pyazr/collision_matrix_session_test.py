#!/usr/bin/env python3
"""pyazr.physics.collision_matrix against a session's own numbers on the test projects T1, T2 and T3, in the three
parametrizations of the levels (a session with its defaults: Brune's; use_brune=False: formal; use_park=True:
Park's).

What is compared.  The arrays of every J^pi block are read off the session and handed to
``pyazr.physics.collision_matrix``; an observable is formed from its U or T in numpy and compared with the model
value the session computes at the same energy (``Session.evaluate``):

- the levels: the session's parameter records at the start vector (``describe()["parameters"]``: ``start`` is the
  level energy in MeV or the reduced-width amplitude in the session's parametrization);
- the channel functions at the energy and S_c at each level energy: ``Session.channel_functions`` (P, S, the
  hard-sphere phase) and the rows ``coulomb_phase``, ``photon_penetrability`` of ``pyazr.physics``;
- the boundary values of the formal levels: B_c = S_c at the lowest level energy of the block
  (``docs/upgrade/level_parameters.md`` Eq. (3.1));
- Park's norms: ``Session.park_norms``.

The observables (``physics/Observables.h``; k from the row ``wave_number``, g = the pair's ``i1i2factor``):

    phase shift       delta = arg(U_cc exp(-2i omega_c)) / 2  (deg; in [0, 180) for identical particles)
    sigma(a -> a')    = (pi / k^2) g sum_J (2J+1) sum_{c in a, c' in a'} |T_cc'|^2
    sigma_total(a)    = (pi / k^2) g sum_J (2J+1) sum_{c in a} 2 Re T_cc          (a neutral pair)

T1: the s-wave phase shift of p + 12C (the 1/2+ block, which has two photon channels), and the two capture cross
sections of a copy of T1 whose photon pairs take no external capture (the multipole mask of their rows set to 0: the
capture is then the R-matrix's photon channels alone, which is what the function computes).  T2: the ratio
sigma(p, alpha_1) / sigma(p, alpha_0) of p + 15N (three coupled 1- levels, two 3- levels, three particle pairs).
T3: the d-wave phase shift of alpha + alpha (identical particles, the hybrid potential's channel functions, a
beta-decay channel in the block), the total and the (n,p) cross section of n + 7Be, and 7Li(p,n).

The function and the session share physics/CollisionMatrix.h's arithmetic of one pathway; what this test adds to
tests/pyazr/collision_matrix_test.py (numpy, no session) is that the arrays mean what the docstring says: a
parametrization, a boundary value, a shift at a level energy or a photon factor read differently from the engine's
moves these numbers.  The agreement measured is printed per comparison.

Run from anywhere:  python3 tests/pyazr/collision_matrix_session_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
UPGRADE = os.path.join(ROOT, "tests", "upgrade")
sys.path.insert(0, ROOT)

failures = []

# The largest relative difference (for a phase shift: the largest difference in degrees over 1 degree) accepted
# between the function and the session.  Both run the same arithmetic on the same channel functions; what is left
# is the rounding of the sums over channels taken in another order.
TOLERANCE = 1e-9

PARAMETRIZATIONS = [("brune", {}), ("formal", {"use_brune": False}), ("park", {"use_park": True})]


def check(what, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(what)


class Model:
    """The blocks of a session at its start vector, and their collision matrices at an excitation energy."""

    def __init__(self, session, parametrization):
        import numpy as np
        self.np = np
        self.s = session
        self.parametrization = parametrization
        described = session.describe()
        self.pairs = {p["key"]: p for p in described["pairs"]}
        records = [r for r in described["parameters"] if r["kind"] in ("energy", "width")]
        norms = np.asarray(session.park_norms(), float)  # one per energy row, in their order
        energy_rows = [r for r in records if r["kind"] == "energy"]
        self.blocks = []
        for group in sorted({r["jgroup"] for r in records}):
            rows = [r for r in records if r["jgroup"] == group]
            levels = [r for r in rows if r["kind"] == "energy"]
            widths = [r for r in rows if r["kind"] == "width"]
            channels = {}
            for r in widths:
                channels.setdefault(r["channel"], dict(pair=r["pair_key"], l=r["L"], s=r["S"],
                                                       photon=self.pairs[r["pair_key"]]["pair_kind"] != "particle"))
            numbers = sorted(channels)
            gamma = np.zeros((len(levels), len(numbers)))
            for r in widths:
                gamma[r["level"] - 1, numbers.index(r["channel"])] = r["start"]
            self.blocks.append(dict(
                J=levels[0]["J"], parity=levels[0]["parity"], channels=[channels[n] for n in numbers],
                E=np.array([r["start"] for r in levels]), gamma=gamma,
                norms=np.array([norms[energy_rows.index(r)] for r in levels])))
        for block in self.blocks:  # S_c at each level energy; B_c = S_c at the lowest level
            block["S_levels"] = np.array([[self.functions(c, e)["S"] for c in block["channels"]] for e in block["E"]])
            block["B"] = block["S_levels"][int(np.argmin(block["E"]))]

    def channel_energy(self, channel, e_x):
        pair = self.pairs[channel["pair"]]
        return e_x - pair["sep_energy"] - pair["excitation"] - pair["light_excitation"]

    def functions(self, channel, e_x):
        """P, S, phi (rad) and omega (rad) of one channel at the excitation energy e_x (MeV)."""
        from pyazr import physics as ph
        pair, e = self.pairs[channel["pair"]], self.channel_energy(channel, e_x)
        zero = dict(P=0.0, S=0.0, phi=0.0, omega=0.0)
        if pair["pair_kind"] == "photon":  # the engine's factor (E_gamma / hbar c)^(2L+1); nothing else
            return dict(zero, P=float(ph.photon_penetrability(e, channel["l"])) if e > 0 else 0.0)
        if pair["pair_kind"] != "particle":  # a beta-decay channel: outside the level matrix, not compared
            return zero
        f = self.s.channel_functions(channel["pair"], e, channel["l"])
        if not f["open"]:
            return dict(zero, S=float(f["S"]))
        omega = 0.0
        if pair["Z1"] * pair["Z2"] != 0:
            eta = ph.sommerfeld_eta(pair["Z1"], pair["Z2"], pair["M1"], pair["M2"], e)
            omega = float(ph.coulomb_phase(channel["l"], eta) - ph.coulomb_phase(0, eta))
        return dict(P=float(f["P"]), S=float(f["S"]), phi=-float(f["delta_hs"]), omega=omega)

    def matrices(self, block, e_x):
        """The dict of collision_matrix for one block at e_x, and omega per channel."""
        from pyazr import physics as ph
        np = self.np
        f = [self.functions(c, e_x) for c in block["channels"]]
        photon = np.array([c["photon"] for c in block["channels"]])
        arguments = dict(penetrability=[x["P"] for x in f], shift=[x["S"] for x in f],
                         hard_sphere_phase=[x["phi"] for x in f], coulomb_phase=[x["omega"] for x in f],
                         photon_channels=photon)
        if self.parametrization == "formal":
            arguments["boundary"] = np.where(photon, 0.0, block["B"])
        else:
            arguments.update(parametrization=self.parametrization, shift_at_levels=block["S_levels"])
            if self.parametrization == "park":
                arguments["park_norms"] = block["norms"]
        return ph.collision_matrix(block["E"], block["gamma"], e_x, **arguments), [x["omega"] for x in f]

    def phase_shift(self, e_x, J, parity, pair, l, s, identical=False):
        """The phase shift (deg) of the channel (pair, l, s) of the block J^parity."""
        np = self.np
        block = next(b for b in self.blocks if b["J"] == J and b["parity"] == parity)
        c = next(i for i, ch in enumerate(block["channels"]) if (ch["pair"], ch["l"], ch["s"]) == (pair, l, s))
        out, omega = self.matrices(block, e_x)
        u = out["U"][c, c] * np.exp(-2j * omega[c])
        delta = np.degrees(np.arctan2(u.imag, u.real)) / 2.0
        return delta + 180.0 if identical and delta < 0.0 else delta

    def sigma(self, e_x, entrance, exit=None):
        """sigma(entrance -> exit) in barns; the total cross section of the entrance pair for exit None."""
        from pyazr import physics as ph
        np = self.np
        pair = self.pairs[entrance]
        e = e_x - pair["sep_energy"] - pair["excitation"] - pair["light_excitation"]
        k = float(ph.wave_number(pair["M1"], pair["M2"], e))
        total = 0.0
        for block in self.blocks:
            ins = [i for i, ch in enumerate(block["channels"]) if ch["pair"] == entrance]
            outs = [i for i, ch in enumerate(block["channels"]) if ch["pair"] == exit]
            if not ins or (exit is not None and not outs):
                continue
            T = self.matrices(block, e_x)[0]["T"]
            if exit is None:
                total += (2.0 * block["J"] + 1.0) * sum(2.0 * T[c, c].real for c in ins)
            else:
                total += (2.0 * block["J"] + 1.0) * sum(abs(T[c, cp]) ** 2 for c in ins for cp in outs)
        return np.pi / k ** 2 * pair["i1i2factor"] * total / 100.0  # fm^2 -> b


def without_external_capture(src, dest):
    """A copy of a project whose photon pairs take no external capture: the multipole mask (the last column of a
    <levels> row, io/LevelRows.h) of every photon row (pair kind 10, the fifth column from the end) set to 0.
    Returns the number of rows changed."""
    changed, inside, lines = 0, False, []
    with open(src) as f:
        for line in f:
            if line.strip() in ("<levels>", "</levels>"):
                inside = line.strip() == "<levels>"
            parts = line.split()
            if inside and len(parts) == 31 and parts[26] == "10" and parts[30] != "0":
                parts[30] = "0"
                line = " ".join(parts) + "\n"
                changed += 1
            lines.append(line)
    with open(dest, "w") as f:
        f.writelines(lines)
    return changed


def compare(what, got, want, absolute=False):
    """One comparison over the points of a data set; prints the measured agreement."""
    import numpy as np
    got, want = np.asarray(got, float), np.asarray(want, float)
    d = np.abs(got - want) if absolute else np.abs(got - want) / np.abs(want)
    worst = float(d.max())
    unit = "deg" if absolute else "relative"
    digits = "all" if worst == 0.0 else f"{-np.log10(worst / (1.0 if not absolute else np.abs(want).max())):.1f}"
    check(f"{what}: {got.size} points, largest difference {worst:.1e} ({unit}; {digits} digits)", worst < TOLERANCE,
          f"limit {TOLERANCE:.0e}")
    return worst


def main():
    try:
        import numpy as np
        import pyazr
        from pyazr import physics as ph  # noqa: F401
    except Exception as err:  # the skip of the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 77
    warnings.simplefilter("ignore", pyazr.EngineWarning)

    work = tempfile.mkdtemp(prefix="collision_matrix_session_")
    try:
        for name in ("T1", "T2", "T3"):
            shutil.copytree(os.path.join(UPGRADE, name), os.path.join(work, name))
        t1_internal = os.path.join(work, "T1", "T1_internal.azr")
        changed = without_external_capture(os.path.join(work, "T1", "T1.azr"), t1_internal)
        check("the copy of T1 without external capture: photon rows changed", changed > 0, f"{changed} rows")

        capture_with_ec = {}
        for parametrization, options in PARAMETRIZATIONS:
            print(f"T1, {parametrization}")
            with pyazr.Session(os.path.join(work, "T1", "T1.azr"), **options) as s:
                model, ev = Model(s, parametrization), s.evaluate(s.parameters.x0)
                data = {d.name: d for d in ev.data}
                d = data["phase_s12"]
                compare("phase shift of p + 12C, 1/2+ l = 0 (two photon channels in the block)",
                        [model.phase_shift(e, 0.5, 1, 1, 0, 0.5) for e in d.excitation_energies], d.model,
                        absolute=True)
                capture_with_ec = {n: np.array(data[n].model) for n in ("cap_g0", "cap_g1")}
            with pyazr.Session(t1_internal, **options) as s:
                model, ev = Model(s, parametrization), s.evaluate(s.parameters.x0)
                data = {d.name: d for d in ev.data}
                for dataset, exit_pair in (("cap_g0", 2), ("cap_g1", 3)):
                    d = data[dataset]
                    moved = np.abs(np.array(d.model) / capture_with_ec[dataset] - 1.0).max()
                    check(f"{dataset}: the copy without external capture differs from T1 (the edit bites)",
                          moved > 1e-3, f"largest relative change {moved:.1e}")
                    compare(f"sigma(p, gamma) to pair {exit_pair} without external capture (photon channels)",
                            [model.sigma(e, 1, exit_pair) for e in d.excitation_energies], d.model)

            print(f"T2, {parametrization}")
            with pyazr.Session(os.path.join(work, "T2", "T2.azr"), **options) as s:
                model, ev = Model(s, parametrization), s.evaluate(s.parameters.x0)
                d = {x.name: x for x in ev.data}["s09_ratio"]
                compare("sigma(p, alpha_1) / sigma(p, alpha_0) of p + 15N",
                        [model.sigma(e, 1, 3) / model.sigma(e, 1, 2) for e in d.excitation_energies], d.model)

            print(f"T3, {parametrization}")
            with pyazr.Session(os.path.join(work, "T3", "T3.azr"), **options) as s:
                model, ev = Model(s, parametrization), s.evaluate(s.parameters.x0)
                data = {x.name: x for x in ev.data}
                d = data["aa_phase_d"]
                compare("phase shift of alpha + alpha, 2+ l = 2 (identical particles, hybrid potential)",
                        [model.phase_shift(e, 2.0, 1, 1, 2, 0.0, identical=True) for e in d.excitation_energies],
                        d.model, absolute=True)
                d = data["n7be_total"]
                compare("total cross section of n + 7Be (optical theorem)",
                        [model.sigma(e, 3) for e in d.excitation_energies], d.model)
                d = data["n7be_np"]
                compare("sigma of 7Be(n,p)7Li", [model.sigma(e, 3, 2) for e in d.excitation_energies], d.model)
                d = data["p7li_pn"]
                compare("sigma of 7Li(p,n)7Be", [model.sigma(e, 2, 3) for e in d.excitation_energies], d.model)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
