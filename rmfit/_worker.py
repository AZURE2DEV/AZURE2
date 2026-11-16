"""One engine session in one process, driven over a pipe.

Every worker is a separate process holding exactly one ``pyazr.Session``.  A hard abort inside the engine (an
error the API layer does not catch) kills the process; the parent notices the broken pipe, restarts the worker
and re-issues the task once.

The numbers come from the session's two objects (pyazr 3.0): ``Session.objective()`` for the least-squares rows
and their Jacobian, ``Session.evaluate(x)`` for chi-squared and the model values.  The objective's rows are the
data rows first, then one row per prior and constraint of the engine; rmfit keeps the data rows only and appends
its own penalty rows in the parent (:mod:`rmfit.engine`), once for all shards.
"""
import os
import sys
import traceback


def worker_main(conn, azr, cwd, flags, seg_map, omp_threads=1):
    os.environ["OMP_NUM_THREADS"] = str(omp_threads)
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here not in sys.path:
        sys.path.insert(0, here)
    import numpy as np
    from rmfit import AZURE2_ROOT  # noqa: F401
    from pyazr import Session
    from rmfit.azr3 import nuisance_rows
    from rmfit.keys import build_keys, jpi_string, ParamKey, level_label, jgroup_nucleus_tags

    try:
        with Session(azr, cwd=cwd, **flags) as m:
            keys = build_keys(m, seg_map=seg_map)
            # level labels gain "@<nucleus>" in a multi-nucleus model (keys.level_label; 2026-10-30)
            tags = jgroup_nucleus_tags(m)
            active = m.active_datasets
            P = m.parameters
            f = m.objective()
            x0 = np.asarray(P.x0, float)
            # Nothing is evaluated here: the engine freezes the model of a ``relvector`` covariance term at the
            # first evaluation of the session, which must be the caller's first vector, not the start vector.
            npts = [int(len(e)) for e in m.energies]
            #: the objective's data rows: [0, ndata); the engine's priors and constraints follow and are dropped
            ndata = int(sum(npts))
            info = dict(
                keys=keys,
                x0=x0,
                npts=npts,
                seg_keys=[int(seg_map[d.key]) if seg_map else int(d.key) for d in active],
                n_rmatrix=len(P.free.of_kind("energy", "width")),
                wigner=[(p.free_index, float(p.wigner_limit))
                        for p in P
                        if p.kind == "width" and not p.fixed and p.free_index is not None
                        and p.wigner_limit is not None],
                # fixed, nonzero amplitudes: not in the free vector, but they
                # pin the sign gauge (see keys.canonical_signs)
                fixed_amplitudes=[
                    ParamKey("width", level=level_label(p.J, p.parity, p.level, p.jgroup, tags),
                             pair=int(p.pair_key), L=int(p.L), S=float(p.S if p.S is not None else 0.0))
                    for p in P
                    if p.kind == "width" and p.fixed and p.start != 0.0],
                level_energy={level_label(p.J, p.parity, p.level, p.jgroup, tags):
                              (0.0 if p.level_energy is None else float(p.level_energy))
                              for p in P.of_kind("energy")},
                # Park mode: the level of every entry of park_norms() (the order of the energy rows), for the
                # J > 0 wall rows (engine.ParkWall)
                park_levels=([level_label(p.J, p.parity, p.level, p.jgroup, tags) for p in P.of_kind("energy")]
                             if flags.get("use_park") else []),
                # measured scattering lengths (engine.ScatLenRows, 2026-11-07): 1/((2j1+1)(2j2+1)) per pair key,
                # the weight of the coherent combination
                pair_i1i2={int(q.key): float(q.i1i2factor) for q in m.pairs
                           if getattr(q, "i1i2factor", None) is not None},
                # the nuisances and the priors of the engine's objective, by parameter row (azr3.nuisance_rows)
                nuisances=nuisance_rows(P),
            )

            #: the nuisance priors are the last rows of the objective, one per free parameter with a prior
            n_prior = sum(1 for p in P if not p.fixed and (getattr(p, "prior_sigma", None) or 0.0) > 0)
            prior_rows = [(int(p.free_index), float(p.prior_centre), float(p.prior_sigma)) for p in P
                          if not p.fixed and (getattr(p, "prior_sigma", None) or 0.0) > 0]

            def engine_chi2(x):
                """The engine's objective without the priors of the nuisances: the data chi-squared
                (``Evaluation.chi2``, which holds the data covariances) with the model's own terms (the rows of
                engine.ModelRows: the J > 0 wall, abundance sums, model data, ...).  The normalization and
                energy-shift penalties are the parent's rows (engine.Penalties)."""
                ev = m.evaluate(x)
                return float(ev.objective) - sum(((x[j] - c) / s) ** 2 for j, c, s in prior_rows)

            def model_rows(x, want_jac):
                """The objective's rows between the data rows and the nuisance priors (engine.ModelRows).  The
                tail is checked to be the nuisance priors (sum of squares), so that a change of the engine's row
                layout cannot pass silently."""
                if want_jac:
                    r, J = f.residuals_and_jacobian(x)
                else:
                    r, J = f.residuals(x), None
                r = np.asarray(r, float)
                end = r.size - n_prior
                if end < ndata:
                    raise RuntimeError(f"objective rows: {r.size} for {ndata} data rows and {n_prior} priors")
                want = sum(((x[j] - c) / s) ** 2 for j, c, s in prior_rows)
                got = float(np.sum(r[end:] ** 2))
                if abs(got - want) > 1e-8 * max(1.0, abs(want)):
                    raise RuntimeError(f"objective rows: the last {n_prior} rows are not the nuisance priors "
                                       f"({got:.10g} vs {want:.10g})")
                return r[ndata:end], (None if J is None else np.asarray(J, float)[ndata:end])

            conn.send(("ready", info))
            while True:
                cmd, payload = conn.recv()
                try:
                    if cmd == "quit":
                        conn.send(("ok", None))
                        break
                    elif cmd == "eval":
                        conn.send(("ok", engine_chi2(np.asarray(payload, float))))
                    elif cmd == "resid":
                        x, want_jac = payload[0], payload[1]
                        skip = payload[2] if len(payload) > 2 else None
                        x = np.asarray(x, float)
                        if want_jac:
                            r, J = f.residuals_and_jacobian(x, skip=(list(skip) if skip else None))
                            conn.send(("ok", (np.asarray(r, float)[:ndata], np.asarray(J, float)[:ndata])))
                        else:
                            conn.send(("ok", (np.asarray(f.residuals(x), float)[:ndata], None)))
                    elif cmd == "model_rows":
                        x, want_jac = payload
                        conn.send(("ok", model_rows(np.asarray(x, float), bool(want_jac))))
                    elif cmd == "fast_status":
                        # pyazr 3.0: the rows without a Jacobian are the objective's own (one forward evaluation)
                        conn.send(("ok", dict(ok=True, checked=True)))
                    elif cmd == "transform":
                        # the physical values of the free parameters, in the order of the fit vector
                        x = np.asarray(payload, float)
                        conn.send(("ok", np.asarray(m.transform(x), float)[np.asarray(P.free_mask, bool)]))
                    elif cmd == "transform_all":
                        # 2026-10-24: physical values of EVERY parameter (fixed ones included, parameter order) for the
                        # free vector ``payload``, as save_fit writes them.  A fixed parameter is fixed in the fit's
                        # basis; its observed value moves with the other amplitudes of its level.
                        x = np.asarray(payload, float)
                        if x.size != len(x0):
                            raise ValueError(f"transform_all: {x.size} values for {len(x0)} free parameters")
                        conn.send(("ok", np.asarray(m.transform(x), float)))
                    elif cmd == "cross":
                        x = np.asarray(payload, float)
                        conn.send(("ok", [np.asarray(d.model, float) for d in m.evaluate(x).data]))
                    elif cmd == "data":
                        conn.send(("ok", dict(
                            energies=[np.asarray(e, float) for e in m.energies],
                            excitation=[np.asarray(e, float) for e in m.excitation_energies],
                            cross=[np.asarray(c, float) for c in m.cross],
                            cross_err=[np.asarray(c, float) for c in m.cross_err],
                            # 2026-09-28: per-point data angles (plots split multi-angle segments by detector)
                            angles=[np.asarray(a, float) for a in m.angles])))
                    elif cmd == "write_output":
                        m.write_outputs(np.asarray(payload, float))
                        conn.send(("ok", None))
                    elif cmd == "park_norms":
                        conn.send(("ok", np.asarray(m.park_norms(np.asarray(payload, float)), float)))
                    elif cmd == "park_jac":
                        # dJ/dx by central differences on the given free columns, (len(cols), nlevels).
                        # J is exactly quadratic in a width amplitude, so the difference is exact there;
                        # an energy column carries the engine's numerical d2S/dE2 (step 1e-6 MeV).
                        x, cols = payload
                        x = np.asarray(x, float)
                        D = []
                        for c in cols:
                            h = 1e-6 if keys[c].kind == "energy" else max(1e-6 * abs(x[c]), 1e-9)
                            xp, xm = x.copy(), x.copy()
                            xp[c] += h
                            xm[c] -= h
                            D.append((np.asarray(m.park_norms(xp), float) - np.asarray(m.park_norms(xm), float))
                                     / (2.0 * h))
                        conn.send(("ok", np.asarray(D, float).reshape(len(cols), -1)))
                    elif cmd == "scatlen":
                        # every pair's s-wave scattering lengths at x (pyazr azure2.scattering_lengths), one row
                        # per (pair, channel spin): pair, s, parity, Re a0, Im a0, has_levels (engine.ScatLenRows)
                        v = m.scattering_lengths(np.asarray(payload, float))
                        conn.send(("ok", np.asarray([[q.pair, q.s, q.parity, q.a0.real, q.a0.imag,
                                                      float(q.has_levels)] for q in v], float).reshape(-1, 6)))
                    elif cmd == "scatlen_jac":
                        # d(Re a0, Im a0)/dx by central differences on the given free columns,
                        # (len(cols), n_channels, 2); steps as park_jac (energy 1e-6 MeV, amplitude 1e-6 relative)
                        x, cols = payload
                        x = np.asarray(x, float)

                        def _a(xx):
                            return np.asarray([[q.a0.real, q.a0.imag] for q in m.scattering_lengths(xx)], float)
                        D = []
                        for c in cols:
                            h = 1e-6 if keys[c].kind == "energy" else max(1e-6 * abs(x[c]), 1e-9)
                            xp, xm = x.copy(), x.copy()
                            xp[c] += h
                            xm[c] -= h
                            D.append((_a(xp) - _a(xm)) / (2.0 * h))
                        conn.send(("ok", np.asarray(D, float).reshape(len(cols), -1, 2)))
                    elif cmd == "wigner_widths":
                        conn.send(("ok", dict(m.wigner_widths(np.asarray(payload, float)))))
                    elif cmd == "snapshot":
                        # the parameter table as plain data (name, kind, fixed, start, the owner fields) with the
                        # start of every parameter in the fit's basis
                        conn.send(("ok", dict(
                            parameters=[dict(name=p.name, kind=p.kind, fixed=bool(p.fixed), free_index=p.free_index,
                                             start=float(p.start), value=p.value, J=p.J, parity=p.parity,
                                             level=p.level, jgroup=p.jgroup, level_label=p.level_label,
                                             pair_key=p.pair_key, L=p.L, S=p.S, segment_key=p.segment_key,
                                             fn_term=p.fn_term) for p in P],
                            pairs=[q.to_dict() for q in m.pairs],
                            names=[p.name for p in P],
                            fixed=[bool(p.fixed) for p in P],
                            basis=str(m.basis))))
                    elif cmd == "save_fit":
                        # a snapshot of the fit as a project and its parameter file (Session.save_fit: the physical
                        # values into the rows they come from, checked by reopening what it wrote)
                        path, x, kw = payload
                        conn.send(("ok", tuple(m.save_fit(path, np.asarray(x, float), **(kw or {})))))
                    else:
                        conn.send(("err", f"unknown command {cmd!r}"))
                except Exception:
                    conn.send(("err", traceback.format_exc()))
    except Exception:
        try:
            conn.send(("fatal", traceback.format_exc()))
        except Exception:
            pass
