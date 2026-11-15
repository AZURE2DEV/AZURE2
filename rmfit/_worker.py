"""One engine session in one process, driven over a pipe.

The engine is not thread-safe and does not scale with OpenMP threads, so every
worker is a separate single-threaded process holding exactly one
``pyazr.azure2`` session.  A hard abort inside the engine (a GSL error the API
layer does not catch) kills the process; the parent notices the broken pipe,
restarts the worker and re-issues the task once.
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
    from pyazr import azure2
    from rmfit.keys import build_keys, jpi_string, ParamKey, level_label, jgroup_nucleus_tags

    try:
        with azure2(azr, cwd=cwd, **flags) as m:
            keys = build_keys(m, seg_map=seg_map)
            # level labels gain "@<nucleus>" in a multi-nucleus model (keys.level_label; 2026-10-30)
            tags = jgroup_nucleus_tags(m)
            active = m.active_datasets
            info = dict(
                keys=keys,
                x0=np.asarray(m.params_rwa, float),
                npts=[int(len(e)) for e in m.energies],
                seg_keys=[int(seg_map[d.key]) if seg_map else int(d.key) for d in active],
                n_rmatrix=int(m.n_rmatrix),
                wigner=[(p.free_index, float(p.wigner_limit))
                        for p in m.parameters
                        if p.kind == "width" and not p.fixed and p.free_index is not None
                        and p.wigner_limit is not None],
                # fixed, nonzero amplitudes: not in the free vector, but they
                # pin the sign gauge (see keys.canonical_signs)
                fixed_amplitudes=[
                    ParamKey("width", level=level_label(p.J, p.parity, p.level, p.jgroup, tags),
                             pair=int({q.number: q.key for q in m.pairs}.get(p.pair, p.pair)),
                             L=int(p.L), S=float(p.S if p.S is not None else 0.0))
                    for p in m.parameters
                    if p.kind == "width" and p.fixed and p.value != 0.0],
                level_energy={level_label(k.J, k.parity, k.level, k.jgroup, tags):
                              (0.0 if k.energy is None else float(k.energy))
                              for k in m.parameters.by_physical_level()},
                # Park mode: the level of every entry of park_norms() (physical_levels() order), for the
                # J > 0 wall rows (engine.ParkWall)
                park_levels=([level_label(k.J, k.parity, k.level, k.jgroup, tags) for k in m.physical_levels()]
                             if flags.get("use_park") else []),
                # measured scattering lengths (engine.ScatLenRows, 2026-11-07): 1/((2j1+1)(2j2+1)) per pair key,
                # the weight of the coherent combination (pyazr azure2.scattering_length_chi2)
                pair_i1i2={int(q.key): float(q.i1i2factor) for q in m.pairs
                           if getattr(q, "i1i2factor", None) is not None},
            )
            # residuals without the Jacobian: (fit - data*norm) / (err*norm)
            # from one forward evaluation, instead of residual_jacobian's
            # adjoint pass (~2-3 evaluations).  Checked against the engine's
            # own residuals at start-up; falls back if they ever disagree.
            data = [np.asarray(c, float) for c in m.cross]
            err = [np.asarray(c, float) for c in m.cross_err]
            norm_free = {p.segment_key: p.free_index for p in m.parameters.norms
                         if not p.fixed and p.free_index is not None}
            seg_norm = [(norm_free.get(d.key), float(d.norm)) for d in active]
            # a normgroup follower takes its group leader's norm (EData::ResolveNormGroups); shard
            # files carry no normgroup, so this matters for single sessions only
            from rmfit.covariance import NormGroups
            groups = NormGroups.from_azr(azr)
            if groups:
                by_key = {d.key: d for d in active}
                seg_norm = [(norm_free.get(groups.owner(d.key)), float(by_key[groups.owner(d.key)].norm))
                            if d.key in groups.leader else sn for d, sn in zip(active, seg_norm)]
            x0 = np.asarray(m.params_rwa, float)

            # 2026-10-22: the engine converts a segment's data to the c.m. frame at its
            # SHIFTED energy, so with an energy shift the data and errors read at start-up are stale (exact only
            # where the lab -> c.m. factor does not depend on energy, e.g. elastic scattering).  Segments with any
            # shift term are re-read from the engine after each evaluation (update_data copies its current points).
            from rmfit.covariance import split_keyword_tail, segment_lines
            _lines = segment_lines(azr)
            def _has_shift(d):
                if d.vary_shift or d.energy_shift != 0.0 or getattr(d, "vary_shift_sqrt", False) \
                        or getattr(d, "energy_shift_sqrt", 0.0) != 0.0 or getattr(d, "shift_fn_params", ()):
                    return True
                t = _lines[d.key - 1] if 0 < d.key <= len(_lines) else []
                return any(kw in ("shiftgroup", "shiftfn", "shiftfngroup", "sqrtshift") for kw, _ in split_keyword_tail(t)[1])
            shifted = [i for i, d in enumerate(active) if _has_shift(d)]

            def fast_residuals(x):
                fits = m.calculate_rwa(x)
                cur_data, cur_err = data, err
                if shifted:
                    m.sess.update_data()
                    cur_data, cur_err = list(data), list(err)
                    for i in shifted:
                        cur_data[i] = np.asarray(m.sess.data_segments(i), float)
                        cur_err[i] = np.asarray(m.sess.data_segments_errors(i), float)
                out = []
                for i, (fi, nom) in enumerate(seg_norm):
                    n = float(x[fi]) if fi is not None else nom
                    f = np.asarray(fits[i], float)
                    e = cur_err[i] * n
                    with np.errstate(divide="ignore", invalid="ignore"):
                        r = np.where(e != 0, (f - cur_data[i] * n) / e, 0.0)
                    out.append(r)
                return np.concatenate(out) if out else np.zeros(0)

            # The fast path is checked against the engine's own residuals the
            # first time a Jacobian is requested (residual_jacobian costs ~14
            # forward evaluations single-threaded on a large model, so
            # it is not run at start-up just for the check).
            fast_state = {"ok": True, "checked": False}
            # residual_jacobian(skip=): held columns not finite-differenced (AZURE2 2026-10-30); older pyazr: all
            import inspect
            try:
                jac_skip = {"ok": "skip" in inspect.signature(m.residual_jacobian).parameters}
            except (TypeError, ValueError):
                jac_skip = {"ok": False}
            info["fast_residuals"] = "unchecked"
            conn.send(("ready", info))
            while True:
                cmd, payload = conn.recv()
                try:
                    if cmd == "quit":
                        conn.send(("ok", None))
                        break
                    elif cmd == "eval":
                        x = np.asarray(payload, float)
                        conn.send(("ok", float(np.sum(m.calculate_chi2_rwa(x)))))
                    elif cmd == "resid":
                        x, want_jac = payload[0], payload[1]
                        skip = payload[2] if len(payload) > 2 else None
                        x = np.asarray(x, float)
                        if want_jac:
                            if skip and jac_skip["ok"]:
                                try:
                                    r, J = m.residual_jacobian(x, skip=skip)
                                except TypeError:          # a .py that knows skip on an .so that does not
                                    jac_skip["ok"] = False
                                    r, J = m.residual_jacobian(x)
                            else:
                                r, J = m.residual_jacobian(x)
                            r = np.asarray(r, float)
                            if not fast_state["checked"]:
                                fast_state["checked"] = True
                                try:
                                    rf = fast_residuals(x)
                                    fast_state["ok"] = (rf.size == r.size
                                                        and np.allclose(rf, r, rtol=1e-9, atol=1e-12))
                                except Exception:
                                    fast_state["ok"] = False
                                if not fast_state["ok"]:
                                    print(f"rmfit worker: fast residuals disagree with the engine on {azr}; "
                                          f"falling back to residual_jacobian", flush=True)
                            conn.send(("ok", (r, np.asarray(J, float))))
                        elif fast_state["ok"]:
                            conn.send(("ok", (fast_residuals(x), None)))
                        else:
                            r = m.residual_jacobian(x)[0]
                            conn.send(("ok", (np.asarray(r, float), None)))
                    elif cmd == "fast_status":
                        conn.send(("ok", dict(fast_state)))
                    elif cmd == "transform":
                        conn.send(("ok", np.asarray(m.transform_rwa(np.asarray(payload, float)), float)))
                    elif cmd == "transform_all":
                        # 2026-10-24: physical values of EVERY parameter (fixed ones included, parameter order) for the
                        # free vector ``payload``, as pyazr's save_fit writes them.  A fixed parameter is fixed in rwa
                        # space; its observed value moves with the other amplitudes of its level.
                        x = np.asarray(payload, float)
                        full = np.asarray(m.sess.params_all_rwa(), float).copy()
                        free = [i for i, p in enumerate(m.parameters) if not p.fixed]
                        if x.size != len(free):
                            raise ValueError(f"transform_all: {x.size} values for {len(free)} free parameters")
                        full[free] = x
                        conn.send(("ok", np.asarray(m.transform_all_rwa(full, free_only=False), float)))
                    elif cmd == "cross":
                        x = np.asarray(payload, float)
                        conn.send(("ok", [np.asarray(c, float) for c in m.calculate_rwa(x)]))
                    elif cmd == "data":
                        conn.send(("ok", dict(
                            energies=[np.asarray(e, float) for e in m.energies],
                            excitation=[np.asarray(e, float) for e in m.excitation_energies],
                            cross=[np.asarray(c, float) for c in m.cross],
                            cross_err=[np.asarray(c, float) for c in m.cross_err],
                            # 2026-09-28: per-point data angles (plots split multi-angle segments by detector)
                            angles=[np.asarray(a, float) for a in getattr(m, "angles", [])])))
                    elif cmd == "write_output":
                        m.write_output_files(np.asarray(payload, float))
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
                        # everything a bake needs, so the file work happens in
                        # the parent and no second session ever shares this process
                        conn.send(("ok", dict(
                            parameters=list(m.parameters),
                            pairs=list(m.pairs),
                            names=[p.name for p in m.parameters],
                            fixed=[bool(f) for f in m.fixed_params],
                            all_rwa=np.asarray(m.sess.params_all_rwa(), float))))
                    else:
                        conn.send(("err", f"unknown command {cmd!r}"))
                except Exception:
                    conn.send(("err", traceback.format_exc()))
    except Exception:
        try:
            conn.send(("fatal", traceback.format_exc()))
        except Exception:
            pass
