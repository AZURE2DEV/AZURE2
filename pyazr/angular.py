"""Angular distributions at energies of your choosing.

AZURE2 computes Legendre coefficients only for segments declared as angular
distributions, and segments live in the ``.azr`` file rather than being
something a running instance can be asked for. So obtaining a distribution at
an arbitrary energy means writing a file that requests it and evaluating that.

:func:`angular_distribution` does exactly that and cleans up after itself. For
coefficients at the grids a model already declares, use
``Session.evaluate(x).data[i].extra["legendre"]`` on a live session instead -- no
temporary file, no second session.
"""

import os
import tempfile

import numpy as np

from .session import Session
from .project import Project


def angular_distribution(azr_file, energies, entrance=1, exit=1, order=4,
                         params=None, lab=True, **azure2_kwargs):
    r"""Legendre coefficients of the angular distribution at given energies.

    Parameters
    ----------
    azr_file : str
        The model to evaluate. It is read, never modified.
    energies : float or sequence of float
        Energies at which the distribution is wanted. **Laboratory** frame by
        default, matching what the ``.azr`` file holds; pass ``lab=False`` to
        give centre-of-mass energies instead.
    entrance, exit : int
        Particle-pair keys of the reaction.
    order : int
        Highest Legendre order to compute.
    params : array_like, optional
        The free vector (the fit vector, ``Session.params_rwa`` of the model's
        session) to evaluate at. Defaults to the model's own values.
    lab : bool
        Whether ``energies`` are laboratory (default) or centre-of-mass.
    **azure2_kwargs
        Passed through to :class:`~pyazr.Session` (``binary``, ``cwd``…).

    Returns
    -------
    energies_cm : ndarray
        Centre-of-mass energies actually calculated, in the requested order.
        AZURE2 returns centre-of-mass regardless of what went in.
    coefficients : ndarray, shape (n_energies, order + 1)
        The :math:`a_k` of :math:`W(\theta) = \sum_k a_k P_k(\cos\theta)`.
        Rows for energies AZURE2 could not evaluate are filled with NaN rather
        than silently dropped, so the rows always line up with the input.

    Notes
    -----
    Every call writes a temporary ``.azr`` and starts an AZURE2 process. For
    many energies, pass them all in one call rather than looping.
    """
    energies = np.atleast_1d(np.asarray(energies, dtype=float))
    if energies.size == 0:
        return np.empty(0), np.empty((0, order + 1))

    frame = "lab" if lab else "cm"
    tmpdir = tempfile.mkdtemp(prefix="pyazr-angdist-")
    tmp_azr = os.path.join(tmpdir, "_angular_distribution.azr")
    # The project with its extrapolations replaced by one single-point grid per energy, so segment i is energy i
    # and no interpolation onto a grid is involved: the Legendre coefficients a_k / a_0 (expansion a_over_a0), the
    # energy in the frame asked for.
    model = Project.open(azr_file)
    for name, r in model.records().items():
        if r.record_kind == "extrapolation":
            model.remove(name)
    reaction = model.reaction(entrance, exit)
    for i, e in enumerate(energies):
        model.add_extrapolation(f"angular_{i + 1}", reaction=reaction, quantity="legendre", unit="1",
                                expansion="a_over_a0", orders=f"0 {int(order)}",
                                energies=f"list {float(e)!r} MeV {frame}")
    write = model.save

    try:
        write(tmp_azr)

        # The .azr keeps its paths relative to itself, so the instance has to
        # run where the original file lives.
        azure2_kwargs.setdefault("cwd", os.path.dirname(os.path.abspath(azr_file)) or ".")

        with Session(tmp_azr, **azure2_kwargs) as m:
            m.extrap_mode()
            grids = m.evaluate(m.params_rwa if params is None else params).data

        out_e = np.full(energies.size, np.nan)
        out_c = np.full((energies.size, order + 1), np.nan)
        for i in range(min(len(grids), energies.size)):
            rows = [r for r in grids[i].extra.get("legendre", []) if r.size]
            if not rows:
                continue
            coeffs = rows[0]
            out_c[i, :min(coeffs.size, order + 1)] = coeffs[:order + 1]
            if len(grids[i].energies):
                out_e[i] = grids[i].energies[0]
        return out_e, out_c
    finally:
        # The temporary model and whatever AZURE2 wrote beside it.
        for name in os.listdir(tmpdir):
            try:
                os.remove(os.path.join(tmpdir, name))
            except OSError:
                pass
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass
