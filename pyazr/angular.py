"""Angular distributions at energies of your choosing.

AZURE2 computes Legendre coefficients only for segments declared as angular
distributions, and segments live in the ``.azr`` file rather than being
something a running instance can be asked for. So obtaining a distribution at
an arbitrary energy means writing a file that requests it and evaluating that.

:func:`angular_distribution` does exactly that and cleans up after itself. For
coefficients at the grids a model already declares, use
``azure2.calculate_angular_dists`` on a live instance instead -- no temporary
file, no second process.
"""

import os
import tempfile

import numpy as np

from .azrfile import AzrModel
from .azure2 import azure2


def angular_distribution(azr_file, energies, entrance=1, exit=1, order=4,
                         params=None, use_rwa=True, lab=True, **azure2_kwargs):
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
        Parameter vector to evaluate at. Defaults to the model's own values.
        Interpreted as reduced-width amplitudes unless ``use_rwa=False``.
    use_rwa : bool
        Whether ``params`` is in reduced-width-amplitude space.
    lab : bool
        Whether ``energies`` are laboratory (default) or centre-of-mass.
    **azure2_kwargs
        Passed through to :class:`~pyazr.azure2.azure2` (``binary``, ``cwd``…).

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

    model = AzrModel.from_file(azr_file)

    # Convert to the lab frame the file expects. The factor is the same one the
    # extrapolation grids use: E_lab = E_cm (m_beam + m_target) / m_target.
    if lab:
        e_lab = energies
    else:
        pair = model.pairs()[entrance] if hasattr(model, "pairs") else None
        if pair is None:
            raise ValueError(
                "centre-of-mass input needs the entrance pair's masses, which "
                "this model does not expose; pass lab energies instead")
        e_lab = energies * (pair["m1"] + pair["m2"]) / pair["m2"]

    # One single-point extrapolation per energy, so segment i is energy i and
    # no interpolation onto a grid is involved.
    model.clear_extrapolations()
    for e in e_lab:
        model.add_extrapolation(entrance=entrance, exit=exit,
                                e_min=float(e), e_max=float(e), e_step=1.0,
                                observable="angular-distribution", order=order)

    tmpdir = tempfile.mkdtemp(prefix="pyazr-angdist-")
    tmp_azr = os.path.join(tmpdir, "_angular_distribution.azr")
    try:
        model.write(tmp_azr)

        # The .azr keeps its paths relative to itself, so the instance has to
        # run where the original file lives.
        azure2_kwargs.setdefault("cwd", os.path.dirname(os.path.abspath(azr_file)) or ".")

        with azure2(tmp_azr, **azure2_kwargs) as m:
            m.extrap_mode()

            # The grid energies do not depend on the parameters, and
            # calculate_energies speaks the physical convention -- handing it a
            # reduced-width vector is what the two conventions being separate
            # means, and AZURE2 does not survive it.
            calc_e = m.calculate_energies(m.params)

            if use_rwa:
                x = m.params_rwa if params is None else params
                dists = m.calculate_angular_dists_rwa(x)
            else:
                x = m.params if params is None else params
                dists = m.calculate_angular_dists(x)

        out_e = np.full(energies.size, np.nan)
        out_c = np.full((energies.size, order + 1), np.nan)
        for i in range(min(len(dists), energies.size)):
            rows = [r for r in dists[i] if r.size]
            if not rows:
                continue
            coeffs = rows[0]
            out_c[i, :min(coeffs.size, order + 1)] = coeffs[:order + 1]
            if i < len(calc_e) and len(calc_e[i]):
                out_e[i] = calc_e[i][0]
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
