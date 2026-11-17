"""The engine's own ensemble sampler on a session (docs/upgrade/pyazr3.md section 12):
:meth:`SessionSampler.sample_with_engine` runs what ``AZURE2 mcmc`` runs -- the affine-invariant ensemble sampler of
the engine, with the command's options as keywords -- on the session's model, in the process, through the engine
function the command calls.  The chain is the command's, number for number.  Mixed into the session class
(``class Session(..., SessionSampler)``).

Example (T3 of ``tests/upgrade``)::

    s = pyazr.Session("T3.azr")
    r = s.sample_with_engine(walkers=20, steps=500, sampler_threads=2, overwrite_samples=True)
    r.chain.shape, r.log_prob.shape            # (500, 20, 8), (500, 20)
    r.chain[200:].reshape(-1, len(r.names))    # the draws after a burn-in of 200 steps

For a sampler of your own (emcee, a NUTS) the rule of the same posterior is ``s.objective().log_prob``.
"""
from dataclasses import dataclass, field
from typing import List

import numpy as np

from . import engine_text
from .errors import EngineError
from .locking import engine_call
from .objective import require_data_mode


@dataclass
class EngineSampleResult:
    """What :meth:`SessionSampler.sample_with_engine` returns.

    ``chain`` has the shape (steps, walkers, parameters): every walker's position after every step, the last axis in
    the order of ``names`` (the session's free parameters), in the units of ``basis`` (``"rwa"``: the fit vector,
    reduced width amplitudes; ``"physical"``: the project's physical parameters, widths in eV).  ``log_prob``,
    ``log_likelihood`` and ``log_prior`` have the shape (steps, walkers); ``log_prob`` is their sum, -1/2 of the
    fit's objective plus the normalisation of the Gaussian priors.  The arrays are ``samples.mcmc`` read back by the
    engine (written with 17 digits, so they are the doubles the sampler held): every step of the file, those of an
    earlier run this one continued included (``steps`` counts them).  ``status`` is ``"completed"`` or ``"stopped"``
    (a stop request ended the run; the files are complete to the last step).  ``samples_file`` and ``walkers_file``
    are the chain and the walkers' last positions (from which a later run continues), in ``output_dir``, where the
    sampler also writes the best walker's output files and ``param.mcmc`` every 10 steps.  ``log`` holds the
    sampler's lines (its table of the free parameters and priors, the acceptance fraction, its warnings).
    """
    chain: np.ndarray
    log_prob: np.ndarray
    log_likelihood: np.ndarray
    log_prior: np.ndarray
    names: List[str]
    basis: str
    walkers: int
    steps: int
    status: str
    samples_file: str
    walkers_file: str
    output_dir: str
    log: List[str] = field(repr=False, default_factory=list)


class SessionSampler:
    """Mixed into :class:`pyazr.Session`: uses its ``_engine``, ``mode`` and ``parameters``."""

    @engine_call
    def sample_with_engine(self, walkers, steps, *, sampler_threads=1, basis="rwa", parameter_spread=5.0,
                           energy_spread=1.0, overwrite_samples=False, x0=None, progress=None, stop=None):
        """The program's ensemble sampler (``AZURE2 mcmc``) on this session, through the engine's own sampler.

        The keywords are the command's options, with its ranges (``docs/source/reference/command_line.rst``); a value
        outside its range is refused with the command line's words (``ValueError``).

        Parameters
        ----------
        walkers : ``--walkers``: the size of the ensemble, 2 to 1000; more than the free parameters.
        steps : ``--steps``: steps per walker, 100 to 10^7.  The chain has ``walkers * steps`` rows.
        sampler_threads : ``--sampler-threads``: the walkers of a half-ensemble are evaluated on this many threads
            (>= 1).  The chain does not depend on it.
        basis : ``--basis``: ``"rwa"`` samples the fit vector (reduced width amplitudes), ``"physical"`` the
            project's physical parameters.
        parameter_spread : ``--parameter-spread``: the start ensemble's spread, percent of each parameter, above 0
            and up to 50 (level energies and the parameters with a prior of their own have their own rule).
        energy_spread : ``--energy-spread``: the start ensemble's spread of a level energy, keV, above 0 and up to 1.
        overwrite_samples : ``--overwrite-samples``: start a new chain.  ``False`` (the command's default): an
            existing ``samples.mcmc`` in the session's output directory is continued from ``walkers.mcmc``, and a
            file that already holds ``steps`` steps is returned as it is.
        x0 : the start, a fit vector (default ``parameters.x0``): what ``--parameter-file`` is to the command.  In
            the physical basis the engine transforms it.
        progress : called after every step with a dict: ``step`` and ``steps`` (of this run), ``log_prob``,
            ``log_likelihood``, ``log_prior`` (the step's best walker; what the program's progress line shows).  It
            must not call the session.  An exception in it ends the run and is raised afterwards.
        stop : a :data:`pyazr.session_fit.StopRequest` or a ``threading.Event``, looked at after every step; once set,
            the run ends after the step, with status ``"stopped"``.

        The command has no seed, burn-in or thinning: the sampler's seed is fixed in the engine (the same call gives
        the same chain; a continued run draws a new stretch), and every step is kept.  The Gaussian priors are the
        engine's (normalizations, energy shifts, ``<parameterSettings>``, a prior covariance) and the rows of the
        project's ``<mcmc>`` section.

        Returns an :class:`EngineSampleResult`.  Writes what the command's sampler writes into the session's output
        directory: ``samples.mcmc``, ``walkers.mcmc`` and, every 10 steps, the best walker's output files and
        ``param.mcmc``.  The session's parameter values are not changed.  A run the sampler refuses (no walker of the
        start ensemble with a finite probability, a row of ``<mcmc>`` that names no parameter, no free parameter)
        raises :class:`pyazr.EngineError` with the sampler's words.
        """
        require_data_mode(self, "sample_with_engine()")
        names = list(self.parameters.free.names)
        start = [] if x0 is None else [float(v) for v in np.asarray(x0, float).ravel()]
        if x0 is not None and len(start) != len(names):
            raise ValueError(f"x0 has {len(start)} values; the session has {len(names)} free parameters.")
        if stop is not None and not callable(getattr(stop, "is_set", None)):
            raise TypeError("stop must be a pyazr.session_fit.StopRequest or a threading.Event (has is_set()).")
        if basis not in ("rwa", "physical"):
            raise ValueError(f"basis = {basis!r}; allowed: 'rwa' (the fit vector) or 'physical' (--basis).")
        for name, value in (("walkers", walkers), ("steps", steps), ("sampler_threads", sampler_threads)):
            if isinstance(value, bool) or int(value) != value:
                raise ValueError(f"{name} = {value!r}; allowed: a whole number (--{name.replace('_', '-')}).")
        with engine_text.held() as lines:
            out = self._engine.sample_with_engine(
                walkers=int(walkers), steps=int(steps), sampler_threads=int(sampler_threads), basis=basis,
                parameter_spread=float(parameter_spread), energy_spread=float(energy_spread),
                overwrite_samples=bool(overwrite_samples), x0=start, progress=progress, stop=stop)
        if out["outcome"] != "sampled":
            said = [line for line in lines if line.lstrip().startswith(("ERROR", "Error", "MCMC"))] or lines[-3:]
            raise EngineError("sample_with_engine: the engine's sampler did not run: " + " ".join(said)
                              + "  The sampler's other lines are in the log 'pyazr.engine'.")
        stopped = bool(out["stop_requested"]) and int(out["steps"]) < int(steps)
        return EngineSampleResult(
            chain=out["chain"], log_prob=out["log_prob"], log_likelihood=out["log_likelihood"],
            log_prior=out["log_prior"], names=names, basis=basis, walkers=int(out["walkers"]),
            steps=int(out["steps"]), status="stopped" if stopped else "completed",
            samples_file=out["samples_file"], walkers_file=out["walkers_file"], output_dir=out["output_dir"],
            log=list(lines))
