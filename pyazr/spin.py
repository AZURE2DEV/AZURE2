"""Literature names of spin observables: the operator a format-3 ``spin`` record takes for a name of the papers.

A spin observable is named in a project by one operator in the engine's grammar (``spin_operator`` of
``docs/source/reference/format3.rst``): ``b:y`` (A_y), ``t:y`` (the target analysing power), ``b:y*t:y`` (C_yy),
``b:yy`` and ``b:T20`` (tensor analysing powers of a spin-1 beam), ``e1:y`` (the outgoing polarization),
``b:y>e1:y`` (a polarization transfer).  The papers' names are shorter and one of them is ambiguous: A_yy is the
tensor analysing power of a spin-1 beam and the spin correlation of two spin-1/2 particles.  :func:`spin_operator`
gives the operator of a name for the spins of the beam and the target, and refuses where the spins do not decide.

The grammar itself is the engine's (it checks an operator against the four particles when a project is read):
this module holds the table of names only.
"""
import re

from .errors import SpinOperatorError

#: The names that mean one operator whatever the spins: name -> operator.
LITERATURE_NAMES = {"Ay": "b:y", "Ay0": "b:y", "A0y": "t:y", "Py": "e1:y", "Kyy": "b:y>e1:y",
                    "iT11": "b:iT11", "T20": "b:T20", "T21": "b:T21", "T22": "b:T22"}

_PAIR = re.compile(r"^([AC])([xyz])([xyz])$")
_FIX = "write the operator itself (its grammar: spin_operator in docs/source/reference/format3.rst)"


def _is(spin, value):
    return spin is not None and abs(float(spin) - value) < 1e-9


def spin_operator(name, beam_spin=None, target_spin=None):
    """The operator of a spin observable given by its name in the literature (or by the operator itself).

    ``name``: ``Ay`` or ``Ay0`` (the beam's A_y), ``A0y`` (the target's), ``Py`` (the outgoing polarization of the
    light exit particle), ``Kyy`` (the transfer from the beam to the light exit particle), ``iT11``, ``T20``,
    ``T21``, ``T22`` (the beam's spherical tensors), and ``Aij`` or ``Cij`` with i, j of x, y, z.  For a spin-1 beam
    ``Aij`` is the tensor analysing power ``b:ij``; for spin 1/2 on spin 1/2 ``Aij`` and ``Cij`` are the spin
    correlation ``b:i*t:j``; ``Cij`` without spins is the correlation too.  A name with a colon is taken as an
    operator and returned unchanged (the engine checks it when the project is read).

    ``beam_spin``, ``target_spin``: the spins of the two particles of the entrance pair (0.5, 1, ...), needed for
    ``Aij`` only.

    Raises :class:`pyazr.SpinOperatorError` (a ``ValueError``) for a name it does not know and for ``Aij`` where the
    spins do not decide: it does not guess.
    """
    name = str(name).strip()
    if ":" in name:
        return name
    if name in LITERATURE_NAMES:
        return LITERATURE_NAMES[name]
    m = _PAIR.match(name)
    if m:
        kind, i, j = m.groups()
        if kind == "A" and _is(beam_spin, 1.0):
            return f"b:{i}{j}"
        if _is(beam_spin, 0.5) and _is(target_spin, 0.5):
            return f"b:{i}*t:{j}"
        if kind == "C" and beam_spin is None and target_spin is None:
            return f"b:{i}*t:{j}"
        what = f"spin observable {name!r} with beam spin {beam_spin} and target spin {target_spin}"
        allowed = (f"b:{i}{j} (the tensor analysing power of a spin-1 beam) or b:{i}*t:{j} (the spin correlation of "
                   "spin 1/2 on spin 1/2)")
        raise SpinOperatorError(f"{what}: the spins do not say which of {allowed} is meant; {_FIX}",
                                what=what, allowed=allowed, fix=_FIX)
    what = f"spin observable {name!r}"
    allowed = ", ".join(sorted(LITERATURE_NAMES)) + ", Aij and Cij with i, j of x, y, z"
    raise SpinOperatorError(f"{what}: not a name of the table ({allowed}); {_FIX}", what=what, allowed=allowed, fix=_FIX)
