"""Spin observables (isDiff 9 segments): token grammar and literature aliases.

A spin-observable segment names its observable with one token (the engine's grammar,
``include/SpinObservable.h``)::

    token := side [ '>' side ]      side := term { '*' term }      term := who ':' op
    who   := b | t | e1 | e2        beam, target, exit particle 1 (light) or 2 (heavy)
    op    := x | y | z              vector, S_i / j (Pauli matrices for spin 1/2)
           | xx yy zz xy xz yz      Cartesian tensor (spin 1 only, Ohlsen)
           | Tkq | iTkq             spherical tensor (Madison), e.g. T20, iT11

Examples: ``b:y`` (A_y), ``t:y`` (A_0y), ``b:y*t:y`` (C_yy), ``b:yy`` (spin-1 A_yy), ``b:T20``,
``e1:y`` (outgoing polarization), ``b:y>e1:y`` (transfer K_y^y').

Literature names are ambiguous between spin-1/2 spin correlations and spin-1 tensor analyzing
powers (A_yy is either), so files carry only the explicit grammar; :func:`resolve` maps a
literature name to a token for given beam and target spins.
"""
import re

_TERM = re.compile(r"^(b|t|e1|e2):([xyz]|[xyz][xyz]|i?T\d\d)$")


def check(token):
    """Raise ValueError unless ``token`` is a valid spin-observable token (grammar only; the
    engine also checks it against the particles' spins)."""
    if not token or any(c.isspace() for c in token):
        raise ValueError(f"spin observable {token!r}: one token without blanks expected")
    sides = token.split(">")
    if len(sides) > 2:
        raise ValueError(f"spin observable {token!r}: more than one '>'")
    seen, n_i = set(), 0
    for n, side in enumerate(sides):
        if not side:
            raise ValueError(f"spin observable {token!r}: empty side")
        whos = []
        for term in side.split("*"):
            m = _TERM.match(term)
            if not m:
                raise ValueError(f"spin observable {token!r}: bad term {term!r}")
            who, op = m.group(1), m.group(2)
            if who in seen:
                raise ValueError(f"spin observable {token!r}: {who!r} appears twice")
            seen.add(who)
            whos.append(who)
            if op.startswith("i"):
                n_i += 1
            if "T" in op:
                k, q = int(op[-2]), int(op[-1])
                if k < 1 or q > k:
                    raise ValueError(f"spin observable {token!r}: {op} needs 1 <= k, 0 <= q <= k")
        entrance = [w in ("b", "t") for w in whos]
        if len(sides) == 2 and n == 0 and not all(entrance):
            raise ValueError(f"spin observable {token!r}: exit particle before '>'")
        if len(sides) == 2 and n == 1 and any(entrance):
            raise ValueError(f"spin observable {token!r}: entrance particle after '>'")
        if len(sides) == 1 and any(entrance) and not all(entrance):
            raise ValueError(f"spin observable {token!r}: write <entrance>'>'<exit>")
    if n_i > 1:
        raise ValueError(f"spin observable {token!r}: at most one iT_kq term")
    return token


def resolve(name, beam_spin=None, target_spin=None):
    """Token for a spin observable given as a token or a literature name.

    Literature names: ``Ay``/``Ay0`` (beam A_y), ``A0y`` (target), ``Py`` (outgoing polarization of
    the light exit particle), ``Kyy`` (transfer beam -> light exit particle), ``iT11``, ``T20``,
    ``T21``, ``T22`` (beam), and ``Axx``, ``Ayy``, ``Azz``, ``Axz``, ``Cxx`` ... : for a spin-1 beam
    ``Aij`` is the tensor analyzing power, for spin-1/2 on spin-1/2 ``Aij``/``Cij`` is the spin
    correlation (beam i, target j); any other case raises, rather than guess.
    """
    if ":" in name:
        return check(name)
    simple = {"Ay": "b:y", "Ay0": "b:y", "A0y": "t:y", "Py": "e1:y", "Kyy": "b:y>e1:y",
              "iT11": "b:iT11", "T20": "b:T20", "T21": "b:T21", "T22": "b:T22"}
    if name in simple:
        return simple[name]
    m = re.match(r"^([AC])([xyz])([xyz])$", name)
    if m:
        kind, i, j = m.groups()
        if kind == "A" and beam_spin is not None and abs(beam_spin - 1.0) < 1e-9:
            return f"b:{i}{j}"
        if (beam_spin is not None and target_spin is not None and abs(beam_spin - 0.5) < 1e-9
                and abs(target_spin - 0.5) < 1e-9):
            return f"b:{i}*t:{j}"
        if kind == "C" and beam_spin is None:
            return f"b:{i}*t:{j}"
        raise ValueError(f"{name!r} is ambiguous for beam spin {beam_spin}, target spin {target_spin}: "
                         "write the token (e.g. 'b:yy' or 'b:y*t:y')")
    raise ValueError(f"unknown spin observable {name!r}")
