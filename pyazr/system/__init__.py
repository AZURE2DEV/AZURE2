"""The built system of a session: its particle pairs, its level scheme as a printable view, and the compound
nuclei a project's pairs form.

- :class:`Pair`, :class:`PairSet` (``pairs.py``): the engine's records of the particle pairs (``Session.pairs``).
- :class:`LevelScheme`, :class:`SchemeLevel`, :class:`SchemeChannel` (``scheme.py``): pairs, J-groups, levels and
  channels grouped as a physicist reads them (``Session.level_scheme``; read-only).
- :class:`Nucleus`, :func:`nucleus_of`, :func:`nucleus_of_pair`, :func:`merge_projects`, :class:`MergeReport`
  (``nuclei.py``): which compound nucleus a pair forms, and several single-nucleus projects merged into one file.
"""
from .pairs import Pair, PairSet
from .scheme import LevelScheme, SchemeLevel, SchemeChannel
from .nuclei import Nucleus, MergeReport, merge_projects, nucleus_of, nucleus_of_pair

__all__ = ["Pair", "PairSet", "LevelScheme", "SchemeLevel", "SchemeChannel",
           "Nucleus", "MergeReport", "merge_projects", "nucleus_of", "nucleus_of_pair"]
