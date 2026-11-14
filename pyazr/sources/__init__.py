"""The libraries pyazr reads nuclear data from (no engine needed; each module needs numpy only):

- :mod:`pyazr.sources.nds` -- IAEA EXFOR (``search_exfor``, ``fetch_exfor``, ``ExforData.to_azr``, a local copy of the
  EXFOR master files as fallback), LiveChart/ENSDF levels, gammas and ground states (``fetch_levels``,
  ``fetch_gammas``, ``fetch_ground_state``), references and DOIs; ``python -m pyazr.sources.nds`` on the command
  line.
- :mod:`pyazr.sources.masses` -- the AME2020 mass table (``pyazr/data/mass_1.mas20``): atomic and nuclear masses,
  separation energies.
- :mod:`pyazr.sources.exfor_coulomb`, :mod:`pyazr.sources.exfor_sfactor` -- the rules by which an EXFOR ratio to
  Rutherford or Mott and an EXFOR S factor become format-3 data sets (``ExforData.to_dataset``).

The functions of ``nds`` are also the top-level names ``pyazr.fetch_exfor``, ``pyazr.fetch_levels``, ...
"""
from . import exfor_coulomb, exfor_sfactor, masses, nds

__all__ = ["nds", "masses", "exfor_coulomb", "exfor_sfactor"]
