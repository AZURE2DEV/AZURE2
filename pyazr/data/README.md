# Atomic mass table (AME2020)

`mass_1.mas20` is the atomic mass table of the 2020 Atomic Mass Evaluation,
unchanged from the Atomic Mass Data Center's distribution:

- W.J. Huang, M. Wang, F.G. Kondev, G. Audi and S. Naimi, "The AME 2020 atomic
  mass evaluation (I)", Chinese Physics C 45, 030002 (2021).
- M. Wang, W.J. Huang, F.G. Kondev, G. Audi and S. Naimi, "The AME 2020 atomic
  mass evaluation (II)", Chinese Physics C 45, 030003 (2021).
- Source: https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
  (MD5 6d28b75833cf53c7cc230223f63da6f6, checked 2026-10-29).

It is the one mass table of AZURE2: `pyazr.masses` reads it, and the GUI
embeds the same file as a Qt resource (`gui/src/AZURESetup.qrc`), so both fill
masses and separation energies from identical numbers.  The masses are
**atomic** masses (electrons included), the convention AZURE2 projects use
throughout; a separation energy is then `S = m1 + m2 - M` of the compound
nucleus, with the electrons cancelling.

To move to a newer evaluation, replace the file (keeping the AME fixed-width
format), update the version string in `pyazr/masses.py` (`AME_VERSION`) and in
`gui/include/MassTable.h` (`MassTable::version()`), and rerun
`tests/pyazr/masses_test.py` and the GUI `mass_table` test, whose reference
values are AME2020's.
