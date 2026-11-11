# SRIM-2013 data files

The two files in this directory are data files of **SRIM-2013**, *The Stopping and Range of Ions in Matter*, by
J. F. Ziegler, J. P. Biersack and M. D. Ziegler (www.srim.org):

| File | SHA-256 | What |
|---|---|---|
| `SCOEF03.DAT` | `a328eaaa04e48efed1d15cfbcdf810849115ff352feaac1ee323c3bf4229ecea` | one row per element Z = 1–92: masses, densities, the electronic stopping fit coefficients for H, He and Li ions, the heavy-ion scaling polynomials, mean excitation energies (solid, gas), shell-correction polynomials, high-energy effective-charge parameters |
| `SNUC03.dat` | `4c7ff553b7be34a1b28a31a6ac89055b800392cbbf0401784da6afc5a693a030` | one row per pair Z1 ≤ Z2 (4278 rows): the four coefficients of SRIM's nuclear stopping fit |

They came to AZURE2 through **LunaSRIM 1.0.0** (J. Skowronski, https://github.com/skowrons94/LunaSRIM, the Python
re-implementation of SRIM-2013 from which `stopping/` is ported; its `lunasrim/data/`), and are copied unchanged:
the SHA-256 above are those of LunaSRIM's files and of the files here (CRLF line ends kept; `.gitattributes` here
stops git from converting them).

**Column layout.** `SCOEF03.DAT`: two quoted names (symbol, element name), then 93 numbers; their meaning is the
table at the top of LunaSRIM's `lunasrim/data_files.py`, and `stopping/include/azure2/stopping/Tables.h`
(`SrimElement`, with the column of each field) with `stopping/src/Tables.cpp` (`ParseSrimData`) is the reader here.
`SNUC03.dat`: `Z1 Z2 c1 c2 c3 c4`, blank-separated; the formula is in `stopping/include/azure2/stopping/Nuclear.h`.

**How the library finds them.** At build time `stopping/CMakeLists.txt` compiles both files into the library
(`src/TableText.cpp.in`), so nothing is read at run time. If the environment variable `AZURE2_SRIM_DATA_DIR` names a
directory, the library reads the two files from there instead (file names matched without regard to case); a caller
can also read any directory with `ReadSrimData`. If this directory is absent when AZURE2 is built, nothing is
compiled in and the variable is required: a user then fetches the two files (for instance from LunaSRIM's
`lunasrim/data/`, or from an installation of SRIM-2013, whose `Data/` folder holds them) and points the variable at
them.

**Terms.** These files are **not covered by AZURE2's licence**. They remain under the terms of their authors.
LunaSRIM's licence says the same (its MIT licence covers its code, not these files), and its README describes SRIM's
terms as free copying and distribution for personal, non-commercial use, commercial use requiring a written licence
from SRIM.org. They are copied here on the decision of the PI, Jakub Skowronski, of **2026-11-11**: "I would use the
SRIM stopping powers from my LunaSRIM repository. They are free to use as long as SRIM is cited." Whether AZURE2's
public repository ships them is the owner's decision; without them the library works with the download step above.

**Citation.** Every result that uses these stopping powers cites SRIM:

> J. F. Ziegler, M. D. Ziegler and J. P. Biersack, *SRIM – The stopping and range of ions in matter (2010)*,
> Nucl. Instrum. Methods Phys. Res. B **268** (2010) 1818–1823, doi:10.1016/j.nimb.2010.02.091.

(the reference LunaSRIM's `CITATION.cff` gives), and, for the implementation, LunaSRIM (its `CITATION.cff`).
