# AZURE2

AZURE2 is a software package for performing multi-channel, multi-level
[R-matrix](https://en.wikipedia.org/wiki/R-matrix) analyses of low-energy
nuclear reaction and scattering data. It provides a Qt graphical setup utility,
a fast C++ calculation engine, parameter fitting via Minuit2, optional Bayesian
(MCMC) sampling, and an in-process Python interface (`pyazr`, built on pybind11)
for scripting and external samplers.

**Documentation:** <https://rdeboer1.github.io/AZURE2/> — user guide, the
physics and conventions, and the C++ API reference at
[`/api/`](https://rdeboer1.github.io/AZURE2/api/).

Upstream project: <https://azure.nd.edu/> · Source:
<https://github.com/rdeboer1/AZURE2>

---

## What's in this repository

| Path           | Description |
|----------------|-------------|
| `src/`         | Core C++ R-matrix engine (`CNuc`, `EData`, cross-section/χ² calculation, output). |
| `include/`     | Public headers for the core engine. |
| `gui/`         | Qt5 graphical setup utility (`AZURESetup`). |
| `api/`         | Python bindings: `AZUREAPI` C++ API + the pybind11 `_azure2` module. |
| `numcmc/`      | Affine-invariant ensemble MCMC sampler (`USE_MCMC`). |
| `pyazr/`       | Python package that drives AZURE2 in-process via the `_azure2` module. |
| `external/`    | Vendored third-party code, built in-tree: `minuit2/` (ROOT Minuit2), `coul/` (Coulomb wave functions), `erya/` (SRIM stopping powers + pugixml). |
| `tests/`       | Physics regression suite — reference evaluations and the runner. |
| `cmake/`       | Custom CMake find-modules (e.g. `FindQwt.cmake`) and toolchains. |
| `packaging/`   | Distribution: `docker/` images; the AppImage, Windows and macOS packaging live in the CI workflow. |
| `scripts/`     | Convenience build scripts (Linux, macOS, Windows, Docker). |
| `examples/`    | Example run scripts (GUI / MCMC via Docker). |
| `docs/`        | Documentation sources: Sphinx under `source/`, Doxygen via `Doxyfile`. |
| `.github/`     | Continuous integration. |
| `.claude/`     | [Claude Code](https://claude.com/claude-code) project skills (see below). |

> **Note:** Minuit2, the Coulomb library, the SRIM utilities and pugixml are
> vendored under `external/` and compiled as part of the build, and the MCMC
> sampler in `numcmc/` is built alongside them. None of these need to be
> downloaded or installed separately, unlike in older versions of AZURE2.

---

## Ready-built downloads

If you only want to *run* AZURE2, take one of these instead of building it.
Each is self-contained: Qt, Qwt, GSL, readline and the OpenMP runtime are
inside, so nothing else has to be installed.

| Platform | Download | How to run it |
|---|---|---|
| **Windows** (x86_64) | [`AZURE2-windows-x86_64.zip`](https://github.com/AZURE2DEV/AZURE2/releases/download/continuous-windows/AZURE2-windows-x86_64.zip) | Unzip anywhere, run `AZURE2.exe`. |
| **macOS** | [`AZURE2-macos.dmg`](https://github.com/AZURE2DEV/AZURE2/releases/download/continuous-windows/AZURE2-macos.dmg) | Open the disk image, drag `AZURE2.app` onto Applications. |
| **Linux** (x86_64) | build artifact `AZURE2-linux` on any [run page](https://github.com/AZURE2DEV/AZURE2/actions) | An AppImage: `chmod +x` it and run it. |

These are rolling builds of the latest `dev` commit, replaced on every push, so
they are for using the code rather than citing it — for a citable version use a
tagged release.

On macOS the first launch is refused with "the developer cannot be verified":
the build is signed ad-hoc rather than notarized with an Apple Developer ID.
Right-click (or control-click) `AZURE2.app`, choose **Open**, and confirm. Only
the first launch needs this. The macOS build is compiled for Intel and runs on
Apple Silicon through Rosetta.

---

## Dependencies

These are the only components you need to install yourself; everything else is
bundled in-tree.

**Build tools**
- A C++ compiler with OpenMP support (GCC or Clang)
- [CMake](https://cmake.org/) ≥ 3.16
- `pkg-config` — CMake locates GSL and Qwt through it, so configuring fails
  without it even when both libraries are installed

**Libraries**
- [GSL](https://www.gnu.org/software/gsl/) — GNU Scientific Library (math routines)
- [Qt5](https://www.qt.io/) — `Core`, `Widgets`, `Svg`, `Script` (for the GUI; `BUILD_GUI`)
- [Qwt](https://qwt.sourceforge.io/) (Qt5 build) — for the in-app plotting tab (`USE_QWT`)
- [Readline](https://tiswww.case.edu/php/chet/readline/rltop.html) — for CLI input (`USE_READLINE`)

**For the Python client (`pyazr`)**
- Python 3 with [NumPy](https://numpy.org/)
- [pybind11](https://pybind11.readthedocs.io/) ≥ 2.6 (to build the `_azure2` module)

### Installing dependencies

**Ubuntu / Debian**
```bash
sudo apt-get update
sudo apt-get install build-essential cmake pkg-config libgsl-dev libreadline-dev \
    qtscript5-dev libqwt-qt5-dev libqt5svg5-dev qtwebengine5-dev \
    python3 python3-numpy
```

**macOS (Homebrew)**
```bash
brew install cmake pkg-config gsl readline libomp qt@5 qwt-qt5
```

Three Homebrew-specific notes, none of which apply to a conda/miniforge
toolchain (which ships its own compiler, OpenMP and Qt in a single prefix):

- **`libomp` is required.** Apple's clang understands no `-fopenmp` and ships
  no OpenMP runtime, so without it CMake stops at `find_package(OpenMP)`.
- **`qt@5` and `qwt-qt5` are keg-only**, meaning Homebrew deliberately leaves
  them off the default search path. CMake asks `brew` for their locations
  automatically, so a plain `cmake -S . -B build` works; you only need
  `-DCMAKE_PREFIX_PATH="$(brew --prefix qt@5);$(brew --prefix qwt-qt5)"` if
  `brew` is not on your `PATH`.
- The plotting tab needs **`qwt-qt5`**, not `qwt` — plain `qwt` is the Qt6
  build and will not link against a Qt5 GUI.

Both Qt5 formulae are deprecated upstream and Homebrew will disable them on
2027-05-19, which is the real deadline for a Qt6 migration.

---

## Building from source

```bash
git clone https://github.com/rdeboer1/AZURE2.git
cd AZURE2
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The resulting executable is `build/src/AZURE2`. `CMAKE_BUILD_TYPE` defaults to
`Release`; set it to `Debug` explicitly if you want an unoptimized build with
symbols.

### Convenience scripts

Platform build scripts that auto-detect Qwt and set sensible options live in
`scripts/`:

```bash
source scripts/build_linux.sh     # Linux       -> build-linux/
source scripts/build_macos.sh     # macOS .app  -> build-macos/
source scripts/build_windows.sh   # Windows cross-compile (MinGW)
source scripts/build_docker.sh    # Docker image
```

### Build options

Pass options to CMake with `-D<OPTION>=ON|OFF` (or edit them interactively with
`ccmake ..`).

| Option                  | Default | Description |
|-------------------------|:-------:|-------------|
| `BUILD_GUI`             | ON      | Build and link the Qt graphical setup utility. |
| `USE_QWT`               | ON      | Include the built-in plotting tab (needs Qwt). |
| `USE_API`               | ON      | Build the `AZUREAPI` library and the pybind11 `_azure2` module (required for `pyazr`). |
| `USE_MCMC`              | ON      | Enable MCMC Bayesian sampling via the bundled `numcmc`. |
| `USE_ERYA`              | ON      | Enable the SRIM stopping-power utilities. |
| `USE_READLINE`          | ON      | Use Readline for console input. |
| `USE_STAT`              | ON      | Use `stat()` for directory checks (turn OFF on Windows). |
| `USE_NLOPT`             | OFF     | Use NLopt as an alternative minimizer (expects an `nlopt/` tree). |
| `BUILD_LIBRARY`         | OFF     | Build AZURE2 as a library. |
| `BUILD_MACOS_BUNDLE`    | OFF     | Produce a macOS `.app` bundle / DMG. |
| `CROSS_COMPILE_WINDOWS` | OFF     | Cross-compile a Windows binary with MinGW. |

If GSL or Qwt are installed in a non-standard location, point CMake at them with
`-DCMAKE_PREFIX_PATH=/path/to/prefix`.

---

## Running AZURE2

AZURE2 reads a `.azr` configuration file describing the compound nucleus,
channels, levels, and data segments.

**Graphical mode** (default — opens the setup utility):
```bash
./build/src/AZURE2 path/to/config.azr
```

**Console mode** (no GUI):
```bash
./build/src/AZURE2 --no-gui path/to/config.azr
```

Useful flags (`AZURE2 --help` for the full list):

| Flag                  | Effect |
|-----------------------|--------|
| `--no-gui`            | Run without the graphical setup utility. |
| `--use-brune`         | Use the alternative level matrix of C. R. Brune. |
| `--gsl-coul`          | Use GSL Coulomb functions (faster, less accurate). |
| `--use-rmc`           | Reich–Moore approximation for capture (neutron capture). |
| `--ignore-externals`  | Ignore external resonant capture when the internal width is zero. |
| `--no-transform`      | Skip the initial parameter transformations. |

---

## Python interface (`pyazr`)

`pyazr` embeds the AZURE2 engine in-process: the R-matrix code is compiled into
a pybind11 extension module (`_azure2`), and an `azure2()` object is a real
`AZUREAPI` session living in the interpreter. No subprocesses, no sockets, no
port bookkeeping — which also makes it easy to plug AZURE2 into external
samplers such as [`emcee`](https://emcee.readthedocs.io/)
or [`brick`](https://github.com/odell/brick).

```python
from pyazr import azure2

# Build the in-process engine bound to a configuration file.
azr = azure2("config.azr")

# Experimental data, grouped by segment.
energies = azr.energies
cross    = azr.cross

# Evaluate the model for a set of (free) parameters.
fit = azr.calculate(azr.params)

# Inspect the fit parameters with full physics metadata.
print(azr.parameters.table())          # readable overview
for w in azr.parameters.widths.free:   # free reduced-width amplitudes
    print(w.name, "J^pi", w.jpi, "L", w.L, "S", w.S)

azr.close()
```

The `parameters` view returns a `ParameterSet` of `Parameter` objects, each
carrying what the parameter *is* — for R-matrix parameters the level it belongs
to (`J`, `parity`, `level_energy`) and for widths the channel (`L`, `S`, `pair`,
`radiation_type`) — plus filtered views (`.free`, `.energies`, `.widths`,
`.norms`, `.shifts`) and lookups (`.by_level(...)`, `.by_name(...)`).

Build the `_azure2` module with CMake (`USE_API=ON`, the default) — it lands in
`pyazr/`, so `import pyazr` picks it up from the repository root; a `pip install`
ships it as package data.

### Fetching experimental data (`pyazr.nds`)

`pyazr.nds` is the package's client for the IAEA data services: **EXFOR** for
measured cross sections, differential cross sections, analyzing powers and
yields, and **LiveChart/ENSDF** for evaluated level schemes and gamma
transitions. It converts what it fetches into AZURE2's own conventions — lab
energies, barns or barns/sr — so a published dataset becomes a data segment
without hand editing:

```python
from pyazr import nds, AzrModel

hits = nds.search_exfor(target="C-13", reaction="p,g", quantity="SIG")
data = nds.fetch_exfor("O2599004")                # S-factor, reported in b·keV
kw   = data.to_azr("run/data", entrance=1, exit=2,   # -> barns, lab energies
                   observable="total-capture")
AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")

nds.fetch_levels("14n")                           # ENSDF level scheme
nds.resolve_doi(nds.reference("O2599004"))        # the paper behind the data
```

Units and frames are handled from the column headers (`EN-CM` → lab, `B*KEV`
S-factors → barns, `NB/SR` → b/sr, ratio-to-Rutherford → b/sr). It needs only
the standard library beyond NumPy, and works from the shell too:

```bash
python -m pyazr.nds search --target C-13 --reaction p,g --quantity SIG
python -m pyazr.nds download O2599004 -o data/skowronski.dat
```

The Qt setup utility has its own EXFOR dialog backed by `gui/src/ExforData.cpp`;
the two are independent implementations of the same Web-API, so a parsing rule
learned by either belongs in both. See `pyazr/examples/exfor_fetch.py`.

---

## Containers

**Docker** — a self-contained image (Ubuntu 22.04 + ROOT + Python tooling) is
defined in `packaging/docker/Dockerfile.azure2`:
```bash
source scripts/build_docker.sh
source examples/run_gui.sh      # run the GUI from the container
```

**Singularity / Apptainer** — for HPC environments, convert the Docker image:
```bash
sudo apptainer build AZURE2.sif docker-daemon://azure2:latest
```
Copy the resulting `AZURE2.sif` to your HPC resource.

---

## Downloads

Every push builds AZURE2 on Linux, macOS and Windows and attaches the results
to the run under **Actions → build → Artifacts**. Tagging a commit `v*`
collects the same artifacts into a draft GitHub release.

| Platform | Artifact | Notes |
|----------|----------|-------|
| Linux    | `AZURE2-x86_64.AppImage` | Self-contained; `chmod +x` and run. No installation, no root. |
| Windows  | `AZURE2.exe` + Qt runtime | Built with MSYS2/mingw-w64. |
| macOS    | `AZURE2` (x86_64) | Not yet a signed `.app` bundle — see `dist/README-packaging.txt` in the artifact. |

---

## Tests

`tests/` holds reference evaluations. Each is run through AZURE2 and its
χ² compared against a recorded result, which catches a change that still
compiles and still runs but moves the physics:

```bash
./tests/run_tests.sh                  # uses build/src/AZURE2
./tests/run_tests.sh path/to/AZURE2   # or an explicit binary
TOL=1e-4 ./tests/run_tests.sh         # tighter tolerance
```

The suite runs on all three platforms in CI. AZURE2 currently agrees to within
about 1×10⁻⁵ across GCC, Clang and MinGW, so the default relative tolerance is
1×10⁻³ — loose enough for that spread, far tighter than any real regression.

**Adding an evaluation.** Create `tests/<name>/` containing `<name>.azr`, its
`data/`, and `expected/chiSquared.out` taken from a run you trust. Cases are
discovered automatically; nothing in the runner needs editing. The comparison
covers the total χ², each segment's χ², and each segment's point count (that
one exactly — a change there means data was dropped or misread).

---

## Documentation

The site at <https://rdeboer1.github.io/AZURE2/> is published by CI from the
default branch (`qt5`); pushes to other branches build the documentation to
catch breakage but do not publish. Sources live in `docs/`; the generated output is not tracked.
To build it locally:

```bash
pip install -r docs/requirements.txt
make -C docs html     # user guide  -> docs/_build/html/index.html
make -C docs api      # C++ internals via Doxygen -> docs/api/html/index.html
```

---

## Claude Code skill

`.claude/skills/azure2-eval/SKILL.md` is a project-scoped
[Claude Code](https://claude.com/claude-code) skill describing how to drive
AZURE2 — the CLI menu modes, the `pyazr` API, adding and removing levels,
decomposing cross sections, and the conventions that are easy to get wrong
(lab vs. centre-of-mass frames, segment indexing, the external-capture integral
caches). It is picked up automatically when Claude Code runs in this
repository; no setup is needed. Editing it is the way to teach Claude something
new about the project.

---

## License

AZURE2 is distributed under the terms of the GNU General Public License v3.
See the upstream project at <https://azure.nd.edu/> for details and citation
information.
