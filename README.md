# AZURE2

AZURE2 is a software package for performing multi-channel, multi-level
[R-matrix](https://en.wikipedia.org/wiki/R-matrix) analyses of low-energy
nuclear reaction and scattering data. It provides a Qt graphical setup utility,
a fast C++ calculation engine, parameter fitting via Minuit2, optional Bayesian
(MCMC) sampling, and a socket API with a Python client (`pyazr`) for scripting
and external samplers.

Upstream project: <https://azure.nd.edu/> · Source:
<https://github.com/rdeboer1/AZURE2>

---

## What's in this repository

| Path           | Description |
|----------------|-------------|
| `src/`         | Core C++ R-matrix engine (`CNuc`, `EData`, cross-section/χ² calculation, output). |
| `include/`     | Public headers for the core engine. |
| `gui/`         | Qt5 graphical setup utility (`AZURESetup`). |
| `api/`         | Socket API server (`--use-api`) used by the Python client. |
| `numcmc/`      | Affine-invariant ensemble MCMC sampler (`USE_MCMC`). |
| `pyazr/`       | Python package that drives AZURE2 over the socket API. |
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

## Dependencies

These are the only components you need to install yourself; everything else is
bundled in-tree.

**Build tools**
- A C++ compiler with OpenMP support (GCC or Clang)
- [CMake](https://cmake.org/) ≥ 3.16

**Libraries**
- [GSL](https://www.gnu.org/software/gsl/) — GNU Scientific Library (math routines)
- [Qt5](https://www.qt.io/) — `Core`, `Widgets`, `Svg`, `Script` (for the GUI; `BUILD_GUI`)
- [Qwt](https://qwt.sourceforge.io/) (Qt5 build) — for the in-app plotting tab (`USE_QWT`)
- [Readline](https://tiswww.case.edu/php/chet/readline/rltop.html) — for CLI input (`USE_READLINE`)

**For the Python client (`pyazr`)**
- Python 3 with [NumPy](https://numpy.org/)

### Installing dependencies

**Ubuntu / Debian**
```bash
sudo apt-get update
sudo apt-get install build-essential cmake libgsl-dev libreadline-dev \
    qtscript5-dev libqwt-qt5-dev libqt5svg5-dev qtwebengine5-dev \
    python3 python3-numpy
```

**macOS (Homebrew)**
```bash
brew install cmake gsl readline qt@5 qwt libomp
```

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
| `USE_API`               | ON      | Build the socket API server (required for `pyazr`). |
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

`pyazr` drives one or more headless AZURE2 processes over the socket API
(`AZURE2 --no-gui --use-api <port> <file>`), making it easy to script fits and
plug AZURE2 into external samplers such as [`emcee`](https://emcee.readthedocs.io/)
or [`brick`](https://github.com/odell/brick).

```python
from pyazr import azure2

# Spawn an AZURE2 instance bound to a configuration file.
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

Set `AZURE2_BINARY` (or pass `binary=...`) if the executable is not at
`build/src/AZURE2`.

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

Sources live in `docs/`; the generated output is not tracked.

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
