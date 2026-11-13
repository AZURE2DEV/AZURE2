// What pyazr.describe() and pyazr.kinds() read of the build and of the engine's tables (pyazr/described.py): the
// build options, the Coulomb table a session would use, and the table of session options.
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "CoulombTable.h"

namespace py = pybind11;

namespace {

// One row per field of RuntimeOptions (pybind_api.cpp): what pyazr.kinds()["options"] says of it (its type and
// default come from a default RuntimeOptions).  `keyword` is the argument of pyazr.azure2; `flag` and
// `negation` are the program's command-line flags that set it on and off ("" where the program has none).  A new
// field takes a row here; tests/pyazr/describe_test.py fails until it has one.
struct SessionOptionRow {
  const char *name;
  const char *keyword;
  const char *flag;
  const char *negation;
  const char *doc;
};
const std::vector<SessionOptionRow> &SessionOptionRows() {
  static const std::vector<SessionOptionRow> rows = {
      {"data_mode", "data_mode", "", "",
       "evaluate the data sets (True) or the extrapolations (False); the program's menu choice 1 or 3"},
      {"use_brune", "use_brune", "--use-brune", "--no-brune",
       "Brune's alternative level matrix (on-shell boundary conditions); off: Lane-Thomas constant boundary "
       "conditions"},
      {"ignore_externals", "ignore_externals", "--ignore-externals", "--no-ignore-externals",
       "ignore the external capture of a photon channel whose input width is exactly zero"},
      {"transform", "transform", "", "--no-transform", "transform the input parameters to the run's basis at start"},
      {"use_long_wavelength", "use_long_wavelength", "", "--no-long-wavelength",
       "the long-wavelength approximation for EL capture"},
      {"use_gsl_coul", "use_gsl_coul", "--gsl-coul", "", "GSL's Coulomb functions instead of AZURE2's own"},
      {"use_rmc", "use_rmc", "--use-rmc", "", "the Reich-Moore approximation for capture (neutron capture only)"},
      {"use_park", "use_park", "--use-park", "",
       "fit the observed reduced width amplitudes (T.-S. Park, PRC 104, 064612); implies use_brune"},
      {"use_coulomb_table", "use_coulomb_table", "--coulomb-table", "--no-coulomb-table",
       "P, S and the hard-sphere phase of charged channels from the universal Coulomb table, exact outside it"},
      {"coulomb_table", "coulomb_table", "--coulomb-table=FILE", "",
       "the Coulomb table file (\"\": $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin)"},
      {"replace_grids", "replace_grids", "--grid-replacement", "--no-grid-replacement",
       "lay the adaptive sub-point grids out again as levels and energy shifts move"},
      {"use_ec_table", "use_ec_table", "--ec-table", "--no-ec-table",
       "external-capture integrals from per-term energy tables, interpolated"},
      {"ec_convection_current", "ec_convection_current", "--ec-convection-current", "--no-ec-convection-current",
       "with use_long_wavelength off: external capture with the convection-current part of the EL operator too"},
      {"base_directory", "cwd", "", "",
       "the directory the project's relative paths are resolved in (\"\": the project file's directory)"},
  };
  return rows;
}

}  // namespace

// What pyazr.describe() and pyazr.kinds() read of the build and of the engine's tables.
void BindDescription(py::module_ &m) {
  m.def(
      "build_info",
      []() {
        py::dict d;
#ifdef AZURE2_VERSION_STRING
        d["version"] = AZURE2_VERSION_STRING;
#else
        d["version"] = "unknown";
#endif
#ifdef AZURE2_GIT_VERSION
        d["git"] = AZURE2_GIT_VERSION;
#else
        d["git"] = "";
#endif
#ifdef USE_MCMC
        d["mcmc"] = true;
#else
        d["mcmc"] = false;
#endif
#ifdef USE_NLOPT
        d["nlopt"] = true;
#else
        d["nlopt"] = false;
#endif
#ifdef NO_READLINE
        d["readline"] = false;
#else
        d["readline"] = true;
#endif
#ifdef _OPENMP
        d["openmp"] = true;
        d["max_threads"] = omp_get_max_threads();
#else
        d["openmp"] = false;
        d["max_threads"] = 1;
#endif
        return d;
      },
      "The build: version, git (git describe at CMake configure time), mcmc, nlopt, readline (the program's), "
      "openmp, max_threads (omp_get_max_threads(): the threads an evaluation uses now).");
  m.def(
      "coulomb_table_info",
      [](const std::string &path) {
        const std::string file = path.empty() ? CoulombTable::DefaultPath() : path;
        const CoulombTable *table = nullptr;
        std::string identity;
        {
          py::gil_scoped_release release;
          table = CoulombTable::Shared(path);  // loaded once per process, as a session loads it
          if (table) identity = table->ContentIdentity();
        }
        py::dict d;
        d["path"] = file;
        d["loaded"] = table != nullptr;
        d["identity"] = identity;
        return d;
      },
      py::arg("path") = "",
      "The Coulomb table a session with this coulomb_table option uses: path (\"\" = the default file), loaded "
      "(False: it cannot be read, and the exact routines are used), identity (CoulombTable::ContentIdentity).");
  m.def(
      "session_options",
      []() {
        py::list out;
        for (const SessionOptionRow &row : SessionOptionRows()) {
          py::dict d;
          d["name"] = row.name;
          d["keyword"] = row.keyword;
          d["flag"] = row.flag;
          d["negation"] = row.negation;
          d["doc"] = row.doc;
          out.append(d);
        }
        return out;
      },
      "The session options (fields of RuntimeOptions), one dict each: name, keyword (of pyazr.azure2), flag and "
      "negation (the program's command-line flags, \"\" where it has none), doc.");
}
