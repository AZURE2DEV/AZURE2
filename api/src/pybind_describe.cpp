// What pyazr.describe() and pyazr.kinds() read of the build and of the engine's tables (pyazr/described.py): the
// build options, the Coulomb table a session would use, and the table of session options.
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "CoulombTable.h"
#include "azure2/io/ConfigSection.h"

namespace py = pybind11;

namespace {

// The options a session takes that are not a property of the project, each with the reason (the fields of
// RuntimeOptions in pybind_api.cpp besides `options`).  The options of the calculation are not listed here: they are the
// option rows of io::ConfigKeys() (<config> keys, docs/upgrade/project_format3.md section 3.12), and
// session_options() below lists both.  tests/pyazr/describe_test.py fails for a field of RuntimeOptions in neither.
struct NotProjectOptionRow {
  const char *name;
  const char *keyword;
  const char *flag;
  const char *type;  ///< "bool" (default true) or "str" (default "")
  const char *doc;
  const char *reason;
};
const std::vector<NotProjectOptionRow> &NotProjectOptionRows() {
  static const std::vector<NotProjectOptionRow> rows = {
      {"data_mode", "data_mode", "", "bool",
       "evaluate the data sets (True) or the extrapolations (False); the program's menu choice 1 or 3",
       "what a run evaluates, chosen per run (the menu, pyazr's set_data / set_extrap)"},
      {"coulomb_table", "coulomb_table", "--coulomb-table=FILE", "str",
       "the Coulomb table file (\"\": $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin)",
       "a file of this machine; whether the table is used is the key use_coulomb_table"},
      {"base_directory", "cwd", "", "str",
       "the directory the project's relative paths are resolved in (\"\": the project file's directory)",
       "a directory of this machine"},
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
        for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
          if (!k.option) continue;
          py::dict d;
          d["name"] = k.name;
          d["keyword"] = k.name;
          d["key"] = k.name;
          d["flag"] = k.flag;
          d["negation"] = k.negation;
          d["type"] = "bool";
          d["default"] = k.defaultOn;
          d["effect"] = k.effect;
          d["reason"] = "";
          d["doc"] = k.doc;
          out.append(d);
        }
        for (const NotProjectOptionRow &row : NotProjectOptionRows()) {
          py::dict d;
          d["name"] = row.name;
          d["keyword"] = row.keyword;
          d["key"] = "";
          d["flag"] = row.flag;
          d["negation"] = "";
          d["type"] = row.type;
          if (std::string(row.type) == "bool")
            d["default"] = true;
          else
            d["default"] = "";
          d["effect"] = "";
          d["reason"] = row.reason;
          d["doc"] = row.doc;
          out.append(d);
        }
        return out;
      },
      "The session options, one dict each: name (in the engine), keyword (of pyazr.azure2), key (of <config>; \"\" "
      "for an option that is not project state, with its reason), flag and negation (the program's, \"\" where it "
      "has none), type, default, effect (physics, convention, numerics), doc.");
}
