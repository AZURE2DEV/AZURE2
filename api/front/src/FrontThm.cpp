// azure2::front's <thm> (Engine.h; docs/upgrade/qt_gui.md section 8): the tables of its keys (ThmKeys.h, read by the
// engine's parser) and the engine's edits of the section's lines (ThmBlockEdits.h), as pyazr's Project edits them
// (api/src/pybind_nuclear_part.cpp: thm_keys, thm_entries, thm_option_edits, thm_experiment_edits).
#include "azure2/front/Engine.h"

#include <sstream>

#include "ThmBlockEdits.h"
#include "ThmKeys.h"

namespace azure2 {
namespace front {

namespace {

FieldInfo FieldOf(const ThmKeyRow &r) {
  FieldInfo f;
  f.key = r.key;
  f.title = r.title.empty() ? r.key : r.title;
  f.type = r.type;
  f.grammar = r.grammar;
  f.doc = r.doc;
  f.defaultText = r.defaultValue;
  f.values = r.values;
  if (!r.unit.empty()) f.units = {r.unit};
  if (!r.presentWith.empty()) f.appliesWhen = r.presentWith + " is given";
  return f;
}

/// The first line of <thm> (1-based), 0 when the file has none.
int ThmLine(const Project &project) {
  std::istringstream text(project.Text());
  std::string line;
  int n = 0;
  while (std::getline(text, line)) {
    n++;
    const std::size_t a = line.find_first_not_of(" \t");
    if (a != std::string::npos && line.compare(a, 5, "<thm>") == 0) return n;
  }
  return 0;
}

}  // namespace

std::vector<KindInfo> ThmKinds() {
  std::vector<KindInfo> out(2);
  out[0].kind = ThmOptionKind().name;
  out[0].title = ThmOptionKind().title;
  out[0].doc = ThmOptionKind().doc;
  out[0].sections = {"thm"};
  for (const ThmOptionKey &k : ThmOptionKeys()) out[0].fields.push_back(FieldOf(k.row));
  out[1].kind = ThmExperimentKind().name;
  out[1].title = ThmExperimentKind().title;
  out[1].doc = ThmExperimentKind().doc;
  out[1].sections = {"thm"};
  for (const ThmExperimentKey &k : ThmExperimentKeys()) out[1].fields.push_back(FieldOf(k.row));
  return out;
}

std::vector<ThmEntry> ThmEntries(const Project &project) {
  std::vector<ThmEntry> out;
  const int first = ThmLine(project);
  for (const ::ThmEntry &e : ThmBlockEntries(project.SectionBody("thm")))
    out.push_back({e.experiment, e.key, e.value, first ? first + 1 + e.line : 0});
  return out;
}

bool SetThmOption(Project &project, const std::string &key, const std::string &value, std::string &why, bool remove) {
  const std::vector<std::string> body = project.SectionBody("thm");
  std::vector<ThmLineEdit> edits;
  if (!ThmOptionEdits(body, key, value, remove ? ThmChange::Remove : ThmChange::Set, edits, why)) return false;
  project.SetSectionBody("thm", ApplyThmLineEdits(body, edits));
  return true;
}

bool SetThmExperimentKey(Project &project, const std::string &name, const std::string &key, const std::string &value,
                         std::string &why, bool remove) {
  const std::vector<std::string> body = project.SectionBody("thm");
  std::vector<ThmLineEdit> edits;
  if (!ThmExperimentEdits(body, name, key, value, remove ? ThmChange::Remove : ThmChange::Set, edits, why)) return false;
  project.SetSectionBody("thm", ApplyThmLineEdits(body, edits));
  return true;
}

}  // namespace front
}  // namespace azure2
