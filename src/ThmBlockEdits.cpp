// The writer of the <thm> block beside its reader (ThmBlockEdits.h).
#include "ThmBlockEdits.h"
#include "ThmKeys.h"

#include <algorithm>

#include "ThmExperiment.h"
#include "azure2/experiment/ThmSectionLines.h"

namespace {

std::string Trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

// The line before its comment, trimmed (as the reader takes it).
std::string Entry(const std::string &line) { return Trim(line.substr(0, line.find('#'))); }

bool IsExperiment(const std::string &entry) { return azure2::experiment::ReadThmSectionLine(entry).experiment; }

std::string ExperimentName(const std::string &entry) { return azure2::experiment::ReadThmSectionLine(entry).name; }

// The names of the experiments that have an active data set (a member).
std::vector<std::string> MemberExperiments(const azure2::experiment::ThmMembershipSpec &thm) {
  std::vector<std::string> names;
  for (const auto &x : thm.experiments)
    if (!x.datasetKeys.empty()) names.push_back(x.name);
  return names;
}

// The key of a global line (text before '=', trimmed), "" for a line that is not one.
std::string OptionKey(const std::string &entry) {
  const std::size_t eq = entry.find('=');
  if (entry.empty() || IsExperiment(entry) || eq == std::string::npos) return std::string();
  return Trim(entry.substr(0, eq));
}

// The tokens key=value of an experiment line after "experiment[<name>]": (first byte, one past the last) in the line.
std::vector<std::pair<std::size_t, std::size_t>> KeyTokens(const std::string &line) {
  std::vector<std::pair<std::size_t, std::size_t>> out;
  const std::size_t end = std::min(line.find('#'), line.size());
  std::size_t i = line.find(']');
  if (i == std::string::npos) return out;
  for (i++; i < end;) {
    while (i < end && (line[i] == ' ' || line[i] == '\t')) i++;
    const std::size_t a = i;
    while (i < end && line[i] != ' ' && line[i] != '\t') i++;
    if (i > a) out.push_back({a, i});
  }
  return out;
}

}  // namespace

std::vector<ThmEntry> ThmBlockEntries(const std::vector<std::string> &lines) {
  std::vector<ThmEntry> out;
  for (std::size_t i = 0; i < lines.size(); i++) {
    const std::string entry = Entry(lines[i]);
    if (IsExperiment(entry)) {
      const std::string name = ExperimentName(entry);
      for (const auto &t : KeyTokens(lines[i])) {
        const std::string token = lines[i].substr(t.first, t.second - t.first);
        const std::size_t eq = token.find('=');
        out.push_back({name, token.substr(0, eq), eq == std::string::npos ? std::string() : token.substr(eq + 1), (int)i});
      }
      continue;
    }
    const std::string key = OptionKey(entry);
    if (!key.empty()) out.push_back({"", key, Trim(entry.substr(entry.find('=') + 1)), (int)i});
  }
  return out;
}

bool ThmOptionEdits(const std::vector<std::string> &lines, const std::string &key, const std::string &value,
                    ThmChange change, std::vector<ThmLineEdit> &edits, std::string &why) {
  edits.clear();
  ThmOptionLine parsed;
  if (key.compare(0, 6, "weight") == 0) {
    why = key + ": in format 3 a THM weight is the thm_weight field of its data set or extrapolation; set that "
                "field of the record";
    return false;
  }
  if (change == ThmChange::Set && !ParseThmOptionLine(key, value, parsed)) {
    why = "<thm> line not understood: '" + key + "=" + value + "' (the options: " + ThmOptionUsage() +
          "); give one of them (pyazr.kinds()[\"thm\"] lists them)";
    return false;
  }
  std::vector<int> own;
  int lastGlobal = -1;
  for (std::size_t i = 0; i < lines.size(); i++) {
    const std::string k = OptionKey(Entry(lines[i]));
    if (k.empty()) continue;
    lastGlobal = (int)i;
    if (k == key) own.push_back((int)i);
  }
  if (change == ThmChange::Remove) {
    if (own.empty()) {
      why = "<thm> has no line " + key + "=: only a line it has can be removed; Project.thm() lists them";
      return false;
    }
    for (auto it = own.rbegin(); it != own.rend(); ++it) edits.push_back({ThmLineEdit::Kind::Remove, *it, ""});
    return true;
  }
  if (own.empty()) {
    edits.push_back({ThmLineEdit::Kind::Insert, lastGlobal + 1, key + "=" + value});
    return true;
  }
  // The value's span in the winning line (after '=', trimmed, before the comment) takes the new value.
  std::string line = lines[own.back()];
  const std::size_t eq = line.find('='), end = std::min(line.find('#'), line.size());
  std::size_t a = line.find_first_not_of(" \t\r", eq + 1), b = a;
  if (a == std::string::npos || a >= end)
    a = b = eq + 1;
  else
    b = line.find_last_not_of(" \t\r", end - 1) + 1;
  line.replace(a, b - a, value);
  edits.push_back({ThmLineEdit::Kind::Replace, own.back(), line});
  for (std::size_t i = own.size() - 1; i-- > 0;) edits.push_back({ThmLineEdit::Kind::Remove, own[i], ""});
  return true;
}

bool ThmExperimentEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &key,
                        const std::string &value, ThmChange change, std::vector<ThmLineEdit> &edits, std::string &why) {
  edits.clear();
  if (key == "segments") {
    why = "segments= is not written in format 3: a data set names its experiment (thm " + name +
          "); set the field thm of the data sets instead";
    return false;
  }
  if (change == ThmChange::Set) {
    // The key and value as the reader takes them, on a line of their own; the reader's message without the
    // experiment it names first (the caller names it).
    std::vector<ThmExperiment> probe;
    why = ParseThmExperimentLine("experiment[" + name + "] " + key + "=" + value, probe);
    if (!why.empty()) {
      const std::string own = "experiment[" + name + "]: ";
      if (why.compare(0, own.size(), own) == 0) why = why.substr(own.size());
      return false;
    }
  }
  std::vector<int> own;
  for (std::size_t i = 0; i < lines.size(); i++) {
    const std::string entry = Entry(lines[i]);
    if (IsExperiment(entry) && ExperimentName(entry) == name) own.push_back((int)i);
  }
  if (change == ThmChange::Remove && key.empty()) {
    if (own.empty()) {
      why = "<thm> has no experiment[" + name + "]: only an experiment it has can be removed; Project.thm() lists "
            "them";
      return false;
    }
    for (auto it = own.rbegin(); it != own.rend(); ++it) edits.push_back({ThmLineEdit::Kind::Remove, *it, ""});
    return true;
  }
  // The line and token of the key (the last one given).
  int at = -1;
  std::pair<std::size_t, std::size_t> span;
  for (int i : own)
    for (const auto &t : KeyTokens(lines[i]))
      if (lines[i].compare(t.first, key.size() + 1, key + "=") == 0) {
        at = i;
        span = t;
      }
  if (change == ThmChange::Remove) {
    if (at < 0) {
      why = "experiment[" + name + "] has no key " + key + ": only a key it has can be removed; Project.thm() "
            "lists them";
      return false;
    }
    std::string line = lines[at];
    std::size_t a = span.first;
    while (a > 0 && (line[a - 1] == ' ' || line[a - 1] == '\t')) a--;
    line.erase(a, span.second - a);
    edits.push_back({ThmLineEdit::Kind::Replace, at, line});
    return true;
  }
  const std::string token = key + "=" + value;
  if (at >= 0) {
    std::string line = lines[at];
    line.replace(span.first, span.second - span.first, token);
    edits.push_back({ThmLineEdit::Kind::Replace, at, line});
  } else if (!own.empty()) {
    std::string line = lines[own.back()];
    const std::size_t hash = line.find('#');
    const std::size_t end = (hash == std::string::npos ? line : line.substr(0, hash)).find_last_not_of(" \t\r") + 1;
    line.insert(end, " " + token);
    edits.push_back({ThmLineEdit::Kind::Replace, own.back(), line});
  } else {
    edits.push_back({ThmLineEdit::Kind::Insert, (int)lines.size(), "experiment[" + name + "] " + token});
  }
  return true;
}

std::vector<std::string> ApplyThmLineEdits(const std::vector<std::string> &lines, const std::vector<ThmLineEdit> &edits) {
  std::vector<std::string> out = lines;
  for (const ThmLineEdit &e : edits) {
    if (e.kind == ThmLineEdit::Kind::Replace)
      out[e.line] = e.text;
    else if (e.kind == ThmLineEdit::Kind::Remove)
      out.erase(out.begin() + e.line);
    else
      out.insert(out.begin() + e.line, e.text);
  }
  return out;
}

std::vector<std::string> ThmExperimentsWithoutMembers(const std::vector<std::string> &section,
                                                      const azure2::experiment::ThmMembershipSpec &thm) {
  const std::vector<std::string> members = MemberExperiments(thm);
  std::vector<std::string> names;
  for (const auto &x : thm.experiments)
    if (x.datasetKeys.empty()) names.push_back(x.name);
  for (const std::string &line : section) {
    const azure2::experiment::ThmSectionLine l = azure2::experiment::ReadThmSectionLine(line);
    if (!l.experiment || l.text.find(']') == std::string::npos) continue;
    if (std::find(members.begin(), members.end(), l.name) == members.end() &&
        std::find(names.begin(), names.end(), l.name) == names.end())
      names.push_back(l.name);
  }
  return names;
}

std::string CheckFormat3ThmMembers(const std::vector<std::string> &lines, const azure2::experiment::ThmMembershipSpec &thm) {
  for (const std::string &name : MemberExperiments(thm)) {
    bool found = false;
    for (std::size_t i = 0; i < lines.size() && !found; i++) {
      const azure2::experiment::ThmSectionLine l = azure2::experiment::ReadThmSectionLine(lines[i]);
      found = l.experiment && l.name == name;
    }
    if (!found)
      return "internal error: <thm> has no experiment[" + name + "], which active data sets name (thm " + name +
             "); the project's reading should have refused it (E_THM_EXPERIMENT_UNKNOWN)";
  }
  return "";
}
