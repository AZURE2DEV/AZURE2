// The <config> section: its reader and the lines of an edited entry (ConfigSection.h).
#include "azure2/io/ConfigSection.h"

#include <algorithm>
#include <cstdlib>

#include "azure2/io/Document.h"

namespace azure2 {
namespace io {

namespace {

std::string Trim(const std::string &s, const char *blanks) {
  const std::size_t a = s.find_first_not_of(blanks);
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(blanks) - a + 1);
}

bool BlankLine(const std::string &s) { return s.find_first_not_of(" \t\r\n\f\v") == std::string::npos; }

// The first token of a line and where it lies (the stream's >> on that line).
bool FirstToken(const std::string &line, std::size_t &a, std::size_t &b) {
  a = line.find_first_not_of(" \t\r\n\f\v");
  if (a == std::string::npos) return false;
  b = line.find_first_of(" \t\r\n\f\v", a);
  if (b == std::string::npos) b = line.size();
  return true;
}

// The positional entries and how each is read: a token (>> then the rest of the line skipped) or a whole line.
struct Positional {
  const char *name;
  bool token;
};
const Positional kPositional[] = {
    {"amatrix", true},           {"output_dir", false},       {"checks_dir", false},
    {"check_compound_nucleus", true}, {"check_boundary_conditions", true}, {"check_data", true},
    {"check_energy_dependence", true}, {"check_legendre", true}, {"check_coulomb_amplitudes", true},
    {"check_pathways", true},    {"check_angular_distributions", true}};

// key=value of a line after the positional ones: false when it has no '=' before its comment.
bool KeyValue(const std::string &line, std::string &key, std::string &value) {
  std::string entry = line;
  const std::size_t hash = entry.find('#');
  if (hash != std::string::npos) entry = entry.substr(0, hash);
  const std::size_t eq = entry.find('=');
  if (eq == std::string::npos) return false;
  key = Trim(entry.substr(0, eq), " \t\r");
  value = Trim(entry.substr(eq + 1), " \t\r");
  return true;
}

// What follows a value that changed length by -`shrink` characters: when a comment follows after blanks, the comment
// keeps its column (blanks added, or removed down to one); otherwise unchanged.
std::string Refit(const std::string &rest, long shrink) {
  const std::size_t hash = rest.find('#');
  if (hash == std::string::npos || rest.find_first_not_of(" \t") != hash || shrink == 0) return rest;
  const long blanks = (long)hash + shrink;
  return std::string((std::size_t)std::max(1L, blanks), ' ') + rest.substr(hash);
}

const ConfigKey *KeyNamed(const std::string &name) {
  for (const ConfigKey &k : ConfigKeys())
    if (name == k.name) return &k;
  return nullptr;
}

bool OneLine(const std::string &value, std::string &why) {
  if (value.find_first_of("#\n\r") != std::string::npos) {
    why = "'" + value + "': a <config> value is one line without '#' (the line is cut there)";
    return false;
  }
  if (Trim(value, " \t") != value) {
    why = "'" + value + "': no blank before or after a <config> value (the engine trims them)";
    return false;
  }
  return true;
}

bool ValueTaken(const ConfigKey &k, const std::string &value, std::string &why) {
  if (!OneLine(value, why)) return false;
  char *end = nullptr;
  switch (k.kind) {
    case ConfigKey::Bool:
      if (value == "0" || value == "1" || value == "true" || value == "false") return true;
      why = std::string(k.name) + "=" + value + ": the value is 0, 1, true or false";
      return false;
    case ConfigKey::Int: {
      const long v = std::strtol(value.c_str(), &end, 10);
      if (!value.empty() && *end == '\0' && v >= k.minimum) return true;
      why = std::string(k.name) + "=" + value + ": the value is a whole number >= " + std::to_string(k.minimum);
      return false;
    }
    case ConfigKey::Number: {
      const double v = std::strtod(value.c_str(), &end);
      if (!value.empty() && *end == '\0' && v >= 0.0) return true;
      why = std::string(k.name) + "=" + value + ": the value is a number >= 0";
      return false;
    }
    case ConfigKey::Path:
      if (!value.empty()) return true;
      why = std::string(k.name) + ": the path is empty (remove the key instead)";
      return false;
    case ConfigKey::Version:
      why = std::string(k.name) + " is the file's format, written by the converter; it is not edited";
      return false;
  }
  return false;
}

}  // namespace

const std::vector<std::string> &ConfigPositionalNames() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> n;
    for (const Positional &p : kPositional) n.push_back(p.name);
    return n;
  }();
  return names;
}

const std::vector<ConfigKey> &ConfigKeys() {
  static const std::vector<ConfigKey> keys = {
      {"includeAllJPi", ConfigKey::Bool, 0, "every J^pi up to maxOrbitalMomentum gets a hard-sphere channel"},
      {"maxOrbitalMomentum", ConfigKey::Int, 0, "the largest l of the hard-sphere channels"},
      {"abundanceSumSigma", ConfigKey::Number, 0, "the width of the penalty on the sum of the abundances (0: none)"},
      {"relativistic", ConfigKey::Bool, 0, "relativistic two-body kinematics"},
      {"relativisticReferencePair", ConfigKey::Int, 1, "the pair whose masses define the relativistic frames"},
      {"priorCovariance", ConfigKey::Path, 0, "the prior parameter covariance file, relative to the project"},
      {"format", ConfigKey::Version, 0, "the project format (3)"}};
  return keys;
}

const ConfigEntry *ConfigSection::Find(const std::string &name) const {
  const ConfigEntry *found = nullptr;
  for (const ConfigEntry &e : entries)
    if (e.name == name) found = &e;
  return found;
}

ConfigSection ReadConfigSection(const std::vector<std::string> &lines) {
  ConfigSection c;
  std::size_t next = 1;  // line 0 is the <config> tag
  for (const Positional &p : kPositional) {
    ConfigEntry e;
    e.name = p.name;
    e.positional = true;
    if (p.token) {
      while (next < lines.size() && BlankLine(lines[next])) next++;
      std::size_t a = 0, b = 0;
      if (next >= lines.size() || !FirstToken(lines[next], a, b)) return c;  // the stream ran out: not closed
      e.value = lines[next].substr(a, b - a);
    } else {
      if (next >= lines.size()) return c;
      std::string temp = lines[next];
      const std::size_t hash = temp.find_last_of('#');
      if (hash != std::string::npos) temp = temp.substr(0, hash);
      e.value = Trim(temp, " \n\t\r");
    }
    e.line = (int)next++;
    c.entries.push_back(e);
  }
  for (; next < lines.size(); next++) {
    if (lines[next] == "</config>") {
      c.closed = true;
      break;
    }
    ConfigEntry e;
    if (!KeyValue(lines[next], e.name, e.value)) continue;
    e.line = (int)next;
    c.entries.push_back(e);
  }
  return c;
}

std::vector<std::string> ConfigLinesOf(const Document &doc, std::size_t &open) {
  std::vector<std::string> lines;
  const Document::Section *s = doc.Find("config");
  if (!s) return lines;
  open = s->open;
  const std::size_t end = s->close == Document::Section::npos ? doc.NumLines() : s->close + 1;
  for (std::size_t i = s->open; i < end; i++) {
    bool closing = false;
    const std::string tag = TagOf(doc.Line(i), closing);
    lines.push_back((i == s->open || (i == s->close && !tag.empty())) ? (closing ? "</" : "<") + tag + ">" : doc.Line(i));
  }
  return lines;
}

bool ConfigEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &value,
                 ConfigChange change, std::vector<ConfigLineEdit> &edits, std::string &why) {
  edits.clear();
  const ConfigSection c = ReadConfigSection(lines);
  if (!c.closed) {
    why = "the <config> section is not complete (eleven entries, then </config>)";
    return false;
  }
  int k = 0;
  for (const Positional &p : kPositional) {
    if (name != p.name) {
      k++;
      continue;
    }
    if (change == ConfigChange::Remove) {
      why = name + " is a positional entry of <config>: it is set, not removed";
      return false;
    }
    if (!OneLine(value, why)) return false;
    if (p.token) {
      const bool amatrix = k == 0;
      if (amatrix ? (value != "true" && value != "false") : (value != "none" && value != "screen" && value != "file")) {
        why = name + "=" + value + ": the value is " + (amatrix ? "true or false" : "none, screen or file");
        return false;
      }
    }
    const ConfigEntry &e = c.entries[k];
    const std::string &old = lines[e.line];
    // The span [a, b) of the value as read; a token entry's value is its first token, a directory's the text before
    // the last '#', trimmed (an empty value stands at the start of the line).
    std::size_t a = 0, b = 0;
    if (p.token) {
      FirstToken(old, a, b);
    } else {
      const std::size_t hash = old.find_last_of('#');
      const std::string before = hash == std::string::npos ? old : old.substr(0, hash);
      const std::size_t first = before.find_first_not_of(" \n\t\r");
      if (first != std::string::npos) {
        a = first;
        b = before.find_last_not_of(" \n\t\r") + 1;
      }
    }
    const std::string text = old.substr(0, a) + value + Refit(old.substr(b), (long)(b - a) - (long)value.size());
    edits.push_back({e.line, text, false, false});
    return true;
  }
  const ConfigKey *key = KeyNamed(name);
  if (!key) {
    std::string list;
    for (const std::string &n : ConfigPositionalNames()) list += (list.empty() ? "" : ", ") + n;
    for (const ConfigKey &x : ConfigKeys()) list += std::string(", ") + x.name;
    why = "'" + name + "' is not an entry of <config> (the entries: " + list + ")";
    return false;
  }
  if (change == ConfigChange::Set && !ValueTaken(*key, value, why)) return false;
  if (change == ConfigChange::Remove && key->kind == ConfigKey::Version) {
    why = name + " is the file's format, written by the converter; it is not removed";
    return false;
  }
  std::vector<int> own;
  for (const ConfigEntry &e : c.entries)
    if (!e.positional && e.name == name) own.push_back(e.line);
  int close = 0;
  for (std::size_t i = 1; i < lines.size(); i++)
    if (lines[i] == "</config>") close = (int)i;
  if (change == ConfigChange::Remove) {
    for (int line : own) edits.push_back({line, "", true, false});
    return true;
  }
  if (own.empty()) {
    edits.push_back({close, name + "=" + value, false, true});
    return true;
  }
  for (std::size_t i = 0; i + 1 < own.size(); i++) edits.push_back({own[i], "", true, false});
  // The value's span in the winning line: after '=', trimmed, before the comment.
  const std::string &old = lines[own.back()];
  const std::size_t eq = old.find('='), end = std::min(old.find('#'), old.size());
  std::size_t a = old.find_first_not_of(" \t\r", eq + 1), b = a;
  if (a == std::string::npos || a > end) {
    a = b = eq + 1;
  } else {
    b = old.find_last_not_of(" \t\r", end - 1) + 1;
  }
  edits.push_back({own.back(), old.substr(0, a) + value + Refit(old.substr(b), (long)(b - a) - (long)value.size()), false,
                   false});
  return true;
}

}  // namespace io
}  // namespace azure2
