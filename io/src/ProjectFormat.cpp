// The format of a project file and the refusal of one that is not format 3 (ProjectFormat.h).
#include "azure2/io/ProjectFormat.h"

#include <vector>

namespace azure2 {
namespace io {

namespace {

std::string Trim(const std::string &s) {
  std::size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t')) a++;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) b--;
  return s.substr(a, b - a);
}

/// The sections of formats 1 and 2 that format 3 replaces (project_format3.md section 3.2), in the order a refusal
/// names them.
const std::vector<std::string> kOldSections = {"segmentsData", "segmentsTest", "targetInt", "abundances", "scatteringLengths"};

}  // namespace

bool ConfigKeyValue(const std::string &line, std::string &key, std::string &value) {
  std::string entry = line.substr(0, line.find('#'));
  const std::size_t eq = entry.find('=');
  if (eq == std::string::npos) return false;
  key = Trim(entry.substr(0, eq));
  value = Trim(entry.substr(eq + 1));
  return true;
}

int FormatOf(const Document &doc) {
  const Document::Section *c = doc.Find("config");
  if (!c) return 1;
  for (const std::string &line : doc.Body(*c)) {
    std::string key, value;
    if (ConfigKeyValue(line, key, value) && key == "format") return value == "2" ? 2 : value == "3" ? 3 : 1;
  }
  return 1;
}

std::string NotFormat3Message(const std::string &text, const std::string &fileName) {
  const Document doc = Document::Parse(text);
  // The name the conversion is offered under: the file's stem with "_f3" (T2.azr -> T2_f3.azr).
  const std::size_t slash = fileName.find_last_of("/\\"), dot = fileName.rfind('.');
  const std::string stem = (dot == std::string::npos || (slash != std::string::npos && dot < slash)) ? fileName : fileName.substr(0, dot);
  const std::string converted = stem + "_f3.azr";
  const std::string howTo =
      ". This AZURE2 reads format 3 only.\n"
      "       Convert it:  AZURE2 --convert " + fileName + " " + converted + "     (or pyazr.convert_azr(\"" + fileName +
      "\", \"" + converted + "\"))\n"
      "       The converter keeps " + fileName + " unchanged and reports every change it makes.";
  if (!doc.Find("config"))
    return fileName + " is not an AZURE2 project: it has no <config> section" +
           (text.find_first_not_of(" \t\r\n") == std::string::npos ? std::string(" (the file is empty)") : std::string()) +
           ". This AZURE2 reads format-3 projects.";
  if (FormatOf(doc) == 3) {
    for (const std::string &name : kOldSections)
      if (const Document::Section *s = doc.Find(name))
        return fileName + " says format=3 but holds <" + name + "> at line " + std::to_string(s->open + 1) +
               ", a section of format 2. This AZURE2 reads format 3 only.\n"
               "       Convert the original (format 1 or 2) instead:  AZURE2 --convert ORIGINAL.azr NEW.azr\n"
               "       A format-3 file states its data in <data>, <extrapolations>, <targets> and <modelData> records.";
    return "";
  }
  std::string where;
  for (const std::string &name : kOldSections)
    if (const Document::Section *s = doc.Find(name)) {
      where = " (<" + name + "> at line " + std::to_string(s->open + 1) + ")";
      break;
    }
  return fileName + " is a format-" + std::to_string(FormatOf(doc)) + " project" + where + howTo;
}

}  // namespace io
}  // namespace azure2
