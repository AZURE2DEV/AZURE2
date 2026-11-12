// The project-file commands of the command line (ProjectCommands.h): --convert (an old project to format 3,
// docs/upgrade/project_format3.md section 8), --convert-params (a parameter file of an old project, section 8.4), --validate (section 6.2), --describe (section 6.3), --schema
// (section 6.1) and --print-predictions (docs/upgrade/data_model.md section 3).
#include "ProjectCommands.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Predictions.h"
#include "TargetStopping.h"
#include "azure2/convert/Convert.h"
#include "azure2/experiment/Describe.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/ProjectFormat.h"
#include "azure2/experiment/Schema.h"

namespace {

/// The whole file as bytes; false when it cannot be read.
bool ReadFileBytes(const std::string &path, std::string &bytes) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  std::ostringstream s;
  s << in.rdbuf();
  bytes = s.str();
  return true;
}

/// Writes a new file; refuses one that exists (the converter never overwrites, docs/upgrade/project_format3.md 8).
bool WriteNewFile(const std::string &path, const std::string &bytes) {
  if (std::ifstream(path)) {
    std::cout << "ERROR: --convert: " << path << " exists; it is not overwritten" << std::endl;
    return false;
  }
  std::ofstream out(path, std::ios::binary);
  out << bytes;
  if (!out) {
    std::cout << "ERROR: --convert: cannot write " << path << std::endl;
    return false;
  }
  return true;
}

/*!
 * --convert OLD NEW [--to 3] [--json-report R] [--answers A] [--exact-grids] [--strict] [--allow-questions]: OLD in
 * format 3 (convert/, docs/upgrade/project_format3.md section 8) written to NEW, with its side files (<f>.f3), and the
 * report; NEW, the side files and R must not exist.  Exit 0 when written, 1 when not (the report says why).
 */
int convertProjectTo3(const std::string &from, const std::string &to, const std::vector<std::string> &options) {
  azure2::convert::ConvertOptions o;
  std::string reportFile, answersFile;
  for (std::size_t i = 0; i < options.size(); i++) {
    const std::string &a = options[i];
    if (a == "--exact-grids") o.exactGrids = true;
    else if (a == "--strict") o.strict = true;
    else if (a == "--allow-questions") o.allowQuestions = true;
    else if (a == "--json-report" && i + 1 < options.size()) reportFile = options[++i];
    else if (a == "--answers" && i + 1 < options.size()) answersFile = options[++i];
    else {
      std::cout << "ERROR: --convert: unknown option " << a << " (see --help)" << std::endl;
      return 1;
    }
  }
  if (from == to) {
    std::cout << "ERROR: --convert: the input and the output are the same file" << std::endl;
    return 1;
  }
  std::string bytes;
  if (!ReadFileBytes(from, bytes)) {
    std::cout << "ERROR: --convert: cannot read " << from << std::endl;
    return 1;
  }
  if (!answersFile.empty()) {
    std::string answers;
    if (!ReadFileBytes(answersFile, answers)) {
      std::cout << "ERROR: --convert: cannot read the answers " << answersFile << std::endl;
      return 1;
    }
    const std::string why = azure2::convert::ReadAnswersJson(answers, o.answers);
    if (!why.empty()) {
      std::cout << "ERROR: --convert: " << answersFile << ": " << why << std::endl;
      return 1;
    }
  }
  const std::size_t slash = from.find_last_of("/\\");
  o.projectDir = slash == std::string::npos ? "." : from.substr(0, slash);
  const azure2::convert::ConvertResult r = azure2::convert::ConvertToFormat3(bytes, o);
  if (!reportFile.empty() && !WriteNewFile(reportFile, r.ReportJson(from, to))) return 1;
  for (const azure2::convert::ConvertNote &n : r.notes)
    std::cout << "  note: " << (n.record.empty() ? "" : n.record + ": ") << n.text << std::endl;
  for (const azure2::convert::ConvertQuestion &q : r.questions)
    std::cout << "  question " << q.key << ": " << q.text << (q.proposal.empty() ? "" : " [" + q.proposal + "]") << std::endl;
  if (!r.ok) {
    for (const azure2::experiment::Diagnostic &d : r.diagnostics) std::cout << "  " << d.ToText() << std::endl;
    std::cout << "ERROR: --convert: " << from << " is not converted: " << r.refusal << std::endl;
    return 1;
  }
  const std::size_t outSlash = to.find_last_of("/\\");
  const std::string outDir = outSlash == std::string::npos ? "" : to.substr(0, outSlash + 1);
  for (const auto &f : r.sideFiles)
    if (!WriteNewFile(outDir + f.first, f.second)) return 1;
  if (!WriteNewFile(to, r.text)) return 1;
  std::cout << "Wrote " << to << " (format 3): " << r.notes.size() << " note(s), " << r.questions.size() << " question(s)"
            << (reportFile.empty() ? "" : "; report " + reportFile) << "." << std::endl;
  return 0;
}

/*!
 * --convert-params OLD.sav NEW.sav --project OLD.azr: a parameter file written for the old project OLD.azr with the names
 * of OLD.azr's format-3 conversion (convert::ConvertParameterFile; docs/upgrade/project_format3.md section 8.4), written
 * to NEW.sav, which must not exist.  Exit 0 when written, 1 when not.
 */
int convertParameters(const std::string &from, const std::string &to, const std::string &project) {
  std::string bytes, projectBytes;
  if (!ReadFileBytes(from, bytes)) {
    std::cout << "ERROR: --convert-params: cannot read " << from << std::endl;
    return 1;
  }
  if (!ReadFileBytes(project, projectBytes)) {
    std::cout << "ERROR: --convert-params: cannot read the project " << project << std::endl;
    return 1;
  }
  const std::size_t slash = project.find_last_of("/\\");
  const azure2::convert::ParameterFileConversion c =
      azure2::convert::ConvertParameterFile(bytes, projectBytes, slash == std::string::npos ? "." : project.substr(0, slash));
  if (!c.refusal.empty()) {
    std::cout << "ERROR: --convert-params: " << project << ": " << c.refusal << std::endl;
    return 1;
  }
  if (!WriteNewFile(to, c.text)) return 1;
  for (const auto &r : c.renamed) std::cout << "  " << r.first << " -> " << r.second << std::endl;
  std::cout << "Wrote " << to << ": " << c.renamed.size() << " name(s) renamed for the format-3 conversion of " << project
            << "." << std::endl;
  return 0;
}

/*!
 * --validate FILE [--json]: the format-3 project FILE checked against the table (docs/upgrade/project_format3.md
 * section 6.2).  Every diagnostic is printed, one per line (file:line: record.field: message) or as one JSON array;
 * exit 0 without errors, 2 with errors, 3 for a file that is not a format-3 project (section 8.4), 1 when FILE cannot be read.
 */
/// The target effects of a built project, data sets' and extrapolations'.
std::vector<azure2::experiment::TargetEffectSpec> TargetEffectSpecs(const azure2::experiment::Format3Project &c) {
  std::vector<azure2::experiment::TargetEffectSpec> out;
  for (const auto &item : c.built.data.targetEffects.items) out.push_back(item.spec);
  for (const auto &item : c.built.extrapolations.targetEffects.items) out.push_back(item.spec);
  return out;
}

int validateProject(const std::string &file, bool json) {
  std::string bytes;
  if (!ReadFileBytes(file, bytes)) {
    std::cout << "ERROR: --validate: cannot read " << file << std::endl;
    return 1;
  }
  const std::string refusal = azure2::io::NotFormat3Message(bytes, file);
  if (!refusal.empty()) {
    std::cout << "ERROR: " << refusal << std::endl;
    return 3;
  }
  azure2::experiment::Format3ReadOptions read;
  read.fileName = file;
        // data files as the engine opens them: relative to the working directory
  azure2::experiment::Format3Project c = azure2::experiment::ReadFormat3Text(bytes, read);
  // the stopping of each target from the library (energy_grid.md section 6): its source, or why it has none
  const StoppingReport stopping = ReportTargetStoppings(c.Ok() ? TargetEffectSpecs(c) : std::vector<azure2::experiment::TargetEffectSpec>());
  for (const azure2::experiment::Diagnostic &d : stopping.refusals) c.diagnostics.push_back(d);
  const std::string note = "note: " + std::to_string(c.buildChecksSkipped) +
                           " check(s) skipped (the build checks need a valid project and its <levels>).";
  if (json) {
    std::cout << "[";
    for (std::size_t i = 0; i < c.diagnostics.size(); i++) std::cout << (i ? ",\n " : "\n ") << c.diagnostics[i].ToJson();
    std::cout << (c.diagnostics.empty() ? "]" : "\n]") << std::endl;
    std::cerr << note << std::endl;
  } else {
    for (const azure2::experiment::Diagnostic &d : c.diagnostics) std::cout << d.ToText() << std::endl;
    std::cout << stopping.text << note << std::endl;
    std::cout << file << ": " << c.diagnostics.size() << " diagnostic(s), " << (c.Ok() ? "valid" : "not valid") << "." << std::endl;
  }
  return c.Ok() ? 0 : 2;
}

/// The expression tree of every data set and extrapolation of `file`, one line each; 0, or 1 (with the set-up's
/// messages) when the project could not be set up or read, 3 for a file that is not a format-3 project (section 8.4).
int printPredictions(const std::string &file) {
  std::string bytes;
  if (!ReadFileBytes(file, bytes)) {
    std::cout << "ERROR: --print-predictions: cannot read " << file << std::endl;
    return 1;
  }
  const std::string refusal = azure2::io::NotFormat3Message(bytes, file);
  if (!refusal.empty()) {
    std::cout << "ERROR: " << refusal << std::endl;
    return 3;
  }
  std::ostringstream messages;
  const std::vector<std::string> lines = azure2::ProjectPredictions(file, messages);
  bool failed = false;
  for (const std::string &line : lines) {
    std::cout << line << std::endl;
    failed = failed || line.find(": set-up failed") != std::string::npos || line.rfind("project: ", 0) == 0;
  }
  if (failed) std::cerr << messages.str();
  return failed ? 1 : 0;
}

/*!
 * --describe FILE [--json] [--record NAME]: what the engine makes of the format-3 project FILE (docs/upgrade/
 * project_format3.md section 6.3; experiment/Describe.h): per data set its quantity, observable, reaction with the pair
 * keys, windows, the data file's resolved path, effects and parameters; with --record only the data set, extrapolation,
 * datum or potential NAME.  A project with an error, or no such record: exit 2; a file that is not a format-3 project: exit 3; a
 * file that cannot be read: exit 1.
 */
int describeProject(const std::string &file, bool json, const std::string &record) {
  std::string bytes;
  if (!ReadFileBytes(file, bytes)) {
    std::cout << "ERROR: --describe: cannot read " << file << std::endl;
    return 1;
  }
  const std::string refusal = azure2::io::NotFormat3Message(bytes, file);
  if (!refusal.empty()) {
    std::cout << "ERROR: " << refusal << std::endl;
    return 3;
  }
  azure2::experiment::Format3ReadOptions read;
  read.fileName = file;  // data files as the engine opens them: relative to the working directory
  const azure2::experiment::Format3Project c = azure2::experiment::ReadFormat3Text(bytes, read);
  for (const azure2::experiment::Diagnostic &d : c.diagnostics)
    (d.severity == azure2::experiment::Severity::Error ? std::cout : std::cerr) << d.ToText() << std::endl;
  if (!c.Ok() || !c.levelsRead) {
    std::cout << file << ": not described (" << (c.levelsRead ? "errors above" : "no readable <levels>") << ")." << std::endl;
    return 2;
  }
  if (!record.empty()) {
    bool described = false;  // the kinds Describe.h describes
    for (const azure2::experiment::ValidatedRecord &r : c.validated.records)
      described = described || (r.name == record && (r.kind == "dataset" || r.kind == "extrapolation" || r.kind == "datum" ||
                                                      r.kind == "potential"));
    if (!described) {
      std::cout << "ERROR: --describe --record " << record << ": " << file << " has no data set, extrapolation, datum or "
                << "potential of that name." << std::endl;
      return 2;
    }
  }
  std::cout << (json ? azure2::experiment::DescribeJson(c, record) : azure2::experiment::DescribeText(c, record)) << std::flush;
  if (json) std::cout << std::endl;
  // the stopping of each target from the library: its source (text form only), or the refusal
  const StoppingReport stopping = ReportTargetStoppings(TargetEffectSpecs(c));
  for (const azure2::experiment::Diagnostic &d : stopping.refusals) std::cout << d.ToText() << std::endl;
  if (!json && record.empty()) std::cout << stopping.text << std::flush;
  return stopping.refusals.empty() ? 0 : 2;
}

/// The command was given; it exits with `value`.
bool Exit(int &exitCode, int value) {
  exitCode = value;
  return true;
}

}  // namespace

bool RunProjectCommand(int argc, char **argv, int &exitCode) {
  // --convert OLD NEW [--to 3] writes OLD (format 1 or 2) in format 3 to NEW and exits (docs/upgrade/project_format3.md
  // section 8): with --convert-params, the only commands that read an old file
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--convert") == 0) {
      if (i + 2 >= argc) {
        std::cout << "ERROR: --convert needs two files: --convert old.azr new.azr [options]" << std::endl;
        return Exit(exitCode, 1);
      }
      std::vector<std::string> rest;
      for (int k = i + 3; k < argc; k++) {
        if (strcmp(argv[k], "--to") == 0 && k + 1 < argc) {
          if (strcmp(argv[k + 1], "3") != 0) {
            std::cout << "ERROR: --convert --to " << argv[k + 1] << ": the converter writes format 3 only" << std::endl;
            return Exit(exitCode, 1);
          }
          k++;
          continue;
        }
        rest.push_back(argv[k]);
      }
      return Exit(exitCode, convertProjectTo3(argv[i + 1], argv[i + 2], rest));
    }

  // --convert-params OLD.sav NEW.sav --project OLD.azr renames a parameter file of an old project and exits
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--convert-params") == 0) {
      std::string project;
      for (int k = 1; k + 1 < argc; k++)
        if (strcmp(argv[k], "--project") == 0) project = argv[k + 1];
      if (i + 2 >= argc || project.empty()) {
        std::cout << "ERROR: --convert-params needs two files and the old project: --convert-params old.sav new.sav "
                     "--project old.azr" << std::endl;
        return Exit(exitCode, 1);
      }
      return Exit(exitCode, convertParameters(argv[i + 1], argv[i + 2], project));
    }

  // --validate FILE [--json] checks a format-3 project and exits (docs/upgrade/project_format3.md section 6.2)
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--validate") == 0) {
      if (i + 1 >= argc) {
        std::cout << "ERROR: --validate needs a file: --validate project.azr [--json]" << std::endl;
        return Exit(exitCode, 1);
      }
      bool json = false;
      for (int k = 1; k < argc; k++) json = json || strcmp(argv[k], "--json") == 0;
      return Exit(exitCode, validateProject(argv[i + 1], json));
    }

  // --print-predictions FILE prints what every data set and extrapolation of FILE predicts, as a tree, and exits
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--print-predictions") == 0) {
      if (i + 1 >= argc) {
        std::cout << "ERROR: --print-predictions needs a file: --print-predictions project.azr" << std::endl;
        return Exit(exitCode, 1);
      }
      return Exit(exitCode, printPredictions(argv[i + 1]));
    }

  // --describe FILE [--json] [--record NAME] prints what the engine makes of a format-3 project, or of its record
  // NAME, and exits (section 6.3)
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--describe") == 0) {
      if (i + 1 >= argc) {
        std::cout << "ERROR: --describe needs a file: --describe project.azr [--json] [--record NAME]" << std::endl;
        return Exit(exitCode, 1);
      }
      bool json = false;
      std::string record;
      for (int k = 1; k < argc; k++) {
        json = json || strcmp(argv[k], "--json") == 0;
        if (strcmp(argv[k], "--record") == 0) {
          if (k + 1 >= argc) {
            std::cout << "ERROR: --record needs the name of a record: --describe project.azr --record NAME" << std::endl;
            return Exit(exitCode, 1);
          }
          record = argv[++k];
        }
      }
      return Exit(exitCode, describeProject(argv[i + 1], json, record));
    }

  // --schema prints the table of project format 3 as JSON and exits (docs/upgrade/project_format3.md section 6.1);
  // tests/schema/schema_v3.json is its committed copy
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "--schema") == 0) {
      std::cout << azure2::experiment::Schema::ToJson() << std::endl;
      return Exit(exitCode, 0);
    }
  return false;
}
