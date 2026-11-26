#ifndef AZURE2_CONVERT_CONVERT_H
#define AZURE2_CONVERT_CONVERT_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "azure2/convert/LineRepair.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/io/Document.h"

namespace azure2 {
namespace convert {

/*!
 * The converter from the old project formats (1 and 2) to format 3 (docs/upgrade/project_format3.md section 8).
 * All knowledge of the old formats lives in this directory: the line parsers the engine still uses until the
 * format-3 reader replaces them (SegLine.h, ExtrapLine.h, ComponentList.h, TargetIntLine.h, AbundanceLine.h,
 * ScatteringLengthLine.h, PriorCentreRow.h, ThmLines.h) and this conversion, which is used only by
 * `AZURE2 --convert ... --to 3` and `pyazr.convert_azr(..., to=3)`.
 *
 * The conversion keeps the meaning the old engine gave a file (section 8.2): where that meaning was implicit, the
 * new text states it with a `# CONVERTED:` comment and a note in the report; where the old meaning is undefined it
 * asks (section 8.3), and a question without an answer stops the conversion unless placeholders are allowed.  The
 * output is validated with the format-3 table before it is returned; a text with errors is never returned as a
 * result to write.
 */
struct ConvertOptions {
  /// The directory of the old project: data files and the prior-covariance file are read relative to it.
  std::string projectDir;
  /// Write every extrapolation grid as the list of the energies (and angles) the old engine accumulated (decision
  /// F10), instead of `range lo hi step s` with E_k = lo + k s.
  bool exactGrids = false;
  /// Ask, instead of writing with a note, wherever the old meaning was implicit (decision F12).
  bool strict = false;
  /// Write '?' for an unanswered question (the file then does not validate until it is answered, F18).
  bool allowQuestions = false;
  /// Answers to questions: "<record>.<field>" -> the field's value text.
  std::map<std::string, std::string> answers;
};

/// A change the converter made without asking, or a value it dropped because the old engine did not use it.
struct ConvertNote {
  std::string construct;  ///< the construct of section 8.1 ("normgroup", "<targetInt> line", ...)
  std::string record;     ///< the format-3 record it concerns ("" for a section)
  std::string text;
};

/// What the converter cannot decide: answered with --answers {"<key>": "<value text>"}.
struct ConvertQuestion {
  std::string construct;
  std::string key;  ///< "<record>.<field>"
  std::string text;
  std::string proposal;  ///< what the converter would write (empty when it cannot propose)
};

/// Per construct of section 8.1: how many old lines (or blocks) it converted, with how many notes and questions.
struct ConstructCount {
  int converted = 0;
  int notes = 0;
  int questions = 0;
};

struct ConvertResult {
  bool ok = false;      ///< `text` is a valid format-3 project (or has only '?' placeholders, with allowQuestions)
  std::string refusal;  ///< why there is no result: "already format 3", a line the old engine refused, ...
  /// When the refusal is about one <segmentsData>, <segmentsTest> or <targetInt> line: where it is, what is expected,
  /// the readings proposed (LineRepair.h); `refusal` then ends with lineProblem.Text().
  LineProblem lineProblem;
  std::string text;     ///< the converted project
  std::vector<std::string> format2Changes;  ///< io::ConvertToFormat2's list
  std::vector<ConvertNote> notes;
  std::vector<ConvertQuestion> questions;
  std::vector<experiment::Diagnostic> diagnostics;  ///< the format-3 validation of `text` (errors only when !ok)
  /// Old parameter name -> new name (segment_3_norm -> s03.norm).
  std::vector<std::pair<std::string, std::string>> parameterNames;
  std::map<std::string, ConstructCount> constructs;
  /// Side files the caller writes beside the output (path relative to the project, contents): `<f>.f3` copies of a
  /// prior-covariance file with renamed parameters.  The converter never writes data files.
  std::vector<std::pair<std::string, std::string>> sideFiles;

  /// The JSON report (--json-report): every list above.
  std::string ReportJson(const std::string &input, const std::string &output) const;
};

/*!
 * Converts a format-1 project to format 2 in place, the converter's first step (docs/upgrade/project.md section 4), and
 * returns what it changed, one line per change (empty for a file already in format 2 or 3):
 *  - <lastRun> removed (the Qt GUI's runtime options, read by nothing; ConvertToFormat3 notes a ticked "use GSL
 *    Coulomb functions" box first);
 *  - <mcmc>'s run settings removed (the CLI asks for them), and <mcmc> with them when its <parameters> are empty;
 *  - useAdaptiveGrid moved from <potential> to <config>, and <potential> removed unless useHybridPotential=1 for
 *    the default or for any pair;
 *  - "format=2" added as the last line of <config>.
 * Every other byte is kept.  Converting twice is converting once.
 */
std::vector<std::string> ConvertToFormat2(io::Document &doc);

/// Converts the text of an old project.  Never throws for bad input; never reads anything but the data files and the
/// prior-covariance file under options.projectDir.
ConvertResult ConvertToFormat3(const std::string &oldText, const ConvertOptions &options);

/// The answers of --answers: a flat JSON object {"<record>.<field>": <string or number>, ...}.  Returns "" or what is
/// wrong with the text.
std::string ReadAnswersJson(const std::string &text, std::map<std::string, std::string> &answers);

/*!
 * A parameter file written for an old project (param.sav, param.par, param.fit: "name value error [fixed]" lines) with
 * the names the project's conversion gives its parameters (ConvertResult::parameterNames, section 8.4: AZURE2
 * --convert-params, pyazr.convert_params): the first token of a line that the map renames (segment_3_norm -> s03.norm)
 * is replaced, the name column kept right-aligned where the new name fits; every other byte is kept.  `renamed` gets
 * the pairs applied, in file order.
 */
std::string RenameParameterFile(const std::string &text, const std::vector<std::pair<std::string, std::string>> &names,
                                std::vector<std::pair<std::string, std::string>> &renamed);

/// A parameter file of an old project, renamed (RenameParameterFile) by the map of that project's conversion.
struct ParameterFileConversion {
  std::string refusal;  ///< "" or why there is no result (the project is format 3, or the converter refuses it)
  std::string text;     ///< the renamed file
  std::vector<std::pair<std::string, std::string>> renamed;  ///< old name -> new name, in file order
};

/*!
 * `parameterText` (a param.sav, param.par or param.fit) renamed for the format-3 conversion of `oldProjectText`, the old
 * project the file was written for (the converter keeps it unchanged): the map is ConvertResult::parameterNames of its
 * conversion (the "parameters" of the JSON report), which a format-3 project alone does not determine (its records
 * may have been renamed or reordered since).  Questions of the conversion do not matter for names: they are allowed.
 * `projectDir` as ConvertOptions::projectDir.  AZURE2 --convert-params and pyazr.convert_params.
 */
ParameterFileConversion ConvertParameterFile(const std::string &parameterText, const std::string &oldProjectText,
                                             const std::string &projectDir);

}  // namespace convert
}  // namespace azure2

#endif
