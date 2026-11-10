#ifndef AZURE2_CONVERT_CONVERT_H
#define AZURE2_CONVERT_CONVERT_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Diagnostics.h"

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

/// Converts the text of an old project.  Never throws for bad input; never reads anything but the data files and the
/// prior-covariance file under options.projectDir.
ConvertResult ConvertToFormat3(const std::string &oldText, const ConvertOptions &options);

/*!
 * The format-3 sections of a project text read with io::ReadRecords and validated with experiment's table:
 * every diagnostic (the reader's and the validator's), and the number of checks the validator leaves to the
 * format-3 reader of step 3b because they need <levels> or the data files (the quantities' refusals marked
 * "build", one per record and refusal).  Used by `AZURE2 --validate`, pyazr.validate and the converter on its own
 * output.
 */
struct Format3Check {
  std::vector<experiment::Diagnostic> diagnostics;
  int buildChecksSkipped = 0;
  bool Ok() const;
};
Format3Check CheckFormat3Text(const std::string &text, const std::string &fileName);

/// The answers of --answers: a flat JSON object {"<record>.<field>": <string or number>, ...}.  Returns "" or what is
/// wrong with the text.
std::string ReadAnswersJson(const std::string &text, std::map<std::string, std::string> &answers);

/// The refusal of section 8.4 for a project that is not format 3 ("" for a format-3 text): it names the first old
/// section and the conversion command.
std::string NotFormat3Message(const std::string &text, const std::string &fileName);

}  // namespace convert
}  // namespace azure2

#endif
