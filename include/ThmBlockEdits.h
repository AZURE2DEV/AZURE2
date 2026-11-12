#ifndef THM_BLOCK_EDITS_H
#define THM_BLOCK_EDITS_H

#include <string>
#include <vector>

#include "ThmOptions.h"
#include "azure2/experiment/SetupSpecs.h"

/*!
 * The writer of the <thm> block beside its reader (ThmExperiment.h): the lines of an edited global option or
 * experiment key, each checked by the reader (ParseThmOptionLine, ParseThmExperimentLine) before it is offered; and
 * the section as the engine reads it (ReadFormat3ThmSection).  Lines an edit does not touch keep their bytes.  pyazr's Project edits <thm> through these
 * (api/src/pybind_thm.cpp, docs/upgrade/pyazr_on_records.md "Edits of the nuclear part").
 *
 * `lines` are the lines strictly between <thm> and </thm>; an edit's `line` indexes them (Insert: before that line,
 * lines.size() = at the end).
 */
struct ThmLineEdit {
  enum class Kind { Replace, Remove, Insert };
  Kind kind = Kind::Replace;
  int line = 0;
  std::string text;
};
enum class ThmChange { Set, Remove };

/// One key=value of the block as the reader takes it: a global option (experiment "") or a key of an experiment.
struct ThmEntry {
  std::string experiment;
  std::string key, value;
  int line = 0;  ///< index into the lines
};
/// Every key=value of the lines, in order (a global line that is not key=value is left out: the reader refuses it).
std::vector<ThmEntry> ThmBlockEntries(const std::vector<std::string> &lines);

/// A global option `key=value` (vertex, kinematics, entranceL, coulombIntegral, spectatorEnergy,
/// spectatorEnergy[<pair key>]): Set rewrites the line that wins (the last of the key; earlier ones are removed) or
/// adds one after the last global line; Remove removes every line of the key.  weight[...] is a data set's field in
/// format 3 (thm_weight) and refused here.  False with `why` (the reader's message) when the line would not read.
bool ThmOptionEdits(const std::vector<std::string> &lines, const std::string &key, const std::string &value,
                    ThmChange change, std::vector<ThmLineEdit> &edits, std::string &why);

/// A key of experiment[`name`]: Set replaces the key's `key=value` token where it stands or appends it to the
/// experiment's last line (a new experiment gets a line of its own at the end); Remove drops the token, or with an
/// empty `key` every line of the experiment.  segments= is refused: in format 3 a data set names its experiment.
/// `why` does not repeat "experiment[<name>]": the caller names the experiment it edits.
bool ThmExperimentEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &key,
                        const std::string &value, ThmChange change, std::vector<ThmLineEdit> &edits, std::string &why);

/// `lines` with the edits applied (in the order given; they are made so that applying them in order is right).
std::vector<std::string> ApplyThmLineEdits(const std::vector<std::string> &lines, const std::vector<ThmLineEdit> &edits);

/*!
 * The <thm> section of a format-3 project read into `thm` (Config::ReadThmBlock; pyazr's check of an edit,
 * check_format3_thm): `section` its lines from the <thm> tag (as Config::SectionStream gives them; empty when the
 * project has none), the global options and experiment[<name>] lines parsed (ParseThmOptionLine,
 * ParseThmExperimentLine); the members of each experiment and the weight tables are those of the records
 * (`members`, ThmMembershipSpec: a data set's fields thm and thm_weight), the tables read relative to `projectFile`
 * when `readFiles`.  An experiment that no active data set names is left out and listed in `idle`; `warnings` gets
 * the Coulomb-consistency warnings.  Returns "" or the refusal (printed after "ERROR: ").
 */
std::string ReadFormat3ThmSection(const std::vector<std::string> &section,
                                  const azure2::experiment::ThmMembershipSpec &members, const std::string &projectFile,
                                  bool readFiles, ThmOptions &thm, std::vector<std::string> &idle,
                                  std::vector<std::string> &warnings);

/// The experiments no active data set names: those only inactive data sets name (in the order of the records), then
/// those of the section's experiment[<name>] lines that no data set names.  They are ignored, with a warning
/// (ReadFormat3ThmSection).
std::vector<std::string> ThmExperimentsWithoutMembers(const std::vector<std::string> &section,
                                                      const azure2::experiment::ThmMembershipSpec &thm);

/// "" or why the block does not hold an experiment an active data set names (`thm <name>`).  An internal guard: the
/// project's reading refuses such a name first (E_THM_EXPERIMENT_UNKNOWN, ProjectReader.cpp), so a project file
/// cannot reach it; a caller that builds a ThmMembershipSpec without the reader could.
std::string CheckFormat3ThmMembers(const std::vector<std::string> &lines, const azure2::experiment::ThmMembershipSpec &thm);

#endif
