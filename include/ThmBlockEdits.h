#ifndef THM_BLOCK_EDITS_H
#define THM_BLOCK_EDITS_H

#include <string>
#include <vector>

#include "azure2/experiment/SetupSpecs.h"

/*!
 * The writer of the <thm> block beside its reader (ThmExperiment.h): the lines of an edited global option or
 * experiment key, each checked by the reader (ParseThmOptionLine, ParseThmExperimentLine) before it is offered; and
 * the block as the engine reads a format-3 one (Format3ThmLines), with which the whole block is checked
 * (CheckThmBlockLines).  Lines an edit does not touch keep their bytes.  pyazr's Project edits <thm> through these
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
bool ThmExperimentEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &key,
                        const std::string &value, ThmChange change, std::vector<ThmLineEdit> &edits, std::string &why);

/// `lines` with the edits applied (in the order given; they are made so that applying them in order is right).
std::vector<std::string> ApplyThmLineEdits(const std::vector<std::string> &lines, const std::vector<ThmLineEdit> &edits);

/*!
 * The <thm> section of a format-3 project (its lines from the <thm> tag to </thm>) in the grammar the reader takes:
 * format 3 names an experiment's data sets and a data set's weight table in its records (project_format3.md section
 * 3.9; Schema::Build, ThmMembershipSpec), so the members are appended to each experiment[<name>] line as
 * segments=<keys> and the weights as weight[<k>]= / weightTest[<k>]= lines before </thm>.
 */
std::vector<std::string> Format3ThmLines(const std::vector<std::string> &section,
                                         const azure2::experiment::ThmMembershipSpec &thm);

/// "" or why the block does not hold an experiment a data set or extrapolation names (`thm <name>`): the reader
/// would not know it.
std::string CheckFormat3ThmMembers(const std::vector<std::string> &lines, const azure2::experiment::ThmMembershipSpec &thm);

#endif
