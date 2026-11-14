#ifndef THM_KEYS_H
#define THM_KEYS_H

#include <string>
#include <vector>

#include "ThmExperiment.h"

/*!
 * The keys of the <thm> block as two tables, the one place each key is read and described: the global options (every
 * line but experiment[...]) and the keys of an experiment[<name>] line.  ParseThmOptionLine and
 * ParseThmExperimentLine (ThmExperiment.h) read a key by its row; the binding serves the rows as the record kinds of
 * a schema (pyazr.kinds()["thm"]), so the web page draws its forms from them.  A new key is one row: its description
 * and its parser.
 */

/// The description of a key, in the shape of a field of a record kind of the schema.
struct ThmKeyRow {
  std::string key;
  std::string type;                 ///< the schema's value type: enum, number, integer or text
  std::vector<std::string> values;  ///< the words of an enum, in the order they are offered
  std::string title;                ///< the title with its unit ("Ebeam (MeV)"); "" for the key itself
  std::string unit;                 ///< the value's unit ("" for none)
  std::string defaultValue;         ///< the value the engine takes when the key is not given ("" for none)
  std::string grammar;              ///< the value's grammar, for a text
  std::string presentWith;          ///< the key that must be given for this one to apply ("" for none)
  std::string doc;
};

/// A global option: its row, whether a per-pair form `key[<pair key>]` exists, the value's grammar as a refusal
/// states it, and its parser (false when the value is not one the engine takes).
struct ThmOptionKey {
  ThmKeyRow row;
  bool perPair = false;
  std::string usage;
  bool (*parse)(const std::string &value, ThmOptionLine &out) = nullptr;
};

/// A key of an experiment line: its row and its parser (sets the experiment's fields; "" or what is wrong).
struct ThmExperimentKey {
  ThmKeyRow row;
  std::string (*parse)(const std::string &value, ThmExperiment &x) = nullptr;
};

/// The two record kinds: their names, titles and docs.
struct ThmKeyKind {
  std::string name, title, doc;
};
const ThmKeyKind &ThmOptionKind();
const ThmKeyKind &ThmExperimentKind();

const std::vector<ThmOptionKey> &ThmOptionKeys();
const std::vector<ThmExperimentKey> &ThmExperimentKeys();

/// "vertex=onshell|constant|perlevel, ..., spectatorEnergy[<pair key>]=<MeV>": the options as a refusal lists them.
std::string ThmOptionUsage();

#endif
