#ifndef RUNMENU_H
#define RUNMENU_H

// The text menu: the reader of a run request that asks (include/RunRequest.h).  It prints the questions of the
// table's rows and reads the answers from the terminal (readline, where built and wanted) or from a piped script.

#include <optional>

#include "RunRequest.h"

class Config;

class RunMenu {
 public:
  /// A menu on `configure`'s output stream; `useReadline`: file names are read through readline (where built).
  RunMenu(Config &configure, bool useReadline) : configure_(configure), useReadline_(useReadline) {}
  /// The command chosen in the menu; none for its Exit line.  Asks until the choice is valid; stops the run (exit 1)
  /// when the input ends.
  std::optional<RunCommand> AskCommand();
  /// Asks `question` and keeps its answer in `request`; asks again an answer out of range (a y/n question from a pipe
  /// stops the run instead, exit 1).
  void Ask(const RunQuestion &question, RunRequest &request);
  /// Says that the answers broke `rule` (the menu then asks again from the rule's first question).
  void Refuse(const RunRule &rule);

 private:
  Config &configure_;
  bool useReadline_;
};

#endif
