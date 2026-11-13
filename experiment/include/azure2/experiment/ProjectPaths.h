#ifndef AZURE2_EXPERIMENT_PROJECTPATHS_H
#define AZURE2_EXPERIMENT_PROJECTPATHS_H

#include <string>

namespace azure2 {
namespace experiment {

/*!
 * The one rule for a path a project names (a data, reference, grid, covariance or table file, the output and checks
 * directories, a file named at a prompt): it is relative to the project's directory, never to the process's working
 * directory, which a library caller (pyazr, the web UI) shares with every other session of the process.
 *
 * The file to open for `path` as written: `path` itself when it is absolute (a leading '/' or '\\') or when no
 * directory is given, `directory` for an empty `path` (an empty output directory is the project's), else
 * `directory`/`path`.  Messages print `path` as written, never this value
 * (docs/upgrade/evaluation.md, section 7.5).  The engine's Config::OpenPath calls it with its baseDirectory.
 */
inline std::string PathToOpen(const std::string &directory, const std::string &path) {
  const bool absolute = !path.empty() && (path[0] == '/' || path[0] == '\\');
  if (absolute || directory.empty()) return path;
  if (path.empty()) return directory;
  return directory.back() == '/' ? directory + path : directory + "/" + path;
}

}  // namespace experiment
}  // namespace azure2

#endif
