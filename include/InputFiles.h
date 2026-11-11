#ifndef INPUT_FILES_H
#define INPUT_FILES_H

#include <istream>
#include <memory>
#include <string>

/*!
 * Reading the input files of a project (the .azr, data, covariance and target-effect files).
 *
 * A pyazr session rebuilds its model from the files on every channel-radius change, nuclear-potential
 * change and mode switch.  Before 2026-11-04 each rebuild read every file again: a free-radius MCMC
 * (3H+d, 16 forked workers, data on /groups) stalled at 1 % CPU on those reads.  With caching on,
 * Open() keeps each file's contents in memory and reads it again only when the file changes on disk
 * (its size, modification time, or on POSIX its inode, device or status-change time); the values a
 * rebuild sees are therefore always those of the file as it now is.  Caching is per process, shared
 * by the sessions in it, and capped in size.
 */
namespace InputFiles {

/// A stream over the file at path: with cached = true (Config::cacheInputFiles) one over the
/// remembered contents, otherwise a plain std::ifstream.  Test the stream as before: a missing or
/// unreadable file gives a stream in the fail state either way.
std::unique_ptr<std::istream> Open(const std::string &path, bool cached);

struct Stats {
  long reads = 0;   ///< files read from disk through the cache
  long served = 0;  ///< opens answered from memory
  long files = 0;   ///< files held now
};
Stats GetStats();
void Clear();

}  // namespace InputFiles

#endif
