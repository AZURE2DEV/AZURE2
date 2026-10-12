#ifndef PATHTOKEN_H
#define PATHTOKEN_H

#include <istream>
#include <string>

/*!
 * Reads a file name that may be wrapped in double quotes.
 *
 * The name sits in the middle of a segment line: optional fields describing
 * composite and UPOS segments follow it, so the rest of the line cannot simply
 * be taken as the name, and a bare ">>" stops at the first space -- which
 * silently truncates any path containing one to its first word.  A quoted name
 * is therefore read as a single unit.  Unquoted names are returned exactly as
 * before, so every existing .azr file keeps parsing identically; only files
 * written with a space in the path use the quoted form.
 */
inline bool ReadPathToken(std::istream &stream, std::string &path) {
  std::string token;
  if (!(stream >> token)) return false;
  if (token.empty() || token[0] != '"') {
    path = token;
    return true;
  }
  path = token.substr(1);
  if (!path.empty() && path[path.size() - 1] == '"') {
    path.erase(path.size() - 1);  // the whole name fitted in one token
    return true;
  }
  // ">>" stopped at the space inside the name; take the remainder up to the
  // closing quote.  It still carries that space, so it is appended as read.
  std::string rest;
  if (std::getline(stream, rest, '"')) path += rest;
  return true;
}

/*!
 * Removes a matching pair of surrounding double quotes, if present.
 */
inline std::string StripPathQuotes(const std::string &path) {
  if (path.size() >= 2 && path[0] == '"' && path[path.size() - 1] == '"')
    return path.substr(1, path.size() - 2);
  return path;
}

#endif
