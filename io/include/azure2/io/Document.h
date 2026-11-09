#ifndef AZURE2_IO_DOCUMENT_H
#define AZURE2_IO_DOCUMENT_H

#include <cstddef>
#include <string>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * An .azr file as a lossless document (docs/upgrade/project.md section 2): its lines, each with its own line ending,
 * and the sections <name> ... </name> found in them.  Write() of an unedited document gives back the input byte for
 * byte, whatever its line endings, a missing final newline, comments, blank lines and text outside the sections.
 *
 * A tag is a line that is <name> or </name> after trimming spaces and tabs (and the line ending): a tag inside a
 * comment or on a line with other text is not a tag (P2).  Sections nest (<mcmc> holds <parameters>); a section
 * opened and not closed is reported by Problems() and has no end.
 */
class Document {
 public:
  struct Section {
    std::string name;
    std::size_t open = 0;   ///< line index of <name>
    std::size_t close = 0;  ///< line index of </name>; == npos if the section is not closed
    int depth = 0;          ///< 0 for a top-level section
    std::size_t parent = npos;  ///< index into Sections() of the enclosing section, or npos
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);
  };

  /// Parse bytes (any line endings: LF, CRLF or CR).
  static Document Parse(const std::string &bytes);

  /// The bytes: every line with its own line ending, in order.
  std::string Write() const;

  /// Line i without its line ending.
  const std::string &Line(std::size_t i) const { return lines_[i].text; }
  /// Line i's ending ("\n", "\r\n", "\r", or "" for a last line without one).
  const std::string &Ending(std::size_t i) const { return lines_[i].ending; }
  std::size_t NumLines() const { return lines_.size(); }

  /// Every section, in order of their opening tags.
  const std::vector<Section> &Sections() const { return sections_; }
  /// The first top-level section of that name, or nullptr.
  const Section *Find(const std::string &name) const;
  /// The first section of that name nested directly in `parent`, or nullptr.
  const Section *FindChild(const Section &parent, const std::string &name) const;
  /// The lines strictly between a section's tags (empty for an unclosed section).
  std::vector<std::string> Body(const Section &s) const;

  /// Structural problems: unclosed sections, a closing tag without its opening tag, tags crossing (<a><b></a></b>),
  /// a top-level section given twice (only the first is read).
  const std::vector<std::string> &Problems() const { return problems_; }

  /*!
   * Replace lines [first, last) by `replacement` (lines without endings; each takes the ending of the first replaced
   * line, or the document's most common ending when none is replaced), then re-find the sections.  The bytes outside
   * the range are untouched.
   */
  void ReplaceLines(std::size_t first, std::size_t last, const std::vector<std::string> &replacement);
  /// Remove a section, its tags included.
  void RemoveSection(const Section &s);

  /// The document's most common line ending ("\n" for a document without one).
  std::string CommonEnding() const;

 private:
  struct LineRec {
    std::string text;
    std::string ending;
  };
  void FindSections();
  std::vector<LineRec> lines_;
  std::vector<Section> sections_;
  std::vector<std::string> problems_;
};

/// The tag name of a line that is "<name>" (opening) or "</name>" (closing) after trimming; empty otherwise.
std::string TagOf(const std::string &line, bool &closing);

}  // namespace io
}  // namespace azure2

#endif
