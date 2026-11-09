// The lossless .azr document (Document.h; docs/upgrade/project.md section 2), on synthetic documents:
//  1. Write() = input for LF, CRLF, CR, mixed endings, no final newline, empty input, comments and text outside sections;
//  2. tags: whole trimmed lines only (not in a comment, not with other text on the line), nested sections, problems
//     (unclosed, closing without opening, crossing);
//  3. edits: ReplaceLines keeps every byte outside its range and gives new lines the replaced line's ending;
//     RemoveSection removes exactly the section's lines.
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/io/Document.h"

using azure2::io::Document;

namespace {

int failures = 0;
int checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

}  // namespace

int main() {
  const std::vector<std::string> inputs = {
      "",
      "<config>\n1\n</config>\n",
      "<config>\r\n1\r\n</config>\r\n",
      "<config>\r1\r</config>\r",
      "<config>\n1\r\n</config>\r",
      "<config>\n1\n</config>",
      "# a comment <levels>\n  <config>  \n\t1 2 3   \n</config>\t\n\ntrailing text\n\n\n",
      "<mcmc>\nnwalkers 50\n<parameters>\nenergy_1 1 2\n</parameters>\n</mcmc>\n",
  };
  for (size_t i = 0; i < inputs.size(); i++)
    check("round trip of input " + std::to_string(i), Document::Parse(inputs[i]).Write() == inputs[i]);

  {
    Document d = Document::Parse("# <levels>\nx <levels>\n  <config>  \n1\n</config>\n");
    check("a tag in a comment or with other text is not a tag", d.Find("levels") == nullptr);
    const Document::Section *c = d.Find("config");
    check("a trimmed tag is a tag", c && c->open == 2 && c->close == 4);
    check("no problems", d.Problems().empty());
    check("body", c && d.Body(*c) == std::vector<std::string>{"1"});
  }
  {
    Document d = Document::Parse("<mcmc>\nnwalkers 50\n<parameters>\na\n</parameters>\n</mcmc>\n");
    const Document::Section *m = d.Find("mcmc");
    check("nested: outer found", m != nullptr && m->close == 5);
    const Document::Section *p = m ? d.FindChild(*m, "parameters") : nullptr;
    check("nested: inner found as child", p && p->open == 2 && p->close == 4 && p->depth == 1);
    check("nested: inner not top level", d.Find("parameters") == nullptr);
  }
  {
    check("unclosed reported", Document::Parse("<levels>\n1\n").Problems().size() == 1);
    check("closing without opening reported", Document::Parse("1\n</levels>\n").Problems().size() == 1);
    Document d = Document::Parse("<a>\n<b>\n</a>\n</b>\n");
    check("crossing reported", !d.Problems().empty());
  }
  {
    const std::string in = "<config>\r\n1\r\n2\r\n</config>\r\n<levels>\r\nL\r\n</levels>\r\n";
    Document d = Document::Parse(in);
    d.ReplaceLines(2, 3, {"two", "three"});
    check("replace keeps the rest, new lines take the ending",
          d.Write() == "<config>\r\n1\r\ntwo\r\nthree\r\n</config>\r\n<levels>\r\nL\r\n</levels>\r\n");
    const Document::Section *l = d.Find("levels");
    check("sections re-found after an edit", l && l->open == 5);
    d.RemoveSection(*d.Find("config"));
    check("remove a section", d.Write() == "<levels>\r\nL\r\n</levels>\r\n");
  }
  {
    Document d = Document::Parse("<a>\nx\n</a>");
    d.ReplaceLines(1, 3, {"y", "</a>"});
    check("a last line without an ending keeps none", d.Write() == "<a>\ny\n</a>");
  }
  std::printf("%s: %d checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
