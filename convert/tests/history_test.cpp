// The forms of the old project file through its history (tests/convert/fixtures/history/README.md): each era's form as
// its setup program wrote it, the variants of a file an engine of its era read, and project files of the repository as
// they were at their revisions.  A form an engine of its era read with one meaning converts with that meaning, checked
// against the same project written in the form of a later era (records equal); a form with no meaning in format 3, or
// two meanings, is refused with its reason.
//   convert_history_test <repository root>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/convert/Convert.h"

using namespace azure2::convert;

namespace {

int failures = 0;
std::string root;

void check(const std::string &what, bool ok) {
  std::printf("%s  %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

std::string Slurp(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

std::string Dir(const std::string &form) { return root + "/tests/convert/fixtures/history/" + form; }

std::string Text(const std::string &form) {
  const std::string name = form.substr(form.find_last_of('/') + 1);
  return Slurp(Dir(form) + "/" + name + ".azr");
}

ConvertResult Convert(const std::string &text, const std::string &form) {
  ConvertOptions o;
  o.projectDir = Dir(form);
  return ConvertToFormat3(text, o);
}

bool HasNote(const ConvertResult &r, const std::string &part) {
  for (const ConvertNote &n : r.notes)
    if (n.text.find(part) != std::string::npos) return true;
  return false;
}

/// The section `name` of a text without its whole-line comments and with CR removed (the records as written).
std::string Section(const std::string &text, const std::string &name) {
  std::istringstream in(text);
  std::string line, out;
  bool inside = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line == "<" + name + ">") inside = true;
    else if (line == "</" + name + ">") inside = false;
    else if (inside && line.find_first_not_of(' ') != std::string::npos && line[line.find_first_not_of(' ')] != '#')
      out += line + "\n";
  }
  return out;
}

/// The <levels> rows as lists of words (the layout of blanks is the writer's).
std::vector<std::vector<std::string>> Rows(const std::string &text) {
  std::vector<std::vector<std::string>> out;
  std::istringstream in(Section(text, "levels"));
  for (std::string line; std::getline(in, line);) {
    std::istringstream w(line);
    std::vector<std::string> row;
    for (std::string t; w >> t;) row.push_back(t);
    if (!row.empty()) out.push_back(row);
  }
  return out;
}

bool SameRecords(const std::string &a, const std::string &b) {
  for (const char *s : {"data", "extrapolations", "targets"})
    if (Section(a, s) != Section(b, s)) return false;
  return true;
}

std::string Field(const std::string &text, const std::string &header, const std::string &key) {
  std::istringstream in(text);
  std::string line;
  bool inside = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!inside) {
      inside = line == header;
      continue;
    }
    if (line.empty() || line[0] == '<' || line[0] != ' ') break;
    std::istringstream w(line);
    std::string k;
    w >> k;
    if (k != key) continue;
    std::string rest;
    std::getline(w, rest);
    return rest.substr(rest.find_first_not_of(' '));
  }
  return "";
}

std::string Replace(std::string text, const std::string &from, const std::string &to) {
  for (std::size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size()))
    text.replace(at, from.size(), to);
  return text;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: convert_history_test <repository root>\n");
    return 2;
  }
  root = argv[1];
  const std::string g2020 = Text("gui_2020");
  const ConvertResult base = Convert(g2020, "gui_2020");
  check("2020 (6ecaee22, origin/qt5): the setup program's file converts", base.ok);

  // ---- before 2020: 30 columns, <externalCapture>
  {
    const ConvertResult r = Convert(Text("pre2020_external_capture"), "pre2020_external_capture");
    check("before 2020, <externalCapture>: converts with a note", r.ok && HasNote(r, "as the 2020 setup program imported it"));
    const auto a = Rows(r.text), b = Rows(base.text);
    check("before 2020, <externalCapture>: the <levels> rows those of the 2020 file (mask 5 on pair 2)", a == b);
    check("before 2020, <externalCapture>: the records those of the 2020 file", SameRecords(r.text, base.text));
    const ConvertResult n = Convert(Text("pre2020_30_columns"), "pre2020_30_columns");
    bool zero = n.ok;
    const auto c = Rows(n.text);
    for (std::size_t i = 0; zero && i < c.size(); i++) {
      zero = c[i].size() == 31 && c[i][30] == "0";
      std::vector<std::string> x = c[i], y = b[i];
      x.pop_back();
      y.pop_back();
      zero = zero && x == y;
    }
    check("before 2020, 30 columns alone: ecMultMask 0 in every row, every other column the 2020 file's", zero && HasNote(n, "30 columns"));
    const std::string inactive = Replace(Text("pre2020_external_capture"), "<externalCapture>\n1 2", "<externalCapture>\n0 2");
    const ConvertResult i = Convert(inactive, "pre2020_external_capture");
    check("before 2020, an inactive <externalCapture> line: refused (the import applied it, the engine is unknown)",
          !i.ok && i.refusal.find("inactive") != std::string::npos);
    const std::string twice = Replace(Text("pre2020_external_capture"), "1 2 0.5 0.5 5\n", "1 2 0.5 0.5 5\n1 2 0.5 0.5 3\n");
    const ConvertResult t = Convert(twice, "pre2020_external_capture");
    check("before 2020, two masks for one pair: refused", !t.ok && t.refusal.find("second mask") != std::string::npos);
  }

  // ---- the 2020 file's variants an engine or the setup program of the era read
  {
    std::string crlf;
    for (char c : g2020) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
    const ConvertResult r = Convert(crlf, "gui_2020");
    check("2020, CRLF line ends: the same records", r.ok && SameRecords(r.text, base.text));
    const ConvertResult o = Convert(Replace(g2020, "</lastRun>\n", ""), "gui_2020");
    check("2020, <lastRun> not closed at the end: converts, the same records", o.ok && SameRecords(o.text, base.text) && HasNote(o, "<lastRun> is not closed"));
    const ConvertResult oc = Convert(Replace(crlf, "</lastRun>\r\n", ""), "gui_2020");
    check("2020, CRLF and <lastRun> not closed: converts", oc.ok && SameRecords(oc.text, base.text));
    const ConvertResult s = Convert(Replace(g2020, "data/artemov.dat", "data/artemov copy.dat"), "gui_2020");
    check("2020, a data file name with a blank, unquoted, the rest of the line (no such file here): read as that name",
          s.ok && Field(s.text, "dataset artemov_copy", "file") == "\"data/artemov copy.dat\"");
    std::string tabs;
    {
      std::istringstream in(g2020);
      for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] != '<') {
          std::string t;
          bool blank = false;
          for (char c : line) {
            if (c == ' ') {
              if (!blank) t += '\t';
              blank = true;
            } else {
              t += c;
              blank = false;
            }
          }
          line = t;
        }
        tabs += line + "\n";
      }
    }
    const ConvertResult tb = Convert(tabs, "gui_2020");
    check("2020, tabs between the fields: the same records", tabs != g2020 && tb.ok && SameRecords(tb.text, base.text));
    const ConvertResult bl = Convert(Replace(g2020, "<segmentsData>\n", "<segmentsData>\n\n   \n"), "gui_2020");
    check("2020, blank lines and blanks inside a section: the same records", bl.ok && SameRecords(bl.text, base.text));
    std::string noLastRun = g2020;
    noLastRun.erase(noLastRun.find("<lastRun>"), noLastRun.find("</lastRun>\n") + 11 - noLastRun.find("<lastRun>"));
    const ConvertResult nl = Convert(noLastRun, "gui_2020");
    check("2020, no <lastRun>: the same records", nl.ok && SameRecords(nl.text, base.text));
    const std::string winOut = Replace(g2020, std::string("output/") + std::string(93, ' '), std::string("C:\\Users\\J S\\output\\") + std::string(81, ' '));
    const ConvertResult wo = Convert(winOut, "gui_2020");
    check("2020, a Windows output directory with a blank: converts, the directory kept",
          winOut != g2020 && wo.ok && wo.text.find("C:\\Users\\J S\\output\\") != std::string::npos);
    const ConvertResult bs = Convert(Replace(g2020, "data/artemov.dat", "data\\artemov.dat"), "gui_2020");
    check("2020, a data file written with a backslash: converts, the name kept", bs.ok && Field(bs.text, "dataset artemov", "file") == "data\\artemov.dat");
    const ConvertResult w = Convert(Replace(g2020, "data/artemov.dat", "C:/Users/J S/data/artemov.dat"), "gui_2020");
    check("2020, an absolute Windows path with a blank: read as that name",
          w.ok && Field(w.text, "dataset artemov", "file") == "\"C:/Users/J S/data/artemov.dat\"");
  }

  // ---- <targetInt> of 2025-01-26 .. 01-30: no equation token
  {
    const ConvertResult r = Convert(Text("target_2025_01_no_equation"), "target_2025_01_no_equation");
    check("target line of 2025-01-26 (17281967), no coefficients: the targets of the 2020 file",
          r.ok && Section(r.text, "targets") == Section(base.text, "targets"));
    const std::string t1 = Text("target_2025_01_coefficients");
    const ConvertResult c = Convert(t1, "target_2025_01_coefficients");
    const std::string t2 = Replace(t1, "              1              2 0.001 0.0005 ",
                                   "              1               \"a0+a1*x\" 2 0.001 0.0005 ");
    const ConvertResult e = Convert(t2, "target_2025_01_coefficients");
    check("target line of 2025-01-26 with c0 + c1 E: converts with a note", c.ok && HasNote(c, "no equation token"));
    check("target line of 2025-01-26 with c0 + c1 E: the records of the same line written with its equation (1c01a007)",
          t2 != t1 && e.ok && SameRecords(c.text, e.text));
    check("target line of 2025-01-26 with c0 + c1 E: sigma written as the polynomial",
          c.text.find("formula \"a0+a1*x\" a0=0.001 a1=0.0005") != std::string::npos);
  }

  // ---- data lines: energy shift (2025-08-01), components (2025-08-30, -12-22, 2026-05-21), secondary decay (2026-03-09)
  {
    const ConvertResult r = Convert(Text("energy_shift_2025_08"), "energy_shift_2025_08");
    check("energy shift (dcfdb344): norm 1 varied with 5 %, shift 0 varied with 0.005 MeV",
          r.ok && Field(r.text, "dataset meyer_84", "norm") == "1 fit prior 1 5%" &&
              Field(r.text, "dataset meyer_84", "shift") == "offset 0 MeV fit prior 0 0.005");
    const char *forms[] = {"components_2025_08", "components_angle_2025_12", "secondary_decay_2026_03", "components_scaling_2026_05"};
    for (const char *f : forms) {
      const ConvertResult c = Convert(Text(f), f);
      check(std::string(f) + ": converts", c.ok);
      check(std::string(f) + ": SUM line = the line's own term + its components 1->1, 1->2",
            Field(c.text, "dataset meyer_84_2", "reaction") == "12C(p,p)12C + 12C(p,p)12C + 12C(p,g)13N");
      check(std::string(f) + ": RATIO line = the line's term over its component 1->2",
            Field(c.text, "dataset meyer_84_3", "reaction") == "12C(p,p)12C / 12C(p,g)13N");
    }
    const ConvertResult a = Convert(Text("components_angle_2025_12"), "components_angle_2025_12");
    check("components with angles (de248934): 30 and 60 deg kept",
          Field(a.text, "dataset meyer_84_4", "reaction") == "12C(p,p)12C + 12C(p,p)12C @ 30 deg lab + 12C(p,g)13N @ 60 deg lab");
    const ConvertResult u = Convert(Text("secondary_decay_2026_03"), "secondary_decay_2026_03");
    check("secondary decay (2c406c7a): L 2, Ic 3/2, delta 0.1",
          u.text.find("secondary_gamma L 2 Ic 3/2 delta 0.1") != std::string::npos);
    const ConvertResult sc = Convert(Text("components_scaling_2026_05"), "components_scaling_2026_05");
    check("components with scalings (9b5084d8): 0.5 and 2, an angle on one",
          Field(sc.text, "dataset meyer_84_4", "reaction") == "12C(p,p)12C + 0.5*12C(p,p)12C + 2*12C(p,g)13N @ 60 deg lab");
    check("components with scalings: a line of the same file without scalings",
          Field(sc.text, "dataset meyer_84_5", "reaction") == "12C(p,p)12C + 12C(p,p)12C");
  }

  // ---- numbers run together by the writer of 2026-09-25 .. 10-07 (two meanings until split): refused, the split proposed
  {
    const ConvertResult d = Convert(Text("glued_data_2026_09"), "glued_data_2026_09");
    check("glued norm and vary flag: refused with the split as the first reading",
          !d.ok && !d.lineProblem.readings.empty() &&
              d.lineProblem.readings[0].line.find("3.659747237e-06 1 5 0 0.005 1") != std::string::npos);
    const ConvertResult l = Convert(Text("glued_levels_2026_09"), "glued_levels_2026_09");
    check("glued J^pi and E_x in a <levels> row: refused, naming the token",
          !l.ok && l.refusal.find("114.918600000000001") != std::string::npos);
  }

  // ---- project files of the repository as they were at their revisions
  {
    const ConvertResult a = Convert(Text("revisions/13N_2026_04_22"), "revisions/13N_2026_04_22");
    check("test/13N.azr at 81243f81 (2026-04-22, the first project file of the repository): converts", a.ok);
    const ConvertResult b = Convert(Text("revisions/12C_pp_ay_2026_08_04"), "revisions/12C_pp_ay_2026_08_04");
    check("12C_pp_ay at df9dc4d5: a target line of inactive lines takes their nucleus (no question)",
          b.ok && HasNote(b, "applies to no active line"));
    const std::string spin = Text("revisions/spin_observables_reaction_2026_10_18");
    const ConvertResult c = Convert(spin, "revisions/spin_observables_reaction_2026_10_18");
    std::string as10;
    {
      std::istringstream in(spin);
      for (std::string line; std::getline(in, line);) {
        std::istringstream w(line);
        std::vector<std::string> t;
        for (std::string x; w >> x;) t.push_back(x);
        if (t.size() > 8 && t[7] == "9" && t[8].find(':') != std::string::npos) {
          t[7] = "10";
          line.clear();
          for (const std::string &x : t) line += (line.empty() ? "" : " ") + x;
        }
        as10 += line + "\n";
      }
    }
    const ConvertResult c10 = Convert(as10, "revisions/spin_observables_reaction_2026_10_18");
    check("spin observables at 8f80abde (type 9 with an operator): converts as type 10 with a note",
          c.ok && HasNote(c, "read as type 10"));
    check("spin observables at 8f80abde: the records of the same lines written as type 10", c10.ok && SameRecords(c.text, c10.text));
    const ConvertResult m = Convert(Text("revisions/mix_data_2026_10_12"), "revisions/mix_data_2026_10_12");
    check("mixture at a38494e8: refused, naming the revert", !m.ok && m.refusal.find("a60577f6") != std::string::npos);
    const ConvertResult p = Convert(Text("revisions/polarization_product_2026_09_17"), "revisions/polarization_product_2026_09_17");
    check("isDiff 8 at 922c0537: refused, naming its removal", !p.ok && p.refusal.find("330dc1eb") != std::string::npos);
  }

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
