// Headless round trip of the optional <thm> block through the GUI.
//
// The GUI rebuilds a .azr from its tabs on every save.  The <thm> block (THM
// options, read by the engine's Config::ReadThmBlock) has no tab, so it has to
// be carried verbatim from the file that was opened to the file that is
// written -- and a project without the block must be written exactly as
// before.  Checked here through AZURESetup's own open/save, on tests/7Li_p_a.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <QTextStream>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

#ifndef AZURE2_SOURCE_DIR
#error "AZURE2_SOURCE_DIR must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll());
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}

// open() then save() through the main window; returns the saved text.
static QString openAndSave(AZURESetup& w, const QString& in, const QString& out) {
  QFile::remove(out);
  QFile::copy(in, out);
  w.open(out);
  w.saveProject();
  return slurp(out);
}

int main(int argc, char** argv) {
  // Keep the recent-file list the GUI writes out of the user's settings.
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("thm_block_test");

  const QString block =
      "kinematics=kf3body   # data divided by the full three-body KF\n"
      "\n"
      "  spectatorEnergy[5] = 0.5\n"
      "entranceL=incoherent\n";
  const QStringList blockLines = QString(block).split('\n').mid(0, 4);

  // 1. The stream helper, on its own.
  {
    QString text = QString("<config>\n</config>\n  <thm>  # options\n") + block +
                   "</thm>   # end\n<levels>\n";
    QTextStream in(&text);
    QStringList lines;
    bool present = false;
    bool good = AZURESetup::readThmBlock(in, lines, present);
    ok("block found anywhere in the file", good && present);
    ok("lines kept verbatim, comments and blanks included", lines == blockLines,
       lines.join("|"));
    QString written;
    QTextStream out(&written);
    AZURESetup::writeThmBlock(out, lines);
    out.flush();
    ok("written back between <thm> and </thm>", written == "<thm>\n" + block + "</thm>\n", written);

    QString none = "<config>\n</config>\n<levels>\n</levels>\n";
    QTextStream in2(&none);
    good = AZURESetup::readThmBlock(in2, lines, present);
    ok("no block: success, nothing stored", good && !present && lines.isEmpty());

    QString open = "<thm>\nvertex=onshell\n<levels>\n";
    QTextStream in3(&open);
    good = AZURESetup::readThmBlock(in3, lines, present);
    ok("unterminated block is an error, as in the engine", !good && !present);
  }

  // 2. Through AZURESetup: open a project, save it, compare.
  QTemporaryDir work;
  const QString src = QString(AZURE2_SOURCE_DIR) + "/tests/7Li_p_a/7Li_p_a.azr";
  QString plain = slurp(src);
  ok("found tests/7Li_p_a", !plain.isEmpty(), src);
  if(!plain.endsWith('\n')) plain += '\n';
  const QString full = "<thm>\n" + block + "</thm>\n";
  spit(work.filePath("plain.in"), plain);
  spit(work.filePath("end.in"), plain + full);
  QString top = plain;
  top.replace("<levels>", full + "<levels>");
  spit(work.filePath("top.in"), top);

  AZURESetup w;
  const QString savedPlain = openAndSave(w, work.filePath("plain.in"), work.filePath("a.azr"));
  const QString savedEnd = openAndSave(w, work.filePath("end.in"), work.filePath("a.azr"));
  const QString savedTop = openAndSave(w, work.filePath("top.in"), work.filePath("a.azr"));

  ok("GUI wrote the project", savedPlain.contains("</targetInt>"));
  ok("no block in, no block out", !savedPlain.contains("<thm>"));
  ok("block at the end is kept, once", savedEnd.count(full) == 1, savedEnd.right(600));
  ok("block before <levels> is kept, once", savedTop.count(full) == 1);
  ok("its position in the input does not matter", savedTop == savedEnd);
  QString stripped = savedEnd;
  stripped.remove(full);
  ok("nothing else changes", stripped == savedPlain);
  ok("written after </targetInt>", savedEnd.contains("</targetInt>\n" + full));

  // A second open/save of the GUI's own output keeps the block in place.  The
  // rest of that file is not a fixed point (the <parameterSettings> rows are
  // regenerated on open), so it is compared with the plain file's second save
  // rather than with the first.
  spit(work.filePath("again.in"), savedEnd);
  spit(work.filePath("againPlain.in"), savedPlain);
  const QString savedAgain = openAndSave(w, work.filePath("again.in"), work.filePath("a.azr"));
  const QString savedAgainPlain =
      openAndSave(w, work.filePath("againPlain.in"), work.filePath("a.azr"));
  ok("second round trip keeps the block, once", savedAgain.count(full) == 1);
  QString strippedAgain = savedAgain;
  strippedAgain.remove(full);
  ok("second round trip: nothing else differs from the plain file", strippedAgain == savedAgainPlain);
  if(savedAgainPlain != savedPlain)
    std::cout << "  note  a second open/save of the plain project changes it too "
                 "(<parameterSettings>, not the <thm> block)" << std::endl;

  // Opening a project without the block after one with it must not carry it over.
  const QString savedAfter = openAndSave(w, work.filePath("plain.in"), work.filePath("a.azr"));
  ok("block does not leak into the next project", savedAfter == savedPlain);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
