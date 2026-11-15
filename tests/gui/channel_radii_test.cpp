// Per-channel radii in the GUI (<channelRadii>, 2026-11-10): ChannelRadiiModel's line format and precedence (the
// engine's: l+s+J > l+J > l+s > l), the Pairs tab's "Per-channel radii" box, and the project file round trip -- a
// block that the GUI dropped on save before 2026-11-10.
//
// Project: tests/hard_sphere/models/no_dummies.azr (p + 12C, pair radius 5.0 fm) with a <channelRadii> block added.
// With AZURE2_CHANNEL_RADII_PNG=<file> (under xvfb-run) the Pairs tab is saved as an image for a layout check.
#include <QApplication>
#include <QFile>
#include <QPushButton>
#include <QTableView>
#include <QTemporaryDir>
#include <QTextStream>
#include <cmath>
#include <iostream>
#include "AZURESetup.h"
#include "ChannelRadiiModel.h"
#include "Config.h"
#include "PairsTab.h"
struct SegPairs {int firstPair; int secondPair;};
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

namespace {
int failures = 0;
void ok(const std::string &what, bool cond, const QString &detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if (!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) failures++;
}
QString readAll(const QString &path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return QString();
  return QTextStream(&f).readAll();
}
void writeAll(const QString &path, const QString &text) {
  QFile f(path);
  f.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream(&f) << text;
}
bool sameRows(const QList<ChannelRadiusRow> &a, const QList<ChannelRadiusRow> &b) {
  if (a.size() != b.size()) return false;
  for (int i = 0; i < a.size(); i++)
    if (a[i].pairKey != b[i].pairKey || a[i].l != b[i].l || a[i].s != b[i].s || a[i].J != b[i].J ||
        a[i].radius != b[i].radius)
      return false;
  return true;
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;

  std::cout << "1. line format" << std::endl;
  double v = 0;
  ok("spin '3/2' -> 1.5", ChannelRadiiModel::parseSpin("3/2", v) && v == 1.5);
  ok("spin '1' -> 1, '' and 'any' -> -1", ChannelRadiiModel::parseSpin("1", v) && v == 1.0 &&
                                              ChannelRadiiModel::parseSpin("", v) && v == -1.0 &&
                                              ChannelRadiiModel::parseSpin("any", v) && v == -1.0);
  ok("spin '1/3' and '-1/2' refused", !ChannelRadiiModel::parseSpin("1/3", v) && !ChannelRadiiModel::parseSpin("-1/2", v));
  ok("spinText 0.5 -> '1/2', 2 -> '2', -1 -> ''", ChannelRadiiModel::spinText(0.5) == "1/2" &&
                                                     ChannelRadiiModel::spinText(2.0) == "2" &&
                                                     ChannelRadiiModel::spinText(-1.0).isEmpty());
  ChannelRadiusRow r;
  ok("'2 1 4.25 s=3/2 J=1/2  # c' parses", ChannelRadiiModel::parseLine("2 1 4.25 s=3/2 J=1/2  # c", r) && r.pairKey == 2 &&
                                                r.l == 1 && r.radius == 4.25 && r.s == 1.5 && r.J == 0.5);
  ok("and writes back as '2 1 4.25 s=3/2 J=1/2'", ChannelRadiiModel::lineFor(r) == "2 1 4.25 s=3/2 J=1/2");
  ok("malformed lines refused", !ChannelRadiiModel::parseLine("1 1", r) && !ChannelRadiiModel::parseLine("1 1 0", r) &&
                                    !ChannelRadiiModel::parseLine("1 -1 4", r) &&
                                    !ChannelRadiiModel::parseLine("1 1 4 t=1/2", r) &&
                                    !ChannelRadiiModel::parseLine("1 1 4 s=1/2 s=3/2", r));

  std::cout << "2. precedence and editing" << std::endl;
  ChannelRadiiModel m;
  m.addRow(1, 1, 4.0);
  const int sRow = m.addRow(1, 1, 4.1);
  m.setData(m.index(sRow, ChannelRadiiModel::SCol), "1/2");
  const int jRow = m.addRow(1, 1, 4.2);
  m.setData(m.index(jRow, ChannelRadiiModel::JCol), "3/2");
  const int sjRow = m.addRow(1, 1, 4.3);
  m.setData(m.index(sjRow, ChannelRadiiModel::SCol), "1/2");
  m.setData(m.index(sjRow, ChannelRadiiModel::JCol), "3/2");
  ok("l+s+J row wins", m.radiusFor(1, 1, 0.5, 1.5, 5.0) == 4.3);
  m.removeRowAt(sjRow);
  ok("then l+J over l+s", m.radiusFor(1, 1, 0.5, 1.5, 5.0) == 4.2);
  ok("l+s for another J", m.radiusFor(1, 1, 0.5, 0.5, 5.0) == 4.1);
  ok("l alone for another s", m.radiusFor(1, 1, 1.5, 0.5, 5.0) == 4.0);
  ok("another l or pair: the pair radius", m.radiusFor(1, 2, 0.5, 1.5, 5.0) == 5.0 && m.radiusFor(2, 1, 0.5, 1.5, 6.0) == 6.0);
  ok("edits refused: radius 0, l -1, spin 1/3, pair 0",
     !m.setData(m.index(0, ChannelRadiiModel::RadiusCol), "0") && !m.setData(m.index(0, ChannelRadiiModel::LCol), "-1") &&
         !m.setData(m.index(0, ChannelRadiiModel::SCol), "1/3") && !m.setData(m.index(0, ChannelRadiiModel::PairCol), "0"));
  ok("display 'any' for an open s", m.data(m.index(0, ChannelRadiiModel::SCol), Qt::DisplayRole).toString() == "any");

  std::cout << "3. project file" << std::endl;
  QString text = readAll(QString(AZURE2_TEST_SOURCE_DIR) + "/tests/hard_sphere/models/no_dummies.azr");
  ok("read no_dummies.azr", text.contains("</levels>"));
  const QString block = "<channelRadii>\n1 1 4 J=3/2\n1 2 4.5 s=1/2   # a comment\n</channelRadii>\n";
  QString withBlock = text;
  withBlock.replace("</levels>\n", "</levels>\n" + block);
  const QString in1 = tmp.filePath("in1.azr");
  writeAll(in1, withBlock);
  {
    AZURESetup setup;
    ok("project with the block opens", setup.readProjectFile(in1));
    ChannelRadiiModel *cr = setup.getPairsTab()->getChannelRadiiModel();
    ok("two rows read", cr->rows().size() == 2 && cr->rows().at(0).J == 1.5 && cr->rows().at(1).s == 0.5 &&
                            cr->rows().at(1).radius == 4.5);
    const QString out1 = tmp.filePath("out1.azr");
    ok("written", setup.writeProjectFile(out1));
    const QString written = readAll(out1);
    ok("block written right after </levels>",
       written.contains("</levels>\n<channelRadii>\n1 1 4 J=3/2\n1 2 4.5 s=1/2\n</channelRadii>\n<segmentsData>"));
    AZURESetup again;
    ok("written file opens", again.readProjectFile(out1));
    ok("same rows after the round trip", sameRows(again.getPairsTab()->getChannelRadiiModel()->rows(), cr->rows()));

    // The Pairs tab lists every particle channel of the active levels (2026-11-15) with the radius it takes; an edit
    // writes a row naming the channel's l, s and J, and "Clear own radius" removes it again.
    PairsTab *tab = setup.getPairsTab();
    QWidget *box = tab->findChild<QWidget *>("channelRadiiBox");
    QTableView *view = tab->findChild<QTableView *>("channelRadiiView");
    QPushButton *clear = tab->findChild<QPushButton *>("resetChannelRadiusButton");
    ok("box, view and button present", box && view && clear);
    tab->setPerChannelRadiiOption(false);
    ok("the box shows for a file with rows even with the option off", box && box->isVisibleTo(tab));
    ChannelRadiiListModel *lm = tab->getChannelRadiiListModel();
    auto rowOf = [&](int l, double s, double J) {
      for (int i = 0; i < lm->channels().size(); i++) {
        const ChannelKey &c = lm->channels().at(i);
        if (c.pairKey == 1 && c.l == l && c.s == s && c.J == J) return i;
      }
      return -1;
    };
    const int pJ = rowOf(1, 0.5, 1.5), sHalf = rowOf(0, 0.5, 0.5);
    std::cout << "    " << lm->channels().size() << " channels listed" << std::endl;
    ok("the channels of the levels are listed (p 3/2- and s 1/2+ of pair 1)", pJ >= 0 && sHalf >= 0);
    if (pJ >= 0 && sHalf >= 0 && view && clear) {
      auto cell = [&](int row, int col) { return lm->data(lm->index(row, col), Qt::DisplayRole).toString(); };
      ok("l = 1, J = 3/2: 4 fm from the file's J row", cell(pJ, ChannelRadiiListModel::RadiusCol) == "4" &&
                                                           cell(pJ, ChannelRadiiListModel::SourceCol) == "row: l = 1, J = 3/2",
         cell(pJ, ChannelRadiiListModel::RadiusCol) + " / " + cell(pJ, ChannelRadiiListModel::SourceCol));
      ok("s 1/2+: the pair's 5 fm", cell(sHalf, ChannelRadiiListModel::RadiusCol) == "5" &&
                                        cell(sHalf, ChannelRadiiListModel::SourceCol) == "pair");
      lm->setData(lm->index(pJ, ChannelRadiiListModel::RadiusCol), "4.2");
      ok("an edit adds a row naming l, s and J", cr->rows().size() == 3 && cr->rows().at(2).s == 0.5 &&
                                                     cr->rows().at(2).J == 1.5 && cr->rows().at(2).radius == 4.2 &&
                                                     cell(pJ, ChannelRadiiListModel::SourceCol) == "this channel");
      lm->setData(lm->index(pJ, ChannelRadiiListModel::RadiusCol), "5");
      ok("the pair's radius under a broader row stays a row of its own", cr->rows().size() == 3 &&
                                                                             cr->rows().at(2).radius == 5.0);
      view->selectRow(pJ);
      clear->click();
      ok("Clear own radius: back to the broader row", cr->rows().size() == 2 && cell(pJ, ChannelRadiiListModel::RadiusCol) == "4");
      lm->setData(lm->index(sHalf, ChannelRadiiListModel::RadiusCol), "4.8");
      ok("a channel without a row gets one", cr->rows().size() == 3 && cr->radiusFor(1, 0, 0.5, 0.5, 5.0) == 4.8);
      lm->setData(lm->index(sHalf, ChannelRadiiListModel::RadiusCol), "5");
      ok("... and the pair's radius removes it", cr->rows().size() == 2);
      ok("radius 0 refused", !lm->setData(lm->index(sHalf, ChannelRadiiListModel::RadiusCol), "0"));
    }
    if (qEnvironmentVariableIsSet("AZURE2_CHANNEL_RADII_PNG")) {
      tab->resize(900, 640);
      tab->show();
      app.processEvents();
      tab->grab().save(qEnvironmentVariable("AZURE2_CHANNEL_RADII_PNG"));
    }
  }
  {
    // without the block nothing is written; the block after an <abundances> block is read too
    const QString in2 = tmp.filePath("in2.azr");
    writeAll(in2, text);
    AZURESetup setup;
    ok("plain project opens", setup.readProjectFile(in2));
    PairsTab *tab = setup.getPairsTab();
    QWidget *box = tab->findChild<QWidget *>("channelRadiiBox");
    tab->setPerChannelRadiiOption(false);
    ok("no rows, option off: the box is hidden", box && !box->isVisibleTo(tab));
    tab->setPerChannelRadiiOption(true);
    ok("option on: shown, with every channel", box && box->isVisibleTo(tab) &&
                                                   !tab->getChannelRadiiListModel()->channels().isEmpty());
    tab->setPerChannelRadiiOption(false);
    const QString out2 = tmp.filePath("out2.azr");
    setup.writeProjectFile(out2);
    ok("no block written for a project without one", !readAll(out2).contains("channelRadii"));
    QString both = text;
    both.replace("</levels>\n", "</levels>\n<abundances>\n</abundances>\n" + block);
    const QString in3 = tmp.filePath("in3.azr");
    writeAll(in3, both);
    AZURESetup setup3;
    ok("block after <abundances> is read", setup3.readProjectFile(in3) &&
                                               setup3.getPairsTab()->getChannelRadiiModel()->rows().size() == 2);
    QString bad = text;
    bad.replace("</levels>\n", "</levels>\n<channelRadii>\n1 1 four\n</channelRadii>\n");
    const QString in4 = tmp.filePath("in4.azr");
    writeAll(in4, bad);
    AZURESetup setup4;
    ok("a malformed line makes the read fail", !setup4.readProjectFile(in4));
  }

  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
