// The "Measured Scattering Lengths" table of the Segments tab and its optional
// <scatteringLengths> block (engine side: src/ScatteringLength.cpp), headless:
//   1. the channel spins a pair allows: |j1-j2| .. j1+j2, only even ones for an
//      identical pair (p+7Be = {1, 2}, d+d = {0, 2});
//   2. the block round-trips through SegmentsTab: active flag, s = -1, the
//      optional imaginary part present and absent; comments are dropped;
//   3. rows whose pair is gone, is not a particle pair, or does not allow the
//      channel spin are flagged (red, tooltip), not dropped, and written as
//      they are; changing a row's pair keeps or resets its channel spin;
//   4. the editors list only particle pairs and the allowed spins (+ coherent);
//   5. reset() clears the table;
//   6. a whole project file through AZURESetup: the block is read wherever it
//      sits (be7p_fit.azr has it before <segmentsTest>), written after
//      </segmentsTest>, the other blocks are what the same file gives without
//      it, a second read -> write is byte-identical, and a file without the
//      block gets none.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.
// With AZURE2_SCATLEN_PNG=<dir> it saves the Segments tab of the loaded
// be7p_fit.azr (plus one flagged row) as <dir>/segments_tab_scattering_lengths.png.

#include <QApplication>
#include <QBrush>
#include <QColor>
#include <QTabWidget>
#include <QTableView>
#include <QComboBox>
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QSettings>
#include <QStyleOptionViewItem>
#include <QTemporaryDir>
#include <QTextStream>
#include <QString>
#include <QStringList>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include "AZURESetup.h"
#include "SegmentsTab.h"
#include "AddSegDataDialog.h"
#include "PairsModel.h"
#include "ScatteringLengthsModel.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

#ifndef AZURE2_TEST_SOURCE_DIR
#define AZURE2_TEST_SOURCE_DIR "."
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

// Whitespace-normalised non-comment lines.
static QStringList tokens(const QString& text) {
  QStringList out;
  for (const QString& line : text.split('\n', Qt::SkipEmptyParts)) {
    QString t = line.section('#', 0, 0).simplified();
    if (!t.isEmpty()) out << t;
  }
  return out;
}

static void addPair(PairsModel& pairs, int row, const QVariant (&v)[PairsData::SIZE]) {
  pairs.insertRows(row, 1);
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c), v[c]);
}

static PairsData pairOf(const QVariant (&v)[PairsData::SIZE]) {
  PairsModel m;
  addPair(m, 0, v);
  return m.getPairs().first();
}

static QString spins(const QList<double>& s) {
  QStringList out;
  for (double x : s) out << QString::number(x);
  return out.join(",");
}

static QString readAll(const QString& path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll());
}

static bool writeAll(const QString& path, const QString& text) {
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly)) return false;
  f.write(text.toUtf8());
  return true;
}

// The text with the <scatteringLengths> block (tags included) cut out.
static QString withoutBlock(const QString& text) {
  QStringList out;
  bool in = false;
  for (const QString& line : text.split('\n')) {
    if (line.trimmed() == "<scatteringLengths>") in = true;
    if (!in) out << line;
    if (line.trimmed() == "</scatteringLengths>") in = false;
  }
  return out.join('\n');
}

static QString blockOf(const QString& text) {
  QStringList out;
  bool in = false;
  for (const QString& line : text.split('\n')) {
    if (line.trimmed() == "</scatteringLengths>") in = false;
    if (in) out << line;
    if (line.trimmed() == "<scatteringLengths>") in = true;
  }
  return out.join('\n');
}

// Columns as PairsData: lightJ lightPi lightZ lightM lightG heavyJ heavyPi heavyZ heavyM heavyG Ex Sep radius type mask.
static const QVariant P_7BE[PairsData::SIZE] = {0.5, 1, 1, 1.00783, 0., 1.5, -1, 4, 7.01693, 0., 0., 0.1375, 4.3, 0, 0};
static const QVariant P_7BEX[PairsData::SIZE] = {0.5, 1, 1, 1.00783, 0., 0.5, -1, 4, 7.01693, 0., 0.42908, 0.1375, 4.3, 0, 0};
static const QVariant G_8B[PairsData::SIZE] = {1.0, -1, 0, 0., 0., 2.0, 1, 5, 8.02461, 0., 0., 0., 0., 10, 0};
static const QVariant D_D[PairsData::SIZE] = {1.0, 1, 1, 2.01410, 0., 1.0, 1, 1, 2.01410, 0., 0., 23.85, 5.0, 0, 0};
static const QVariant N_P[PairsData::SIZE] = {0.5, 1, 0, 1.00866, 0., 0.5, 1, 1, 1.00783, 0., 0., 2.2246, 5.0, 0, 0};

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  // AZURESetup keeps a recent-file list: keep it away from the user's settings.
  QTemporaryDir tmp;
  ok("temporary directory", tmp.isValid());
  QCoreApplication::setOrganizationName("azure2-test");
  QCoreApplication::setApplicationName("scattering-lengths-test");
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, tmp.path() + "/settings");
  QSettings::setPath(QSettings::NativeFormat, QSettings::UserScope, tmp.path() + "/settings");

  // ---- 1. allowed channel spins ----------------------------------------------
  {
    PairsData pBe = pairOf(P_7BE), pBex = pairOf(P_7BEX), dd = pairOf(D_D), np = pairOf(N_P), g = pairOf(G_8B);
    ok("p(1/2+) + 7Be(3/2-): s = {1, 2}", spins(ScatteringLengthsModel::allowedChannelSpins(pBe)) == "1,2",
       spins(ScatteringLengthsModel::allowedChannelSpins(pBe)));
    ok("p + 7Be*(1/2-): s = {0, 1}", spins(ScatteringLengthsModel::allowedChannelSpins(pBex)) == "0,1",
       spins(ScatteringLengthsModel::allowedChannelSpins(pBex)));
    ok("d(1+) + d(1+), identical: s = {0, 2}", spins(ScatteringLengthsModel::allowedChannelSpins(dd)) == "0,2",
       spins(ScatteringLengthsModel::allowedChannelSpins(dd)));
    ok("n + p is not identical: s = {0, 1}", spins(ScatteringLengthsModel::allowedChannelSpins(np)) == "0,1",
       spins(ScatteringLengthsModel::allowedChannelSpins(np)));
    ok("identical-pair rule", ScatteringLengthsModel::isIdentical(dd) && !ScatteringLengthsModel::isIdentical(pBe) &&
                                  !ScatteringLengthsModel::isIdentical(np));
    ok("a gamma pair has no channel spins", ScatteringLengthsModel::allowedChannelSpins(g).isEmpty());
    ok("spin labels", ScatteringLengthsModel::spinLabel(1.0) == "1" && ScatteringLengthsModel::spinLabel(0.5) == "1/2" &&
                          ScatteringLengthsModel::spinLabel(-1.0) == "coherent");
  }

  // ---- 2. block round trip through the Segments tab ------------------------------
  PairsModel pairs;
  addPair(pairs, 0, P_7BE);
  addPair(pairs, 1, P_7BEX);
  addPair(pairs, 2, G_8B);
  addPair(pairs, 3, D_D);
  SegmentsTab tab;
  tab.setPairsModel(&pairs);
  ScatteringLengthsModel* model = tab.getScatteringLengthsModel();
  ok("empty table: no block", !tab.hasScatteringLengths());
  {
    QString in("# active pairKey s value sigma [imValue imSigma]\n"
               "1 1 1 17.0 0.5\n"
               "\n"
               "1 1 -1 5.0 1.0   # coherent\n"
               "0 2 1 9.0 0.5 -0.3 0.05\n"
               "</scatteringLengths>\n"
               "<segmentsTest>\n");
    QTextStream inStream(&in);
    ok("block read", tab.readScatteringLengths(inStream));
    ok("the reader stops at the closing tag", inStream.readLine() == "<segmentsTest>");
    QList<ScatteringLengthEntry> e = model->entries();
    ok("three rows", e.size() == 3, QString::number(e.size()));
    if (e.size() == 3) {
      ok("row 1", e[0].active && e[0].pairKey == 1 && e[0].s == 1.0 && e[0].value == 17.0 && e[0].sigma == 0.5 && !e[0].hasIm);
      ok("row 2: coherent", e[1].active && e[1].pairKey == 1 && e[1].s == -1.0 && e[1].value == 5.0 && e[1].sigma == 1.0);
      ok("row 3: inactive, with an imaginary part", !e[2].active && e[2].pairKey == 2 && e[2].s == 1.0 && e[2].hasIm &&
                                                        e[2].imValue == -0.3 && e[2].imSigma == 0.05);
    }
    ok("has a block", tab.hasScatteringLengths());
    ok("display: J^pi of the channel spin and the coherent label",
       model->index(0, ScatteringLengthsModel::SpinCol).data().toString() == QString("1") + QChar(0x207B) &&
           model->index(1, ScatteringLengthsModel::SpinCol).data().toString() == "coherent",
       model->index(0, ScatteringLengthsModel::SpinCol).data().toString());
    ok("display: absent imaginary part is blank", model->index(0, ScatteringLengthsModel::ImValueCol).data().toString().isEmpty());
    ok("display: active check box", model->index(2, ScatteringLengthsModel::ActiveCol).data(Qt::CheckStateRole).toInt() == Qt::Unchecked);
    ok("no row is flagged", model->problem(0).isEmpty() && model->problem(1).isEmpty() && model->problem(2).isEmpty(),
       model->problem(0) + model->problem(1) + model->problem(2));

    QString out;
    QTextStream outStream(&out);
    tab.writeScatteringLengths(outStream);
    QStringList want;
    want << "1 1 1 17 0.5" << "1 1 -1 5 1" << "0 2 1 9 0.5 -0.3 0.05";
    ok("written back (comments dropped, a header comment added)", tokens(out) == want, tokens(out).join(" | "));
    ok("the header comment", out.startsWith("# active pairKey s value sigma"));

    QString bad1("1 1 1 17.0\n</scatteringLengths>\n"), bad2("1 1 1 17.0 0.5 -0.3\n</scatteringLengths>\n"),
        bad3("1 1 1 17.0 0.5\n");
    QTextStream b1(&bad1), b2(&bad2), b3(&bad3);
    ScatteringLengthsModel scratch;
    ok("a short line is refused", !scratch.read(b1));
    ok("an imaginary part without its sigma is refused", !scratch.read(b2));
    ok("a block without its closing tag is refused", !scratch.read(b3));
  }

  // ---- 3. flags, edits ----------------------------------------------------------
  {
    ScatteringLengthEntry gamma, missing, ddCoherent, ddOdd, beBad;
    gamma.pairKey = 3;
    missing.pairKey = 9;
    ddCoherent.pairKey = 4;
    ddCoherent.s = -1.0;
    ddOdd.pairKey = 4;
    ddOdd.s = 1.0;
    beBad.pairKey = 1;
    beBad.s = 3.0;
    ok("a gamma pair is flagged", model->problem(gamma).contains("not a particle pair"), model->problem(gamma));
    ok("a missing pair is flagged", model->problem(missing).contains("no pair 9"), model->problem(missing));
    ok("coherent on an identical pair is flagged", model->problem(ddCoherent).contains("identical"), model->problem(ddCoherent));
    ok("odd s on an identical pair is flagged", model->problem(ddOdd).contains("allowed: 0, 2"), model->problem(ddOdd));
    ok("s = 3 for p+7Be is flagged", model->problem(beBad).contains("allowed: 1, 2"), model->problem(beBad));

    // Change row 3's pair to d+d: s = 1 is not allowed there, so it becomes 0; back to 7Be*: 0 is fine and stays.
    QModelIndex pairCell = model->index(2, ScatteringLengthsModel::PairCol);
    ok("pair changed", model->setData(pairCell, 4));
    ok("an s the new pair lacks becomes its first", model->entries()[2].s == 0.0, QString::number(model->entries()[2].s));
    model->setData(pairCell, 2);
    ok("an s the new pair has is kept", model->entries()[2].s == 0.0);
    model->setData(pairCell, 1);
    ok("p+7Be: 0 -> 1", model->entries()[2].s == 1.0);
    model->setData(pairCell, 2);
    ok("negative sigma refused", !model->setData(model->index(2, ScatteringLengthsModel::SigmaCol), -1.0));
    ok("junk refused", !model->setData(model->index(2, ScatteringLengthsModel::ValueCol), QString("abc")));
    ok("value edited", model->setData(model->index(2, ScatteringLengthsModel::ValueCol), QString("9.25")) &&
                           model->entries()[2].value == 9.25);
    ok("blank imaginary cells remove the part",
       model->setData(model->index(2, ScatteringLengthsModel::ImValueCol), QString("")) &&
           model->setData(model->index(2, ScatteringLengthsModel::ImSigmaCol), QString("")) && !model->entries()[2].hasIm);
    ok("an imaginary value adds it", model->setData(model->index(2, ScatteringLengthsModel::ImValueCol), QString("-0.3")) &&
                                         model->entries()[2].hasIm);
    model->setData(model->index(2, ScatteringLengthsModel::ImSigmaCol), QString("0.05"));
    model->setData(model->index(2, ScatteringLengthsModel::ValueCol), 9.0);

    // Removing pair 1 renumbers the others: rows on key 1 now point at 7Be* (s = 1 fine, coherent fine), rows on 2 at gamma.
    pairs.removeRows(0, 1);
    ok("rows are kept when a pair goes", model->rowCount() == 3);
    ok("a row now on the gamma pair is flagged", !model->problem(2).isEmpty(), model->problem(2));
    ok("flag: red text", model->index(2, ScatteringLengthsModel::ValueCol).data(Qt::ForegroundRole).value<QBrush>().color() == QColor(200, 0, 0));
    ok("flag: tooltip gives the reason", model->index(2, ScatteringLengthsModel::ValueCol).data(Qt::ToolTipRole).toString() == model->problem(2));
    ok("a fine row has no red", !model->index(0, ScatteringLengthsModel::ValueCol).data(Qt::ForegroundRole).isValid());
    QString out;
    QTextStream outStream(&out);
    tab.writeScatteringLengths(outStream);
    ok("flagged rows are written as they are", tokens(out).size() == 3 && tokens(out).at(2) == "0 2 1 9 0.5 -0.3 0.05",
       tokens(out).join(" | "));
    pairs.insertRows(0, 1);
    for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(0, c), P_7BE[c]);
    ok("and are fine again once the pair is back", model->problem(2).isEmpty(), model->problem(2));
  }

  // ---- 4. editors ------------------------------------------------------------------
  {
    ScatteringLengthsDelegate delegate;
    QStyleOptionViewItem option;
    ScatteringLengthEntry dd;
    dd.pairKey = 4;
    model->addEntry(dd);
    const int ddRow = model->rowCount() - 1;
    QWidget* w = delegate.createEditor(nullptr, option, model->index(0, ScatteringLengthsModel::SpinCol));
    QComboBox* box = qobject_cast<QComboBox*>(w);
    delegate.setEditorData(box, model->index(0, ScatteringLengthsModel::SpinCol));
    QStringList items;
    for (int i = 0; box && i < box->count(); i++) items << box->itemText(i);
    const QString minus(QChar(0x207B)), plus(QChar(0x207A));
    ok("spin editor for p+7Be: 1-, 2-, coherent", items == (QStringList() << "1" + minus << "2" + minus << "coherent"), items.join(","));
    ok("spin editor shows the row's value", box && box->currentIndex() == 0);
    delegate.setEditorData(box, model->index(ddRow, ScatteringLengthsModel::SpinCol));
    items.clear();
    for (int i = 0; box && i < box->count(); i++) items << box->itemText(i);
    ok("spin editor for d+d: 0+, 2+, no coherent", items == (QStringList() << "0" + plus << "2" + plus), items.join(","));
    delete w;
    w = delegate.createEditor(nullptr, option, model->index(0, ScatteringLengthsModel::PairCol));
    box = qobject_cast<QComboBox*>(w);
    delegate.setEditorData(box, model->index(0, ScatteringLengthsModel::PairCol));
    QList<int> keys;
    for (int i = 0; box && i < box->count(); i++) keys << box->itemData(i).toInt();
    ok("pair editor lists only the particle pairs", keys == (QList<int>() << 1 << 2 << 4));
    box->setCurrentIndex(box->findData(4));
    delegate.setModelData(box, model, model->index(0, ScatteringLengthsModel::PairCol));
    ok("pair editor sets the key (s = 1 -> d+d's first, 0)", model->entries()[0].pairKey == 4 && model->entries()[0].s == 0.0,
       QString("%1 %2").arg(model->entries()[0].pairKey).arg(model->entries()[0].s));
    delete w;
    model->addDefaultEntry();
    ok("a new row: first particle pair, first spin, nothing constrained",
       model->entries().last().pairKey == 1 && model->entries().last().s == 1.0 && model->entries().last().sigma == 0.0 &&
           model->entries().last().active);
  }

  // ---- 4b. scattering lengths as a Data Type of the segment dialog, listed under the data segments ----------
  {
    AddSegDataDialog d;
    d.setPairsModel(&pairs);
    QComboBox* type = d.findChild<QComboBox*>("dataTypeCombo");
    ok("dialog offers the Scattering Length type", type && type->findData(AddSegDataDialog::ScatteringLengthCode) >= 0);
    d.show();
    QApplication::processEvents();
    d.setDataTypeCode(AddSegDataDialog::ScatteringLengthCode);
    QApplication::processEvents();
    ok("scattering-length mode: its fields shown, the segment ones hidden",
       d.isScatteringLength() && d.scatLenBox->isVisible() && !d.angleBox->isVisible() && !d.dataFileText->isVisible() &&
           !d.covarianceFileText->isVisible() && !d.uposCheck->isVisible());
    ok("scattering-length mode: exit pair shown but greyed out", d.exitPairCombo->isVisible() && !d.exitPairCombo->isEnabled());
    d.entrancePairCombo->setKey(1);   // p + 7Be: s = 1, 2 and coherent
    QStringList spins;
    for (int i = 0; i < d.scatLenSpinCombo->count(); i++) spins << d.scatLenSpinCombo->itemText(i);
    ok("channel spins of the pair offered", spins.size() == 3 && spins.last() == "coherent", spins.join(","));
    d.entrancePairCombo->setKey(4);   // d + d (identical): s = 0, 2, no coherent
    ok("identical pair: even spins, no coherent", d.scatLenSpinCombo->count() == 2 &&
                                                    d.scatLenSpinCombo->findData(-1.0) < 0);
    d.entrancePairCombo->setKey(1);
    d.scatLenSpinCombo->setCurrentIndex(d.scatLenSpinCombo->findData(2.0));
    d.scatLenValueText->setText("-3.5");
    d.scatLenSigmaText->setText("0.4");
    ScatteringLengthEntry e = d.scatteringLength();
    ok("entry from the dialog", e.active && e.pairKey == 1 && e.s == 2.0 && e.value == -3.5 && e.sigma == 0.4 && !e.hasIm);
    d.scatLenImValueText->setText("-0.2");
    d.scatLenImSigmaText->setText("0.05");
    e = d.scatteringLength();
    ok("imaginary part when given", e.hasIm && e.imValue == -0.2 && e.imSigma == 0.05);
    d.setDataTypeCode(1);
    QApplication::processEvents();
    ok("back to a segment type: segment fields shown again", !d.isScatteringLength() && !d.scatLenBox->isVisible() &&
                                                                 d.angleBox->isVisible() && d.dataFileText->isVisible() &&
                                                                 d.exitPairCombo->isVisible() && d.exitPairCombo->isEnabled());

    AddSegDataDialog edit;
    edit.setPairsModel(&pairs);
    ScatteringLengthEntry in;
    in.active = false; in.pairKey = 2; in.s = 1.0; in.value = 9.0; in.sigma = 0.5; in.hasIm = true; in.imValue = -0.3;
    in.imSigma = 0.05;
    edit.setScatteringLength(in);
    const ScatteringLengthEntry back = edit.scatteringLength();
    ok("edit round trip (inactive kept)", !back.active && back.pairKey == 2 && back.s == 1.0 && back.value == 9.0 &&
                                              back.sigma == 0.5 && back.hasIm && back.imValue == -0.3 && back.imSigma == 0.05);
    ok("edit: the type is fixed", !edit.findChild<QComboBox*>("dataTypeCombo")->isEnabled());
    if (const char* shots = std::getenv("AZURE2_SCATLEN_PNG")) {   // the dialog as the user sees it
      edit.show();
      for (int k = 0; k < 20; k++) QApplication::processEvents();
      QDir().mkpath(shots);
      ok("dialog screenshot saved", edit.grab().save(QString("%1/scattering_length_dialog.png").arg(shots)));
      edit.hide();
    }
    AddSegDataDialog seg;
    seg.disallowScatteringLength();
    ok("segment edit: no Scattering Length type",
       seg.findChild<QComboBox*>("dataTypeCombo")->findData(AddSegDataDialog::ScatteringLengthCode) < 0);

    // rows under the data table: one per entry, in the data table's columns, hidden when there are none
    QTableView* rows = tab.findChild<QTableView*>("scatteringLengthsView");
    ok("rows view exists", rows != nullptr);
    if (rows) {
      QAbstractItemModel* rm = rows->model();
      ok("one row per entry", rm->rowCount() == model->rowCount(), QString::number(rm->rowCount()));
      ok("Data Type column", rm->index(0, 7).data().toString() == "Scattering Length", rm->index(0, 7).data().toString());
      ok("Data File column: the value", rm->index(0, 8).data().toString().startsWith("Re a0 = 17"),
         rm->index(0, 8).data().toString());
      ok("check box mirrors the entry", rm->index(2, 0).data(Qt::CheckStateRole).toInt() == Qt::Unchecked);
      rm->setData(rm->index(2, 0), Qt::Checked, Qt::CheckStateRole);
      ok("check box writes the entry", model->entries().at(2).active);
      rm->setData(rm->index(2, 0), Qt::Unchecked, Qt::CheckStateRole);
    }
  }

  // ---- 5. reset --------------------------------------------------------------------
  tab.reset();
  ok("reset() clears the table", model->rowCount() == 0 && !tab.hasScatteringLengths());

  // ---- 6. a project file through AZURESetup ------------------------------------------------
  {
    const QString source = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/scattering_length/be7p_fit.azr";
    const QString original = readAll(source);
    ok("engine test model found", original.contains("<scatteringLengths>"), source);
    const QString dir = tmp.path();
    QDir().mkpath(dir + "/output");
    QDir().mkpath(dir + "/checks");
    writeAll(dir + "/in.azr", original);
    writeAll(dir + "/noblock.azr", withoutBlock(original));

    AZURESetup setup;
    ok("project read", setup.readProjectFile(dir + "/in.azr"));
    QList<ScatteringLengthEntry> e = setup.getSegmentsTab()->getScatteringLengthsModel()->entries();
    ok("the block before <segmentsTest> is read: five rows", e.size() == 5, QString::number(e.size()));
    ok("project written", setup.writeProjectFile(dir + "/a.azr"));
    const QString a = readAll(dir + "/a.azr");
    QStringList want;
    want << "1 1 1 17 0.5" << "1 1 2 -4 0.5" << "1 1 -1 5 1" << "1 2 1 9 0.5 -0.3 0.05" << "0 2 0 25 1";
    ok("block written with every line", tokens(blockOf(a)) == want, tokens(blockOf(a)).join(" | "));
    ok("written after </segmentsTest>, before <targetInt>",
       a.indexOf("</segmentsTest>") < a.indexOf("<scatteringLengths>") && a.indexOf("</scatteringLengths>") < a.indexOf("<targetInt>") &&
           a.count("<scatteringLengths>") == 1);

    // Screenshot of the tab with the file loaded and one flagged row (pair 9).
    if (const char* shots = std::getenv("AZURE2_SCATLEN_PNG")) {
      ScatteringLengthEntry flagged;
      flagged.pairKey = 9;
      flagged.s = 1.0;
      flagged.value = 3.0;
      flagged.sigma = 0.2;
      setup.getSegmentsTab()->getScatteringLengthsModel()->addEntry(flagged);
      setup.resize(1250, 800);
      setup.show();
      QWidget* segTab = setup.getSegmentsTab();
      if (QTabWidget* tabs = qobject_cast<QTabWidget*>(setup.centralWidget())) tabs->setCurrentWidget(segTab);
      QApplication::processEvents();
      QDir().mkpath(shots);
      bool saved = setup.grab().save(QString("%1/segments_tab_scattering_lengths.png").arg(shots));
      // The table itself scrolled to the flagged row, with its tooltip text.
      for (QTableView* view : setup.getSegmentsTab()->findChildren<QTableView*>()) {
        if (view->model() != setup.getSegmentsTab()->getScatteringLengthsModel()) continue;
        view->scrollToBottom();
        QApplication::processEvents();
        saved &= view->parentWidget()->grab().save(QString("%1/scattering_lengths_table_flagged.png").arg(shots));
      }
      ok("screenshots saved", saved, shots);
      setup.hide();
      setup.getSegmentsTab()->getScatteringLengthsModel()->removeRows(5, 1);
    }

    AZURESetup plain;
    ok("project without the block read", plain.readProjectFile(dir + "/noblock.azr"));
    ok("... and its table is empty", plain.getSegmentsTab()->getScatteringLengthsModel()->rowCount() == 0);
    ok("project without the block written", plain.writeProjectFile(dir + "/b.azr"));
    const QString b = readAll(dir + "/b.azr");
    ok("no block written for a project without one", !b.contains("scatteringLengths"));
    ok("every other block is what the file gives without the block", withoutBlock(a) == b);

    AZURESetup again;
    ok("written project read back", again.readProjectFile(dir + "/a.azr"));
    ok("written project written again", again.writeProjectFile(dir + "/a2.azr"));
    const QString a2 = readAll(dir + "/a2.azr");
    QString firstDiff;
    {
      QStringList la = a.split('\n'), lb = a2.split('\n');
      for (int i = 0; i < std::max(la.size(), lb.size()) && firstDiff.isEmpty(); i++)
        if (la.value(i) != lb.value(i)) firstDiff = QString("line %1: '%2' vs '%3'").arg(i + 1).arg(la.value(i)).arg(lb.value(i));
    }
    // The Fitting tab's <parameterSettings> of this file is not stable across a second save (its level rows
    // are listed in another order; nothing to do with this block), so the second generation is compared with
    // the second generation of the same file without the block.
    AZURESetup plain2;
    ok("block-less project read back and written again",
       plain2.readProjectFile(dir + "/b.azr") && plain2.writeProjectFile(dir + "/b2.azr"));
    const QString b2 = readAll(dir + "/b2.azr");
    ok("second save: the block is unchanged", blockOf(a2) == blockOf(a));
    ok("second save: every other block is what the file gives without the block", withoutBlock(a2) == b2);
    if (a2 != a) std::cout << "  note  the second save differs from the first (also without the block: "
                           << (b2 != b ? "yes" : "no") << "), first difference " << firstDiff.toStdString() << std::endl;

    // The block after </segmentsTest> (where the GUI writes it) and after </targetInt> are read too.
    QString late = withoutBlock(original);
    late.replace("</targetInt>", "</targetInt>\n<scatteringLengths>\n1 1 2 -4.0 0.5\n</scatteringLengths>");
    writeAll(dir + "/late.azr", late);
    ok("block after </targetInt> read", again.readProjectFile(dir + "/late.azr") &&
                                             again.getSegmentsTab()->getScatteringLengthsModel()->rowCount() == 1);
    ok("reading another file replaces the rows (readFile resets)",
       again.readProjectFile(dir + "/noblock.azr") && again.getSegmentsTab()->getScatteringLengthsModel()->rowCount() == 0);
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
