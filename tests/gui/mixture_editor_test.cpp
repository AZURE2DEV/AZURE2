// Headless checks on the mixture-block editor of the segment dialogs.
//
// A segment line's trailing "mixture <w_self> <n> {"file" <w>}..." block is
// carried on the row as verbatim text beside any "sqrtshift" block.  The
// editor must read it, leave the other block alone, write it back last (as the
// engine's writers do), keep a block it cannot parse verbatim, and refuse the
// observables the engine refuses.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.
// With AZURE2_MIXTURE_EDITOR_PNG=<stem> it also saves the data dialog as
// <stem>_data.png for a layout review.

#include <QApplication>
#include <QString>
#include <QStringList>
#include <QTableWidgetItem>
#include <QScreen>
#include <QRect>
#include <cstdlib>
#include <iostream>
#include "MixtureEditor.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"
#include "PairsModel.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if (!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) fails++;
}

static void addPair(PairsModel& pairs, int row, const QVariant* v) {
  pairs.insertRows(row, 1, QModelIndex());
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c, QModelIndex()), v[c], Qt::EditRole);
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  const QString junghans =
      "mixture 0.9048 2 \"tables/22Ne_ntot_endfb8_tendl2015.dat\" 0.0925 \"tables/21Ne_ntot_endfb8_tendl2015.dat\" 0.0027";

  // ---- tokens and blocks ---------------------------------------------------
  std::cout << "1. token helpers" << std::endl;
  {
    QStringList t = MixtureEditor::splitTokens("a \"b c\" d  \"e\"");
    ok("quoted token kept whole", t.size() == 4 && t[1] == "b c" && t[3] == "e", t.join("|"));

    bool present; double w; QList<MixtureEditor::Table> tables; QString other, err;
    ok("no block parses", MixtureEditor::parseBlock("sqrtshift 0.003 0.001 0", present, w, tables, other) && !present
                              && other == "sqrtshift 0.003 0.001 0", other);
    ok("the Junghans block parses", MixtureEditor::parseBlock(junghans, present, w, tables, other) && present
                                        && w == 0.9048 && tables.size() == 2 && other.isEmpty()
                                        && tables[1].file == "tables/21Ne_ntot_endfb8_tendl2015.dat" && tables[1].weight == 0.0027,
       other);
    ok("formatted back verbatim", MixtureEditor::formatBlock(w, tables) == junghans, MixtureEditor::formatBlock(w, tables));
    ok("sqrtshift before the block is kept aside",
       MixtureEditor::parseBlock("sqrtshift 0.001 0.002 1 " + junghans, present, w, tables, other) && present
           && other == "sqrtshift 0.001 0.002 1" && tables.size() == 2, other);
    ok("sqrtshift after the block is kept aside",
       MixtureEditor::parseBlock(junghans + " sqrtshift 0.001 0.002 1", present, w, tables, other) && present
           && other == "sqrtshift 0.001 0.002 1" && tables.size() == 2, other);
    ok("a path with a space survives", MixtureEditor::parseBlock("mixture 0.5 1 \"my tables/a.dat\" 0.5", present, w, tables, other)
                                            && tables.size() == 1 && tables[0].file == "my tables/a.dat");
    ok("a short block is refused with a reason",
       !MixtureEditor::parseBlock("mixture 0.9 2 \"a.dat\" 0.1", present, w, tables, other, &err) && !err.isEmpty(), err);
    ok("a non-numeric weight is refused", !MixtureEditor::parseBlock("mixture 0.9 1 \"a.dat\" x", present, w, tables, other, &err), err);
  }

  // ---- the widget ----------------------------------------------------------
  std::cout << "\n2. the editor" << std::endl;
  {
    MixtureEditor ed;
    ok("starts disabled with no block", !ed.isMixtureEnabled() && ed.tokens().isEmpty(), ed.tokens());
    ed.setTokens("sqrtshift 0.003 0.001 0 " + junghans);
    ok("enabled by a block", ed.isMixtureEnabled());
    ok("two table rows", ed.tablesTable->rowCount() == 2);
    ok("w_self shown", ed.selfWeightText->text() == "0.9048", ed.selfWeightText->text());
    ok("sum of weights shown", ed.sumLabel->text().contains("1"), ed.sumLabel->text());
    ok("round trip puts the mixture block last", ed.tokens() == "sqrtshift 0.003 0.001 0 " + junghans, ed.tokens());
    ok("valid as loaded", ed.validate());

    ed.enableCheck->setChecked(false);
    ok("unchecked: only the other block remains", ed.tokens() == "sqrtshift 0.003 0.001 0", ed.tokens());
    ed.enableCheck->setChecked(true);
    ed.selfWeightText->setText("0.5");
    ed.tablesTable->item(1, 1)->setText("0.25");
    ok("edits are written", ed.tokens() == "sqrtshift 0.003 0.001 0 mixture 0.5 2 \"tables/22Ne_ntot_endfb8_tendl2015.dat\" 0.0925 "
                                           "\"tables/21Ne_ntot_endfb8_tendl2015.dat\" 0.25", ed.tokens());
    ed.tablesTable->item(0, 0)->setText("");
    QString why;
    ok("an empty file is refused", !ed.validate(&why) && why.contains("file"), why);
    ed.tablesTable->item(0, 0)->setText("t.dat");
    ed.selfWeightText->setText("abc");
    ok("a non-numeric w_self is refused", !ed.validate(&why), why);

    ed.setTokens("");
    ok("cleared by an empty row", !ed.isMixtureEnabled() && ed.tokens().isEmpty() && ed.tablesTable->rowCount() == 0, ed.tokens());

    ed.setTokens("mixture 0.9 x");
    ok("an unreadable block is kept verbatim and locked",
       ed.tokens() == "mixture 0.9 x" && !ed.isMixtureEnabled() && !ed.enableCheck->isEnabled() && ed.rawLabel->isVisibleTo(&ed),
       ed.tokens());
    ed.setTokens(junghans);
    ok("a readable block unlocks it again", ed.isMixtureEnabled() && ed.enableCheck->isEnabled());
  }

  // ---- in the dialogs ------------------------------------------------------
  std::cout << "\n3. the segment dialogs" << std::endl;
  {
    PairsModel pairs;
    const QVariant p1[PairsData::SIZE] = {0.0, 1, 2, 4.0026, 0., 2.5, 1, 8, 16.9991, 0., 0., 7.347, 1.0, 0, 0};
    const QVariant p2[PairsData::SIZE] = {0.5, 1, 0, 1.00866, 0., 0.0, 1, 10, 19.9924, 0., 0., 7.347, 1.0, 0, 0};
    addPair(pairs, 0, p1); addPair(pairs, 1, p2);

    AddSegDataDialog d;
    d.setPairsModel(&pairs);
    d.mixtureEditor->setTokens(junghans);
    d.setDataTypeCode(0);
    QString why;
    ok("data dialog: angle-integrated accepts a mixture", d.mixtureCompatible(&why), why);
    d.setDataTypeCode(4);
    ok("data dialog: c.m. differential accepts a mixture", d.mixtureCompatible(&why), why);
    d.setDataTypeCode(2);
    ok("data dialog: a phase shift refuses it", !d.mixtureCompatible(&why) && !why.isEmpty(), why);
    d.setDataTypeCode(0);
    d.advancedModeCheck->setChecked(true);
    d.operationCombo->setCurrentIndex(1);
    ok("data dialog: a ratio refuses it", !d.mixtureCompatible(&why) && why.contains("ratio"), why);
    d.operationCombo->setCurrentIndex(0);
    ok("data dialog: a sum accepts it", d.mixtureCompatible(&why), why);
    ok("data dialog hands the block back", d.mixtureEditor->tokens() == junghans, d.mixtureEditor->tokens());

    if (const char* stem = std::getenv("AZURE2_MIXTURE_EDITOR_PNG")) {
      d.advancedModeCheck->setChecked(false);
      d.show();   // sized to its content within the screen by its showEvent
      QApplication::processEvents();
      d.grab().save(QString::fromLocal8Bit(stem) + "_data.png");
      {
        const QRect a = d.screen() ? d.screen()->availableGeometry() : QRect();
        std::cout << "  dialog " << d.width() << "x" << d.height() << " (max " << d.maximumWidth() << "x" << d.maximumHeight()
                  << "), screen available " << a.width() << "x" << a.height() << std::endl;
      }
      d.hide();
    }

    AddSegTestDialog t;
    t.setPairsModel(&pairs);
    t.mixtureEditor->setTokens(junghans);
    t.setDataTypeCode(5);
    ok("test dialog: c.m. differential (5) accepts a mixture", t.mixtureCompatible(&why), why);
    t.setDataTypeCode(4);
    ok("test dialog: total capture (4) refuses it", !t.mixtureCompatible(&why), why);
    t.setDataTypeCode(0);
    ok("test dialog hands the block back", t.mixtureEditor->tokens() == junghans, t.mixtureEditor->tokens());
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
