// The light particle's excitation energy (column 32 of a <levels> line, 2026-11-05), GUI side (headless):
//   1. PairsData::lightExcitation is model column 15 (PairsModel::LightExcitationColumn), the compound-nucleus
//      column moved to 16; header text; setData/data round trip;
//   2. labels: an excited light particle is marked "*(E)" in the pair and reaction labels;
//   3. identical particles (PairsModel::isIdentical, used by the scattering-lengths table and the level
//      channels): all per-particle inputs equal, including excitation -- 12C(2+)+12C(2+) entered as e1 = e2 is,
//      e1 = 0 / e2 = 2E* is not; a level of an identical pair is offered l + s even channels only;
//   4. LevelsTab reads a 32nd column into the pair and writes it back only when it is nonzero: a project
//      without it is written with 31 columns, one with it keeps it on that pair's lines only;
//   5. the Add/Edit Pair dialog has the field (disabled for gamma and beta pairs), and PairsTab::editPair
//      stores it -- and the EC multipolarity mask, which before 2026-11-05 was compared against the pair type;
//   6. the mass-table check of a pair does not depend on the excitation.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QLineEdit>
#include <QComboBox>
#include <QTextStream>
#include <QString>
#include <QStringList>
#include <cmath>
#include <iostream>
#include "AddPairDialog.h"
#include "PairsTab.h"
#include "PairsModel.h"
#include "LevelsTab.h"
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "MassTable.h"
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

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QList<QStringList> levelLines(const QString& text) {
  QList<QStringList> out;
  for (const QString& line : text.split('\n', Qt::SkipEmptyParts)) {
    QStringList t = line.simplified().split(' ', Qt::SkipEmptyParts);
    if (t.size() >= 31) out << t;
  }
  return out;
}

static PairsData pair(double j1, int pi1, int z1, double m1, double j2, int pi2, int z2, double m2, double e2, double S,
                      double e1, int type = 0) {
  PairsData p{};
  p.lightJ = j1; p.lightPi = pi1; p.lightZ = z1; p.lightM = m1;
  p.heavyJ = j2; p.heavyPi = pi2; p.heavyZ = z2; p.heavyM = m2;
  p.excitationEnergy = e2; p.seperationEnergy = S; p.channelRadius = 5.43; p.pairType = type; p.lightExcitation = e1;
  return p;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  const double E2P = 4.43891;

  // ---- 1. the model column -------------------------------------------------
  PairsTab pairsTab;
  PairsModel& pairs = *pairsTab.getPairsModel();
  ok("PairsData::SIZE is 16, light excitation column 15, nucleus column 16",
     PairsData::SIZE == 16 && PairsModel::LightExcitationColumn == 15 && PairsModel::NucleusColumn == 16);
  pairsTab.addPair(pair(0, 1, 2, 4.0026, 0, 1, 6, 12., 0., 7.16191, 0.), 0, true);
  pairsTab.addPair(pair(2, 1, 2, 4.0026, 0, 1, 6, 12., 0., 7.16191, E2P), 1, true);
  ok("columnCount includes the new column", pairs.columnCount(QModelIndex()) == 17);
  ok("header", pairs.headerData(PairsModel::LightExcitationColumn, Qt::Horizontal, Qt::DisplayRole).toString() == "Light\nExcitation" &&
     pairs.headerData(10, Qt::Horizontal, Qt::DisplayRole).toString() == "Heavy\nExcitation");
  ok("addPair(PairsData) stores the light excitation",
     pairs.getPairs().at(1).lightExcitation == E2P && pairs.getPairs().at(0).lightExcitation == 0.,
     QString::number(pairs.getPairs().at(1).lightExcitation));
  ok("data(EditRole) and data(DisplayRole) of column 15",
     pairs.data(pairs.index(1, 15), Qt::EditRole).toDouble() == E2P && pairs.data(pairs.index(1, 15), Qt::DisplayRole).toDouble() == E2P);
  ok("isPair tells the two apart", pairs.isPair(pair(2, 1, 2, 4.0026, 0, 1, 6, 12., 0., 7.16191, E2P)) == 1 &&
     pairs.isPair(pair(2, 1, 2, 4.0026, 0, 1, 6, 12., 0., 7.16191, 1.0)) == -1);

  // ---- 2. labels ---------------------------------------------------------------
  {
    QString ground = pairs.getPairLabel(pairs.getPairs().at(0));
    QString excited = pairs.getPairLabel(pairs.getPairs().at(1));
    ok("ground-state label unchanged", ground == QString::fromUtf8("¹²C+α [0.000 MeV]"), ground);
    ok("excited light particle marked", excited == QString::fromUtf8("¹²C+α*(4.439) [0.000 MeV]"), excited);
    QString reaction = PairsModel::plainText(pairs.getReactionLabel(pairs.getPairs().at(0), pairs.getPairs().at(1)));
    ok("reaction label marks the excited ejectile", reaction == QString::fromUtf8("¹²C(α,α*(4.439))¹²C [0.000 MeV]"), reaction);
    QString single = PairsModel::plainText(pairs.getParticleLabel(pairs.getPairs().at(1), 0));
    ok("light-particle column label marked", single == QString::fromUtf8("α*(4.439)"), single);
  }

  // ---- 3. identical particles -----------------------------------------------
  {
    PairsData mutual = pair(2, 1, 6, 12., 2, 1, 6, 12., E2P, 13.933, E2P);
    PairsData oldWay = pair(2, 1, 6, 12., 2, 1, 6, 12., 2 * E2P, 13.933, 0.);
    PairsData ground = pair(0, 1, 6, 12., 0, 1, 6, 12., 0., 13.933, 0.);
    PairsData oneExcited = pair(0, 1, 6, 12., 2, 1, 6, 12., E2P, 13.933, 0.);
    PairsData gDiff = ground; gDiff.lightG = 1.0;
    ok("12C(2+)+12C(2+) with e1 = e2 is identical", PairsModel::isIdentical(mutual) && ScatteringLengthsModel::isIdentical(mutual));
    ok("the old way (e1 = 0, e2 = 2E*) is not", !PairsModel::isIdentical(oldWay));
    ok("12C+12C is, 12C+12C(2+) is not", PairsModel::isIdentical(ground) && !PairsModel::isIdentical(oneExcited));
    ok("different g-factors are not", !PairsModel::isIdentical(gDiff));
    ok("identical spin-2 bosons: s = 0, 2, 4 for l = 0", ScatteringLengthsModel::allowedChannelSpins(mutual) == (QList<double>() << 0. << 2. << 4.));
  }

  // ---- 3b. the channels the GUI offers for a level of an identical pair ----------
  {
    PairsTab tab;
    LevelsTab levelsTab;
    levelsTab.setPairsModel(tab.getPairsModel());
    tab.addPair(pair(2, 1, 6, 12., 2, 1, 6, 12., E2P, 13.933, E2P), 0, true);
    LevelsData four = {1, 0, 4.0, 1, 26.0, 12, 24};
    QList<ChannelsData> chans = levelsTab.calculateChannels(four, -1);
    bool even = !chans.isEmpty();
    for (const ChannelsData& c : chans) even &= (c.lValue + int(c.sValue + 0.5)) % 2 == 0;
    // 4+ of two spin-2 bosons, l <= maxL: s = 0 (l 4), s = 2 (l 2, 4, 6), s = 4 (l 0, 2, 4, 6, 8) with l <= maxL
    ok("a 4+ level of 12C(2+)+12C(2+) is offered only l + s even channels", even, QString::number(chans.size()));
    PairsTab tab2;
    LevelsTab levels2;
    levels2.setPairsModel(tab2.getPairsModel());
    tab2.addPair(pair(2, 1, 6, 12., 2, 1, 6, 12., 2 * E2P, 13.933, 0.), 0, true);
    bool odd = false;
    for (const ChannelsData& c : levels2.calculateChannels(four, -1)) odd |= (c.lValue + int(c.sValue + 0.5)) % 2 != 0;
    ok("entered the old way (e1 = 0, e2 = 2E*) it is offered l + s odd ones too", odd);
  }

  // ---- 4. the 32nd column through LevelsTab -----------------------------------
  {
    PairsTab tab;
    LevelsTab levelsTab;
    levelsTab.setPairsModel(tab.getPairsModel());
    QObject::connect(&levelsTab, SIGNAL(readNewPair(PairsData, int, bool)), &tab, SLOT(addPair(PairsData, int, bool)));
    QObject::connect(&levelsTab, SIGNAL(readExistingPair(PairsData, int, bool)), &tab, SLOT(editPair(PairsData, int, bool)));
    const QString block =
        "2 1 13 1 1 1 0 4 1 1 0 40000 0 1 0 1 0 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0\n"
        "2 1 13 1 1 2 4 0 1 1 0 20000 2 1 0 1 0 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0 4.43891\n"
        "\n"
        "1 -1 13.6 1 1 1 0 2 2 1 0 60000 0 1 0 1 0 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0 0\n"
        "1 -1 13.6 1 1 2 4 2 2 1 0 15000 2 1 0 1 0 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0 4.43891\n"
        "</levels>\n";
    QString in(block);
    QTextStream inStream(&in);
    ok("a block with 32-column lines is read", levelsTab.readNuclearFile(inStream));
    QList<PairsData> read = tab.getPairsModel()->getPairs();
    ok("two pairs, pair 2 with the light excitation, pair 1 without (an explicit 0 is 0)",
       read.size() == 2 && read.at(0).lightExcitation == 0. && read.at(1).lightExcitation == E2P,
       QString("%1 pairs").arg(read.size()));
    QString out;
    QTextStream outStream(&out);
    levelsTab.writeNuclearFile(outStream);
    QList<QStringList> lines = levelLines(out);
    bool layout = lines.size() == 4;
    for (const QStringList& t : lines) {
      if (t.at(5) == "1") layout &= t.size() == 31;
      else layout &= t.size() == 32 && std::fabs(t.at(31).toDouble() - E2P) < 1e-12;
    }
    ok("written back: 32 columns on pair 2's lines only (an explicit 0 is not written)", layout, out);

    PairsTab plainTab;
    LevelsTab plainLevels;
    plainLevels.setPairsModel(plainTab.getPairsModel());
    QObject::connect(&plainLevels, SIGNAL(readNewPair(PairsData, int, bool)), &plainTab, SLOT(addPair(PairsData, int, bool)));
    QObject::connect(&plainLevels, SIGNAL(readExistingPair(PairsData, int, bool)), &plainTab, SLOT(editPair(PairsData, int, bool)));
    QString in2 =
        "2 1 13 1 1 1 0 4 1 1 0 40000 0 1 0 1 0 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0\n"
        "2 1 13 1 1 2 4 0 1 1 0 20000 0 1 2 1 4.43891 4.0026 12 2 6 7.16191 7.16191 0 0 0 0 5.43 0 0 0\n"
        "</levels>\n";
    QTextStream inStream2(&in2);
    ok("a 31-column block is read", plainLevels.readNuclearFile(inStream2));
    QString out2;
    QTextStream outStream2(&out2);
    plainLevels.writeNuclearFile(outStream2);
    bool all31 = true;
    for (const QStringList& t : levelLines(out2)) all31 &= t.size() == 31;
    ok("and written with 31 columns", all31 && levelLines(out2).size() == 2, out2);
  }

  // ---- 5. the pair dialog and editPair ---------------------------------------------
  {
    AddPairDialog dialog;
    ok("the dialog has the light-excitation field, enabled for a particle pair",
       dialog.lightExcitationText != nullptr && dialog.lightExcitationText->isEnabled());
    dialog.lightExcitationText->setText("1.5");
    dialog.pairTypeCombo->setCurrentIndex(1);
    ok("a gamma pair disables and clears it", !dialog.lightExcitationText->isEnabled() &&
       dialog.lightExcitationText->text().toDouble() == 0.);
    dialog.pairTypeCombo->setCurrentIndex(2);
    ok("a beta pair too", !dialog.lightExcitationText->isEnabled());
    dialog.pairTypeCombo->setCurrentIndex(0);
    ok("a particle pair enables it again", dialog.lightExcitationText->isEnabled());

    PairsData edited = pairs.getPairs().at(0);
    edited.lightExcitation = 2.5;
    edited.ecMultMask = 5;
    pairsTab.editPair(edited, 0, true);
    ok("editPair stores the light excitation", pairs.getPairs().at(0).lightExcitation == 2.5);
    ok("editPair stores a new EC mask on a type-0 pair whose old mask was 0", pairs.getPairs().at(0).ecMultMask == 5,
       QString::number(pairs.getPairs().at(0).ecMultMask));
  }

  // ---- 6. the mass table -----------------------------------------------------------
  {
    const MassTable& mt = MassTable::instance();
    if (!mt.isLoaded()) {
      std::cout << "  skip  mass table not loaded" << std::endl;
    } else {
      PairsData g = pair(0, 1, 2, 4.0026032541, 0, 1, 6, 12., 0., 7.16191, 0.);
      PairsData x = g;
      x.lightExcitation = E2P;
      const PairMassCheck a = mt.check(g), b = mt.check(x);
      ok("the mass-table check ignores the excitation (S is between ground states)",
         a.hasS && b.hasS && a.sTable == b.sTable && a.dSkeV == b.dSkeV && a.massOffset == b.massOffset);
    }
  }

  std::cout << (fails ? "FAIL" : "PASS") << ": light_excitation (" << fails << " failures)" << std::endl;
  return fails ? 1 : 0;
}
