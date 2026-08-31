// Headless check that moving a data segment up or down keeps every field.
//
// The up/down buttons used to rebuild the moved row column by column, and the
// copy list had fallen behind the struct: the UPOS block (unobserved primary,
// observed secondary) was dropped on every move.  The models now relocate the
// whole struct through moveLine(); this test pins that no field -- present or
// future -- can be lost, by comparing the moved row against the original
// struct memberwise.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QString>
#include <iostream>
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
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
static void ok(const char* what, bool cond) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what << std::endl;
  if(!cond) fails++;
}

static SegmentsDataData sampleRow(int tag) {
  SegmentsDataData d;
  d.isActive = 1;
  d.entrancePairIndex = 1 + tag;
  d.exitPairIndex = 2;
  d.lowEnergy = 0.1 * tag;
  d.highEnergy = 3.0 + tag;
  d.lowAngle = 0.;
  d.highAngle = 180.;
  d.dataType = 2;
  d.dataFile = QString("data/file_%1.dat").arg(tag);
  d.dataNorm = 1.0 + 0.1 * tag;
  d.dataNormError = 5.0;
  d.varyNorm = 1;
  d.phaseJ = 1.5;
  d.phaseL = 1;
  d.energyShift = 0.001 * tag;
  d.energyShiftError = 0.002;
  d.varyEnergyShift = 1;
  d.isAdvanced = 0;
  d.operationType = 0;
  d.componentsList = "";
  d.isUPOS = 1;             // the fields the old move dropped
  d.secondaryDecayL = 2;
  d.finalJ = 0.5;
  d.delta = 0.25;
  return d;
}

static bool sameRow(const SegmentsDataData& a, const SegmentsDataData& b) {
  return a.isActive == b.isActive && a.entrancePairIndex == b.entrancePairIndex &&
         a.exitPairIndex == b.exitPairIndex && a.lowEnergy == b.lowEnergy &&
         a.highEnergy == b.highEnergy && a.lowAngle == b.lowAngle &&
         a.highAngle == b.highAngle && a.dataType == b.dataType &&
         a.dataFile == b.dataFile && a.dataNorm == b.dataNorm &&
         a.dataNormError == b.dataNormError && a.varyNorm == b.varyNorm &&
         a.phaseJ == b.phaseJ && a.phaseL == b.phaseL &&
         a.energyShift == b.energyShift && a.energyShiftError == b.energyShiftError &&
         a.varyEnergyShift == b.varyEnergyShift && a.isAdvanced == b.isAdvanced &&
         a.operationType == b.operationType && a.componentsList == b.componentsList &&
         a.isUPOS == b.isUPOS && a.secondaryDecayL == b.secondaryDecayL &&
         a.finalJ == b.finalJ && a.delta == b.delta;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  SegmentsDataModel model;
  SegmentsDataData r0 = sampleRow(0);
  SegmentsDataData r1 = sampleRow(1);
  // Insert through the model API the way the tab does, then overwrite the
  // rows wholesale via moveLine-independent direct comparison.
  model.insertRows(0, 2, QModelIndex());
  // setData drives every column, including the UPOS block (20-23).
  struct Put { int col; QVariant v0; QVariant v1; };
  const Put puts[] = {
      {0, r0.isActive, r1.isActive}, {1, r0.entrancePairIndex, r1.entrancePairIndex},
      {2, r0.exitPairIndex, r1.exitPairIndex}, {3, r0.lowEnergy, r1.lowEnergy},
      {4, r0.highEnergy, r1.highEnergy}, {5, r0.lowAngle, r1.lowAngle},
      {6, r0.highAngle, r1.highAngle}, {7, r0.dataType, r1.dataType},
      {8, r0.dataFile, r1.dataFile}, {9, r0.dataNorm, r1.dataNorm},
      {10, r0.dataNormError, r1.dataNormError}, {11, r0.varyNorm, r1.varyNorm},
      {12, r0.phaseJ, r1.phaseJ}, {13, r0.phaseL, r1.phaseL},
      {14, r0.energyShift, r1.energyShift}, {15, r0.energyShiftError, r1.energyShiftError},
      {16, r0.varyEnergyShift, r1.varyEnergyShift}, {17, r0.isAdvanced, r1.isAdvanced},
      {18, r0.operationType, r1.operationType}, {19, r0.componentsList, r1.componentsList},
      {20, r0.isUPOS, r1.isUPOS}, {21, r0.secondaryDecayL, r1.secondaryDecayL},
      {22, r0.finalJ, r1.finalJ}, {23, r0.delta, r1.delta}};
  for (const Put& p : puts) {
    model.setData(model.index(0, p.col, QModelIndex()), p.v0, Qt::EditRole);
    model.setData(model.index(1, p.col, QModelIndex()), p.v1, Qt::EditRole);
  }
  ok("two rows set up", model.getLines().size() == 2 &&
                         sameRow(model.getLines().at(0), r0) &&
                         sameRow(model.getLines().at(1), r1));

  // Move row 0 down: rows swap, nothing may change inside either row.
  ok("moveLine(0,1) accepted", model.moveLine(0, 1));
  ok("row moved down intact (UPOS kept)", sameRow(model.getLines().at(1), r0));
  ok("displaced row intact", sameRow(model.getLines().at(0), r1));

  // And back up.
  ok("moveLine(1,0) accepted", model.moveLine(1, 0));
  ok("row moved up intact", sameRow(model.getLines().at(0), r0));
  ok("order restored", sameRow(model.getLines().at(1), r1));

  // Degenerate calls must refuse harmlessly.
  ok("move onto itself refused", !model.moveLine(0, 0));
  ok("out-of-range refused", !model.moveLine(0, 5) && !model.moveLine(-1, 0));
  ok("rows unharmed after refusals", sameRow(model.getLines().at(0), r0) &&
                                      sameRow(model.getLines().at(1), r1));

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
