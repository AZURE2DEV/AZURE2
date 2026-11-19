#include "PairKeyCombo.h"
#include "PairsModel.h"

PairKeyCombo::PairKeyCombo(Mode mode, QWidget *parent) : QComboBox(parent), mode(mode) {
  setSizeAdjustPolicy(QComboBox::AdjustToContents);
#if QT_VERSION >= QT_VERSION_CHECK(5, 15, 0)
  setPlaceholderText(tr("(choose a pair)"));
#endif
  connect(this, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) { emit keyChanged(key()); });
}

void PairKeyCombo::setPairsModel(PairsModel *model) {
  pairsModel = model;
  refresh();
}

void PairKeyCombo::setEntranceKey(int key) {
  if (entranceKey == key) return;
  entranceKey = key;
  if (mode == Exit) refresh();
}

bool PairKeyCombo::isGammaPair(int key) const {
  if (!pairsModel || key < 1 || key > pairsModel->getPairs().size()) return false;
  return pairsModel->getPairs().at(key - 1).pairType == 10;
}

bool PairKeyCombo::hasGammaPair() const {
  if (!pairsModel) return false;
  for (const PairsData &pair : pairsModel->getPairs())
    if (pair.pairType == 10) return true;
  return false;
}

void PairKeyCombo::setTotalCaptureOffered(bool offered) {
  if (totalCaptureOffered == offered) return;
  totalCaptureOffered = offered;
  refresh();
}

QString PairKeyCombo::labelFor(int key) const {
  if (key == TotalCapture) return tr("Total capture (all γ pairs)");
  QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  if (key < 1 || key > pairs.size()) return tr("%1: (no such pair)").arg(key);
  const PairsData &pair = pairs.at(key - 1);
  // With several compound nuclei in the project, say which one the pair forms.
  QString nucleus;
  if (pairsModel->nuclei().size() > 1) nucleus = QString("  -> %1").arg(PairsModel::nucleusLabel(PairsModel::nucleusOf(pair), false));
  if (mode == Exit && entranceKey >= 1 && entranceKey <= pairs.size())
    return QString("%1: %2%3").arg(key).arg(PairsModel::plainText(pairsModel->getReactionLabel(pairs.at(entranceKey - 1), pair))).arg(nucleus);
  return QString("%1: %2%3").arg(key).arg(pairsModel->getPairLabel(pair)).arg(nucleus);
}

void PairKeyCombo::refresh() {
  int current = key();
  blockSignals(true);
  clear();
  int n = pairsModel ? pairsModel->getPairs().size() : 0;
  for (int k = 1; k <= n; k++) addItem(labelFor(k), k);
  if (current > n) addItem(labelFor(current), current);
  const bool withTotal = mode == Exit && totalCaptureOffered && (hasGammaPair() || current == TotalCapture);
  if (withTotal) addItem(labelFor(TotalCapture), TotalCapture);
  setCurrentIndex(current >= 1 || (withTotal && current == TotalCapture) ? findData(current) : -1);
  blockSignals(false);
}

void PairKeyCombo::showPopup() {
  refresh();
  QComboBox::showPopup();
}

void PairKeyCombo::setKey(int key) {
  if (key == TotalCapture && mode == Exit && totalCaptureOffered) {
    if (findData(TotalCapture) == -1) addItem(labelFor(TotalCapture), TotalCapture);
    setCurrentIndex(findData(TotalCapture));
    return;
  }
  if (key < 1) {
    setCurrentIndex(-1);
    return;
  }
  if (findData(key) == -1) {
    if (pairsModel && key <= pairsModel->getPairs().size())
      refresh();
    else
      addItem(labelFor(key), key);
  }
  setCurrentIndex(findData(key));
}

int PairKeyCombo::key() const { return currentIndex() >= 0 ? currentData().toInt() : 0; }
