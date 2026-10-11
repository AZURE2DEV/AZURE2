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

QString PairKeyCombo::labelFor(int key) const {
  QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  if (key < 1 || key > pairs.size()) return tr("%1: (no such pair)").arg(key);
  const PairsData &pair = pairs.at(key - 1);
  if (mode == Exit && entranceKey >= 1 && entranceKey <= pairs.size())
    return QString("%1: %2").arg(key).arg(PairsModel::plainText(pairsModel->getReactionLabel(pairs.at(entranceKey - 1), pair)));
  return QString("%1: %2").arg(key).arg(pairsModel->getPairLabel(pair));
}

void PairKeyCombo::refresh() {
  int current = key();
  blockSignals(true);
  clear();
  int n = pairsModel ? pairsModel->getPairs().size() : 0;
  for (int k = 1; k <= n; k++) addItem(labelFor(k), k);
  if (current > n) addItem(labelFor(current), current);
  setCurrentIndex(current >= 1 ? findData(current) : -1);
  blockSignals(false);
}

void PairKeyCombo::showPopup() {
  refresh();
  QComboBox::showPopup();
}

void PairKeyCombo::setKey(int key) {
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
