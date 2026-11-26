// The Levels tab's check of the partial waves without a level and its two remedies (PartialWavesDialog.h,
// docs/upgrade/qt_gui.md section 20.3).
#include <algorithm>
#include <cmath>

#include <QDialogButtonBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QRegularExpression>
#include <QVBoxLayout>

#include "AddLevelDialog.h"
#include "LevelsTab.h"
#include "PartialWavesDialog.h"

namespace front = azure2::front;

std::vector<front::MissingPartialWave> LevelsTab::partialWavesWithoutLevel(int *keyMaxL) {
  int key = -1;
  const int maxOrbital = maxOrbitalMomentumSource_ ? maxOrbitalMomentumSource_() : -1;
  std::vector<front::MissingPartialWave> out = front::PartialWavesWithoutLevel(levelRows(), maxL(), maxOrbital, key);
  if (keyMaxL) *keyMaxL = key;
  return out;
}

int LevelsTab::addBlankLevels(const std::vector<front::MissingPartialWave> &waves) {
  if (waves.empty()) return 0;
  // far from everything: 10 MeV above the highest level and the highest threshold, rounded up to whole MeV.  A level
  // with every width 0 adds no resonance at any energy; the energy only has to be a level's (and not another's of its
  // J^pi within 1 keV, which the engine would join).
  double top = 0.0;
  for (const LevelsData &l : levelsModel->getLevels()) top = std::max(top, l.energy);
  if (pairsModel)
    for (const PairsData &p : pairsModel->getPairs()) top = std::max(top, p.seperationEnergy + p.excitationEnergy + p.lightExcitation);
  const double energy = std::ceil(top + 10.0);
  int added = 0;
  undoText_ = tr("blank levels");
  for (const front::MissingPartialWave &w : waves) {
    LevelsData level{1, 1, w.twoJ / 2.0, w.parity, energy, w.nucleusZ, w.nucleusA};
    if (levelsModel->isLevel(level) != -1) continue;
    addLevel(level, false);
    const int row = (int)levelsModel->getLevels().size() - 1;
    const QList<ChannelsData> channels = channelsModel->getChannels();
    for (int i = 0; i < channels.size(); i++)
      if (channels.at(i).levelIndex == row) {
        channelsModel->setData(channelsModel->index(i, 6), 0.0, Qt::EditRole);
        channelsModel->setData(channelsModel->index(i, 0), 1, Qt::EditRole);
      }
    added++;
  }
  setStatus(tr("Added %n blank level(s) at E_x = %1 MeV: every width 0, energy and widths fixed; they carry their "
               "Jπ into the calculation and add no resonance.", "", added).arg(energy));
  updateMissingSpinParities();
  return added;
}

void LevelsTab::showPartialWaves() {
  if (!partialWaves_) partialWaves_ = new PartialWavesDialog(this, this);
  partialWaves_->refresh();
  partialWaves_->show();
  partialWaves_->raise();
}

PartialWavesDialog::PartialWavesDialog(LevelsTab *tab, QWidget *parent) : QDialog(parent), tab_(tab) {
  setWindowTitle(tr("Partial waves without a level"));
  setObjectName("partialWavesDialog");
  const QString jpi = QString("J<sup>%1</sup>").arg(QChar(0x03C0));
  text_ = new QLabel(
      tr("<p>AZURE2 builds its R-matrix one %1 group at a time, and it makes a group only for the %1 of an <b>active "
         "level</b>.  A %1 listed below has none, so the engine makes no channel for it: its partial waves are <b>left "
         "out of the calculation</b> &mdash; no hard-sphere (Coulomb or potential) scattering in them, no background of a "
         "pair's potential, no external capture from them.  Elastic scattering then misses those partial waves without "
         "any message.</p>"
         "<p>Two ways to include them:</p><ul>"
         "<li><b>Hard-sphere scattering in every %1</b> (the &lt;config&gt; key <code>includeAllJPi=1</code>; off when the "
         "file does not set it): for every %1 without a group the engine adds one with every particle channel up to "
         "l = <code>maxOrbitalMomentum</code> (absent: the largest l of a particle channel of the levels) and no level.  "
         "Switched on here, the tab's maximum l is written as <code>maxOrbitalMomentum</code> so that it reaches every "
         "%1 below.</li>"
         "<li><b>A blank level</b> of the %1: every width 0, its energy and widths fixed, far from the data.  It carries "
         "the group into the calculation and adds no resonance.</li></ul>").arg(jpi));
  text_->setWordWrap(true);
  text_->setTextFormat(Qt::RichText);
  state_ = new QLabel;
  state_->setWordWrap(true);
  state_->setTextFormat(Qt::RichText);
  list_ = new QListWidget;
  list_->setObjectName("partialWavesList");
  switchOn_ = new QPushButton(tr("Switch hard-sphere scattering in every Jπ on"));
  switchOn_->setObjectName("partialWavesSwitchOn");
  connect(switchOn_, &QPushButton::clicked, this, [this]() {
    tab_->setIncludeAllJPi(true);
    refresh();
  });
  addBlank_ = new QPushButton(tr("Add blank levels for the ticked"));
  addBlank_->setObjectName("partialWavesAddBlank");
  connect(addBlank_, &QPushButton::clicked, this, [this]() {
    int key = -1;
    const std::vector<front::MissingPartialWave> all = tab_->partialWavesWithoutLevel(&key);
    std::vector<front::MissingPartialWave> ticked;
    for (int i = 0; i < list_->count() && i < (int)all.size(); i++)
      if (list_->item(i)->checkState() == Qt::Checked) ticked.push_back(all[i]);
    tab_->addBlankLevels(ticked);
    refresh();
  });
  all_ = new QPushButton(tr("Tick all"));
  connect(all_, &QPushButton::clicked, this, [this]() {
    for (int i = 0; i < list_->count(); i++) list_->item(i)->setCheckState(Qt::Checked);
  });
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(all_);
  buttons->addStretch(1);
  buttons->addWidget(switchOn_);
  buttons->addWidget(addBlank_);
  QDialogButtonBox *close = new QDialogButtonBox(QDialogButtonBox::Close);
  connect(close, &QDialogButtonBox::rejected, this, &QDialog::hide);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addWidget(text_);
  layout->addWidget(state_);
  layout->addWidget(list_, 1);
  layout->addLayout(buttons);
  layout->addWidget(close);
  resize(720, 560);
}

void PartialWavesDialog::refresh() {
  int key = -1;
  const std::vector<front::MissingPartialWave> missing = tab_->partialWavesWithoutLevel(&key);
  list_->clear();
  for (const front::MissingPartialWave &m : missing) {
    QString text = tr("%1   (lowest l = %2, pair %3)").arg(AddLevelDialog::spinParityLabel(m.twoJ / 2.0, m.parity).remove(QRegularExpression("<[^>]*>")))
                       .arg(m.l)
                       .arg(m.pairKey);
    if (tab_->includeAllJPi()) text += m.key ? tr(" — scatters as a hard sphere (includeAllJPi)") : tr(" — beyond the key's l = %1: left out").arg(key);
    else text += tr(" — left out of the calculation");
    QListWidgetItem *item = new QListWidgetItem(text);
    item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable | Qt::ItemIsSelectable);
    item->setCheckState(Qt::Unchecked);
    list_->addItem(item);
  }
  const QString jpi = QString("J<sup>%1</sup>").arg(QChar(0x03C0));
  if (missing.empty())
    state_->setText(tr("<b>Every %1 the pairs form up to l = %2 has an active level.</b>").arg(jpi).arg(tab_->maxL()));
  else
    state_->setText(tr("<b>%1 %2 without an active level</b> (up to l = %3, the tab's maximum l).  Hard-sphere scattering "
                       "in every %2 is <b>%4</b>%5.")
                        .arg(missing.size())
                        .arg(jpi)
                        .arg(tab_->maxL())
                        .arg(tab_->includeAllJPi() ? tr("on") : tr("off"))
                        .arg(tab_->includeAllJPi() ? tr(", up to l = %1").arg(key) : QString()));
  switchOn_->setEnabled(!tab_->includeAllJPi() && !missing.empty());
  addBlank_->setEnabled(!missing.empty());
  all_->setEnabled(!missing.empty());
}
