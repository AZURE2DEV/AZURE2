#include "ThmChannelsPage.h"

#include <QCheckBox>
#include <QDoubleValidator>
#include <QGridLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QTableWidget>
#include <QVBoxLayout>
#include <cmath>
#include <sstream>

#include "ChannelsModel.h"
#include "LevelsModel.h"
#include "PairsModel.h"
#include "RichTextDelegate.h"
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"

namespace {

// Shortest text that reads back as the same double.
QString numberText(double x) { return QString::number(x, 'g', QLocale::FloatingPointShortest); }

// A whole token as a number, read as the engine reads a <levels> column.
bool readWholeDouble(const QString &text, double &x) {
  std::istringstream s(text.trimmed().toStdString());
  std::string rest;
  return !!(s >> x) && !(s >> rest) && std::isfinite(x);
}

const char *kDocs =
    "docs/source/theory/thm_implementation.rst; online: "
    "<a href=\"https://rdeboer1.github.io/AZURE2/\">rdeboer1.github.io/AZURE2</a> "
    "(Theory &gt; Trojan Horse (HOES) Observable).";

}  // namespace

QSet<int> ThmChannelsPage::thmEntrancePairs(SegmentsDataModel *segmentsData, SegmentsTestModel *segmentsTest) {
  QSet<int> keys;
  for (const SegmentsDataData &s : segmentsData->getLines())
    if (s.isTHM) keys.insert(s.entrancePairIndex);
  for (const SegmentsTestData &s : segmentsTest->getLines())
    if (s.isTHM) keys.insert(s.entrancePairIndex);
  return keys;
}

ThmChannelsPage::ThmChannelsPage(PairsModel *pairs, LevelsModel *levels, ChannelsModel *channels,
                                 SegmentsDataModel *segmentsData, SegmentsTestModel *segmentsTest, QWidget *parent) :
  QWidget(parent),
  pairs_(pairs),
  channels_(channels) {
  const QSet<int> thm = thmEntrancePairs(segmentsData, segmentsTest);
  const QList<PairsData> pairList = pairs->getPairs();
  const QList<LevelsData> levelList = levels->getLevels();
  const QList<ChannelsData> channelList = channels->getChannels();

  // Which segments use each pair as their THM entrance.
  QMap<int, QStringList> users;
  const QList<SegmentsDataData> data = segmentsData->getLines();
  for (int k = 0; k < data.size(); k++)
    if (data.at(k).isTHM) users[data.at(k).entrancePairIndex] << tr("data %1").arg(k + 1);
  const QList<SegmentsTestData> test = segmentsTest->getLines();
  for (int k = 0; k < test.size(); k++)
    if (test.at(k).isTHM) users[test.at(k).entrancePairIndex] << tr("test %1").arg(k + 1);

  pairTable = new QTableWidget(0, 4);
  pairTable->setHorizontalHeaderLabels(QStringList() << tr("Pair") << tr("Nuclei") << tr("THM segments")
                                                     << tr("B (MeV)"));
  pairTable->horizontalHeader()->setStretchLastSection(true);
  pairTable->verticalHeader()->hide();
  pairTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  RichTextDelegate *rich = new RichTextDelegate;
  rich->setParent(pairTable);
  pairTable->setItemDelegateForColumn(1, rich);
  pairTable->setToolTip(
      tr("B: binding energy (MeV) of the transferred particle x in the Trojan horse a = x + s; the "
         "half-off-shell momentum is p = sqrt(2 mu (E + B)). Field 32 of the pair's <levels> lines."));
  for (int key = 1; key <= pairList.size(); key++) {
    const PairsData &p = pairList.at(key - 1);
    if (p.pairType != 0 || !(thm.contains(key) || p.bindingEnergy != 0.0)) continue;
    int row = pairTable->rowCount();
    pairTable->insertRow(row);
    pairRows_ << key;
    pairTable->setItem(row, 0, new QTableWidgetItem(QString::number(key)));
    QString nuclei = pairs->getParticleLabel(p, 0) + " + " + pairs->getParticleLabel(p, 1);
    nuclei.remove("<center>").remove("</center>");  // one line: the delegate renders the rest
    pairTable->setItem(row, 1, new QTableWidgetItem(nuclei));
    pairTable->setItem(row, 2, new QTableWidgetItem(users.contains(key) ? users[key].join(", ")
                                                                        : tr("none (B kept as in the file)")));
    QLineEdit *edit = new QLineEdit(numberText(p.bindingEnergy));
    QDoubleValidator *v = new QDoubleValidator(edit);
    v->setLocale(QLocale::c());
    edit->setValidator(v);
    edit->setProperty("shown", edit->text());
    pairTable->setCellWidget(row, 3, edit);
  }
  pairTable->resizeColumnsToContents();
  pairTable->horizontalHeader()->setStretchLastSection(true);

  channelTable = new QTableWidget(0, 6);
  channelTable->setHorizontalHeaderLabels(QStringList() << tr("Level") << tr("Pair") << tr("l") << tr("s")
                                                        << tr("Width") << tr("Amplitude"));
  channelTable->horizontalHeaderItem(4)->setToolTip(tr("The width column of the Levels tab, as entered"));
  channelTable->horizontalHeaderItem(5)->setToolTip(tr("Entered as a reduced width amplitude, MeV^(1/2)"));
  channelTable->horizontalHeader()->setStretchLastSection(true);
  channelTable->verticalHeader()->hide();
  channelTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  channelTable->setToolTip(
      tr("Field 33: ticked, the width column of this channel is a reduced width amplitude gamma (MeV^(1/2)) "
         "instead of a partial width (eV) or ANC, e.g. for a level known only from THM data. The value is not "
         "converted when the flag changes: re-enter it in the Levels tab."));
  for (int c = 0; c < channelList.size(); c++) {
    const ChannelsData &ch = channelList.at(c);
    if (ch.radType != 'P') continue;
    if (!(thm.contains(ch.pairIndex + 1) || ch.gammaIsRWA == 1)) continue;
    if (ch.levelIndex < 0 || ch.levelIndex >= levelList.size()) continue;
    const LevelsData &level = levelList.at(ch.levelIndex);
    int row = channelTable->rowCount();
    channelTable->insertRow(row);
    channelRows_ << c;
    channelTable->setItem(row, 0, new QTableWidgetItem(levels->getSpinLabel(level) + ", " +
                                                        numberText(level.energy) + " MeV"));
    channelTable->setItem(row, 1, new QTableWidgetItem(QString::number(ch.pairIndex + 1)));
    channelTable->setItem(row, 2, new QTableWidgetItem(QString::number(ch.lValue)));
    channelTable->setItem(row, 3, new QTableWidgetItem(channels->getSpinLabel(ch)));
    channelTable->setItem(row, 4, new QTableWidgetItem(numberText(ch.reducedWidth)));
    QTableWidgetItem *flag = new QTableWidgetItem;
    flag->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    flag->setCheckState(ch.gammaIsRWA == 1 ? Qt::Checked : Qt::Unchecked);
    flag->setData(Qt::UserRole, flag->checkState() == Qt::Checked);  // as read
    channelTable->setItem(row, 5, flag);
  }
  channelTable->resizeColumnsToContents();
  channelTable->horizontalHeader()->setStretchLastSection(true);

  QGroupBox *pairBox = new QGroupBox(tr("Binding energy of the transferred particle"));
  QVBoxLayout *pl = new QVBoxLayout;
  pl->addWidget(pairTable);
  pairBox->setLayout(pl);
  QGroupBox *channelBox = new QGroupBox(tr("Width input convention"));
  QVBoxLayout *cl = new QVBoxLayout;
  cl->addWidget(channelTable);
  channelBox->setLayout(cl);

  QLabel *docs = new QLabel(tr("Pairs that are the entrance of a THM segment, and their particle channels. "
                               "Documentation: ") +
                            QString(kDocs));
  docs->setWordWrap(true);
  docs->setOpenExternalLinks(true);
  docs->setTextFormat(Qt::RichText);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(pairBox, 1);
  mainLayout->addWidget(channelBox, 3);
  mainLayout->addWidget(docs);
  setLayout(mainLayout);
}

void ThmChannelsPage::setBindingText(int pairKey, const QString &text) {
  int row = pairRow(pairKey);
  if (row >= 0) qobject_cast<QLineEdit *>(pairTable->cellWidget(row, 3))->setText(text);
}

void ThmChannelsPage::setReducedWidthFlag(int channelIndex, bool rwa) {
  int row = channelRow(channelIndex);
  if (row >= 0) channelTable->item(row, 5)->setCheckState(rwa ? Qt::Checked : Qt::Unchecked);
}

QString ThmChannelsPage::check() const {
  for (int row = 0; row < pairRows_.size(); row++) {
    const QString text = qobject_cast<QLineEdit *>(pairTable->cellWidget(row, 3))->text();
    double x;
    if (!readWholeDouble(text, x))
      return tr("Pair %1: the binding energy '%2' is not a number (MeV).").arg(pairRows_.at(row)).arg(text);
  }
  return QString();
}

void ThmChannelsPage::apply() {
  for (int row = 0; row < pairRows_.size(); row++) {
    QLineEdit *edit = qobject_cast<QLineEdit *>(pairTable->cellWidget(row, 3));
    if (edit->text().trimmed() == edit->property("shown").toString()) continue;  // exactly as read
    double x;
    if (!readWholeDouble(edit->text(), x)) continue;
    pairs_->setData(pairs_->index(pairRows_.at(row) - 1, 15), x, Qt::EditRole);
  }
  const QList<ChannelsData> channelList = channels_->getChannels();
  for (int row = 0; row < channelRows_.size(); row++) {
    const int c = channelRows_.at(row);
    const QTableWidgetItem *item = channelTable->item(row, 5);
    const bool rwa = item->checkState() == Qt::Checked;
    if (rwa == item->data(Qt::UserRole).toBool() || c >= channelList.size()) continue;  // as read
    channels_->setData(channels_->index(c, 7), rwa ? 1 : 0, Qt::EditRole);
  }
}
