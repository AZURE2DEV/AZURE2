// A tab's content and the records panels below it (RecordExtras.h).
#include "RecordExtras.h"

#include <QSettings>
#include <QTabWidget>

#include "ExpertMode.h"
#include "RecordTableModel.h"
#include "RecordsPanel.h"

RecordExtras::RecordExtras(QWidget *content, const QString &name, QWidget *parent) : QSplitter(Qt::Vertical, parent) {
  setObjectName(name);
  expert_ = expert::Enabled();
  tabs_ = new QTabWidget;
  tabs_->setObjectName(name + "Tabs");
  addWidget(content);
  addWidget(tabs_);
  setCollapsible(0, false);
  setCollapsible(1, true);
  setStretchFactor(0, 3);
  setStretchFactor(1, 1);
  tabs_->setVisible(false);
  connect(this, &QSplitter::splitterMoved, this, [this]() { QSettings().setValue("layout/" + objectName(), saveState()); });
}

void RecordExtras::addPanel(const QString &title, RecordsPanel *panel) {
  panels_ << panel;
  titles_ << title;
  tabs_->addTab(panel, title);
  connect(panel->model(), &QAbstractItemModel::modelReset, this, &RecordExtras::refresh);
  connect(panel->model(), &QAbstractItemModel::rowsInserted, this, &RecordExtras::refresh);
  connect(panel->model(), &QAbstractItemModel::rowsRemoved, this, &RecordExtras::refresh);
  refresh();
}

void RecordExtras::addWidgetPanel(const QString &title, QWidget *panel, bool first) {
  widgetPanels_ << panel;
  if (first) {
    tabs_->insertTab(0, panel, title);
    tabs_->setCurrentIndex(0);
  } else {
    tabs_->addTab(panel, title);
  }
  refresh();
}

void RecordExtras::setExpertMode(bool on) {
  expert_ = on;
  for (RecordsPanel *p : panels_) p->setExpertMode(on);
  refresh();
}

bool RecordExtras::panelShown(RecordsPanel *panel) const {
  const int i = tabs_->indexOf(panel);
  return i >= 0 && tabs_->isTabVisible(i) && tabs_->isVisible();
}

void RecordExtras::refresh() {
  bool any = !widgetPanels_.isEmpty();  // a panel without records is no expert option: always shown
  for (int i = 0; i < panels_.size(); i++) {
    const int n = panels_[i]->recordCount();
    const int at = tabs_->indexOf(panels_[i]);
    const bool shown = expert_ || n > 0;  // a panel that holds records is shown in either view
    tabs_->setTabVisible(at, shown);
    tabs_->setTabText(at, n > 0 ? QString("%1 (%2)").arg(titles_[i]).arg(n) : titles_[i]);
    tabs_->setTabToolTip(at, !expert_ && n > 0 ? expert::ShownBecauseSet() : QString());
    any = any || shown;
  }
  tabs_->setVisible(any);
}

void RecordExtras::showEvent(QShowEvent *event) {
  QSplitter::showEvent(event);
  if (sized_ || height() < 200) return;
  sized_ = true;
  const QByteArray saved = QSettings().value("layout/" + objectName()).toByteArray();
  if (!saved.isEmpty() && restoreState(saved)) return;
  // the content first: the panels below get a quarter of the height (an empty one less), dragged larger when needed
  // (a panel of its own, the Potential panel, needs the room: 5 parts of 7, as before the splitter)
  const int lower = widgetPanels_.isEmpty() ? std::min(height() / 4, 220) : height() * 5 / 7;
  setSizes({height() - lower, lower});
}
