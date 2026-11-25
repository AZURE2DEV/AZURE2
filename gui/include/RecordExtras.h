#ifndef RECORDEXTRAS_H
#define RECORDEXTRAS_H

#include <QList>
#include <QSplitter>

class RecordsPanel;

QT_BEGIN_NAMESPACE
class QTabWidget;
QT_END_NAMESPACE

/*!
 * A tab's main content above and the records panels that belong to it below (the Pairs tab's per-channel radii, the
 * Levels tab's gamma decays and isospin): docs/upgrade/qt_gui.md sections 19.3 and 19.4.  The lower part is a page of
 * sub-tabs in a splitter, small at first (the sizes a user sets are kept), its tabs titled with their record counts.  In
 * standard view a panel without records is hidden (an expert option), one with records is shown; the lower part is
 * hidden when no panel is shown.
 */
class RecordExtras : public QSplitter {
  Q_OBJECT

 public:
  /// `content` is the tab's own widget; `name` keys the saved sizes.
  RecordExtras(QWidget *content, const QString &name, QWidget *parent = nullptr);
  void addPanel(const QString &title, RecordsPanel *panel);
  /// A panel that holds no records (the Pairs tab's Potential panel, qt_gui.md section 20): shown in either view, and
  /// the lower part then gets most of the height at first; `first`: as the first and current sub-tab.
  void addWidgetPanel(const QString &title, QWidget *panel, bool first);
  void setExpertMode(bool on);
  QTabWidget *tabs() const { return tabs_; }
  /// Whether the panel's sub-tab is shown (the tests).
  bool panelShown(RecordsPanel *panel) const;

 protected:
  void showEvent(QShowEvent *event) override;

 private:
  void refresh();
  QTabWidget *tabs_ = nullptr;
  QList<RecordsPanel *> panels_;
  QStringList titles_;
  QList<QWidget *> widgetPanels_;
  bool expert_ = false;
  bool sized_ = false;
};

#endif
