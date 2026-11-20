#ifndef PARTIALWAVESDIALOG_H
#define PARTIALWAVESDIALOG_H

#include <QDialog>

class LevelsTab;

QT_BEGIN_NAMESPACE
class QLabel;
class QListWidget;
class QPushButton;
QT_END_NAMESPACE

/*!
 * The Levels tab's check of the partial waves (docs/upgrade/qt_gui.md section 19.3): every J^pi the particle pairs form
 * up to the tab's maximum l that no active level carries, with what that means (the engine makes no J^pi group for
 * it: its partial waves are left out of the calculation unless includeAllJPi is on) and two ways out: switch
 * includeAllJPi on (the tab's "Hard-sphere scattering in every J^pi"), or add a blank level -- zero widths, energy and
 * widths fixed -- for each J^pi ticked.  The list comes from the facade (front::PartialWavesWithoutLevel, the engine's
 * rule); the levels are added through the tab (one undo step).
 */
class PartialWavesDialog : public QDialog {
  Q_OBJECT

 public:
  explicit PartialWavesDialog(LevelsTab *tab, QWidget *parent = nullptr);
  QListWidget *list() const { return list_; }
  QLabel *explanation() const { return text_; }
  QPushButton *switchOnButton() const { return switchOn_; }
  QPushButton *addBlankButton() const { return addBlank_; }

 public slots:
  /// The list again (after an edit of the tab).
  void refresh();

 private:
  LevelsTab *tab_;
  QLabel *text_, *state_;
  QListWidget *list_;
  QPushButton *switchOn_, *addBlank_, *all_;
};

#endif
