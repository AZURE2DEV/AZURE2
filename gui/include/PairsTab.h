#ifndef PAIRSTAB_H
#define PAIRSTAB_H

#include <QWidget>
#include <QItemSelection>
#include <QSignalMapper>
#include <QTableView>
#include <QPushButton>
#include <QSignalMapper>
#include <QPointer>

#include "PairsModel.h"
#include <QGroupBox>
#include <QLabel>
#include <QShowEvent>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QTextStream>
#include "AddPairDialog.h"

QT_BEGIN_NAMESPACE

class QPushButton;
class QTabWidget;
class QUndoStack;

QT_END_NAMESPACE

class InfoDialog;
class RecordsPanel;
class RecordExtras;

/*!
 * The Particle Pairs tab: the pairs of the compound nucleus, and the source every other tab takes its pair list from.
 */
class PairsTab : public QWidget {
  Q_OBJECT

 public:
  PairsTab(QWidget *parent = 0);
  PairsModel *getPairsModel();
  /// A panel shown below the pairs (the per-channel radii of <channelRadii>, AZURESetup).
  void addExtraPanel(const QString &title, RecordsPanel *panel);
  /// Standard or expert view (qt_gui.md section 19.4): the per-channel radii hidden in standard view unless they hold records.
  void setExpertMode(bool on);
  RecordExtras *extraPanels() const { return extras; }
  /*! A new file: no pairs. */
  void reset();
  /// Select pair `row` and remove it as the Remove button does (the tests).
  void removePairAt(int row);
  /// The undo stack shared with the Levels tab (LevelsTab::setPairsModel): Ctrl+Z, Ctrl+Shift+Z here too.
  void setUndoStack(QUndoStack *stack);
  /// A number of pair `row` edited in place (column of PairsModel: 10, 11, 12, 15): false with the reason.
  bool setPairNumber(int row, int column, const QString &text, QString &why);
  /// What goes with pair `row` (its levels' channels, the data sets and extrapolations that name it), in words.
  QStringList pairConsequences(int row) const;
  QTableView *getPairsView() const { return pairsView; }
  /// Ask before a removal (false: the tests' removePairAt removes at once).
  void setConfirmRemoval(bool on) { confirmRemoval_ = on; }

 public slots:
  void addPair();
  /// The selected pair's copy in the pair dialog, added as a new pair when it differs (Ctrl+D).
  void duplicatePair();
  void addPair(PairsData pair, int pairIndex, bool fromFile);
  void editPair();
  void editPair(PairsData pair, int pairIndex, bool fromFile);
  void removePair();
  void updateButtons(const QItemSelection &selection);
  void showInfo(int which = 0, QString title = "");
  void updateNucleusColumn();
  /*! "Masses from AME2020...": compare every pair with the mass table and fill the ticked ones (MassCheckDialog). */
  void checkMasses();

 signals:
  void pairAdded(int);
  void pairRemoved(int);
  void pairEdited(int);

 private:
  PairsModel *pairsModel;
  QTableView *pairsView;
  RecordExtras *extras = nullptr;
  bool confirmRemoval_ = true;

 private:
  QPushButton *addButton;
  QPushButton *deleteButton;
  QPushButton *duplicateButton = nullptr;
  QPushButton *massButton;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
