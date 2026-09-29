#ifndef THMMODELPAGE_H
#define THMMODELPAGE_H

#include <QString>
#include <QStringList>
#include <QWidget>

#include "ThmSettings.h"

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLineEdit;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * The Model page of the THM workspace: the options of the <thm> block, i.e.
 * the choices in the Trojan Horse (HOES) observable of THM segments
 * (isDiff >= 10).  See docs/source/theory/thm_implementation.rst.
 */
class ThmModelPage : public QWidget {
  Q_OBJECT

 public:
  ThmModelPage(const ThmSettings &settings, const QString &projectDir, QWidget *parent = 0);

  /// The settings as the widgets hold them (experiment lines as given).
  ThmSettings settings() const;
  /// Put settings into the widgets.
  void setSettings(const ThmSettings &settings);
  /// What would be refused (duplicate rows, ThmSettings::validate); "" if none.
  QString check() const;

  /// Table rows, for the tests and the Add buttons.
  void addSpectatorRow(int pair, double energy);
  void addWeightRow(bool test, int segment, const QString &file);
  void clearRows();

  QComboBox *entranceLCombo;
  QComboBox *vertexCombo;
  QComboBox *kinematicsCombo;
  QCheckBox *coulombIntegralCheck;
  QDoubleSpinBox *spectatorEnergySpin;
  QTableWidget *spectatorTable;
  QTableWidget *weightTable;

 private slots:
  void chooseWeightFile(QLineEdit *lineEdit);
  void removeSelectedSpectatorRows();
  void removeSelectedWeightRows();

 private:
  QString projectDir_;
  QStringList experimentLines_;  // carried through unchanged (ThmSettings::experimentLines)
};

#endif
