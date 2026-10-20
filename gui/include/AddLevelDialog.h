#ifndef ADDLEVELDIALOG_H
#define ADDLEVELDIALOG_H

#include <QDialog>
#include <QtWidgets>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;

QT_END_NAMESPACE

/*!
 * One J^pi the level dialog offers, with the number of channels the current pairs and limits give it.
 */
struct SpinParityOption {
  double jValue;
  int piValue;
  int numChannels;
};

/*!
 * Dialog for adding or editing one level.
 *
 * J^pi is picked from a list rather than typed: the Levels tab passes the values that give the level at least one
 * channel under the current pairs and limits, so a J^pi with the wrong integer/half-integer character cannot be entered.
 */
class AddLevelDialog : public QDialog {
  Q_OBJECT

 public:
  AddLevelDialog(QWidget *parent = 0);
  /*!
   * Fills the J^pi list. An empty list (no pairs yet) falls back to J = 0 to 10 in steps of 1/2, both parities. When
   * haveCurrent is set, the current J^pi is selected, and added (marked as having no channels) if it is not in the list.
   */
  void setSpinParityOptions(const QList<SpinParityOption> &options, double currentJ = 0., int currentPi = 1,
                            bool haveCurrent = false);
  double jValue() const;
  int piValue() const;
  static QString spinParityLabel(double jValue, int piValue);
  /*!
   * A project whose pairs form several compound nuclei: a drop-down picks the level's nucleus and the J^pi list
   * follows it (one option list per nucleus, same order as the labels). Hidden with fewer than two nuclei; on an edit
   * (editable = false) the nucleus is shown but cannot be changed, since the channels would have to be rebuilt.
   */
  void setNucleusOptions(const QStringList &labels, const QList<QList<SpinParityOption> > &optionsPerNucleus,
                         int current = 0, bool editable = true);
  /*! Index into the labels passed to setNucleusOptions (0 when none were). */
  int nucleusIndex() const;

  QLineEdit *energyText;

 private slots:
  void updateOkButton();
  void nucleusChanged(int index);

 private:
  void addOption(double jValue, int piValue, const QString &suffix);

  QLabel *spinParityLabelWidget;
  QLabel *energyLabel;
  QComboBox *spinParityCombo;
  QLabel *nucleusLabelWidget;
  QComboBox *nucleusCombo;
  QList<QList<SpinParityOption> > nucleusOptions;
  QLabel *noteLabel;

  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
