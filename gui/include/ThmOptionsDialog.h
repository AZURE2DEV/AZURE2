#ifndef THMOPTIONSDIALOG_H
#define THMOPTIONSDIALOG_H

#include <QDialog>
#include <QMap>
#include <QString>
#include <QStringList>

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLineEdit;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * The options of the optional <thm> block (Config::ThmOptions), as the GUI
 * edits them.  Values are canonical: vertex=real is read as perlevel, flags as
 * bool.  Defaults are those of the engine, so a default-constructed object is
 * a project without the block.
 */
struct ThmSettings {
  QString entranceL = "incoherent";  // incoherent | coherent
  QString vertex = "constant";       // constant | perlevel | onshell
  QString kinematics = "lacognata";  // lacognata | triple | kf3body | lambda32
  bool coulombIntegral = false;
  double spectatorEnergy = 0.0;       // MeV, every THM entrance pair
  QMap<int, double> spectatorByPair;  // spectatorEnergy[<pair key>]
  QMap<int, QString> weight;          // weight[<segmentsData line>]
  QMap<int, QString> weightTest;      // weightTest[<segmentsTest line>]

  bool operator==(const ThmSettings &o) const;
  bool operator!=(const ThmSettings &o) const { return !(*this == o); }
  bool isDefault() const { return *this == ThmSettings(); }

  /// The non-default keys, canonical key -> canonical value, in the order
  /// they are written.
  QList<QPair<QString, QString>> keyValues() const;

  /*! Parses the lines between <thm> and </thm> with the rules of
      Config::ReadThmBlock ('#' comments, key=value, same keys and values).
      Returns false and a message naming the line on anything the engine
      would refuse.  Weight files are not opened here (see validate). */
  static bool parse(const QStringList &lines, ThmSettings &out, QString *error = nullptr);

  /*! Parses one line; returns false on a line the engine refuses.  A blank or
      comment-only line gives an empty key. */
  static bool parseLine(const QString &line, QString &key, QString &value, ThmSettings &into);

  /*! What the engine will refuse at startup that can be seen without the
      data: every weight file must be readable, two columns, strictly
      increasing E and w > 0 (ThmWeightTable::Read, the engine's own reader).
      A relative path is taken from projectDir, as the engine takes it from
      the directory of the .azr.  Returns an empty string when all is well. */
  QString validate(const QString &projectDir) const;

  /*! The block lines for these settings, given the lines of the block that
      was read (possibly none).  Comment and blank lines are kept in place; a
      key line whose value did not change is kept verbatim (inline comment
      included), a changed one is rewritten keeping its inline comment, a key
      that went back to its default is dropped, and new keys are appended. */
  QStringList compose(const QStringList &oldLines) const;
};

/*!
 * Editor of the <thm> block (Configure > THM Options...): the choices in the
 * Trojan Horse (HOES) observable of THM segments (isDiff >= 10).  See
 * docs/source/theory/thm_implementation.rst.
 */
class ThmOptionsDialog : public QDialog {
  Q_OBJECT

 public:
  ThmOptionsDialog(const ThmSettings &settings, const QString &projectDir, QWidget *parent = 0);

  /// The settings as the widgets hold them.
  ThmSettings settings() const;
  /// Put settings into the widgets.
  void setSettings(const ThmSettings &settings);

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

 public slots:
  void accept() override;

 private slots:
  void chooseWeightFile(QLineEdit *lineEdit);
  void removeSelectedSpectatorRows();
  void removeSelectedWeightRows();

 private:
  QString projectDir_;
};

#endif
