#ifndef MIXTUREEDITOR_H
#define MIXTUREEDITOR_H

#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QLineEdit>
#include <QTableWidget>
#include <QPushButton>
#include <QLabel>
#include <QList>
#include <QString>
#include <QStringList>

/*!
 * Editor for a segment line's trailing "mixture" block -- a natural-abundance
 * or compound target whose other nuclei are frozen as cross-section tables:
 *
 *   mixture <w_self> <n> {"<file>" <w_i>}...
 *
 * The block rides on the row as verbatim text (SegmentsDataData::extraTokens)
 * beside any other keyword block ("sqrtshift ..."), which this editor leaves
 * untouched: setTokens() splits the mixture block out, tokens() puts the row's
 * text back together with the mixture block last, as the engine's writers do.
 * A block it cannot parse is kept verbatim and shown, never dropped.
 */
class MixtureEditor : public QWidget {
  Q_OBJECT

 public:
  struct Table {
    QString file;
    double weight;
  };

  explicit MixtureEditor(QWidget *parent = 0);

  /*! Load the row's keyword text; the mixture block fills the editor, the rest is kept aside. */
  void setTokens(const QString &extraTokens);
  /*! The row's keyword text with the editor's mixture block (if enabled) appended last. */
  QString tokens() const;
  /*! Is a mixture block being edited (the check box)? */
  bool isMixtureEnabled() const;
  /*! Every field parses and every table has a file; the reason otherwise. */
  bool validate(QString *error = 0) const;
  /*! The editor's current block, independent of the row's other text ("" when disabled). */
  QString block() const;

  /*! Split keyword text into tokens, a double-quoted token (a path with spaces) kept whole. */
  static QStringList splitTokens(const QString &text);
  /*! Find and decode the mixture block; \a other receives the remaining text. */
  static bool parseBlock(const QString &extraTokens, bool &present, double &selfWeight,
                         QList<Table> &tables, QString &other, QString *error = 0);
  /*! Encode a block in the engine's syntax, every file quoted. */
  static QString formatBlock(double selfWeight, const QList<Table> &tables);

  QCheckBox *enableCheck;
  QLineEdit *selfWeightText;
  QTableWidget *tablesTable;
  QPushButton *addTableButton;
  QPushButton *removeTableButton;
  QLabel *sumLabel;
  QLabel *rawLabel;

 public slots:
  void enableChanged(int state);
  void addTable();
  void removeTable();
  void updateSum();

 private:
  QGroupBox *box_;
  QString other_;   // the row's keyword text without the mixture block
  QString raw_;     // a mixture block that did not parse, kept verbatim
  QList<Table> currentTables(bool *ok = 0) const;
};

#endif
