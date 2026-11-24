#ifndef VALUEEDITOR_H
#define VALUEEDITOR_H

#include <QList>
#include <QString>
#include <QStringList>
#include <QWidget>
#include <functional>

#include "azure2/front/Engine.h"

class ProjectHolder;

QT_BEGIN_NAMESPACE
class QCheckBox;
class QLabel;
class QLineEdit;
class QStackedWidget;
class QToolButton;
QT_END_NAMESPACE

/*!
 * The editors of a record's values drawn from the schema (docs/upgrade/qt_gui.md section 11): every choice a widget
 * offers comes from the facade (ValueInfo: units, enum words, alternatives, grid kinds, the kinds a reference names;
 * the project's records and pairs), never from a word spelled here.
 *
 * A PieceEditor edits one value of a ValueInfo as its pieces (azure2::front::ValueParts): a nuisance as a number, a
 * unit, "Fit" and a prior; a structured value as its alternative and parts; a grid, the columns, a reaction ...  A
 * FieldEditor is one field of a record: a PieceEditor, an "Edit as text" toggle, the parser's message under it, and
 * the round-trip rule: until a widget is changed, text() is the text it was given, byte for byte.
 */

/// What an editor needs of the record it edits.
struct EditorContext {
  ProjectHolder *holder = nullptr;
  std::string kind;                     ///< the record's kind
  std::function<QString()> quantity;    ///< the record's quantity as the form holds it ("" none)
  /// The value an optional field or part takes when it is left out ("" none): shown greyed as the empty choice
  /// ("yes (default)") or as a line's placeholder, never written.
  std::string defaultText;
};

class PieceEditor : public QWidget {
  Q_OBJECT
 public:
  explicit PieceEditor(QWidget *parent = nullptr) : QWidget(parent) {}
  /// Show the pieces; false when this editor cannot show them (the field then opens as text).
  virtual bool set(const azure2::front::ValueParts &parts) = 0;
  virtual azure2::front::ValueParts get() const = 0;
  /// Read the project again (record names, pairs, the quantity's units) keeping the choice.
  virtual void refresh() {}
  /// Choose sensible first values (a new record): the first word, unit, record or reaction; a dimensionless
  /// nuisance 1, fixed.
  virtual void prefill() {}
 signals:
  void changed();
};

/// The editor of a value of shape `info` (a field's value or a part of a structured one).
PieceEditor *MakePieceEditor(const azure2::front::ValueInfo &info, const EditorContext &context, bool optional,
                             QWidget *parent = nullptr);

class FieldEditor : public QWidget {
  Q_OBJECT
 public:
  FieldEditor(const azure2::front::FieldInfo &field, const EditorContext &context, QWidget *parent = nullptr);

  const azure2::front::FieldInfo &field() const { return field_; }
  /// Show a value as written ("" for a field the record does not give).  Not an edit.
  void setText(const QString &text);
  /// The value to write: the text given until a widget changed, else the text of the pieces; "" when the field is not
  /// to be written (an optional field left empty or switched off).
  QString text() const;
  bool isEdited() const { return edited_; }
  /// The field is given (an optional structured field has a "use" switch).
  bool isUsed() const;
  bool textMode() const { return textMode_; }
  void setTextMode(bool on);
  /// Type `text` into the field's text line, as the user would (an edit; the field shows as text).
  void typeText(const QString &text);
  /// The message under the field ("" none): the parser's, or the validator's after Apply.
  void showMessage(const QString &message, bool error = true);
  QString message() const;
  /// Check the text with the engine's parser now (the editor does it on every change).
  bool check();
  PieceEditor *piece() const { return piece_; }
  /// Read the project's records, pairs and the quantity's units again.
  void refresh();
  /// Take the widgets' first choices as the value (a new record): an edit.
  void prefill();
  /// Give or leave out an optional structured field (its "use" switch); an edit.
  void setUsed(bool on);
  /// Whether the field applies to the record as the form holds it (its schema condition, "kind is layer or ..."): a
  /// required field that does not apply is not marked missing (RecordsPanel::applyConditions).
  void setApplies(bool on);
  bool applies() const { return applies_; }

 signals:
  void changed();

 private:
  void pieceChanged();
  void useToggled(bool on);

  azure2::front::FieldInfo field_;
  EditorContext context_;
  PieceEditor *piece_ = nullptr;
  QLineEdit *textEdit_;
  QToolButton *textToggle_;
  QCheckBox *use_ = nullptr;
  QLabel *message_;
  QString original_;
  bool edited_ = false;
  bool textMode_ = false;
  bool loading_ = false;
  bool applies_ = true;
};

/// Is `field` the record's quantity (an enum whose words are the quantities of its kind)?
bool IsQuantityField(const azure2::front::FieldInfo &field, const std::string &kind);
/// The names of the records of `kinds` in the project (each kind read in the sections the schema puts it).
QStringList RecordNamesOfKinds(const ProjectHolder *holder, const std::vector<std::string> &kinds);
/// A value's grammar as a hint a user reads (qt_gui.md section 19.9): "<A><Symbol>" as "e.g. 20Ne", a single
/// "<word>" as the word, anything with alternatives or optional parts ("fixed|fit [prior ...]") as "" (the grammar
/// stays in the expert view's tips and the text mode).
QString PlainHint(const std::string &grammar);
/// The text with each number written with more than 12 significant digits shown with 12 (shortest form: "5.5" for
/// "5.4999999999999997", "0.3" for "0.30000000000000004"); other tokens as they are.  Display only (qt_gui.md 19.9).
QString RoundedNumbers(const QString &text);
/// The value with one blank between tokens outside quotes.
QString OneBlank(const QString &text);

#endif
