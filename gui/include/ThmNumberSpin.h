#ifndef THMNUMBERSPIN_H
#define THMNUMBERSPIN_H

#include <QDoubleSpinBox>
#include <QLocale>
#include <QRegularExpression>
#include <QValidator>
#include <cmath>
#include <limits>
#include <sstream>

/*!
 * A numeric field of the THM workspace: a spin box with its unit as a suffix,
 * right-aligned, that shows the shortest text of its value and remembers the
 * text it was loaded with.  writtenText() gives that text back as long as the
 * value has not been changed, so a number written by hand ("0.0", "40.00")
 * reads back verbatim; a changed value, or one loaded from a text that is no
 * number in range, gives the shortest text of the value shown.  With a
 * special value text, the minimum stands for "not given" and gives "".  What
 * is typed is read as the engine reads a number (C locale, exponents allowed:
 * "1e-05", the way the field shows small values), and the value is kept to
 * the last bit, not rounded to a number of decimals.
 * No Q_OBJECT: it adds no signal or slot.
 */
class ThmNumberSpin : public QDoubleSpinBox {
 public:
  explicit ThmNumberSpin(const QString &suffix = QString(), double lo = 0.0, double hi = 1.0e6, double step = 1.0,
                         QWidget *parent = nullptr) :
    QDoubleSpinBox(parent) {
    setLocale(QLocale::c());
    // QDoubleSpinBox rounds every value to this many decimals: the most it
    // allows, so that the value keeps what is typed (1e-12 stays 1e-12); the
    // text is the shortest (textFromValue).
    setDecimals(323);
    setRange(lo, hi);
    setSingleStep(step);
    if (!suffix.isEmpty()) setSuffix(suffix);
    setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    setAccelerated(true);
  }

  /// Shows `text` (as read from the file); unreadable or empty: the minimum.
  void setWrittenText(const QString &text) {
    double x = minimum();
    std::istringstream s(text.trimmed().toStdString());
    std::string rest;
    written_ = text.trimmed();
    readable_ = !written_.isEmpty() && (s >> x) && !(s >> rest) && std::isfinite(x) && x >= minimum() && x <= maximum();
    setValue(readable_ ? x : minimum());
    writtenValue_ = value();
  }
  /*! The text to write: as loaded while the value is unchanged (and the text
      is the number shown), else the shortest text of the value; "" for the
      minimum under a special value text (not given). */
  QString writtenText() const {
    if (readable_ && value() == writtenValue_) return written_;
    if (!specialValueText().isEmpty() && value() == minimum()) return QString();
    return shortest(value());
  }
  static QString shortest(double x) { return QString::number(x, 'g', QLocale::FloatingPointShortest); }

  /// Acceptable: a number in range; Intermediate: on the way to one ("", "-", "1e", "2.", out of range).
  QValidator::State validate(QString &text, int &pos) const override {
    (void)pos;
    const QString t = bare(text);
    if (t.isEmpty()) return QValidator::Intermediate;
    bool ok = false;
    const double x = QLocale::c().toDouble(t, &ok);
    if (ok)
      return std::isfinite(x) && x >= minimum() && x <= maximum() ? QValidator::Acceptable : QValidator::Intermediate;
    static const QRegularExpression partial("^[+-]?[0-9]*\\.?[0-9]*([eE][+-]?[0-9]*)?$");
    return partial.match(t).hasMatch() ? QValidator::Intermediate : QValidator::Invalid;
  }
  double valueFromText(const QString &text) const override {
    bool ok = false;
    const double x = QLocale::c().toDouble(bare(text), &ok);
    return ok ? x : value();
  }

 protected:
  QString textFromValue(double value) const override { return shortest(value); }

 private:
  /// The text without the prefix and the suffix, trimmed.
  QString bare(const QString &text) const {
    QString t = text;
    if (!prefix().isEmpty() && t.startsWith(prefix())) t.remove(0, prefix().size());
    if (!suffix().isEmpty() && t.endsWith(suffix())) t.chop(suffix().size());
    return t.trimmed();
  }
  QString written_;
  double writtenValue_ = std::numeric_limits<double>::quiet_NaN();
  bool readable_ = false;
};

#endif
