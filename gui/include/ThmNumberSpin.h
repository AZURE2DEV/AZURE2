#ifndef THMNUMBERSPIN_H
#define THMNUMBERSPIN_H

#include <QDoubleSpinBox>
#include <QLocale>
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
 * special value text, the minimum stands for "not given" and gives "".
 * No Q_OBJECT: it adds no signal or slot.
 */
class ThmNumberSpin : public QDoubleSpinBox {
 public:
  explicit ThmNumberSpin(const QString &suffix = QString(), double lo = 0.0, double hi = 1.0e6, double step = 1.0,
                         QWidget *parent = nullptr) :
    QDoubleSpinBox(parent) {
    setLocale(QLocale::c());
    setDecimals(10);  // the value keeps what is typed; the text is the shortest
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

 protected:
  QString textFromValue(double value) const override { return shortest(value); }

 private:
  QString written_;
  double writtenValue_ = std::numeric_limits<double>::quiet_NaN();
  bool readable_ = false;
};

#endif
