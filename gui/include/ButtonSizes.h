#ifndef BUTTONSIZES_H
#define BUTTONSIZES_H

#include <QAbstractButton>
#include <QFontMetrics>
#include <algorithm>

/*!
 * Button sizes from the font, not in pixels (docs/upgrade/qt_gui.md section 19.11): a button sized for one font
 * (setMaximumSize(50, 28)) cuts its label under a larger one ("NDS" at a 16-18 pt application font).
 */
namespace buttons {

/// The height of a row of small buttons: 28 px at the usual fonts, the style's height of a text button at larger ones.
inline int RowHeight(const QAbstractButton *b) { return std::max(28, b->sizeHint().height()); }

/// The width the button's label needs, with a margin of one line height (the style's padding).
inline int TextWidth(const QAbstractButton *b) {
  const QFontMetrics fm = b->fontMetrics();
  return fm.horizontalAdvance(b->text()) + fm.height();
}

/// A text button of a row: at least as wide as its label, as high as the row.
inline void FitText(QAbstractButton *b, int rowHeight) {
  b->setMinimumWidth(TextWidth(b));
  b->setFixedHeight(rowHeight);
}

/// A one-glyph button ("+", "-"): a square of the row's height, wide enough for its glyph.
inline void FitGlyph(QAbstractButton *b, int rowHeight) { b->setFixedSize(std::max(rowHeight, TextWidth(b)), rowHeight); }

}  // namespace buttons

#endif
