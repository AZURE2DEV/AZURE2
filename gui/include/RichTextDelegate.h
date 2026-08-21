#ifndef RICHTEXTDELEGATE_H
#define RICHTEXTDELEGATE_H

#include <QStyledItemDelegate>

class QPainter;

/*!
 * Item delegate that renders a cell's rich text -- the superscripts and subscripts the pair and channel labels use.
 */
class RichTextDelegate : public QStyledItemDelegate {
 protected:
  void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;
  QSize sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const;
};

#endif
