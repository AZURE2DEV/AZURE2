#ifndef RICHTEXTDELEGATE_H
#define RICHTEXTDELEGATE_H

#include <QStyledItemDelegate>

class QPainter;

class RichTextDelegate : public QStyledItemDelegate {

 protected:
  void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;
  QSize sizeHint ( const QStyleOptionViewItem & option, const QModelIndex & index ) const;

};

#endif
