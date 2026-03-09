#include "RichTextDelegate.h"
#include <QTextDocument>
#include <QPainter>
#include <QAbstractTextDocumentLayout>
#include <QApplication>

void RichTextDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const {
  QTextDocument document;
  QString value = index.data(Qt::DisplayRole).toString();
  
  document.setHtml(value);
  QStyleOptionViewItemV4 opt(option);
  initStyleOption(&opt, index);

  document.setTextWidth(opt.rect.width());

  opt.text=QString();
  QStyle *style = opt.widget ? opt.widget->style() : QApplication::style();
  style->drawControl(QStyle::CE_ItemViewItem, &opt, painter);

  QAbstractTextDocumentLayout::PaintContext ctx;
  if (opt.state & QStyle::State_Selected) {
    if(opt.state & QStyle::State_Active)
      ctx.palette.setColor(QPalette::Text, opt.palette.color(QPalette::Active, QPalette::HighlightedText));
    else ctx.palette.setColor(QPalette::Text, opt.palette.color(QPalette::Inactive, QPalette::HighlightedText));
  }
  
  painter->translate(opt.rect.topLeft());
  document.documentLayout()->draw(painter,ctx);
  painter->translate(-opt.rect.topLeft());
}

QSize RichTextDelegate::sizeHint(const QStyleOptionViewItem & option, const QModelIndex & index ) const {
  QTextDocument document;
  QString value = index.data(Qt::DisplayRole).toString();
  document.setHtml(value);
  document.setTextWidth(option.rect.width());
  return QSize(option.rect.width(),document.size().height());
} 

