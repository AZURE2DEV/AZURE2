#ifndef CONFIGDIALOGS_H
#define CONFIGDIALOGS_H

#include <QDialog>
#include <QJsonObject>
#include <QMap>
#include <QStringList>

#include "ProjectHolder.h"

/*!
 * The Configure menu's dialogs (Checks, Directories, Runtime Options): the <config> entries of the given groups, one
 * row each, drawn from the engine's table (azure2::front::ConfigTableJson: the old program's titles, the docs, the
 * defaults).  OK writes each changed entry with the engine's edit (ProjectHolder::setConfig): the line keeps its
 * comment, an option set to its default stays written.  An option that is absent shows its default; it is written
 * only when changed.
 */
class ConfigDialog : public QDialog {
  Q_OBJECT

 public:
  /// `groups`: the table's groups shown (empty: every group but those of `exclude`).
  ConfigDialog(ProjectHolder *holder, const QString &title, const QStringList &groups, const QStringList &exclude,
               QWidget *parent = nullptr);
  /// Write the changed entries; false with the engine's reason (the dialog stays open).
  bool write(QString &why);
  /// The editor of an entry, for the tests.
  QWidget *editor(const QString &name) const { return editors_.value(name, nullptr); }
  /// A row that is no <config> entry (a setting of this computer, as the Python interpreter of the Tools menu),
  /// below the entries; the caller reads it after exec().
  void addMachineRow(const QString &label, QWidget *editor);

 public slots:
  void accept() override;

 private:
  ProjectHolder *holder_;
  QMap<QString, QWidget *> editors_;
  QMap<QString, QString> initial_;
  QMap<QString, QJsonObject> rows_;
  class QFormLayout *machine_ = nullptr;
};

#endif
