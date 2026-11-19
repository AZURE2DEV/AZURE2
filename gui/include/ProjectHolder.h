#ifndef PROJECTHOLDER_H
#define PROJECTHOLDER_H

#include <QObject>
#include <QString>
#include <QStringList>

#include "azure2/front/Engine.h"

/*!
 * The open project of the window: one azure2::front::Project (the file's lossless document, edited with the engine's
 * grammars; docs/upgrade/qt_gui.md section 2) and the schema's tables, shared by the tabs.  Lives in the GUI thread;
 * a run reads the saved file, never this object.
 */
class ProjectHolder : public QObject {
  Q_OBJECT

 public:
  explicit ProjectHolder(QObject *parent = nullptr);

  azure2::front::Project &project() { return project_; }
  const azure2::front::Project &project() const { return project_; }
  /// Replace the project (after Open or New) and tell the tabs.
  void setProject(azure2::front::Project &&project);
  QString path() const { return QString::fromStdString(project_.Path()); }
  QString directory() const { return QString::fromStdString(project_.Directory()); }

  /// The schema's record kinds (in its order) and one kind by name (nullptr when unknown).
  const std::vector<azure2::front::KindInfo> &kinds() const { return kinds_; }
  const azure2::front::KindInfo *kind(const std::string &name) const;
  /// The kinds a section of records holds.
  QStringList kindsOfSection(const QString &section) const;

  /// Records of a section; set/remove/rename one (false with the engine's reason in `why`), each followed by
  /// sectionEdited(section).
  std::vector<azure2::front::Record> records(const QString &section) const;
  bool setRecord(const QString &section, const azure2::front::Record &record, QString &why);
  bool removeRecord(const QString &section, const QString &name, QString &why);
  bool renameRecord(const QString &section, const QString &from, const QString &to, QString &why);

  QString configValue(const QString &name, const QString &absent = QString()) const;
  bool setConfig(const QString &name, const QString &value, QString &why, bool remove = false);

  QStringList sectionBody(const QString &name) const;
  void setSectionBody(const QString &name, const QStringList &lines);

  /// The engine's diagnostics of the text as it stands (AZURE2 --validate).
  azure2::front::Validation validate() const { return project_.Validate(); }

 signals:
  /// A new project is in place: every tab reads it again.
  void opened();
  /// A section was edited by a tab (another tab that shows it reads it again).
  void sectionEdited(const QString &section);

 private:
  azure2::front::Project project_;
  std::vector<azure2::front::KindInfo> kinds_;
  std::vector<std::pair<std::string, std::vector<std::string>>> sections_;
};

#endif
