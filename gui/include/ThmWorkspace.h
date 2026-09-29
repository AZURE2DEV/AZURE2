#ifndef THMWORKSPACE_H
#define THMWORKSPACE_H

#include <QDialog>

#include "ThmSettings.h"

QT_BEGIN_NAMESPACE
class QLabel;
class QPushButton;
class QTabWidget;
QT_END_NAMESPACE

class AZURESetup;
class ThmChannelsPage;
class ThmExperimentsPage;
class ThmModelPage;

/*!
 * The THM workspace (Configure > THM Workspace...): everything that belongs
 * to the Trojan Horse (HOES) observable, in one window, apart from the classic
 * editor.  Pages:
 *   Model        the options of the <thm> block (ThmModelPage);
 *   Experiments  its experiment[<name>] lines (ThmExperimentsPage);
 *   Channels     the THM columns of the <levels> lines: binding energy B
 *                (field 32) and width input flag (field 33) (ThmChannelsPage).
 * It writes nothing else, and nothing at all when its values are left as
 * they were read: such a project saves byte for byte as before.  Without a
 * THM segment the pages are disabled under a short explanation.
 */
class ThmWorkspace : public QDialog {
  Q_OBJECT

 public:
  /// `settings`: the project's <thm> block as read (AZURESetup::thmSettings).
  ThmWorkspace(AZURESetup *setup, const ThmSettings &settings, QWidget *parent = 0);

  /// True if the project has a THM data or test segment (the pages are enabled).
  bool hasThmSegments() const { return hasThm_; }
  /// The first thing AZURE2 would refuse, page by page ("" if none); shows that page.
  QString validate();
  /// Installs the pages' values into the project (no check; see validate).
  void apply();

  QLabel *noThmLabel;
  QTabWidget *pages;
  ThmModelPage *modelPage;
  ThmExperimentsPage *experimentsPage;
  ThmChannelsPage *channelsPage;
  QPushButton *acceptButton;

 public slots:
  void accept() override;

 private:
  AZURESetup *setup_;
  bool hasThm_;
};

#endif
