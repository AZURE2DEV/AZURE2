#include "AzureTheme.h"

#include <QApplication>
#include <QPalette>
#include <QStyleFactory>
#include <cstdlib>

void applyAzureTheme(QApplication &app) {
  const char *choice = std::getenv("AZURE2_THEME");
  if (choice && QString(choice).compare("system", Qt::CaseInsensitive) == 0) return;

  // The icon's "(A)" blue is about #2457a8; the accent is a touch lighter so highlighted text stays readable.
  const QColor azure("#2a6fc9"), azureDark("#1f4f8f");
  const QColor window("#e8f0f8"), panel("#f6f9fc"), base("#ffffff"), alternate("#f2f7fc");
  const QColor button("#d6e4f3"), line("#b7cbe0"), text("#1b2a3a"), disabled("#8a98a8");

  if (QStyle *fusion = QStyleFactory::create("Fusion")) app.setStyle(fusion);

  QPalette p;
  p.setColor(QPalette::Window, window);
  p.setColor(QPalette::WindowText, text);
  // Tables, lists and text areas take their background from the palette, slightly translucent so the
  // watermark behind the tab shows where they are empty. (A translucent background in the stylesheet
  // instead breaks the style's check-box indicators in item views.)
  p.setColor(QPalette::Base, QColor(base.red(), base.green(), base.blue(), 165));
  p.setColor(QPalette::AlternateBase, QColor(alternate.red(), alternate.green(), alternate.blue(), 165));
  p.setColor(QPalette::Text, text);
  p.setColor(QPalette::Button, button);
  p.setColor(QPalette::ButtonText, text);
  p.setColor(QPalette::Highlight, azure);
  p.setColor(QPalette::HighlightedText, Qt::white);
  p.setColor(QPalette::Link, azureDark);
  p.setColor(QPalette::ToolTipBase, QColor("#fffbe6"));
  p.setColor(QPalette::ToolTipText, text);
  p.setColor(QPalette::Light, QColor("#ffffff"));
  p.setColor(QPalette::Midlight, QColor("#e3edf8"));
  p.setColor(QPalette::Mid, line);
  p.setColor(QPalette::Dark, QColor("#9fb9d6"));
  p.setColor(QPalette::Shadow, QColor("#6f8aa8"));
  for (QPalette::ColorRole role : {QPalette::WindowText, QPalette::Text, QPalette::ButtonText})
    p.setColor(QPalette::Disabled, role, disabled);
  p.setColor(QPalette::Disabled, QPalette::Base, QColor("#eef2f6"));
  p.setColor(QPalette::Disabled, QPalette::Button, QColor("#e3eaf1"));
  app.setPalette(p);

  // Accents the palette alone cannot give. Combo and spin boxes are left to the style so their arrows survive.
  app.setStyleSheet(QString(R"(
    QTabWidget::pane { border: 1px solid %LINE%; top: -1px; background-color: %PANEL%;
                       background-image: url(:/azuma-background.png); background-repeat: no-repeat;
                       background-position: top left; }
    QTabBar::tab { background: %BUTTON%; color: %TEXT%; padding: 6px 11px; margin-right: 2px;
                   border: 1px solid %LINE%; border-bottom: none;
                   border-top-left-radius: 4px; border-top-right-radius: 4px; }
    QTabBar::tab:selected { background: %PANEL%; color: %AZUREDARK%; border-top: 3px solid %AZURE%; }
    QTabBar::tab:hover:!selected { background: #e3edf8; }
    QGroupBox { border: 1px solid %LINE%; border-radius: 5px; margin-top: 1.1em; padding-top: 0.3em; background: rgba(246, 249, 252, 0.55); }
    QGroupBox[title=""] { margin-top: 0.3em; }
    QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; color: %AZUREDARK%; }
    QHeaderView::section { background: %BUTTON%; color: %TEXT%; border: 1px solid %LINE%; padding: 3px; }
    QTableView, QTableWidget, QListWidget, QTreeView { gridline-color: #d5e2ef;
                   selection-background-color: %AZURE%; selection-color: white; border: 1px solid %LINE%; }
    QPushButton { background: #e4eef8; border: 1px solid #9fb9d6; border-radius: 4px; padding: 3px 6px; }
    QPushButton:hover { background: %BUTTON%; }
    QPushButton:pressed { background: #c3d6ec; }
    QPushButton:default { border: 2px solid %AZURE%; }
    QPushButton:checked { background: %AZURE%; color: white; }
    QPushButton:disabled { color: %DISABLED%; background: #eef2f6; border-color: #c9d6e3; }
    QLineEdit, QTextEdit, QPlainTextEdit { border: 1px solid %LINE%; border-radius: 3px; padding: 2px;
                   selection-background-color: %AZURE%; selection-color: white; }
    QLineEdit { background: %BASE%; }
    QLineEdit:disabled, QTextEdit:disabled { background: #eef2f6; color: %DISABLED%; }
    QLineEdit:focus, QTextEdit:focus { border: 1px solid %AZURE%; }
    QMenuBar { background: #dbe7f3; }
    QMenuBar::item:selected { background: %AZURE%; color: white; }
    QMenu::item:selected { background: %AZURE%; color: white; }
    QProgressBar { border: 1px solid %LINE%; border-radius: 3px; text-align: center; background: %BASE%; }
    QProgressBar::chunk { background: %AZURE%; }
    QToolTip { background: #fffbe6; color: %TEXT%; border: 1px solid #d9c98a; padding: 3px; }
    QStatusBar { background: #dbe7f3; }
  )")
                         .replace("%LINE%", line.name()).replace("%PANEL%", panel.name()).replace("%BUTTON%", button.name())
                         .replace("%TEXT%", text.name()).replace("%AZUREDARK%", azureDark.name()).replace("%AZURE%", azure.name())
                         .replace("%ALT%", alternate.name()).replace("%BASE%", base.name()).replace("%DISABLED%", disabled.name()));
}
