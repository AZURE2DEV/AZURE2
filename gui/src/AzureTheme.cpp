// The program's look: Fusion, a light and a dark palette, one style sheet over the palette (AzureTheme.h).
#include "AzureTheme.h"

#include <QApplication>
#include <QPalette>
#include <QSettings>
#include <QStyle>
#include <QStyleFactory>
#include <QStyleHints>
#include <QWidget>
#include <cstdlib>

namespace {

/// One look's colours, by role.  The light one keeps the azure of the AZURE2 icon (#2457a8) as its accent, a
/// touch lighter so highlighted text stays readable; the dark one the same hue, lighter still.
struct Colours {
  const char *window, *base, *alternate, *text, *button, *highlight, *highlightedText, *light, *midlight, *mid, *dark,
      *shadow, *toolTipBase, *toolTipText, *link, *placeholder, *disabledText, *disabledBase, *disabledButton, *indicator;
};

constexpr Colours kLight = {"#eef3f8", "#ffffff", "#f3f7fb", "#1b2a3a", "#e2ebf5", "#2a6fc9", "#ffffff",
                            "#ffffff", "#e3edf8", "#b7cbe0", "#9fb9d6", "#6f8aa8", "#fffbe6", "#1b2a3a",
                            "#1f4f8f", "#7a8796", "#8a98a8", "#eef2f6", "#e3eaf1", "#6f86a0"};
constexpr Colours kDark = {"#1f242b", "#161a1f", "#1c2127", "#e3e8ee", "#2c333c", "#3d8be0", "#ffffff",
                           "#3a434d", "#323a43", "#46515d", "#121519", "#0a0c0e", "#2c333c", "#e3e8ee",
                           "#7fb6ff", "#8592a0", "#6b7682", "#1b1f24", "#252b32", "#8a97a6"};

/// The style sheet: palette roles (palette(...); the three roles a style sheet cannot name as %TOKENS%, filled from
/// the look's table), lengths in em, borders in px.  Combo and spin boxes keep the
/// style's own frame so that their arrows survive.
const char *kStyleSheet = R"(
  QMainWindow::separator { width: 4px; height: 4px; }
  QTabWidget::pane { border: 1px solid palette(mid); top: -1px; background-color: palette(window); }
  QTabBar::tab { background: palette(button); color: palette(button-text); padding: 0.4em 0.9em;
                 margin-right: 2px; border: 1px solid palette(mid); border-bottom: none;
                 border-top-left-radius: 4px; border-top-right-radius: 4px; }
  QTabBar::tab:selected { background: palette(window); color: palette(link); border-top: 3px solid palette(highlight); }
  QTabBar::tab:hover:!selected { background: palette(midlight); }
  QGroupBox { border: 1px solid palette(mid); border-radius: 5px; margin-top: 1.1em; padding: 0.4em 0.3em 0.3em 0.3em; }
  QGroupBox[title=""] { margin-top: 0.3em; }
  QGroupBox::title { subcontrol-origin: margin; left: 0.7em; padding: 0 0.3em; color: palette(link); }
  QHeaderView::section { background: palette(button); color: palette(button-text); padding: 0.3em 0.5em;
                         border: none; border-right: 1px solid palette(mid); border-bottom: 1px solid palette(mid);
                         font-weight: 600; }
  QAbstractItemView { qproperty-alternatingRowColors: 1; }
  QTableView, QTreeView, QListView { gridline-color: palette(midlight); border: 1px solid palette(mid);
                 selection-background-color: palette(highlight); selection-color: palette(highlighted-text); }
  QPushButton { background: palette(button); color: palette(button-text); border: 1px solid palette(mid);
                border-radius: 4px; padding: 0.3em 0.8em; min-height: 1.3em; }
  QPushButton:hover { background: palette(midlight); border-color: palette(highlight); }
  QPushButton:pressed { background: palette(mid); }
  QPushButton:focus { border: 1px solid palette(highlight); }
  QPushButton:default { border: 2px solid palette(highlight); }
  QPushButton:checked { background: palette(highlight); color: palette(highlighted-text); }
  QPushButton:disabled { color: %PLACEHOLDER%; background: palette(window); border-color: palette(midlight); }
  QPushButton[primary="true"] { background: palette(highlight); color: palette(highlighted-text);
                                border: 1px solid palette(highlight); font-weight: 600; }
  QPushButton[primary="true"]:hover { background: palette(link); }
  QPushButton[primary="true"]:disabled { background: palette(window); color: %PLACEHOLDER%;
                                         border-color: palette(midlight); font-weight: normal; }
  QPushButton[danger="true"]:enabled { background: #c62828; color: white; border: 1px solid #a51d1d; font-weight: 600; }
  QPushButton[danger="true"]:enabled:hover { background: #a51d1d; }
  QToolBar { spacing: 0.2em; padding: 0.15em 0.3em; border: none; border-bottom: 1px solid palette(mid); }
  QToolButton { padding: 0.2em 0.35em; border: 1px solid transparent; border-radius: 4px; }
  QToolButton:hover { background: palette(midlight); border-color: palette(mid); }
  QToolButton:pressed, QToolButton:checked { background: palette(mid); }
  QLineEdit, QTextEdit, QPlainTextEdit { border: 1px solid palette(mid); border-radius: 3px; padding: 0.15em;
                 selection-background-color: palette(highlight); selection-color: palette(highlighted-text); }
  QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus { border: 1px solid palette(highlight); }
  QLineEdit:disabled, QTextEdit:disabled, QPlainTextEdit:disabled { color: %PLACEHOLDER%; }
  QMenuBar { background: palette(window); border-bottom: 1px solid palette(midlight); }
  QMenuBar::item { padding: 0.25em 0.6em; }
  QMenuBar::item:selected, QMenu::item:selected { background: palette(highlight); color: palette(highlighted-text); }
  QProgressBar { border: 1px solid palette(mid); border-radius: 3px; text-align: center; background: palette(base);
                 min-height: 1.1em; }
  QProgressBar::chunk { background: palette(highlight); border-radius: 2px; }
  QToolTip { background: %TIPBASE%; color: %TIPTEXT%; border: 1px solid palette(mid);
             padding: 0.25em; }
  QStatusBar { background: palette(window); border-top: 1px solid palette(mid); }
  QStatusBar::item { border: none; }
  QDockWidget::title { background: palette(button); padding: 0.25em 0.5em; border-top: 1px solid palette(mid); }
  QLabel[role="hint"] { color: %PLACEHOLDER%; }
  QCheckBox::indicator, QGroupBox::indicator, QAbstractItemView::indicator {
                 width: 0.85em; height: 0.85em; border: 1px solid %INDICATOR%; border-radius: 3px; background: palette(base); }
  QCheckBox::indicator:hover, QAbstractItemView::indicator:hover { border-color: palette(highlight); }
  QCheckBox::indicator:checked, QGroupBox::indicator:checked, QAbstractItemView::indicator:checked {
                 background: palette(highlight); border-color: palette(highlight); image: url(:/actions/check.svg); }
  QCheckBox::indicator:indeterminate, QAbstractItemView::indicator:indeterminate {
                 background: palette(highlight); border-color: palette(highlight); image: url(:/actions/partial.svg); }
  QCheckBox::indicator:disabled, QAbstractItemView::indicator:disabled { border-color: palette(midlight); background: palette(window); }
  QCheckBox::indicator:checked:disabled { background: palette(mid); border-color: palette(mid); }
  QRadioButton::indicator { width: 0.8em; height: 0.8em; border: 1px solid %INDICATOR%; border-radius: 0.45em;
                 background: palette(base); }
  QRadioButton::indicator:hover { border-color: palette(highlight); }
  QRadioButton::indicator:checked { border: 1px solid palette(highlight);
                 background: qradialgradient(cx: 0.5, cy: 0.5, radius: 0.5, fx: 0.5, fy: 0.5, stop: 0 %ACCENT%,
                                             stop: 0.5 %ACCENT%, stop: 0.6 %BASE%, stop: 1 %BASE%); }
  QLabel[role="title"] { color: palette(link); }
  QFrame[role="card"] { background: palette(base); border: 1px solid palette(mid); border-radius: 6px; }
)";

/// The watermark of the light look (the old program's), under every tab.
const char *kWatermark = R"(
  QTabWidget::pane { background-image: url(:/azuma-background.png); background-repeat: no-repeat;
                     background-position: bottom left; }
)";

bool PlatformIsDark(const QApplication &app) {
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
  const Qt::ColorScheme scheme = QGuiApplication::styleHints()->colorScheme();
  if (scheme != Qt::ColorScheme::Unknown) return scheme == Qt::ColorScheme::Dark;
#endif
  // the platform's own window colour, kept when the theme was first applied
  return app.property("azure2PlatformDark").toBool();
}

void Apply(QApplication &app, ThemeChoice choice) {
  const bool dark = choice == ThemeChoice::Dark || (choice == ThemeChoice::System && PlatformIsDark(app));
  if (QStyle *fusion = QStyleFactory::create("Fusion")) app.setStyle(fusion);
  app.setPalette(azurePalette(dark));
  const Colours &c = dark ? kDark : kLight;
  app.setStyleSheet((QString(kStyleSheet) + (dark ? QString() : QString(kWatermark)))
                        .replace("%PLACEHOLDER%", c.placeholder)
                        .replace("%TIPBASE%", c.toolTipBase)
                        .replace("%TIPTEXT%", c.toolTipText)
                        .replace("%INDICATOR%", c.indicator)
                        .replace("%ACCENT%", c.highlight)
                        .replace("%BASE%", c.base));
  app.setProperty("azure2ThemeDark", dark);
  app.setProperty("azure2ThemeChoice", themeChoiceName(choice));
}

}  // namespace

QString themeChoiceName(ThemeChoice choice) {
  return choice == ThemeChoice::Light ? "light" : choice == ThemeChoice::Dark ? "dark" : "system";
}

ThemeChoice themeChoiceFromName(const QString &name) {
  if (name.compare("light", Qt::CaseInsensitive) == 0) return ThemeChoice::Light;
  if (name.compare("dark", Qt::CaseInsensitive) == 0) return ThemeChoice::Dark;
  return ThemeChoice::System;
}

QPalette azurePalette(bool dark) {
  const Colours &c = dark ? kDark : kLight;
  QPalette p;
  p.setColor(QPalette::Window, QColor(c.window));
  p.setColor(QPalette::WindowText, QColor(c.text));
  // the light look's tables are slightly translucent so the watermark shows where they are empty (a translucent
  // background in the style sheet instead breaks the check-box indicators of item views)
  QColor base(c.base), alternate(c.alternate);
  if (!dark) {
    base.setAlpha(200);
    alternate.setAlpha(200);
  }
  p.setColor(QPalette::Base, base);
  p.setColor(QPalette::AlternateBase, alternate);
  p.setColor(QPalette::Text, QColor(c.text));
  p.setColor(QPalette::Button, QColor(c.button));
  p.setColor(QPalette::ButtonText, QColor(c.text));
  p.setColor(QPalette::BrightText, QColor("#ff5252"));
  p.setColor(QPalette::Highlight, QColor(c.highlight));
  p.setColor(QPalette::HighlightedText, QColor(c.highlightedText));
  p.setColor(QPalette::Link, QColor(c.link));
  p.setColor(QPalette::LinkVisited, QColor(c.link));
  p.setColor(QPalette::ToolTipBase, QColor(c.toolTipBase));
  p.setColor(QPalette::ToolTipText, QColor(c.toolTipText));
  p.setColor(QPalette::PlaceholderText, QColor(c.placeholder));
  p.setColor(QPalette::Light, QColor(c.light));
  p.setColor(QPalette::Midlight, QColor(c.midlight));
  p.setColor(QPalette::Mid, QColor(c.mid));
  p.setColor(QPalette::Dark, QColor(c.dark));
  p.setColor(QPalette::Shadow, QColor(c.shadow));
  for (QPalette::ColorRole role : {QPalette::WindowText, QPalette::Text, QPalette::ButtonText})
    p.setColor(QPalette::Disabled, role, QColor(c.disabledText));
  p.setColor(QPalette::Disabled, QPalette::Base, QColor(c.disabledBase));
  p.setColor(QPalette::Disabled, QPalette::Button, QColor(c.disabledButton));
  p.setColor(QPalette::Disabled, QPalette::Highlight, QColor(c.mid));
  return p;
}

void applyAzureTheme(QApplication &app) {
  const QString env = QString::fromLocal8Bit(std::getenv("AZURE2_THEME"));
  if (env.compare("native", Qt::CaseInsensitive) == 0 || env.compare("system", Qt::CaseInsensitive) == 0) return;
  if (!app.property("azure2PlatformDark").isValid()) {
    // before the first palette of ours: whether the platform's own window is dark
    app.setProperty("azure2PlatformDark", QGuiApplication::palette().color(QPalette::Window).lightness() < 128);
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
    QObject::connect(QGuiApplication::styleHints(), &QStyleHints::colorSchemeChanged, &app, [&app]() {
      if (azureThemeChoice() == ThemeChoice::System) Apply(app, ThemeChoice::System);
    });
#endif
  }
  ThemeChoice choice = themeChoiceFromName(QSettings().value("view/theme", "system").toString());
  if (env.compare("light", Qt::CaseInsensitive) == 0) choice = ThemeChoice::Light;
  if (env.compare("dark", Qt::CaseInsensitive) == 0) choice = ThemeChoice::Dark;
  Apply(app, choice);
}

void setAzureTheme(ThemeChoice choice) {
  QSettings().setValue("view/theme", themeChoiceName(choice));
  if (QApplication *app = qobject_cast<QApplication *>(QCoreApplication::instance())) Apply(*app, choice);
}

ThemeChoice azureThemeChoice() {
  return themeChoiceFromName(qApp ? qApp->property("azure2ThemeChoice").toString() : QString());
}

bool azureThemeIsDark() { return qApp && qApp->property("azure2ThemeDark").toBool(); }
