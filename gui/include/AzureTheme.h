#ifndef AZURETHEME_H
#define AZURETHEME_H

#include <QPalette>
#include <QString>

class QApplication;

/*!
 * The program's look (docs/upgrade/qt_gui.md section 11): the Fusion style on every platform, a light and a dark
 * palette built from one table of roles, and one restrained style sheet over the palette (spacing in em, readable
 * table headers, alternating rows in every item view, a visible focus frame, primary buttons in the accent colour,
 * a Stop button in red).  The style sheet names palette roles, never a colour of its own, so one sheet serves both
 * palettes and reaches every tab's widgets.  In the light palette the watermark of R.E. Azuma (azure.nd.edu) sits at
 * the bottom left of every tab, showing through the translucent tables where they are empty.
 *
 * The choice (System, Light, Dark) is the View menu's, kept in QSettings ("view/theme").  System follows the
 * platform's colour scheme: QStyleHints::colorScheme where Qt >= 6.5 has it, else the lightness of the platform's own
 * window colour.  AZURE2_THEME=native in the environment leaves the platform style and palette alone (the old
 * AZURE2_THEME=system does the same).
 *
 * A widget asks for an accent with a dynamic property: primary=true (Run, Save and Run, the start page's buttons),
 * danger=true (Stop: red while enabled).
 */
enum class ThemeChoice { System, Light, Dark };

/// Apply the saved choice (or AZURE2_THEME) to the application; called once in main() and by the tests.
void applyAzureTheme(QApplication &app);
/// Apply and keep a choice (the View menu).
void setAzureTheme(ThemeChoice choice);
ThemeChoice azureThemeChoice();
/// The palette in use is the dark one.
bool azureThemeIsDark();
/// The palette of the light or the dark look.
QPalette azurePalette(bool dark);
/// "system", "light", "dark" and back.
QString themeChoiceName(ThemeChoice choice);
ThemeChoice themeChoiceFromName(const QString &name);

#endif
