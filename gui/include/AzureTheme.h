#ifndef AZURETHEME_H
#define AZURETHEME_H

class QApplication;

/*!
 * The GUI's colour scheme: the Fusion style with an azure palette and a few accents (tab bar, group-box
 * titles, headers, buttons) in the blue of the AZURE2 icon, in place of the platform's default grey. The
 * watermark of R.E. Azuma from the AZURE website (azure.nd.edu, azuma_background2.gif) sits at the bottom
 * left of every tab, showing through the translucent tables and text areas where they are empty.
 *
 * Fusion is used on every platform so the colours look the same everywhere. AZURE2_THEME=system in the
 * environment leaves the platform style and palette alone.
 */
void applyAzureTheme(QApplication &app);

#endif
