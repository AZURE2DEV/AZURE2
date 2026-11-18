#ifndef FILEDIALOGPOLICY_H
#define FILEDIALOGPOLICY_H

/*!
 * Which file dialogs the GUI uses: the platform's own (Finder, GTK, Windows Explorer) or Qt's.
 *
 * On macOS the native dialogs are switched off for the whole application (Qt::AA_DontUseNativeDialogs), so every
 * "Choose..." and Open/Save dialog is Qt's own. With the Homebrew Qt 5.15 the .dmg is built from, the native
 * NSOpenPanel did not work for users on recent macOS (report of 2026-11-05): the "Choose..." buttons failed until
 * each call was given QFileDialog::DontUseNativeDialog, one button at a time. The application attribute covers all
 * of them, and any dialog added later, in one place. Linux and Windows keep their native dialogs.
 *
 * AZURE2_NATIVE_DIALOGS=1 (or 0) in the environment forces native (or Qt) dialogs on any platform: for a Mac user
 * who wants the Finder dialog back, and for the tests.
 */
bool useNativeFileDialogs();

/*!
 * Sets or clears Qt::AA_DontUseNativeDialogs to match useNativeFileDialogs(). Call it before the QApplication is
 * constructed (main.cpp does); it can also be called again later, and affects dialogs created after the call.
 */
void applyFileDialogPolicy();

#endif
