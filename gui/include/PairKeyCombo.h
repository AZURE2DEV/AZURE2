#ifndef PAIRKEYCOMBO_H
#define PAIRKEYCOMBO_H

#include <QComboBox>

class PairsModel;

/*!
 * A drop-down that stands in for a typed particle-pair key.
 *
 * The .azr and the engine identify a pair by its 1-based position in the pairs list, and the Particle Pairs tab
 * shows that position only as a row number, which is easy to miss. This box lists the pairs by name instead and
 * hands back the key. In Entrance mode an item reads "2: 14N + n"; in Exit mode it reads as the reaction the
 * chosen entrance pair makes with it, "2: 11B(a,n)14N [0.000 MeV]", so two exit pairs that differ only in the
 * excitation energy of the residual are told apart at a glance. The items are rebuilt from the pairs model each
 * time the list drops down, so a pair added or edited on its tab shows up without any other wiring.
 *
 * An exit box can also offer a last item, "Total capture (all γ pairs)", whose key is TotalCapture (-1): the
 * exit key a total-capture segment line has always carried. It is listed only while the pairs include a gamma
 * pair (pair type 10), or while it is the selected key (an edited line keeps its value).
 */
class PairKeyCombo : public QComboBox {
  Q_OBJECT

 public:
  enum Mode { Entrance, Exit };
  /*! The key of the total-capture item: the exit key a total-capture segment line is written with. */
  static const int TotalCapture = -1;
  PairKeyCombo(Mode mode, QWidget *parent = 0);

  void setPairsModel(PairsModel *model);
  /*! 1-based pair key; 0 (or any key below 1) clears the selection. A key beyond the pairs list is kept, marked.
   *  TotalCapture selects the total-capture item when the box offers it (and clears the selection otherwise). */
  void setKey(int key);
  int key() const;
  /*! Exit mode: list the "Total capture (all γ pairs)" item (key TotalCapture) when the pairs include a gamma pair. */
  void setTotalCaptureOffered(bool offered);
  /*! True when the pair with this key is a gamma (capture) pair, pair type 10. */
  bool isGammaPair(int key) const;
  /*! True when the pairs model holds at least one gamma pair. */
  bool hasGammaPair() const;
  /*! Rebuilds the items from the pairs model, keeping the current key. */
  void refresh();
  void showPopup() override;

 public slots:
  /*! Exit mode only: the entrance pair the reaction labels are written against (0 = unknown, pair names are used). */
  void setEntranceKey(int key);

 signals:
  void keyChanged(int key);

 private:
  QString labelFor(int key) const;

  Mode mode;
  PairsModel *pairsModel = nullptr;
  int entranceKey = 0;
  bool totalCaptureOffered = false;
};

#endif
