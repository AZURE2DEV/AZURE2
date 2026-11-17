// The energy-dependent tables of a point and their shared form: the channel rows (EDependentTables), what a
// point at the same energy can take from another (EDependentShare), the THM vertex pieces (ThmVertexRows) and the
// THM spectator-window table (ThmPsTable).  Unit names them as its own types (Unit::EDependentTables, ...).
#ifndef EPOINTTABLES_H
#define EPOINTTABLES_H

#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

#include "Constants.h"
#include "SharedTable.h"

struct ThmPsTable;  // defined below, named by EDependentShare

/// The energy-dependent channel tables of a point (one row per J group); see Site::Values::edep (include/Site.h).
/// Packed: one flat array per kind and one row index, five allocations a point
/// instead of four per J group (natural neon: ~60 a sub-point).  Row j (from 1)
/// holds the channels of J group j, at slots [rowStart[j-1], rowStart[j]).
/// They depend on the compound energy E_x alone (and, through a beta-decay channel's Fermi function, on the exit
/// pair's charge), never on the point's entrance pair or route: every site at one E_x holds one table
/// (WorkingCopy::FillChannelRows), the THM vertex pieces apart (ThmVertexRows).
struct EDependentTables {
  std::vector<int> rowStart{0};
  vector_c lo_elements;
  vector_r penetrabilities;
  vector_c coulombphase;
  vector_c hardspherephase;
  /// Slot of (J group, channel); rows and channels are filled in order and extended here.
  int Slot(int jGroupNum, int channelNum) {
    while ((int)rowStart.size() < jGroupNum + 1) rowStart.push_back(rowStart.back());
    const int slot = rowStart[jGroupNum - 1] + channelNum - 1;
    if (slot >= rowStart[jGroupNum]) {
      assert(jGroupNum + 1 == (int)rowStart.size() && slot == rowStart[jGroupNum]);
      rowStart[jGroupNum] = slot + 1;
    }
    return slot;
  }
  int At(int jGroupNum, int channelNum) const { return rowStart[jGroupNum - 1] + channelNum - 1; }
  template <class V, class T>
  static void Put(V &v, int slot, const T &x) {
    if ((int)v.size() <= slot) v.resize(slot + 1);
    v[slot] = x;
  }
  void Clear() {
    rowStart.assign(1, 0);
    lo_elements.clear();
    penetrabilities.clear();
    coulombphase.clear();
    hardspherephase.clear();
  }
  void ShrinkToFit() {
    rowStart.shrink_to_fit();
    lo_elements.shrink_to_fit();
    penetrabilities.shrink_to_fit();
    coulombphase.shrink_to_fit();
    hardspherephase.shrink_to_fit();
  }
};
/// The THM (HOES) entrance-vertex pieces of a THM point (Site::GetThmFormFactor), for the channels of its entrance
/// pair only (no other channel's vertex is read): j_l(rho) and rho dj_l/drho, boundary-independent (M_l is assembled
/// at the vertex with the per-level or channel boundary), and the external Coulomb term (0 when not used).  The
/// entries of J group j (from 1) are [groupStart[j-1], groupStart[j]), in channel order, with their channel numbers in
/// `channel`.  Immutable once built, one per point (shared by its copies).
struct ThmVertexRows {
  std::vector<int> groupStart{0};
  std::vector<int> channel;
  vector_r jl, rhodjl;
  vector_c coul;
  /// eps = E + B_xs + T_s (MeV) of the entrance pair at this point: p^2/2mu of the vertex's half-off-shell momentum,
  /// which the full form factor's volume term reads (vertexModel=ff, ThmFormFactor.h); 0 below it.
  double eps = 0.0;
  /// The entry of (J group, channel), both 1-based, or -1 when the channel is not the entrance pair's.
  int Find(int jGroupNum, int channelNum) const {
    if (jGroupNum < 1 || jGroupNum >= (int)groupStart.size()) return -1;
    for (int k = groupStart[jGroupNum - 1]; k < groupStart[jGroupNum]; k++)
      if (channel[k] == channelNum) return k;
    return -1;
  }
  size_t Bytes() const {
    return sizeof(*this) + (groupStart.capacity() + channel.capacity()) * sizeof(int) +
           (jl.capacity() + rhodjl.capacity()) * sizeof(double) + coul.capacity() * sizeof(complex);
  }
};

/// Everything CalcEDependentValues and CalcLegendreP compute for a point, by
/// reference: the tables are shared, not copied.  These depend only on the
/// point's c.m. energy, angle and pairs, never on a fit parameter, so a
/// re-placed sub-point at an unchanged node takes them from the sub-point it
/// replaces (WorkingCopy::ReplaceSubPointGrids).
struct EDependentShare {
  double energy = 0.0;
  double angle = 0.0;
  int entranceKey = 0, exitKey = 0;
  double geofactor = 0.0, sfactorconv = 0.0, screen_k = 0.0;
  SharedTable<EDependentTables> tables;
  SharedTable<vector_r> legendre;
  std::shared_ptr<const ThmPsTable> thmPs;  ///< THM spectator-window table (immutable), empty for a classic point
  std::shared_ptr<const ThmVertexRows> thmVertex;  ///< THM vertex pieces (immutable), empty for a classic point
};

/// With a spectator-momentum window (ps=...): the same three pieces per
/// node of the window at the point's energy (the accepted directions,
/// ThmSpectatorWindow::Nodes) with its weight, for the entrance-pair channels only
/// (every other channel is 0).  One flat, immutable table per point,
/// shared by its mapped points; the single-node matrices above are then
/// unused (zero).  It used to be [node][jGroup][channel] vectors of every
/// channel: ~117 MB per node for the 19F THM model (16 nodes did not fit).
struct ThmPsTable {
  int nodes = 0;
  int nSlots = 0;                 ///< entrance channels stored (not `slots`: a Qt keyword in the GUI)
  std::vector<int> offset;        ///< J group j -> first flat index of its channels (j = 1..N, N+1 entries)
  std::vector<int> slot;          ///< flat (J group, channel) -> slot, -1 when not stored
  std::vector<double> jl;         ///< [slot * nodes + node]
  std::vector<double> rhodjl;
  std::vector<complex> coul;      ///< empty unless an external Coulomb term was computed
  std::vector<double> weight;     ///< [node], normalized
  std::vector<double> eps;        ///< [node]: E + B_xs + T_k (MeV), as ThmVertexRows::eps; 0 below it
  size_t Bytes() const {
    return sizeof(*this) + offset.capacity() * sizeof(int) + slot.capacity() * sizeof(int) +
           (jl.capacity() + rhodjl.capacity() + weight.capacity() + eps.capacity()) * sizeof(double) +
           coul.capacity() * sizeof(complex);
  }
};

#endif
