#ifndef AZURE2_EXPERIMENT_CHANNELRADII_H
#define AZURE2_EXPERIMENT_CHANNELRADII_H

#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * One `channel_radius` record of <channelRadii>: the radius of the channels (l, and where given s and J) of one
 * particle pair, in place of the pair's radius (column 28 of <levels>).  Every other channel keeps its pair's.
 */
struct ChannelRadiusEntry {
  int pair_key = 0;        ///< the pair key of <levels>
  int l = 0;               ///< orbital angular momentum
  double radius_fm = 0.0;  ///< the channels' radius (fm)
  bool has_s = false;      ///< the record names one channel spin
  int two_s = 0;           ///< 2s
  bool has_J = false;      ///< the record names one total J
  int two_J = 0;           ///< 2J
  std::string record;      ///< the record's name ("" for an entry of the API)
  int line = 0;            ///< its line in the project file (0 for an entry of the API)

  /// How many of s and J the entry names: the most specific entry that matches a channel gives its radius.
  int Specificity() const { return (has_s ? 1 : 0) + (has_J ? 1 : 0); }
  /// The entry names a channel (l, 2s) of a group 2J of pair `key`.
  bool Matches(int key, int channelL, int twoS, int twoJ) const {
    return key == pair_key && channelL == l && (!has_s || twoS == two_s) && (!has_J || twoJ == two_J);
  }
  /// Two entries name the same channels (same pair, l, and the same s and J or absence of them).
  bool SameSelection(const ChannelRadiusEntry &o) const {
    return pair_key == o.pair_key && l == o.l && has_s == o.has_s && has_J == o.has_J && (!has_s || two_s == o.two_s) &&
           (!has_J || two_J == o.two_J);
  }
  /// "pair 1, l = 2, J = 5/2" for messages.
  std::string Text() const;
};

/*!
 * The entry that gives channel (l, 2s) of a group 2J of pair `key` its radius: the most specific that matches
 * (l+s+J, then l+J or l+s, then l); nullptr when none does.  `ambiguous` is set when an entry with s only and an
 * entry with J only both match and none names both: the channel's radius is then not defined, and the set-up refuses.
 */
const ChannelRadiusEntry *ChannelRadiusFor(const std::vector<ChannelRadiusEntry> &entries, int key, int l, int twoS,
                                           int twoJ, bool *ambiguous);

}  // namespace experiment
}  // namespace azure2

#endif
