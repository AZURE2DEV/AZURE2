#ifndef NODETABLECACHE_H
#define NODETABLECACHE_H

#include "EPoint.h"

#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <unordered_map>

/*!
 * The energy-dependent tables of re-placed sub-points, by node, shared by an
 * EData and all its copies (EData::Clone copies the pointer).
 *
 * Every pooled copy of a fit re-places the sub-point grids for its own
 * parameters (EData::ReplaceSubPointGrids).  A node that moved gets new tables,
 * and a copy used to keep its own for every node that no longer coincided with
 * one of the data it was cloned from: memory grew towards a full copy per
 * pooled copy as the levels drifted (natural neon: 11 -> 17 GB in 50 min).  But
 * the copies evaluate parameter sets that differ in one parameter, and a grid's
 * nodes depend only on the levels inside its window, so most copies place most
 * nodes at exactly the same energies.  The first copy to compute a node
 * publishes it here; the others take its tables (EPoint::AdoptEDependentValues
 * re-checks energy, angle and pairs exactly).
 *
 * Keyed by (segment key, c.m. energy, c.m. angle), bit for bit: the Legendre
 * polynomials depend on the segment's target effect (Q coefficients).  Holds
 * the tables weakly: it keeps nothing alive, and an entry whose tables no copy
 * uses any more is a miss and is swept out.  Sharded, one mutex per shard.
 */
class NodeTableCache {
 public:
  bool Find(int segmentKey, double energy, double angle, EPoint::EDependentShare &out) const {
    const Key key = MakeKey(segmentKey, energy, angle);
    const Shard &s = shards_[ShardOf(key)];
    std::lock_guard<std::mutex> lock(s.mutex);
    auto it = s.map.find(key);
    if (it == s.map.end()) return false;
    out = it->second.share;
    out.tables = it->second.tables.lock();
    out.legendre = it->second.legendre.lock();
    return out.tables && out.legendre;
  }
  void Store(int segmentKey, const EPoint::EDependentShare &share) {
    if (!share.tables || !share.legendre) return;
    const Key key = MakeKey(segmentKey, share.energy, share.angle);
    Entry e;
    e.share = share;
    e.share.tables.reset();
    e.share.legendre.reset();
    e.tables = share.tables;
    e.legendre = share.legendre;
    Shard &s = shards_[ShardOf(key)];
    std::lock_guard<std::mutex> lock(s.mutex);
    s.map[key] = std::move(e);
    if (s.map.size() > s.sweepAt) {  // drop the entries no copy uses any more
      for (auto it = s.map.begin(); it != s.map.end();)
        it = it->second.tables.expired() ? s.map.erase(it) : std::next(it);
      s.sweepAt = 2 * s.map.size() + 1024;
    }
  }

 private:
  struct Key {
    int segment;
    std::uint64_t energy, angle;
    bool operator==(const Key &o) const { return segment == o.segment && energy == o.energy && angle == o.angle; }
  };
  struct KeyHash {
    size_t operator()(const Key &k) const {
      std::uint64_t h = k.energy * 0x9E3779B97F4A7C15ULL ^ (k.angle + 0x632BE59BD9B4E019ULL + (h0(k.segment) << 6));
      return (size_t)(h ^ (h >> 29));
    }
    static std::uint64_t h0(int s) { return (std::uint64_t)(unsigned)s * 0xBF58476D1CE4E5B9ULL; }
  };
  struct Entry {
    EPoint::EDependentShare share;  ///< without its pointers
    std::weak_ptr<const EPoint::EDependentTables> tables;
    std::weak_ptr<const vector_r> legendre;
  };
  struct Shard {
    mutable std::mutex mutex;
    std::unordered_map<Key, Entry, KeyHash> map;
    size_t sweepAt = 1024;
  };
  static constexpr int kShards = 16;
  static Key MakeKey(int segmentKey, double energy, double angle) {
    Key k;
    k.segment = segmentKey;
    std::memcpy(&k.energy, &energy, sizeof(double));
    std::memcpy(&k.angle, &angle, sizeof(double));
    return k;
  }
  static int ShardOf(const Key &k) { return (int)((KeyHash()(k) >> 7) % kShards); }
  Shard shards_[kShards];
};

#endif
