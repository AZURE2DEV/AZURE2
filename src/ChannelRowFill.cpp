// The channel rows of the set-up (docs/upgrade/channel_table.md, docs/upgrade/data_model.md section 9): the stage's
// channel table and the fill of every site's rows, one table per compound energy.
#include <algorithm>
#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "ChannelTable.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "GSLException.h"
#include "Segment.h"
#include "Unit.h"
#include "WorkingCopy.h"
#include "azure2/support/Environment.h"
#include "azure2/support/FloatingPointPrecision.h"

/*!
 * The channel table of the rows the given sites will compute, built before the threads that fill them start
 * (docs/upgrade/channel_table.md).  Rows the current table holds are copied from it.
 */
std::shared_ptr<const ChannelTable> WorkingCopy::BuildChannelTable(const CompoundSystem *compound,
                                                                   const Config &configure,
                                                                   const std::vector<Site *> &sites) {
  // Each site's keys are made unique before they are added: its slots repeat (pair, l, E_c) once per channel spin
  // and J^pi block, and the full list was the peak memory of the build (natural_neon: 11 M keys).
  std::vector<ChannelTable::Key> keys, own;
  std::size_t requested = 0;
  for (Site *site : sites) {
    own.clear();
    site->ChannelKeys(this, compound, own);
    requested += own.size();
    std::sort(own.begin(), own.end());
    own.erase(std::unique(own.begin(), own.end()), own.end());
    keys.insert(keys.end(), own.begin(), own.end());
  }
  std::shared_ptr<const ChannelTable> table =
      ChannelTable::Build(compound, configure, std::move(keys), channelTable_.get());
  if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kProfile))
    configure.outStream << "Channel table: " << requested << " rows asked for, " << table->Size() << " distinct, "
                        << table->Computed() << " computed" << std::endl;
  return table;
}

int WorkingCopy::FillChannelRows(const std::vector<std::vector<Unit *>> &batches, const CompoundSystem *compound,
                                 const Config &configure, const ChannelRowFailure &failed) {
  // Every site of every batch in data order (a point's own site, then the slots it writes), its point, and the site
  // whose rows it holds: the first site at its ChannelRowsKey, which computes them.
  struct Entry {
    Site *site;
    Unit *point;
    int batch;
    int rowsOf;  // index of the site that computes the rows, itself for that site
  };
  std::vector<Entry> sites;
  std::map<std::tuple<int, std::uint64_t, int>, int> first;
  std::vector<Site *> computing;
  for (int b = 0; b < (int)batches.size(); b++)
    for (Unit *point : batches[b]) {
      auto add = [&](Site *site) {
        const int k = (int)sites.size();
        const auto at = first.emplace(site->ChannelRowsKey(compound), k);
        sites.push_back({site, point, b, at.first->second});
        if (at.second) computing.push_back(site);
      };
      add(point);
      for (int i = 1; i <= point->NumSubPoints(); i++)
        if (WritesSlot(*point, i)) add(point->GetSubPoint(*this, i));
    }
  if (!computing.empty() || !channelTable_) {
    try {
      channelTable_ = BuildChannelTable(compound, configure, computing);
    } catch (std::exception &e) {
      configure.outStream << "Channel functions: " << e.what() << std::endl;
      return -1;
    }
  }
  ChannelRowSource rows;
  rows.tables.push_back(channelTable_.get());
  // Batch by batch, as the points were filled one data set at a time: first the sites that compute their rows, then
  // the sites that take a table computed in this batch or an earlier one.
  std::vector<int> computes, takes;
  std::size_t at = 0;
  for (int b = 0; b < (int)batches.size(); b++) {
    computes.clear(), takes.clear();
    for (; at < sites.size() && sites[at].batch == b; at++)
      (sites[at].rowsOf == (int)at ? computes : takes).push_back((int)at);
    bool stop = false;
    for (const std::vector<int> *part : {&computes, &takes}) {
      const int n = (int)part->size();
#pragma omp parallel for schedule(dynamic, 16) shared(stop, configure)
      for (int q = 0; q < n; q++) {
        const azure2::support::ExtendedPrecisionScope fpPrecision;
        if (configure.stopFlag || stop) continue;
        const Entry &e = sites[(*part)[q]];
        try {
          e.site->Site::CalcEDependentValues(this, compound, configure, &rows,
                                             e.rowsOf == (*part)[q] ? nullptr : sites[e.rowsOf].site);
        } catch (GSLException err) {
#pragma omp critical
          {
            failed(*e.point, err);
            stop = true;
          }
        }
      }
      if (configure.stopFlag || stop) return -1;
    }
  }
  return 0;
}

int WorkingCopy::CalcEDependentValues(const CompoundSystem *compound, const Config &configure) {
  std::vector<std::vector<Unit *>> batches;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    batches.emplace_back();
    for (int i = 1; i <= segment->NumPoints(); i++) batches.back().push_back(segment->GetPoint(*this, i));
  }
  return FillChannelRows(batches, compound, configure, [&configure](const Unit &point, const GSLException &e) {
    std::cout << point.GetLabEnergy() << "\t" << point.GetLabAngle() << std::endl;
    configure.outStream << e.what() << std::endl;
  });
}

void WorkingCopy::ShareChannelRows(const std::vector<Site *> &sites, const CompoundSystem *compound) {
  std::vector<std::pair<std::tuple<int, std::uint64_t, int>, int>> keys;
  keys.reserve(sites.size());
  for (int k = 0; k < (int)sites.size(); k++) keys.push_back({sites[k]->ChannelRowsKey(compound), k});
  std::sort(keys.begin(), keys.end());  // by key, then by position: the first of each key holds its table
  for (std::size_t k = 1, first = 0; k < keys.size(); k++) {
    if (keys[k].first != keys[first].first)
      first = k;
    else
      sites[keys[k].second]->TakeChannelRowsOf(*sites[keys[first].second]);
  }
}
