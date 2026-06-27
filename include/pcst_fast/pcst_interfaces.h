#pragma once

#include "pcst_fast/pcst_types.h"
#include "pcst_fast/statistics.h"
#include "pcst_fast/logger.h"
#include "pcst_fast/pcst_core_internals.h"

#include <vector>
#include <utility>
#include <span>
#include <optional>
#include <memory>

namespace cluster_approx {

struct InactiveMergeEvent {
    ClusterId active_cluster_index = kInvalidClusterId;
    ClusterId inactive_cluster_index = kInvalidClusterId;
    NodeId active_cluster_node = kInvalidNodeId;
    NodeId inactive_cluster_node = kInvalidNodeId;
};

struct GraphData {
    std::span<const std::pair<NodeId, NodeId>> edges;
    std::span<const double> prizes;
    std::span<const double> costs;
    NodeId root = kInvalidNodeId;
};

struct CoreAlgorithmResult {
    std::vector<EdgeId> phase1_edges;
    std::vector<uint8_t> initial_node_filter;
    std::vector<EventId> edge_inactive_merge_event_ids;
    std::vector<InactiveMergeEvent> inactive_merge_events;
    std::vector<Cluster> final_cluster_state;
    Statistics statistics;
};

struct PruningInput {
    const GraphData& graph;
    CoreAlgorithmResult& core_result;
    Logger* logger;
};

struct PruningResult {
    std::vector<NodeId> nodes;
    std::vector<EdgeId> edges;
};

class IPruner {
  public:
    virtual ~IPruner() = default;
    [[nodiscard]] virtual PruningResult prune(const PruningInput& input) = 0;
};

}