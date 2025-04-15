#pragma once

#include "pcst_fast/pcst_types.h"
#include "pcst_fast/statistics.h"
#include "pcst_fast/logger.h"

#include <vector>
#include <utility>
#include <span>
#include <optional>

#include "pcst_fast/pcst_core_internals.h"

namespace cluster_approx {
struct Cluster;
}

namespace cluster_approx {

/**
 * @brief Structure storing details about a merge event involving one active and one inactive cluster.
 * Used during GW pruning to decide whether to keep the corresponding edge.
 * Definition moved here from pcst_core_internals.h to allow use in CoreAlgorithmResult.
 */
struct InactiveMergeEvent {
    ClusterId active_cluster_index = kInvalidClusterId;
    ClusterId inactive_cluster_index = kInvalidClusterId;
    NodeId active_cluster_node = kInvalidNodeId;
    NodeId inactive_cluster_node = kInvalidNodeId;
};

/**
 * @brief Structure holding the input graph data.
 * Uses std::span for non-owning views of the data provided by the caller.
 */
struct GraphData {
    std::span<const std::pair<NodeId, NodeId>> edges;
    std::span<const double> prizes;
    std::span<const double> costs;
    NodeId root = kInvalidNodeId;
};

/**
 * @brief Structure holding the intermediate result from the core Goemans-Williamson algorithm phase.
 */
struct CoreAlgorithmResult {

    std::vector<EdgeId> phase1_edges;

    std::vector<bool> initial_node_filter;

    std::vector<EventId> edge_inactive_merge_event_ids;

    std::vector<InactiveMergeEvent> inactive_merge_events;

    std::vector<Cluster> final_cluster_state;

    Statistics statistics;
};

/**
 * @brief Structure holding the input required by a pruning strategy.
 */
struct PruningInput {
    const GraphData& graph;
    const CoreAlgorithmResult& core_result;
    Logger* logger;
};

/**
 * @brief Structure holding the final result of the PCST algorithm after pruning.
 */
struct PruningResult {
    std::vector<NodeId> nodes;
    std::vector<EdgeId> edges;
};

/**
 * @brief Interface for different pruning strategies applied after the core PCST algorithm.
 */
class IPruner {
  public:
    virtual ~IPruner() = default;

    /**
     * @brief Applies the specific pruning strategy to the intermediate result.
     * @param input Data structure containing the graph and the result from the core algorithm.
     * @return A PruningResult containing the final set of nodes and edges.
     * @throws std::exception or derived classes on error.
     */
    [[nodiscard]] virtual PruningResult prune(const PruningInput& input) = 0;
};

}