#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"

#include <vector>
#include <utility>

namespace cluster_approx {
namespace pruning {

/**
 * @brief Builds the final set of nodes based on filtering flags.
 * Includes nodes that passed the initial filter and were not subsequently deleted by pruning.
 * @param num_nodes The total number of nodes in the original graph.
 * @param node_deleted_filter A vector indicating if a node was marked for deletion by the pruning process.
 * @param initial_node_filter A vector indicating if a node passed the core algorithm's initial filtering (e.g., node_good).
 * @return A vector containing the NodeIds of the final selected nodes.
 */
[[nodiscard]] std::vector<NodeId> build_final_node_set(
    size_t num_nodes,
    const std::vector<bool>& node_deleted_filter,
    const std::vector<bool>& initial_node_filter);


/**
 * @brief Builds an adjacency list representation of a graph subset defined by selected edges.
 * @param edges A vector of EdgeIds representing the subgraph.
 * @param graph The original graph data containing edge endpoints and costs.
 * @return An adjacency list where each element is a vector of pairs {neighbor_node_id, edge_cost}.
 * The size of the outer vector is determined by the highest node index encountered.
 */
[[nodiscard]] std::vector<std::vector<std::pair<NodeId, double>>> build_adjacency_list(
    const std::vector<EdgeId>& edges,
    const GraphData& graph);


}
}