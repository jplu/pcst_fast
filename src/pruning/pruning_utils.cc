#include "pcst_fast/pruning/pruning_utils.h"
#include "pcst_fast/logger.h"

#include <vector>
#include <algorithm>
#include <cassert>

namespace cluster_approx {
namespace pruning {

std::vector<NodeId> build_final_node_set(
    size_t num_nodes,
    const std::vector<bool>& node_deleted_filter,
    const std::vector<bool>& initial_node_filter) {
    assert(node_deleted_filter.size() == num_nodes);
    assert(initial_node_filter.size() == num_nodes);

    std::vector<NodeId> final_nodes;
    final_nodes.reserve(num_nodes);

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes); ++i) {

        if (initial_node_filter[i] && !node_deleted_filter[i]) {
            final_nodes.push_back(i);
        }
    }

    return final_nodes;
}

std::vector<std::vector<std::pair<NodeId, double>>> build_adjacency_list(
    const std::vector<EdgeId>& edges,
    const GraphData& graph) {
    NodeId max_node_id = -1;
    for (EdgeId edge_idx : edges) {
        assert(static_cast<size_t>(edge_idx) < graph.edges.size());
        const auto& edge = graph.edges[edge_idx];
        max_node_id = std::max({max_node_id, edge.first, edge.second});
    }

    size_t adj_list_size = (max_node_id == -1) ? 0 : static_cast<size_t>(max_node_id + 1);

    adj_list_size = std::max(adj_list_size, graph.prizes.size());

    std::vector<std::vector<std::pair<NodeId, double>>> adj_list(adj_list_size);

    for (EdgeId edge_idx : edges) {
        assert(static_cast<size_t>(edge_idx) < graph.edges.size());
        assert(static_cast<size_t>(edge_idx) < graph.costs.size());
        const auto& edge = graph.edges[edge_idx];
        double cost = graph.costs[edge_idx];
        NodeId u = edge.first;
        NodeId v = edge.second;

        assert(static_cast<size_t>(u) < adj_list.size());
        assert(static_cast<size_t>(v) < adj_list.size());

        adj_list[u].emplace_back(v, cost);
        adj_list[v].emplace_back(u, cost);
    }

    return adj_list;
}

}
}