#include "pcst_fast/pruning/simple_pruner.h"
#include "pcst_fast/pruning/pruning_utils.h"

#include <vector>
#include <cassert>

namespace cluster_approx {
namespace pruning {

PruningResult SimplePruner::prune(const PruningInput& input) {
    assert(input.logger != nullptr);
    input.logger->log(LogLevel::INFO, "Applying SimplePruning strategy.");

    PruningResult result;
    size_t num_nodes = input.graph.prizes.size();

    result.edges.reserve(input.core_result.phase1_edges.size());
    for (EdgeId edge_idx : input.core_result.phase1_edges) {
        assert(static_cast<size_t>(edge_idx) < input.graph.edges.size());
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;
        assert(static_cast<size_t>(u) < num_nodes && static_cast<size_t>(v) < num_nodes);
        assert(static_cast<size_t>(u) < input.core_result.initial_node_filter.size());
        assert(static_cast<size_t>(v) < input.core_result.initial_node_filter.size());

        if (input.core_result.initial_node_filter[u] && input.core_result.initial_node_filter[v]) {
            result.edges.push_back(edge_idx);
        }
    }
    input.logger->log(LogLevel::DEBUG, "SimplePruning: Filtered phase1 edges down to {} intermediate edges.", result.edges.size());

    std::vector<bool> node_deleted(num_nodes, false);
    result.nodes = build_final_node_set(num_nodes, node_deleted, input.core_result.initial_node_filter);

    input.logger->log(LogLevel::DEBUG, "SimplePruning: Derived {} nodes.", result.nodes.size());

    return result;
}

}
}