#include "pcst_fast/pruning/no_pruner.h"

#include <vector>
#include <cassert>

namespace cluster_approx {
namespace pruning {

PruningResult NoPruner::prune(const PruningInput& input) {
    assert(input.logger != nullptr);
    input.logger->log(LogLevel::INFO, "Applying NoPruning strategy.");

    PruningResult result;

    result.edges = input.core_result.phase1_edges;
    input.logger->log(LogLevel::DEBUG, "NoPruning: Returning {} phase1 edges.", result.edges.size());

    size_t num_nodes = input.graph.prizes.size();
    std::vector<bool> included(num_nodes, false);
    result.nodes.reserve(num_nodes);

    for (EdgeId edge_idx : result.edges) {
        assert(static_cast<size_t>(edge_idx) < input.graph.edges.size());
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;
        assert(static_cast<size_t>(u) < num_nodes && static_cast<size_t>(v) < num_nodes);

        if (!included[u]) {
            included[u] = true;
            result.nodes.push_back(u);
        }
        if (!included[v]) {
            included[v] = true;
            result.nodes.push_back(v);
        }
    }

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes); ++i) {
        assert(static_cast<size_t>(i) < input.core_result.initial_node_filter.size());
        if (input.core_result.initial_node_filter[i] && !included[i]) {
            result.nodes.push_back(i);
        }
    }

    input.logger->log(LogLevel::DEBUG, "NoPruning: Derived {} nodes using original logic.", result.nodes.size());

    return result;
}

}
}