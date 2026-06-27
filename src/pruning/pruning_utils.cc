#include "pcst_fast/pruning/pruning_utils.h"
#include "pcst_fast/logger.h"

#include <vector>
#include <algorithm>
#include <cassert>

namespace cluster_approx {
namespace pruning {

std::vector<NodeId> build_final_node_set(
    size_t num_nodes,
    const std::vector<uint8_t>& node_deleted_filter,
    const std::vector<uint8_t>& initial_node_filter) {
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

CSRGraph build_adjacency_list_csr(
    const std::vector<EdgeId>& edges,
    const GraphData& graph) {
    
    size_t num_nodes = graph.prizes.size();
    CSRGraph csr;
    csr.row_ptr.assign(num_nodes + 1, 0);

    // Step 1: Count degrees
    std::vector<uint32_t> degrees(num_nodes, 0);
    for (EdgeId edge_idx : edges) {
        const auto& edge = graph.edges[edge_idx];
        degrees[edge.first]++;
        degrees[edge.second]++;
    }

    // Step 2: Establish base row index pointers leveraging degrees array as tracking offsets implicitly
    uint32_t sum = 0;
    for (size_t i = 0; i < num_nodes; ++i) {
        csr.row_ptr[i] = sum;
        uint32_t deg = degrees[i];
        degrees[i] = sum;
        sum += deg;
    }
    csr.row_ptr[num_nodes] = sum;

    csr.col_indices.resize(sum);
    csr.edge_costs.resize(sum);

    // Step 3: Contiguously write neighbor indices directly bypassing vector offset copies
    for (EdgeId edge_idx : edges) {
        const auto& edge = graph.edges[edge_idx];
        double cost = graph.costs[edge_idx];
        NodeId u = edge.first;
        NodeId v = edge.second;

        uint32_t u_pos = degrees[u]++;
        uint32_t v_pos = degrees[v]++;

        csr.col_indices[u_pos] = v;
        csr.edge_costs[u_pos] = cost;

        csr.col_indices[v_pos] = u;
        csr.edge_costs[v_pos] = cost;
    }

    return csr;
}

}
}