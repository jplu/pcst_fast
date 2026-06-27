#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"

#include <vector>
#include <utility>
#include <cstdint>

namespace cluster_approx {
namespace pruning {

struct CSRGraph {
    std::vector<uint32_t> row_ptr;
    std::vector<NodeId> col_indices;
    std::vector<double> edge_costs;

    struct NeighborsRange {
        struct Iterator {
            const NodeId* col_ptr;
            const double* cost_ptr;

            bool operator!=(const Iterator& other) const noexcept { return col_ptr != other.col_ptr; }
            void operator++() noexcept { ++col_ptr; ++cost_ptr; }
            std::pair<NodeId, double> operator*() const noexcept { return {*col_ptr, *cost_ptr}; }
        };

        Iterator begin_it;
        Iterator end_it;

        Iterator begin() const noexcept { return begin_it; }
        Iterator end() const noexcept { return end_it; }
    };

    [[nodiscard]] NeighborsRange get_neighbors(NodeId node) const noexcept {
        uint32_t start = row_ptr[node];
        uint32_t end = row_ptr[node + 1];
        return NeighborsRange{
            .begin_it = {&col_indices[start], &edge_costs[start]},
            .end_it = {&col_indices[end], &edge_costs[end]}
        };
    }
};

[[nodiscard]] std::vector<NodeId> build_final_node_set(
    size_t num_nodes,
    const std::vector<uint8_t>& node_deleted_filter,
    const std::vector<uint8_t>& initial_node_filter);

[[nodiscard]] CSRGraph build_adjacency_list_csr(
    const std::vector<EdgeId>& edges,
    const GraphData& graph);

}
}