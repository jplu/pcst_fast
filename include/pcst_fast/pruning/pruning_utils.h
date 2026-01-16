#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"

#include <vector>
#include <utility>

namespace cluster_approx {
namespace pruning {

[[nodiscard]] std::vector<NodeId> build_final_node_set(
    size_t num_nodes,
    const std::vector<bool>& node_deleted_filter,
    const std::vector<bool>& initial_node_filter);

[[nodiscard]] std::vector<std::vector<std::pair<NodeId, double>>> build_adjacency_list(
    const std::vector<EdgeId>& edges,
    const GraphData& graph);

}
}