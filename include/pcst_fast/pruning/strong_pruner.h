#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pruning/pruning_utils.h"
#include <vector>
#include <utility>

namespace cluster_approx {
namespace pruning {

class StrongPruner final : public IPruner {
  public:
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;

  private:
    const PruningInput* input_ = nullptr;
    size_t num_nodes_ = 0;
    Logger* logger_ = nullptr;

    CSRGraph neighbors_;

    std::vector<bool> node_deleted_;
    std::vector<ClusterId> final_component_label_;
    std::vector<std::vector<NodeId>> final_components_;
    ClusterId root_component_index_ = kInvalidClusterId;

    std::vector<std::pair<NodeId, double>> strong_pruning_parent_;
    std::vector<double> strong_pruning_payoff_;

    std::vector<std::pair<bool, NodeId>> dfs_stack_;
    std::vector<NodeId> dfs_stack2_;
    std::vector<NodeId> node_queue_;

    void label_final_component(NodeId start_node_index, ClusterId component_index);
    void strong_pruning_dfs(NodeId start_node_index, bool mark_as_deleted);
    [[nodiscard]] NodeId find_best_component_root(ClusterId component_index);
    void mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index);
};

}
}