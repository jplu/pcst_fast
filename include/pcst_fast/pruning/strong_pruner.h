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

    std::vector<uint8_t> node_deleted_;
    std::vector<ClusterId> final_component_label_;
    std::vector<std::vector<NodeId>> final_components_;
    ClusterId root_component_index_ = kInvalidClusterId;

    std::vector<std::pair<NodeId, double>> strong_pruning_parent_;
    std::vector<double> strong_pruning_payoff_;

    void label_final_component(NodeId start_node_index, ClusterId component_index, std::vector<NodeId>& dfs_stack2);
    
    // Adapted specifically for OpenMP threaded stack contexts
    void strong_pruning_dfs(NodeId start_node_index, bool mark_as_deleted,
                            std::vector<std::pair<bool, NodeId>>& local_dfs_stack,
                            std::vector<NodeId>& local_node_queue);
    
    [[nodiscard]] NodeId find_best_component_root(ClusterId component_index,
                                                  std::vector<std::pair<bool, NodeId>>& local_dfs_stack,
                                                  std::vector<NodeId>& local_dfs_stack2);
                                                  
    void mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index,
                               std::vector<NodeId>& local_node_queue);
};

}
}