#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include <vector>
#include <utility>

namespace cluster_approx {
namespace pruning {

/**
 * @brief Implements the "strong pruning" strategy.
 *
 * This method identifies connected components in the graph formed by the intermediate
 * edges. For each component, it finds an optimal root (maximizing prize minus cost)
 * and performs a traversal to prune subtrees whose total prize does not cover the
 * cost of the edge connecting them to the rest of the component.
 */
class StrongPruner final : public IPruner {
  public:
    /**
     * @brief Applies the strong pruning strategy.
     * @param input Data structure containing the graph, core algorithm result, and logger.
     * @return A PruningResult containing the final filtered set of nodes and edges.
     * @throws std::exception or derived classes on error.
     */
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;

  private:

    const PruningInput* input_ = nullptr;
    size_t num_nodes_ = 0;
    Logger* logger_ = nullptr;

    std::vector<std::vector<std::pair<NodeId, double>>> neighbors_;

    std::vector<bool> node_deleted_;
    std::vector<ClusterId> final_component_label_;
    std::vector<std::vector<NodeId>> final_components_;
    ClusterId root_component_index_ = kInvalidClusterId;

    std::vector<std::pair<NodeId, double>> strong_pruning_parent_;
    std::vector<double> strong_pruning_payoff_;

    std::vector<std::pair<bool, NodeId>> dfs_stack_;
    std::vector<NodeId> dfs_stack2_;

    /**
     * @brief Performs a Depth First Search (DFS) starting from a node to label a connected component.
     * Populates `final_component_label_` and `final_components_`. Identifies `root_component_index_`.
     * @param start_node_index The node to start the DFS from.
     * @param component_index The integer label to assign to this component.
     */
    void label_final_component(NodeId start_node_index, ClusterId component_index);

    /**
     * @brief Performs the core strong pruning logic using DFS starting from a given node.
     * Calculates subtree payoffs and marks subtrees with non-positive contributions as deleted.
     * Uses iterative DFS with `dfs_stack_`.
     * @param start_node_index The node to start the DFS from (component root).
     * @param mark_as_deleted If true, actually marks nodes in `node_deleted_`; if false, only calculates payoffs.
     */
    void strong_pruning_dfs(NodeId start_node_index, bool mark_as_deleted);

    /**
     * @brief Finds the best root node within a given component that maximizes the pruned subtree value.
     * Uses results from a previous `strong_pruning_dfs` call (with mark_as_deleted=false) and performs additional traversals.
     * Uses iterative DFS with `dfs_stack2_`.
     * @param component_index The index of the component to analyze.
     * @return The NodeId of the best root found for this component.
     */
    [[nodiscard]] NodeId find_best_component_root(ClusterId component_index);

    /**
     * @brief Marks a node and its subtree (excluding the parent) as deleted. (Identical to GWPruner's version).
     * Uses BFS traversal on the `neighbors_` graph.
     * @param start_node_index The node to start the deletion marking from.
     * @param parent_node_index The node from which start_node_index was reached (to avoid traversing back up).
     */
    void mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index);

    std::vector<NodeId> node_queue_;

};

}
}