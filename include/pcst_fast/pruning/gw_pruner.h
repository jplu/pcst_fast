#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include <vector>
#include <utility>

namespace cluster_approx {

struct Cluster;

namespace pruning {

/**
 * @brief Implements the Goemans-Williamson style pruning strategy.
 *
 * This pruner examines the merge events from the core algorithm run. It traverses
 * the intermediate edge set in reverse order of selection. Edges corresponding to
 * active-active merges are kept. Edges from active-inactive merges are kept only
 * if the specific *intermediate cluster* corresponding to the inactive side
 * (at the time of the merge) is later deemed necessary for connectivity.
 * Nodes connected only through discarded active-inactive merge edges are removed.
 * NOTE: This implementation modifies the 'necessary' flags within the cluster
 * state passed from the core algorithm result to replicate original behavior.
 */
class GWPruner final : public IPruner {
  public:
    /**
     * @brief Applies the GW pruning strategy.
     * @param input Data structure containing the graph, core algorithm result (including merge events and final cluster state), and logger.
     * @return A PruningResult containing the final filtered set of nodes and edges.
     * @throws std::exception or derived classes on error.
     */
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;

  private:

    const PruningInput* input_ = nullptr;
    size_t num_nodes_ = 0;
    Logger* logger_ = nullptr;

    std::vector<bool> node_deleted_;

    std::vector<Cluster>* clusters_ptr_ = nullptr;

    std::vector<std::vector<std::pair<NodeId, double>>> neighbors_;

    std::vector<NodeId> node_queue_;

    /**
     * @brief Marks a cluster (corresponding to an original node) and its ancestors
     * in the merge tree as necessary by setting the 'necessary' flag on the Cluster structs.
     * Operates on the clusters_ptr_ vector, modifying the state passed in PruningInput.
     * @param start_node_index The original node ID whose corresponding cluster needs marking.
     */
    void mark_clusters_as_necessary_from_node(NodeId start_node_index);

    /**
     * @brief Marks a node and its subtree (excluding the parent) as deleted.
     * Uses BFS traversal on the `neighbors_` graph.
     * @param start_node_index The node to start the deletion marking from.
     * @param parent_node_index The node from which start_node_index was reached (to avoid traversing back up).
     */
    void mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index);

};

}
}