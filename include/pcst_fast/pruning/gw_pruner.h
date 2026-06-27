#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pruning/pruning_utils.h"
#include <vector>
#include <utility>
#include <cstdint>

namespace cluster_approx {

struct Cluster;

namespace pruning {

class GWPruner final : public IPruner {
  public:
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;

  private:
    const PruningInput* input_ = nullptr;
    size_t num_nodes_ = 0;
    Logger* logger_ = nullptr;

    std::vector<uint8_t> node_deleted_;
    std::vector<Cluster>* clusters_ptr_ = nullptr;
    CSRGraph neighbors_;
    std::vector<NodeId> node_queue_;

    void mark_clusters_as_necessary_from_node(NodeId start_node_index);
    void mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index);
};

}
}