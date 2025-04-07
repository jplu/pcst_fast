#pragma once

#include "ipruner.h"
#include <vector>
#include <utility>
#include <set>
#include <queue>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

        class ConnectFinalPruner : public IPruner {
            private:
                struct DijkstraState {
                    PCSTFast::ValueType distance;
                    PCSTFast::IndexType node;
                    PCSTFast::IndexType edge_to_parent;

                    bool operator>(const DijkstraState& other) const {
                        return distance > other.distance;
                    }
                };

                void build_full_adjacency(const PruningContext& context,
                                          std::vector<std::vector<std::pair<PCSTFast::IndexType,
                                          PCSTFast::ValueType>>>& adj) const;

                void run_steiner_approximation(const PruningContext& context,
                                               const std::set<PCSTFast::IndexType>& target_nodes,
                                               const std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>>& adj,
                                               std::set<PCSTFast::IndexType>& steiner_nodes_out,
                                               std::set<PCSTFast::IndexType>& steiner_edges_out
                                               );

            public:
                void prune(
                    const PruningContext& context,
                    std::vector<PCSTFast::IndexType>& result_nodes,
                    std::vector<PCSTFast::IndexType>& result_edges
                ) override;
        };

    }
}
