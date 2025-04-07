#pragma once

#include "ipruner.h"
#include <vector>
#include <utility>
#include <set>
#include <queue>
#include <tuple>
#include "pcst_fast.h"
#include "logger.h"

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

                using AdjacencyList = std::vector<std::vector<std::tuple<PCSTFast::IndexType,
                                                                        PCSTFast::ValueType,
                                                                        PCSTFast::IndexType>>>;

                void build_full_adjacency(const PruningContext& context,
                                          AdjacencyList& adj) const;

                void run_steiner_approximation(const PruningContext& context,
                                               const std::set<PCSTFast::IndexType>& target_nodes,
                                               const AdjacencyList& adj,
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
