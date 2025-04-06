#pragma once

#include <vector>
#include <utility>
#include <memory>
#include <queue>
#include <set>

#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {
        struct PruningContext {
            const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges;
            const std::vector<PCSTFast::ValueType>& prizes;
            const std::vector<PCSTFast::ValueType>& costs;
            const PCSTFast::IndexType root;
            const std::vector<PCSTFast::Cluster>& clusters;
            const std::vector<PCSTFast::EdgeInfo>& edge_info;
            const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events;
            const std::vector<uint8_t>& node_good;
            const std::vector<int>& phase1_result;
            const PCSTFast::Logger& logger;
            int verbosity_level;

            PruningContext(
                const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges_ref,
                const std::vector<PCSTFast::ValueType>& prizes_ref,
                const std::vector<PCSTFast::ValueType>& costs_ref,
                PCSTFast::IndexType root_val,
                const std::vector<PCSTFast::Cluster>& clusters_ref,
                const std::vector<PCSTFast::EdgeInfo>& edge_info_ref,
                const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events_ref,
                const std::vector<uint8_t>& node_good_ref,
                const std::vector<int>& phase1_result_ref,
                const PCSTFast::Logger& logger_ref,
                int verbosity
            ) : edges(edges_ref),
                prizes(prizes_ref),
                costs(costs_ref),
                root(root_val),
                clusters(clusters_ref),
                edge_info(edge_info_ref),
                inactive_merge_events(inactive_merge_events_ref),
                node_good(node_good_ref),
                phase1_result(phase1_result_ref),
                logger(logger_ref),
                verbosity_level(verbosity)
                {}

            void log(int level, const char* format, ...) const
                #if defined(__GNUC__) || defined(__clang__)
                    __attribute__ ((format (printf, 3, 4)))
                #endif
                ;
        };

        class IPruner {
            public:
                virtual ~IPruner() = default;
                virtual void prune(
                    const PruningContext& context,
                    std::vector<PCSTFast::IndexType>& result_nodes,
                    std::vector<PCSTFast::IndexType>& result_edges
                ) = 0;

            protected:
                void build_node_set_base(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set);
            };

        class AdvancedPrunerBase : public IPruner {
            protected:
                std::vector<PCSTFast::IndexType> phase2_result_local_;
                std::vector<PCSTFast::IndexType> phase3_result_local_;
                std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>> phase3_neighbors_;
                std::vector<uint8_t> node_deleted_;

                void build_phase2(const PruningContext& context);
                void build_phase3_adjacency(const PruningContext& context);
                void build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set);

            public:
                virtual ~AdvancedPrunerBase() = default;
                virtual void setup(const PruningContext& context);
        };


        class NoPruner : public IPruner {
            public:
                void prune(const PruningContext& context,
                        std::vector<PCSTFast::IndexType>& result_nodes,
                        std::vector<PCSTFast::IndexType>& result_edges) override;
            };

        class SimplePruner : public IPruner {
            private:
                std::vector<PCSTFast::IndexType> phase2_result_local_;
                void build_phase2(const PruningContext& context);
            public:
                void prune(const PruningContext& context,
                        std::vector<PCSTFast::IndexType>& result_nodes,
                        std::vector<PCSTFast::IndexType>& result_edges) override;
            };

        class GWPruner : public AdvancedPrunerBase {
            private:
                std::vector<bool> cluster_necessary_local_;

                void mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index);
                void mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index);
                void run_gw_pruning(const PruningContext& context);
            public:
                void prune(const PruningContext& context,
                        std::vector<PCSTFast::IndexType>& result_nodes,
                        std::vector<PCSTFast::IndexType>& result_edges) override;
            };

        class StrongPruner : public AdvancedPrunerBase {
            private:
                std::vector<PCSTFast::IndexType> final_component_label_;
                std::vector<std::vector<PCSTFast::IndexType>> final_components_;
                PCSTFast::IndexType root_component_index_ = PCSTFast::kInvalidIndex;
                std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>> strong_pruning_parent_;
                std::vector<PCSTFast::ValueType> strong_pruning_payoff_;

                void label_final_components(const PruningContext& context);
                void label_component_bfs(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id);
                PCSTFast::IndexType find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx);
                void propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root, PCSTFast::ValueType& best_value);
                void strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted);
                void mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node);
                void run_strong_pruning(const PruningContext& context);
            public:
                void prune(const PruningContext& context,
                           std::vector<PCSTFast::IndexType>& result_nodes,
                           std::vector<PCSTFast::IndexType>& result_edges) override;
                void setup(const PruningContext& context) override;
            };

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
                                               std::set<PCSTFast::IndexType>& steiner_edges_out);

            public:
                void prune(
                    const PruningContext& context,
                    std::vector<PCSTFast::IndexType>& result_nodes,
                    std::vector<PCSTFast::IndexType>& result_edges
                ) override;
            };

        std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method);

    }
}
