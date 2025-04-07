#include "advanced_pruner_base.h"
#include "pruning_context.h"
#include "logger.h"

#include <vector>
#include <algorithm>
#include <utility>
#include <cassert>

namespace cluster_approx {
    namespace internal {

        using internal::LogLevel;

        void AdvancedPrunerBase::build_phase2(const PruningContext& context) {
            context.logger.log(LogLevel::INFO, "AdvancedPrunerBase::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
            phase2_result_local_.clear();
            phase2_result_local_.reserve(context.phase1_result.size());

            for (int edge_idx_int : context.phase1_result) {
                PCSTFast::IndexType edge_idx = static_cast<PCSTFast::IndexType>(edge_idx_int);

                if (static_cast<size_t>(edge_idx) < context.edges.size()) {
                    const auto& edge = context.edges[edge_idx];
                    bool u_good = static_cast<size_t>(edge.first) < context.node_good.size() && context.node_good[edge.first];
                    bool v_good = static_cast<size_t>(edge.second) < context.node_good.size() && context.node_good[edge.second];

                    if (u_good && v_good) {
                        phase2_result_local_.push_back(edge_idx);
                    } else {
                        context.logger.log(LogLevel::DEBUG, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                    }
                } else {
                     context.logger.log(LogLevel::WARNING, "Warning: Invalid edge index %d in AdvancedPrunerBase::build_phase2.\n", edge_idx);
                }
            }
            context.logger.log(LogLevel::INFO, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
        }

        void AdvancedPrunerBase::build_phase3_adjacency(const PruningContext& context) {
            context.logger.log(LogLevel::INFO, "AdvancedPrunerBase::build_phase3_adjacency Entry (Using %zu phase 2 edges).\n", phase2_result_local_.size());
            size_t num_nodes = context.prizes.size();

            if (phase3_neighbors_.size() != num_nodes) {
                phase3_neighbors_.resize(num_nodes);
            }
            for (auto& neighbors : phase3_neighbors_) { neighbors.clear(); }

            int edges_added = 0;
            for (PCSTFast::IndexType edge_idx : phase2_result_local_) {
                if (static_cast<size_t>(edge_idx) < context.edges.size() && static_cast<size_t>(edge_idx) < context.costs.size()) {
                    const auto& edge = context.edges[edge_idx];
                    PCSTFast::ValueType cost = context.costs[edge_idx];

                    if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                        phase3_neighbors_[edge.first].emplace_back(edge.second, cost);
                        phase3_neighbors_[edge.second].emplace_back(edge.first, cost);
                        edges_added++;
                    } else {
                        context.logger.log(LogLevel::FATAL, "Error: Invalid node index in edge %d while building phase 3 adjacency list.\n", edge_idx);
                        assert(false && "Invalid node index in build_phase3_adjacency");
                    }
                } else {
                    context.logger.log(LogLevel::WARNING, "Warning: Invalid edge index %d found in phase2_result_ during adjacency build.\n", edge_idx);
                }
            }
            context.logger.log(LogLevel::INFO, "AdvancedPrunerBase::build_phase3_adjacency Exit. Added %d edges (x2) to lists.\n", edges_added);
        }

        void AdvancedPrunerBase::build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set) {
            context.logger.log(LogLevel::INFO, "AdvancedPrunerBase::build_pruned_node_set Entry.\n");
            node_set.clear();
            size_t num_nodes = context.prizes.size();
            node_set.reserve(num_nodes);

            int nodes_included = 0;
            for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
                bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
                bool is_deleted = (static_cast<size_t>(ii) < node_deleted_.size() && node_deleted_[ii]);

                if (is_good && !is_deleted) {
                    node_set.push_back(ii);
                    nodes_included++;
                }
            }

            std::sort(node_set.begin(), node_set.end());
            context.logger.log(LogLevel::INFO, "AdvancedPrunerBase::build_pruned_node_set Exit. Final node set size: %zu\n", static_cast<size_t>(nodes_included));
        }

        void AdvancedPrunerBase::setup(const PruningContext& context) {
            context.logger.log(LogLevel::DEBUG, "AdvancedPrunerBase::setup Entry.\n");
            phase2_result_local_.reserve(context.phase1_result.size());

            size_t num_nodes = context.prizes.size();
            node_deleted_.assign(num_nodes, 0);
            phase3_neighbors_.resize(num_nodes);

            build_phase2(context);
            build_phase3_adjacency(context);
            context.logger.log(LogLevel::DEBUG, "AdvancedPrunerBase::setup Exit.\n");
        }

    }
}
