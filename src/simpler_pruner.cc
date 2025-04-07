#include "simple_pruner.h"
#include "pruning_context.h"

#include <vector>

namespace cluster_approx {
    namespace internal {

        void SimplePruner::build_phase2(const PruningContext& context) {
             context.log(3, "SimplePruner::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
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
                        context.log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                    }
                } else {
                    context.log(2, "Warning: Invalid edge index %d in SimplePruner::build_phase2.\n", edge_idx);
                }
            }
            context.log(3, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
        }

        void SimplePruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.log(3, "Pruning: Simple. Running Phase 2 filtering.\n");
            build_phase2(context);
            result_edges = phase2_result_local_;
            build_node_set_base(context, result_edges, result_nodes);
            context.log(3, "Final Result (Simple Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
