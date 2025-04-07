#include "no_pruner.h"
#include "pruning_context.h"

#include <vector>
#include <numeric>

namespace cluster_approx {
    namespace internal {

        void NoPruner::prune(const PruningContext& context,
                    std::vector<PCSTFast::IndexType>& result_nodes,
                    std::vector<PCSTFast::IndexType>& result_edges) {
            context.log(3, "Pruning: None. Using Phase 1 result directly.\n");
            result_edges.assign(context.phase1_result.begin(), context.phase1_result.end());
            build_node_set_base(context, result_edges, result_nodes);
            context.log(3, "Final Result (No Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
