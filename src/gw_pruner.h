#pragma once

#include "advanced_pruner_base.h"
#include <vector>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

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

    }
}
