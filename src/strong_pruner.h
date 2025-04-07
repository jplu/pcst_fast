#pragma once

#include "advanced_pruner_base.h"
#include <vector>
#include <utility>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

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

    }
}
