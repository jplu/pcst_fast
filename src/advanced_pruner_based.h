#pragma once

#include "ipruner.h"
#include <vector>
#include <utility>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

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

    }
}
