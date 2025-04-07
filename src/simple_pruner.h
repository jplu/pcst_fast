#pragma once

#include "ipruner.h"
#include <vector>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

        class SimplePruner : public IPruner {
            private:
                std::vector<PCSTFast::IndexType> phase2_result_local_;
                void build_phase2(const PruningContext& context);
            public:
                void prune(const PruningContext& context,
                        std::vector<PCSTFast::IndexType>& result_nodes,
                        std::vector<PCSTFast::IndexType>& result_edges) override;
        };

    }
}
