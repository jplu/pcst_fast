#pragma once

#include "ipruner.h"
#include <vector>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

        class NoPruner : public IPruner {
            public:
                void prune(const PruningContext& context,
                        std::vector<PCSTFast::IndexType>& result_nodes,
                        std::vector<PCSTFast::IndexType>& result_edges) override;
        };

    }
}
