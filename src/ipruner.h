#pragma once

#include <vector>
#include <memory>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext;

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

    }
}
