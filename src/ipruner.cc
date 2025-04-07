#include "ipruner.h"
#include "pruning_context.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <set>
#include <map>
#include <limits>

namespace cluster_approx {
    namespace internal {

        void IPruner::build_node_set_base(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
            context.log(4, "IPruner::build_node_set_base Entry (using %zu edges).\n", edge_set.size());
            node_set.clear();
            size_t num_nodes = context.prizes.size();
            node_set.reserve(num_nodes);

            std::vector<uint8_t> included_nodes_local(num_nodes, 0);

            for (PCSTFast::IndexType edge_idx : edge_set) {
                if (static_cast<size_t>(edge_idx) >= context.edges.size()) {
                    context.log(2, "Warning: Invalid edge index %d in build_node_set_base.\n", edge_idx);
                    continue;
                }
                const auto& edge = context.edges[edge_idx];
                PCSTFast::IndexType uu = edge.first;
                PCSTFast::IndexType vv = edge.second;

                if (static_cast<size_t>(uu) < included_nodes_local.size() && !included_nodes_local[uu]) {
                    included_nodes_local[uu] = 1;
                    node_set.push_back(uu);
                }
                if (static_cast<size_t>(vv) < included_nodes_local.size() && !included_nodes_local[vv]) {
                    included_nodes_local[vv] = 1;
                    node_set.push_back(vv);
                }
            }

            for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
                bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
                bool is_included = (static_cast<size_t>(ii) < included_nodes_local.size() && included_nodes_local[ii]);

                if (is_good && !is_included) {
                    context.log(4, "  Adding isolated good node %d.\n", ii);
                    node_set.push_back(ii);
                }
            }

            std::sort(node_set.begin(), node_set.end());
            context.log(4, "IPruner::build_node_set_base Exit. Final node set size: %zu\n", node_set.size());
        }

    }
}
