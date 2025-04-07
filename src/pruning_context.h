#pragma once

#include <vector>
#include <utility>
#include <memory>
#include <string>
#include <cstdarg>

#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        struct PruningContext {
            const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges;
            const std::vector<PCSTFast::ValueType>& prizes;
            const std::vector<PCSTFast::ValueType>& costs;
            const PCSTFast::IndexType root;
            const std::vector<PCSTFast::Cluster>& clusters;
            const std::vector<PCSTFast::EdgeInfo>& edge_info;
            const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events;
            const std::vector<uint8_t>& node_good;
            const std::vector<int>& phase1_result;
            const PCSTFast::Logger& logger;
            int verbosity_level;

            PruningContext(
                const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges_ref,
                const std::vector<PCSTFast::ValueType>& prizes_ref,
                const std::vector<PCSTFast::ValueType>& costs_ref,
                PCSTFast::IndexType root_val,
                const std::vector<PCSTFast::Cluster>& clusters_ref,
                const std::vector<PCSTFast::EdgeInfo>& edge_info_ref,
                const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events_ref,
                const std::vector<uint8_t>& node_good_ref,
                const std::vector<int>& phase1_result_ref,
                const PCSTFast::Logger& logger_ref,
                int verbosity
            ) : edges(edges_ref),
                prizes(prizes_ref),
                costs(costs_ref),
                root(root_val),
                clusters(clusters_ref),
                edge_info(edge_info_ref),
                inactive_merge_events(inactive_merge_events_ref),
                node_good(node_good_ref),
                phase1_result(phase1_result_ref),
                logger(logger_ref),
                verbosity_level(verbosity)
                {}

            void log(int level, const char* format, ...) const
                #if defined(__GNUC__) || defined(__clang__)
                    __attribute__ ((format (printf, 3, 4)))
                #endif
                ;
        };

    }
}
