#pragma once

#include <vector>
#include <utility>
#include <functional>
#include <memory> // For unique_ptr
#include <string> // For context logger helper
#include <cstdio> // For context logger helper
#include <cstdarg> // For context logger helper
#include <stdexcept> // For exceptions

#include "pcst_fast.h" // Include full header for enum, types, structs

namespace cluster_approx {
namespace internal {

    // Context structure to pass necessary data from PCSTFast to the Pruner
    struct PruningContext {
        const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges;
        const std::vector<PCSTFast::ValueType>& prizes;
        const std::vector<PCSTFast::ValueType>& costs;
        const PCSTFast::IndexType root;
        const std::vector<PCSTFast::Cluster>& clusters;
        const std::vector<PCSTFast::EdgeInfo>& edge_info;
        const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events;
        const std::vector<uint8_t>& node_good;
        const std::vector<int>& phase1_result; // phase1 is edge indices (int)

        const PCSTFast::Logger& logger;
        int verbosity_level;

        // Constructor
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

        // Logger helper using C++ parameter pack directly with snprintf
         template<typename... Args>
         void log(int level, const char* format, Args... args) const {
             if (level <= verbosity_level && logger) {
                 // ***** CHANGED HERE: Handle zero-argument case *****
                 if constexpr (sizeof...(args) == 0) {
                     // If there are no arguments to format, pass the format string directly
                     logger(level, std::string(format));
                 } else {
                     // If there ARE arguments, format them using snprintf
                     constexpr int BUFFER_SIZE = 1024;
                     char buffer[BUFFER_SIZE];
                     int needed = snprintf(buffer, BUFFER_SIZE, format, args...);

                     if (needed < 0) {
                         logger(level, "[Pruner Logging Error: snprintf failed]");
                     } else if (static_cast<size_t>(needed) >= BUFFER_SIZE) {
                          buffer[BUFFER_SIZE - 1] = '\0';
                          if (BUFFER_SIZE > 4) {
                               snprintf(buffer + BUFFER_SIZE - 5, 5, "...");
                          }
                          logger(level, std::string(buffer) + "(TRUNCATED)");
                     } else {
                         logger(level, std::string(buffer));
                     }
                 }
                 // ***** END CHANGE *****
             }
         }
    };

    // Interface for different pruning strategies
    class IPruner {
    public:
        virtual ~IPruner() = default;
        virtual void prune(
            const PruningContext& context,
            std::vector<PCSTFast::IndexType>& result_nodes, // Use IndexType
            std::vector<PCSTFast::IndexType>& result_edges  // Use IndexType (edge indices are int)
        ) = 0;
    };


    // Factory function declaration
    std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method);

} // namespace internal
} // namespace cluster_approx
