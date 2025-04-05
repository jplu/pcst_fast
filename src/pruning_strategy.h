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

    // Context structure to pass necessary read-only data and configuration
    // from PCSTFast to the Pruner strategy object.
    struct PruningContext {
        // Input graph data (references to PCSTFast's data)
        const std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>>& edges;
        const std::vector<PCSTFast::ValueType>& prizes;
        const std::vector<PCSTFast::ValueType>& costs;
        const PCSTFast::IndexType root; // kNoRoot if unrooted

        // State from the main algorithm phase
        const std::vector<PCSTFast::Cluster>& clusters;             // Full cluster history
        const std::vector<PCSTFast::EdgeInfo>& edge_info;           // Extra info per edge (e.g., inactive merge event link)
        const std::vector<PCSTFast::InactiveMergeEvent>& inactive_merge_events; // Log of A-I merges
        const std::vector<uint8_t>& node_good;                      // Nodes connected to final active/root clusters
        const std::vector<int>& phase1_result;                      // Edge indices selected by main loop

        // Configuration
        const PCSTFast::Logger& logger;         // Logger function object
        int verbosity_level;                    // Verbosity level for logging

        // Constructor (initializes all references)
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

        // Helper function within the context for easy logging by pruners
         template<typename... Args>
         void log(int level, const char* format, Args... args) const {
             // Check level and if logger is valid before formatting
             if (level <= verbosity_level && logger) {
                 // Handle case where no format arguments are provided
                 if constexpr (sizeof...(args) == 0) {
                     logger(level, std::string(format));
                 } else {
                     // Format with arguments using snprintf
                     constexpr int BUFFER_SIZE = 1024; // Match PCSTFast buffer size
                     char buffer[BUFFER_SIZE];
                     int needed = snprintf(buffer, BUFFER_SIZE, format, args...);

                     if (needed < 0) {
                         // Formatting error
                         logger(level, "[Pruner Logging Error: snprintf failed]");
                     } else if (static_cast<size_t>(needed) >= BUFFER_SIZE) {
                          // Truncation occurred
                          buffer[BUFFER_SIZE - 1] = '\0';
                          if (BUFFER_SIZE > 4) { // Ensure space for ellipsis
                               snprintf(buffer + BUFFER_SIZE - 5, 5, "...");
                          }
                          logger(level, std::string(buffer) + "(TRUNCATED)");
                     } else {
                         // Message fits
                         logger(level, std::string(buffer));
                     }
                 }
             }
         }
    };

    // Interface (Abstract Base Class) for different pruning strategies
    class IPruner {
    public:
        virtual ~IPruner() = default; // Ensure virtual destructor for base class

        /**
         * @brief Applies the specific pruning strategy.
         * @param context Provides read-only access to necessary data from PCSTFast.
         * @param result_nodes Output vector to store the final node indices. Will be cleared by the method.
         * @param result_edges Output vector to store the final edge indices. Will be cleared by the method.
         */
        virtual void prune(
            const PruningContext& context,
            std::vector<PCSTFast::IndexType>& result_nodes, // Use consistent IndexType
            std::vector<PCSTFast::IndexType>& result_edges  // Edge indices are ints, matching phase1_result
        ) = 0;
    };


    // --- Base class for GW and Strong Pruning (contains common members/methods) ---
    // This helps reduce code duplication between GWPruner and StrongPruner.
    class AdvancedPrunerBase : public IPruner {
    protected:
        // Common state needed during the pruning phase by derived classes
        std::vector<PCSTFast::IndexType> phase2_result_local_; // Edges after basic connectivity check
        std::vector<PCSTFast::IndexType> phase3_result_local_; // Edges after GW/Strong specific logic
        std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>> phase3_neighbors_; // Adjacency list of phase 2 graph
        std::vector<uint8_t> node_deleted_;                    // Flags nodes deleted by the specific pruning logic

        // No queue members here - they are local to methods needing them

        // Common setup steps implemented once
        void build_phase2(const PruningContext& context); // Filters phase1_result based on node_good
        void build_phase3_adjacency(const PruningContext& context); // Builds phase3_neighbors_ from phase2_result_local_
        // Common helper to generate final node list from final edge list and node_good/node_deleted flags
        void build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set);

    public:
        virtual ~AdvancedPrunerBase() = default; // Virtual destructor

        // Common setup called by derived classes before their specific logic
        virtual void setup(const PruningContext& context);
    };

    // --- Concrete Pruner Class Declarations (interfaces defined here) ---

    // No pruning, just uses phase 1 result.
    class NoPruner : public IPruner {
    public:
        void prune(const PruningContext& context,
                   std::vector<PCSTFast::IndexType>& result_nodes,
                   std::vector<PCSTFast::IndexType>& result_edges) override;
    private:
        // Helper specific to NoPruner for building node set
        void build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set);
    };

    // Simple pruning based on connectivity to 'good' nodes (Phase 2).
    class SimplePruner : public IPruner {
    private:
        std::vector<PCSTFast::IndexType> phase2_result_local_; // Stores intermediate result
        void build_phase2(const PruningContext& context);
        // Helper specific to SimplePruner for building node set
        void build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set);
    public:
        void prune(const PruningContext& context,
                   std::vector<PCSTFast::IndexType>& result_nodes,
                   std::vector<PCSTFast::IndexType>& result_edges) override;
    };

    // Goemans-Williamson style pruning.
    class GWPruner : public AdvancedPrunerBase {
    private:
        std::vector<bool> cluster_necessary_local_; // Flags used during GW reverse pass

        void mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index);
        void mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index); // Uses local queue
        void run_gw_pruning(const PruningContext& context); // Main GW logic
    public:
        void prune(const PruningContext& context,
                   std::vector<PCSTFast::IndexType>& result_nodes,
                   std::vector<PCSTFast::IndexType>& result_edges) override;
    };

    // Strong pruning based on subtree payoffs.
    class StrongPruner : public AdvancedPrunerBase {
    private:
        // Strong pruning specific state (members needed across methods)
        std::vector<PCSTFast::IndexType> final_component_label_;
        std::vector<std::vector<PCSTFast::IndexType>> final_components_;
        PCSTFast::IndexType root_component_index_ = PCSTFast::kInvalidIndex;
        // State used by DFS passes
        std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>> strong_pruning_parent_;
        std::vector<PCSTFast::ValueType> strong_pruning_payoff_;

        // No stack members here - they are local to methods

        // Strong pruning specific methods
        void label_final_components(const PruningContext& context); // Uses local queue
        void label_component_recursive(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id); // Uses local queue
        PCSTFast::IndexType find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx); // Uses local stack
        void propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root, PCSTFast::ValueType& best_value); // Uses local stack
        void strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted); // Uses local stack
        void mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node); // Uses local queue
        void run_strong_pruning(const PruningContext& context); // Main strong pruning logic orchestration
    public:
        void prune(const PruningContext& context,
                   std::vector<PCSTFast::IndexType>& result_nodes,
                   std::vector<PCSTFast::IndexType>& result_edges) override;
        // Override setup to resize strong-specific members
        void setup(const PruningContext& context) override;
    };


    // Factory function declaration (creates the appropriate pruner instance)
    std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method);

} // namespace internal
} // namespace cluster_approx
