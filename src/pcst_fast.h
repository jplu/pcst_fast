#pragma once

#include <string_view>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>   // For std::string in parse_pruning_method
#include <limits>  // For numeric_limits
#include <optional> // Potentially useful alternative return types

#include "pairing_heap.h"
#include "priority_queue.h"

namespace cluster_approx {

    /**
     * @brief Implements the "Fast" Prize-Collecting Steiner Tree (PCST) algorithm.
     *
     * Uses a growth-based strategy (like GW) with pairing heaps and priority queues
     * to efficiently find approximate solutions to the PCST problem. Supports
     * various pruning methods.
     *
     * @note This class manages significant internal state and its methods (especially `run`)
     *       are NOT thread-safe for concurrent calls on the *same* instance.
     *       However, multiple independent `PCSTFast` instances can be created and run
     *       concurrently in different threads on different problems.
     */
    class PCSTFast {
    public:
        // --- Public Types ---

        enum class PruningMethod {
            kNoPruning = 0,
            kSimplePruning, // Phase 1 + connectivity check
            kGWPruning,     // Goemans-Williamson style pruning
            kStrongPruning, // Stronger pruning based on subtree payoffs
            kUnknownPruning,
        };

        struct Statistics {
            long long total_num_edge_events = 0;
            long long num_deleted_edge_events = 0;
            long long num_merged_edge_events = 0; // Edges connecting already merged clusters
            long long total_num_merge_events = 0; // Events that merge two clusters
            long long num_active_active_merge_events = 0;
            long long num_active_inactive_merge_events = 0;
            long long total_num_edge_growth_events = 0; // Events that advance edge processing
            long long num_active_active_edge_growth_events = 0;
            long long num_active_inactive_edge_growth_events = 0;
            long long num_cluster_events = 0; // Deactivation events

            // Reset stats to zero
            void reset() {
                total_num_edge_events = 0;
                num_deleted_edge_events = 0;
                num_merged_edge_events = 0;
                total_num_merge_events = 0;
                num_active_active_merge_events = 0;
                num_active_inactive_merge_events = 0;
                total_num_edge_growth_events = 0;
                num_active_active_edge_growth_events = 0;
                num_active_inactive_edge_growth_events = 0;
                num_cluster_events = 0;
            }
        };

        // --- Constants ---
        static constexpr int kNoRoot = -1;
        static constexpr int kInvalidIndex = -1; // General invalid index/pointer

    private:
        // --- Internal Types ---
        // Note: Using double for values/costs. Consider a custom fixed-point type
        //       or careful comparison with epsilon if precision issues arise.
        using ValueType = double;
        using PayloadType = int; // Often edge part index or cluster index
        using IndexType = int;   // Node indices, cluster indices
        using PairingHeapType = PairingHeap<ValueType, PayloadType>;
        using PriorityQueueType = PriorityQueue<ValueType, IndexType>; // Using default std::less (min-heap)

        struct EdgeInfo {
            // Index into inactive_merge_events_ if this edge caused such an event during GW phase.
            int inactive_merge_event = kInvalidIndex;
        };

        // Represents one 'half' of an edge incident to a cluster.
        struct EdgePart {
            ValueType next_event_val = 0.0; // Time (or accumulated cost) for next event
            bool deleted = false;           // If this edge part is logically deleted
            PairingHeapType::ItemHandle heap_node = PairingHeapType::kInvalidHandle; // Handle in cluster's heap
        };

        // Stores information about a merge event involving an inactive cluster (for GW pruning).
        struct InactiveMergeEvent {
            IndexType active_cluster_index = kInvalidIndex;   // Representative cluster index (at time of merge)
            IndexType inactive_cluster_index = kInvalidIndex; // Representative cluster index (at time of merge)
            IndexType active_cluster_node = kInvalidIndex;    // Original node on the active side
            IndexType inactive_cluster_node = kInvalidIndex;  // Original node on the inactive side
        };

        // Represents a cluster during the algorithm's execution.
        struct Cluster {
            PairingHeapType edge_parts; // Min-heap of edge parts incident to this cluster (by event time)
            bool active = false;
            ValueType active_start_time = 0.0; // Time when cluster became active
            ValueType active_end_time = -1.0; // Time when cluster became inactive (-1 if still active)
            IndexType merged_into = kInvalidIndex; // Index of the cluster this one merged into
            ValueType prize_sum = 0.0;         // Sum of prizes of original nodes in this cluster
            ValueType subcluster_moat_sum = 0.0; // Sum of moats of subclusters merged into this one
            ValueType moat = 0.0;             // Cost accumulated while active (active_end_time - active_start_time)
            bool contains_root = false;       // True if the designated root node is in this cluster
            // Path compression for get_sum_on_edge_part
            IndexType skip_up = kInvalidIndex;
            ValueType skip_up_sum = 0.0;
            // Information for reconstructing the tree / GW pruning
            IndexType merged_along = kInvalidIndex;  // Edge index used to merge children into this cluster
            IndexType child_cluster_1 = kInvalidIndex; // Child cluster index (representative at time of merge)
            IndexType child_cluster_2 = kInvalidIndex; // Child cluster index (representative at time of merge)
            bool necessary = false;            // Used during GW/Strong pruning phases

            // Cluster doesn't need external buffer passed in anymore
            Cluster() = default; // Rely on PairingHeap default constructor
        };

        // --- Member Variables ---

        // Input Data (references to external data)
        const std::vector<std::pair<int, int> >& edges_;
        const std::vector<ValueType>& prizes_;
        const std::vector<ValueType>& costs_;
        const IndexType root_; // kNoRoot if unrooted
        const int target_num_active_clusters_; // Target for unrooted variant termination
        const PruningMethod pruning_;
        const int verbosity_level_;
        void (*output_function_)(const char*); // Optional logging callback

        // Algorithm State
        Statistics stats_;
        std::vector<EdgePart> edge_parts_; // Size 2 * num_edges
        std::vector<EdgeInfo> edge_info_;  // Size num_edges
        std::vector<Cluster> clusters_;    // Stores all clusters (initial nodes + merged clusters)
        std::vector<InactiveMergeEvent> inactive_merge_events_; // Logged for GW pruning
        PriorityQueueType clusters_deactivation_; // Min-heap: {deactivation_time, cluster_index}
        PriorityQueueType clusters_next_edge_event_; // Min-heap: {min_edge_event_time, cluster_index}
        ValueType current_time_ = 0.0;

        // Pruning/Output Generation State
        std::vector<uint8_t> node_good_;     // Node included in the result before phase 3 pruning? (size num_nodes)
        std::vector<uint8_t> node_deleted_;  // Node pruned by GW/Strong pruning? (size num_nodes)
        std::vector<int> phase1_result_;    // Edges selected by main GW loop
        std::vector<int> phase2_result_;    // Edges after simple pruning (connectivity to 'good' nodes)
        std::vector<int> phase3_result_;    // Edges after GW/Strong pruning

        // Temporary Buffers (used by specific methods to avoid repeated allocations)
        std::vector<std::pair<IndexType, ValueType>> path_compression_visited_; // For get_sum_on_edge_part
        std::vector<IndexType> cluster_queue_; // General purpose queue for BFS/DFS traversals
        std::vector<std::vector<std::pair<IndexType, ValueType>>> phase3_neighbors_; // Adjacency list for pruning phases
        std::vector<uint8_t> build_phase1_included_nodes_; // Temp for build_phase1_node_set
        // Strong Pruning Specific State
        std::vector<IndexType> final_component_label_; // Maps node index -> component index
        std::vector<std::vector<IndexType>> final_components_; // Stores nodes for each component
        IndexType root_component_index_ = kInvalidIndex;
        std::vector<std::pair<IndexType, ValueType>> strong_pruning_parent_; // {parent_idx, edge_cost} in DFS tree
        std::vector<ValueType> strong_pruning_payoff_; // Accumulated payoff in DFS tree
        std::vector<std::pair<bool, IndexType>> strong_pruning_stack_; // DFS stack: {is_pre_order, node_idx}
        std::vector<IndexType> strong_pruning_stack2_; // Stack for find_best_component_root

        // Configuration
        static constexpr ValueType kEpsilon = 1e-9; // Tolerance for float comparisons
        static constexpr int kOutputBufferSize = 1024; // Size for local logging buffer

        // --- Private Methods ---

        // Initialization Helpers
        void initialize_clusters();
        void initialize_edges();

        // Main Loop Event Handling
        struct NextEvent {
             ValueType time = std::numeric_limits<ValueType>::infinity();
             enum class Type { EDGE, CLUSTER_DEACTIVATION, NONE } type = Type::NONE;
             IndexType cluster_index = kInvalidIndex;
             PayloadType edge_part_index = kInvalidIndex; // Only for edge events
        };
        [[nodiscard]] NextEvent get_next_event() const;
        void process_event(const NextEvent& event);
        void handle_edge_event(PayloadType edge_part_index);
        void handle_cluster_deactivation_event(IndexType cluster_index);

        // Edge Event Sub-logic
        struct EdgeProcessingInfo {
            ValueType sum_current = 0.0;
            ValueType sum_other = 0.0;
            ValueType finished_moat_current = 0.0;
            ValueType finished_moat_other = 0.0;
            IndexType current_cluster_idx = kInvalidIndex;
            IndexType other_cluster_idx = kInvalidIndex;
            IndexType edge_idx = kInvalidIndex;
            PayloadType current_part_idx = kInvalidIndex;
            PayloadType other_part_idx = kInvalidIndex;
            ValueType cost = 0.0;
            ValueType remainder = 0.0;
        };
        [[nodiscard]] std::optional<EdgeProcessingInfo> get_edge_processing_info(PayloadType edge_part_index);
        void merge_clusters(const EdgeProcessingInfo& info);
        void handle_active_active_growth(const EdgeProcessingInfo& info);
        void handle_active_inactive_growth(const EdgeProcessingInfo& info);
        void update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index);
        void update_cluster_queues_post_growth(IndexType cluster_idx);
        void log_inactive_merge_event(const EdgeProcessingInfo& info);


        // Path Compression Helper
        void get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                                  ValueType* finished_moat_sum, IndexType* current_cluster_index);
        [[nodiscard]] IndexType find_representative_and_compress(IndexType start_node, ValueType& path_sum_out);


        // Pruning Phase Helpers
        void select_initial_active_clusters(std::vector<uint8_t>& node_good_flag); // Renamed from mark_nodes_as_good
        void build_phase2_result();
        void build_phase3_adjacency();

        // GW Pruning
        void prune_gw();
        void mark_clusters_as_necessary_gw(IndexType start_cluster_index); // Traverses merge tree
        void mark_nodes_as_deleted_gw(IndexType start_node_index, IndexType parent_node_index); // Traverses phase3 graph


        // Strong Pruning
        void prune_strong();
        void label_final_components(); // Builds final_components_ and final_component_label_
        void label_component_recursive(IndexType start_node_index, IndexType component_id); // DFS/BFS for labeling
        [[nodiscard]] IndexType find_best_component_root(IndexType component_index); // Finds best root for payoff calc
        ValueType calculate_initial_payoffs(IndexType start_node_index); // First pass DFS for payoff
        void propagate_payoffs_and_find_best(IndexType initial_root, IndexType& best_root_out, ValueType& best_value_out); // Second pass for rerooting
        void strong_pruning_dfs(IndexType start_node_index, bool mark_deleted); // Final pruning DFS
        void mark_nodes_as_deleted_strong(IndexType start_node_index, IndexType parent_node_index); // Subtree deletion


        // Output Generation Helpers
        void build_node_set_from_edges(const std::vector<int>& edge_set, std::vector<int>* node_set); // Common helper
        void build_final_node_set(std::vector<int>* node_set) const; // Uses node_good_ and node_deleted_

        // Utility Helpers
        void log_message(int level, const char* format, ...) const;
        [[nodiscard]] static constexpr int get_other_edge_part_index(int edge_part_index) {
            return edge_part_index ^ 1;
        }

    public:
        // --- Public Interface ---

        [[nodiscard]] static PruningMethod parse_pruning_method(std::string_view input);

        PCSTFast(const std::vector<std::pair<int, int> >& edges,
                 const std::vector<ValueType>& prizes,
                 const std::vector<ValueType>& costs,
                 IndexType root, // Use IndexType consistently
                 int target_num_active_clusters,
                 PruningMethod pruning,
                 int verbosity_level,
                 void (*output_function)(const char*) = nullptr); // Default logger to null

        ~PCSTFast() = default; // Rely on member destructors

        // Delete copy/move operations - instances manage significant state and aren't meant to be copied/moved easily.
        PCSTFast(const PCSTFast&) = delete;
        PCSTFast& operator=(const PCSTFast&) = delete;
        PCSTFast(PCSTFast&&) = delete;
        PCSTFast& operator=(PCSTFast&&) = delete;


        /**
         * @brief Runs the PCST algorithm.
         * @param result_nodes Pointer to a vector to store the indices of nodes in the result Steiner tree.
         * @param result_edges Pointer to a vector to store the indices of edges in the result Steiner tree.
         * @return True if the algorithm runs successfully, false on configuration errors.
         * @throws std::invalid_argument if result pointers are null.
         */
        [[nodiscard]] bool run(std::vector<int>* result_nodes, std::vector<int>* result_edges);

        /**
         * @brief Gets the statistics collected during the last run.
         * @param s Pointer to a Statistics struct to fill. Must not be null.
         */
        void get_statistics(Statistics* s) const;

    }; // class PCSTFast

} // namespace cluster_approx
