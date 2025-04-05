#pragma once

#include <string_view>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <limits>
#include <optional>
#include <functional> // For std::function (Logger)
#include <memory>     // For std::unique_ptr (IPruner)
#include <cstdarg>    // For va_list etc. in log_message
#include <cstdio>     // For vsnprintf in log_message
#include <stdexcept>  // For exceptions

#include "pairing_heap.h"
#include "priority_queue.h"

namespace cluster_approx {

// Forward declaration for internal types used by PCSTFast public interface
namespace internal {
    class IPruner; // Forward declare the interface
}


/**
 * @brief Implements the "Fast" Prize-Collecting Steiner Tree (PCST) algorithm.
 *
 * Uses a growth-based strategy (like GW) with pairing heaps and priority queues
 * to efficiently find approximate solutions to the PCST problem. Supports
 * various pruning methods delegated to a Pruner strategy object.
 *
 * @note This class manages significant internal state and its methods (especially `run`)
 *       are NOT thread-safe for concurrent calls on the *same* instance.
 *       However, multiple independent `PCSTFast` instances can be created and run
 *       concurrently in different threads on different problems. Input data referenced
 *       via const& must not be modified concurrently by other threads. Output vectors
 *       passed by pointer must be distinct per thread or protected by locks externally.
 *       The provided Logger must be thread-safe if used concurrently.
 */
class PCSTFast {
public:
    // --- Public Types ---
    using ValueType = double; // Consistent alias
    using IndexType = int;    // Consistent alias
    using PayloadType = int;  // Consistent alias (e.g., edge part index, cluster index)

    // Logger function type definition
    // Parameters: LogLevel, Message string
    using Logger = std::function<void(int level, const std::string& message)>;

    enum class PruningMethod {
        kNoPruning = 0,
        kSimplePruning, // Phase 1 + connectivity check
        kGWPruning,     // Goemans-Williamson style pruning
        kStrongPruning, // Stronger pruning based on subtree payoffs
        kUnknownPruning,
    };

    // --- Public Algorithm Structures ---
    // (These are defined here so the external PruningContext can use them)
    struct Statistics {
        long long total_num_edge_events = 0;
        long long num_deleted_edge_events = 0;
        long long num_merged_edge_events = 0;
        long long total_num_merge_events = 0;
        long long num_active_active_merge_events = 0;
        long long num_active_inactive_merge_events = 0;
        long long total_num_edge_growth_events = 0;
        long long num_active_active_edge_growth_events = 0;
        long long num_active_inactive_edge_growth_events = 0;
        long long num_cluster_events = 0;

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

    // Information attached to original edges, primarily for pruning
    struct EdgeInfo {
        IndexType inactive_merge_event = kInvalidIndex; // Index into inactive_merge_events_ if applicable
    };

    // Information logged during an Active-Inactive merge event for GW pruning
    struct InactiveMergeEvent {
        IndexType active_cluster_index = kInvalidIndex;   // Representative active cluster index (at time of merge)
        IndexType inactive_cluster_index = kInvalidIndex; // Representative inactive cluster index (at time of merge)
        IndexType active_cluster_node = kInvalidIndex;    // Original node on the active side of merging edge
        IndexType inactive_cluster_node = kInvalidIndex;  // Original node on the inactive side of merging edge
    };

    // Represents a cluster (initially single nodes, later merged entities)
    struct Cluster {
        PairingHeap<ValueType, PayloadType> edge_parts; // Min-heap of incident edge parts
        bool active = false;                            // Is the cluster currently growing?
        ValueType active_start_time = 0.0;              // Time it became active
        ValueType active_end_time = -1.0;               // Time it became inactive (-1 if still active/never active)
        IndexType merged_into = kInvalidIndex;          // Index of the parent cluster in the merge tree
        ValueType prize_sum = 0.0;                      // Sum of original node prizes in this cluster
        ValueType subcluster_moat_sum = 0.0;            // Sum of final moats of children merged into this
        ValueType moat = 0.0;                           // Cost accumulated while active (active_end_time - active_start_time)
        bool contains_root = false;                     // Does this cluster contain the designated root node?
        // Path compression fields
        IndexType skip_up = kInvalidIndex;              // Shortcut pointer for find_representative
        ValueType skip_up_sum = 0.0;                    // Accumulated moat sum along skip path
        // Merge tree structure fields
        IndexType merged_along = kInvalidIndex;         // Edge index used to form this cluster by merging children
        IndexType child_cluster_1 = kInvalidIndex;      // First child cluster index (representative at merge time)
        IndexType child_cluster_2 = kInvalidIndex;      // Second child cluster index (representative at merge time)
        // Pruning related flags
        bool necessary = false;                         // Used by GW pruner logic via PruningContext

        Cluster() = default; // Default constructor is sufficient
    };

    // --- Constants ---
    static constexpr IndexType kNoRoot = -1;     // Value indicating an unrooted problem
    static constexpr IndexType kInvalidIndex = -1; // General invalid index marker

private:
    // --- Private Internal Types ---
    using PairingHeapType = PairingHeap<ValueType, PayloadType>;
    using PriorityQueueType = PriorityQueue<ValueType, IndexType>; // Min-heap by default

    // Represents one 'half' of an edge incident to a cluster
    // Kept private as its manipulation is internal to PCSTFast logic
    struct EdgePart {
        ValueType next_event_val = 0.0; // Time/cost for the next event associated with this part
        bool deleted = false;           // Flag indicating if this edge part is logically removed
        PairingHeapType::ItemHandle heap_node = PairingHeapType::kInvalidHandle; // Handle in the cluster's heap
    };

    // --- Member Variables ---

    // Input Data (references passed in constructor)
    const std::vector<std::pair<IndexType, IndexType>>& edges_; // Graph edges {u, v}
    const std::vector<ValueType>& prizes_;                      // Node prizes
    const std::vector<ValueType>& costs_;                       // Edge costs
    const IndexType root_;                                      // Designated root node (kNoRoot if unrooted)
    const int target_num_active_clusters_;                      // Termination condition for unrooted variant

    // Core Algorithm State (managed during run)
    Statistics stats_;                                         // Performance counters
    std::vector<EdgePart> edge_parts_;                         // Stores state for each half-edge (size 2 * num_edges)
    std::vector<EdgeInfo> edge_info_;                          // Extra info per original edge (size num_edges)
    std::vector<Cluster> clusters_;                            // Stores all clusters (nodes + merged clusters)
    std::vector<InactiveMergeEvent> inactive_merge_events_;    // Log for GW pruning
    PriorityQueueType clusters_deactivation_;                  // PQ for cluster deactivation events {time, cluster_idx}
    PriorityQueueType clusters_next_edge_event_;               // PQ for next edge events {min_edge_time, cluster_idx}
    ValueType current_time_ = 0.0;                             // Current time in the simulation

    // Pruning Strategy (chosen at construction)
    std::unique_ptr<internal::IPruner> pruner_;                // Polymorphic pruner object

    // Logging Configuration
    Logger logger_;                                            // Logger function object
    const int verbosity_level_;                                // Controls level of detail in logging

    // State shared with Pruning Strategy via PruningContext
    std::vector<uint8_t> node_good_;                           // Flags indicating nodes connected to the final active/root set before pruning
    std::vector<int> phase1_result_;                           // Indices of edges selected during the main growth phase

    // Algorithm Configuration Constants
    static constexpr ValueType kEpsilon = 1e-9;                // Tolerance for floating-point comparisons
    static constexpr int kOutputBufferSize = 1024;             // Buffer size for log formatting

    // --- Private Methods ---

    // Initialization Helpers
    void initialize_clusters();
    void initialize_edges();

    // Main Loop Event Handling
    struct NextEvent {
         ValueType time = std::numeric_limits<ValueType>::infinity();
         enum class Type { EDGE, CLUSTER_DEACTIVATION, NONE } type = Type::NONE;
         IndexType cluster_index = kInvalidIndex;
         PayloadType edge_part_index = kInvalidIndex; // Index into edge_parts_ vector
    };
    [[nodiscard]] NextEvent get_next_event() const;
    void process_event(const NextEvent& event);
    void handle_edge_event(PayloadType edge_part_index); // Handles logic when an edge event occurs
    void handle_cluster_deactivation_event(IndexType cluster_index); // Handles logic when a cluster deactivates

    // Edge Event Sub-logic Types/Helpers
    // Structure to hold information needed during edge event processing
    struct EdgeProcessingInfo {
        ValueType sum_current = 0.0;              // Accumulated cost on the current edge part side
        ValueType sum_other = 0.0;                // Accumulated cost on the other edge part side
        ValueType finished_moat_current = 0.0;    // Sum of completed moats up to current representative
        ValueType finished_moat_other = 0.0;      // Sum of completed moats up to other representative
        IndexType current_cluster_idx = kInvalidIndex; // Representative cluster for current side
        IndexType other_cluster_idx = kInvalidIndex;   // Representative cluster for other side
        IndexType edge_idx = kInvalidIndex;        // Index into original edges_ vector
        PayloadType current_part_idx = kInvalidIndex; // Index into edge_parts_ vector for current side
        PayloadType other_part_idx = kInvalidIndex;   // Index into edge_parts_ vector for other side
        ValueType cost = 0.0;                     // Original cost of the edge
        ValueType remainder = 0.0;                // Remaining cost to be covered (cost - sum_current - sum_other)
    };
    [[nodiscard]] std::optional<EdgeProcessingInfo> get_edge_processing_info(PayloadType edge_part_index); // Gathers info for edge event
    void merge_clusters(const EdgeProcessingInfo& info);                     // Logic for merging two clusters
    void handle_active_active_growth(const EdgeProcessingInfo& info);      // Logic for edge growth between two active clusters
    void handle_active_inactive_growth(const EdgeProcessingInfo& info);    // Logic for edge growth between active and inactive clusters
    void update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index); // Updates PQs after merge
    void update_cluster_queues_post_growth(IndexType cluster_idx);        // Updates PQs after edge growth event for a cluster
    void log_inactive_merge_event(const EdgeProcessingInfo& info);         // Records event for GW pruning

    // Path Compression Helper
    // Calculates accumulated cost/moats and finds representative cluster for an edge part
    void get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                              ValueType* finished_moat_sum, IndexType* current_cluster_index);
    // Core path compression logic: finds representative and updates skip pointers
    [[nodiscard]] IndexType find_representative_and_compress(IndexType start_node, ValueType& path_sum_out);

    // Helper needed before pruning starts
    // Traverses merge tree to identify nodes connected to final active/root clusters
    void select_initial_active_clusters(std::vector<uint8_t>& node_good_flag);

    // Utility Helpers
    // Formats a log message and passes it to the configured logger_
    void format_and_log(int level, const char* format, ...) const;

    // Helper to get the index of the other part of an edge
    [[nodiscard]] static constexpr PayloadType get_other_edge_part_index(PayloadType edge_part_index) {
        return edge_part_index ^ 1;
    }

public:
    // --- Public Interface ---

    /**
     * @brief Parses a string representation of a pruning method name.
     * @param input The string name (e.g., "none", "gw", "strong"). Case-insensitive.
     * @return The corresponding PruningMethod enum value, or kUnknownPruning if not recognized.
     */
    [[nodiscard]] static PruningMethod parse_pruning_method(std::string_view input);

    /**
     * @brief Constructs a PCSTFast instance.
     * @param edges Vector of edges, each represented as a pair of node indices.
     * @param prizes Vector of node prizes. Must be non-negative. Size defines number of nodes.
     * @param costs Vector of edge costs corresponding to the edges vector. Must be non-negative.
     * @param root The index of the root node, or kNoRoot for the unrooted variant.
     * @param target_num_active_clusters Target number of clusters for termination in unrooted variant (ignored if root >= 0).
     * @param pruning The pruning method to use after the main algorithm phase.
     * @param verbosity_level Controls the level of detail for logging (0=silent).
     * @param logger Optional logger function (defaults to no-op). If provided, must be thread-safe if used concurrently.
     * @throws std::invalid_argument If input parameters are invalid (e.g., negative costs/prizes, bad root index, unknown pruning).
     */
    PCSTFast(const std::vector<std::pair<IndexType, IndexType>>& edges,
             const std::vector<ValueType>& prizes,
             const std::vector<ValueType>& costs,
             IndexType root,
             int target_num_active_clusters,
             PruningMethod pruning,
             int verbosity_level,
             Logger logger = nullptr
             );

    /**
     * @brief Destructor.
     */
    ~PCSTFast(); // Defined in .cc needed because of unique_ptr to forward-declared type

    // Delete copy/move operations to prevent accidental slicing or state mismanagement.
    PCSTFast(const PCSTFast&) = delete;
    PCSTFast& operator=(const PCSTFast&) = delete;
    PCSTFast(PCSTFast&&) = delete;
    PCSTFast& operator=(PCSTFast&&) = delete;

    /**
     * @brief Runs the PCST algorithm and pruning.
     * @param result_nodes Pointer to a vector where the indices of nodes in the final solution will be stored. Cleared before use.
     * @param result_edges Pointer to a vector where the indices of edges in the final solution will be stored. Cleared before use.
     * @return True if the algorithm runs successfully. False might indicate internal errors (though exceptions are preferred).
     * @throws std::invalid_argument if result pointers are null.
     */
    [[nodiscard]] bool run(std::vector<IndexType>* result_nodes, std::vector<IndexType>* result_edges);

    /**
     * @brief Gets the statistics collected during the last run.
     * @param s Pointer to a Statistics struct to fill. Must not be null.
     */
    void get_statistics(Statistics* s) const;

}; // class PCSTFast

} // namespace cluster_approx
