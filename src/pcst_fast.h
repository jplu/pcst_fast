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


class PCSTFast {
public:
    // --- Public Types ---
    using ValueType = double; // Consistent alias
    using IndexType = int;    // Consistent alias
    using PayloadType = int;  // Consistent alias

    // Logger function type definition
    using Logger = std::function<void(int level, const std::string& message)>;

    enum class PruningMethod {
        kNoPruning = 0,
        kSimplePruning,
        kGWPruning,
        kStrongPruning,
        kUnknownPruning,
    };

    // --- Public Algorithm Structures ---
    // Moved here so PruningContext can access them
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

    struct EdgeInfo {
        IndexType inactive_merge_event = kInvalidIndex; // Use IndexType
    };

    struct InactiveMergeEvent {
        IndexType active_cluster_index = kInvalidIndex;
        IndexType inactive_cluster_index = kInvalidIndex;
        IndexType active_cluster_node = kInvalidIndex;
        IndexType inactive_cluster_node = kInvalidIndex;
    };

    struct Cluster {
        PairingHeap<ValueType, PayloadType> edge_parts; // Use aliases
        bool active = false;
        ValueType active_start_time = 0.0;
        ValueType active_end_time = -1.0;
        IndexType merged_into = kInvalidIndex;
        ValueType prize_sum = 0.0;
        ValueType subcluster_moat_sum = 0.0;
        ValueType moat = 0.0;
        bool contains_root = false;
        IndexType skip_up = kInvalidIndex;
        ValueType skip_up_sum = 0.0;
        IndexType merged_along = kInvalidIndex;
        IndexType child_cluster_1 = kInvalidIndex;
        IndexType child_cluster_2 = kInvalidIndex;
        bool necessary = false; // Still used internally by GWPruner logic via context
        Cluster() = default;
    };

    // --- Constants ---
    static constexpr IndexType kNoRoot = -1;
    static constexpr IndexType kInvalidIndex = -1;


private:
    // --- Private Internal Types ---
    using PairingHeapType = PairingHeap<ValueType, PayloadType>;
    using PriorityQueueType = PriorityQueue<ValueType, IndexType>; // Min-heap

    // EdgePart can remain private as it's only manipulated internally by PCSTFast
    struct EdgePart {
        ValueType next_event_val = 0.0;
        bool deleted = false;
        PairingHeapType::ItemHandle heap_node = PairingHeapType::kInvalidHandle;
    };

    // --- Member Variables ---

    // Input Data (references)
    const std::vector<std::pair<IndexType, IndexType>>& edges_; // Use IndexType
    const std::vector<ValueType>& prizes_;
    const std::vector<ValueType>& costs_;
    const IndexType root_;
    const int target_num_active_clusters_; // Keep as int if it represents a count

    // Core Algorithm State (Order matches initializer list)
    Statistics stats_;
    std::vector<EdgePart> edge_parts_; // Size 2 * num_edges
    std::vector<EdgeInfo> edge_info_;  // Size num_edges
    std::vector<Cluster> clusters_;    // Stores all clusters
    std::vector<InactiveMergeEvent> inactive_merge_events_; // Logged for GW pruning
    PriorityQueueType clusters_deactivation_;
    PriorityQueueType clusters_next_edge_event_;
    ValueType current_time_ = 0.0;

    // Pruning Strategy (Polymorphic)
    std::unique_ptr<internal::IPruner> pruner_; // Holds the chosen pruning strategy

    // Logging (Order matches initializer list)
    Logger logger_; // Use std::function for logger
    const int verbosity_level_;

    // State passed to Pruners
    std::vector<uint8_t> node_good_;     // Node included before phase 3 pruning (size num_nodes)
    std::vector<int> phase1_result_;    // Edges selected by main GW loop (Keep as int if indices are int)

    // Configuration
    static constexpr ValueType kEpsilon = 1e-9;
    static constexpr int kOutputBufferSize = 1024; // For log_message formatting

    // --- Private Methods ---

    // Initialization Helpers
    void initialize_clusters();
    void initialize_edges();

    // Main Loop Event Handling
    struct NextEvent {
         ValueType time = std::numeric_limits<ValueType>::infinity();
         enum class Type { EDGE, CLUSTER_DEACTIVATION, NONE } type = Type::NONE;
         IndexType cluster_index = kInvalidIndex;
         PayloadType edge_part_index = kInvalidIndex; // Use PayloadType
    };
    [[nodiscard]] NextEvent get_next_event() const;
    void process_event(const NextEvent& event);
    void handle_edge_event(PayloadType edge_part_index); // Use PayloadType
    void handle_cluster_deactivation_event(IndexType cluster_index);

    // Edge Event Sub-logic Types/Helpers
    struct EdgeProcessingInfo {
        ValueType sum_current = 0.0;
        ValueType sum_other = 0.0;
        ValueType finished_moat_current = 0.0;
        ValueType finished_moat_other = 0.0;
        IndexType current_cluster_idx = kInvalidIndex;
        IndexType other_cluster_idx = kInvalidIndex;
        IndexType edge_idx = kInvalidIndex; // Index into original edges_ vector
        PayloadType current_part_idx = kInvalidIndex; // Index into edge_parts_ vector
        PayloadType other_part_idx = kInvalidIndex;   // Index into edge_parts_ vector
        ValueType cost = 0.0;
        ValueType remainder = 0.0;
    };
    [[nodiscard]] std::optional<EdgeProcessingInfo> get_edge_processing_info(PayloadType edge_part_index); // Use PayloadType
    void merge_clusters(const EdgeProcessingInfo& info);
    void handle_active_active_growth(const EdgeProcessingInfo& info);
    void handle_active_inactive_growth(const EdgeProcessingInfo& info);
    void update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index);
    void update_cluster_queues_post_growth(IndexType cluster_idx);
    void log_inactive_merge_event(const EdgeProcessingInfo& info);

    // Path Compression Helper
    void get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum, // Use PayloadType
                              ValueType* finished_moat_sum, IndexType* current_cluster_index);
    [[nodiscard]] IndexType find_representative_and_compress(IndexType start_node, ValueType& path_sum_out);

    // Helper needed before pruning starts
    void select_initial_active_clusters(std::vector<uint8_t>& node_good_flag);

    // Utility Helpers
    void format_and_log(int level, const char* format, ...) const; // Formats then calls logger_

    [[nodiscard]] static constexpr PayloadType get_other_edge_part_index(PayloadType edge_part_index) { // Use PayloadType
        return edge_part_index ^ 1;
    }

public:
    // --- Public Interface ---

    [[nodiscard]] static PruningMethod parse_pruning_method(std::string_view input);

    PCSTFast(const std::vector<std::pair<IndexType, IndexType>>& edges, // Use IndexType
             const std::vector<ValueType>& prizes,
             const std::vector<ValueType>& costs,
             IndexType root,
             int target_num_active_clusters,
             PruningMethod pruning,
             int verbosity_level,
             Logger logger = nullptr // Make logger optional, default to null/no-op
             );

    ~PCSTFast(); // Defined in .cc

    PCSTFast(const PCSTFast&) = delete;
    PCSTFast& operator=(const PCSTFast&) = delete;
    PCSTFast(PCSTFast&&) = delete;
    PCSTFast& operator=(PCSTFast&&) = delete;

    [[nodiscard]] bool run(std::vector<IndexType>* result_nodes, std::vector<IndexType>* result_edges); // Use IndexType

    void get_statistics(Statistics* s) const;

}; // class PCSTFast

} // namespace cluster_approx
