#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/datastructures/priority_queue.h"

#include <vector>
#include <limits>
#include <optional>
#include <cstdint>

namespace cluster_approx {

struct InactiveMergeEvent;
struct SkipUpNode {
    ClusterId skip_up = kInvalidClusterId;
    double sum = 0.0;
};

class PCSTCoreAlgorithm {
  public:
    PCSTCoreAlgorithm(const GraphData& graph,
                      int target_num_active_clusters,
                      Logger* logger);

    ~PCSTCoreAlgorithm();

    PCSTCoreAlgorithm(const PCSTCoreAlgorithm&) = delete;
    PCSTCoreAlgorithm& operator=(const PCSTCoreAlgorithm&) = delete;

    PCSTCoreAlgorithm(PCSTCoreAlgorithm&&) noexcept = default;
    PCSTCoreAlgorithm& operator=(PCSTCoreAlgorithm&&) noexcept = default;

    [[nodiscard]] CoreAlgorithmResult run();

  private:
    using PriorityQueueType = PriorityQueue<double, ClusterId>;

    void initialize();
    void handle_edge_event(double edge_event_time, EdgePartId edge_part_index);
    void handle_cluster_event(double cluster_event_time, ClusterId cluster_index);
    ClusterId merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder);
    void remove_next_edge_event(ClusterId cluster_index);
    void mark_nodes_as_good(ClusterId start_cluster_index);
    [[nodiscard]] CoreAlgorithmResult build_core_result();

    FORCE_INLINE void get_sum_on_edge_part(EdgePartId edge_part_index,
      double* total_sum,
      double* finished_moat_sum,
      ClusterId* current_cluster_index);

    FORCE_INLINE int32_t ph_link(int32_t node1, int32_t node2);
    FORCE_INLINE int32_t ph_insert(int32_t root, int32_t node, double value);
    FORCE_INLINE void ph_add_to_heap(int32_t root, double value);
    FORCE_INLINE int32_t ph_decrease_key(int32_t root, int32_t node, double from_value, double to_value);
    FORCE_INLINE int32_t ph_delete_min(int32_t root, double* out_value, int32_t* out_node);
    FORCE_INLINE int32_t ph_meld(int32_t root1, int32_t root2);
    FORCE_INLINE bool ph_get_min(int32_t root, double* out_val, int32_t* out_node) const;
    FORCE_INLINE void update_cluster_edge_pq(ClusterId cluster_idx);

    [[nodiscard]] static FORCE_INLINE constexpr EdgePartId get_other_edge_part_index(EdgePartId edge_part_index) noexcept {
      return edge_part_index ^ 1;
    }

    const GraphData& graph_;
    int target_num_active_clusters_;
    Logger* logger_;

    double current_time_ = 0.0;
    double eps_ = 1e-6;
    int num_active_clusters_ = 0;
    int32_t num_created_clusters_ = 0;

    // Data-Oriented Design (SoA) for extreme cache throughput
    std::vector<int32_t> cluster_edge_parts_root_;
    std::vector<uint8_t> cluster_active_;
    std::vector<double> cluster_active_start_time_;
    std::vector<double> cluster_active_end_time_;
    std::vector<ClusterId> cluster_merged_into_;
    std::vector<double> cluster_prize_sum_;
    std::vector<double> cluster_subcluster_moat_sum_;
    std::vector<double> cluster_moat_;
    std::vector<uint8_t> cluster_contains_root_;
    std::vector<SkipUpNode> cluster_skip_up_;
    std::vector<EdgeId> cluster_merged_along_;
    std::vector<ClusterId> cluster_child_cluster_1_;
    std::vector<ClusterId> cluster_child_cluster_2_;

    std::vector<EdgePart> edge_parts_;
    std::vector<EventId> edge_inactive_merge_event_ids_;
    std::vector<InactiveMergeEvent> inactive_merge_events_;
    
    std::vector<PairingHeapNode> ph_nodes_;
    std::vector<int32_t> pairing_heap_buffer_;

    PriorityQueueType clusters_deactivation_;
    PriorityQueueType clusters_next_edge_event_;

    Statistics stats_;

    std::vector<uint8_t> node_good_;
    std::vector<EdgeId> phase1_result_edges_;

    std::vector<ClusterId> cluster_queue_;
};

}