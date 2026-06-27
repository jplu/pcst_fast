#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/datastructures/priority_queue.h"
#include "pcst_fast/datastructures/pairing_heap.h"

#include <vector>
#include <memory>
#include <limits>
#include <optional>
#include <cstdint>

namespace cluster_approx {

struct Cluster;
struct EdgePart;
struct EdgeInfo;
struct InactiveMergeEvent;

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
    using PairingHeapType = PairingHeap<double, EdgePartId>;
    using PriorityQueueType = PriorityQueue<double, ClusterId>;

    void initialize();
    void handle_edge_event(double edge_event_time, EdgePartId edge_part_index);
    void handle_cluster_event(double cluster_event_time, ClusterId cluster_index);
    ClusterId merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder);
    [[nodiscard]] std::optional<std::pair<double, std::pair<ClusterId, EdgePartId>>> get_next_edge_event();
    void remove_next_edge_event(ClusterId cluster_index);
    [[nodiscard]] std::optional<std::pair<double, ClusterId>> get_next_cluster_event();
    void remove_next_cluster_event();
    void get_sum_on_edge_part(EdgePartId edge_part_index,
                              double* total_sum,
                              double* finished_moat_sum,
                              ClusterId* current_cluster_index);
    void mark_nodes_as_good(ClusterId start_cluster_index);
    [[nodiscard]] CoreAlgorithmResult build_core_result();

    [[nodiscard]] static constexpr EdgePartId get_other_edge_part_index(EdgePartId edge_part_index) noexcept {
        return (edge_part_index % 2 == 0) ? (edge_part_index + 1) : (edge_part_index - 1);
    }

    const GraphData& graph_;
    int target_num_active_clusters_;
    Logger* logger_;

    double current_time_ = 0.0;
    double eps_ = 1e-6;
    int num_active_clusters_ = 0;

    // Structural cluster elements used for tracking hierarchy and local pairing heaps
    std::vector<Cluster> clusters_;
    
    // Flat parallel vectors for hot path loop variables
    std::vector<ClusterId> cluster_merged_into_;
    std::vector<ClusterId> cluster_skip_up_;
    std::vector<double> cluster_skip_up_sum_;
    std::vector<double> cluster_moat_;
    std::vector<uint8_t> cluster_active_;
    std::vector<double> cluster_active_start_time_;
    std::vector<double> cluster_active_end_time_;

    std::vector<EdgePart> edge_parts_;
    std::vector<EdgeInfo> edge_info_;
    std::vector<InactiveMergeEvent> inactive_merge_events_;
    
    std::unique_ptr<PairingHeapType::AllocatorType> heap_node_allocator_;
    std::vector<PairingHeapType::ItemHandle> pairing_heap_buffer_;

    PriorityQueueType clusters_deactivation_;
    PriorityQueueType clusters_next_edge_event_;

    Statistics stats_;

    std::vector<bool> node_good_;
    std::vector<EdgeId> phase1_result_edges_;

    std::vector<ClusterId> cluster_queue_;
};

}