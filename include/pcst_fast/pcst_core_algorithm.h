#pragma once

#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/datastructures/priority_queue.h"
#include "pcst_fast/datastructures/pairing_heap.h"

#include <vector>
#include <memory>
#include <limits>
#include <optional>

namespace cluster_approx {

struct Cluster;
struct EdgePart;
struct EdgeInfo;
struct InactiveMergeEvent;

/**
 * @brief Implements the core Goemans-Williamson based clustering algorithm for PCST.
 *
 * This class runs the main event loop, managing clusters, edge parts, and events
 * (merges, growths, deactivations) until the target number of active clusters is reached
 * or no more events can occur. It produces an intermediate result containing selected edges,
 * node filters, and event information required for subsequent pruning steps.
 */
class PCSTCoreAlgorithm {
  public:
    /**
     * @brief Constructs the core algorithm runner.
     * @param graph The input graph data (edges, prizes, costs, root). References held must remain valid for the lifetime of this object or until run() completes.
     * @param target_num_active_clusters The desired number of active clusters remaining at the end. Must be 0 for rooted problems.
     * @param logger A non-owning pointer to the logger instance. Must not be null and remain valid.
     * @throws std::invalid_argument if inputs are inconsistent (e.g., negative prizes/costs, invalid root, incorrect target_num_active_clusters for rooted).
     */
    PCSTCoreAlgorithm(const GraphData& graph,
                      int target_num_active_clusters,
                      Logger* logger);

    ~PCSTCoreAlgorithm();

    PCSTCoreAlgorithm(const PCSTCoreAlgorithm&) = delete;
    PCSTCoreAlgorithm& operator=(const PCSTCoreAlgorithm&) = delete;

    PCSTCoreAlgorithm(PCSTCoreAlgorithm&&) noexcept = default;
    PCSTCoreAlgorithm& operator=(PCSTCoreAlgorithm&&) noexcept = default;

    /**
     * @brief Executes the core Goemans-Williamson algorithm.
     * @return A CoreAlgorithmResult struct containing intermediate edges, node filters, event info, and statistics.
     * @throws std::runtime_error if an unexpected state occurs during execution.
     */
    [[nodiscard]] CoreAlgorithmResult run();

  private:
    using PairingHeapType = PairingHeap<double, EdgePartId>;
    using PriorityQueueType = PriorityQueue<double, ClusterId>;

    /**
     * @brief Initializes internal data structures (clusters, edge parts, heaps, queues).
     */
    void initialize();

    /**
     * @brief Handles the next edge event (growth or merge).
     * @param edge_event_time The time of the event.
     * @param edge_part_index The index of the edge part triggering the event.
     */
    void handle_edge_event(double edge_event_time, EdgePartId edge_part_index);

    /**
     * @brief Handles the next cluster deactivation event.
     * @param cluster_event_time The time of the event.
     * @param cluster_index The index of the cluster becoming inactive.
     */
    void handle_cluster_event(double cluster_event_time, ClusterId cluster_index);

    /**
     * @brief Merges two clusters into a new one.
     * @param cluster1_idx Index of the first cluster.
     * @param cluster2_idx Index of the second cluster.
     * @param merge_edge_idx Index of the edge causing the merge.
     * @param event_time Time of the merge event.
     * @param remainder Any remaining cost difference affecting the merge time/state.
     * @return Index of the newly created cluster.
     */
    ClusterId merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder);

    /**
     * @brief Retrieves the next edge event from the priority queue.
     * @return An optional pair containing {event_time, {cluster_index, edge_part_index}}, or nullopt if no events exist.
     */
    [[nodiscard]] std::optional<std::pair<double, std::pair<ClusterId, EdgePartId>>> get_next_edge_event();

    /**
     * @brief Removes the top edge event associated with a cluster from the queues/heaps.
     * Re-inserts the next event for that cluster if one exists.
     * @param cluster_index The cluster whose top event should be removed.
     */
    void remove_next_edge_event(ClusterId cluster_index);

    /**
     * @brief Retrieves the next cluster deactivation event from the priority queue.
     * @return An optional pair containing {event_time, cluster_index}, or nullopt if no events exist.
     */
    [[nodiscard]] std::optional<std::pair<double, ClusterId>> get_next_cluster_event();

    /**
     * @brief Removes the top cluster deactivation event from the priority queue.
     */
    void remove_next_cluster_event();

    /**
     * @brief Calculates the accumulated moat sum along the path from an edge part's endpoint up to its current cluster root.
     * Implements path compression optimization.
     * @param edge_part_index The index of the edge part.
     * @param[out] total_sum The total accumulated sum (moat sum + potentially active time).
     * @param[out] finished_moat_sum The sum including only completed moats (up to the root's start time if active).
     * @param[out] current_cluster_index The index of the representative cluster root found.
     */
    void get_sum_on_edge_part(EdgePartId edge_part_index,
                              double* total_sum,
                              double* finished_moat_sum,
                              ClusterId* current_cluster_index);

    /**
     * @brief Traverses the cluster merge tree starting from a final cluster to mark all underlying original nodes as 'good'.
     * Uses BFS/DFS internally.
     * @param start_cluster_index The index of the final (root or active) cluster to start from.
     */
    void mark_nodes_as_good(ClusterId start_cluster_index);

    /**
     * @brief Constructs the final CoreAlgorithmResult object after the main loop finishes.
     * Filters phase1 edges based on `node_good_` and populates the result struct.
     * @return The completed CoreAlgorithmResult.
     */
    [[nodiscard]] CoreAlgorithmResult build_core_result();

    /**
     * @brief Helper to get the other edge part index for a given edge part index.
     * @param edge_part_index The input edge part index (0 to 2*num_edges - 1).
     * @return The index of the corresponding edge part for the other endpoint.
     */
    [[nodiscard]] static constexpr EdgePartId get_other_edge_part_index(EdgePartId edge_part_index) noexcept {

        return (edge_part_index % 2 == 0) ? (edge_part_index + 1) : (edge_part_index - 1);
    }

    const GraphData& graph_;
    int target_num_active_clusters_;
    Logger* logger_;

    double current_time_ = 0.0;
    double eps_ = 1e-6;
    int num_active_clusters_ = 0;

    std::vector<Cluster> clusters_;
    std::vector<EdgePart> edge_parts_;
    std::vector<EdgeInfo> edge_info_;
    std::vector<InactiveMergeEvent> inactive_merge_events_;

    std::vector<PairingHeapType::ItemHandle> pairing_heap_buffer_;

    PriorityQueueType clusters_deactivation_;
    PriorityQueueType clusters_next_edge_event_;

    Statistics stats_;

    std::vector<bool> node_good_;
    std::vector<EdgeId> phase1_result_edges_;

    std::vector<std::pair<ClusterId, double>> path_compression_visited_;
    std::vector<ClusterId> cluster_queue_;
};

}