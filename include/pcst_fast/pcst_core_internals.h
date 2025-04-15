#pragma once

#include "pcst_fast/pcst_types.h"
#include "pcst_fast/datastructures/pairing_heap.h"

#include <vector>
#include <limits>

namespace cluster_approx {

using PairingHeapType = PairingHeap<double, EdgePartId>;

/**
 * @brief Internal structure holding information about an edge relevant to the core algorithm.
 * Specifically used to track inactive merge events for GW pruning.
 */
struct EdgeInfo {
    EventId inactive_merge_event = kInvalidEventId;
};

/**
 * @brief Internal structure representing one "half" of an edge incident to a cluster.
 * Tracks the state of this edge part within its cluster's pairing heap.
 */
struct EdgePart {
    double next_event_val = std::numeric_limits<double>::infinity();
    bool deleted = false;
    PairingHeapType::ItemHandle heap_node = nullptr;
};

/**
 * @brief Internal structure representing a cluster during the Goemans-Williamson algorithm execution.
 */
struct Cluster {

    PairingHeapType edge_parts;
    bool active = false;
    double active_start_time = 0.0;

    double active_end_time = -1.0;
    ClusterId merged_into = kInvalidClusterId;
    double prize_sum = 0.0;
    double subcluster_moat_sum = 0.0;
    double moat = 0.0;
    bool contains_root = false;

    ClusterId skip_up = kInvalidClusterId;
    double skip_up_sum = 0.0;

    EdgeId merged_along = kInvalidEdgeId;
    ClusterId child_cluster_1 = kInvalidClusterId;
    ClusterId child_cluster_2 = kInvalidClusterId;

    bool necessary = false;

    /**
     * @brief Constructor requiring the shared buffer for the PairingHeap.
     */
    explicit Cluster(std::vector<PairingHeapType::ItemHandle>* heap_buffer)
        : edge_parts(heap_buffer) {}

    Cluster(Cluster&& other) noexcept = default;
    Cluster& operator=(Cluster&& other) noexcept = default;

    Cluster(const Cluster&) = delete;
    Cluster& operator=(const Cluster&) = delete;
};

}