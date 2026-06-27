#pragma once

#include "pcst_fast/pcst_types.h"
#include "pcst_fast/datastructures/pairing_heap.h"
#include <vector>
#include <limits>

namespace cluster_approx {

using PairingHeapType = PairingHeap<double, EdgePartId>;

struct EdgeInfo {
    EventId inactive_merge_event = kInvalidEventId;
};

struct EdgePart {
    double next_event_val = std::numeric_limits<double>::infinity();
    bool deleted = false;
    // Replacing pointer-based handles with 32-bit index handles.
    // -1 represents a null/invalid handle.
    PairingHeapType::ItemHandle heap_node = -1;
};

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

    // Constructor updated to take the index-based allocator and workspace buffer
    Cluster(PairingHeapType::AllocatorType* allocator, std::vector<PairingHeapType::ItemHandle>* heap_buffer)
        : edge_parts(allocator, heap_buffer) {}

    Cluster(Cluster&& other) noexcept = default;
    Cluster& operator=(Cluster&& other) noexcept = default;

    Cluster(const Cluster&) = delete;
    Cluster& operator=(const Cluster&) = delete;
};

}