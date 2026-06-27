#pragma once

#include "pcst_fast/pcst_types.h"
#include <vector>
#include <limits>

namespace cluster_approx {

struct alignas(32) PairingHeapNode {
    double value;
    double child_offset;
    int32_t sibling;
    int32_t child;
    int32_t left_up;
};

struct EdgePart {
    double next_event_val = std::numeric_limits<double>::infinity();
    NodeId endpoint_node = kInvalidNodeId;
    bool deleted = false;
    bool in_heap = false;
};

struct Cluster {
    int32_t edge_parts_root = -1;
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
};

}