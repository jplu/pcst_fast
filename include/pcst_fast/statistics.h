#pragma once

#include <cstdint>

namespace cluster_approx {

struct Statistics {
    int64_t total_num_edge_events = 0;
    int64_t num_deleted_edge_events = 0;
    int64_t num_merged_edge_events = 0;
    int64_t total_num_merge_events = 0;
    int64_t num_active_active_merge_events = 0;
    int64_t num_active_inactive_merge_events = 0;
    int64_t total_num_edge_growth_events = 0;
    int64_t num_active_active_edge_growth_events = 0;
    int64_t num_active_inactive_edge_growth_events = 0;
    int64_t num_cluster_events = 0;

    Statistics();
};

}