#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pcst_core_internals.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>
#include <format>

namespace cluster_approx {

PCSTCoreAlgorithm::PCSTCoreAlgorithm(const GraphData& graph,
                                     int target_num_active_clusters,
                                     Logger* logger)
    : graph_(graph),
      target_num_active_clusters_(target_num_active_clusters),
      logger_(logger) {
    assert(logger_ != nullptr && "Logger cannot be null.");
    const size_t num_nodes = graph_.prizes.size();
    const size_t num_edges = graph_.edges.size();

    if (graph_.root != kInvalidNodeId && target_num_active_clusters != 0) {
        logger_->log(LogLevel::ERROR,
                     "Target number of active clusters ({}) must be 0 for rooted problems (root = {}).",
                     target_num_active_clusters, graph_.root);
        throw std::invalid_argument(std::format(
                                        "Target number of active clusters ({}) must be 0 for rooted problems (root = {}).",
                                        target_num_active_clusters, graph_.root));
    }
    if (target_num_active_clusters < 0) {
        logger_->log(LogLevel::ERROR,
                     "Target number of active clusters ({}) cannot be negative.",
                     target_num_active_clusters);
        throw std::invalid_argument(std::format(
                                        "Target number of active clusters ({}) cannot be negative.",
                                        target_num_active_clusters));
    }

    if (graph_.prizes.empty()) {
        logger_->log(LogLevel::ERROR, "Prizes data cannot be empty.");
        throw std::invalid_argument("Prizes data cannot be empty.");
    }
    if (graph_.edges.size() != graph_.costs.size()) {
        logger_->log(LogLevel::ERROR,
                     "Number of edges ({}) does not match number of costs ({}).",
                     graph_.edges.size(), graph_.costs.size());
        throw std::invalid_argument(std::format(
                                        "Number of edges ({}) does not match number of costs ({}).",
                                        graph_.edges.size(), graph_.costs.size()));
    }

    for (size_t i = 0; i < num_nodes; ++i) {
        if (graph_.prizes[i] < 0.0) {
            logger_->log(LogLevel::ERROR, "Prize for node {} ({}) is negative.", i, graph_.prizes[i]);
            throw std::invalid_argument(std::format("Prize for node {} ({}) is negative.", i, graph_.prizes[i]));
        }
    }

    for (size_t i = 0; i < num_edges; ++i) {
        const double cost = graph_.costs[i];
        if (cost < 0.0) {
            logger_->log(LogLevel::ERROR, "Cost for edge {} ({}) is negative.", i, cost);
            throw std::invalid_argument(std::format("Cost for edge {} ({}) is negative.", i, cost));
        }
        
        const NodeId u = graph_.edges[i].first;
        const NodeId v = graph_.edges[i].second;
        if (u < 0 || static_cast<size_t>(u) >= num_nodes || v < 0 || static_cast<size_t>(v) >= num_nodes) {
            logger_->log(LogLevel::ERROR, "Edge {} ({}, {}) endpoint out of range [0, {}).", i, u, v, num_nodes);
            throw std::invalid_argument(std::format("Edge {} ({}, {}) endpoint out of range [0, {}).", i, u, v, num_nodes));
        }
    }

    logger_->log(LogLevel::INFO, "PCSTCoreAlgorithm initialized. Target clusters: {}.", target_num_active_clusters_);
}

PCSTCoreAlgorithm::~PCSTCoreAlgorithm() {}

void PCSTCoreAlgorithm::initialize() {
    current_time_ = 0.0;
    num_active_clusters_ = 0;
    stats_ = Statistics();
    phase1_result_edges_.clear();
    inactive_merge_events_.clear();

    const size_t num_nodes = graph_.prizes.size();
    const size_t num_edges = graph_.edges.size();
    const size_t max_clusters = num_nodes + (num_nodes > 0 ? num_nodes - 1 : 0);
    
    ph_nodes_.assign(2 * num_edges, PairingHeapNode{0.0, 0.0, -1, -1, -1});
    pairing_heap_buffer_.assign(2 * num_edges, 0);

    // Structure of Arrays setup for extreme cache coherence
    cluster_edge_parts_root_.assign(max_clusters, -1);
    cluster_active_.assign(max_clusters, 0);
    cluster_active_start_time_.assign(max_clusters, 0.0);
    cluster_active_end_time_.assign(max_clusters, -1.0);
    cluster_merged_into_.assign(max_clusters, kInvalidClusterId);
    cluster_prize_sum_.assign(max_clusters, 0.0);
    cluster_subcluster_moat_sum_.assign(max_clusters, 0.0);
    cluster_moat_.assign(max_clusters, 0.0);
    cluster_contains_root_.assign(max_clusters, 0);
    cluster_skip_up_.assign(max_clusters, SkipUpNode{kInvalidClusterId, 0.0});
    cluster_merged_along_.assign(max_clusters, kInvalidEdgeId);
    cluster_child_cluster_1_.assign(max_clusters, kInvalidClusterId);
    cluster_child_cluster_2_.assign(max_clusters, kInvalidClusterId);
    
    num_created_clusters_ = num_nodes;

    edge_parts_.assign(2 * num_edges, EdgePart{});
    edge_inactive_merge_event_ids_.assign(num_edges, kInvalidEventId);
    inactive_merge_events_.reserve(num_nodes);
    
    node_good_.assign(num_nodes, 0);

    clusters_deactivation_.reserve(max_clusters, max_clusters);
    clusters_next_edge_event_.reserve(max_clusters, max_clusters);

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes); ++i) {
        cluster_active_[i] = (i != graph_.root);
        cluster_active_end_time_[i] = (i == graph_.root) ? 0.0 : -1.0;
        cluster_prize_sum_[i] = graph_.prizes[i];
        cluster_contains_root_[i] = (i == graph_.root);

        if (cluster_active_[i]) {
            num_active_clusters_++;
            clusters_deactivation_.push_back_fast(cluster_prize_sum_[i], i);
        }
    }
    
    clusters_deactivation_.build_heap();

    for (int i = 0; i < static_cast<int>(num_edges); ++i) {
        const NodeId u = graph_.edges[i].first;
        const NodeId v = graph_.edges[i].second;
        const double cost = graph_.costs[i];

        EdgePart& u_part = edge_parts_[2 * i];
        EdgePart& v_part = edge_parts_[2 * i + 1];

        u_part.endpoint_node = u;
        v_part.endpoint_node = v;

        if (u == v) {
            u_part.deleted = true;
            v_part.deleted = true;
            continue;
        }

        bool u_active = (u != graph_.root);
        bool v_active = (v != graph_.root);

        if (u_active && v_active) {
            double event_val = cost * 0.5;
            u_part.next_event_val = event_val;
            v_part.next_event_val = event_val;
        } else if (u_active) {
            u_part.next_event_val = cost;
            v_part.next_event_val = 0.0;
        } else if (v_active) {
            u_part.next_event_val = 0.0;
            v_part.next_event_val = cost;
        } else {
            u_part.next_event_val = 0.0;
            v_part.next_event_val = 0.0;
        }

        u_part.in_heap = u_active;
        v_part.in_heap = v_active;

        if (u_active) {
            cluster_edge_parts_root_[u] = ph_insert(cluster_edge_parts_root_[u], 2 * i, u_part.next_event_val);
        }
        if (v_active) {
            cluster_edge_parts_root_[v] = ph_insert(cluster_edge_parts_root_[v], 2 * i + 1, v_part.next_event_val);
        }
    }

    for (ClusterId i = 0; i < static_cast<ClusterId>(num_nodes); ++i) {
        if (cluster_active_[i]) {
            double min_val;
            int32_t min_edge;
            if (ph_get_min(cluster_edge_parts_root_[i], &min_val, &min_edge)) {
                clusters_next_edge_event_.push_back_fast(min_val, i);
            }
        }
    }
    
    clusters_next_edge_event_.build_heap();
}

CoreAlgorithmResult PCSTCoreAlgorithm::run() {
    initialize();

    while (num_active_clusters_ > target_num_active_clusters_) {
        double edge_event_time = std::numeric_limits<double>::infinity();
        EdgePartId edge_part_idx = kInvalidEdgePartId;

        auto min_edge = clusters_next_edge_event_.get_min();
        if (min_edge) {
            ClusterId edge_cluster_idx = min_edge->second;
            ph_get_min(cluster_edge_parts_root_[edge_cluster_idx], &edge_event_time, &edge_part_idx);
        }

        double cluster_event_time = std::numeric_limits<double>::infinity();
        ClusterId cluster_idx = kInvalidClusterId;
        auto min_cluster = clusters_deactivation_.get_min();
        if (min_cluster) {
            cluster_event_time = min_cluster->first;
            cluster_idx = min_cluster->second;
        }

        if (edge_event_time == std::numeric_limits<double>::infinity() &&
            cluster_event_time == std::numeric_limits<double>::infinity()) {
            break;
        }

        if (edge_event_time <= cluster_event_time + eps_) {
            stats_.total_num_edge_events++;
            current_time_ = edge_event_time;

            ClusterId triggering_cluster_idx = min_edge->second;
            remove_next_edge_event(triggering_cluster_idx);
            handle_edge_event(current_time_, edge_part_idx);
        } else {
            stats_.num_cluster_events++;
            current_time_ = cluster_event_time;
            clusters_deactivation_.delete_min();
            handle_cluster_event(current_time_, cluster_idx);
        }
    }

    node_good_.assign(graph_.prizes.size(), 0);
    if (graph_.root != kInvalidNodeId) {
        ClusterId final_root_cluster = kInvalidClusterId;
        for(ClusterId i = 0; i < num_created_clusters_; ++i) {
            if (cluster_contains_root_[i] && cluster_merged_into_[i] == kInvalidClusterId) {
                final_root_cluster = i;
                break;
            }
        }
        if(final_root_cluster != kInvalidClusterId) mark_nodes_as_good(final_root_cluster);
        else if(graph_.root >= 0 && static_cast<size_t>(graph_.root) < node_good_.size()) node_good_[graph_.root] = 1;
    } else {
        for (ClusterId i = 0; i < num_created_clusters_; ++i) {
            if (cluster_active_[i] && cluster_merged_into_[i] == kInvalidClusterId) {
                mark_nodes_as_good(i);
            }
        }
    }

    return build_core_result();
}

void PCSTCoreAlgorithm::handle_edge_event(double event_time, EdgePartId edge_part_index) {
    if (edge_parts_[edge_part_index].deleted) {
        stats_.num_deleted_edge_events++;
        return;
    }

    EdgePartId other_edge_part_index = get_other_edge_part_index(edge_part_index);
    EdgeId edge_index = edge_part_index / 2;
    double current_edge_cost = graph_.costs[edge_index];

    double sum_current, finished_moat_current;
    ClusterId cluster_idx_current;
    get_sum_on_edge_part(edge_part_index, &sum_current, &finished_moat_current, &cluster_idx_current);

    double sum_other, finished_moat_other;
    ClusterId cluster_idx_other;
    get_sum_on_edge_part(other_edge_part_index, &sum_other, &finished_moat_other, &cluster_idx_other);

    if (cluster_idx_current == cluster_idx_other) {
        stats_.num_merged_edge_events++;
        edge_parts_[edge_part_index].deleted = true;
        edge_parts_[other_edge_part_index].deleted = true;
        return;
    }

    if (edge_parts_[other_edge_part_index].deleted) {
        stats_.num_deleted_edge_events++;
        edge_parts_[edge_part_index].deleted = true;
        return;
    }

    double remainder = current_edge_cost - sum_current - sum_other;

    if (remainder <= eps_ * current_edge_cost || std::fabs(remainder) < eps_) {
        stats_.total_num_merge_events++;
        phase1_result_edges_.push_back(edge_index);

        edge_parts_[other_edge_part_index].deleted = true;
        merge_clusters(cluster_idx_current, cluster_idx_other, edge_index, event_time, std::max(0.0, remainder));
    } else {
        bool current_cluster_active = cluster_active_[cluster_idx_current];
        bool other_cluster_active = cluster_active_[cluster_idx_other];
        EdgePart& current_edge_part_ref = edge_parts_[edge_part_index];
        EdgePart& other_edge_part_ref = edge_parts_[other_edge_part_index];

        if (current_cluster_active && other_cluster_active) {
            stats_.total_num_edge_growth_events++;
            stats_.num_active_active_edge_growth_events++;

            double time_to_meet = event_time + remainder * 0.5;
            double val_at_meet_current = sum_current + remainder * 0.5;
            double val_at_meet_other = sum_other + remainder * 0.5;

            current_edge_part_ref.next_event_val = val_at_meet_current;
            cluster_edge_parts_root_[cluster_idx_current] = ph_insert(
                cluster_edge_parts_root_[cluster_idx_current], edge_part_index, time_to_meet);
            current_edge_part_ref.in_heap = true;
            update_cluster_edge_pq(cluster_idx_current);

            double old_event_time_other = cluster_active_start_time_[cluster_idx_other] + other_edge_part_ref.next_event_val - finished_moat_other;

            if (other_edge_part_ref.in_heap) {
                cluster_edge_parts_root_[cluster_idx_other] = ph_decrease_key(
                    cluster_edge_parts_root_[cluster_idx_other], other_edge_part_index, old_event_time_other, time_to_meet);
                other_edge_part_ref.next_event_val = val_at_meet_other;
                update_cluster_edge_pq(cluster_idx_other);
            } else {
                other_edge_part_ref.next_event_val = val_at_meet_other;
            }
        } else {
            stats_.total_num_edge_growth_events++;
            stats_.num_active_inactive_edge_growth_events++;

            EdgePart& active_part = current_cluster_active ? current_edge_part_ref : other_edge_part_ref;
            EdgePart& inactive_part = current_cluster_active ? other_edge_part_ref : current_edge_part_ref;
            EdgePartId active_part_idx = current_cluster_active ? edge_part_index : other_edge_part_index;
            EdgePartId inactive_part_idx = current_cluster_active ? other_edge_part_index : edge_part_index;
            ClusterId active_cluster_idx = current_cluster_active ? cluster_idx_current : cluster_idx_other;
            ClusterId inactive_cluster_idx = current_cluster_active ? cluster_idx_other : cluster_idx_current;

            double finished_moat_inactive = current_cluster_active ? finished_moat_other : finished_moat_current;
            double time_to_meet = event_time + remainder;
            double val_at_meet_active = current_edge_cost - finished_moat_inactive;

            active_part.next_event_val = val_at_meet_active;
            cluster_edge_parts_root_[active_cluster_idx] = ph_insert(
                cluster_edge_parts_root_[active_cluster_idx], active_part_idx, time_to_meet);
            active_part.in_heap = true;
            update_cluster_edge_pq(active_cluster_idx);

            double inactive_deactivation_time = cluster_active_end_time_[inactive_cluster_idx];

            if (inactive_part.in_heap) {
                double old_event_time_inactive = inactive_deactivation_time + inactive_part.next_event_val - finished_moat_inactive;
                cluster_edge_parts_root_[inactive_cluster_idx] = ph_decrease_key(
                    cluster_edge_parts_root_[inactive_cluster_idx], inactive_part_idx, old_event_time_inactive, inactive_deactivation_time);
                inactive_part.next_event_val = finished_moat_inactive;
            } else {
                inactive_part.next_event_val = finished_moat_inactive;
            }
        }
    }
}

void PCSTCoreAlgorithm::handle_cluster_event(double event_time, ClusterId cluster_index) {
    cluster_active_[cluster_index] = 0;
    cluster_active_end_time_[cluster_index] = event_time;
    cluster_moat_[cluster_index] = event_time - cluster_active_start_time_[cluster_index];
    num_active_clusters_--;

    clusters_next_edge_event_.delete_element(cluster_index);
}

ClusterId PCSTCoreAlgorithm::merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder) {
    ClusterId new_cluster_idx = num_created_clusters_++;
    
    bool cluster1_active = cluster_active_[cluster1_idx];
    bool cluster2_active = cluster_active_[cluster2_idx];

    if (cluster1_active && cluster2_active) {
        stats_.num_active_active_merge_events++;
    } else {
        stats_.num_active_inactive_merge_events++;

        ClusterId active_original_cluster_idx = cluster1_active ? cluster1_idx : cluster2_idx;
        ClusterId inactive_original_cluster_idx = cluster1_active ? cluster2_idx : cluster1_idx;
        int32_t inactive_edge_parts_root = cluster_edge_parts_root_[inactive_original_cluster_idx];

        NodeId u_node = graph_.edges[merge_edge_idx].first;
        NodeId v_node = graph_.edges[merge_edge_idx].second;
        double temp_sum, temp_moat;
        ClusterId u_repr_cluster, v_repr_cluster;
        get_sum_on_edge_part(2 * merge_edge_idx, &temp_sum, &temp_moat, &u_repr_cluster);
        get_sum_on_edge_part(2 * merge_edge_idx + 1, &temp_sum, &temp_moat, &v_repr_cluster);

        NodeId active_node = (cluster1_active ? u_node : v_node);
        NodeId inactive_node = (cluster1_active ? v_node : u_node);
        
        if (u_repr_cluster == active_original_cluster_idx && v_repr_cluster == inactive_original_cluster_idx) {
            active_node = u_node; inactive_node = v_node;
        } else if (v_repr_cluster == active_original_cluster_idx && u_repr_cluster == inactive_original_cluster_idx) {
            active_node = v_node; inactive_node = u_node;
        }

        inactive_merge_events_.push_back({active_original_cluster_idx, inactive_original_cluster_idx, active_node, inactive_node});
        edge_inactive_merge_event_ids_[merge_edge_idx] = inactive_merge_events_.size() - 1;

        if (inactive_edge_parts_root != -1) {
            double inactive_deact = cluster_active_end_time_[inactive_original_cluster_idx];
            double time_diff = (event_time + remainder) - inactive_deact;
            if (time_diff > 0.0) {
                ph_add_to_heap(inactive_edge_parts_root, time_diff);
            }
        }
    }

    if (cluster1_active) {
        cluster_active_[cluster1_idx] = 0;
        cluster_active_end_time_[cluster1_idx] = event_time + remainder;
        cluster_moat_[cluster1_idx] = cluster_active_end_time_[cluster1_idx] - cluster_active_start_time_[cluster1_idx];
        clusters_deactivation_.delete_element(cluster1_idx);
        clusters_next_edge_event_.delete_element(cluster1_idx);
        num_active_clusters_--;
    }
    cluster_merged_into_[cluster1_idx] = new_cluster_idx;
    cluster_skip_up_[cluster1_idx] = SkipUpNode{new_cluster_idx, cluster_moat_[cluster1_idx]};

    if (cluster2_active) {
        cluster_active_[cluster2_idx] = 0;
        cluster_active_end_time_[cluster2_idx] = event_time + remainder;
        cluster_moat_[cluster2_idx] = cluster_active_end_time_[cluster2_idx] - cluster_active_start_time_[cluster2_idx];
        clusters_deactivation_.delete_element(cluster2_idx);
        clusters_next_edge_event_.delete_element(cluster2_idx);
        num_active_clusters_--;
    }
    cluster_merged_into_[cluster2_idx] = new_cluster_idx;
    cluster_skip_up_[cluster2_idx] = SkipUpNode{new_cluster_idx, cluster_moat_[cluster2_idx]};

    cluster_prize_sum_[new_cluster_idx] = cluster_prize_sum_[cluster1_idx] + cluster_prize_sum_[cluster2_idx];
    cluster_subcluster_moat_sum_[new_cluster_idx] = cluster_subcluster_moat_sum_[cluster1_idx] + cluster_subcluster_moat_sum_[cluster2_idx] + cluster_moat_[cluster1_idx] + cluster_moat_[cluster2_idx];
    cluster_contains_root_[new_cluster_idx] = cluster_contains_root_[cluster1_idx] || cluster_contains_root_[cluster2_idx];
    
    bool new_cluster_active = !cluster_contains_root_[new_cluster_idx];
    cluster_active_[new_cluster_idx] = new_cluster_active;

    cluster_merged_along_[new_cluster_idx] = merge_edge_idx;
    cluster_child_cluster_1_[new_cluster_idx] = cluster1_idx;
    cluster_child_cluster_2_[new_cluster_idx] = cluster2_idx;

    cluster_edge_parts_root_[new_cluster_idx] = ph_meld(cluster_edge_parts_root_[cluster1_idx], cluster_edge_parts_root_[cluster2_idx]);

    if (new_cluster_active) {
        cluster_active_start_time_[new_cluster_idx] = event_time + remainder;
        num_active_clusters_++;

        double potential_deactivation_time = cluster_active_start_time_[new_cluster_idx] + cluster_prize_sum_[new_cluster_idx] - cluster_subcluster_moat_sum_[new_cluster_idx];
        if (potential_deactivation_time < cluster_active_start_time_[new_cluster_idx] - eps_) {
            potential_deactivation_time = cluster_active_start_time_[new_cluster_idx];
        }

        clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);
        update_cluster_edge_pq(new_cluster_idx);
    }  else {
        cluster_active_end_time_[new_cluster_idx] = event_time + remainder;
    }

    return new_cluster_idx;
}

void PCSTCoreAlgorithm::remove_next_edge_event(ClusterId cluster_index) {
    double tmp_value = 0.0;
    int32_t tmp_edge_part = -1;
    
    cluster_edge_parts_root_[cluster_index] = ph_delete_min(cluster_edge_parts_root_[cluster_index], &tmp_value, &tmp_edge_part);
    
    if (tmp_edge_part != -1) {
        edge_parts_[tmp_edge_part].in_heap = false;
    }
    
    update_cluster_edge_pq(cluster_index);
}

FORCE_INLINE void PCSTCoreAlgorithm::get_sum_on_edge_part(EdgePartId edge_part_index,
    double* total_sum,
    double* finished_moat_sum,
    ClusterId* current_cluster_index) {

    ClusterId curr = edge_parts_[edge_part_index].endpoint_node;
    double total_sum_val = 0.0;

    while (true) {
        ClusterId parent = cluster_skip_up_[curr].skip_up;
        if (parent == kInvalidClusterId) {
            break;
        }
        
        double parent_moat = cluster_skip_up_[curr].sum;
        ClusterId grandparent = cluster_skip_up_[parent].skip_up;
        
        if (grandparent != kInvalidClusterId) {
            double grandparent_moat = cluster_skip_up_[parent].sum;
            cluster_skip_up_[curr] = SkipUpNode{grandparent, parent_moat + grandparent_moat};
            
            total_sum_val += parent_moat + grandparent_moat;
            curr = grandparent;
        } else {
            total_sum_val += parent_moat;
            curr = parent;
            break; // Break immediately to avoid re-reading the root's skip_up pointer
        }
    }

    if (cluster_active_[curr]) {
        *finished_moat_sum = total_sum_val;
        total_sum_val += current_time_ - cluster_active_start_time_[curr];
    } else {
        total_sum_val += cluster_moat_[curr];
        *finished_moat_sum = total_sum_val;
    }

    *total_sum = total_sum_val;
    *current_cluster_index = curr;
}

void PCSTCoreAlgorithm::mark_nodes_as_good(ClusterId start_cluster_index) {
    cluster_queue_.clear();
    cluster_queue_.push_back(start_cluster_index);

    std::vector<uint8_t> visited_clusters(num_created_clusters_, 0);
    visited_clusters[start_cluster_index] = 1;

    size_t queue_index = 0;
    while (queue_index < cluster_queue_.size()) {
        ClusterId current_cluster_idx = cluster_queue_[queue_index++];

        if (cluster_merged_along_[current_cluster_idx] == kInvalidEdgeId) {
            if (current_cluster_idx >= 0 && static_cast<size_t>(current_cluster_idx) < node_good_.size()) {
                node_good_[current_cluster_idx] = 1;
            }
        } else {
            if (cluster_child_cluster_1_[current_cluster_idx] != kInvalidClusterId && !visited_clusters[cluster_child_cluster_1_[current_cluster_idx]]) {
                visited_clusters[cluster_child_cluster_1_[current_cluster_idx]] = 1;
                cluster_queue_.push_back(cluster_child_cluster_1_[current_cluster_idx]);
            }
            if (cluster_child_cluster_2_[current_cluster_idx] != kInvalidClusterId && !visited_clusters[cluster_child_cluster_2_[current_cluster_idx]]) {
                visited_clusters[cluster_child_cluster_2_[current_cluster_idx]] = 1;
                cluster_queue_.push_back(cluster_child_cluster_2_[current_cluster_idx]);
            }
        }
    }
}

CoreAlgorithmResult PCSTCoreAlgorithm::build_core_result() {
    CoreAlgorithmResult result;
    result.statistics = stats_;
    result.phase1_edges = std::move(phase1_result_edges_);
    result.initial_node_filter = std::move(node_good_);
    result.edge_inactive_merge_event_ids = std::move(edge_inactive_merge_event_ids_);
    result.inactive_merge_events = std::move(inactive_merge_events_);
    
    std::vector<Cluster> final_state(num_created_clusters_);
    for(int i = 0; i < num_created_clusters_; ++i) {
        final_state[i].edge_parts_root = cluster_edge_parts_root_[i];
        final_state[i].active = cluster_active_[i];
        final_state[i].active_start_time = cluster_active_start_time_[i];
        final_state[i].active_end_time = cluster_active_end_time_[i];
        final_state[i].merged_into = cluster_merged_into_[i];
        final_state[i].prize_sum = cluster_prize_sum_[i];
        final_state[i].subcluster_moat_sum = cluster_subcluster_moat_sum_[i];
        final_state[i].moat = cluster_moat_[i];
        final_state[i].contains_root = cluster_contains_root_[i];
        final_state[i].skip_up = cluster_skip_up_[i].skip_up;
        final_state[i].skip_up_sum = cluster_skip_up_[i].sum;
        final_state[i].merged_along = cluster_merged_along_[i];
        final_state[i].child_cluster_1 = cluster_child_cluster_1_[i];
        final_state[i].child_cluster_2 = cluster_child_cluster_2_[i];
    }
    
    result.final_cluster_state = std::move(final_state);
    return result;
}

FORCE_INLINE int32_t PCSTCoreAlgorithm::ph_link(int32_t node1, int32_t node2) {
    if (node1 == -1) return node2;
    if (node2 == -1) return node1;

    auto* s_node = &ph_nodes_[node1];
    auto* l_node = &ph_nodes_[node2];
    int32_t smaller_node = node1;
    int32_t larger_node = node2;

    if (l_node->value < s_node->value) {
        std::swap(smaller_node, larger_node);
        std::swap(s_node, l_node);
    }

    l_node->sibling = s_node->child;
    if (s_node->child != -1) ph_nodes_[s_node->child].left_up = larger_node;
    l_node->left_up = smaller_node;
    s_node->child = larger_node;

    l_node->value -= s_node->child_offset;
    l_node->child_offset -= s_node->child_offset;

    return smaller_node;
}

FORCE_INLINE int32_t PCSTCoreAlgorithm::ph_insert(int32_t root, int32_t node, double value) {
    auto& n = ph_nodes_[node];
    n.value = value;
    n.child_offset = 0.0;
    n.sibling = -1;
    n.child = -1;
    n.left_up = -1;
    return ph_link(root, node);
}

FORCE_INLINE void PCSTCoreAlgorithm::ph_add_to_heap(int32_t root, double value) {
    if (root != -1 && value > 0.0) {
        ph_nodes_[root].value += value;
        ph_nodes_[root].child_offset += value;
    }
}

FORCE_INLINE int32_t PCSTCoreAlgorithm::ph_decrease_key(int32_t root, int32_t node, double from_value, double to_value) {
    auto& n = ph_nodes_[node];
    n.child_offset += (from_value - n.value);
    n.value = to_value;

    if (node == root) return root;

    if (n.left_up != -1) {
        int32_t p_or_l = n.left_up;
        auto& p = ph_nodes_[p_or_l];
        if (p.child == node) p.child = n.sibling;
        else p.sibling = n.sibling;

        if (n.sibling != -1) ph_nodes_[n.sibling].left_up = p_or_l;

        n.left_up = -1;
        n.sibling = -1;

        return ph_link(root, node);
    }
    return root;
}

FORCE_INLINE int32_t PCSTCoreAlgorithm::ph_delete_min(int32_t root, double* out_value, int32_t* out_node) {
    if (root == -1) return -1;

    auto& old_root = ph_nodes_[root];
    *out_value = old_root.value;
    *out_node = root;

    double root_offset = old_root.child_offset;

    // Fast-path for exactly 1 or 2 children
    if (old_root.child != -1) {
        int32_t child1 = old_root.child;
        int32_t child2 = ph_nodes_[child1].sibling;
        
        if (child2 == -1) {
            auto& c1_node = ph_nodes_[child1];
            c1_node.value += root_offset;
            c1_node.child_offset += root_offset;
            c1_node.left_up = -1;
            c1_node.sibling = -1;
            return child1;
        } else if (ph_nodes_[child2].sibling == -1) {
            auto& c1_node = ph_nodes_[child1];
            c1_node.value += root_offset;
            c1_node.child_offset += root_offset;
            c1_node.left_up = -1;
            c1_node.sibling = -1;

            auto& c2_node = ph_nodes_[child2];
            c2_node.value += root_offset;
            c2_node.child_offset += root_offset;
            c2_node.left_up = -1;
            c2_node.sibling = -1;

            return ph_link(child1, child2);
        }
    }

    // Generic two-pass merge for 3 or more children
    int ph_buf_sz = 0;
    int32_t current_child = old_root.child;
    
    while (current_child != -1) {
        auto& child_node = ph_nodes_[current_child];
        int32_t next_sibling = child_node.sibling;
        child_node.value += root_offset;
        child_node.child_offset += root_offset;
        child_node.left_up = -1;
        child_node.sibling = -1;
        pairing_heap_buffer_[ph_buf_sz++] = current_child;
        current_child = next_sibling;
    }

    if (ph_buf_sz == 0) return -1;

    int write_idx = 0;
    int merged = 0;
    
    while (merged + 1 < ph_buf_sz) {
        pairing_heap_buffer_[write_idx] = ph_link(pairing_heap_buffer_[merged], pairing_heap_buffer_[merged + 1]);
        merged += 2;
        write_idx++;
    }
    
    if (merged < ph_buf_sz) {
        pairing_heap_buffer_[write_idx++] = pairing_heap_buffer_[merged];
    }
    
    int32_t new_root = pairing_heap_buffer_[write_idx - 1];
    for (int i = write_idx - 2; i >= 0; --i) {
        new_root = ph_link(new_root, pairing_heap_buffer_[i]);
    }
    return new_root;
}

FORCE_INLINE int32_t PCSTCoreAlgorithm::ph_meld(int32_t root1, int32_t root2) {
    return ph_link(root1, root2);
}

FORCE_INLINE bool PCSTCoreAlgorithm::ph_get_min(int32_t root, double* out_val, int32_t* out_node) const {
    if (root != -1) {
        *out_val = ph_nodes_[root].value;
        *out_node = root;
        return true;
    }
    return false;
}

FORCE_INLINE void PCSTCoreAlgorithm::update_cluster_edge_pq(ClusterId cluster_idx) {
    double min_val;
    int32_t min_edge;
    if (ph_get_min(cluster_edge_parts_root_[cluster_idx], &min_val, &min_edge)) {
        clusters_next_edge_event_.insert_or_update(min_val, cluster_idx);
    } else {
        clusters_next_edge_event_.delete_element(cluster_idx);
    }
}

}