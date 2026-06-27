#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pcst_core_internals.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>
#include <format>

#ifdef _OPENMP
#include <omp.h>
#endif

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

PCSTCoreAlgorithm::~PCSTCoreAlgorithm() {
    logger_->log(LogLevel::DEBUG, "PCSTCoreAlgorithm destructor called.");
}

void PCSTCoreAlgorithm::initialize() {
    logger_->log(LogLevel::DEBUG, "Initializing core algorithm state.");
    current_time_ = 0.0;
    num_active_clusters_ = 0;
    stats_ = Statistics();
    phase1_result_edges_.clear();
    inactive_merge_events_.clear();
    pairing_heap_buffer_.clear();

    const size_t num_nodes = graph_.prizes.size();
    const size_t num_edges = graph_.edges.size();
    
    // Allocate single contiguous heap node pool
    heap_node_allocator_ = std::make_unique<PairingHeapType::AllocatorType>(2 * num_edges);

    clusters_.clear();
    const size_t max_clusters = num_nodes + (num_nodes > 0 ? num_nodes - 1 : 0);
    clusters_.reserve(max_clusters);

    // Resize parallel vectors and pre-allocate capacity to completely prevent reallocation in the loop
    cluster_merged_into_.assign(num_nodes, kInvalidClusterId);
    cluster_skip_up_.assign(num_nodes, kInvalidClusterId);
    cluster_skip_up_sum_.assign(num_nodes, 0.0);
    cluster_moat_.assign(num_nodes, 0.0);
    cluster_active_.assign(num_nodes, 0);
    cluster_active_start_time_.assign(num_nodes, 0.0);
    cluster_active_end_time_.assign(num_nodes, -1.0);

    cluster_merged_into_.reserve(max_clusters);
    cluster_skip_up_.reserve(max_clusters);
    cluster_skip_up_sum_.reserve(max_clusters);
    cluster_moat_.reserve(max_clusters);
    cluster_active_.reserve(max_clusters);
    cluster_active_start_time_.reserve(max_clusters);
    cluster_active_end_time_.reserve(max_clusters);

    edge_parts_.assign(2 * num_edges, EdgePart{});
    edge_info_.assign(num_edges, EdgeInfo{});
    node_good_.assign(num_nodes, false);

    clusters_deactivation_ = PriorityQueueType();
    clusters_next_edge_event_ = PriorityQueueType();

    // Parallel initialization of cluster objects and parallel vectors
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(num_nodes); ++i) {
        cluster_active_[i] = (i != graph_.root) ? 1 : 0;
        cluster_active_start_time_[i] = 0.0;
        cluster_active_end_time_[i] = (i == graph_.root) ? 0.0 : -1.0;
        cluster_merged_into_[i] = kInvalidClusterId;
        cluster_moat_[i] = 0.0;
        cluster_skip_up_[i] = kInvalidClusterId;
        cluster_skip_up_sum_[i] = 0.0;
    }

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes); ++i) {
        clusters_.emplace_back(heap_node_allocator_.get(), &pairing_heap_buffer_);
        Cluster& cluster = clusters_.back();

        cluster.active = (cluster_active_[i] != 0);
        cluster.active_start_time = cluster_active_start_time_[i];
        cluster.active_end_time = cluster_active_end_time_[i];
        cluster.merged_into = cluster_merged_into_[i];
        cluster.prize_sum = graph_.prizes[i];
        cluster.subcluster_moat_sum = 0.0;
        cluster.moat = cluster_moat_[i];
        cluster.contains_root = (i == graph_.root);
        cluster.skip_up = cluster_skip_up_[i];
        cluster.skip_up_sum = cluster_skip_up_sum_[i];
        cluster.merged_along = kInvalidEdgeId;
        cluster.child_cluster_1 = kInvalidClusterId;
        cluster.child_cluster_2 = kInvalidClusterId;
        cluster.necessary = false;

        if (cluster.active) {
            num_active_clusters_++;
            clusters_deactivation_.insert(cluster.prize_sum, i);
        }
    }

    // Parallel edge setup configuration
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(num_edges); ++i) {
        const NodeId u = graph_.edges[i].first;
        const NodeId v = graph_.edges[i].second;
        const double cost = graph_.costs[i];

        if (u == v) {
            edge_parts_[2 * i].deleted = true;
            edge_parts_[2 * i + 1].deleted = true;
            continue;
        }

        EdgePart& u_part = edge_parts_[2 * i];
        EdgePart& v_part = edge_parts_[2 * i + 1];

        u_part.deleted = false;
        v_part.deleted = false;

        bool u_active = (u != graph_.root);
        bool v_active = (v != graph_.root);

        if (u_active && v_active) {
            double event_val = cost / 2.0;
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
    }

    // Sequentially insert edge parts to preserve exact allocation/heap order
    for (EdgeId i = 0; i < static_cast<EdgeId>(num_edges); ++i) {
        const NodeId u = graph_.edges[i].first;
        const NodeId v = graph_.edges[i].second;
        if (u == v) continue;

        EdgePart& u_part = edge_parts_[2 * i];
        EdgePart& v_part = edge_parts_[2 * i + 1];

        if (cluster_active_[u]) {
            u_part.heap_node = clusters_[u].edge_parts.insert(u_part.next_event_val, 2 * i);
        } else {
            u_part.heap_node = -1;
        }
        if (cluster_active_[v]) {
            v_part.heap_node = clusters_[v].edge_parts.insert(v_part.next_event_val, 2 * i + 1);
        } else {
            v_part.heap_node = -1;
        }
    }

    for (ClusterId i = 0; i < static_cast<ClusterId>(num_nodes); ++i) {
        if (cluster_active_[i] && !clusters_[i].edge_parts.is_empty()) {
            double min_val;
            EdgePartId min_edge_part;
            [[maybe_unused]] bool success = clusters_[i].edge_parts.get_min(&min_val, &min_edge_part);
            assert(success);
            clusters_next_edge_event_.insert(min_val, i);
        }
    }
}

CoreAlgorithmResult PCSTCoreAlgorithm::run() {
    initialize();

    logger_->log(LogLevel::INFO, "Starting core algorithm run. Initial active clusters: {}", num_active_clusters_);

    while (num_active_clusters_ > target_num_active_clusters_) {
        auto next_edge_event = get_next_edge_event();
        auto next_cluster_event = get_next_cluster_event();

        double edge_event_time = std::numeric_limits<double>::infinity();
        EdgePartId edge_part_idx = kInvalidEdgePartId;

        if (next_edge_event) {
            edge_event_time = next_edge_event->first;
            edge_part_idx = next_edge_event->second.second;
        }

        double cluster_event_time = std::numeric_limits<double>::infinity();
        ClusterId cluster_idx = kInvalidClusterId;
        if (next_cluster_event) {
            cluster_event_time = next_cluster_event->first;
            cluster_idx = next_cluster_event->second;
        }

        if (edge_event_time == std::numeric_limits<double>::infinity() &&
                cluster_event_time == std::numeric_limits<double>::infinity()) {
            break;
        }

        double time_delta = std::min(edge_event_time, cluster_event_time) - current_time_;
        
        if (time_delta < -eps_) {
            throw std::runtime_error(std::format(
                                         "Negative time delta detected! Next event time {} < current time {}. Aborting.",
                                         std::min(edge_event_time, cluster_event_time), current_time_));
        }

        if (edge_event_time <= cluster_event_time + eps_) {
            stats_.total_num_edge_events++;
            current_time_ = edge_event_time;

            assert(next_edge_event);
            ClusterId triggering_cluster_idx = next_edge_event->second.first;
            remove_next_edge_event(triggering_cluster_idx);
            handle_edge_event(current_time_, edge_part_idx);
        } else {
            stats_.num_cluster_events++;
            current_time_ = cluster_event_time;
            assert(next_cluster_event);
            remove_next_cluster_event();
            handle_cluster_event(current_time_, cluster_idx);
        }
    }

    node_good_.assign(graph_.prizes.size(), false);
    if (graph_.root != kInvalidNodeId) {
        ClusterId final_root_cluster = kInvalidClusterId;
        for(ClusterId i = 0; i < static_cast<ClusterId>(clusters_.size()); ++i) {
            if (clusters_[i].contains_root && cluster_merged_into_[i] == kInvalidClusterId) {
                final_root_cluster = i;
                break;
            }
        }
        if(final_root_cluster != kInvalidClusterId) {
            mark_nodes_as_good(final_root_cluster);
        } else {
            if(graph_.root >= 0 && static_cast<size_t>(graph_.root) < node_good_.size()) {
                node_good_[graph_.root] = true;
            }
        }
    } else {
        for (ClusterId i = 0; i < static_cast<ClusterId>(clusters_.size()); ++i) {
            if (cluster_active_[i] && cluster_merged_into_[i] == kInvalidClusterId) {
                mark_nodes_as_good(i);
            }
        }
    }

    return build_core_result();
}

void PCSTCoreAlgorithm::handle_edge_event(double event_time, EdgePartId edge_part_index) {
    assert(static_cast<size_t>(edge_part_index) < edge_parts_.size());

    if (edge_parts_[edge_part_index].deleted) {
        stats_.num_deleted_edge_events++;
        return;
    }

    EdgePartId other_edge_part_index = get_other_edge_part_index(edge_part_index);
    assert(static_cast<size_t>(other_edge_part_index) < edge_parts_.size());
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
        bool current_cluster_active = (cluster_active_[cluster_idx_current] != 0);
        bool other_cluster_active = (cluster_active_[cluster_idx_other] != 0);
        Cluster& current_cluster = clusters_[cluster_idx_current];
        Cluster& other_cluster = clusters_[cluster_idx_other];
        EdgePart& current_edge_part_ref = edge_parts_[edge_part_index];
        EdgePart& other_edge_part_ref = edge_parts_[other_edge_part_index];

        if (current_cluster_active && other_cluster_active) {
            stats_.total_num_edge_growth_events++;
            stats_.num_active_active_edge_growth_events++;

            assert(remainder > 0.0 && "Remainder should be positive here.");
            double time_to_meet = event_time + remainder / 2.0;
            double val_at_meet_current = sum_current + remainder / 2.0;
            double val_at_meet_other = sum_other + remainder / 2.0;

            current_edge_part_ref.next_event_val = val_at_meet_current;
            current_edge_part_ref.heap_node = current_cluster.edge_parts.insert(time_to_meet, edge_part_index);

            if (!current_cluster.edge_parts.is_empty()) {
                double min_val;
                EdgePartId min_part;
                [[maybe_unused]] bool success = current_cluster.edge_parts.get_min(&min_val, &min_part);
                assert(success);
                clusters_next_edge_event_.insert(min_val, cluster_idx_current);
            }

            double old_event_time_other = cluster_active_start_time_[cluster_idx_other] + other_edge_part_ref.next_event_val - finished_moat_other;

            if (other_edge_part_ref.heap_node != -1) {
                clusters_next_edge_event_.delete_element(cluster_idx_other);
                other_cluster.edge_parts.decrease_key(other_edge_part_ref.heap_node, old_event_time_other, time_to_meet);
                other_edge_part_ref.next_event_val = val_at_meet_other;

                if (!other_cluster.edge_parts.is_empty()) {
                    double min_val;
                    EdgePartId min_part;
                    [[maybe_unused]] bool success = other_cluster.edge_parts.get_min(&min_val, &min_part);
                    assert(success);
                    clusters_next_edge_event_.insert(min_val, cluster_idx_other);
                }
            } else {
                other_edge_part_ref.next_event_val = val_at_meet_other;
            }
        } else {
            assert(current_cluster_active != other_cluster_active);
            assert(remainder > 0.0);
            stats_.total_num_edge_growth_events++;
            stats_.num_active_inactive_edge_growth_events++;

            Cluster& active_cluster = current_cluster_active ? current_cluster : other_cluster;
            Cluster& inactive_cluster = current_cluster_active ? other_cluster : current_cluster;
            EdgePart& active_part = current_cluster_active ? current_edge_part_ref : other_edge_part_ref;
            EdgePart& inactive_part = current_cluster_active ? other_edge_part_ref : current_edge_part_ref;
            EdgePartId active_part_idx = current_cluster_active ? edge_part_index : other_edge_part_index;
            //EdgePartId inactive_part_idx = current_cluster_active ? other_edge_part_index : edge_part_index;
            ClusterId active_cluster_idx = current_cluster_active ? cluster_idx_current : cluster_idx_other;
            ClusterId inactive_cluster_idx = current_cluster_active ? cluster_idx_other : cluster_idx_current;

            double finished_moat_inactive = current_cluster_active ? finished_moat_other : finished_moat_current;
            double time_to_meet = event_time + remainder;
            double val_at_meet_active = current_edge_cost - finished_moat_inactive;

            active_part.next_event_val = val_at_meet_active;
            active_part.heap_node = active_cluster.edge_parts.insert(time_to_meet, active_part_idx);

            if (!active_cluster.edge_parts.is_empty()) {
                double min_val;
                EdgePartId min_part;
                [[maybe_unused]] bool success = active_cluster.edge_parts.get_min(&min_val, &min_part);
                assert(success);
                clusters_next_edge_event_.insert(min_val, active_cluster_idx);
            }

            double inactive_deactivation_time = cluster_active_end_time_[inactive_cluster_idx];
            assert(inactive_deactivation_time >= 0.0 && "Inactive cluster must have a valid end time.");

            if (inactive_part.heap_node != -1) {
                double old_event_time_inactive = inactive_deactivation_time + inactive_part.next_event_val - finished_moat_inactive;
                inactive_cluster.edge_parts.decrease_key(inactive_part.heap_node, old_event_time_inactive, inactive_deactivation_time);
                inactive_part.next_event_val = finished_moat_inactive;
            } else {
                inactive_part.next_event_val = finished_moat_inactive;
            }
        }
    }
}

void PCSTCoreAlgorithm::handle_cluster_event(double event_time, ClusterId cluster_index) {
    assert(static_cast<size_t>(cluster_index) < clusters_.size());
    assert(cluster_active_[cluster_index] && "Cluster deactivation event for an already inactive cluster!");
    if (!cluster_active_[cluster_index]) return;

    cluster_active_[cluster_index] = 0;
    cluster_active_end_time_[cluster_index] = event_time;
    assert(event_time >= cluster_active_start_time_[cluster_index] && "Deactivation time cannot be before start time.");
    cluster_moat_[cluster_index] = cluster_active_end_time_[cluster_index] - cluster_active_start_time_[cluster_index];
    num_active_clusters_--;

    if (!clusters_[cluster_index].edge_parts.is_empty()) {
        clusters_next_edge_event_.delete_element(cluster_index);
    }
}

ClusterId PCSTCoreAlgorithm::merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder) {
    assert(static_cast<size_t>(cluster1_idx) < clusters_.size());
    assert(static_cast<size_t>(cluster2_idx) < clusters_.size());
    assert(cluster1_idx != cluster2_idx);

    clusters_.emplace_back(heap_node_allocator_.get(), &pairing_heap_buffer_);
    ClusterId new_cluster_idx = clusters_.size() - 1;

    // Maintain flat parallel vector boundaries
    cluster_merged_into_.push_back(kInvalidClusterId);
    cluster_skip_up_.push_back(kInvalidClusterId);
    cluster_skip_up_sum_.push_back(0.0);
    cluster_moat_.push_back(0.0);
    cluster_active_.push_back(0);
    cluster_active_start_time_.push_back(0.0);
    cluster_active_end_time_.push_back(-1.0);

    Cluster& cluster1 = clusters_[cluster1_idx];
    Cluster& cluster2 = clusters_[cluster2_idx];
    Cluster& new_cluster = clusters_[new_cluster_idx];

    bool cluster1_active = (cluster_active_[cluster1_idx] != 0);
    bool cluster2_active = (cluster_active_[cluster2_idx] != 0);

    if (cluster1_active && cluster2_active) {
        stats_.num_active_active_merge_events++;
    } else {
        assert((cluster1_active != cluster2_active) && "Cannot merge two inactive clusters.");
        stats_.num_active_inactive_merge_events++;

        Cluster& inactive_cluster_ref = cluster1_active ? cluster2 : cluster1;
        ClusterId active_original_cluster_idx = cluster1_active ? cluster1_idx : cluster2_idx;
        ClusterId inactive_original_cluster_idx = cluster1_active ? cluster2_idx : cluster1_idx;

        NodeId u_node = graph_.edges[merge_edge_idx].first;
        NodeId v_node = graph_.edges[merge_edge_idx].second;
        double temp_sum, temp_moat;
        ClusterId u_repr_cluster, v_repr_cluster;
        get_sum_on_edge_part(2 * merge_edge_idx, &temp_sum, &temp_moat, &u_repr_cluster);
        get_sum_on_edge_part(2 * merge_edge_idx + 1, &temp_sum, &temp_moat, &v_repr_cluster);

        NodeId active_node = kInvalidNodeId;
        NodeId inactive_node = kInvalidNodeId;

        if (u_repr_cluster == active_original_cluster_idx && v_repr_cluster == inactive_original_cluster_idx) {
            active_node = u_node;
            inactive_node = v_node;
        } else if (v_repr_cluster == active_original_cluster_idx && u_repr_cluster == inactive_original_cluster_idx) {
            active_node = v_node;
            inactive_node = u_node;
        } else {
            active_node = (cluster1_active ? u_node : v_node);
            inactive_node = (cluster1_active ? v_node : u_node);
            assert(false && "Logic error determining active/inactive nodes in merge.");
        }

        inactive_merge_events_.push_back({active_original_cluster_idx, inactive_original_cluster_idx, active_node, inactive_node});
        edge_info_[merge_edge_idx].inactive_merge_event = inactive_merge_events_.size() - 1;

        if (!inactive_cluster_ref.edge_parts.is_empty()) {
            double inactive_deact = cluster_active_end_time_[inactive_original_cluster_idx];
            assert(inactive_deact >= 0.0);
            double time_diff = (event_time + remainder) - inactive_deact;

            if (time_diff < -eps_) {
                time_diff = 0.0;
            }
            inactive_cluster_ref.edge_parts.add_to_heap(std::max(0.0, time_diff));
        }
    }

    if (cluster1_active) {
        cluster_active_[cluster1_idx] = 0;
        cluster_active_end_time_[cluster1_idx] = event_time + remainder;
        assert(cluster_active_end_time_[cluster1_idx] >= cluster_active_start_time_[cluster1_idx]);
        cluster_moat_[cluster1_idx] = cluster_active_end_time_[cluster1_idx] - cluster_active_start_time_[cluster1_idx];
        clusters_deactivation_.delete_element(cluster1_idx);
        if (!cluster1.edge_parts.is_empty()) {
            clusters_next_edge_event_.delete_element(cluster1_idx);
        }
        num_active_clusters_--;
    }
    cluster_merged_into_[cluster1_idx] = new_cluster_idx;

    if (cluster2_active) {
        cluster_active_[cluster2_idx] = 0;
        cluster_active_end_time_[cluster2_idx] = event_time + remainder;
        assert(cluster_active_end_time_[cluster2_idx] >= cluster_active_start_time_[cluster2_idx]);
        cluster_moat_[cluster2_idx] = cluster_active_end_time_[cluster2_idx] - cluster_active_start_time_[cluster2_idx];
        clusters_deactivation_.delete_element(cluster2_idx);
        if (!cluster2.edge_parts.is_empty()) {
            clusters_next_edge_event_.delete_element(cluster2_idx);
        }
        num_active_clusters_--;
    }
    cluster_merged_into_[cluster2_idx] = new_cluster_idx;

    new_cluster.prize_sum = cluster1.prize_sum + cluster2.prize_sum;
    new_cluster.subcluster_moat_sum = cluster1.subcluster_moat_sum + cluster2.subcluster_moat_sum
                                      + cluster_moat_[cluster1_idx] + cluster_moat_[cluster2_idx];
    new_cluster.contains_root = cluster1.contains_root || cluster2.contains_root;
    
    bool new_cluster_active = !new_cluster.contains_root;
    cluster_active_[new_cluster_idx] = new_cluster_active ? 1 : 0;

    new_cluster.merged_along = merge_edge_idx;
    new_cluster.child_cluster_1 = cluster1_idx;
    new_cluster.child_cluster_2 = cluster2_idx;

    new_cluster.merged_into = kInvalidClusterId;
    new_cluster.moat = 0.0;
    new_cluster.skip_up = kInvalidClusterId;
    new_cluster.skip_up_sum = 0.0;
    new_cluster.active_end_time = -1.0;

    new_cluster.edge_parts = PairingHeapType::meld(&cluster1.edge_parts, &cluster2.edge_parts);

    if (new_cluster_active) {
        cluster_active_start_time_[new_cluster_idx] = event_time + remainder;
        num_active_clusters_++;

        double potential_deactivation_time = cluster_active_start_time_[new_cluster_idx]
                                             + new_cluster.prize_sum
                                             - new_cluster.subcluster_moat_sum;

        if (potential_deactivation_time < cluster_active_start_time_[new_cluster_idx] - eps_) {
            potential_deactivation_time = cluster_active_start_time_[new_cluster_idx];
        }

        clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);

        if (!new_cluster.edge_parts.is_empty()) {
            double min_val;
            EdgePartId min_part;
            [[maybe_unused]] bool success = new_cluster.edge_parts.get_min(&min_val, &min_part);
            assert(success);
            clusters_next_edge_event_.insert(min_val, new_cluster_idx);
        }
    }

    return new_cluster_idx;
}

std::optional<std::pair<double, std::pair<ClusterId, EdgePartId>>> PCSTCoreAlgorithm::get_next_edge_event() {
    auto min_cluster_event = clusters_next_edge_event_.get_min();
    if (!min_cluster_event) {
        return std::nullopt;
    }

    //double global_event_time = min_cluster_event->first;
    ClusterId cluster_index = min_cluster_event->second;

    assert(static_cast<size_t>(cluster_index) < clusters_.size());

    if (clusters_[cluster_index].edge_parts.is_empty()) {
        clusters_next_edge_event_.delete_element(cluster_index);

        while(true) {
            min_cluster_event = clusters_next_edge_event_.get_min();
            if (!min_cluster_event) {
                return std::nullopt;
            }
            //global_event_time = min_cluster_event->first;
            cluster_index = min_cluster_event->second;
            assert(static_cast<size_t>(cluster_index) < clusters_.size());
            if (!clusters_[cluster_index].edge_parts.is_empty()) {
                break;
            }
            clusters_next_edge_event_.delete_element(cluster_index);
        }
    }

    double actual_heap_min_val = std::numeric_limits<double>::infinity();
    EdgePartId edge_part_index = kInvalidEdgePartId;

    bool success = clusters_[cluster_index].edge_parts.get_min(&actual_heap_min_val, &edge_part_index);

    if (!success) {
        throw std::runtime_error(std::format("Internal Error: Failed get_min for cluster {}", cluster_index));
    }

    double event_time = actual_heap_min_val;
    return std::make_pair(event_time, std::make_pair(cluster_index, edge_part_index));
}

void PCSTCoreAlgorithm::remove_next_edge_event(ClusterId cluster_index) {
    assert(static_cast<size_t>(cluster_index) < clusters_.size());

    clusters_next_edge_event_.delete_element(cluster_index);

    double tmp_value;
    EdgePartId tmp_edge_part;
    [[maybe_unused]] bool deleted = clusters_[cluster_index].edge_parts.delete_min(&tmp_value, &tmp_edge_part);
    assert(deleted);

    if (!clusters_[cluster_index].edge_parts.is_empty()) {
        [[maybe_unused]] bool success = clusters_[cluster_index].edge_parts.get_min(&tmp_value, &tmp_edge_part);
        assert(success);
        clusters_next_edge_event_.insert(tmp_value, cluster_index);
    }
}

std::optional<std::pair<double, ClusterId>> PCSTCoreAlgorithm::get_next_cluster_event() {
    return clusters_deactivation_.get_min();
}

void PCSTCoreAlgorithm::remove_next_cluster_event() {
    [[maybe_unused]] auto deleted_event = clusters_deactivation_.delete_min();
    assert(deleted_event.has_value());
}

/**
 * @brief Zero-allocation path compression utilizing a two-pass iterative traversal.
 */
void PCSTCoreAlgorithm::get_sum_on_edge_part(EdgePartId edge_part_index,
        double* total_sum,
        double* finished_moat_sum,
        ClusterId* current_cluster_index) {
    assert(total_sum != nullptr && finished_moat_sum != nullptr && current_cluster_index != nullptr);
    assert(static_cast<size_t>(edge_part_index / 2) < graph_.edges.size());

    NodeId endpoint_node = (edge_part_index % 2 == 0)
                           ? graph_.edges[edge_part_index / 2].first
                           : graph_.edges[edge_part_index / 2].second;

    assert(static_cast<size_t>(endpoint_node) < clusters_.size());

    double total_sum_val = 0.0;
    ClusterId curr = endpoint_node;

    // Pass 1: Trace to active root and calculate total intermediate moat sums
    while (cluster_merged_into_[curr] != kInvalidClusterId) {
        if (cluster_skip_up_[curr] != kInvalidClusterId) {
            total_sum_val += cluster_skip_up_sum_[curr];
            curr = cluster_skip_up_[curr];
        } else {
            total_sum_val += cluster_moat_[curr];
            curr = cluster_merged_into_[curr];
        }
        assert(static_cast<size_t>(curr) < clusters_.size());
    }

    ClusterId root = curr;

    // Pass 2: Re-traverse and apply zero-allocation iterative path compression
    double sum_before_x = 0.0;
    curr = endpoint_node;
    while (curr != root) {
        ClusterId next_node;
        double edge_moat;
        if (cluster_skip_up_[curr] != kInvalidClusterId) {
            next_node = cluster_skip_up_[curr];
            edge_moat = cluster_skip_up_sum_[curr];
        } else {
            next_node = cluster_merged_into_[curr];
            edge_moat = cluster_moat_[curr];
        }

        cluster_skip_up_[curr] = root;
        cluster_skip_up_sum_[curr] = total_sum_val - sum_before_x;

        sum_before_x += edge_moat;
        curr = next_node;
    }

    if (cluster_active_[root]) {
        *finished_moat_sum = total_sum_val;
        total_sum_val += current_time_ - cluster_active_start_time_[root];
    } else {
        total_sum_val += cluster_moat_[root];
        *finished_moat_sum = total_sum_val;
    }

    *total_sum = total_sum_val;
    *current_cluster_index = root;
}

void PCSTCoreAlgorithm::mark_nodes_as_good(ClusterId start_cluster_index) {
    assert(static_cast<size_t>(start_cluster_index) < clusters_.size());

    cluster_queue_.clear();
    cluster_queue_.push_back(start_cluster_index);

    std::vector<bool> visited_clusters(clusters_.size(), false);
    visited_clusters[start_cluster_index] = true;

    size_t queue_index = 0;
    while (queue_index < cluster_queue_.size()) {
        ClusterId current_cluster_idx = cluster_queue_[queue_index++];
        const Cluster& cluster = clusters_[current_cluster_idx];

        if (cluster.merged_along == kInvalidEdgeId) {
            assert(current_cluster_idx < static_cast<ClusterId>(graph_.prizes.size()));
            if (current_cluster_idx >= 0 && static_cast<size_t>(current_cluster_idx) < node_good_.size()) {
                if (!node_good_[current_cluster_idx]) {
                    node_good_[current_cluster_idx] = true;
                }
            }
        } else {
            if (cluster.child_cluster_1 != kInvalidClusterId && !visited_clusters[cluster.child_cluster_1]) {
                visited_clusters[cluster.child_cluster_1] = true;
                cluster_queue_.push_back(cluster.child_cluster_1);
            }
            if (cluster.child_cluster_2 != kInvalidClusterId && !visited_clusters[cluster.child_cluster_2]) {
                visited_clusters[cluster.child_cluster_2] = true;
                cluster_queue_.push_back(cluster.child_cluster_2);
            }
        }
    }
}

/**
 * @brief Constructs the output results, reconstructing structural Cluster components from fast flat parallel vectors.
 */
CoreAlgorithmResult PCSTCoreAlgorithm::build_core_result() {
    logger_->log(LogLevel::INFO, "Building core algorithm result.");
    CoreAlgorithmResult result;
    result.statistics = stats_;

    result.phase1_edges = phase1_result_edges_;
    logger_->log(LogLevel::DEBUG, "Phase 1 selected edges (unfiltered): {}.", result.phase1_edges.size());

    result.initial_node_filter = node_good_;

    result.edge_inactive_merge_event_ids.resize(edge_info_.size());
    for(size_t i=0; i < edge_info_.size(); ++i) {
        result.edge_inactive_merge_event_ids[i] = edge_info_[i].inactive_merge_event;
    }

    result.inactive_merge_events = std::move(inactive_merge_events_);
    
    // Copy the hot parameters into the Cluster array for downstream pruners
    size_t num_clusters_total = clusters_.size();
    for (size_t i = 0; i < num_clusters_total; ++i) {
        clusters_[i].active = (cluster_active_[i] != 0);
        clusters_[i].active_start_time = cluster_active_start_time_[i];
        clusters_[i].active_end_time = cluster_active_end_time_[i];
        clusters_[i].merged_into = cluster_merged_into_[i];
        clusters_[i].moat = cluster_moat_[i];
        clusters_[i].skip_up = cluster_skip_up_[i];
        clusters_[i].skip_up_sum = cluster_skip_up_sum_[i];
    }

    result.final_cluster_state = std::move(clusters_);
    result.heap_node_allocator = std::move(heap_node_allocator_);
    
    return result;
}

}