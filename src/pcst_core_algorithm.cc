#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pcst_core_internals.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>

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

    clusters_.clear();

    clusters_.reserve(num_nodes + (num_nodes > 0 ? num_nodes - 1 : 0));
    edge_parts_.assign(2 * num_edges, EdgePart{});
    edge_info_.assign(num_edges, EdgeInfo{});
    node_good_.assign(num_nodes, false);

    clusters_deactivation_ = PriorityQueueType();
    clusters_next_edge_event_ = PriorityQueueType();

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes); ++i) {

        clusters_.emplace_back(&pairing_heap_buffer_);
        Cluster& cluster = clusters_.back();

        cluster.active = (i != graph_.root);
        cluster.active_start_time = 0.0;
        cluster.active_end_time = (i == graph_.root) ? 0.0 : -1.0;
        cluster.merged_into = kInvalidClusterId;
        cluster.prize_sum = graph_.prizes[i];
        cluster.subcluster_moat_sum = 0.0;
        cluster.moat = 0.0;
        cluster.contains_root = (i == graph_.root);
        cluster.skip_up = kInvalidClusterId;
        cluster.skip_up_sum = 0.0;
        cluster.merged_along = kInvalidEdgeId;
        cluster.child_cluster_1 = kInvalidClusterId;
        cluster.child_cluster_2 = kInvalidClusterId;
        cluster.necessary = false;

        if (cluster.active) {
            num_active_clusters_++;
            clusters_deactivation_.insert(cluster.prize_sum, i);
            logger_->log(LogLevel::TRACE, "Node {} initialized as active cluster {}. Prize: {}. Added to deactivation queue.", i, i, cluster.prize_sum);
        } else {
            logger_->log(LogLevel::TRACE, "Node {} initialized as inactive cluster {}.", i, i);
        }
    }
    assert(static_cast<int>(clusters_.size()) == static_cast<int>(num_nodes));

    for (EdgeId i = 0; i < static_cast<EdgeId>(num_edges); ++i) {
        const NodeId u = graph_.edges[i].first;
        const NodeId v = graph_.edges[i].second;
        const double cost = graph_.costs[i];

        if (u == v) {
            logger_->log(LogLevel::WARNING, "Ignoring self-loop edge {} ({}, {}) with cost {}.", i, u, v, cost);

            if (static_cast<size_t>(2*i + 1) < edge_parts_.size()) {
                edge_parts_[2 * i].deleted = true;
                edge_parts_[2 * i + 1].deleted = true;
            }
            continue;
        }

        assert(static_cast<size_t>(u) < clusters_.size() && static_cast<size_t>(v) < clusters_.size());

        EdgePart& u_part = edge_parts_[2 * i];
        EdgePart& v_part = edge_parts_[2 * i + 1];
        Cluster& u_cluster = clusters_[u];
        Cluster& v_cluster = clusters_[v];

        u_part.deleted = false;
        v_part.deleted = false;

        if (u_cluster.active && v_cluster.active) {
            double event_val = cost / 2.0;
            u_part.next_event_val = event_val;
            v_part.next_event_val = event_val;
            logger_->log(LogLevel::TRACE, "Edge {}({},{}): Active-Active. Initial event val: {}", i, u, v, event_val);
        } else if (u_cluster.active) {
            u_part.next_event_val = cost;
            v_part.next_event_val = 0.0;
            logger_->log(LogLevel::TRACE, "Edge {}({},{}): Active(u)-Inactive(v). Initial event val u: {}, v: {}", i, u, v, cost, 0.0);
        } else if (v_cluster.active) {
            u_part.next_event_val = 0.0;
            v_part.next_event_val = cost;
            logger_->log(LogLevel::TRACE, "Edge {}({},{}): Inactive(u)-Active(v). Initial event val u: {}, v: {}", i, u, v, 0.0, cost);
        } else {
            u_part.next_event_val = 0.0;
            v_part.next_event_val = 0.0;
            logger_->log(LogLevel::TRACE, "Edge {}({},{}): Inactive-Inactive. Initial event val u: {}, v: {}", i, u, v, 0.0, 0.0);
        }

        if (u_cluster.active) {
            u_part.heap_node = u_cluster.edge_parts.insert(u_part.next_event_val, 2 * i);
        } else {
            u_part.heap_node = nullptr;
        }
        if (v_cluster.active) {
            v_part.heap_node = v_cluster.edge_parts.insert(v_part.next_event_val, 2 * i + 1);
        } else {
            v_part.heap_node = nullptr;
        }
    }

    for (ClusterId i = 0; i < static_cast<ClusterId>(num_nodes); ++i) {
        if (clusters_[i].active && !clusters_[i].edge_parts.is_empty()) {
            double min_val;
            EdgePartId min_edge_part;
            [[maybe_unused]] bool success = clusters_[i].edge_parts.get_min(&min_val, &min_edge_part);
            assert(success);
            clusters_next_edge_event_.insert(min_val, i);
            logger_->log(LogLevel::TRACE, "Cluster {}: Initial min edge event at time {} from part {}.", i, min_val, min_edge_part);
        }
    }
    logger_->log(LogLevel::INFO, "Initialization complete. {} nodes, {} edges. {} active clusters.", num_nodes, num_edges, num_active_clusters_);
}

CoreAlgorithmResult PCSTCoreAlgorithm::run() {
    initialize();

    logger_->log(LogLevel::INFO, "Starting core algorithm run. Initial active clusters: {}", num_active_clusters_);

    while (num_active_clusters_ > target_num_active_clusters_) {
        logger_->log(LogLevel::TRACE, "-----------------------------------------");
        logger_->log(LogLevel::DEBUG, "Start main loop iteration. Current time: {}, Active clusters: {}", current_time_, num_active_clusters_);

        auto next_edge_event = get_next_edge_event();
        auto next_cluster_event = get_next_cluster_event();

        double edge_event_time = std::numeric_limits<double>::infinity();
        EdgePartId edge_part_idx = kInvalidEdgePartId;

        if (next_edge_event) {
            edge_event_time = next_edge_event->first;
            edge_part_idx = next_edge_event->second.second;
            logger_->log(LogLevel::TRACE, "Next edge event: Time={}, Cluster={}, Part={}",
                         edge_event_time, next_edge_event->second.first, edge_part_idx);
        } else {
            logger_->log(LogLevel::TRACE, "No further edge events.");
        }

        double cluster_event_time = std::numeric_limits<double>::infinity();
        ClusterId cluster_idx = kInvalidClusterId;
        if (next_cluster_event) {
            cluster_event_time = next_cluster_event->first;
            cluster_idx = next_cluster_event->second;
            logger_->log(LogLevel::TRACE, "Next cluster event: Time={}, Cluster={}", cluster_event_time, cluster_idx);
        } else {
            logger_->log(LogLevel::TRACE, "No further cluster events.");
        }

        if (edge_event_time == std::numeric_limits<double>::infinity() &&
                cluster_event_time == std::numeric_limits<double>::infinity()) {
            logger_->log(LogLevel::WARNING, "No more events, but target active clusters ({}) not reached ({} remaining). Stopping early.",
                         target_num_active_clusters_, num_active_clusters_);
            break;
        }

        double time_delta = std::min(edge_event_time, cluster_event_time) - current_time_;
        if (time_delta < -eps_) {
            logger_->log(LogLevel::ERROR,
                         "Negative time delta detected! Next event time {} < current time {}. Aborting.",
                         std::min(edge_event_time, cluster_event_time), current_time_);
            throw std::runtime_error(std::format(
                                         "Negative time delta detected! Next event time {} < current time {}. Aborting.",
                                         std::min(edge_event_time, cluster_event_time), current_time_));
        }

        if (edge_event_time <= cluster_event_time + eps_) {
            stats_.total_num_edge_events++;
            current_time_ = edge_event_time;

            assert(next_edge_event);
            ClusterId triggering_cluster_idx = next_edge_event->second.first;
            logger_->log(LogLevel::DEBUG, "Processing edge event for part {} from cluster {} at time {}", edge_part_idx, triggering_cluster_idx, current_time_);
            remove_next_edge_event(triggering_cluster_idx);
            handle_edge_event(current_time_, edge_part_idx);
        } else {
            stats_.num_cluster_events++;
            current_time_ = cluster_event_time;
            assert(next_cluster_event);
            logger_->log(LogLevel::DEBUG, "Processing cluster event for cluster {} at time {}", cluster_idx, current_time_);
            remove_next_cluster_event();
            handle_cluster_event(current_time_, cluster_idx);
        }
    }

    logger_->log(LogLevel::INFO, "Finished core algorithm loop. Final time: {}, Active clusters: {}", current_time_, num_active_clusters_);

    node_good_.assign(graph_.prizes.size(), false);
    logger_->log(LogLevel::DEBUG, "Marking 'good' nodes reachable from final clusters.");
    if (graph_.root != kInvalidNodeId) {

        ClusterId final_root_cluster = kInvalidClusterId;
        for(ClusterId i = 0; i < static_cast<ClusterId>(clusters_.size()); ++i) {
            if (clusters_[i].contains_root && clusters_[i].merged_into == kInvalidClusterId) {
                final_root_cluster = i;
                break;
            }
        }
        if(final_root_cluster != kInvalidClusterId) {
            logger_->log(LogLevel::DEBUG, "Rooted case: Marking nodes from final root cluster {}.", final_root_cluster);
            mark_nodes_as_good(final_root_cluster);
        } else {
            logger_->log(LogLevel::WARNING, "Rooted case: Could not find the final cluster containing root {}. No nodes marked good.", graph_.root);

            if(graph_.root >= 0 && static_cast<size_t>(graph_.root) < node_good_.size()) {
                node_good_[graph_.root] = true;
            }
        }
    } else {

        logger_->log(LogLevel::DEBUG, "Unrooted case: Marking nodes from {} remaining active clusters.", num_active_clusters_);
        for (ClusterId i = 0; i < static_cast<ClusterId>(clusters_.size()); ++i) {
            if (clusters_[i].active && clusters_[i].merged_into == kInvalidClusterId) {
                logger_->log(LogLevel::TRACE, "Marking nodes starting from active cluster {}.", i);
                mark_nodes_as_good(i);
            }
        }
    }

    return build_core_result();
}

void PCSTCoreAlgorithm::handle_edge_event(double event_time, EdgePartId edge_part_index) {
    logger_->log(LogLevel::TRACE, "Entering handle_edge_event for part {}", edge_part_index);
    assert(static_cast<size_t>(edge_part_index) < edge_parts_.size());

    if (edge_parts_[edge_part_index].deleted) {
        stats_.num_deleted_edge_events++;
        logger_->log(LogLevel::TRACE, "Edge part {} (triggering) already deleted, skipping.", edge_part_index);
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

    logger_->log(LogLevel::TRACE, "Edge event details: Edge={}, Cost={:.4f}", edge_index, current_edge_cost);
    logger_->log(LogLevel::TRACE, "  Part {}: Cluster={}, Sum={:.4f}, FinishedMoat={:.4f}", edge_part_index, cluster_idx_current, sum_current, finished_moat_current);
    logger_->log(LogLevel::TRACE, "  Part {}: Cluster={}, Sum={:.4f}, FinishedMoat={:.4f}", other_edge_part_index, cluster_idx_other, sum_other, finished_moat_other);

    if (cluster_idx_current == cluster_idx_other) {
        stats_.num_merged_edge_events++;
        logger_->log(LogLevel::DEBUG, "Edge part {} connects already merged clusters ({}), ignoring.", edge_part_index, cluster_idx_current);

        edge_parts_[edge_part_index].deleted = true;
        edge_parts_[other_edge_part_index].deleted = true;
        return;
    }

    if (edge_parts_[other_edge_part_index].deleted) {
        stats_.num_deleted_edge_events++;
        logger_->log(LogLevel::TRACE,"Other edge part {} was deleted, skipping event for part {}.", other_edge_part_index, edge_part_index);

        edge_parts_[edge_part_index].deleted = true;
        return;
    }

    double remainder = current_edge_cost - sum_current - sum_other;
    logger_->log(LogLevel::TRACE, "  Remainder: {:.4f}", remainder);

    if (remainder <= eps_ * current_edge_cost || std::fabs(remainder) < eps_) {
        logger_->log(LogLevel::DEBUG, "Edge {} covered (remainder {:.4f}). Merging clusters {} and {}.",
                     edge_index, remainder, cluster_idx_current, cluster_idx_other);
        stats_.total_num_merge_events++;
        phase1_result_edges_.push_back(edge_index);

        edge_parts_[other_edge_part_index].deleted = true;
        logger_->log(LogLevel::TRACE, "Marking other edge part {} as deleted pre-merge.", other_edge_part_index);

        merge_clusters(cluster_idx_current, cluster_idx_other, edge_index, event_time, std::max(0.0, remainder));

    } else {

        Cluster& current_cluster = clusters_[cluster_idx_current];
        Cluster& other_cluster = clusters_[cluster_idx_other];
        EdgePart& current_edge_part_ref = edge_parts_[edge_part_index];
        EdgePart& other_edge_part_ref = edge_parts_[other_edge_part_index];

        if (current_cluster.active && other_cluster.active) {
            logger_->log(LogLevel::DEBUG, "Edge {} growth (Active-Active). Remainder: {:.4f}", edge_index, remainder);
            stats_.total_num_edge_growth_events++;
            stats_.num_active_active_edge_growth_events++;

            assert(remainder > 0.0 && "Remainder should be positive here.");
            double time_to_meet = event_time + remainder / 2.0;
            double val_at_meet_current = sum_current + remainder / 2.0;
            double val_at_meet_other = sum_other + remainder / 2.0;

            logger_->log(LogLevel::TRACE, "  Updating part {}: New event time={:.4f}, New val={:.4f}", edge_part_index, time_to_meet, val_at_meet_current);
            current_edge_part_ref.next_event_val = val_at_meet_current;

            current_edge_part_ref.heap_node = current_cluster.edge_parts.insert(time_to_meet, edge_part_index);

            if (!current_cluster.edge_parts.is_empty()) {
                double min_val;
                EdgePartId min_part;
                [[maybe_unused]] bool success = current_cluster.edge_parts.get_min(&min_val, &min_part);
                assert(success);
                clusters_next_edge_event_.insert(min_val, cluster_idx_current);
            }

            logger_->log(LogLevel::TRACE, "  Updating part {}: Decrease key to time={:.4f}, New val={:.4f}", other_edge_part_index, time_to_meet, val_at_meet_other);
            double old_event_time_other = other_cluster.active_start_time + other_edge_part_ref.next_event_val - finished_moat_other;

            if (other_edge_part_ref.heap_node != nullptr) {
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
                logger_->log(LogLevel::WARNING, "Other edge part {} heap node is null, cannot decrease key.", other_edge_part_index);

                other_edge_part_ref.next_event_val = val_at_meet_other;
            }
        }

        else {
            logger_->log(LogLevel::DEBUG, "Edge {} growth (Active-Inactive). Remainder: {:.4f}", edge_index, remainder);
            assert(current_cluster.active != other_cluster.active);
            assert(remainder > 0.0);
            stats_.total_num_edge_growth_events++;
            stats_.num_active_inactive_edge_growth_events++;

            Cluster& active_cluster = current_cluster.active ? current_cluster : other_cluster;
            Cluster& inactive_cluster = current_cluster.active ? other_cluster : current_cluster;
            EdgePart& active_part = current_cluster.active ? current_edge_part_ref : other_edge_part_ref;
            EdgePart& inactive_part = current_cluster.active ? other_edge_part_ref : current_edge_part_ref;
            EdgePartId active_part_idx = current_cluster.active ? edge_part_index : other_edge_part_index;
            EdgePartId inactive_part_idx = current_cluster.active ? other_edge_part_index : edge_part_index;
            ClusterId active_cluster_idx = current_cluster.active ? cluster_idx_current : cluster_idx_other;

            double finished_moat_inactive = current_cluster.active ? finished_moat_other : finished_moat_current;

            double time_to_meet = event_time + remainder;

            double val_at_meet_active = current_edge_cost - finished_moat_inactive;

            logger_->log(LogLevel::TRACE, "  Updating active part {}: New event time={:.4f}, New val={:.4f}", active_part_idx, time_to_meet, val_at_meet_active);
            active_part.next_event_val = val_at_meet_active;

            active_part.heap_node = active_cluster.edge_parts.insert(time_to_meet, active_part_idx);

            if (!active_cluster.edge_parts.is_empty()) {
                double min_val;
                EdgePartId min_part;
                [[maybe_unused]] bool success = active_cluster.edge_parts.get_min(&min_val, &min_part);
                assert(success);
                clusters_next_edge_event_.insert(min_val, active_cluster_idx);
            }

            double inactive_deactivation_time = inactive_cluster.active_end_time;
            assert(inactive_deactivation_time >= 0.0 && "Inactive cluster must have a valid end time.");
            logger_->log(LogLevel::TRACE, "  Updating inactive part {}: Decrease key to time={:.4f}, New val={:.4f}", inactive_part_idx, inactive_deactivation_time, finished_moat_inactive);

            if (inactive_part.heap_node != nullptr) {

                double old_event_time_inactive = inactive_cluster.active_end_time + inactive_part.next_event_val - finished_moat_inactive;

                inactive_cluster.edge_parts.decrease_key(inactive_part.heap_node,
                                           old_event_time_inactive,
                                           inactive_deactivation_time);
                inactive_part.next_event_val = finished_moat_inactive;
            } else {

                logger_->log(LogLevel::TRACE, "  Inactive part {} has no heap node. Just updating value.", inactive_part_idx);
                inactive_part.next_event_val = finished_moat_inactive;
            }
        }
    }
    logger_->log(LogLevel::TRACE, "Exiting handle_edge_event for part {}", edge_part_index);
}

void PCSTCoreAlgorithm::handle_cluster_event(double event_time, ClusterId cluster_index) {
    logger_->log(LogLevel::TRACE, "Entering handle_cluster_event for cluster {}", cluster_index);
    assert(static_cast<size_t>(cluster_index) < clusters_.size());
    Cluster& cluster = clusters_[cluster_index];

    assert(cluster.active && "Cluster deactivation event for an already inactive cluster!");
    if (!cluster.active) return;

    cluster.active = false;
    cluster.active_end_time = event_time;
    assert(event_time >= cluster.active_start_time && "Deactivation time cannot be before start time.");
    cluster.moat = cluster.active_end_time - cluster.active_start_time;
    num_active_clusters_--;

    logger_->log(LogLevel::DEBUG, "Cluster {} deactivated at time {:.4f}. Moat size: {:.4f}. Active clusters remaining: {}",
                 cluster_index, event_time, cluster.moat, num_active_clusters_);

    if (!cluster.edge_parts.is_empty()) {
        clusters_next_edge_event_.delete_element(cluster_index);
        logger_->log(LogLevel::TRACE, "Removed cluster {} from edge event queue.", cluster_index);
    }

    logger_->log(LogLevel::TRACE, "Exiting handle_cluster_event for cluster {}", cluster_index);
}

ClusterId PCSTCoreAlgorithm::merge_clusters(ClusterId cluster1_idx, ClusterId cluster2_idx, EdgeId merge_edge_idx, double event_time, double remainder) {
    assert(static_cast<size_t>(cluster1_idx) < clusters_.size());
    assert(static_cast<size_t>(cluster2_idx) < clusters_.size());
    assert(cluster1_idx != cluster2_idx);

    clusters_.emplace_back(&pairing_heap_buffer_);
    ClusterId new_cluster_idx = clusters_.size() - 1;
    logger_->log(LogLevel::DEBUG, "Merging clusters {} and {} into new cluster {} along edge {} at time {:.4f}",
                 cluster1_idx, cluster2_idx, new_cluster_idx, merge_edge_idx, event_time);

    Cluster& cluster1 = clusters_[cluster1_idx];
    Cluster& cluster2 = clusters_[cluster2_idx];
    Cluster& new_cluster = clusters_[new_cluster_idx];

    if (cluster1.active && cluster2.active) {
        stats_.num_active_active_merge_events++;
        logger_->log(LogLevel::TRACE, "  Merge type: Active-Active");
    } else {
        assert(cluster1.active != cluster2.active && "Cannot merge two inactive clusters.");
        stats_.num_active_inactive_merge_events++;
        logger_->log(LogLevel::TRACE, "  Merge type: Active-Inactive");

        Cluster& inactive_cluster_ref = cluster1.active ? cluster2 : cluster1;
        ClusterId active_original_cluster_idx = cluster1.active ? cluster1_idx : cluster2_idx;
        ClusterId inactive_original_cluster_idx = cluster1.active ? cluster2_idx : cluster1_idx;

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

            logger_->log(LogLevel::ERROR, "Could not reliably determine active/inactive nodes for merge edge {}. Repr clusters: u={}, v={}. Original clusters: {}, {}",
                         merge_edge_idx, u_repr_cluster, v_repr_cluster, active_original_cluster_idx, inactive_original_cluster_idx);

            active_node = (cluster1.active ? u_node : v_node);
            inactive_node = (cluster1.active ? v_node : u_node);
            assert(false && "Logic error determining active/inactive nodes in merge.");
        }

        inactive_merge_events_.push_back({active_original_cluster_idx, inactive_original_cluster_idx, active_node, inactive_node});
        edge_info_[merge_edge_idx].inactive_merge_event = inactive_merge_events_.size() - 1;
        logger_->log(LogLevel::TRACE, "  Recorded inactive merge event {}: active_cluster={}, inactive_cluster={}, active_node={}, inactive_node={}",
                     edge_info_[merge_edge_idx].inactive_merge_event, active_original_cluster_idx, inactive_original_cluster_idx, active_node, inactive_node);

        if (!inactive_cluster_ref.edge_parts.is_empty()) {
            assert(inactive_cluster_ref.active_end_time >= 0.0);
            double time_diff = (event_time + remainder) - inactive_cluster_ref.active_end_time;
            logger_->log(LogLevel::TRACE, "  Adding offset {} to inactive cluster {} heap.", time_diff, inactive_original_cluster_idx);

            if (time_diff < -eps_) {
                logger_->log(LogLevel::WARNING, "Negative time diff ({}) when updating inactive heap {}. Clamping to 0.", time_diff, inactive_original_cluster_idx);
                time_diff = 0.0;
            }
            inactive_cluster_ref.edge_parts.add_to_heap(std::max(0.0, time_diff));
        }
    }

    if (cluster1.active) {
        cluster1.active = false;
        cluster1.active_end_time = event_time + remainder;
        assert(cluster1.active_end_time >= cluster1.active_start_time);
        cluster1.moat = cluster1.active_end_time - cluster1.active_start_time;
        logger_->log(LogLevel::TRACE, "  Deactivating cluster {} at time {:.4f}. Moat: {:.4f}", cluster1_idx, cluster1.active_end_time, cluster1.moat);
        clusters_deactivation_.delete_element(cluster1_idx);
        if (!cluster1.edge_parts.is_empty()) {
            clusters_next_edge_event_.delete_element(cluster1_idx);
        }
        num_active_clusters_--;
    } else {

        logger_->log(LogLevel::TRACE, "  Cluster {} was already inactive.", cluster1_idx);
    }
    cluster1.merged_into = new_cluster_idx;

    if (cluster2.active) {
        cluster2.active = false;
        cluster2.active_end_time = event_time + remainder;
        assert(cluster2.active_end_time >= cluster2.active_start_time);
        cluster2.moat = cluster2.active_end_time - cluster2.active_start_time;
        logger_->log(LogLevel::TRACE, "  Deactivating cluster {} at time {:.4f}. Moat: {:.4f}", cluster2_idx, cluster2.active_end_time, cluster2.moat);
        clusters_deactivation_.delete_element(cluster2_idx);
        if (!cluster2.edge_parts.is_empty()) {
            clusters_next_edge_event_.delete_element(cluster2_idx);
        }
        num_active_clusters_--;
    } else {
        logger_->log(LogLevel::TRACE, "  Cluster {} was already inactive.", cluster2_idx);
    }
    cluster2.merged_into = new_cluster_idx;

    new_cluster.prize_sum = cluster1.prize_sum + cluster2.prize_sum;
    new_cluster.subcluster_moat_sum = cluster1.subcluster_moat_sum + cluster2.subcluster_moat_sum
                                      + cluster1.moat + cluster2.moat;
    new_cluster.contains_root = cluster1.contains_root || cluster2.contains_root;
    new_cluster.active = !new_cluster.contains_root;
    new_cluster.merged_along = merge_edge_idx;
    new_cluster.child_cluster_1 = cluster1_idx;
    new_cluster.child_cluster_2 = cluster2_idx;

    new_cluster.merged_into = kInvalidClusterId;
    new_cluster.moat = 0.0;
    new_cluster.skip_up = kInvalidClusterId;
    new_cluster.skip_up_sum = 0.0;
    new_cluster.active_end_time = -1.0;

    new_cluster.edge_parts = PairingHeapType::meld(&cluster1.edge_parts, &cluster2.edge_parts);

    if (new_cluster.active) {
        new_cluster.active_start_time = event_time + remainder;
        num_active_clusters_++;

        double potential_deactivation_time = new_cluster.active_start_time
                                             + new_cluster.prize_sum
                                             - new_cluster.subcluster_moat_sum;
        logger_->log(LogLevel::TRACE, "  New cluster {} activated at time {:.4f}. Active clusters: {}", new_cluster_idx, new_cluster.active_start_time, num_active_clusters_);
        logger_->log(LogLevel::TRACE, "    PrizeSum={:.4f}, SubMoatSum={:.4f}, DeactivationTime={:.4f}",
                     new_cluster.prize_sum, new_cluster.subcluster_moat_sum, potential_deactivation_time);

        if (potential_deactivation_time < new_cluster.active_start_time - eps_) {
            logger_->log(LogLevel::WARNING, "  Potential deactivation time ({:.4f}) is before start time ({:.4f}) for new cluster {}. Clamping.",
                         potential_deactivation_time, new_cluster.active_start_time, new_cluster_idx);
            potential_deactivation_time = new_cluster.active_start_time;
        }

        clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);

        if (!new_cluster.edge_parts.is_empty()) {
            double min_val;
            EdgePartId min_part;
            [[maybe_unused]] bool success = new_cluster.edge_parts.get_min(&min_val, &min_part);
            assert(success);
            logger_->log(LogLevel::TRACE, "  New cluster {} added to edge event queue. Min time: {:.4f}", new_cluster_idx, min_val);
            clusters_next_edge_event_.insert(min_val, new_cluster_idx);
        }
    } else {
        logger_->log(LogLevel::TRACE, "  New cluster {} contains root, remains inactive.", new_cluster_idx);
    }

    return new_cluster_idx;
}

std::optional<std::pair<double, std::pair<ClusterId, EdgePartId>>> PCSTCoreAlgorithm::get_next_edge_event() {
    auto min_cluster_event = clusters_next_edge_event_.get_min();
    if (!min_cluster_event) {
        logger_->log(LogLevel::TRACE,"Global edge event queue is empty.");
        return std::nullopt;
    }

    double global_event_time = min_cluster_event->first;
    ClusterId cluster_index = min_cluster_event->second;

    assert(static_cast<size_t>(cluster_index) < clusters_.size());

    if (clusters_[cluster_index].edge_parts.is_empty()) {
        logger_->log(LogLevel::ERROR, "Mismatch: Global edge queue has event for cluster {} but its local heap is empty! Removing stale global entry.", cluster_index);
        clusters_next_edge_event_.delete_element(cluster_index);

        while(true) {
            min_cluster_event = clusters_next_edge_event_.get_min();
            if (!min_cluster_event) {
                logger_->log(LogLevel::TRACE,"Global edge event queue became empty after removing stale entries.");
                return std::nullopt;
            }
            global_event_time = min_cluster_event->first;
            cluster_index = min_cluster_event->second;
            assert(static_cast<size_t>(cluster_index) < clusters_.size());
            if (!clusters_[cluster_index].edge_parts.is_empty()) {

                break;
            }

            logger_->log(LogLevel::ERROR, "Mismatch: Global edge queue has event for cluster {} but its local heap is empty! Removing stale global entry.", cluster_index);
            clusters_next_edge_event_.delete_element(cluster_index);
        }

    }

    double actual_heap_min_val = std::numeric_limits<double>::infinity();
    EdgePartId edge_part_index = kInvalidEdgePartId;

    bool success = clusters_[cluster_index].edge_parts.get_min(&actual_heap_min_val, &edge_part_index);

    if (!success) {

        logger_->log(LogLevel::FATAL, "Internal Error: Failed to get_min from supposedly non-empty heap for cluster {}!", cluster_index);

        throw std::runtime_error(std::format("Internal Error: Failed get_min for cluster {}", cluster_index));
    }

    double event_time = actual_heap_min_val;

    if (std::fabs(global_event_time - actual_heap_min_val) > eps_ * std::fabs(global_event_time)) {
        logger_->log(LogLevel::WARNING, "Mismatch between global edge event time ({:.6f}) and cluster {} heap min time ({:.6f}). Using heap min.",
                     global_event_time, cluster_index, actual_heap_min_val);

    }

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
        logger_->log(LogLevel::TRACE, "Re-inserted cluster {} into edge event queue with new min time {:.4f}", cluster_index, tmp_value);
    } else {
        logger_->log(LogLevel::TRACE, "Cluster {} edge heap is now empty.", cluster_index);
    }
}

std::optional<std::pair<double, ClusterId>> PCSTCoreAlgorithm::get_next_cluster_event() {
    return clusters_deactivation_.get_min();
}

void PCSTCoreAlgorithm::remove_next_cluster_event() {
    [[maybe_unused]] auto deleted_event = clusters_deactivation_.delete_min();
    assert(deleted_event.has_value());
}

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

    *total_sum = 0.0;
    *current_cluster_index = endpoint_node;

    path_compression_visited_.clear();

    while (clusters_[*current_cluster_index].merged_into != kInvalidClusterId) {
        ClusterId cluster_id = *current_cluster_index;
        path_compression_visited_.push_back({cluster_id, *total_sum});

        if (clusters_[cluster_id].skip_up != kInvalidClusterId) {
            *total_sum += clusters_[cluster_id].skip_up_sum;
            *current_cluster_index = clusters_[cluster_id].skip_up;
            logger_->log(LogLevel::TRACE, "Path compression: Skipping from {} to {} (sum {})", cluster_id, *current_cluster_index, clusters_[cluster_id].skip_up_sum);
        } else {

            *total_sum += clusters_[cluster_id].moat;
            *current_cluster_index = clusters_[cluster_id].merged_into;
            logger_->log(LogLevel::TRACE, "Path traversal: Moving from {} to {} (added moat {})", cluster_id, *current_cluster_index, clusters_[cluster_id].moat);
        }
        assert(static_cast<size_t>(*current_cluster_index) < clusters_.size());
    }

    if (!path_compression_visited_.empty()) {
        logger_->log(LogLevel::TRACE, "Applying path compression for {} visited nodes. Final root: {}", path_compression_visited_.size(), *current_cluster_index);
        for (const auto& visited : path_compression_visited_) {
            ClusterId visited_cluster_index = visited.first;
            double visited_sum = visited.second;
            clusters_[visited_cluster_index].skip_up = *current_cluster_index;
            clusters_[visited_cluster_index].skip_up_sum = *total_sum - visited_sum;
            assert(clusters_[visited_cluster_index].skip_up_sum >= -eps_ && "Skip up sum should be non-negative");
        }
    }

    Cluster& root_cluster = clusters_[*current_cluster_index];
    if (root_cluster.active) {
        *finished_moat_sum = *total_sum;
        *total_sum += current_time_ - root_cluster.active_start_time;
        logger_->log(LogLevel::TRACE, "Root cluster {} is active. FinishedMoat={:.4f}, TotalSum={:.4f} (CurrentTime={:.4f}, StartTime={:.4f})",
                     *current_cluster_index, *finished_moat_sum, *total_sum, current_time_, root_cluster.active_start_time);
        assert(*total_sum >= *finished_moat_sum - eps_);
    } else {

        *total_sum += root_cluster.moat;
        *finished_moat_sum = *total_sum;
        logger_->log(LogLevel::TRACE, "Root cluster {} is inactive. Moat={:.4f}, FinishedMoat=TotalSum={:.4f}",
                     *current_cluster_index, root_cluster.moat, *total_sum);
    }
    assert(*total_sum >= -eps_ && *finished_moat_sum >= -eps_);
}

void PCSTCoreAlgorithm::mark_nodes_as_good(ClusterId start_cluster_index) {
    logger_->log(LogLevel::TRACE, "Entering mark_nodes_as_good from cluster {}", start_cluster_index);
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
                    logger_->log(LogLevel::TRACE, "Marked original node {} as good.", current_cluster_idx);
                }
            } else {
                logger_->log(LogLevel::ERROR, "Cluster {} appears to be original node but index is out of range [0, {}).", current_cluster_idx, node_good_.size());
                assert(false);
            }
        } else {

            logger_->log(LogLevel::TRACE, "Exploring children ({}, {}) of merged cluster {}", cluster.child_cluster_1, cluster.child_cluster_2, current_cluster_idx);
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
    logger_->log(LogLevel::TRACE, "Exiting mark_nodes_as_good from cluster {}", start_cluster_index);
}

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
    result.final_cluster_state = std::move(clusters_);
    logger_->log(LogLevel::DEBUG, "Moved {} inactive merge events to result.", result.inactive_merge_events.size());

    return result;
}

}