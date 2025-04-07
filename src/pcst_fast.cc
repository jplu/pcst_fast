#include "pcst_fast.h"
#include "pruning_context.h"
#include "pruner_factory.h"
#include "ipruner.h"
#include "logger.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdarg>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <optional>
#include <memory>
#include <cassert>


namespace cluster_approx {

    using internal::LogLevel;

    PCSTFast::PruningMethod PCSTFast::parse_pruning_method(std::string_view input) {
        std::string lower_input;
        lower_input.reserve(input.length());
        std::transform(input.begin(), input.end(), std::back_inserter(lower_input),
                       [](unsigned char c){ return std::tolower(c); });

        if (lower_input == "none") return PruningMethod::kNoPruning;
        if (lower_input == "simple") return PruningMethod::kSimplePruning;
        if (lower_input == "gw") return PruningMethod::kGWPruning;
        if (lower_input == "strong") return PruningMethod::kStrongPruning;
        if (lower_input == "connectfinal") return PruningMethod::kConnectFinalComponents;

        return PruningMethod::kUnknownPruning;
    }

    PCSTFast::PCSTFast(const std::vector<std::pair<IndexType, IndexType>>& edges,
                       const std::vector<ValueType>& prizes,
                       const std::vector<ValueType>& costs,
                       IndexType root,
                       int target_num_active_clusters,
                       PruningMethod pruning,
                       internal::Logger& logger
                       )
        : edges_(edges),
          prizes_(prizes),
          costs_(costs),
          root_(root),
          target_num_active_clusters_( (root >= 0 && target_num_active_clusters > 0) ? 0 : target_num_active_clusters),
          clusters_deactivation_(),
          clusters_next_edge_event_(),
          logger_(logger)
    {
        logger_.log(LogLevel::INFO, "PCSTFast Constructor Entry.\n");

        const size_t num_nodes = prizes_.size();
        const size_t num_edges = edges_.size();

        logger_.log(LogLevel::INFO, "Input sizes: Nodes=%zu, Edges=%zu\n", num_nodes, num_edges);
        logger_.log(LogLevel::INFO, "Parameters: Root=%d, TargetClusters=%d, Pruning=%d, Verbosity=%d\n",
                    root_, target_num_active_clusters_, static_cast<int>(pruning), static_cast<int>(logger_.get_level()));

        if (edges_.size() != costs_.size()) {
             logger_.log(LogLevel::FATAL, "Error: Mismatch between number of edges (%zu) and costs (%zu).\n", edges_.size(), costs_.size());
             throw std::invalid_argument("Edges and costs vectors must have the same size.");
        }

        if (root_ >= static_cast<IndexType>(num_nodes)) {
            logger_.log(LogLevel::FATAL, "Error: Root node index %d out of range [0, %zu).\n", root_, num_nodes);
            throw std::invalid_argument("Root node index out of range.");
        }

        if (target_num_active_clusters_ < 0) {
            logger_.log(LogLevel::FATAL, "Error: Target number of active clusters %d cannot be negative.\n", target_num_active_clusters_);
            throw std::invalid_argument("Target number of active clusters cannot be negative.");
        }

        if (root_ >= 0 && target_num_active_clusters > 0) {
            logger_.log(LogLevel::WARNING, "Warning: target_num_active_clusters > 0 is ignored in the rooted case (using 0 internally).\n");
        }

        try {
            pruner_ = internal::create_pruner(pruning);
            std::string pruner_name = "Unknown";
             if (pruning == PruningMethod::kNoPruning) pruner_name="No";
             else if (pruning == PruningMethod::kSimplePruning) pruner_name="Simple";
             else if (pruning == PruningMethod::kGWPruning) pruner_name="GW";
             else if (pruning == PruningMethod::kStrongPruning) pruner_name="Strong";
             else if (pruning == PruningMethod::kConnectFinalComponents) pruner_name="ConnectFinal";

            logger_.log(LogLevel::INFO, "Instantiated %s Pruner.\n", pruner_name.c_str());
        } catch (const std::invalid_argument& e) {
            logger_.log(LogLevel::FATAL, "Error creating pruner: %s\n", e.what());
            throw;
        }

        logger_.log(LogLevel::DEBUG, "Resizing internal vectors...\n");
         edge_parts_.resize(2 * num_edges);
         edge_info_.resize(num_edges);
         node_good_.resize(num_nodes, 0);
         clusters_.reserve(2 * num_nodes);
         inactive_merge_events_.reserve(num_edges);
         phase1_result_.reserve(num_edges);

        logger_.log(LogLevel::DEBUG, "Initializing algorithm state...\n");
        current_time_ = 0.0;
        stats_.reset();
        initialize_clusters();
        initialize_edges();

        logger_.log(LogLevel::DEBUG, "Initializing cluster_next_edge_event PQ...\n");
        for (IndexType ii = 0; ii < static_cast<IndexType>(clusters_.size()) && ii < static_cast<IndexType>(prizes_.size()) ; ++ii) {
            if (clusters_[ii].active && !clusters_[ii].edge_parts.is_empty()) {
                ValueType val;
                PayloadType edge_part;
                if (clusters_[ii].edge_parts.get_min(&val, &edge_part)) {
                    logger_.log(LogLevel::DEBUG, "  Cluster %d initial min edge event at time %.9g (part %d).\n", ii, val, edge_part);
                    clusters_next_edge_event_.insert(val, ii);
                } else {
                     logger_.log(LogLevel::FATAL, "Error: Failed to get min edge for active cluster %d during init, though heap reported not empty.\n", ii);
                     assert(false && "Pairing heap get_min failed after is_empty check");
                }
            }
        }

        logger_.log(LogLevel::INFO, "PCSTFast Constructor Exit.\n");
    }

    PCSTFast::~PCSTFast() {
         logger_.log(LogLevel::INFO, "PCSTFast Destructor.\n");
    }

    void PCSTFast::initialize_clusters() {
        logger_.log(LogLevel::INFO, "initialize_clusters Entry.\n");
        const size_t num_nodes = prizes_.size();
        clusters_.clear();

        for (size_t ii = 0; ii < num_nodes; ++ii) {
            IndexType current_node_index = static_cast<IndexType>(ii);
            logger_.log(LogLevel::DEBUG, "Initializing cluster %d (Node %d).\n", current_node_index, current_node_index);

            if (prizes_[ii] < 0.0) {
                 logger_.log(LogLevel::FATAL, "Error: Prize for node %d is negative (%.9g).\n", current_node_index, prizes_[ii]);
                throw std::invalid_argument("Node prize cannot be negative.");
            }

            clusters_.emplace_back();
            Cluster& current_cluster = clusters_.back();
            current_cluster.active = (current_node_index != root_);
            current_cluster.active_start_time = 0.0;
            current_cluster.active_end_time = (current_cluster.active) ? -1.0 : 0.0;
            current_cluster.prize_sum = prizes_[ii];
            current_cluster.contains_root = (current_node_index == root_);
            current_cluster.moat = 0.0;

             logger_.log(LogLevel::DEBUG, "  Node %d: Prize=%.9g, Active=%d, ContainsRoot=%d\n",
                            current_node_index, current_cluster.prize_sum, current_cluster.active, current_cluster.contains_root);

            if (current_cluster.active) {
                if (prizes_[ii] > kEpsilon) {
                    ValueType deactivation_time = prizes_[ii];
                    logger_.log(LogLevel::DEBUG, "  Node %d: Scheduling deactivation at time %.9g.\n", current_node_index, deactivation_time);
                    clusters_deactivation_.insert(deactivation_time, current_node_index);
                } else {
                    logger_.log(LogLevel::DEBUG, "  Node %d: Prize is negligible, immediately inactive.\n", current_node_index);
                    current_cluster.active = false;
                    current_cluster.active_end_time = 0.0;
                }
            }
        }
        logger_.log(LogLevel::INFO, "initialize_clusters Exit.\n");
    }

    void PCSTFast::initialize_edges() {
        logger_.log(LogLevel::INFO, "initialize_edges Entry.\n");
        const size_t num_edges = edges_.size();
        const size_t num_nodes = prizes_.size();

        for (size_t ii = 0; ii < num_edges; ++ii) {
            const auto& edge = edges_[ii];
            IndexType uu = edge.first;
            IndexType vv = edge.second;
            ValueType cost = costs_[ii];
            IndexType edge_index = static_cast<IndexType>(ii);

            logger_.log(LogLevel::DEBUG, "Initializing edge %d: (%d <-> %d), Cost=%.9g\n", edge_index, uu, vv, cost);

            if (uu < 0 || vv < 0 || static_cast<size_t>(uu) >= num_nodes || static_cast<size_t>(vv) >= num_nodes) {
                logger_.log(LogLevel::FATAL, "Error: Edge %d endpoint index (%d or %d) out of range [0, %zu).\n", edge_index, uu, vv, num_nodes);
                throw std::invalid_argument("Edge endpoint index out of range.");
            }
            if (cost < 0.0) {
                logger_.log(LogLevel::FATAL, "Error: Edge %d cost is negative (%.9g).\n", edge_index, cost);
                throw std::invalid_argument("Edge cost cannot be negative.");
            }
            if (uu == vv) {
                logger_.log(LogLevel::WARNING, "Warning: Ignoring self-loop edge %d (%d -> %d).\n", edge_index, uu, vv);
                assert(static_cast<size_t>(2 * ii + 1) < edge_parts_.size());
                edge_parts_[2 * ii].deleted = true;
                edge_parts_[2 * ii + 1].deleted = true;
                continue;
            }

            assert(static_cast<size_t>(2 * ii + 1) < edge_parts_.size() && "Edge parts index out of bounds");
            EdgePart& uu_part = edge_parts_[2 * ii];
            EdgePart& vv_part = edge_parts_[2 * ii + 1];

            assert(static_cast<size_t>(uu) < clusters_.size() && static_cast<size_t>(vv) < clusters_.size() && "Cluster index out of bounds");
            Cluster& uu_cluster = clusters_[uu];
            Cluster& vv_cluster = clusters_[vv];
            uu_part.deleted = false;
            vv_part.deleted = false;
            ValueType event_time_u, event_time_v;

            if (uu_cluster.active && vv_cluster.active) {
                event_time_u = cost * 0.5;
                event_time_v = cost * 0.5;
                logger_.log(LogLevel::DEBUG, "  Edge %d: Active-Active. Event times = %.9g\n", edge_index, event_time_u);
            } else if (uu_cluster.active) {
                event_time_u = cost;
                event_time_v = 0.0;
                logger_.log(LogLevel::DEBUG, "  Edge %d: Active(U)-Inactive(V). Event time U=%.9g, V=%.9g\n", edge_index, event_time_u, event_time_v);
            } else if (vv_cluster.active) {
                event_time_u = 0.0;
                event_time_v = cost;
                logger_.log(LogLevel::DEBUG, "  Edge %d: Inactive(U)-Active(V). Event time U=%.9g, V=%.9g\n", edge_index, event_time_u, event_time_v);
            } else {
                event_time_u = 0.0;
                event_time_v = 0.0;
                logger_.log(LogLevel::DEBUG, "  Edge %d: Inactive-Inactive. Event times = 0.0\n", edge_index);
            }

            PayloadType part_payload_u = static_cast<PayloadType>(2 * ii);
            PayloadType part_payload_v = static_cast<PayloadType>(2 * ii + 1);
            uu_part.next_event_val = event_time_u;

            if (uu_cluster.active || vv_cluster.active) {
                 logger_.log(LogLevel::DEBUG, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", part_payload_u, uu, uu, event_time_u);
                 uu_part.heap_node = uu_cluster.edge_parts.insert(uu_part.next_event_val, part_payload_u);
            } else {
                 logger_.log(LogLevel::DEBUG, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", part_payload_u, uu);
            }

            vv_part.next_event_val = event_time_v;
            if (uu_cluster.active || vv_cluster.active) {
                logger_.log(LogLevel::DEBUG, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", part_payload_v, vv, vv, event_time_v);
                vv_part.heap_node = vv_cluster.edge_parts.insert(vv_part.next_event_val, part_payload_v);
            } else {
                logger_.log(LogLevel::DEBUG, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", part_payload_v, vv);
            }
        }
        logger_.log(LogLevel::INFO, "initialize_edges Exit.\n");
    }

    bool PCSTFast::run(std::vector<IndexType>* result_nodes, std::vector<IndexType>* result_edges) {
        logger_.log(LogLevel::INFO, "PCSTFast::run Entry.\n");

        if (!result_nodes || !result_edges) {
             logger_.log(LogLevel::FATAL, "Error: Result pointers cannot be null in run().\n");
            throw std::invalid_argument("Result pointers cannot be null.");
        }

        result_nodes->clear();
        result_edges->clear();
        phase1_result_.clear();
        stats_.reset();

        int current_num_active_clusters = 0;
        for(const auto& cluster : clusters_) {
            if (cluster.active && cluster.merged_into == kInvalidIndex) {
                current_num_active_clusters++;
            }
        }
        logger_.log(LogLevel::INFO, "Starting GW phase. Initial active representative clusters: %d (Target: %d)\n",
                       current_num_active_clusters, target_num_active_clusters_);

        int event_count = 0;
        while (current_num_active_clusters > target_num_active_clusters_) {
            event_count++;
            logger_.log(LogLevel::INFO, "--- GW Loop Iteration %d (Active Rep Clusters: %d) ---\n", event_count, current_num_active_clusters);

            NextEvent event = get_next_event();
            logger_.log(LogLevel::INFO, "get_next_event returned: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                           static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);

            if (event.type == NextEvent::Type::NONE) {
                logger_.log(LogLevel::WARNING, "No more events possible. Terminating GW loop prematurely.\n");
                break;
            }

            logger_.log(LogLevel::DEBUG, "-----------------------------------------\n");
            if (event.time < current_time_ - kEpsilon) {
                logger_.log(LogLevel::FATAL, "CRITICAL Error: Time decreased! Current=%.9g, Event=%.9g. Algorithm state likely corrupt.\n", current_time_, event.time);
                assert(false && "Time decreased in main loop");
                return false;
            }
            current_time_ = event.time;
            logger_.log(LogLevel::DEBUG, "Advanced time to %.9g\n", current_time_);

            process_event(event);

            int previous_active_count = current_num_active_clusters;
            current_num_active_clusters = 0;
            for(const auto& cluster : clusters_) {
                if (cluster.active && cluster.merged_into == kInvalidIndex) {
                    current_num_active_clusters++;
                }
            }
            logger_.log(LogLevel::INFO, "End of event processing iteration %d. Active representative clusters: %d (was %d)\n", event_count, current_num_active_clusters, previous_active_count);
        }

        logger_.log(LogLevel::INFO, "Finished GW main loop: Final time %.9g, Total edge events processed %lld, Final active representative clusters %d\n",
                       current_time_, stats_.total_num_edge_events, current_num_active_clusters);
        logger_.log(LogLevel::INFO, "--- Starting Pruning and Output Generation ---\n");
        logger_.log(LogLevel::INFO, "Selecting initial 'good' nodes based on final active/root clusters...\n");

        node_good_.assign(prizes_.size(), 0);
        select_initial_active_clusters(node_good_);

        int good_node_count = 0;
        for(uint8_t flag : node_good_) good_node_count += flag;
        logger_.log(LogLevel::INFO, "Found %d 'good' nodes connected to final structure.\n", good_node_count);

        logger_.log(LogLevel::INFO, "Creating Pruning Context...\n");
        internal::PruningContext context(
            edges_, prizes_, costs_, root_, clusters_, edge_info_,
            inactive_merge_events_, node_good_, phase1_result_,
            logger_
        );

        logger_.log(LogLevel::INFO, "Calling selected pruner strategy's prune method...\n");
        assert(pruner_ && "Pruner strategy object is null!");
        try {
            pruner_->prune(context, *result_nodes, *result_edges);
        } catch (const std::exception& e) {
            logger_.log(LogLevel::FATAL, "Error during pruning phase: %s\n", e.what());
            return false;
        } catch (...) {
            logger_.log(LogLevel::FATAL, "Unknown error during pruning phase.\n");
            return false;
        }

        logger_.log(LogLevel::INFO, "Pruning complete. Final result: Nodes=%zu, Edges=%zu\n", result_nodes->size(), result_edges->size());
        logger_.log(LogLevel::INFO, "PCSTFast::run Exit (Success).\n");
        return true;
    }

    PCSTFast::NextEvent PCSTFast::get_next_event() const {
        logger_.log(LogLevel::DEBUG, "get_next_event Entry.\n");
        NextEvent event;
        ValueType next_edge_time = std::numeric_limits<ValueType>::infinity();
        IndexType next_edge_cluster_index = kInvalidIndex;
        PayloadType next_edge_part_index = kInvalidIndex;
        bool edge_event_found = false;
        ValueType next_cluster_time = std::numeric_limits<ValueType>::infinity();
        IndexType next_cluster_deactivation_index = kInvalidIndex;
        bool cluster_event_found = false;
        ValueType temp_pq_edge_time;
        IndexType temp_pq_edge_cluster;
        ValueType temp_cluster_min_time;
        PayloadType temp_cluster_min_part;
        ValueType temp_pq_cluster_time;
        IndexType temp_pq_cluster_index;

        logger_.log(LogLevel::DEBUG, "Checking clusters_next_edge_event PQ...\n");
        if (clusters_next_edge_event_.get_min(&temp_pq_edge_time, &temp_pq_edge_cluster)) {
            logger_.log(LogLevel::DEBUG, "  Global Edge PQ min: Time=%.9g, Cluster=%d\n", temp_pq_edge_time, temp_pq_edge_cluster);
            if (temp_pq_edge_cluster != kInvalidIndex && static_cast<size_t>(temp_pq_edge_cluster) < clusters_.size()) {
                logger_.log(LogLevel::DEBUG, "  Checking cluster %d's internal heap...\n", temp_pq_edge_cluster);
                const Cluster& cluster_ref = clusters_[temp_pq_edge_cluster];
                if (cluster_ref.edge_parts.get_min(&temp_cluster_min_time, &temp_cluster_min_part)) {
                    logger_.log(LogLevel::DEBUG, "  Cluster %d internal min: Time=%.9g, Part=%d\n", temp_pq_edge_cluster, temp_cluster_min_time, temp_cluster_min_part);
                    next_edge_time = temp_cluster_min_time;
                    next_edge_cluster_index = temp_pq_edge_cluster;
                    next_edge_part_index = temp_cluster_min_part;
                    edge_event_found = true;
                } else {
                    logger_.log(LogLevel::FATAL, "Error: Inconsistency in get_next_event. Global edge PQ has cluster %d, but its heap is empty.\n", temp_pq_edge_cluster);
                    assert(false && "Global Edge PQ inconsistent with cluster heap");
                }
            } else {
                logger_.log(LogLevel::FATAL, "Error: Invalid cluster index %d retrieved from global edge PQ in get_next_event.\n", temp_pq_edge_cluster);
            }
        } else {
            logger_.log(LogLevel::DEBUG, "  Global Edge PQ is empty.\n");
        }

        logger_.log(LogLevel::DEBUG, "Checking clusters_deactivation PQ...\n");
        if (clusters_deactivation_.get_min(&temp_pq_cluster_time, &temp_pq_cluster_index)) {
            logger_.log(LogLevel::DEBUG, "  Deactivation PQ min: Time=%.9g, Cluster=%d\n", temp_pq_cluster_time, temp_pq_cluster_index);
            if (temp_pq_cluster_index != kInvalidIndex) {
                next_cluster_time = temp_pq_cluster_time;
                next_cluster_deactivation_index = temp_pq_cluster_index;
                cluster_event_found = true;
            } else {
                logger_.log(LogLevel::FATAL, "Error: Invalid cluster index %d retrieved from deactivation PQ in get_next_event.\n", temp_pq_cluster_index);
            }
        } else {
            logger_.log(LogLevel::DEBUG, "  Deactivation PQ is empty.\n");
        }

        logger_.log(LogLevel::DEBUG, "Determining next event: EdgeFound=%d (Time=%.9g), ClusterFound=%d (Time=%.9g)\n",
                       edge_event_found, next_edge_time, cluster_event_found, next_cluster_time);
        if (edge_event_found && cluster_event_found) {
            if (next_edge_time <= next_cluster_time + kEpsilon * std::max(1.0, std::abs(next_cluster_time))) {
                logger_.log(LogLevel::TRACE,"  Edge event is next (or tie).\n");
                event.type = NextEvent::Type::EDGE;
                event.time = next_edge_time;
                event.cluster_index = next_edge_cluster_index;
                event.edge_part_index = next_edge_part_index;
            } else {
                logger_.log(LogLevel::TRACE,"  Cluster event is next.\n");
                event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
                event.time = next_cluster_time;
                event.cluster_index = next_cluster_deactivation_index;
            }
        } else if (edge_event_found) {
            logger_.log(LogLevel::TRACE,"  Only Edge event found.\n");
            event.type = NextEvent::Type::EDGE;
            event.time = next_edge_time;
            event.cluster_index = next_edge_cluster_index;
            event.edge_part_index = next_edge_part_index;
        } else if (cluster_event_found) {
            logger_.log(LogLevel::TRACE,"  Only Cluster event found.\n");
            event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
            event.time = next_cluster_time;
            event.cluster_index = next_cluster_deactivation_index;
        } else {
            logger_.log(LogLevel::DEBUG,"  No events found.\n");
        }
        logger_.log(LogLevel::DEBUG, "get_next_event Exit. Returning Type=%d, Time=%.9g\n", static_cast<int>(event.type), event.time);
        return event;
    }

    void PCSTFast::process_event(const NextEvent& event) {
        logger_.log(LogLevel::INFO, "process_event Entry: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                       static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);
         if (event.type == NextEvent::Type::EDGE) {
             stats_.total_num_edge_events++;
             logger_.log(LogLevel::DEBUG, "Processing EDGE event for Cluster %d, EdgePart %d\n", event.cluster_index, event.edge_part_index);

             if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size()) {
                 logger_.log(LogLevel::WARNING, "Warning: Edge event cluster index %d invalid. Skipping.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index);
                 return;
             }
             if (!clusters_[event.cluster_index].active) {
                 logger_.log(LogLevel::DEBUG, "Edge event for inactive cluster %d. Removing from global PQ. Skipping logic.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index);
                 return;
             }
             if (event.edge_part_index < 0 || static_cast<size_t>(event.edge_part_index) >= edge_parts_.size()) {
                 logger_.log(LogLevel::FATAL, "Error: Invalid edge part index %d in edge event. Skipping logic.\n", event.edge_part_index);
                 assert(false && "Invalid edge part index in event processing");
                 return;
             }
             logger_.log(LogLevel::DEBUG, "Deleting min edge part from cluster %d's heap...\n", event.cluster_index);

             ValueType deleted_val; PayloadType deleted_part;
             if (!clusters_[event.cluster_index].edge_parts.delete_min(&deleted_val, &deleted_part)) {
                 logger_.log(LogLevel::FATAL, "Error: Failed to delete min edge from cluster %d's heap (event time %.9g).\n", event.cluster_index, event.time);
                 clusters_next_edge_event_.delete_element(event.cluster_index);
                 assert(clusters_[event.cluster_index].edge_parts.is_empty() && "delete_min failed but heap not empty");
                 return;
             }
             logger_.log(LogLevel::TRACE, "  Deleted Part=%d, Value=%.9g. Expected Part=%d, Value=%.9g.\n", deleted_part, deleted_val, event.edge_part_index, event.time);
             if (deleted_part != event.edge_part_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                 logger_.log(LogLevel::FATAL, "Error: Mismatch in edge event deletion for cluster %d! Deleted {%d, %.9g}, Expected {%d, %.9g}\n",
                                event.cluster_index, deleted_part, deleted_val, event.edge_part_index, event.time);
                 assert(false && "Mismatch during edge event deletion");
             }
             logger_.log(LogLevel::DEBUG, "Updating global edge PQ for cluster %d after deletion.\n", event.cluster_index);
             update_cluster_queues_post_growth(event.cluster_index);

             if (edge_parts_[event.edge_part_index].deleted) {
                 stats_.num_deleted_edge_events++;
                 logger_.log(LogLevel::DEBUG, "Edge part %d already deleted (flag=true), skipping event logic.\n", event.edge_part_index);
                 return;
             }
             logger_.log(LogLevel::DEBUG, "Calling handle_edge_event for EdgePart %d.\n", event.edge_part_index);
             handle_edge_event(event.edge_part_index);

          } else if (event.type == NextEvent::Type::CLUSTER_DEACTIVATION) {
             stats_.num_cluster_events++;
             logger_.log(LogLevel::DEBUG, "Processing CLUSTER_DEACTIVATION event for Cluster %d.\n", event.cluster_index);
             logger_.log(LogLevel::DEBUG, "Deleting min from cluster deactivation PQ...\n");

             ValueType deleted_val; IndexType deleted_cluster;
             if (!clusters_deactivation_.delete_min(&deleted_val, &deleted_cluster)) {
                 logger_.log(LogLevel::WARNING, "Warning: Failed to delete cluster deactivation event from PQ (expected cluster %d). Maybe already merged?\n", event.cluster_index);
                 if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size() || !clusters_[event.cluster_index].active) {
                     logger_.log(LogLevel::DEBUG, "Skipping cluster deactivation logic as cluster %d is invalid or already inactive.\n", event.cluster_index);
                     return;
                 }
                 logger_.log(LogLevel::FATAL, "Error: PQ delete failed but cluster %d is still active!\n", event.cluster_index);
                 assert(false && "Deactivation PQ delete failed but cluster still active");
             } else {
                 logger_.log(LogLevel::TRACE, "  Deleted Cluster=%d, Value=%.9g. Expected Cluster=%d, Value=%.9g.\n", deleted_cluster, deleted_val, event.cluster_index, event.time);
                 if (deleted_cluster != event.cluster_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                     logger_.log(LogLevel::WARNING, "Warning: Mismatch in cluster deactivation event deletion (possibly due to prior merge). Proceeding with target cluster %d.\n", event.cluster_index);
                 }
                 if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size() || !clusters_[event.cluster_index].active) {
                     logger_.log(LogLevel::DEBUG, "Skipping cluster deactivation logic as target cluster %d is invalid or already inactive.\n", event.cluster_index);
                     return;
                 }
             }
             logger_.log(LogLevel::DEBUG, "Calling handle_cluster_deactivation_event for Cluster %d.\n", event.cluster_index);
             handle_cluster_deactivation_event(event.cluster_index);
         } else {
             logger_.log(LogLevel::FATAL, "Error: Processed an event of type NONE.\n");
             assert(false && "Processed event of type NONE");
         }
        logger_.log(LogLevel::INFO, "process_event Exit.\n");
    }

    void PCSTFast::handle_edge_event(PayloadType edge_part_index) {
        logger_.log(LogLevel::INFO, "handle_edge_event Entry: EdgePart=%d\n", edge_part_index);
        logger_.log(LogLevel::DEBUG, "Getting edge processing info for part %d...\n", edge_part_index);
        std::optional<EdgeProcessingInfo> opt_info = get_edge_processing_info(edge_part_index);

        if (!opt_info) {
            logger_.log(LogLevel::WARNING, "Warning: Could not get valid edge processing info for part %d. Skipping.\n", edge_part_index);
            return;
        }
        EdgeProcessingInfo info = *opt_info;

        logger_.log(LogLevel::DEBUG, "Edge Event Details: Edge %d (%d <-> %d), Cost %.9g\n",
                       info.edge_idx, edges_[info.edge_idx].first, edges_[info.edge_idx].second, info.cost);
        logger_.log(LogLevel::DEBUG, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                       info.current_part_idx, info.current_cluster_idx, info.sum_current, info.finished_moat_current);
        logger_.log(LogLevel::DEBUG, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                       info.other_part_idx, info.other_cluster_idx, info.sum_other, info.finished_moat_other);
        logger_.log(LogLevel::DEBUG, "  Remainder: %.9g\n", info.remainder);

        assert(info.current_cluster_idx != kInvalidIndex && static_cast<size_t>(info.current_cluster_idx) < clusters_.size() && "Invalid current cluster index in handle_edge_event");
        assert(info.other_cluster_idx != kInvalidIndex && static_cast<size_t>(info.other_cluster_idx) < clusters_.size() && "Invalid other cluster index in handle_edge_event");
        assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size() && "Invalid current part index");
        assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size() && "Invalid other part index");

        if (info.current_cluster_idx == info.other_cluster_idx) {
            stats_.num_merged_edge_events++;
            logger_.log(LogLevel::DEBUG, "Clusters %d already merged. Marking other edge part %d as deleted.\n",
                           info.current_cluster_idx, info.other_part_idx);
            edge_parts_[info.other_part_idx].deleted = true;
            logger_.log(LogLevel::INFO, "handle_edge_event Exit (Clusters already merged).\n");
            return;
        }

        bool merge_condition = (info.remainder <= kEpsilon * std::max(1.0, info.cost) || info.remainder <= kEpsilon );
        logger_.log(LogLevel::DEBUG, "Merge condition check: Remainder=%.9g, Epsilon*Max(1,Cost)=%.9g, Epsilon=%.9g -> Merge=%d\n",
                       info.remainder, kEpsilon * std::max(1.0, info.cost), kEpsilon, merge_condition);

        if (merge_condition) {
            stats_.total_num_merge_events++;
            logger_.log(LogLevel::DEBUG, "Merge condition met for edge %d.\n", info.edge_idx);
            logger_.log(LogLevel::INFO, "Adding edge %d to phase1_result.\n", info.edge_idx);
            phase1_result_.push_back(info.edge_idx);
            logger_.log(LogLevel::INFO, "Marking other edge part %d as deleted.\n", info.other_part_idx);
            edge_parts_[info.other_part_idx].deleted = true;
            logger_.log(LogLevel::DEBUG, "Calling merge_clusters...\n");
            merge_clusters(info);
        } else {
            stats_.total_num_edge_growth_events++;
            logger_.log(LogLevel::DEBUG, "Growth condition for edge %d.\n", info.edge_idx);
            assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size() && "Invalid other cluster index before growth check");
            if (clusters_[info.other_cluster_idx].active) {
                logger_.log(LogLevel::DEBUG, "Calling handle_active_active_growth...\n");
                handle_active_active_growth(info);
            } else {
                logger_.log(LogLevel::DEBUG, "Calling handle_active_inactive_growth...\n");
                handle_active_inactive_growth(info);
            }
        }
        logger_.log(LogLevel::INFO, "handle_edge_event Exit.\n");
    }

    std::optional<PCSTFast::EdgeProcessingInfo> PCSTFast::get_edge_processing_info(PayloadType edge_part_index) {
        logger_.log(LogLevel::DEBUG, "get_edge_processing_info Entry: EdgePart=%d\n", edge_part_index);
        EdgeProcessingInfo info;
        info.current_part_idx = edge_part_index;
        info.other_part_idx = get_other_edge_part_index(edge_part_index);
        info.edge_idx = edge_part_index / 2;

        logger_.log(LogLevel::TRACE, "  Derived OtherPart=%d, EdgeIdx=%d\n", info.other_part_idx, info.edge_idx);
        if (static_cast<size_t>(info.edge_idx) >= costs_.size() ||
            static_cast<size_t>(info.other_part_idx) >= edge_parts_.size())
        {
            logger_.log(LogLevel::FATAL, "Error: Invalid edge index %d or part index %d derived from %d.\n",
                           info.edge_idx, info.other_part_idx, edge_part_index);
            assert(false && "Invalid index derived in get_edge_processing_info");
            return std::nullopt;
        }
        info.cost = costs_[info.edge_idx];

        logger_.log(LogLevel::TRACE, "  EdgeCost=%.9g\n", info.cost);
        logger_.log(LogLevel::DEBUG, "  Getting sum for current part %d...\n", info.current_part_idx);
        get_sum_on_edge_part(info.current_part_idx, &info.sum_current, &info.finished_moat_current, &info.current_cluster_idx);
        logger_.log(LogLevel::DEBUG, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.current_cluster_idx, info.sum_current, info.finished_moat_current);
        logger_.log(LogLevel::DEBUG, "  Getting sum for other part %d...\n", info.other_part_idx);
        get_sum_on_edge_part(info.other_part_idx, &info.sum_other, &info.finished_moat_other, &info.other_cluster_idx);
        logger_.log(LogLevel::DEBUG, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.other_cluster_idx, info.sum_other, info.finished_moat_other);

        if (info.current_cluster_idx == kInvalidIndex || info.other_cluster_idx == kInvalidIndex ||
            static_cast<size_t>(info.current_cluster_idx) >= clusters_.size() ||
            static_cast<size_t>(info.other_cluster_idx) >= clusters_.size())
        {
            logger_.log(LogLevel::WARNING, "Warning: Failed to find valid representatives in get_edge_processing_info.\n");
            return std::nullopt;
        }
        info.remainder = info.cost - info.sum_current - info.sum_other;

        logger_.log(LogLevel::DEBUG, "  Calculated Remainder = %.9g - %.9g - %.9g = %.9g\n", info.cost, info.sum_current, info.sum_other, info.remainder);
        logger_.log(LogLevel::DEBUG, "get_edge_processing_info Exit (Success).\n");
        return info;
    }

    void PCSTFast::merge_clusters(const EdgeProcessingInfo& info) {
        logger_.log(LogLevel::INFO, "merge_clusters Entry: Merging Cluster %d and Cluster %d via Edge %d.\n",
                       info.current_cluster_idx, info.other_cluster_idx, info.edge_idx);
        IndexType cluster1_idx = info.current_cluster_idx;
        IndexType cluster2_idx = info.other_cluster_idx;

        assert(static_cast<size_t>(cluster1_idx) < clusters_.size() && "Invalid cluster1 index in merge");
        assert(static_cast<size_t>(cluster2_idx) < clusters_.size() && "Invalid cluster2 index in merge");

        Cluster& cluster1 = clusters_[cluster1_idx];
        Cluster& cluster2 = clusters_[cluster2_idx];
        IndexType new_cluster_idx = clusters_.size();

        logger_.log(LogLevel::DEBUG, "Creating new cluster at index %d.\n", new_cluster_idx);
        clusters_.emplace_back();
        Cluster& new_cluster = clusters_.back();
        new_cluster.prize_sum = cluster1.prize_sum + cluster2.prize_sum;
        new_cluster.subcluster_moat_sum = cluster1.subcluster_moat_sum + cluster2.subcluster_moat_sum;
        new_cluster.contains_root = cluster1.contains_root || cluster2.contains_root;
        new_cluster.active = !new_cluster.contains_root;
        new_cluster.merged_along = info.edge_idx;
        new_cluster.child_cluster_1 = cluster1_idx;
        new_cluster.child_cluster_2 = cluster2_idx;

        logger_.log(LogLevel::DEBUG, "  New cluster %d: PrizeSum=%.9g, BaseMoatSum=%.9g, ContainsRoot=%d, Active=%d\n",
                       new_cluster_idx, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, new_cluster.contains_root, new_cluster.active);
        logger_.log(LogLevel::DEBUG, "Updating state for child cluster %d (was Active=%d).\n", cluster1_idx, cluster1.active);
        cluster1.merged_into = new_cluster_idx;
        if (cluster1.active) {
            cluster1.active = false;
            cluster1.active_end_time = current_time_;
            cluster1.moat = cluster1.active_end_time - cluster1.active_start_time;
            if (cluster1.moat < 0.0) { logger_.log(LogLevel::WARNING,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster1_idx, cluster1.moat); cluster1.moat = 0.0; }
            logger_.log(LogLevel::DEBUG, "  Child %d became inactive. Moat=%.9g\n", cluster1_idx, cluster1.moat);
            new_cluster.subcluster_moat_sum += cluster1.moat;
        } else {
            logger_.log(LogLevel::DEBUG, "  Child %d already inactive. Adding previous moat %.9g.\n", cluster1_idx, cluster1.moat);
            new_cluster.subcluster_moat_sum += cluster1.moat;
        }
        logger_.log(LogLevel::DEBUG, "Updating state for child cluster %d (was Active=%d).\n", cluster2_idx, cluster2.active);
        cluster2.merged_into = new_cluster_idx;
        if (cluster2.active) {
            stats_.num_active_active_merge_events++;
            cluster2.active = false;
            cluster2.active_end_time = current_time_;
            cluster2.moat = cluster2.active_end_time - cluster2.active_start_time;
            if (cluster2.moat < 0.0) { logger_.log(LogLevel::WARNING,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster2_idx, cluster2.moat); cluster2.moat = 0.0; }
            logger_.log(LogLevel::DEBUG, "  Child %d (Active) became inactive. Moat=%.9g (Active-Active Merge)\n", cluster2_idx, cluster2.moat);
            new_cluster.subcluster_moat_sum += cluster2.moat;
        } else {
            stats_.num_active_inactive_merge_events++;
            logger_.log(LogLevel::DEBUG, "  Child %d already inactive. Adding previous moat %.9g (Active-Inactive Merge)\n", cluster2_idx, cluster2.moat);
            new_cluster.subcluster_moat_sum += cluster2.moat;
            if (!cluster2.contains_root) {
                logger_.log(LogLevel::DEBUG, "  Logging inactive merge event for edge %d (Inactive child %d).\n", info.edge_idx, cluster2_idx);
                log_inactive_merge_event(info);
                ValueType edge_event_update_shift = current_time_ - cluster2.active_end_time;
                if (edge_event_update_shift > kEpsilon) {
                    logger_.log(LogLevel::INFO, "  Shifting edge event times for inactive child cluster %d by %.9g.\n", cluster2_idx, edge_event_update_shift);
                    cluster2.edge_parts.add_to_heap(edge_event_update_shift);
                }
            } else {
                logger_.log(LogLevel::DEBUG, "  Skipping inactive merge event log/shift as inactive child %d contains root.\n", cluster2_idx);
            }
        }
        logger_.log(LogLevel::DEBUG, "  New cluster %d updated SubclusterMoatSum=%.9g\n", new_cluster_idx, new_cluster.subcluster_moat_sum);

        logger_.log(LogLevel::DEBUG, "Melding heaps from cluster %d and %d into %d.\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        new_cluster.edge_parts = PairingHeapType::meld(std::move(cluster1.edge_parts), std::move(cluster2.edge_parts));
        logger_.log(LogLevel::DEBUG, "  Meld complete. New cluster %d heap is%s empty.\n", new_cluster_idx, new_cluster.edge_parts.is_empty() ? "" : " NOT");

        if (new_cluster.active) {
            new_cluster.active_start_time = current_time_;
            ValueType potential_deactivation_time = new_cluster.active_start_time + new_cluster.prize_sum - new_cluster.subcluster_moat_sum;
            logger_.log(LogLevel::DEBUG, "  New cluster %d is active. Potential deactivation time = %.9g + %.9g - %.9g = %.9g\n",
                           new_cluster_idx, new_cluster.active_start_time, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, potential_deactivation_time);
            if (potential_deactivation_time > current_time_ + kEpsilon) {
                logger_.log(LogLevel::INFO, "  Scheduling deactivation for new cluster %d at time %.9g.\n", new_cluster_idx, potential_deactivation_time);
                clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);
            } else {
                logger_.log(LogLevel::INFO, "  New cluster %d immediately inactive (potential time %.9g <= current time %.9g).\n", new_cluster_idx, potential_deactivation_time, current_time_);
                new_cluster.active = false;
                new_cluster.active_end_time = current_time_;
                new_cluster.moat = 0.0;
                new_cluster.subcluster_moat_sum += new_cluster.moat;
            }
        } else {
            new_cluster.active_start_time = current_time_;
            new_cluster.active_end_time = current_time_;
            new_cluster.moat = 0.0;
            new_cluster.subcluster_moat_sum += new_cluster.moat;
            logger_.log(LogLevel::INFO, "  New cluster %d contains root or became immediately inactive, setting Active=false.\n", new_cluster_idx);
        }
        logger_.log(LogLevel::DEBUG, "Updating PQs after merge (removing children %d, %d; adding new %d if needed).\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        update_cluster_queues_post_merge(new_cluster_idx, cluster1_idx, cluster2_idx);
        logger_.log(LogLevel::INFO, "merge_clusters Exit.\n");
    }

     void PCSTFast::log_inactive_merge_event(const EdgeProcessingInfo& info) {
        logger_.log(LogLevel::DEBUG, "log_inactive_merge_event Entry for Edge %d.\n", info.edge_idx);
        inactive_merge_events_.emplace_back();
        InactiveMergeEvent& merge_event = inactive_merge_events_.back();
        merge_event.active_cluster_index = info.current_cluster_idx;
        merge_event.inactive_cluster_index = info.other_cluster_idx;

        assert(static_cast<size_t>(info.edge_idx) < edges_.size() && "Invalid edge index in log_inactive");
        const auto& edge = edges_[info.edge_idx];
        IndexType node_on_current_side = (info.current_part_idx % 2 == 0) ? edge.first : edge.second;
        IndexType node_on_other_side = (info.other_part_idx % 2 == 0) ? edge.first : edge.second;
        merge_event.active_cluster_node = node_on_current_side;
        merge_event.inactive_cluster_node = node_on_other_side;

        if (static_cast<size_t>(info.edge_idx) < edge_info_.size()) {
            edge_info_[info.edge_idx].inactive_merge_event = std::ssize(inactive_merge_events_) - 1;
            logger_.log(LogLevel::INFO, "Logged InactiveMergeEvent %zd for edge %d (Active: C%d/N%d, Inactive: C%d/N%d).\n",
                           inactive_merge_events_.size() - 1, info.edge_idx,
                           merge_event.active_cluster_index, merge_event.active_cluster_node,
                           merge_event.inactive_cluster_index, merge_event.inactive_cluster_node);
        } else {
                logger_.log(LogLevel::FATAL, "Error: Edge index %d out of bounds for edge_info_ in log_inactive_merge_event.\n", info.edge_idx);
                assert(false && "Edge index out of bounds for edge_info_");
        }
        logger_.log(LogLevel::DEBUG, "log_inactive_merge_event Exit.\n");
    }

    void PCSTFast::handle_active_active_growth(const EdgeProcessingInfo& info) {
        stats_.num_active_active_edge_growth_events++;
        ValueType next_event_time = current_time_ + info.remainder * 0.5;
        logger_.log(LogLevel::INFO, "handle_active_active_growth Entry: Edge=%d, Clusters=(%d, %d), Remainder=%.9g, NextEventTime=%.9g\n",
                       info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);

        assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.current_cluster_idx) < clusters_.size());
        assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size());

        EdgePart& current_part = edge_parts_[info.current_part_idx];
        EdgePart& other_part = edge_parts_[info.other_part_idx];
        Cluster& current_cluster = clusters_[info.current_cluster_idx];
        Cluster& other_cluster = clusters_[info.other_cluster_idx];

        logger_.log(LogLevel::DEBUG, "Updating current part %d (cluster %d): Inserting into heap with time %.9g.\n",
                       info.current_part_idx, info.current_cluster_idx, next_event_time);
        current_part.next_event_val = next_event_time;
        current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
        update_cluster_queues_post_growth(info.current_cluster_idx);

        logger_.log(LogLevel::DEBUG, "Updating other part %d (cluster %d): Decreasing key in heap to time %.9g.\n",
                       info.other_part_idx, info.other_cluster_idx, next_event_time);
        other_part.next_event_val = next_event_time;
        logger_.log(LogLevel::DEBUG, "  Deleting cluster %d from global edge PQ before decrease_key.\n", info.other_cluster_idx);
        clusters_next_edge_event_.delete_element(info.other_cluster_idx);
        other_cluster.edge_parts.decrease_key(other_part.heap_node, next_event_time);
        update_cluster_queues_post_growth(info.other_cluster_idx);

        logger_.log(LogLevel::INFO, "handle_active_active_growth Exit.\n");
    }

    void PCSTFast::handle_active_inactive_growth(const EdgeProcessingInfo& info) {
        stats_.num_active_inactive_edge_growth_events++;
        ValueType next_event_time = current_time_ + info.remainder;
        logger_.log(LogLevel::INFO, "handle_active_inactive_growth Entry: Edge=%d, Clusters=(Active:%d, Inactive:%d), Remainder=%.9g, NextEventTime=%.9g\n",
                       info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);

        assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.current_cluster_idx) < clusters_.size());
        assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size());
        assert(!clusters_[info.other_cluster_idx].active && "Other cluster should be inactive");

        EdgePart& current_part = edge_parts_[info.current_part_idx];
        EdgePart& other_part = edge_parts_[info.other_part_idx];
        Cluster& current_cluster = clusters_[info.current_cluster_idx];
        Cluster& other_cluster = clusters_[info.other_cluster_idx];

        logger_.log(LogLevel::DEBUG, "Updating active part %d (cluster %d): Inserting into heap with time %.9g.\n",
                       info.current_part_idx, info.current_cluster_idx, next_event_time);
        current_part.next_event_val = next_event_time;
        current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
        update_cluster_queues_post_growth(info.current_cluster_idx);

        ValueType inactive_reset_time = other_cluster.active_end_time;
        logger_.log(LogLevel::DEBUG, "Updating inactive part %d (cluster %d): Decreasing key in heap to time %.9g (cluster active_end_time).\n",
                       info.other_part_idx, info.other_cluster_idx, inactive_reset_time);
        other_part.next_event_val = inactive_reset_time;
        other_cluster.edge_parts.decrease_key(other_part.heap_node, inactive_reset_time);

        logger_.log(LogLevel::INFO, "handle_active_inactive_growth Exit.\n");
    }

    void PCSTFast::handle_cluster_deactivation_event(IndexType cluster_index) {
        logger_.log(LogLevel::INFO, "handle_cluster_deactivation_event Entry: Cluster=%d\n", cluster_index);
        if (cluster_index < 0 || static_cast<size_t>(cluster_index) >= clusters_.size()) {
            logger_.log(LogLevel::FATAL,"Error: Invalid cluster index %d in deactivation event.\n", cluster_index);
            assert(false && "Invalid cluster index in deactivation event");
            return;
        }
        Cluster& cluster = clusters_[cluster_index];

        if (!cluster.active) {
            logger_.log(LogLevel::DEBUG, "Cluster %d already inactive, skipping deactivation event logic.\n", cluster_index);
            return;
        }
        cluster.active = false;
        cluster.active_end_time = current_time_;
        cluster.moat = cluster.active_end_time - cluster.active_start_time;
        if (cluster.moat < 0.0) {
            logger_.log(LogLevel::WARNING,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster_index, cluster.moat);
            cluster.moat = 0.0;
        }
        logger_.log(LogLevel::INFO, "Cluster Deactivation: Cluster %d at time %.9g (Moat size %.9g).\n",
                       cluster_index, current_time_, cluster.moat);

        if (!cluster.edge_parts.is_empty()) {
            logger_.log(LogLevel::DEBUG, "Removing cluster %d from next edge event PQ as it became inactive.\n", cluster_index);
            clusters_next_edge_event_.delete_element(cluster_index);
        } else {
            logger_.log(LogLevel::DEBUG, "Cluster %d had no edges, no need to remove from edge PQ.\n", cluster_index);
        }
        logger_.log(LogLevel::INFO, "handle_cluster_deactivation_event Exit.\n");
    }

    void PCSTFast::update_cluster_queues_post_growth(IndexType cluster_idx) {
        logger_.log(LogLevel::DEBUG, "update_cluster_queues_post_growth Entry: Cluster=%d\n", cluster_idx);
        if (cluster_idx == kInvalidIndex || static_cast<size_t>(cluster_idx) >= clusters_.size() || !clusters_[cluster_idx].active) {
            logger_.log(LogLevel::DEBUG, "  Skipping PQ update (cluster invalid or inactive).\n");
            return;
        }
        Cluster& cluster = clusters_[cluster_idx];
        logger_.log(LogLevel::DEBUG, "  Deleting cluster %d from global edge PQ (if present).\n", cluster_idx);
        clusters_next_edge_event_.delete_element(cluster_idx);

        if (!cluster.edge_parts.is_empty()) {
            ValueType min_val; PayloadType min_part;
            if (cluster.edge_parts.get_min(&min_val, &min_part)) {
                logger_.log(LogLevel::DEBUG, "  Cluster %d heap not empty. Inserting into edge PQ with Time=%.9g, Part=%d.\n", cluster_idx, min_val, min_part);
                clusters_next_edge_event_.insert(min_val, cluster_idx);
            } else {
                logger_.log(LogLevel::FATAL, "Error: Cluster %d heap not empty but get_min failed in update_queues.\n", cluster_idx);
                assert(false && "get_min failed on non-empty heap");
            }
        } else {
            logger_.log(LogLevel::DEBUG, "  Cluster %d heap is empty, not inserting into edge PQ.\n", cluster_idx);
        }
        logger_.log(LogLevel::DEBUG, "update_cluster_queues_post_growth Exit.\n");
    }

    void PCSTFast::update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index) {
        logger_.log(LogLevel::DEBUG, "update_cluster_queues_post_merge Entry: New=%d, Children=(%d, %d)\n", new_cluster_index, child1_index, child2_index);
        logger_.log(LogLevel::DEBUG, "  Deleting child %d from PQs (if present).\n", child1_index);
        clusters_deactivation_.delete_element(child1_index);
        clusters_next_edge_event_.delete_element(child1_index);
        logger_.log(LogLevel::DEBUG, "  Deleting child %d from PQs (if present).\n", child2_index);
        clusters_deactivation_.delete_element(child2_index);
        clusters_next_edge_event_.delete_element(child2_index);

        if (new_cluster_index != kInvalidIndex && static_cast<size_t>(new_cluster_index) < clusters_.size()) {
            Cluster& new_cluster = clusters_[new_cluster_index];
            if (new_cluster.active && !new_cluster.edge_parts.is_empty()) {
                ValueType min_val; PayloadType min_part;
                if (new_cluster.edge_parts.get_min(&min_val, &min_part)) {
                    logger_.log(LogLevel::DEBUG, "  New cluster %d is active with edges. Inserting into edge PQ with Time=%.9g, Part=%d.\n", new_cluster_index, min_val, min_part);
                    clusters_next_edge_event_.insert(min_val, new_cluster_index);
                } else {
                    logger_.log(LogLevel::FATAL, "Error: New cluster %d heap not empty but get_min failed in update_queues_post_merge.\n", new_cluster_index);
                    assert(false && "get_min failed on non-empty heap post-merge");
                }
            } else {
                logger_.log(LogLevel::DEBUG, "  New cluster %d is inactive or has no edges, not inserting into edge PQ.\n", new_cluster_index);
            }
        } else {
            logger_.log(LogLevel::DEBUG, "  New cluster index %d is invalid, skipping edge PQ update.\n", new_cluster_index);
        }
        logger_.log(LogLevel::DEBUG, "update_cluster_queues_post_merge Exit.\n");
    }

     void PCSTFast::get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                                       ValueType* finished_moat_sum, IndexType* current_cluster_index) {
        logger_.log(LogLevel::DEBUG, "get_sum_on_edge_part Entry: EdgePart=%d\n", edge_part_index);
        IndexType edge_idx = edge_part_index / 2;
        assert(static_cast<size_t>(edge_idx) < edges_.size() && "Invalid edge index in get_sum");

        const auto& edge = edges_[edge_idx];
        IndexType start_node = (edge_part_index % 2 == 0) ? edge.first : edge.second;
        logger_.log(LogLevel::TRACE, "  Start node derived from edge part: %d\n", start_node);

        *total_sum = 0.0;
        *finished_moat_sum = 0.0;
        *current_cluster_index = kInvalidIndex;
         ValueType path_sum = 0.0;

        logger_.log(LogLevel::DEBUG, "  Calling find_representative_and_compress for start node %d...\n", start_node);
        IndexType representative_cluster = find_representative_and_compress(start_node, path_sum);
        logger_.log(LogLevel::DEBUG, "  find_representative result: RepCluster=%d, PathSum=%.9g\n", representative_cluster, path_sum);

        if (representative_cluster == kInvalidIndex || static_cast<size_t>(representative_cluster) >= clusters_.size()) {
            logger_.log(LogLevel::FATAL, "Error: Path compression failed to find representative for node %d (edge part %d).\n", start_node, edge_part_index);
            logger_.log(LogLevel::DEBUG, "get_sum_on_edge_part Exit (Failure).\n");
            assert(false && "Path compression failed");
            return;
        }
        *current_cluster_index = representative_cluster;
        const Cluster& rep_cluster = clusters_[representative_cluster];
        logger_.log(LogLevel::DEBUG, "  Representative cluster %d is Active=%d (StartTime=%.9g, EndTime=%.9g, Moat=%.9g)\n",
                       representative_cluster, rep_cluster.active, rep_cluster.active_start_time, rep_cluster.active_end_time, rep_cluster.moat);

        if (rep_cluster.active) {
            *finished_moat_sum = path_sum;
            *total_sum = path_sum + (current_time_ - rep_cluster.active_start_time);
            logger_.log(LogLevel::TRACE, "  Rep is active. FinishedMoatSum=%.9g, TotalSum = %.9g + (%.9g - %.9g) = %.9g\n",
                           *finished_moat_sum, path_sum, current_time_, rep_cluster.active_start_time, *total_sum);
        } else {
            *total_sum = path_sum + rep_cluster.moat;
            *finished_moat_sum = *total_sum;
            logger_.log(LogLevel::TRACE, "  Rep is inactive. FinishedMoatSum=TotalSum = %.9g + %.9g = %.9g\n",
                           path_sum, rep_cluster.moat, *total_sum);
        }
        logger_.log(LogLevel::DEBUG, "get_sum_on_edge_part Exit (Success). Cluster=%d, TotalSum=%.9g, FinishedMoat=%.9g\n",
                       *current_cluster_index, *total_sum, *finished_moat_sum);
    }

    PCSTFast::IndexType PCSTFast::find_representative_and_compress(IndexType start_node, ValueType& path_sum_out) {
        logger_.log(LogLevel::DEBUG, "find_representative_and_compress Entry: StartNode=%d\n", start_node);
        std::vector<std::pair<IndexType, ValueType>> path_compression_visited_local;
        path_compression_visited_local.reserve(clusters_.size());
         path_sum_out = 0.0;
         IndexType current_node = start_node;
         int steps = 0;

        while (static_cast<size_t>(current_node) < clusters_.size() && clusters_[current_node].merged_into != kInvalidIndex) {
            steps++;
            logger_.log(LogLevel::TRACE, "  Step %d: CurrentNode=%d, PathSum=%.9g. Storing locally.\n", steps, current_node, path_sum_out);
            path_compression_visited_local.emplace_back(current_node, path_sum_out);

            ValueType step_sum;
            IndexType next_node;
            bool used_skip = false;

            assert(static_cast<size_t>(current_node) < clusters_.size() && "Current node index out of bounds in find_rep");
            const Cluster& curr_cluster_ref = clusters_[current_node];

            if (curr_cluster_ref.skip_up >= 0) {
                step_sum = curr_cluster_ref.skip_up_sum;
                next_node = curr_cluster_ref.skip_up;
                used_skip = true;
            } else {
                step_sum = curr_cluster_ref.moat;
                next_node = curr_cluster_ref.merged_into;
            }
            logger_.log(LogLevel::TRACE, "  Step %d: NextNode=%d, StepSum=%.9g (UsedSkip=%d)\n", steps, next_node, step_sum, used_skip);
            path_sum_out += step_sum;
            current_node = next_node;
        }
        logger_.log(LogLevel::DEBUG, "  Traversal finished after %d steps. FinalNode=%d, FinalPathSum=%.9g\n", steps, current_node, path_sum_out);

        if (static_cast<size_t>(current_node) >= clusters_.size()) {
            logger_.log(LogLevel::FATAL, "Error: Path compression traversed out of bounds from start node %d.\n", start_node);
            assert(false && "Path compression out of bounds");
            return kInvalidIndex;
        }
        IndexType representative_cluster = current_node;

        logger_.log(LogLevel::DEBUG, "  Representative cluster found: %d\n", representative_cluster);
        logger_.log(LogLevel::DEBUG, "  Applying path compression for %zu visited nodes (updating member clusters_)...\n", path_compression_visited_local.size());
        for (const auto& [visited_node_idx, sum_at_visited] : path_compression_visited_local) {
            assert(static_cast<size_t>(visited_node_idx) < clusters_.size() && "Invalid visited node index during compression");
            ValueType new_skip_sum = path_sum_out - sum_at_visited;
            logger_.log(LogLevel::TRACE, "    Updating skip pointer for Node %d: SkipUp=%d, SkipSum=%.9g\n",
                           visited_node_idx, representative_cluster, new_skip_sum);
            clusters_[visited_node_idx].skip_up = representative_cluster;
            clusters_[visited_node_idx].skip_up_sum = new_skip_sum;
        }
        logger_.log(LogLevel::DEBUG, "find_representative_and_compress Exit. Rep=%d, PathSum=%.9g\n", representative_cluster, path_sum_out);
        return representative_cluster;
    }

    void PCSTFast::select_initial_active_clusters(std::vector<uint8_t>& node_good_flag) {
        logger_.log(LogLevel::INFO, "select_initial_active_clusters Entry.\n");
        std::vector<IndexType> cluster_queue_local;
        cluster_queue_local.reserve(clusters_.size());

        if (root_ >= 0) {
            logger_.log(LogLevel::INFO, "  Rooted case (Root=%d).\n", root_);
            IndexType root_cluster_final_index = root_;
            logger_.log(LogLevel::DEBUG, "  Finding final representative for root %d...\n", root_);

            int path_len = 0;
            while(static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].merged_into != kInvalidIndex) {
                assert(path_len < static_cast<int>(clusters_.size()) && "Infinite loop detected in merge path?");
                root_cluster_final_index = clusters_[root_cluster_final_index].merged_into;
                path_len++;
            }
            logger_.log(LogLevel::DEBUG, "  Final representative for root %d is cluster %d (path length %d).\n", root_, root_cluster_final_index, path_len);

            if (static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].contains_root) {
                logger_.log(LogLevel::INFO, "  Rooted case: Starting traversal from final root cluster %d.\n", root_cluster_final_index);
                cluster_queue_local.push_back(root_cluster_final_index);
            } else {
                logger_.log(LogLevel::WARNING, "Warning: Could not find final cluster containing root %d. No nodes selected initially.\n", root_);
            }
        } else {
            logger_.log(LogLevel::INFO, "  Unrooted case.\n");
            logger_.log(LogLevel::INFO, "  Unrooted case: Starting traversal from all final active clusters.\n");
            for (size_t ii = 0; ii < clusters_.size(); ++ii) {
                if (clusters_[ii].merged_into == kInvalidIndex && clusters_[ii].active) {
                    logger_.log(LogLevel::INFO, "  Adding final active cluster %zu to initial queue.\n", ii);
                    cluster_queue_local.push_back(static_cast<IndexType>(ii));
                }
            }
        }
        logger_.log(LogLevel::INFO, "  Starting traversal from %zu initial clusters to mark good nodes...\n", cluster_queue_local.size());

        size_t queue_index = 0;
        int nodes_marked_good = 0;
        while (queue_index < cluster_queue_local.size()) {
            IndexType current_cluster_index = cluster_queue_local[queue_index++];
            logger_.log(LogLevel::DEBUG, "    Processing cluster %d from queue (index %zu).\n", current_cluster_index, queue_index - 1);

            if (static_cast<size_t>(current_cluster_index) >= clusters_.size()) {
                logger_.log(LogLevel::WARNING, "Warning: Invalid cluster index %d encountered during good node traversal.\n", current_cluster_index);
                assert(false && "Invalid cluster index in good node traversal");
                continue;
            }
            const Cluster& current_cluster = clusters_[current_cluster_index];

            if (current_cluster.merged_along >= 0) {
                logger_.log(LogLevel::DEBUG, "    Cluster %d is merged. Adding children %d, %d to queue.\n",
                               current_cluster_index, current_cluster.child_cluster_1, current_cluster.child_cluster_2);
                if(current_cluster.child_cluster_1 != kInvalidIndex) cluster_queue_local.push_back(current_cluster.child_cluster_1);
                if(current_cluster.child_cluster_2 != kInvalidIndex) cluster_queue_local.push_back(current_cluster.child_cluster_2);
              } else {
                logger_.log(LogLevel::DEBUG,"    Cluster %d is original node cluster.\n", current_cluster_index);
                if (static_cast<size_t>(current_cluster_index) < node_good_flag.size()) {
                    if(node_good_flag[current_cluster_index] == 0) {
                        logger_.log(LogLevel::INFO, "    Marking node %d as good.\n", current_cluster_index);
                        node_good_flag[current_cluster_index] = 1;
                        nodes_marked_good++;
                    } else {
                        logger_.log(LogLevel::DEBUG, "    Node %d already marked good.\n", current_cluster_index);
                    }
                } else {
                    logger_.log(LogLevel::FATAL, "Error: Original cluster index %d out of bounds for node_good_flag (%zu).\n", current_cluster_index, node_good_flag.size());
                    assert(false && "node_good_flag index out of bounds");
                }
            }
        }
        logger_.log(LogLevel::INFO, "select_initial_active_clusters Exit. Marked %d nodes as good.\n", nodes_marked_good);
    }

    void PCSTFast::get_statistics(Statistics* s) const {
        logger_.log(LogLevel::INFO, "get_statistics Entry.\n");
        if(s) {
            *s = stats_;
            logger_.log(LogLevel::INFO, "  Statistics copied.\n");
        } else {
            logger_.log(LogLevel::WARNING, "Warning: Null pointer passed to get_statistics.\n");
        }
        logger_.log(LogLevel::INFO, "get_statistics Exit.\n");
    }
}
