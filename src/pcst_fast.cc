#include "pcst_fast.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio> // For snprintf
#include <cstdarg> // For va_list, etc.
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <numeric> // For std::iota potentially
#include <optional>

namespace cluster_approx {

    // --- Utility Functions ---

    PCSTFast::PruningMethod PCSTFast::parse_pruning_method(std::string_view input) {
        std::string lower_input;
        lower_input.reserve(input.length());
        // Convert to lower case for case-insensitive comparison
        std::transform(input.begin(), input.end(), std::back_inserter(lower_input),
                       [](unsigned char c){ return std::tolower(c); });

        if (lower_input == "none") return PruningMethod::kNoPruning;
        if (lower_input == "simple") return PruningMethod::kSimplePruning;
        if (lower_input == "gw") return PruningMethod::kGWPruning;
        if (lower_input == "strong") return PruningMethod::kStrongPruning;

        return PruningMethod::kUnknownPruning; // Indicate failure
    }

    // Safe logging helper
    void PCSTFast::log_message(int level, const char* format, ...) const {
        if (output_function_ && verbosity_level_ >= level) {
            char buffer[kOutputBufferSize]; // Use local buffer
            va_list args;
            va_start(args, format);
            // vsnprintf returns the number of characters that *would have* been written
            int needed = vsnprintf(buffer, kOutputBufferSize, format, args);
            va_end(args);

             // Handle potential truncation (though kOutputBufferSize is large)
             if (needed < 0) {
                 // Encoding error
                 output_function_("Logging error: vsnprintf failed.\n");
             } else if (static_cast<size_t>(needed) >= kOutputBufferSize) {
                 // Output was truncated
                  buffer[kOutputBufferSize - 1] = '\0'; // Ensure null termination
                  snprintf(buffer + kOutputBufferSize - 5, 5, "..."); // Add ellipsis safely
                  output_function_(buffer);
                  output_function_("(TRUNCATED)\n");
             } else {
                output_function_(buffer);
             }
        }
    }


    // --- Constructor and Initialization ---

    PCSTFast::PCSTFast(const std::vector<std::pair<int, int> >& edges,
                       const std::vector<ValueType>& prizes,
                       const std::vector<ValueType>& costs,
                       IndexType root,
                       int target_num_active_clusters,
                       PruningMethod pruning,
                       int verbosity_level,
                       void (*output_function)(const char*))
        : edges_(edges), prizes_(prizes), costs_(costs), root_(root),
          target_num_active_clusters_( (root >= 0 && target_num_active_clusters > 0) ? 0 : target_num_active_clusters),
          pruning_(pruning), verbosity_level_(verbosity_level),
          output_function_(output_function),
          clusters_deactivation_(),
          clusters_next_edge_event_()
    {
        log_message(3, "PCSTFast Constructor Entry.\n");
        const size_t num_nodes = prizes_.size();
        const size_t num_edges = edges_.size();
        log_message(3, "Input sizes: Nodes=%zu, Edges=%zu\n", num_nodes, num_edges);
        log_message(3, "Parameters: Root=%d, TargetClusters=%d, Pruning=%d, Verbosity=%d\n",
                    root_, target_num_active_clusters_, static_cast<int>(pruning_), verbosity_level_);


        // --- Input Validation ---
        if (root_ >= static_cast<int>(num_nodes)) {
            log_message(0, "Error: Root node index %d out of range (Num nodes: %zu).\n", root_, num_nodes);
            throw std::invalid_argument("Root node index out of range.");
        }
        if (target_num_active_clusters_ < 0) {
             log_message(0, "Error: Target number of active clusters %d cannot be negative.\n", target_num_active_clusters_);
            throw std::invalid_argument("Target number of active clusters cannot be negative.");
        }
        if (root_ >= 0 && target_num_active_clusters > 0) { // Check original parameter here
            log_message(1, "Warning: target_num_active_clusters > 0 is ignored in the rooted case (using 0 internally).\n");
        }
        if (pruning_ == PruningMethod::kUnknownPruning) {
             log_message(0, "Error: Unknown pruning method specified.\n");
             throw std::invalid_argument("Unknown pruning method specified.");
        }

        // --- Preallocate / Resize Member Vectors ---
        log_message(4, "Resizing internal vectors...\n");
        edge_parts_.resize(2 * num_edges);
        edge_info_.resize(num_edges);
        node_deleted_.resize(num_nodes, 0);
        node_good_.resize(num_nodes, 0);
        build_phase1_included_nodes_.resize(num_nodes, 0);

        clusters_.reserve(2 * num_nodes);
        inactive_merge_events_.reserve(num_edges);
        phase1_result_.reserve(num_edges);
        phase2_result_.reserve(num_edges);
        phase3_result_.reserve(num_edges);
        path_compression_visited_.reserve(num_nodes);
        cluster_queue_.reserve(num_nodes);

        if (pruning_ == PruningMethod::kStrongPruning || pruning_ == PruningMethod::kGWPruning) {
            log_message(4, "Resizing phase3_neighbors_ for GW/Strong pruning.\n");
            phase3_neighbors_.resize(num_nodes);
        }

        if (pruning_ == PruningMethod::kStrongPruning) {
            log_message(4, "Resizing vectors for Strong pruning.\n");
            final_component_label_.resize(num_nodes);
            strong_pruning_parent_.resize(num_nodes);
            strong_pruning_payoff_.resize(num_nodes);
            final_components_.reserve(num_nodes);
            strong_pruning_stack_.reserve(num_nodes);
            strong_pruning_stack2_.reserve(num_nodes);
        }

        // --- Initialize State ---
        log_message(4, "Initializing algorithm state...\n");
        current_time_ = 0.0;
        stats_.reset();

        initialize_clusters();
        initialize_edges();

        log_message(4, "Initializing cluster_next_edge_event PQ...\n");
        for (int ii = 0; ii < std::ssize(clusters_) && ii < std::ssize(prizes_) ; ++ii) {
            if (clusters_[ii].active && !clusters_[ii].edge_parts.is_empty()) {
                ValueType val;
                PayloadType edge_part;
                if (clusters_[ii].edge_parts.get_min(&val, &edge_part)) {
                    log_message(4, "  Cluster %d initial min edge event at time %.9g (part %d).\n", ii, val, edge_part);
                    clusters_next_edge_event_.insert(val, ii);
                } else {
                     log_message(0, "Error: Failed to get min edge for active cluster %d during init.\n", ii);
                }
            }
        }
        log_message(3, "PCSTFast Constructor Exit.\n");
    }

    void PCSTFast::initialize_clusters() {
        log_message(3, "initialize_clusters Entry.\n");
        const size_t num_nodes = prizes_.size();
        clusters_.clear(); // Ensure clean state

        for (size_t ii = 0; ii < num_nodes; ++ii) {
            log_message(4, "Initializing cluster %zu (Node %zu).\n", ii, ii);
            if (prizes_[ii] < 0.0) {
                 log_message(0, "Error: Prize for node %zu is negative (%.9g).\n", ii, prizes_[ii]);
                throw std::invalid_argument("Node prize cannot be negative.");
            }

            clusters_.emplace_back();
            Cluster& current_cluster = clusters_.back();

            current_cluster.active = (static_cast<int>(ii) != root_);
            current_cluster.active_start_time = 0.0;
            current_cluster.active_end_time = (current_cluster.active) ? -1.0 : 0.0;
            current_cluster.prize_sum = prizes_[ii];
            current_cluster.contains_root = (static_cast<int>(ii) == root_);
            current_cluster.moat = 0.0;

             log_message(4, "  Node %zu: Prize=%.9g, Active=%d, ContainsRoot=%d\n",
                         ii, current_cluster.prize_sum, current_cluster.active, current_cluster.contains_root);


            if (current_cluster.active) {
                if (prizes_[ii] > kEpsilon) {
                    log_message(4, "  Node %zu: Scheduling deactivation at time %.9g.\n", ii, prizes_[ii]);
                    clusters_deactivation_.insert(prizes_[ii], ii);
                } else {
                    log_message(4, "  Node %zu: Prize is negligible, immediately inactive.\n", ii);
                     current_cluster.active = false;
                     current_cluster.active_end_time = 0.0;
                }
            }
        }
        log_message(3, "initialize_clusters Exit.\n");
    }

    void PCSTFast::initialize_edges() {
        log_message(3, "initialize_edges Entry.\n");
        const size_t num_edges = edges_.size();
        const size_t num_nodes = prizes_.size();

        for (size_t ii = 0; ii < num_edges; ++ii) {
            const auto& edge = edges_[ii];
            int uu = edge.first;
            int vv = edge.second;
            ValueType cost = costs_[ii];
            log_message(4, "Initializing edge %zu: (%d <-> %d), Cost=%.9g\n", ii, uu, vv, cost);

            // --- Edge Validation ---
            if (uu < 0 || vv < 0 || static_cast<size_t>(uu) >= num_nodes || static_cast<size_t>(vv) >= num_nodes) {
                log_message(0, "Error: Edge %zu endpoint index (%d or %d) out of range (Num nodes: %zu).\n", ii, uu, vv, num_nodes);
                throw std::invalid_argument("Edge endpoint index out of range.");
            }
            if (cost < 0.0) {
                log_message(0, "Error: Edge %zu cost is negative (%.9g).\n", ii, cost);
                throw std::invalid_argument("Edge cost cannot be negative.");
            }
            if (uu == vv) {
                log_message(1, "Warning: Ignoring self-loop edge %zu (%d -> %d).\n", ii, uu, vv);
                edge_parts_[2 * ii].deleted = true;
                edge_parts_[2 * ii + 1].deleted = true;
                continue;
            }

            // --- Initialize Edge Parts ---
            EdgePart& uu_part = edge_parts_[2 * ii];
            EdgePart& vv_part = edge_parts_[2 * ii + 1];
            Cluster& uu_cluster = clusters_[uu];
            Cluster& vv_cluster = clusters_[vv];

            uu_part.deleted = false;
            vv_part.deleted = false;

            ValueType event_time_u, event_time_v;
            if (uu_cluster.active && vv_cluster.active) {
                event_time_u = cost * 0.5;
                event_time_v = cost * 0.5;
                 log_message(4, "  Edge %zu: Active-Active. Event times = %.9g\n", ii, event_time_u);
            } else if (uu_cluster.active) {
                event_time_u = cost;
                event_time_v = 0.0;
                 log_message(4, "  Edge %zu: Active(U)-Inactive(V). Event time U=%.9g, V=%.9g\n", ii, event_time_u, event_time_v);
            } else if (vv_cluster.active) {
                event_time_u = 0.0;
                event_time_v = cost;
                 log_message(4, "  Edge %zu: Inactive(U)-Active(V). Event time U=%.9g, V=%.9g\n", ii, event_time_u, event_time_v);
            } else {
                event_time_u = 0.0;
                event_time_v = 0.0;
                 log_message(4, "  Edge %zu: Inactive-Inactive. Event times = 0.0\n", ii);
            }

            uu_part.next_event_val = event_time_u;
            if (uu_cluster.active || vv_cluster.active) {
                 log_message(4, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", 2 * (int)ii, uu, uu, event_time_u);
                 uu_part.heap_node = uu_cluster.edge_parts.insert(uu_part.next_event_val, 2 * ii);
            } else {
                 log_message(4, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", 2 * (int)ii, uu);
            }


            vv_part.next_event_val = event_time_v;
             if (uu_cluster.active || vv_cluster.active) {
                 log_message(4, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", 2 * (int)ii + 1, vv, vv, event_time_v);
                 vv_part.heap_node = vv_cluster.edge_parts.insert(vv_part.next_event_val, 2 * ii + 1);
             } else {
                  log_message(4, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", 2 * (int)ii + 1, vv);
             }
        }
        log_message(3, "initialize_edges Exit.\n");
    }


    // --- Main Algorithm Loop ---

    bool PCSTFast::run(std::vector<int>* result_nodes, std::vector<int>* result_edges) {
        log_message(3, "PCSTFast::run Entry.\n");
        if (!result_nodes || !result_edges) {
             log_message(0, "Error: Result pointers cannot be null in run().\n");
            throw std::invalid_argument("Result pointers cannot be null.");
        }

        result_nodes->clear();
        result_edges->clear();
        phase1_result_.clear();
        phase2_result_.clear();
        phase3_result_.clear();
        stats_.reset(); // Ensure stats are clean for this run

        int current_num_active_clusters = 0;
        for(const auto& cluster : clusters_) {
            if (cluster.active && cluster.merged_into == kInvalidIndex) { // Count only active representatives
                current_num_active_clusters++;
            }
        }
        log_message(1, "Starting GW phase. Initial active clusters: %d (Target: %d)\n",
                    current_num_active_clusters, target_num_active_clusters_);

        // --- Main Event Loop ---
        int event_count = 0;
        while (current_num_active_clusters > target_num_active_clusters_) {
            event_count++;
            log_message(3, "--- GW Loop Iteration %d (Active Clusters: %d) ---\n", event_count, current_num_active_clusters);

            NextEvent event = get_next_event();
             log_message(3, "get_next_event returned: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                        static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);


            if (event.type == NextEvent::Type::NONE) {
                log_message(1, "No more events possible. Terminating GW loop.\n");
                break;
            }

            log_message(2, "-----------------------------------------\n");
             if (event.time < current_time_ - kEpsilon) {
                  log_message(0, "Error: Time decreased! Current=%.9g, Event=%.9g\n", current_time_, event.time);
                  // Optionally throw or break, as this indicates a major issue. Let's log and continue for now.
             }
            current_time_ = event.time; // Advance time
             log_message(2, "Advanced time to %.9g\n", current_time_);

            process_event(event);

            // Recalculate active cluster count after processing
            int previous_active_count = current_num_active_clusters;
            current_num_active_clusters = 0;
            for(const auto& cluster : clusters_) {
                if (cluster.active && cluster.merged_into == kInvalidIndex) {
                    current_num_active_clusters++;
                }
            }
             log_message(3, "End of event processing. Active clusters: %d (was %d)\n", current_num_active_clusters, previous_active_count);

        } // End while loop

        log_message(1, "Finished GW main loop: Final time %.9g, Total edge events processed %lld, Final active clusters %d\n",
                    current_time_, stats_.total_num_edge_events, current_num_active_clusters);

        // --- Pruning and Output Generation ---
         log_message(2, "--- Starting Pruning and Output Generation ---\n");

        // Phase 1: Determine 'good' nodes
         log_message(2, "Phase 1: Selecting initial good nodes...\n");
        node_good_.assign(prizes_.size(), 0);
        select_initial_active_clusters(node_good_);
         int good_node_count = 0;
         for(uint8_t flag : node_good_) good_node_count += flag;
         log_message(2, "Phase 1: Found %d good nodes.\n", good_node_count);

        if (pruning_ == PruningMethod::kNoPruning) {
            log_message(1, "Pruning: None. Using Phase 1 result directly.\n");
            build_node_set_from_edges(phase1_result_, result_nodes);
            *result_edges = phase1_result_;
            log_message(1, "Final Result (No Pruning): Nodes=%zu, Edges=%zu\n", result_nodes->size(), result_edges->size());
            log_message(3, "PCSTFast::run Exit (Success).\n");
            return true;
        }

        // Phase 2: Filter phase 1 edges
         log_message(2, "Phase 2: Building result based on good node connectivity...\n");
        build_phase2_result();
        log_message(1, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_.size());

        if (pruning_ == PruningMethod::kSimplePruning) {
            log_message(1, "Pruning: Simple. Using Phase 2 result.\n");
            build_node_set_from_edges(phase2_result_, result_nodes);
            *result_edges = phase2_result_;
            log_message(1, "Final Result (Simple Pruning): Nodes=%zu, Edges=%zu\n", result_nodes->size(), result_edges->size());
            log_message(3, "PCSTFast::run Exit (Success).\n");
            return true;
        }

        // Phase 3 (GW or Strong)
        log_message(2, "Phase 3: Building adjacency list for GW/Strong pruning (Size %zu edges)...\n", phase2_result_.size());
        build_phase3_adjacency();
        log_message(4, "Resetting node_deleted_ flags.\n");
        node_deleted_.assign(prizes_.size(), 0);

        if (pruning_ == PruningMethod::kGWPruning) {
            log_message(1, "Pruning: Starting GW Pruning.\n");
            prune_gw();
            *result_edges = phase3_result_;
            log_message(2, "GW pruning complete. Building final node set...\n");
            build_final_node_set(result_nodes);
            log_message(1, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes->size(), result_edges->size());

        } else if (pruning_ == PruningMethod::kStrongPruning) {
             log_message(1, "Pruning: Starting Strong Pruning.\n");
             prune_strong();
             log_message(2, "Strong pruning complete. Filtering phase 2 edges...\n");
             phase3_result_.clear();
             for(int edge_idx : phase2_result_) {
                 const auto& edge = edges_[edge_idx];
                 // Check bounds before accessing node_deleted_
                 bool u_del = (static_cast<size_t>(edge.first) >= node_deleted_.size() || node_deleted_[edge.first]);
                 bool v_del = (static_cast<size_t>(edge.second) >= node_deleted_.size() || node_deleted_[edge.second]);
                 if (!u_del && !v_del) {
                     phase3_result_.push_back(edge_idx);
                 } else {
                     log_message(4, "Strong: Edge %d (%d,%d) removed due to deleted endpoint.\n", edge_idx, edge.first, edge.second);
                 }
             }
            *result_edges = phase3_result_;
            log_message(2, "Strong pruning edge filtering done. Building final node set...\n");
            build_final_node_set(result_nodes);
            log_message(1, "Final Result (Strong Pruning): Nodes=%zu, Edges=%zu\n", result_nodes->size(), result_edges->size());
        }

        log_message(3, "PCSTFast::run Exit (Success).\n");
        return true;
    }

    // --- Event Handling ---

    PCSTFast::NextEvent PCSTFast::get_next_event() const {
        // Function body as provided in the previous response (with correct variable usage)
        // Adding entry/exit logs
        log_message(4, "get_next_event Entry.\n");

        NextEvent event; // Default type is NONE, time is infinity

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

        log_message(4, "Checking clusters_next_edge_event PQ...\n");
        if (clusters_next_edge_event_.get_min(&temp_pq_edge_time, &temp_pq_edge_cluster)) {
            log_message(4, "  Global Edge PQ min: Time=%.9g, Cluster=%d\n", temp_pq_edge_time, temp_pq_edge_cluster);
            if (temp_pq_edge_cluster != kInvalidIndex && static_cast<size_t>(temp_pq_edge_cluster) < clusters_.size()) {
                log_message(4, "  Checking cluster %d's internal heap...\n", temp_pq_edge_cluster);
                if (clusters_[temp_pq_edge_cluster].edge_parts.get_min(&temp_cluster_min_time, &temp_cluster_min_part)) {
                     log_message(4, "  Cluster %d internal min: Time=%.9g, Part=%d\n", temp_pq_edge_cluster, temp_cluster_min_time, temp_cluster_min_part);
                     next_edge_time = temp_cluster_min_time;
                     next_edge_cluster_index = temp_pq_edge_cluster;
                     next_edge_part_index = temp_cluster_min_part;
                     edge_event_found = true;
                } else {
                     log_message(0, "Error: Inconsistency in get_next_event. Global edge PQ has cluster %d, but its heap is empty.\n", temp_pq_edge_cluster);
                }
            } else {
                 log_message(0, "Error: Invalid cluster index %d retrieved from global edge PQ in get_next_event.\n", temp_pq_edge_cluster);
            }
        } else {
             log_message(4, "  Global Edge PQ is empty.\n");
        }

        log_message(4, "Checking clusters_deactivation PQ...\n");
        if (clusters_deactivation_.get_min(&temp_pq_cluster_time, &temp_pq_cluster_index)) {
            log_message(4, "  Deactivation PQ min: Time=%.9g, Cluster=%d\n", temp_pq_cluster_time, temp_pq_cluster_index);
             if (temp_pq_cluster_index != kInvalidIndex) {
                 next_cluster_time = temp_pq_cluster_time;
                 next_cluster_deactivation_index = temp_pq_cluster_index;
                 cluster_event_found = true;
             } else {
                 log_message(0, "Error: Invalid cluster index %d retrieved from deactivation PQ in get_next_event.\n", temp_pq_cluster_index);
             }
        } else {
             log_message(4, "  Deactivation PQ is empty.\n");
        }

        log_message(4, "Determining next event: EdgeFound=%d (Time=%.9g), ClusterFound=%d (Time=%.9g)\n",
                    edge_event_found, next_edge_time, cluster_event_found, next_cluster_time);
        if (edge_event_found && cluster_event_found) {
            if (next_edge_time <= next_cluster_time + kEpsilon * std::max(1.0, std::abs(next_cluster_time))) {
                log_message(4,"  Edge event is next (or tie).\n");
                event.type = NextEvent::Type::EDGE;
                event.time = next_edge_time;
                event.cluster_index = next_edge_cluster_index;
                event.edge_part_index = next_edge_part_index;
            } else {
                 log_message(4,"  Cluster event is next.\n");
                event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
                event.time = next_cluster_time;
                event.cluster_index = next_cluster_deactivation_index;
            }
        } else if (edge_event_found) {
             log_message(4,"  Only Edge event found.\n");
            event.type = NextEvent::Type::EDGE;
            event.time = next_edge_time;
            event.cluster_index = next_edge_cluster_index;
            event.edge_part_index = next_edge_part_index;
        } else if (cluster_event_found) {
             log_message(4,"  Only Cluster event found.\n");
            event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
            event.time = next_cluster_time;
            event.cluster_index = next_cluster_deactivation_index;
        } else {
             log_message(4,"  No events found.\n");
        }

        log_message(4, "get_next_event Exit. Returning Type=%d, Time=%.9g\n", static_cast<int>(event.type), event.time);
        return event;
    }

    void PCSTFast::process_event(const NextEvent& event) {
        log_message(3, "process_event Entry: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                    static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);

        if (event.type == NextEvent::Type::EDGE) {
            stats_.total_num_edge_events++;
            log_message(4, "Processing EDGE event for Cluster %d, EdgePart %d\n", event.cluster_index, event.edge_part_index);

            // --- Pre-checks ---
            if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size()) {
                 log_message(1, "Warning: Edge event cluster index %d invalid. Skipping.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index); // Attempt cleanup
                 return;
            }
             if (!clusters_[event.cluster_index].active) {
                 log_message(2, "Edge event for inactive cluster %d. Removing from global PQ. Skipping logic.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index);
                 // No need to delete from cluster heap, it shouldn't have been in global PQ.
                 return;
            }
            if (event.edge_part_index < 0 || static_cast<size_t>(event.edge_part_index) >= edge_parts_.size()) {
                 log_message(0, "Error: Invalid edge part index %d in edge event. Skipping logic.\n", event.edge_part_index);
                 return; // Cannot proceed
            }

            // --- Delete event from originating cluster's heap ---
             log_message(4, "Deleting min edge part from cluster %d's heap...\n", event.cluster_index);
             ValueType deleted_val; PayloadType deleted_part;
             if (!clusters_[event.cluster_index].edge_parts.delete_min(&deleted_val, &deleted_part)) {
                  log_message(0, "Error: Failed to delete min edge from cluster %d's heap (event time %.9g).\n", event.cluster_index, event.time);
                  clusters_next_edge_event_.delete_element(event.cluster_index); // Clean up global PQ
                  return;
             }
              log_message(4, "  Deleted Part=%d, Value=%.9g. Expected Part=%d, Value=%.9g.\n", deleted_part, deleted_val, event.edge_part_index, event.time);
             // Verify consistency (allow small float differences)
             if (deleted_part != event.edge_part_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                  log_message(0, "Error: Mismatch in edge event deletion for cluster %d!\n", event.cluster_index);
                   // Inconsistent state. Try to recover by updating global PQ anyway.
             }

            // --- Update global PQ for the cluster ---
             log_message(4, "Updating global edge PQ for cluster %d after deletion.\n", event.cluster_index);
             update_cluster_queues_post_growth(event.cluster_index); // Re-inserts cluster with new min edge time, if any

            // --- Check if specific edge part was logically deleted ---
            if (edge_parts_[event.edge_part_index].deleted) {
                stats_.num_deleted_edge_events++;
                log_message(2, "Edge part %d already deleted (flag=true), skipping event logic.\n", event.edge_part_index);
                return;
            }

            // --- Handle the core logic ---
            log_message(4, "Calling handle_edge_event for EdgePart %d.\n", event.edge_part_index);
            handle_edge_event(event.edge_part_index);

        } else if (event.type == NextEvent::Type::CLUSTER_DEACTIVATION) {
            stats_.num_cluster_events++;
             log_message(4, "Processing CLUSTER_DEACTIVATION event for Cluster %d.\n", event.cluster_index);

            // --- Remove the event from the deactivation PQ ---
             log_message(4, "Deleting min from cluster deactivation PQ...\n");
             ValueType deleted_val; IndexType deleted_cluster;
             if(!clusters_deactivation_.delete_min(&deleted_val, &deleted_cluster)) {
                 log_message(0, "Error: Failed to delete cluster deactivation event from PQ (expected cluster %d).\n", event.cluster_index);
                 // Event might have been removed by a merge earlier. Don't process if cluster is already inactive.
                 if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size() || !clusters_[event.cluster_index].active) {
                      log_message(2, "Skipping cluster deactivation logic as cluster %d is invalid or already inactive.\n", event.cluster_index);
                      return;
                 }
                 // Otherwise, proceed, but log the PQ inconsistency.
             }
              log_message(4, "  Deleted Cluster=%d, Value=%.9g. Expected Cluster=%d, Value=%.9g.\n", deleted_cluster, deleted_val, event.cluster_index, event.time);
             // Check consistency
             if(deleted_cluster != event.cluster_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                  log_message(1, "Warning: Mismatch in cluster deactivation event deletion (possibly due to prior merge). Proceeding with cluster %d.\n", event.cluster_index);
             }

            // --- Handle the core logic ---
            log_message(4, "Calling handle_cluster_deactivation_event for Cluster %d.\n", event.cluster_index);
            handle_cluster_deactivation_event(event.cluster_index);
        } else {
            log_message(0, "Error: Processed an event of type NONE.\n");
        }
         log_message(3, "process_event Exit.\n");
    }

    // --- Edge Event Sub-logic Implementation ---

    void PCSTFast::handle_edge_event(PayloadType edge_part_index) {
        log_message(3, "handle_edge_event Entry: EdgePart=%d\n", edge_part_index);

         // 1. Get representative clusters and sums
         log_message(4, "Getting edge processing info for part %d...\n", edge_part_index);
         std::optional<EdgeProcessingInfo> opt_info = get_edge_processing_info(edge_part_index);
         if (!opt_info) {
             log_message(1, "Warning: Could not get valid edge processing info for part %d. Skipping.\n", edge_part_index);
             return;
         }
         EdgeProcessingInfo info = *opt_info;

         log_message(2, "Edge Event Details: Edge %d (%d <-> %d), Cost %.9g\n",
                     info.edge_idx, edges_[info.edge_idx].first, edges_[info.edge_idx].second, info.cost);
         log_message(2, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                     info.current_part_idx, info.current_cluster_idx, info.sum_current, info.finished_moat_current);
         log_message(2, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                     info.other_part_idx, info.other_cluster_idx, info.sum_other, info.finished_moat_other);
         log_message(2, "  Remainder: %.9g\n", info.remainder);

         // 2. Check if clusters are already merged
         if (info.current_cluster_idx == info.other_cluster_idx) {
             stats_.num_merged_edge_events++;
             log_message(2, "Clusters %d already merged. Marking other edge part %d as deleted.\n",
                         info.current_cluster_idx, info.other_part_idx);
             edge_parts_[info.other_part_idx].deleted = true;
             log_message(3, "handle_edge_event Exit (Clusters already merged).\n");
             return;
         }

         // 3. Check merge condition
         // Use relative epsilon check too? Original just used absolute eps. Sticking to original for now.
         bool merge_condition = (info.remainder <= kEpsilon * info.cost || info.remainder <= kEpsilon );
         log_message(4, "Merge condition check: Remainder=%.9g, Epsilon*Cost=%.9g, Epsilon=%.9g -> Merge=%d\n",
                     info.remainder, kEpsilon * info.cost, kEpsilon, merge_condition);

         if (merge_condition) {
             stats_.total_num_merge_events++;
             log_message(2, "Merge condition met for edge %d.\n", info.edge_idx);
             log_message(3, "Adding edge %d to phase1_result.\n", info.edge_idx);
             phase1_result_.push_back(info.edge_idx);
             log_message(3, "Marking other edge part %d as deleted.\n", info.other_part_idx);
             edge_parts_[info.other_part_idx].deleted = true;

             log_message(4, "Calling merge_clusters...\n");
             merge_clusters(info);

         } else { // 4. Edge growth condition
             stats_.total_num_edge_growth_events++;
             log_message(2, "Growth condition for edge %d.\n", info.edge_idx);

             if (clusters_[info.other_cluster_idx].active) {
                 log_message(4, "Calling handle_active_active_growth...\n");
                 handle_active_active_growth(info);
             } else {
                 log_message(4, "Calling handle_active_inactive_growth...\n");
                 handle_active_inactive_growth(info);
             }
         }
         log_message(3, "handle_edge_event Exit.\n");
    }


     std::optional<PCSTFast::EdgeProcessingInfo> PCSTFast::get_edge_processing_info(PayloadType edge_part_index) {
         log_message(4, "get_edge_processing_info Entry: EdgePart=%d\n", edge_part_index);
        EdgeProcessingInfo info;
        info.current_part_idx = edge_part_index;
        info.other_part_idx = get_other_edge_part_index(edge_part_index);
        info.edge_idx = edge_part_index / 2;
         log_message(4, "  Derived OtherPart=%d, EdgeIdx=%d\n", info.other_part_idx, info.edge_idx);


        if (static_cast<size_t>(info.edge_idx) >= costs_.size() ||
            static_cast<size_t>(info.other_part_idx) >= edge_parts_.size())
        {
             log_message(0, "Error: Invalid edge index %d or part index %d derived from %d.\n",
                         info.edge_idx, info.other_part_idx, edge_part_index);
            return std::nullopt;
        }
        info.cost = costs_[info.edge_idx];
         log_message(4, "  EdgeCost=%.9g\n", info.cost);

        log_message(4, "  Getting sum for current part %d...\n", info.current_part_idx);
        get_sum_on_edge_part(info.current_part_idx, &info.sum_current, &info.finished_moat_current, &info.current_cluster_idx);
         log_message(4, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.current_cluster_idx, info.sum_current, info.finished_moat_current);

        log_message(4, "  Getting sum for other part %d...\n", info.other_part_idx);
        get_sum_on_edge_part(info.other_part_idx, &info.sum_other, &info.finished_moat_other, &info.other_cluster_idx);
         log_message(4, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.other_cluster_idx, info.sum_other, info.finished_moat_other);


        if (info.current_cluster_idx == kInvalidIndex || info.other_cluster_idx == kInvalidIndex ||
            static_cast<size_t>(info.current_cluster_idx) >= clusters_.size() ||
            static_cast<size_t>(info.other_cluster_idx) >= clusters_.size())
        {
            log_message(1, "Warning: Failed to find valid representatives in get_edge_processing_info.\n");
            return std::nullopt;
        }

        info.remainder = info.cost - info.sum_current - info.sum_other;
         log_message(4, "  Calculated Remainder = %.9g - %.9g - %.9g = %.9g\n", info.cost, info.sum_current, info.sum_other, info.remainder);

         log_message(4, "get_edge_processing_info Exit (Success).\n");
        return info;
     }


    void PCSTFast::merge_clusters(const EdgeProcessingInfo& info) {
         log_message(3, "merge_clusters Entry: Merging Cluster %d and Cluster %d via Edge %d.\n",
                     info.current_cluster_idx, info.other_cluster_idx, info.edge_idx);
         IndexType cluster1_idx = info.current_cluster_idx;
         IndexType cluster2_idx = info.other_cluster_idx;
         Cluster& cluster1 = clusters_[cluster1_idx];
         Cluster& cluster2 = clusters_[cluster2_idx];

         // --- Create New Merged Cluster ---
         IndexType new_cluster_idx = clusters_.size();
         log_message(4, "Creating new cluster at index %d.\n", new_cluster_idx);
         clusters_.emplace_back();
         Cluster& new_cluster = clusters_.back();

         // --- Update New Cluster State ---
         new_cluster.prize_sum = cluster1.prize_sum + cluster2.prize_sum;
         new_cluster.subcluster_moat_sum = cluster1.subcluster_moat_sum + cluster2.subcluster_moat_sum; // Base sum
         new_cluster.contains_root = cluster1.contains_root || cluster2.contains_root;
         new_cluster.active = !new_cluster.contains_root;
         new_cluster.merged_along = info.edge_idx;
         new_cluster.child_cluster_1 = cluster1_idx;
         new_cluster.child_cluster_2 = cluster2_idx;
          log_message(4, "  New cluster %d: PrizeSum=%.9g, BaseMoatSum=%.9g, ContainsRoot=%d, Active=%d\n",
                      new_cluster_idx, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, new_cluster.contains_root, new_cluster.active);

         // --- Update Child Cluster 1 State ---
         log_message(4, "Updating state for child cluster %d (was Active=%d).\n", cluster1_idx, cluster1.active);
         cluster1.merged_into = new_cluster_idx;
         if (cluster1.active) {
             cluster1.active = false;
             cluster1.active_end_time = current_time_;
             cluster1.moat = cluster1.active_end_time - cluster1.active_start_time;
             if (cluster1.moat < 0.0) {
                  log_message(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster1_idx, cluster1.moat);
                  cluster1.moat = 0.0;
             }
              log_message(4, "  Child %d became inactive. Moat=%.9g\n", cluster1_idx, cluster1.moat);
             new_cluster.subcluster_moat_sum += cluster1.moat; // Add calculated moat
         } else {
             log_message(4, "  Child %d already inactive. Adding previous moat %.9g.\n", cluster1_idx, cluster1.moat);
             new_cluster.subcluster_moat_sum += cluster1.moat; // Add previously calculated moat
         }

         // --- Update Child Cluster 2 State ---
         log_message(4, "Updating state for child cluster %d (was Active=%d).\n", cluster2_idx, cluster2.active);
         cluster2.merged_into = new_cluster_idx;
         if (cluster2.active) {
             stats_.num_active_active_merge_events++;
             cluster2.active = false;
             cluster2.active_end_time = current_time_;
             cluster2.moat = cluster2.active_end_time - cluster2.active_start_time;
             if (cluster2.moat < 0.0) {
                  log_message(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster2_idx, cluster2.moat);
                  cluster2.moat = 0.0;
             }
              log_message(4, "  Child %d (Active) became inactive. Moat=%.9g (Active-Active Merge)\n", cluster2_idx, cluster2.moat);
             new_cluster.subcluster_moat_sum += cluster2.moat;
         } else {
             stats_.num_active_inactive_merge_events++;
             log_message(4, "  Child %d already inactive. Adding previous moat %.9g (Active-Inactive Merge)\n", cluster2_idx, cluster2.moat);
             new_cluster.subcluster_moat_sum += cluster2.moat;
             if (!cluster2.contains_root) {
                  log_message(4, "  Logging inactive merge event for edge %d.\n", info.edge_idx);
                  log_inactive_merge_event(info);
                  ValueType edge_event_update_shift = current_time_ - cluster2.active_end_time;
                  if (edge_event_update_shift > kEpsilon) {
                      log_message(3, "  Shifting edge event times for inactive child cluster %d by %.9g.\n", cluster2_idx, edge_event_update_shift);
                      cluster2.edge_parts.add_to_heap(edge_event_update_shift);
                  }
             } else {
                  log_message(4, "  Skipping inactive merge event log/shift as inactive child %d contains root.\n", cluster2_idx);
             }
         }
          log_message(4, "  New cluster %d updated SubclusterMoatSum=%.9g\n", new_cluster_idx, new_cluster.subcluster_moat_sum);

        // --- Meld Pairing Heaps ---
        log_message(4, "Melding heaps from cluster %d and %d into %d.\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        new_cluster.edge_parts = PairingHeapType::meld(std::move(cluster1.edge_parts), std::move(cluster2.edge_parts));
         log_message(4, "  Meld complete. New cluster %d heap is%s empty.\n", new_cluster_idx, new_cluster.edge_parts.is_empty() ? "" : " NOT");


        // --- Update New Cluster Activity and PQs ---
        if (new_cluster.active) {
            new_cluster.active_start_time = current_time_;
            ValueType potential_deactivation_time = new_cluster.active_start_time + new_cluster.prize_sum - new_cluster.subcluster_moat_sum;
            log_message(4, "  New cluster %d is active. Potential deactivation time = %.9g + %.9g - %.9g = %.9g\n",
                        new_cluster_idx, new_cluster.active_start_time, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, potential_deactivation_time);
            if (potential_deactivation_time > current_time_ + kEpsilon) {
                 log_message(3, "  Scheduling deactivation for new cluster %d at time %.9g.\n", new_cluster_idx, potential_deactivation_time);
                 clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);
            } else {
                 log_message(3, "  New cluster %d immediately inactive (potential time %.9g <= current time %.9g).\n", new_cluster_idx, potential_deactivation_time, current_time_);
                 new_cluster.active = false;
                 new_cluster.active_end_time = current_time_;
                 new_cluster.moat = 0.0;
                 new_cluster.subcluster_moat_sum += new_cluster.moat; // Should be zero anyway
            }
        } else {
            new_cluster.active_start_time = current_time_;
            new_cluster.active_end_time = current_time_;
            new_cluster.moat = 0.0;
            new_cluster.subcluster_moat_sum += new_cluster.moat;
             log_message(3, "  New cluster %d contains root or became immediately inactive, setting Active=false.\n", new_cluster_idx);
        }

         log_message(4, "Updating PQs after merge (removing children %d, %d; adding new %d if needed).\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        update_cluster_queues_post_merge(new_cluster_idx, cluster1_idx, cluster2_idx);

         log_message(3, "merge_clusters Exit.\n");
    }

     void PCSTFast::log_inactive_merge_event(const EdgeProcessingInfo& info) {
         log_message(4, "log_inactive_merge_event Entry for Edge %d.\n", info.edge_idx);
         inactive_merge_events_.emplace_back();
         InactiveMergeEvent& merge_event = inactive_merge_events_.back();

         // Assume info.current corresponds to the active cluster side at event time
         merge_event.active_cluster_index = info.current_cluster_idx;
         merge_event.inactive_cluster_index = info.other_cluster_idx;

         const auto& edge = edges_[info.edge_idx];
         IndexType node_on_current_side = (info.current_part_idx % 2 == 0) ? edge.first : edge.second;
         IndexType node_on_other_side = (info.other_part_idx % 2 == 0) ? edge.first : edge.second;

         merge_event.active_cluster_node = node_on_current_side;
         merge_event.inactive_cluster_node = node_on_other_side;

         if (static_cast<size_t>(info.edge_idx) < edge_info_.size()) {
             edge_info_[info.edge_idx].inactive_merge_event = std::ssize(inactive_merge_events_) - 1;
             log_message(3, "Logged InactiveMergeEvent %zd for edge %d (Active: C%d/N%d, Inactive: C%d/N%d).\n",
                         inactive_merge_events_.size() - 1, info.edge_idx,
                         merge_event.active_cluster_index, merge_event.active_cluster_node,
                         merge_event.inactive_cluster_index, merge_event.inactive_cluster_node);
         } else {
              log_message(0, "Error: Edge index %d out of bounds for edge_info_ in log_inactive_merge_event.\n", info.edge_idx);
         }
          log_message(4, "log_inactive_merge_event Exit.\n");
     }


    void PCSTFast::handle_active_active_growth(const EdgeProcessingInfo& info) {
        stats_.num_active_active_edge_growth_events++;
        ValueType next_event_time = current_time_ + info.remainder * 0.5;
        log_message(3, "handle_active_active_growth Entry: Edge=%d, Clusters=(%d, %d), Remainder=%.9g, NextEventTime=%.9g\n",
                    info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);


        EdgePart& current_part = edge_parts_[info.current_part_idx];
        EdgePart& other_part = edge_parts_[info.other_part_idx];
        Cluster& current_cluster = clusters_[info.current_cluster_idx];
        Cluster& other_cluster = clusters_[info.other_cluster_idx];

        log_message(4, "Updating current part %d (cluster %d): Inserting into heap with time %.9g.\n",
                    info.current_part_idx, info.current_cluster_idx, next_event_time);
        current_part.next_event_val = next_event_time;
        current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
        update_cluster_queues_post_growth(info.current_cluster_idx);


        log_message(4, "Updating other part %d (cluster %d): Decreasing key in heap to time %.9g.\n",
                    info.other_part_idx, info.other_cluster_idx, next_event_time);
        other_part.next_event_val = next_event_time;
        log_message(4, "  Deleting cluster %d from global edge PQ before decrease_key.\n", info.other_cluster_idx);
        clusters_next_edge_event_.delete_element(info.other_cluster_idx);
        other_cluster.edge_parts.decrease_key(other_part.heap_node, next_event_time);
        update_cluster_queues_post_growth(info.other_cluster_idx);

         log_message(3, "handle_active_active_growth Exit.\n");
    }


    void PCSTFast::handle_active_inactive_growth(const EdgeProcessingInfo& info) {
         stats_.num_active_inactive_edge_growth_events++;
         ValueType next_event_time = current_time_ + info.remainder;
          log_message(3, "handle_active_inactive_growth Entry: Edge=%d, Clusters=(Active:%d, Inactive:%d), Remainder=%.9g, NextEventTime=%.9g\n",
                     info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);


         EdgePart& current_part = edge_parts_[info.current_part_idx];
         EdgePart& other_part = edge_parts_[info.other_part_idx];
         Cluster& current_cluster = clusters_[info.current_cluster_idx];
         Cluster& other_cluster = clusters_[info.other_cluster_idx];

         log_message(4, "Updating active part %d (cluster %d): Inserting into heap with time %.9g.\n",
                     info.current_part_idx, info.current_cluster_idx, next_event_time);
         current_part.next_event_val = next_event_time;
         current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
         update_cluster_queues_post_growth(info.current_cluster_idx);


         ValueType inactive_reset_time = other_cluster.active_end_time;
         log_message(4, "Updating inactive part %d (cluster %d): Decreasing key in heap to time %.9g (cluster active_end_time).\n",
                     info.other_part_idx, info.other_cluster_idx, inactive_reset_time);
         other_part.next_event_val = inactive_reset_time;
         other_cluster.edge_parts.decrease_key(other_part.heap_node, inactive_reset_time);
         // No global PQ update for inactive cluster

          log_message(3, "handle_active_inactive_growth Exit.\n");
    }

    // --- Cluster Deactivation Event ---

    void PCSTFast::handle_cluster_deactivation_event(IndexType cluster_index) {
         log_message(3, "handle_cluster_deactivation_event Entry: Cluster=%d\n", cluster_index);
        if (cluster_index < 0 || static_cast<size_t>(cluster_index) >= clusters_.size()) {
            log_message(0,"Error: Invalid cluster index %d in deactivation event.\n", cluster_index);
            return;
        }
        Cluster& cluster = clusters_[cluster_index];

        if (!cluster.active) {
            log_message(2, "Cluster %d already inactive, skipping deactivation event logic.\n", cluster_index);
            return; // Already deactivated
        }

        // --- Deactivate Cluster ---
        cluster.active = false;
        cluster.active_end_time = current_time_;
        cluster.moat = cluster.active_end_time - cluster.active_start_time;
        if (cluster.moat < 0.0) {
             log_message(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster_index, cluster.moat);
             cluster.moat = 0.0;
        }

        log_message(2, "Cluster Deactivation: Cluster %d at time %.9g (Moat size %.9g).\n",
                    cluster_index, current_time_, cluster.moat);

        // --- Remove from Edge Event PQ ---
        if (!cluster.edge_parts.is_empty()) {
             log_message(4, "Removing cluster %d from next edge event PQ as it became inactive.\n", cluster_index);
            clusters_next_edge_event_.delete_element(cluster_index);
        } else {
             log_message(4, "Cluster %d had no edges, no need to remove from edge PQ.\n", cluster_index);
        }

         log_message(3, "handle_cluster_deactivation_event Exit.\n");
    }

    // --- Priority Queue Updates ---

    void PCSTFast::update_cluster_queues_post_growth(IndexType cluster_idx) {
         log_message(4, "update_cluster_queues_post_growth Entry: Cluster=%d\n", cluster_idx);
        if (cluster_idx == kInvalidIndex || static_cast<size_t>(cluster_idx) >= clusters_.size() || !clusters_[cluster_idx].active) {
             log_message(4, "  Skipping PQ update (cluster invalid or inactive).\n");
            return;
        }
        Cluster& cluster = clusters_[cluster_idx];

        log_message(4, "  Deleting cluster %d from global edge PQ (if present).\n", cluster_idx);
        clusters_next_edge_event_.delete_element(cluster_idx);

        if (!cluster.edge_parts.is_empty()) {
            ValueType min_val; PayloadType min_part;
            if (cluster.edge_parts.get_min(&min_val, &min_part)) {
                 log_message(4, "  Cluster %d heap not empty. Inserting into edge PQ with Time=%.9g, Part=%d.\n", cluster_idx, min_val, min_part);
                clusters_next_edge_event_.insert(min_val, cluster_idx);
            } else {
                 log_message(0, "Error: Cluster %d heap not empty but get_min failed in update_queues.\n", cluster_idx);
            }
        } else {
             log_message(4, "  Cluster %d heap is empty, not inserting into edge PQ.\n", cluster_idx);
        }
         log_message(4, "update_cluster_queues_post_growth Exit.\n");
    }

     void PCSTFast::update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index) {
          log_message(4, "update_cluster_queues_post_merge Entry: New=%d, Children=(%d, %d)\n", new_cluster_index, child1_index, child2_index);

         log_message(4, "  Deleting child %d from PQs (if present).\n", child1_index);
         clusters_deactivation_.delete_element(child1_index);
         clusters_next_edge_event_.delete_element(child1_index);
          log_message(4, "  Deleting child %d from PQs (if present).\n", child2_index);
         clusters_deactivation_.delete_element(child2_index);
         clusters_next_edge_event_.delete_element(child2_index);

         if (new_cluster_index != kInvalidIndex && static_cast<size_t>(new_cluster_index) < clusters_.size()) {
             Cluster& new_cluster = clusters_[new_cluster_index];
             if (new_cluster.active && !new_cluster.edge_parts.is_empty()) {
                 ValueType min_val; PayloadType min_part;
                 if (new_cluster.edge_parts.get_min(&min_val, &min_part)) {
                     log_message(4, "  New cluster %d is active with edges. Inserting into edge PQ with Time=%.9g, Part=%d.\n", new_cluster_index, min_val, min_part);
                     clusters_next_edge_event_.insert(min_val, new_cluster_index);
                 } else {
                     log_message(0, "Error: New cluster %d heap not empty but get_min failed in update_queues_post_merge.\n", new_cluster_index);
                 }
             } else {
                 log_message(4, "  New cluster %d is inactive or has no edges, not inserting into edge PQ.\n", new_cluster_index);
             }
         } else {
              log_message(4, "  New cluster index %d is invalid, skipping edge PQ update.\n", new_cluster_index);
         }
           log_message(4, "update_cluster_queues_post_merge Exit.\n");
     }


    // --- Path Compression ---

    void PCSTFast::get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                                      ValueType* finished_moat_sum, IndexType* current_cluster_index)
    {
         log_message(4, "get_sum_on_edge_part Entry: EdgePart=%d\n", edge_part_index);
        const auto& edge = edges_[edge_part_index / 2];
        IndexType start_node = (edge_part_index % 2 == 0) ? edge.first : edge.second;
         log_message(4, "  Start node derived from edge part: %d\n", start_node);

        *total_sum = 0.0;
        *finished_moat_sum = 0.0;
        *current_cluster_index = kInvalidIndex;

        ValueType path_sum = 0.0;
         log_message(4, "  Calling find_representative_and_compress for start node %d...\n", start_node);
        IndexType representative_cluster = find_representative_and_compress(start_node, path_sum);
         log_message(4, "  find_representative result: RepCluster=%d, PathSum=%.9g\n", representative_cluster, path_sum);

        if (representative_cluster == kInvalidIndex || static_cast<size_t>(representative_cluster) >= clusters_.size()) {
            log_message(0, "Error: Path compression failed to find representative for node %d (edge part %d).\n", start_node, edge_part_index);
             log_message(4, "get_sum_on_edge_part Exit (Failure).\n");
            return;
        }

        *current_cluster_index = representative_cluster;
        const Cluster& rep_cluster = clusters_[representative_cluster];
         log_message(4, "  Representative cluster %d is Active=%d (StartTime=%.9g, EndTime=%.9g, Moat=%.9g)\n",
                     representative_cluster, rep_cluster.active, rep_cluster.active_start_time, rep_cluster.active_end_time, rep_cluster.moat);

        if (rep_cluster.active) {
            *finished_moat_sum = path_sum;
            *total_sum = path_sum + (current_time_ - rep_cluster.active_start_time);
             log_message(4, "  Rep is active. FinishedMoatSum=%.9g, TotalSum = %.9g + (%.9g - %.9g) = %.9g\n",
                         *finished_moat_sum, path_sum, current_time_, rep_cluster.active_start_time, *total_sum);
        } else {
            *total_sum = path_sum + rep_cluster.moat;
            *finished_moat_sum = *total_sum;
             log_message(4, "  Rep is inactive. FinishedMoatSum=TotalSum = %.9g + %.9g = %.9g\n",
                         path_sum, rep_cluster.moat, *total_sum);
        }
          log_message(4, "get_sum_on_edge_part Exit (Success). Cluster=%d, TotalSum=%.9g, FinishedMoat=%.9g\n",
                      *current_cluster_index, *total_sum, *finished_moat_sum);
    }


     PCSTFast::IndexType PCSTFast::find_representative_and_compress(IndexType start_node, ValueType& path_sum_out) {
          log_message(4, "find_representative_and_compress Entry: StartNode=%d\n", start_node);
         path_compression_visited_.clear();
         path_sum_out = 0.0;
         IndexType current_node = start_node;
         int steps = 0;

         while (static_cast<size_t>(current_node) < clusters_.size() && clusters_[current_node].merged_into != kInvalidIndex) {
             steps++;
             log_message(4, "  Step %d: CurrentNode=%d, PathSum=%.9g. Storing for compression.\n", steps, current_node, path_sum_out);
             path_compression_visited_.emplace_back(current_node, path_sum_out);

             ValueType step_sum;
             IndexType next_node;
             bool used_skip = false;

             if (clusters_[current_node].skip_up >= 0) {
                 step_sum = clusters_[current_node].skip_up_sum;
                 next_node = clusters_[current_node].skip_up;
                 used_skip = true;
             } else {
                 step_sum = clusters_[current_node].moat;
                 next_node = clusters_[current_node].merged_into;
             }
              log_message(4, "  Step %d: NextNode=%d, StepSum=%.9g (UsedSkip=%d)\n", steps, next_node, step_sum, used_skip);

             path_sum_out += step_sum;
             current_node = next_node;
         }
          log_message(4, "  Traversal finished after %d steps. FinalNode=%d, FinalPathSum=%.9g\n", steps, current_node, path_sum_out);

         if (static_cast<size_t>(current_node) >= clusters_.size()) {
              log_message(0, "Error: Path compression traversed out of bounds from start node %d.\n", start_node);
              return kInvalidIndex;
         }

         IndexType representative_cluster = current_node;
          log_message(4, "  Representative cluster found: %d\n", representative_cluster);
          log_message(4, "  Applying path compression for %zu visited nodes...\n", path_compression_visited_.size());

         for (const auto& [visited_node_idx, sum_at_visited] : path_compression_visited_) {
              if (static_cast<size_t>(visited_node_idx) < clusters_.size()) {
                  ValueType new_skip_sum = path_sum_out - sum_at_visited;
                  log_message(4, "    Updating skip pointer for Node %d: SkipUp=%d, SkipSum=%.9g\n",
                              visited_node_idx, representative_cluster, new_skip_sum);
                 clusters_[visited_node_idx].skip_up = representative_cluster;
                 clusters_[visited_node_idx].skip_up_sum = new_skip_sum;
              }
         }

          log_message(4, "find_representative_and_compress Exit. Rep=%d, PathSum=%.9g\n", representative_cluster, path_sum_out);
         return representative_cluster;
     }


    // --- Pruning Phase 1: Select Initial Good Nodes ---

    void PCSTFast::select_initial_active_clusters(std::vector<uint8_t>& node_good_flag) {
         log_message(3, "select_initial_active_clusters Entry.\n");
        cluster_queue_.clear();

        if (root_ >= 0) {
             log_message(3, "  Rooted case (Root=%d).\n", root_);
            IndexType root_cluster_final_index = root_;
             log_message(4, "  Finding final representative for root %d...\n", root_);
             int path_len = 0;
             while(static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].merged_into != kInvalidIndex) {
                 root_cluster_final_index = clusters_[root_cluster_final_index].merged_into;
                 path_len++;
             }
              log_message(4, "  Final representative for root %d is cluster %d (path length %d).\n", root_, root_cluster_final_index, path_len);

             if (static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].contains_root) {
                 log_message(2, "  Rooted case: Selecting nodes connected to final root cluster %d.\n", root_cluster_final_index);
                 cluster_queue_.push_back(root_cluster_final_index);
             } else {
                 log_message(0, "Warning: Could not find final cluster containing root %d. No nodes selected initially.\n", root_);
             }
        } else {
             log_message(3, "  Unrooted case.\n");
             log_message(2, "  Unrooted case: Selecting nodes from final active clusters.\n");
             for (size_t ii = 0; ii < clusters_.size(); ++ii) {
                 if (clusters_[ii].merged_into == kInvalidIndex && clusters_[ii].active) {
                      log_message(3, "  Adding final active cluster %zu to initial queue.\n", ii);
                     cluster_queue_.push_back(ii);
                 }
             }
        }

         log_message(3, "  Starting traversal from %zu initial clusters to mark good nodes...\n", cluster_queue_.size());
        size_t queue_index = 0;
        int nodes_marked_good = 0;
        while (queue_index < cluster_queue_.size()) {
            IndexType current_cluster_index = cluster_queue_[queue_index++];
            log_message(4, "    Processing cluster %d from queue (index %zu).\n", current_cluster_index, queue_index - 1);

            if (static_cast<size_t>(current_cluster_index) >= clusters_.size()) {
                 log_message(1, "Warning: Invalid cluster index %d encountered during good node traversal.\n", current_cluster_index);
                 continue;
            }

            const Cluster& current_cluster = clusters_[current_cluster_index];

            if (current_cluster.merged_along >= 0) {
                 log_message(4, "    Cluster %d is merged. Adding children %d, %d to queue.\n",
                             current_cluster_index, current_cluster.child_cluster_1, current_cluster.child_cluster_2);
                if(current_cluster.child_cluster_1 != kInvalidIndex) cluster_queue_.push_back(current_cluster.child_cluster_1);
                if(current_cluster.child_cluster_2 != kInvalidIndex) cluster_queue_.push_back(current_cluster.child_cluster_2);
            } else {
                 log_message(4,"    Cluster %d is original node cluster.\n", current_cluster_index);
                 if (static_cast<size_t>(current_cluster_index) < node_good_flag.size()) {
                     if(node_good_flag[current_cluster_index] == 0) {
                          log_message(3, "    Marking node %d as good.\n", current_cluster_index);
                          node_good_flag[current_cluster_index] = 1;
                          nodes_marked_good++;
                     } else {
                          log_message(4, "    Node %d already marked good.\n", current_cluster_index);
                     }
                 } else {
                      log_message(0, "Error: Original cluster index %d out of bounds for node_good_flag (%zu).\n", current_cluster_index, node_good_flag.size());
                 }
            }
        }
          log_message(3, "select_initial_active_clusters Exit. Marked %d nodes as good.\n", nodes_marked_good);
    }

    // --- Pruning Phase 2: Filter Edges ---

    void PCSTFast::build_phase2_result() {
         log_message(3, "build_phase2_result Entry (Filtering %zu phase 1 edges).\n", phase1_result_.size());
         phase2_result_.clear();
         int edges_kept = 0;
         int edges_removed = 0;
         for (int edge_idx : phase1_result_) {
             bool keep_edge = false;
             if (static_cast<size_t>(edge_idx) < edges_.size()) {
                 const auto& edge = edges_[edge_idx];
                 bool u_good = static_cast<size_t>(edge.first) < node_good_.size() && node_good_[edge.first];
                 bool v_good = static_cast<size_t>(edge.second) < node_good_.size() && node_good_[edge.second];
                  log_message(4, "  Checking Phase 1 edge %d (%d <-> %d): Node %d Good=%d, Node %d Good=%d\n",
                              edge_idx, edge.first, edge.second, edge.first, u_good, edge.second, v_good);
                 if (u_good && v_good) {
                      log_message(4, "    Keeping edge %d.\n", edge_idx);
                     phase2_result_.push_back(edge_idx);
                     edges_kept++;
                     keep_edge = true;
                 }
             } else {
                  log_message(1, "Warning: Invalid edge index %d found in phase1_result_.\n", edge_idx);
             }
             if (!keep_edge) {
                  edges_removed++;
                   log_message(4, "    Removing edge %d.\n", edge_idx);
             }
         }
          log_message(3, "build_phase2_result Exit. Kept %d edges, Removed %d edges.\n", edges_kept, edges_removed);
    }


    // --- Pruning Phase 3 Setup ---

    void PCSTFast::build_phase3_adjacency() {
         log_message(3, "build_phase3_adjacency Entry (Using %zu phase 2 edges).\n", phase2_result_.size());
         const size_t num_nodes = prizes_.size();
         if (phase3_neighbors_.size() != num_nodes) {
              log_message(3, "  Resizing phase3_neighbors_ to %zu.\n", num_nodes);
              phase3_neighbors_.resize(num_nodes);
         }
          log_message(4, "  Clearing existing adjacency lists.\n");
         for (auto& neighbors : phase3_neighbors_) { neighbors.clear(); }

         int edges_added = 0;
         for (int edge_idx : phase2_result_) {
             if (static_cast<size_t>(edge_idx) < edges_.size() && static_cast<size_t>(edge_idx) < costs_.size()) {
                 const auto& edge = edges_[edge_idx];
                 ValueType cost = costs_[edge_idx];
                 if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                      log_message(4, "  Adding edge %d (%d <-> %d, Cost=%.9g) to adjacency list.\n",
                                  edge_idx, edge.first, edge.second, cost);
                     phase3_neighbors_[edge.first].emplace_back(edge.second, cost);
                     phase3_neighbors_[edge.second].emplace_back(edge.first, cost);
                     edges_added++;
                 } else {
                      log_message(0, "Error: Invalid node index in edge %d while building phase 3 adjacency list.\n", edge_idx);
                 }
             } else {
                 log_message(1, "Warning: Invalid edge index %d found in phase2_result_ during adjacency build.\n", edge_idx);
             }
         }
          log_message(3, "build_phase3_adjacency Exit. Added %d edges to lists.\n", edges_added);
    }


    // --- Pruning Phase 3: GW ---

    void PCSTFast::prune_gw() {
         log_message(3, "prune_gw Entry.\n");
        phase3_result_.clear();
         log_message(4, "Resetting cluster necessary flags.\n");
        for(auto& cluster : clusters_) { cluster.necessary = false; }

        log_message(2, "Starting GW pruning traversal (processing %zu edges in reverse order).\n", phase2_result_.size());
        int edges_kept = 0;
        int edges_discarded = 0;

        for (int ii = std::ssize(phase2_result_) - 1; ii >= 0; --ii) {
            int current_edge_index = phase2_result_[ii];
             log_message(4, "Processing edge %d (index %d from end).\n", current_edge_index, (int)phase2_result_.size()-1-ii);

            if(static_cast<size_t>(current_edge_index) >= edges_.size() ||
               static_cast<size_t>(current_edge_index) >= edge_info_.size()) {
                log_message(0,"Error: Invalid edge index %d during GW pruning.\n", current_edge_index);
                continue;
            }

            const auto& edge = edges_[current_edge_index];
            int uu = edge.first;
            int vv = edge.second;
             log_message(4, "  Edge %d connects nodes (%d, %d).\n", current_edge_index, uu, vv);

            bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
            bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);
             log_message(4, "  Node status: %d Deleted=%d, %d Deleted=%d.\n", uu, u_deleted, vv, v_deleted);

            if (u_deleted && v_deleted) {
                log_message(3, "GW: Both endpoints (%d, %d) of edge %d deleted. Skipping.\n", uu, vv, current_edge_index);
                edges_discarded++;
                continue;
            }

            int inactive_merge_idx = edge_info_[current_edge_index].inactive_merge_event;
             log_message(4, "  Inactive merge event index for edge %d: %d.\n", current_edge_index, inactive_merge_idx);


            if (inactive_merge_idx == kInvalidIndex) {
                 log_message(3, "GW: Edge %d (%d, %d) was Active-Active. Keeping.\n", current_edge_index, uu, vv);
                 phase3_result_.push_back(current_edge_index);
                 edges_kept++;
                  log_message(4, "  Marking clusters necessary from nodes %d and %d (if not deleted).\n", uu, vv);
                 if (!u_deleted) mark_clusters_as_necessary_gw(uu);
                 if (!v_deleted) mark_clusters_as_necessary_gw(vv);
            } else {
                 if (static_cast<size_t>(inactive_merge_idx) >= inactive_merge_events_.size()) {
                     log_message(0,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, current_edge_index);
                     continue;
                 }
                 const InactiveMergeEvent& merge_event = inactive_merge_events_[inactive_merge_idx];
                 IndexType inactive_cluster_rep = merge_event.inactive_cluster_index;
                 IndexType active_node = merge_event.active_cluster_node;
                 IndexType inactive_node = merge_event.inactive_cluster_node;
                  log_message(4, "  Edge %d was Active-Inactive. ActiveNode=%d, InactiveNode=%d, InactiveRep=%d.\n",
                              current_edge_index, active_node, inactive_node, inactive_cluster_rep);


                 bool is_inactive_necessary = (static_cast<size_t>(inactive_cluster_rep) < clusters_.size() &&
                                               clusters_[inactive_cluster_rep].necessary);
                  log_message(4, "  Checking necessity of inactive representative cluster %d: Necessary=%d\n",
                              inactive_cluster_rep, is_inactive_necessary);

                 if (is_inactive_necessary) {
                      log_message(3, "GW: Edge %d (%d, %d) was A-I, Inactive side cluster %d is necessary. Keeping.\n", current_edge_index, uu, vv, inactive_cluster_rep);
                      phase3_result_.push_back(current_edge_index);
                      edges_kept++;
                       log_message(4, "  Marking clusters necessary from nodes %d and %d.\n", active_node, inactive_node);
                      mark_clusters_as_necessary_gw(active_node);
                      mark_clusters_as_necessary_gw(inactive_node);
                 } else {
                      log_message(3, "GW: Edge %d (%d, %d) was A-I, Inactive side cluster %d not necessary. Pruning inactive side starting at node %d.\n", current_edge_index, uu, vv, inactive_cluster_rep, inactive_node);
                      edges_discarded++;
                      mark_nodes_as_deleted_gw(inactive_node, active_node);
                 }
            }
        }

        log_message(2, "GW pruning reverse pass complete. Reversing phase 3 result.\n");
        std::reverse(phase3_result_.begin(), phase3_result_.end());
         log_message(3, "prune_gw Exit. Kept %d edges, Discarded initially %d edges.\n", edges_kept, edges_discarded);
    }


     void PCSTFast::mark_clusters_as_necessary_gw(IndexType start_node_index) {
          log_message(4, "mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);
         IndexType current_cluster_index = start_node_index;
         int steps = 0;
         while (current_cluster_index != kInvalidIndex && static_cast<size_t>(current_cluster_index) < clusters_.size()) {
             steps++;
             if (clusters_[current_cluster_index].necessary) {
                  log_message(4, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
                 break;
             }
              log_message(4, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);
             clusters_[current_cluster_index].necessary = true;
             current_cluster_index = clusters_[current_cluster_index].merged_into;
              log_message(4, "  Moving up to parent cluster %d.\n", current_cluster_index);
         }
           log_message(4, "mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
     }

     void PCSTFast::mark_nodes_as_deleted_gw(IndexType start_node_index, IndexType parent_node_index) {
          log_message(4, "mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);
         if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
              log_message(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
             return;
         }

         cluster_queue_.clear();
         cluster_queue_.push_back(start_node_index);
         node_deleted_[start_node_index] = 1;
          log_message(3, "  Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);

         size_t queue_index = 0;
         int nodes_deleted_count = 1;
         while (queue_index < cluster_queue_.size()) {
             IndexType current_node_index = cluster_queue_[queue_index++];
             log_message(4,"    Processing node %d from deletion queue (index %zu).\n", current_node_index, queue_index-1);


             if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                 log_message(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
                 for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                     log_message(4,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);
                     if (neighbor_node_index == parent_node_index) {
                          log_message(4,"          Neighbor is parent, skipping.\n");
                         continue;
                     }
                     if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                         node_deleted_[neighbor_node_index] = 1;
                         cluster_queue_.push_back(neighbor_node_index);
                         nodes_deleted_count++;
                         log_message(3, "        Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                     } else {
                           log_message(4,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
                     }
                 }
             } else {
                  log_message(4, "      Node %d has no neighbors in phase 3 graph (or index out of bounds).\n", current_node_index);
             }
         }
         log_message(4, "mark_nodes_as_deleted_gw Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
     }


    // --- Pruning Phase 3: Strong ---

    void PCSTFast::prune_strong() {
         log_message(3, "prune_strong Entry.\n");
         log_message(2, "Strong Pruning: Labeling final components...\n");
         label_final_components();
         log_message(2, "Strong Pruning: Found %zu connected components.\n", final_components_.size());

         for (size_t comp_idx = 0; comp_idx < final_components_.size(); ++comp_idx) {
             if (final_components_[comp_idx].empty()) {
                 log_message(3, "Skipping empty component %zu.\n", comp_idx);
                 continue;
             }

             log_message(2, "Strong Pruning: Processing component %zu (size %zu).\n", comp_idx, final_components_[comp_idx].size());

             IndexType component_root_node;
             if (static_cast<int>(comp_idx) == root_component_index_) {
                 component_root_node = root_;
                 log_message(3, "  Component contains designated root %d. Using as pruning root.\n", root_);
             } else {
                 log_message(3, "  Finding best root for component %zu...\n", comp_idx);
                 component_root_node = find_best_component_root(comp_idx);
                 log_message(3, "  Best root for component %zu is node %d. Using as pruning root.\n", comp_idx, component_root_node);
             }

             if (component_root_node != kInvalidIndex && static_cast<size_t>(component_root_node) < prizes_.size()) {
                  log_message(3, "  Running strong pruning DFS from root %d (MarkDeleted=true).\n", component_root_node);
                 // Reset parent/payoff state *specifically* for this component before DFS
                  log_message(4, "  Resetting strong pruning parent/payoff for component %zu nodes.\n", comp_idx);
                 for(IndexType node_idx : final_components_[comp_idx]) {
                     if(static_cast<size_t>(node_idx) < strong_pruning_parent_.size()) strong_pruning_parent_[node_idx] = {kInvalidIndex, 0.0};
                     if(static_cast<size_t>(node_idx) < strong_pruning_payoff_.size()) strong_pruning_payoff_[node_idx] = 0.0;
                 }
                 strong_pruning_dfs(component_root_node, true); // True = mark deleted nodes
             } else {
                 log_message(1, "Warning: Could not determine valid root for component %zu. Skipping strong pruning for it.\n", comp_idx);
             }
         }
         log_message(3, "prune_strong Exit.\n");
    }


     void PCSTFast::label_final_components() {
         log_message(3, "label_final_components Entry.\n");
         const size_t num_nodes = prizes_.size();
         if (final_component_label_.size() != num_nodes) final_component_label_.resize(num_nodes);
         std::fill(final_component_label_.begin(), final_component_label_.end(), kInvalidIndex);
         final_components_.clear();
         root_component_index_ = kInvalidIndex;
         int components_found = 0;

         for (IndexType start_node = 0; start_node < static_cast<IndexType>(num_nodes); ++start_node) {
             bool is_in_phase2_graph = (static_cast<size_t>(start_node) < phase3_neighbors_.size() && !phase3_neighbors_[start_node].empty())
                                      || (static_cast<size_t>(start_node) < node_good_.size() && node_good_[start_node]
                                          && static_cast<size_t>(start_node) < phase3_neighbors_.size() && phase3_neighbors_[start_node].empty());

             if (is_in_phase2_graph && final_component_label_[start_node] == kInvalidIndex) {
                 IndexType new_component_id = final_components_.size();
                 components_found++;
                 log_message(3, "  Found new component %d starting from node %d.\n", new_component_id, start_node);
                 final_components_.emplace_back();
                 label_component_recursive(start_node, new_component_id);
                 log_message(3, "  Component %d labeled. Size=%zu.\n", new_component_id, final_components_[new_component_id].size());
             }
         }
          log_message(3, "label_final_components Exit. Found %d components.\n", components_found);
     }

     void PCSTFast::label_component_recursive(IndexType start_node_index, IndexType component_id) {
          log_message(4, "label_component_recursive Entry: StartNode=%d, ComponentID=%d\n", start_node_index, component_id);
         cluster_queue_.clear();
         cluster_queue_.push_back(start_node_index);
         final_component_label_[start_node_index] = component_id;
         final_components_[component_id].push_back(start_node_index);
         log_message(4, "  Added start node %d to component %d.\n", start_node_index, component_id);
         if (start_node_index == root_) {
             root_component_index_ = component_id;
             log_message(4, "  Start node is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
         }

         size_t queue_next = 0;
         while (queue_next < cluster_queue_.size()) {
             IndexType current_node_index = cluster_queue_[queue_next++];
              log_message(4,"    Processing node %d from component labeling queue (index %zu).\n", current_node_index, queue_next-1);

             if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                 log_message(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
                 for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                     if (static_cast<size_t>(neighbor_node_index) < final_component_label_.size() &&
                         final_component_label_[neighbor_node_index] == kInvalidIndex)
                     {
                         log_message(4,"        Labeling neighbor %d with component %d and adding to queue.\n", neighbor_node_index, component_id);
                         final_component_label_[neighbor_node_index] = component_id;
                         final_components_[component_id].push_back(neighbor_node_index);
                         cluster_queue_.push_back(neighbor_node_index);
                         if (neighbor_node_index == root_) {
                              root_component_index_ = component_id;
                               log_message(4, "        Neighbor is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
                         }
                     } else {
                           log_message(4,"        Neighbor %d already labeled or invalid index.\n", neighbor_node_index);
                     }
                 }
             }
         }
          log_message(4, "label_component_recursive Exit for Component %d.\n", component_id);
     }


     PCSTFast::IndexType PCSTFast::find_best_component_root(IndexType component_index) {
          log_message(3, "find_best_component_root Entry: Component=%d\n", component_index);
         if (static_cast<size_t>(component_index) >= final_components_.size() || final_components_[component_index].empty()) {
              log_message(1, "Warning: Invalid or empty component index %d in find_best_component_root.\n", component_index);
             return kInvalidIndex;
         }
         const auto& component_nodes = final_components_[component_index];
         IndexType initial_root = component_nodes[0];
         log_message(3, "  Using initial root %d for component %d.\n", initial_root, component_index);

          log_message(4, "  Resetting strong pruning parent/payoff for component %d nodes.\n", component_index);
         for (IndexType node_idx : component_nodes) {
             if (static_cast<size_t>(node_idx) < strong_pruning_parent_.size()) strong_pruning_parent_[node_idx] = {kInvalidIndex, 0.0};
             if (static_cast<size_t>(node_idx) < strong_pruning_payoff_.size()) strong_pruning_payoff_[node_idx] = 0.0;
         }

          log_message(4, "  Running initial strong_pruning_dfs from %d (MarkDeleted=false).\n", initial_root);
         strong_pruning_dfs(initial_root, false);

         if (static_cast<size_t>(initial_root) >= strong_pruning_payoff_.size()) {
               log_message(0,"Error: Initial root index %d out of bounds for payoff after DFS.\n", initial_root);
               return kInvalidIndex;
         }

         IndexType current_best_root = initial_root;
         ValueType current_best_value = strong_pruning_payoff_[initial_root];
         log_message(3, "  Initial Payoff at root %d = %.9g\n", current_best_root, current_best_value);

          log_message(4, "  Starting payoff propagation (rerooting) from initial root %d.\n", initial_root);
         propagate_payoffs_and_find_best(initial_root, current_best_root, current_best_value);

          log_message(3, "find_best_component_root Exit. Best Root=%d, Best Payoff=%.9g\n", current_best_root, current_best_value);
         return current_best_root;
     }

     void PCSTFast::propagate_payoffs_and_find_best(IndexType initial_root, IndexType& best_root_out, ValueType& best_value_out) {
         log_message(4, "propagate_payoffs_and_find_best Entry: InitialRoot=%d\n", initial_root);
         strong_pruning_stack2_.clear();

         if (static_cast<size_t>(initial_root) < phase3_neighbors_.size()) {
              log_message(4, "  Adding children of initial root %d to propagation stack...\n", initial_root);
             for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[initial_root]) {
                 if (static_cast<size_t>(neighbor_node_index) < strong_pruning_parent_.size() &&
                     strong_pruning_parent_[neighbor_node_index].first == initial_root)
                 {
                      log_message(4, "    Adding child %d to stack.\n", neighbor_node_index);
                     strong_pruning_stack2_.push_back(neighbor_node_index);
                 }
             }
         }

         int nodes_processed = 0;
         while (!strong_pruning_stack2_.empty()) {
             IndexType current_node_index = strong_pruning_stack2_.back();
             strong_pruning_stack2_.pop_back();
             nodes_processed++;
             log_message(4,"    Processing node %d from propagation stack (item %d).\n", current_node_index, nodes_processed);


             if (static_cast<size_t>(current_node_index) >= prizes_.size() ||
                 static_cast<size_t>(current_node_index) >= strong_pruning_parent_.size() ||
                 static_cast<size_t>(current_node_index) >= strong_pruning_payoff_.size()) {
                      log_message(1,"Warning: Invalid index %d during payoff propagation.\n", current_node_index);
                      continue;
                 }

             IndexType parent_index = strong_pruning_parent_[current_node_index].first;
             ValueType parent_edge_cost = strong_pruning_parent_[current_node_index].second;
             log_message(4,"      Parent=%d, EdgeCost=%.9g\n", parent_index, parent_edge_cost);


             if (parent_index == kInvalidIndex || static_cast<size_t>(parent_index) >= strong_pruning_payoff_.size()) {
                 log_message(1,"Warning: Invalid parent index %d for node %d during payoff propagation.\n", parent_index, current_node_index);
                 continue;
             }

             ValueType payoff_subtree_current = strong_pruning_payoff_[current_node_index];
             ValueType payoff_parent_rooted = strong_pruning_payoff_[parent_index];
              log_message(4,"      Payoffs: Subtree(Current)=%.9g, Parent(Rooted)=%.9g\n", payoff_subtree_current, payoff_parent_rooted);


             ValueType current_contribution_to_parent = 0.0;
             if (payoff_subtree_current > parent_edge_cost) {
                 current_contribution_to_parent = payoff_subtree_current - parent_edge_cost;
             }
              log_message(4,"      Current contribution to parent = max(0, %.9g - %.9g) = %.9g\n", payoff_subtree_current, parent_edge_cost, current_contribution_to_parent);

             ValueType payoff_parent_without_current = payoff_parent_rooted - current_contribution_to_parent;
              log_message(4,"      Parent payoff without current = %.9g - %.9g = %.9g\n", payoff_parent_rooted, current_contribution_to_parent, payoff_parent_without_current);

             ValueType parent_contribution_to_current = 0.0;
             if (payoff_parent_without_current > parent_edge_cost) {
                 parent_contribution_to_current = payoff_parent_without_current - parent_edge_cost;
             }
               log_message(4,"      Parent contribution to current = max(0, %.9g - %.9g) = %.9g\n", payoff_parent_without_current, parent_edge_cost, parent_contribution_to_current);


             ValueType current_node_total_payoff_if_root = payoff_subtree_current + parent_contribution_to_current;
              log_message(4,"      Total payoff if node %d is root = %.9g + %.9g = %.9g\n", current_node_index, payoff_subtree_current, parent_contribution_to_current, current_node_total_payoff_if_root);


             if (current_node_total_payoff_if_root > best_value_out) {
                  log_message(3, "      New best root found: Node=%d, Payoff=%.9g (OldBest: Node=%d, Payoff=%.9g)\n",
                              current_node_index, current_node_total_payoff_if_root, best_root_out, best_value_out);
                 best_root_out = current_node_index;
                 best_value_out = current_node_total_payoff_if_root;
             }

             strong_pruning_payoff_[current_node_index] = current_node_total_payoff_if_root; // Update payoff for propagation

             if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                  log_message(4,"      Adding children of node %d to propagation stack...\n", current_node_index);
                 for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                     if (static_cast<size_t>(neighbor_node_index) < strong_pruning_parent_.size() &&
                         strong_pruning_parent_[neighbor_node_index].first == current_node_index)
                     {
                           log_message(4,"        Adding child %d to stack.\n", neighbor_node_index);
                         strong_pruning_stack2_.push_back(neighbor_node_index);
                     }
                 }
             }
         }
          log_message(4, "propagate_payoffs_and_find_best Exit. Processed %d nodes.\n", nodes_processed);
     }


     void PCSTFast::strong_pruning_dfs(IndexType start_node_index, bool mark_deleted) {
          log_message(3, "strong_pruning_dfs Entry: StartNode=%d, MarkDeleted=%d\n", start_node_index, mark_deleted);
         strong_pruning_stack_.clear();
         const size_t num_nodes = prizes_.size();

         if (static_cast<size_t>(start_node_index) >= num_nodes ||
             static_cast<size_t>(start_node_index) >= strong_pruning_parent_.size()) {
              log_message(1,"Warning: Invalid start node %d for strong_pruning_dfs.\n", start_node_index);
              return;
         }

         log_message(4,"  Initializing DFS root %d: Parent=Invalid, Pushing PreOrder.\n", start_node_index);
         strong_pruning_parent_[start_node_index] = {kInvalidIndex, 0.0};
         strong_pruning_stack_.emplace_back(true, start_node_index);

         int nodes_visited_pre = 0;
         int nodes_visited_post = 0;
         while (!strong_pruning_stack_.empty()) {
             bool is_pre_order = strong_pruning_stack_.back().first;
             IndexType current_node_index = strong_pruning_stack_.back().second;
             strong_pruning_stack_.pop_back();

             if (static_cast<size_t>(current_node_index) >= num_nodes) {
                  log_message(1,"Warning: Invalid node index %d popped from DFS stack.\n", current_node_index);
                  continue;
             }

             if (is_pre_order) {
                 nodes_visited_pre++;
                 log_message(4,"  DFS PreOrder Visit: Node=%d (Visit %d)\n", current_node_index, nodes_visited_pre);
                 log_message(4,"    Pushing PostOrder visit for node %d.\n", current_node_index);
                 strong_pruning_stack_.emplace_back(false, current_node_index);
                 strong_pruning_payoff_[current_node_index] = prizes_[current_node_index];
                  log_message(4,"    Initialized Payoff[%d] = Prize = %.9g\n", current_node_index, strong_pruning_payoff_[current_node_index]);

                 if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                      log_message(4,"    Exploring neighbors of %d...\n", current_node_index);
                     for (const auto& [neighbor_node_index, neighbor_cost] : phase3_neighbors_[current_node_index]) {
                           log_message(4,"      Neighbor %d (Parent is %d).\n", neighbor_node_index, strong_pruning_parent_[current_node_index].first);
                          if (neighbor_node_index == strong_pruning_parent_[current_node_index].first) {
                              log_message(4,"        Is parent, skipping.\n");
                              continue;
                          }
                          if (static_cast<size_t>(neighbor_node_index) >= num_nodes) {
                               log_message(1,"Warning: Invalid neighbor index %d for node %d.\n", neighbor_node_index, current_node_index);
                               continue;
                          }

                           log_message(4,"      Setting Parent[%d] = %d (Cost=%.9g). Pushing PreOrder.\n",
                                       neighbor_node_index, current_node_index, neighbor_cost);
                          strong_pruning_parent_[neighbor_node_index] = {current_node_index, neighbor_cost};
                          strong_pruning_stack_.emplace_back(true, neighbor_node_index);
                     }
                 }
             } else {
                 nodes_visited_post++;
                 log_message(4,"  DFS PostOrder Visit: Node=%d (Visit %d)\n", current_node_index, nodes_visited_post);

                 if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                      log_message(4,"    Aggregating child payoffs for node %d...\n", current_node_index);
                     for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                         // Check if neighbor is a child in the DFS tree
                         if (static_cast<size_t>(neighbor_node_index) < strong_pruning_parent_.size() &&
                             strong_pruning_parent_[neighbor_node_index].first == current_node_index)
                         {
                             log_message(4,"      Processing child %d...\n", neighbor_node_index);
                             ValueType child_edge_cost = strong_pruning_parent_[neighbor_node_index].second;
                             ValueType child_subtree_payoff = strong_pruning_payoff_[neighbor_node_index] - child_edge_cost;
                              log_message(4,"        ChildPayoff=%.9g, EdgeCost=%.9g -> NetChildPayoff=%.9g\n",
                                          strong_pruning_payoff_[neighbor_node_index], child_edge_cost, child_subtree_payoff);

                             if (child_subtree_payoff <= 0.0) {
                                  log_message(3,"        Child %d subtree payoff %.9g <= 0.\n", neighbor_node_index, child_subtree_payoff);
                                 if (mark_deleted) {
                                     log_message(3, "        Strong Pruning: Pruning subtree at %d.\n", neighbor_node_index);
                                     mark_nodes_as_deleted_strong(neighbor_node_index, current_node_index);
                                 } else {
                                      log_message(4,"        (Not marking deleted as MarkDeleted=false)\n");
                                 }
                             } else {
                                 log_message(4,"        Adding child %d payoff %.9g to Payoff[%d] (current=%.9g).\n",
                                             neighbor_node_index, child_subtree_payoff, current_node_index, strong_pruning_payoff_[current_node_index]);
                                 strong_pruning_payoff_[current_node_index] += child_subtree_payoff;
                             }
                         }
                     }
                      log_message(4,"    Final Payoff[%d] = %.9g after aggregating children.\n", current_node_index, strong_pruning_payoff_[current_node_index]);
                 }
             }
         }
          log_message(3, "strong_pruning_dfs Exit. Visited %d nodes (pre), %d nodes (post).\n", nodes_visited_pre, nodes_visited_post);
     }


     void PCSTFast::mark_nodes_as_deleted_strong(IndexType start_node_index, IndexType parent_node_index) {
          log_message(4, "mark_nodes_as_deleted_strong Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);
         // This function is identical to mark_nodes_as_deleted_gw.
         // Consider consolidating later if they remain identical.
         if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
              log_message(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
             return;
         }

         cluster_queue_.clear();
         cluster_queue_.push_back(start_node_index);
         node_deleted_[start_node_index] = 1;
          log_message(3, "  Strong: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);


         size_t queue_index = 0;
         int nodes_deleted_count = 1;
         while (queue_index < cluster_queue_.size()) {
             IndexType current_node_index = cluster_queue_[queue_index++];
              log_message(4,"    Processing node %d from strong deletion queue (index %zu).\n", current_node_index, queue_index-1);

             if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                 log_message(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
                 for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                      log_message(4,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);
                     if (neighbor_node_index == parent_node_index) {
                          log_message(4,"          Neighbor is parent, skipping.\n");
                         continue;
                     }
                     if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                         node_deleted_[neighbor_node_index] = 1;
                         cluster_queue_.push_back(neighbor_node_index);
                         nodes_deleted_count++;
                          log_message(3, "        Strong: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                     } else {
                           log_message(4,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
                     }
                 }
             } else {
                  log_message(4, "      Node %d has no neighbors in phase 3 graph (or index out of bounds).\n", current_node_index);
             }
         }
          log_message(4, "mark_nodes_as_deleted_strong Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
     }


    // --- Output Generation Helpers ---

    void PCSTFast::build_node_set_from_edges(const std::vector<int>& edge_set, std::vector<int>* node_set) {
         log_message(3, "build_node_set_from_edges Entry (using %zu edges).\n", edge_set.size());
        if (!node_set) {
             log_message(1, "Warning: node_set pointer is null in build_node_set_from_edges.\n");
             return;
        }
        node_set->clear();
        node_set->reserve(prizes_.size());

        if (build_phase1_included_nodes_.size() != prizes_.size()) {
             log_message(4, "  Resizing build_phase1_included_nodes_ to %zu.\n", prizes_.size());
            build_phase1_included_nodes_.assign(prizes_.size(), 0);
        } else {
             log_message(4, "  Clearing build_phase1_included_nodes_.\n");
            std::fill(build_phase1_included_nodes_.begin(), build_phase1_included_nodes_.end(), 0);
        }

        int nodes_added_from_edges = 0;
        log_message(4, "  Adding nodes from edge endpoints...\n");
        for (int edge_idx : edge_set) {
            if (static_cast<size_t>(edge_idx) >= edges_.size()) {
                 log_message(1, "Warning: Invalid edge index %d in build_node_set.\n", edge_idx);
                 continue;
            }
            const auto& edge = edges_[edge_idx];
            int uu = edge.first;
            int vv = edge.second;

            if (static_cast<size_t>(uu) < build_phase1_included_nodes_.size() && !build_phase1_included_nodes_[uu]) {
                build_phase1_included_nodes_[uu] = 1;
                node_set->push_back(uu);
                nodes_added_from_edges++;
                 log_message(4, "    Added node %d from edge %d.\n", uu, edge_idx);
            }
            if (static_cast<size_t>(vv) < build_phase1_included_nodes_.size() && !build_phase1_included_nodes_[vv]) {
                build_phase1_included_nodes_[vv] = 1;
                node_set->push_back(vv);
                nodes_added_from_edges++;
                 log_message(4, "    Added node %d from edge %d.\n", vv, edge_idx);
            }
        }
         log_message(4, "  Added %d unique nodes from edge endpoints.\n", nodes_added_from_edges);


        int isolated_good_nodes_added = 0;
        log_message(4, "  Adding isolated 'good' nodes...\n");
        for (int ii = 0; ii < std::ssize(prizes_); ++ii) {
            bool is_good = (static_cast<size_t>(ii) < node_good_.size() && node_good_[ii]);
            bool is_included = (static_cast<size_t>(ii) < build_phase1_included_nodes_.size() && build_phase1_included_nodes_[ii]);
            if (is_good && !is_included) {
                 log_message(4, "    Adding isolated good node %d.\n", ii);
                node_set->push_back(ii);
                isolated_good_nodes_added++;
            }
        }
         log_message(4, "  Added %d isolated good nodes.\n", isolated_good_nodes_added);

         log_message(4, "  Sorting final node set (size %zu).\n", node_set->size());
         std::sort(node_set->begin(), node_set->end());
         log_message(3, "build_node_set_from_edges Exit. Final node set size: %zu\n", node_set->size());
    }

    void PCSTFast::build_final_node_set(std::vector<int>* node_set) const {
         log_message(3, "build_final_node_set Entry.\n");
        if (!node_set) {
             log_message(1, "Warning: node_set pointer is null in build_final_node_set.\n");
             return;
        }
        node_set->clear();
        node_set->reserve(prizes_.size());
        int nodes_included = 0;

        log_message(4, "  Iterating through nodes to build final set (based on good & not deleted)...\n");
        for (int ii = 0; ii < std::ssize(prizes_); ++ii) {
            bool is_good = (static_cast<size_t>(ii) < node_good_.size() && node_good_[ii]);
            bool is_deleted = (static_cast<size_t>(ii) < node_deleted_.size() && node_deleted_[ii]);
             log_message(4, "    Node %d: Good=%d, Deleted=%d\n", ii, is_good, is_deleted);

            if (is_good && !is_deleted) {
                 log_message(4, "      Including node %d.\n", ii);
                node_set->push_back(ii);
                nodes_included++;
            }
        }
         log_message(4, "  Sorting final node set (size %zu).\n", node_set->size());
         std::sort(node_set->begin(), node_set->end());
         log_message(3, "build_final_node_set Exit. Final node set size: %zu\n", nodes_included);
    }


    // --- Statistics ---

    void PCSTFast::get_statistics(Statistics* s) const {
        log_message(3, "get_statistics Entry.\n");
        if(s) {
            *s = stats_;
             log_message(3, "  Statistics copied.\n");
        } else {
             log_message(1, "Warning: Null pointer passed to get_statistics.\n");
        }
         log_message(3, "get_statistics Exit.\n");
    }

} // namespace cluster_approx
