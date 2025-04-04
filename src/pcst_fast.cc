#include "pcst_fast.h"
#include "pruning_strategy.h" // Include the new header with Pruner definitions

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
#include <numeric>
#include <optional>
#include <memory> // For unique_ptr from factory
#include <cassert> // For assertions


// Default logger implementation (writes to stderr)
// Defined outside namespace or in an anonymous namespace if preferred
void DefaultStderrLogger( [[maybe_unused]] int level, const std::string& message) {
    fprintf(stderr, "%s", message.c_str());
    fflush(stderr); // Flush stream
}

// No-op logger
// Defined outside namespace or in an anonymous namespace if preferred
void NoOpLogger( [[maybe_unused]] int level, [[maybe_unused]] const std::string& message) {
    // Do nothing
}


namespace cluster_approx {

    // --- Utility Functions ---
    PCSTFast::PruningMethod PCSTFast::parse_pruning_method(std::string_view input) {
        std::string lower_input;
        lower_input.reserve(input.length());
        std::transform(input.begin(), input.end(), std::back_inserter(lower_input),
                       [](unsigned char c){ return std::tolower(c); });
        if (lower_input == "none") return PruningMethod::kNoPruning;
        if (lower_input == "simple") return PruningMethod::kSimplePruning;
        if (lower_input == "gw") return PruningMethod::kGWPruning;
        if (lower_input == "strong") return PruningMethod::kStrongPruning;
        return PruningMethod::kUnknownPruning;
    }

    // Updated logging helper using std::function logger_
    void PCSTFast::format_and_log(int level, const char* format, ...) const {
        // Only format if logging is enabled and level is sufficient
        if (logger_ && verbosity_level_ >= level) {
            char buffer[kOutputBufferSize];
            va_list args;
            va_start(args, format);
            // Careful: vsnprintf is safe against buffer overflows
            int needed = vsnprintf(buffer, kOutputBufferSize, format, args);
            va_end(args);

            if (needed < 0) {
                // Handle vsnprintf error (e.g., encoding issues)
                logger_(level, "[PCST Logging Error: vsnprintf failed]");
            } else if (static_cast<size_t>(needed) >= kOutputBufferSize) {
                // Output was truncated
                buffer[kOutputBufferSize - 1] = '\0'; // Ensure null termination
                // Add ellipsis safely within buffer bounds
                 if (kOutputBufferSize > 4) {
                    snprintf(buffer + kOutputBufferSize - 5, 5, "...");
                 }
                logger_(level, std::string(buffer) + "(TRUNCATED)");
            } else {
                // Output fits completely
                logger_(level, std::string(buffer));
            }
        }
    }


    // --- Constructor and Initialization ---

    PCSTFast::PCSTFast(const std::vector<std::pair<IndexType, IndexType>>& edges,
                       const std::vector<ValueType>& prizes,
                       const std::vector<ValueType>& costs,
                       IndexType root,
                       int target_num_active_clusters,
                       PruningMethod pruning,
                       int verbosity_level,
                       Logger logger
                       )
        : edges_(edges),
          prizes_(prizes),
          costs_(costs),
          root_(root),
          target_num_active_clusters_( (root >= 0 && target_num_active_clusters > 0) ? 0 : target_num_active_clusters),
          // stats_ default initialized
          // edge_parts_, edge_info_, clusters_, inactive_merge_events_ resized later
          clusters_deactivation_(), // Initialized before logger/verbosity
          clusters_next_edge_event_(),
          // current_time_ initialized inline to 0.0
          // pruner_ initialized below
          logger_(logger ? logger : NoOpLogger), // Use provided logger or NoOpLogger
          verbosity_level_(verbosity_level)     // Initialized last
          // node_good_, phase1_result_ resized later
    {
        format_and_log(3, "PCSTFast Constructor Entry.\n");
        const size_t num_nodes = prizes_.size();
        const size_t num_edges = edges_.size();
        format_and_log(3, "Input sizes: Nodes=%zu, Edges=%zu\n", num_nodes, num_edges);
        format_and_log(3, "Parameters: Root=%d, TargetClusters=%d, Pruning=%d, Verbosity=%d\n",
                    root_, target_num_active_clusters_, static_cast<int>(pruning), verbosity_level_);

        // --- Input Validation ---
        if (root_ >= static_cast<IndexType>(num_nodes)) {
            format_and_log(0, "Error: Root node index %d out of range (Num nodes: %zu).\n", root_, num_nodes);
            throw std::invalid_argument("Root node index out of range.");
        }
        if (target_num_active_clusters_ < 0) {
             format_and_log(0, "Error: Target number of active clusters %d cannot be negative.\n", target_num_active_clusters_);
            throw std::invalid_argument("Target number of active clusters cannot be negative.");
        }
        if (root_ >= 0 && target_num_active_clusters > 0) { // Check original parameter
            format_and_log(1, "Warning: target_num_active_clusters > 0 is ignored in the rooted case (using 0 internally).\n");
        }

        // --- Instantiate Pruner ---
        try {
            pruner_ = internal::create_pruner(pruning); // Use factory function
             format_and_log(3, "Instantiated %s Pruner.\n",
                           pruning == PruningMethod::kNoPruning ? "No" :
                           pruning == PruningMethod::kSimplePruning ? "Simple" :
                           pruning == PruningMethod::kGWPruning ? "GW" :
                           pruning == PruningMethod::kStrongPruning ? "Strong" : "Unknown");
        } catch (const std::invalid_argument& e) {
             format_and_log(0, "Error creating pruner: %s\n", e.what());
             // Allow kUnknownPruning? No, the factory should handle it.
             throw; // Re-throw exception from factory
        }

        // --- Preallocate / Resize Member Vectors ---
        format_and_log(4, "Resizing internal vectors...\n");
        edge_parts_.resize(2 * num_edges);
        edge_info_.resize(num_edges);
        node_good_.resize(num_nodes, 0);

        clusters_.reserve(2 * num_nodes); // Estimate: nodes + potential merge clusters
        inactive_merge_events_.reserve(num_edges); // Upper bound
        phase1_result_.reserve(num_edges); // Upper bound

        // --- Initialize State ---
        format_and_log(4, "Initializing algorithm state...\n");
        current_time_ = 0.0;
        stats_.reset(); // Ensure stats start clean

        initialize_clusters();
        initialize_edges();

        format_and_log(4, "Initializing cluster_next_edge_event PQ...\n");
        // Initialize the cluster next edge event priority queue
        for (IndexType ii = 0; ii < static_cast<IndexType>(clusters_.size()) && ii < static_cast<IndexType>(prizes_.size()) ; ++ii) {
            if (clusters_[ii].active && !clusters_[ii].edge_parts.is_empty()) {
                ValueType val;
                PayloadType edge_part;
                // Check return value of get_min
                if (clusters_[ii].edge_parts.get_min(&val, &edge_part)) {
                    format_and_log(4, "  Cluster %d initial min edge event at time %.9g (part %d).\n", ii, val, edge_part);
                    clusters_next_edge_event_.insert(val, ii);
                } else {
                     format_and_log(0, "Error: Failed to get min edge for active cluster %d during init, though heap reported not empty.\n", ii);
                     // This indicates an internal error in PairingHeap or state management
                     assert(false && "Pairing heap get_min failed after is_empty check");
                }
            }
        }
        format_and_log(3, "PCSTFast Constructor Exit.\n");
    }

    // Define destructor here where internal::IPruner is defined (via includes)
    PCSTFast::~PCSTFast() {
         format_and_log(3, "PCSTFast Destructor.\n");
         // unique_ptr pruner_ will be automatically deleted.
    }

    void PCSTFast::initialize_clusters() {
        format_and_log(3, "initialize_clusters Entry.\n");
        const size_t num_nodes = prizes_.size();
        clusters_.clear(); // Ensure clean state if constructor were reused

        for (size_t ii = 0; ii < num_nodes; ++ii) {
            IndexType current_node_index = static_cast<IndexType>(ii); // Use alias
            format_and_log(4, "Initializing cluster %d (Node %d).\n", current_node_index, current_node_index);
            if (prizes_[ii] < 0.0) {
                 format_and_log(0, "Error: Prize for node %d is negative (%.9g).\n", current_node_index, prizes_[ii]);
                throw std::invalid_argument("Node prize cannot be negative.");
            }

            clusters_.emplace_back(); // Create cluster using default constructor
            Cluster& current_cluster = clusters_.back();

            current_cluster.active = (current_node_index != root_);
            current_cluster.active_start_time = 0.0;
            current_cluster.active_end_time = (current_cluster.active) ? -1.0 : 0.0; // Inactive root starts/ends at t=0
            current_cluster.prize_sum = prizes_[ii];
            current_cluster.contains_root = (current_node_index == root_);
            current_cluster.moat = 0.0; // Initialize moat
            // Other members default/zero initialized

             format_and_log(4, "  Node %d: Prize=%.9g, Active=%d, ContainsRoot=%d\n",
                         current_node_index, current_cluster.prize_sum, current_cluster.active, current_cluster.contains_root);


            if (current_cluster.active) {
                // Insert into deactivation queue only if prize > epsilon (deactivation time > 0)
                if (prizes_[ii] > kEpsilon) {
                    format_and_log(4, "  Node %d: Scheduling deactivation at time %.9g.\n", current_node_index, prizes_[ii]);
                    clusters_deactivation_.insert(prizes_[ii], current_node_index); // Event time = prize
                } else {
                    // If prize is zero or negligible, deactivate immediately
                    format_and_log(4, "  Node %d: Prize is negligible, immediately inactive.\n", current_node_index);
                     current_cluster.active = false;
                     current_cluster.active_end_time = 0.0;
                }
            }
        }
        format_and_log(3, "initialize_clusters Exit.\n");
    }

    void PCSTFast::initialize_edges() {
        format_and_log(3, "initialize_edges Entry.\n");
        const size_t num_edges = edges_.size();
        const size_t num_nodes = prizes_.size();

        for (size_t ii = 0; ii < num_edges; ++ii) {
            const auto& edge = edges_[ii];
            IndexType uu = edge.first; // Use alias
            IndexType vv = edge.second; // Use alias
            ValueType cost = costs_[ii]; // Use alias
            IndexType edge_index = static_cast<IndexType>(ii); // Use alias if needed later

            format_and_log(4, "Initializing edge %d: (%d <-> %d), Cost=%.9g\n", edge_index, uu, vv, cost);

            // --- Edge Validation ---
            if (uu < 0 || vv < 0 || static_cast<size_t>(uu) >= num_nodes || static_cast<size_t>(vv) >= num_nodes) {
                format_and_log(0, "Error: Edge %d endpoint index (%d or %d) out of range (Num nodes: %zu).\n", edge_index, uu, vv, num_nodes);
                throw std::invalid_argument("Edge endpoint index out of range.");
            }
            if (cost < 0.0) {
                format_and_log(0, "Error: Edge %d cost is negative (%.9g).\n", edge_index, cost);
                throw std::invalid_argument("Edge cost cannot be negative.");
            }
            if (uu == vv) {
                format_and_log(1, "Warning: Ignoring self-loop edge %d (%d -> %d).\n", edge_index, uu, vv);
                // Mark edge parts as deleted? Or just skip initialization. Skipping is simpler.
                edge_parts_[2 * ii].deleted = true;
                edge_parts_[2 * ii + 1].deleted = true;
                continue;
            }

            // --- Initialize Edge Parts ---
            assert(static_cast<size_t>(2 * ii + 1) < edge_parts_.size() && "Edge parts index out of bounds");
            EdgePart& uu_part = edge_parts_[2 * ii];
            EdgePart& vv_part = edge_parts_[2 * ii + 1];
            // Initial clusters map 1:1 with nodes
            assert(static_cast<size_t>(uu) < clusters_.size() && static_cast<size_t>(vv) < clusters_.size() && "Cluster index out of bounds");
            Cluster& uu_cluster = clusters_[uu];
            Cluster& vv_cluster = clusters_[vv];

            uu_part.deleted = false;
            vv_part.deleted = false;

            // Determine initial event times based on cluster activity
            ValueType event_time_u, event_time_v;
            if (uu_cluster.active && vv_cluster.active) {
                event_time_u = cost * 0.5;
                event_time_v = cost * 0.5;
                 format_and_log(4, "  Edge %d: Active-Active. Event times = %.9g\n", edge_index, event_time_u);
            } else if (uu_cluster.active) { // Only u is active
                event_time_u = cost;
                event_time_v = 0.0; // Inactive side processes immediately (effectively)
                 format_and_log(4, "  Edge %d: Active(U)-Inactive(V). Event time U=%.9g, V=%.9g\n", edge_index, event_time_u, event_time_v);
            } else if (vv_cluster.active) { // Only v is active
                event_time_u = 0.0;
                event_time_v = cost;
                 format_and_log(4, "  Edge %d: Inactive(U)-Active(V). Event time U=%.9g, V=%.9g\n", edge_index, event_time_u, event_time_v);
            } else { // Both inactive (e.g., root connected to inactive)
                event_time_u = 0.0;
                event_time_v = 0.0;
                 format_and_log(4, "  Edge %d: Inactive-Inactive. Event times = 0.0\n", edge_index);
            }

            // Use PayloadType consistently for heap payload
            PayloadType part_payload_u = static_cast<PayloadType>(2 * ii);
            PayloadType part_payload_v = static_cast<PayloadType>(2 * ii + 1);

            // Insert edge parts into their respective cluster pairing heaps only if edge can participate
            uu_part.next_event_val = event_time_u;
            if (uu_cluster.active || vv_cluster.active) {
                 format_and_log(4, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", part_payload_u, uu, uu, event_time_u);
                 uu_part.heap_node = uu_cluster.edge_parts.insert(uu_part.next_event_val, part_payload_u);
            } else {
                 format_and_log(4, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", part_payload_u, uu);
            }

            vv_part.next_event_val = event_time_v;
             if (uu_cluster.active || vv_cluster.active) {
                 format_and_log(4, "  Inserting edge part %d (Node %d) into cluster %d heap with value %.9g.\n", part_payload_v, vv, vv, event_time_v);
                 vv_part.heap_node = vv_cluster.edge_parts.insert(vv_part.next_event_val, part_payload_v);
             } else {
                  format_and_log(4, "  Skipping heap insert for edge part %d (Node %d) as both clusters inactive.\n", part_payload_v, vv);
             }
        }
        format_and_log(3, "initialize_edges Exit.\n");
    }


    // --- Main Algorithm Loop ---

    bool PCSTFast::run(std::vector<IndexType>* result_nodes, std::vector<IndexType>* result_edges) {
        format_and_log(3, "PCSTFast::run Entry.\n");
        if (!result_nodes || !result_edges) {
             format_and_log(0, "Error: Result pointers cannot be null in run().\n");
            throw std::invalid_argument("Result pointers cannot be null.");
        }

        result_nodes->clear();
        result_edges->clear();
        phase1_result_.clear(); // Clear phase 1 result from previous runs if any
        stats_.reset(); // Ensure stats are clean for this run

        int current_num_active_clusters = 0;
        for(const auto& cluster : clusters_) {
            // Count only active *representative* clusters
            if (cluster.active && cluster.merged_into == kInvalidIndex) {
                current_num_active_clusters++;
            }
        }
        format_and_log(1, "Starting GW phase. Initial active clusters: %d (Target: %d)\n",
                    current_num_active_clusters, target_num_active_clusters_);


        // --- Main Event Loop ---
        int event_count = 0;
        while (current_num_active_clusters > target_num_active_clusters_) {
            event_count++;
            format_and_log(3, "--- GW Loop Iteration %d (Active Clusters: %d) ---\n", event_count, current_num_active_clusters);

            NextEvent event = get_next_event();
             format_and_log(3, "get_next_event returned: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                        static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);


            if (event.type == NextEvent::Type::NONE) {
                format_and_log(1, "No more events possible. Terminating GW loop.\n");
                break;
            }

            format_and_log(2, "-----------------------------------------\n");
             // Check for time going backwards (indicates a serious bug)
             if (event.time < current_time_ - kEpsilon) {
                  format_and_log(0, "Error: Time decreased! Current=%.9g, Event=%.9g\n", current_time_, event.time);
                  assert(false && "Time decreased in main loop"); // Assert in debug builds
                  break; // Stop processing if time decreases
             }
            current_time_ = event.time; // Advance time
             format_and_log(2, "Advanced time to %.9g\n", current_time_);

            process_event(event);

            // Recalculate active cluster count after processing the event
            int previous_active_count = current_num_active_clusters;
            current_num_active_clusters = 0;
            for(const auto& cluster : clusters_) {
                // Count only active *representative* clusters
                if (cluster.active && cluster.merged_into == kInvalidIndex) {
                    current_num_active_clusters++;
                }
            }
             format_and_log(3, "End of event processing. Active clusters: %d (was %d)\n", current_num_active_clusters, previous_active_count);

        } // End while loop

        format_and_log(1, "Finished GW main loop: Final time %.9g, Total edge events processed %lld, Final active clusters %d\n",
                    current_time_, stats_.total_num_edge_events, current_num_active_clusters);

        // --- Pruning and Output Generation ---
         format_and_log(2, "--- Starting Pruning and Output Generation ---\n");

        // Phase 1: Determine 'good' nodes based on connectivity to final active/root clusters
         format_and_log(2, "Selecting initial good nodes...\n");
        node_good_.assign(prizes_.size(), 0); // Reset flags
        select_initial_active_clusters(node_good_);
         int good_node_count = 0;
         for(uint8_t flag : node_good_) good_node_count += flag;
         format_and_log(2, "Found %d good nodes initially.\n", good_node_count);

        format_and_log(3, "Creating Pruning Context...\n");
        // Create the context to pass state to the pruner
        internal::PruningContext context(
            edges_, prizes_, costs_, root_, clusters_, edge_info_,
            inactive_merge_events_, node_good_, phase1_result_,
            logger_, verbosity_level_
        );

        format_and_log(2, "Calling pruner strategy prune method...\n");
        assert(pruner_ && "Pruner should have been initialized in constructor");
        // Delegate the actual pruning and result finalization to the strategy object
        pruner_->prune(context, *result_nodes, *result_edges);

        format_and_log(3, "PCSTFast::run Exit (Success).\n");
        return true;
    }

    // --- Event Handling ---
    PCSTFast::NextEvent PCSTFast::get_next_event() const {
        format_and_log(4, "get_next_event Entry.\n");
        NextEvent event; // Default type is NONE, time is infinity

        // --- Declare local variables using original naming conventions ---
        ValueType next_edge_time = std::numeric_limits<ValueType>::infinity();
        IndexType next_edge_cluster_index = kInvalidIndex;
        PayloadType next_edge_part_index = kInvalidIndex;
        bool edge_event_found = false;

        ValueType next_cluster_time = std::numeric_limits<ValueType>::infinity();
        IndexType next_cluster_deactivation_index = kInvalidIndex;
        bool cluster_event_found = false;

        // --- Temporary variables needed for retrieving values from PQs ---
        ValueType temp_pq_edge_time;
        IndexType temp_pq_edge_cluster;
        ValueType temp_cluster_min_time; // Renamed to avoid conflict
        PayloadType temp_cluster_min_part; // Renamed to avoid conflict

        ValueType temp_pq_cluster_time;
        IndexType temp_pq_cluster_index;

        // --- Get potential next edge event ---
        format_and_log(4, "Checking clusters_next_edge_event PQ...\n");
        // Check return value of the *global* PQ get_min
        if (clusters_next_edge_event_.get_min(&temp_pq_edge_time, &temp_pq_edge_cluster)) {
            format_and_log(4, "  Global Edge PQ min: Time=%.9g, Cluster=%d\n", temp_pq_edge_time, temp_pq_edge_cluster);
            // Verify the cluster index and get the *actual* min from the cluster's heap
            if (temp_pq_edge_cluster != kInvalidIndex && static_cast<size_t>(temp_pq_edge_cluster) < clusters_.size()) {
                format_and_log(4, "  Checking cluster %d's internal heap...\n", temp_pq_edge_cluster);
                // Check return value of the *cluster's* heap get_min
                if (clusters_[temp_pq_edge_cluster].edge_parts.get_min(&temp_cluster_min_time, &temp_cluster_min_part)) {
                     format_and_log(4, "  Cluster %d internal min: Time=%.9g, Part=%d\n", temp_pq_edge_cluster, temp_cluster_min_time, temp_cluster_min_part);
                     // Valid event found - store details in the main variables
                     next_edge_time = temp_cluster_min_time;
                     next_edge_cluster_index = temp_pq_edge_cluster;
                     next_edge_part_index = temp_cluster_min_part;
                     edge_event_found = true;
                } else {
                     // Inconsistency: Global PQ claimed cluster had edges, but its heap was empty. Log this.
                     format_and_log(0, "Error: Inconsistency in get_next_event. Global edge PQ has cluster %d, but its heap is empty.\n", temp_pq_edge_cluster);
                     // This indicates a potential bug elsewhere where the global PQ wasn't updated correctly.
                     assert(false && "Global Edge PQ inconsistent with cluster heap");
                }
            } else {
                 // Invalid cluster index retrieved from global PQ. Log this.
                 format_and_log(0, "Error: Invalid cluster index %d retrieved from global edge PQ in get_next_event.\n", temp_pq_edge_cluster);
            }
        } else {
             format_and_log(4, "  Global Edge PQ is empty.\n");
        }


        // --- Get potential next cluster deactivation event ---
        format_and_log(4, "Checking clusters_deactivation PQ...\n");
        // Check return value of the deactivation PQ get_min
        if (clusters_deactivation_.get_min(&temp_pq_cluster_time, &temp_pq_cluster_index)) {
            format_and_log(4, "  Deactivation PQ min: Time=%.9g, Cluster=%d\n", temp_pq_cluster_time, temp_pq_cluster_index);
             if (temp_pq_cluster_index != kInvalidIndex) {
                 // Valid event found - store details in the main variables
                 next_cluster_time = temp_pq_cluster_time;
                 next_cluster_deactivation_index = temp_pq_cluster_index;
                 cluster_event_found = true;
             } else {
                 // Log if an invalid index was somehow stored and retrieved.
                 format_and_log(0, "Error: Invalid cluster index %d retrieved from deactivation PQ in get_next_event.\n", temp_pq_cluster_index);
                 // Do not mark cluster_event_found = true if index is invalid.
             }
        } else {
             format_and_log(4, "  Deactivation PQ is empty.\n");
        }


        // --- Determine the actual next event based on existence and time ---
        format_and_log(4, "Determining next event: EdgeFound=%d (Time=%.9g), ClusterFound=%d (Time=%.9g)\n",
                    edge_event_found, next_edge_time, cluster_event_found, next_cluster_time);
        if (edge_event_found && cluster_event_found) {
            // Both types of events exist, compare times
            // Use the local variables populated above
            if (next_edge_time <= next_cluster_time + kEpsilon * std::max(1.0, std::abs(next_cluster_time))) {
                // Edge event is next (or ties go to edge)
                format_and_log(4,"  Edge event is next (or tie).\n");
                event.type = NextEvent::Type::EDGE;
                event.time = next_edge_time;
                event.cluster_index = next_edge_cluster_index; // Use declared local variable
                event.edge_part_index = next_edge_part_index; // Use declared local variable
            } else {
                // Cluster event is next
                 format_and_log(4,"  Cluster event is next.\n");
                event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
                event.time = next_cluster_time;
                event.cluster_index = next_cluster_deactivation_index; // Use declared local variable
                // event.edge_part_index remains kInvalidIndex
            }
        } else if (edge_event_found) {
            // Only edge event exists
             format_and_log(4,"  Only Edge event found.\n");
            event.type = NextEvent::Type::EDGE;
            event.time = next_edge_time;
            event.cluster_index = next_edge_cluster_index; // Use declared local variable
            event.edge_part_index = next_edge_part_index; // Use declared local variable
        } else if (cluster_event_found) {
            // Only cluster event exists
             format_and_log(4,"  Only Cluster event found.\n");
            event.type = NextEvent::Type::CLUSTER_DEACTIVATION;
            event.time = next_cluster_time;
            event.cluster_index = next_cluster_deactivation_index; // Use declared local variable
            // event.edge_part_index remains kInvalidIndex
        } else {
            // Neither event exists, event remains Type::NONE with time=infinity
             format_and_log(4,"  No events found.\n");
        }

        format_and_log(4, "get_next_event Exit. Returning Type=%d, Time=%.9g\n", static_cast<int>(event.type), event.time);
        return event;
    }

    void PCSTFast::process_event(const NextEvent& event) {
        format_and_log(3, "process_event Entry: Type=%d, Time=%.9g, Cluster=%d, EdgePart=%d\n",
                    static_cast<int>(event.type), event.time, event.cluster_index, event.edge_part_index);

        if (event.type == NextEvent::Type::EDGE) {
            stats_.total_num_edge_events++;
            format_and_log(4, "Processing EDGE event for Cluster %d, EdgePart %d\n", event.cluster_index, event.edge_part_index);

            // --- Pre-checks ---
            if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size()) {
                 format_and_log(1, "Warning: Edge event cluster index %d invalid. Skipping.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index); // Attempt cleanup
                 return;
            }
             if (!clusters_[event.cluster_index].active) {
                 format_and_log(2, "Edge event for inactive cluster %d. Removing from global PQ. Skipping logic.\n", event.cluster_index);
                 clusters_next_edge_event_.delete_element(event.cluster_index);
                 // No need to delete from cluster heap, it shouldn't have been in global PQ if inactive.
                 return;
            }
            if (event.edge_part_index < 0 || static_cast<size_t>(event.edge_part_index) >= edge_parts_.size()) {
                 format_and_log(0, "Error: Invalid edge part index %d in edge event. Skipping logic.\n", event.edge_part_index);
                 assert(false && "Invalid edge part index in event processing");
                 return; // Cannot proceed
            }

            // --- Delete event from originating cluster's heap ---
             format_and_log(4, "Deleting min edge part from cluster %d's heap...\n", event.cluster_index);
             ValueType deleted_val; PayloadType deleted_part;
             // Use return value of delete_min
             if (!clusters_[event.cluster_index].edge_parts.delete_min(&deleted_val, &deleted_part)) {
                  format_and_log(0, "Error: Failed to delete min edge from cluster %d's heap (event time %.9g).\n", event.cluster_index, event.time);
                  clusters_next_edge_event_.delete_element(event.cluster_index); // Clean up global PQ
                  // If delete_min failed, the cluster's heap might be empty, despite get_min succeeding earlier.
                  assert(clusters_[event.cluster_index].edge_parts.is_empty() && "delete_min failed but heap not empty");
                  return;
             }
              format_and_log(4, "  Deleted Part=%d, Value=%.9g. Expected Part=%d, Value=%.9g.\n", deleted_part, deleted_val, event.edge_part_index, event.time);
             // Verify consistency (allow small float differences)
             if (deleted_part != event.edge_part_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                  format_and_log(0, "Error: Mismatch in edge event deletion for cluster %d! Deleted {%d, %.9g}, Expected {%d, %.9g}\n",
                                 event.cluster_index, deleted_part, deleted_val, event.edge_part_index, event.time);
                   // This suggests the global PQ might be out of sync with the cluster heap.
                   assert(false && "Mismatch during edge event deletion");
                   // Continue processing, but state might be corrupt.
             }

            // --- Update global PQ for the cluster ---
             format_and_log(4, "Updating global edge PQ for cluster %d after deletion.\n", event.cluster_index);
             update_cluster_queues_post_growth(event.cluster_index); // Re-inserts cluster with new min edge time, if any

            // --- Check if specific edge part was logically deleted ---
            if (edge_parts_[event.edge_part_index].deleted) {
                stats_.num_deleted_edge_events++;
                format_and_log(2, "Edge part %d already deleted (flag=true), skipping event logic.\n", event.edge_part_index);
                return; // Already processed/deleted, skip
            }

            // --- Handle the core logic ---
            format_and_log(4, "Calling handle_edge_event for EdgePart %d.\n", event.edge_part_index);
            handle_edge_event(event.edge_part_index);

        } else if (event.type == NextEvent::Type::CLUSTER_DEACTIVATION) {
            stats_.num_cluster_events++;
             format_and_log(4, "Processing CLUSTER_DEACTIVATION event for Cluster %d.\n", event.cluster_index);

            // --- Remove the event from the deactivation PQ ---
             format_and_log(4, "Deleting min from cluster deactivation PQ...\n");
             ValueType deleted_val; IndexType deleted_cluster;
             // Use return value of delete_min
             if(!clusters_deactivation_.delete_min(&deleted_val, &deleted_cluster)) {
                 format_and_log(1, "Warning: Failed to delete cluster deactivation event from PQ (expected cluster %d). Maybe already merged?\n", event.cluster_index);
                 // Event might have been removed by a merge earlier.
                 // We only need to proceed if the cluster is actually still active.
                 if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size() || !clusters_[event.cluster_index].active) {
                      format_and_log(2, "Skipping cluster deactivation logic as cluster %d is invalid or already inactive.\n", event.cluster_index);
                      return; // Don't process if already inactive
                 }
                 // If it *is* still active despite PQ delete failing, log inconsistency and proceed.
                  format_and_log(0, "Error: PQ delete failed but cluster %d is still active!\n", event.cluster_index);
                  assert(false && "Deactivation PQ delete failed but cluster still active");
             } else {
                 // Successfully deleted *something* from PQ. Check if it matches.
                 format_and_log(4, "  Deleted Cluster=%d, Value=%.9g. Expected Cluster=%d, Value=%.9g.\n", deleted_cluster, deleted_val, event.cluster_index, event.time);
                 // Check consistency
                 if(deleted_cluster != event.cluster_index || std::abs(deleted_val - event.time) > kEpsilon * std::max(1.0, std::abs(event.time))) {
                      format_and_log(1, "Warning: Mismatch in cluster deactivation event deletion (possibly due to prior merge). Proceeding with target cluster %d.\n", event.cluster_index);
                      // Ensure we are processing the *intended* cluster, even if PQ was slightly off.
                 }
                 // Check if the *intended* cluster is actually still active
                  if (event.cluster_index < 0 || static_cast<size_t>(event.cluster_index) >= clusters_.size() || !clusters_[event.cluster_index].active) {
                       format_and_log(2, "Skipping cluster deactivation logic as target cluster %d is invalid or already inactive.\n", event.cluster_index);
                       return;
                  }
             }

            // --- Handle the core logic ---
            format_and_log(4, "Calling handle_cluster_deactivation_event for Cluster %d.\n", event.cluster_index);
            handle_cluster_deactivation_event(event.cluster_index); // Process the intended cluster
        } else {
            format_and_log(0, "Error: Processed an event of type NONE.\n");
             assert(false && "Processed event of type NONE");
        }
         format_and_log(3, "process_event Exit.\n");
    }

    // --- Edge Event Sub-logic Implementation ---
    void PCSTFast::handle_edge_event(PayloadType edge_part_index) {
         format_and_log(3, "handle_edge_event Entry: EdgePart=%d\n", edge_part_index);
         format_and_log(4, "Getting edge processing info for part %d...\n", edge_part_index);
         std::optional<EdgeProcessingInfo> opt_info = get_edge_processing_info(edge_part_index);
         if (!opt_info) {
             format_and_log(1, "Warning: Could not get valid edge processing info for part %d. Skipping.\n", edge_part_index);
             return; // Cannot proceed if info is invalid
         }
         EdgeProcessingInfo info = *opt_info; // Copy the info

         // Log detailed info obtained
         format_and_log(2, "Edge Event Details: Edge %d (%d <-> %d), Cost %.9g\n",
                     info.edge_idx, edges_[info.edge_idx].first, edges_[info.edge_idx].second, info.cost);
         format_and_log(2, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                     info.current_part_idx, info.current_cluster_idx, info.sum_current, info.finished_moat_current);
         format_and_log(2, "  Part %d: Cluster %d, Sum %.9g, FinishedMoat %.9g\n",
                     info.other_part_idx, info.other_cluster_idx, info.sum_other, info.finished_moat_other);
         format_and_log(2, "  Remainder: %.9g\n", info.remainder);

         // Basic assertions for sanity
         assert(info.current_cluster_idx != kInvalidIndex && static_cast<size_t>(info.current_cluster_idx) < clusters_.size() && "Invalid current cluster index in handle_edge_event");
         assert(info.other_cluster_idx != kInvalidIndex && static_cast<size_t>(info.other_cluster_idx) < clusters_.size() && "Invalid other cluster index in handle_edge_event");
         assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size() && "Invalid current part index");
         assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size() && "Invalid other part index");


         // Check if clusters are already merged (representatives are the same)
         if (info.current_cluster_idx == info.other_cluster_idx) {
             stats_.num_merged_edge_events++;
             format_and_log(2, "Clusters %d already merged. Marking other edge part %d as deleted.\n",
                         info.current_cluster_idx, info.other_part_idx);
             edge_parts_[info.other_part_idx].deleted = true; // Mark the *other* side as deleted
             // Note: The current side was already removed from its heap in process_event.
             // The other side might still be in its heap, but this 'deleted' flag prevents processing it later.
             format_and_log(3, "handle_edge_event Exit (Clusters already merged).\n");
             return;
         }

         // Check merge condition: remainder is negligible compared to cost or absolutely negligible
         bool merge_condition = (info.remainder <= kEpsilon * std::max(1.0, info.cost) || info.remainder <= kEpsilon );
         format_and_log(4, "Merge condition check: Remainder=%.9g, Epsilon*Max(1,Cost)=%.9g, Epsilon=%.9g -> Merge=%d\n",
                     info.remainder, kEpsilon * std::max(1.0, info.cost), kEpsilon, merge_condition);

         if (merge_condition) {
             // --- Merge Clusters ---
             stats_.total_num_merge_events++;
             format_and_log(2, "Merge condition met for edge %d.\n", info.edge_idx);
             format_and_log(3, "Adding edge %d to phase1_result.\n", info.edge_idx);
             phase1_result_.push_back(info.edge_idx); // Record the merging edge
             format_and_log(3, "Marking other edge part %d as deleted.\n", info.other_part_idx);
             edge_parts_[info.other_part_idx].deleted = true; // Mark other side logically deleted

             format_and_log(4, "Calling merge_clusters...\n");
             merge_clusters(info); // Perform the merge operation

         } else {
             // --- Edge Growth ---
             stats_.total_num_edge_growth_events++;
             format_and_log(2, "Growth condition for edge %d.\n", info.edge_idx);

             assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size() && "Invalid other cluster index before growth check");
             if (clusters_[info.other_cluster_idx].active) {
                 // Active-Active Growth
                 format_and_log(4, "Calling handle_active_active_growth...\n");
                 handle_active_active_growth(info);
             } else {
                 // Active-Inactive Growth
                 format_and_log(4, "Calling handle_active_inactive_growth...\n");
                 handle_active_inactive_growth(info);
             }
         }
         format_and_log(3, "handle_edge_event Exit.\n");
    }

    std::optional<PCSTFast::EdgeProcessingInfo> PCSTFast::get_edge_processing_info(PayloadType edge_part_index) {
        format_and_log(4, "get_edge_processing_info Entry: EdgePart=%d\n", edge_part_index);
        EdgeProcessingInfo info;
        info.current_part_idx = edge_part_index;
        // Calculate related indices
        info.other_part_idx = get_other_edge_part_index(edge_part_index);
        info.edge_idx = edge_part_index / 2;
        format_and_log(4, "  Derived OtherPart=%d, EdgeIdx=%d\n", info.other_part_idx, info.edge_idx);

        // Validate calculated indices
        if (static_cast<size_t>(info.edge_idx) >= costs_.size() ||
            static_cast<size_t>(info.other_part_idx) >= edge_parts_.size())
        {
             format_and_log(0, "Error: Invalid edge index %d or part index %d derived from %d.\n",
                         info.edge_idx, info.other_part_idx, edge_part_index);
              assert(false && "Invalid index derived in get_edge_processing_info");
            return std::nullopt; // Return empty optional on failure
        }
        info.cost = costs_[info.edge_idx];
        format_and_log(4, "  EdgeCost=%.9g\n", info.cost);

        // Get sum and representative for the current side
        format_and_log(4, "  Getting sum for current part %d...\n", info.current_part_idx);
        get_sum_on_edge_part(info.current_part_idx, &info.sum_current, &info.finished_moat_current, &info.current_cluster_idx);
        format_and_log(4, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.current_cluster_idx, info.sum_current, info.finished_moat_current);

        // Get sum and representative for the other side
        format_and_log(4, "  Getting sum for other part %d...\n", info.other_part_idx);
        get_sum_on_edge_part(info.other_part_idx, &info.sum_other, &info.finished_moat_other, &info.other_cluster_idx);
        format_and_log(4, "  Result: Cluster=%d, Sum=%.9g, FinishedMoat=%.9g\n", info.other_cluster_idx, info.sum_other, info.finished_moat_other);

        // Validate the representatives found
        if (info.current_cluster_idx == kInvalidIndex || info.other_cluster_idx == kInvalidIndex ||
            static_cast<size_t>(info.current_cluster_idx) >= clusters_.size() ||
            static_cast<size_t>(info.other_cluster_idx) >= clusters_.size())
        {
            format_and_log(1, "Warning: Failed to find valid representatives in get_edge_processing_info.\n");
            return std::nullopt; // Return empty optional on failure
        }

        // Calculate remaining cost
        info.remainder = info.cost - info.sum_current - info.sum_other;
        format_and_log(4, "  Calculated Remainder = %.9g - %.9g - %.9g = %.9g\n", info.cost, info.sum_current, info.sum_other, info.remainder);

        format_and_log(4, "get_edge_processing_info Exit (Success).\n");
        return info; // Return the populated struct
     }

    void PCSTFast::merge_clusters(const EdgeProcessingInfo& info) {
        format_and_log(3, "merge_clusters Entry: Merging Cluster %d and Cluster %d via Edge %d.\n",
                    info.current_cluster_idx, info.other_cluster_idx, info.edge_idx);
        IndexType cluster1_idx = info.current_cluster_idx;
        IndexType cluster2_idx = info.other_cluster_idx;
        // Assertions to catch errors early in debug builds
        assert(static_cast<size_t>(cluster1_idx) < clusters_.size() && "Invalid cluster1 index in merge");
        assert(static_cast<size_t>(cluster2_idx) < clusters_.size() && "Invalid cluster2 index in merge");
        Cluster& cluster1 = clusters_[cluster1_idx]; // Get references
        Cluster& cluster2 = clusters_[cluster2_idx];

        // --- Create New Merged Cluster ---
        IndexType new_cluster_idx = clusters_.size(); // Index will be the current size
        format_and_log(4, "Creating new cluster at index %d.\n", new_cluster_idx);
        clusters_.emplace_back(); // Add new cluster element
        Cluster& new_cluster = clusters_.back(); // Get reference to the new cluster

        // --- Update New Cluster State ---
        new_cluster.prize_sum = cluster1.prize_sum + cluster2.prize_sum;
        // Base sum only includes *previously finished* moats
        new_cluster.subcluster_moat_sum = cluster1.subcluster_moat_sum + cluster2.subcluster_moat_sum;
        new_cluster.contains_root = cluster1.contains_root || cluster2.contains_root;
        new_cluster.active = !new_cluster.contains_root; // Active only if root not present
        new_cluster.merged_along = info.edge_idx; // Record merging edge
        new_cluster.child_cluster_1 = cluster1_idx; // Record children (representatives at merge time)
        new_cluster.child_cluster_2 = cluster2_idx;
        format_and_log(4, "  New cluster %d: PrizeSum=%.9g, BaseMoatSum=%.9g, ContainsRoot=%d, Active=%d\n",
                     new_cluster_idx, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, new_cluster.contains_root, new_cluster.active);

        // --- Update Child Cluster 1 State ---
        format_and_log(4, "Updating state for child cluster %d (was Active=%d).\n", cluster1_idx, cluster1.active);
        cluster1.merged_into = new_cluster_idx; // Point child to new parent
        if (cluster1.active) {
            // If child was active, finalize its state
            cluster1.active = false;
            cluster1.active_end_time = current_time_;
            cluster1.moat = cluster1.active_end_time - cluster1.active_start_time;
            // Clamp moat if time precision issues occur
            if (cluster1.moat < 0.0) { format_and_log(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster1_idx, cluster1.moat); cluster1.moat = 0.0; }
            format_and_log(4, "  Child %d became inactive. Moat=%.9g\n", cluster1_idx, cluster1.moat);
            new_cluster.subcluster_moat_sum += cluster1.moat; // Add finalized moat to parent's sum
        } else {
            // If child was already inactive, just add its previously calculated final moat
            format_and_log(4, "  Child %d already inactive. Adding previous moat %.9g.\n", cluster1_idx, cluster1.moat);
            new_cluster.subcluster_moat_sum += cluster1.moat;
        }

        // --- Update Child Cluster 2 State ---
        format_and_log(4, "Updating state for child cluster %d (was Active=%d).\n", cluster2_idx, cluster2.active);
        cluster2.merged_into = new_cluster_idx;
        if (cluster2.active) {
            // Active-Active merge case
            stats_.num_active_active_merge_events++;
            cluster2.active = false;
            cluster2.active_end_time = current_time_;
            cluster2.moat = cluster2.active_end_time - cluster2.active_start_time;
            if (cluster2.moat < 0.0) { format_and_log(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster2_idx, cluster2.moat); cluster2.moat = 0.0; }
            format_and_log(4, "  Child %d (Active) became inactive. Moat=%.9g (Active-Active Merge)\n", cluster2_idx, cluster2.moat);
            new_cluster.subcluster_moat_sum += cluster2.moat;
        } else {
            // Active-Inactive merge case
            stats_.num_active_inactive_merge_events++;
            format_and_log(4, "  Child %d already inactive. Adding previous moat %.9g (Active-Inactive Merge)\n", cluster2_idx, cluster2.moat);
            new_cluster.subcluster_moat_sum += cluster2.moat;
            // Potentially log inactive merge event for GW pruning
            if (!cluster2.contains_root) { // Only log if inactive side is not the root component
                 format_and_log(4, "  Logging inactive merge event for edge %d (Inactive child %d).\n", info.edge_idx, cluster2_idx);
                 log_inactive_merge_event(info);
                 // Adjust heap keys for inactive cluster's edges based on time shift
                 ValueType edge_event_update_shift = current_time_ - cluster2.active_end_time;
                 if (edge_event_update_shift > kEpsilon) {
                     format_and_log(3, "  Shifting edge event times for inactive child cluster %d by %.9g.\n", cluster2_idx, edge_event_update_shift);
                     // This add_to_heap propagates the time shift to edge parts within the inactive cluster's heap
                     cluster2.edge_parts.add_to_heap(edge_event_update_shift);
                 }
            } else {
                 format_and_log(4, "  Skipping inactive merge event log/shift as inactive child %d contains root.\n", cluster2_idx);
            }
        }
        format_and_log(4, "  New cluster %d updated SubclusterMoatSum=%.9g\n", new_cluster_idx, new_cluster.subcluster_moat_sum);

        // --- Meld Pairing Heaps ---
        format_and_log(4, "Melding heaps from cluster %d and %d into %d.\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        // Move heaps from children into the new cluster
        new_cluster.edge_parts = PairingHeapType::meld(std::move(cluster1.edge_parts), std::move(cluster2.edge_parts));
        format_and_log(4, "  Meld complete. New cluster %d heap is%s empty.\n", new_cluster_idx, new_cluster.edge_parts.is_empty() ? "" : " NOT");

        // --- Update New Cluster Activity and PQs ---
        if (new_cluster.active) {
            new_cluster.active_start_time = current_time_;
            // Calculate potential deactivation time based on total prize and accumulated costs (moats)
            ValueType potential_deactivation_time = new_cluster.active_start_time + new_cluster.prize_sum - new_cluster.subcluster_moat_sum;
            format_and_log(4, "  New cluster %d is active. Potential deactivation time = %.9g + %.9g - %.9g = %.9g\n",
                        new_cluster_idx, new_cluster.active_start_time, new_cluster.prize_sum, new_cluster.subcluster_moat_sum, potential_deactivation_time);
            // Only insert if deactivation time is meaningfully after current time
            if (potential_deactivation_time > current_time_ + kEpsilon) {
                 format_and_log(3, "  Scheduling deactivation for new cluster %d at time %.9g.\n", new_cluster_idx, potential_deactivation_time);
                 clusters_deactivation_.insert(potential_deactivation_time, new_cluster_idx);
            } else {
                 // If potential deactivation is now or in the past, make it immediately inactive
                 format_and_log(3, "  New cluster %d immediately inactive (potential time %.9g <= current time %.9g).\n", new_cluster_idx, potential_deactivation_time, current_time_);
                 new_cluster.active = false;
                 new_cluster.active_end_time = current_time_;
                 new_cluster.moat = 0.0; // Immediately inactive means zero moat duration
                 new_cluster.subcluster_moat_sum += new_cluster.moat; // Add zero moat (effectively no change)
            }
        } else { // Not active (contains root or became immediately inactive)
            new_cluster.active_start_time = current_time_; // Still note the merge time
            new_cluster.active_end_time = current_time_;
            new_cluster.moat = 0.0;
            new_cluster.subcluster_moat_sum += new_cluster.moat; // Add zero moat
             format_and_log(3, "  New cluster %d contains root or became immediately inactive, setting Active=false.\n", new_cluster_idx);
        }

         // Update PQs (remove children, add new cluster if needed)
         format_and_log(4, "Updating PQs after merge (removing children %d, %d; adding new %d if needed).\n", cluster1_idx, cluster2_idx, new_cluster_idx);
        update_cluster_queues_post_merge(new_cluster_idx, cluster1_idx, cluster2_idx);

         format_and_log(3, "merge_clusters Exit.\n");
    }

    void PCSTFast::log_inactive_merge_event(const EdgeProcessingInfo& info) {
         format_and_log(4, "log_inactive_merge_event Entry for Edge %d.\n", info.edge_idx);
         inactive_merge_events_.emplace_back(); // Add new event entry
         InactiveMergeEvent& merge_event = inactive_merge_events_.back(); // Get reference

         // Store representative cluster indices *at the time of the merge*
         // Assume info.current corresponds to the cluster that triggered the event (the active one)
         merge_event.active_cluster_index = info.current_cluster_idx;
         merge_event.inactive_cluster_index = info.other_cluster_idx;

         // Store original node indices from the edge causing the merge
         assert(static_cast<size_t>(info.edge_idx) < edges_.size() && "Invalid edge index in log_inactive");
         const auto& edge = edges_[info.edge_idx];
         // Determine which original node corresponds to which cluster side based on edge part index
         IndexType node_on_current_side = (info.current_part_idx % 2 == 0) ? edge.first : edge.second;
         IndexType node_on_other_side = (info.other_part_idx % 2 == 0) ? edge.first : edge.second;

         merge_event.active_cluster_node = node_on_current_side;
         merge_event.inactive_cluster_node = node_on_other_side;

         // Link the edge index back to this event in edge_info_
         if (static_cast<size_t>(info.edge_idx) < edge_info_.size()) {
             edge_info_[info.edge_idx].inactive_merge_event = std::ssize(inactive_merge_events_) - 1; // Store index of the event we just added
             format_and_log(3, "Logged InactiveMergeEvent %zd for edge %d (Active: C%d/N%d, Inactive: C%d/N%d).\n",
                         inactive_merge_events_.size() - 1, info.edge_idx,
                         merge_event.active_cluster_index, merge_event.active_cluster_node,
                         merge_event.inactive_cluster_index, merge_event.inactive_cluster_node);
         } else {
              format_and_log(0, "Error: Edge index %d out of bounds for edge_info_ in log_inactive_merge_event.\n", info.edge_idx);
               assert(false && "Edge index out of bounds for edge_info_");
         }
          format_and_log(4, "log_inactive_merge_event Exit.\n");
     }

    void PCSTFast::handle_active_active_growth(const EdgeProcessingInfo& info) {
        stats_.num_active_active_edge_growth_events++;
        // Calculate the time when the edge *would* fully merge if both continue growing at equal rates
        ValueType next_event_time = current_time_ + info.remainder * 0.5;
        format_and_log(3, "handle_active_active_growth Entry: Edge=%d, Clusters=(%d, %d), Remainder=%.9g, NextEventTime=%.9g\n",
                    info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);

        // Basic assertions
        assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size());
        assert(static_cast<size_t>(info.current_cluster_idx) < clusters_.size());
        assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size());

        EdgePart& current_part = edge_parts_[info.current_part_idx];
        EdgePart& other_part = edge_parts_[info.other_part_idx];
        Cluster& current_cluster = clusters_[info.current_cluster_idx];
        Cluster& other_cluster = clusters_[info.other_cluster_idx];

        // Update and re-insert/decrease key for the current part (which triggered the event)
        format_and_log(4, "Updating current part %d (cluster %d): Inserting into heap with time %.9g.\n",
                    info.current_part_idx, info.current_cluster_idx, next_event_time);
        current_part.next_event_val = next_event_time;
        // Insert the edge part back into the heap with the new event time
        current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
        // Update the global PQ for this cluster, as its min edge might have changed
        update_cluster_queues_post_growth(info.current_cluster_idx);

        // Update and decrease key for the other part
        format_and_log(4, "Updating other part %d (cluster %d): Decreasing key in heap to time %.9g.\n",
                    info.other_part_idx, info.other_cluster_idx, next_event_time);
        other_part.next_event_val = next_event_time;
        // Need to remove old entry from global PQ *before* decrease_key,
        // as decrease_key might change the cluster's minimum.
        format_and_log(4, "  Deleting cluster %d from global edge PQ before decrease_key.\n", info.other_cluster_idx);
        clusters_next_edge_event_.delete_element(info.other_cluster_idx);
        // Decrease the key in the other cluster's heap
        other_cluster.edge_parts.decrease_key(other_part.heap_node, next_event_time);
        // Update the global PQ for the other cluster with its potentially new minimum
        update_cluster_queues_post_growth(info.other_cluster_idx);

         format_and_log(3, "handle_active_active_growth Exit.\n");
    }

    void PCSTFast::handle_active_inactive_growth(const EdgeProcessingInfo& info) {
         stats_.num_active_inactive_edge_growth_events++;
         // Calculate the time when the edge *would* fully merge if the active side covers the remainder
         ValueType next_event_time = current_time_ + info.remainder;
          format_and_log(3, "handle_active_inactive_growth Entry: Edge=%d, Clusters=(Active:%d, Inactive:%d), Remainder=%.9g, NextEventTime=%.9g\n",
                     info.edge_idx, info.current_cluster_idx, info.other_cluster_idx, info.remainder, next_event_time);

         // Basic assertions
         assert(static_cast<size_t>(info.current_part_idx) < edge_parts_.size());
         assert(static_cast<size_t>(info.other_part_idx) < edge_parts_.size());
         assert(static_cast<size_t>(info.current_cluster_idx) < clusters_.size());
         assert(static_cast<size_t>(info.other_cluster_idx) < clusters_.size());
         assert(!clusters_[info.other_cluster_idx].active && "Other cluster should be inactive");


         EdgePart& current_part = edge_parts_[info.current_part_idx];
         EdgePart& other_part = edge_parts_[info.other_part_idx];
         Cluster& current_cluster = clusters_[info.current_cluster_idx]; // Active cluster
         Cluster& other_cluster = clusters_[info.other_cluster_idx]; // Inactive cluster

         // Update and re-insert for the active side
         format_and_log(4, "Updating active part %d (cluster %d): Inserting into heap with time %.9g.\n",
                     info.current_part_idx, info.current_cluster_idx, next_event_time);
         current_part.next_event_val = next_event_time;
         current_part.heap_node = current_cluster.edge_parts.insert(next_event_time, info.current_part_idx);
         update_cluster_queues_post_growth(info.current_cluster_idx);

         // Reset the inactive side's event time in its heap and update its state
         // The effective event time becomes the time the inactive cluster *finished* being active.
         ValueType inactive_reset_time = other_cluster.active_end_time;
         format_and_log(4, "Updating inactive part %d (cluster %d): Decreasing key in heap to time %.9g (cluster active_end_time).\n",
                     info.other_part_idx, info.other_cluster_idx, inactive_reset_time);
         other_part.next_event_val = inactive_reset_time;
         // Decrease key in the inactive cluster's heap. Note: This heap might have been time-shifted in merge_clusters.
         other_cluster.edge_parts.decrease_key(other_part.heap_node, inactive_reset_time);
         // No need to update global edge event PQ for the inactive cluster as it's not tracked there.

          format_and_log(3, "handle_active_inactive_growth Exit.\n");
    }

    // --- Cluster Deactivation Event ---
    void PCSTFast::handle_cluster_deactivation_event(IndexType cluster_index) {
         format_and_log(3, "handle_cluster_deactivation_event Entry: Cluster=%d\n", cluster_index);
        // Validate index
        if (cluster_index < 0 || static_cast<size_t>(cluster_index) >= clusters_.size()) {
            format_and_log(0,"Error: Invalid cluster index %d in deactivation event.\n", cluster_index);
             assert(false && "Invalid cluster index in deactivation event");
            return;
        }
        Cluster& cluster = clusters_[cluster_index];

        // Check if already inactive (e.g., due to a merge event happening at the same time)
        if (!cluster.active) {
            format_and_log(2, "Cluster %d already inactive, skipping deactivation event logic.\n", cluster_index);
            return; // Already deactivated, nothing to do
        }

        // --- Deactivate Cluster ---
        cluster.active = false;
        cluster.active_end_time = current_time_; // Mark end time
        // Calculate the duration it was active (moat size)
        cluster.moat = cluster.active_end_time - cluster.active_start_time;
        // Clamp moat if time precision issues occur
        if (cluster.moat < 0.0) {
             format_and_log(1,"Warning: Negative moat calculated for cluster %d (%.9g). Clamping to 0.\n", cluster_index, cluster.moat);
             cluster.moat = 0.0;
        }

        format_and_log(2, "Cluster Deactivation: Cluster %d at time %.9g (Moat size %.9g).\n",
                    cluster_index, current_time_, cluster.moat);

        // --- Remove from Edge Event PQ ---
        // If the cluster had active edges, remove its entry from the global edge event PQ,
        // as its minimum edge event time is no longer relevant for scheduling.
        if (!cluster.edge_parts.is_empty()) {
             format_and_log(4, "Removing cluster %d from next edge event PQ as it became inactive.\n", cluster_index);
            clusters_next_edge_event_.delete_element(cluster_index);
        } else {
             format_and_log(4, "Cluster %d had no edges, no need to remove from edge PQ.\n", cluster_index);
        }

        // Note: Cluster was already removed from the *deactivation* PQ by process_event.

         format_and_log(3, "handle_cluster_deactivation_event Exit.\n");
    }

    // --- Priority Queue Updates ---
    void PCSTFast::update_cluster_queues_post_growth(IndexType cluster_idx) {
         format_and_log(4, "update_cluster_queues_post_growth Entry: Cluster=%d\n", cluster_idx);
        // Only update for valid, *active* clusters
        if (cluster_idx == kInvalidIndex || static_cast<size_t>(cluster_idx) >= clusters_.size() || !clusters_[cluster_idx].active) {
             format_and_log(4, "  Skipping PQ update (cluster invalid or inactive).\n");
            return;
        }
        Cluster& cluster = clusters_[cluster_idx];

        // Remove existing entry first (if any) to handle potential key changes
        format_and_log(4, "  Deleting cluster %d from global edge PQ (if present).\n", cluster_idx);
        clusters_next_edge_event_.delete_element(cluster_idx);

        // Insert new entry if the cluster still has edges
        if (!cluster.edge_parts.is_empty()) {
            ValueType min_val; PayloadType min_part;
            // Check return value of get_min
            if (cluster.edge_parts.get_min(&min_val, &min_part)) {
                 format_and_log(4, "  Cluster %d heap not empty. Inserting into edge PQ with Time=%.9g, Part=%d.\n", cluster_idx, min_val, min_part);
                clusters_next_edge_event_.insert(min_val, cluster_idx); // Insert {new_min_time, cluster_idx}
            } else {
                 format_and_log(0, "Error: Cluster %d heap not empty but get_min failed in update_queues.\n", cluster_idx);
                 assert(false && "get_min failed on non-empty heap"); // Should not happen
            }
        } else {
             format_and_log(4, "  Cluster %d heap is empty, not inserting into edge PQ.\n", cluster_idx);
        }
         format_and_log(4, "update_cluster_queues_post_growth Exit.\n");
    }

     void PCSTFast::update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index) {
          format_and_log(4, "update_cluster_queues_post_merge Entry: New=%d, Children=(%d, %d)\n", new_cluster_index, child1_index, child2_index);

         // Remove children from both PQs (they are no longer active representatives)
         format_and_log(4, "  Deleting child %d from PQs (if present).\n", child1_index);
         clusters_deactivation_.delete_element(child1_index);
         clusters_next_edge_event_.delete_element(child1_index);
          format_and_log(4, "  Deleting child %d from PQs (if present).\n", child2_index);
         clusters_deactivation_.delete_element(child2_index);
         clusters_next_edge_event_.delete_element(child2_index);

         // Add the new cluster to the edge event PQ *if* it's active and has edges
         if (new_cluster_index != kInvalidIndex && static_cast<size_t>(new_cluster_index) < clusters_.size()) {
             Cluster& new_cluster = clusters_[new_cluster_index];
             if (new_cluster.active && !new_cluster.edge_parts.is_empty()) {
                 ValueType min_val; PayloadType min_part;
                 if (new_cluster.edge_parts.get_min(&min_val, &min_part)) {
                     format_and_log(4, "  New cluster %d is active with edges. Inserting into edge PQ with Time=%.9g, Part=%d.\n", new_cluster_index, min_val, min_part);
                     clusters_next_edge_event_.insert(min_val, new_cluster_index);
                 } else {
                     format_and_log(0, "Error: New cluster %d heap not empty but get_min failed in update_queues_post_merge.\n", new_cluster_index);
                     assert(false && "get_min failed on non-empty heap post-merge");
                 }
             } else {
                 format_and_log(4, "  New cluster %d is inactive or has no edges, not inserting into edge PQ.\n", new_cluster_index);
             }
             // Note: Deactivation event for new_cluster was potentially added within merge_clusters logic.
         } else {
              format_and_log(4, "  New cluster index %d is invalid, skipping edge PQ update.\n", new_cluster_index);
         }
           format_and_log(4, "update_cluster_queues_post_merge Exit.\n");
     }

    // --- Path Compression ---
    void PCSTFast::get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                                      ValueType* finished_moat_sum, IndexType* current_cluster_index) {
        format_and_log(4, "get_sum_on_edge_part Entry: EdgePart=%d\n", edge_part_index);
        // Determine the starting node index from the edge part index
        IndexType edge_idx = edge_part_index / 2;
        assert(static_cast<size_t>(edge_idx) < edges_.size() && "Invalid edge index in get_sum");
        const auto& edge = edges_[edge_idx];
        IndexType start_node = (edge_part_index % 2 == 0) ? edge.first : edge.second;
        format_and_log(4, "  Start node derived from edge part: %d\n", start_node);

        // Initialize outputs
        *total_sum = 0.0;
        *finished_moat_sum = 0.0;
        *current_cluster_index = kInvalidIndex; // Assume failure initially

        // Find the representative cluster and perform path compression
        ValueType path_sum = 0.0;
        format_and_log(4, "  Calling find_representative_and_compress for start node %d...\n", start_node);
        IndexType representative_cluster = find_representative_and_compress(start_node, path_sum);
        format_and_log(4, "  find_representative result: RepCluster=%d, PathSum=%.9g\n", representative_cluster, path_sum);

        // Check if traversal ended successfully
        if (representative_cluster == kInvalidIndex || static_cast<size_t>(representative_cluster) >= clusters_.size()) {
            format_and_log(0, "Error: Path compression failed to find representative for node %d (edge part %d).\n", start_node, edge_part_index);
            // Leave outputs initialized to failure state
            format_and_log(4, "get_sum_on_edge_part Exit (Failure).\n");
            return; // Failure case
        }

        *current_cluster_index = representative_cluster;
        const Cluster& rep_cluster = clusters_[representative_cluster];
        format_and_log(4, "  Representative cluster %d is Active=%d (StartTime=%.9g, EndTime=%.9g, Moat=%.9g)\n",
                     representative_cluster, rep_cluster.active, rep_cluster.active_start_time, rep_cluster.active_end_time, rep_cluster.moat);

        // Calculate total sum based on whether the representative is active
        if (rep_cluster.active) {
            *finished_moat_sum = path_sum; // Sum only includes moats up to the active cluster
            // Add current active duration
            *total_sum = path_sum + (current_time_ - rep_cluster.active_start_time);
             format_and_log(4, "  Rep is active. FinishedMoatSum=%.9g, TotalSum = %.9g + (%.9g - %.9g) = %.9g\n",
                         *finished_moat_sum, path_sum, current_time_, rep_cluster.active_start_time, *total_sum);
        } else {
            // Add final moat of the inactive representative
            *total_sum = path_sum + rep_cluster.moat;
            *finished_moat_sum = *total_sum; // For inactive, finished moat is total sum
             format_and_log(4, "  Rep is inactive. FinishedMoatSum=TotalSum = %.9g + %.9g = %.9g\n",
                         path_sum, rep_cluster.moat, *total_sum);
        }
          format_and_log(4, "get_sum_on_edge_part Exit (Success). Cluster=%d, TotalSum=%.9g, FinishedMoat=%.9g\n",
                      *current_cluster_index, *total_sum, *finished_moat_sum);
    }

    PCSTFast::IndexType PCSTFast::find_representative_and_compress(IndexType start_node, ValueType& path_sum_out) {
        format_and_log(4, "find_representative_and_compress Entry: StartNode=%d\n", start_node);
        // Use a local buffer for storing visited nodes during traversal
        std::vector<std::pair<IndexType, ValueType>> path_compression_visited_local;
        // Reserve capacity based on a reasonable estimate (number of clusters)
        path_compression_visited_local.reserve(clusters_.size());

        path_sum_out = 0.0; // Initialize accumulated path sum
        IndexType current_node = start_node;
        int steps = 0;

        // Traverse up the merge tree until a representative (merged_into == kInvalidIndex) is found
        while (static_cast<size_t>(current_node) < clusters_.size() && clusters_[current_node].merged_into != kInvalidIndex) {
            steps++;
            format_and_log(4, "  Step %d: CurrentNode=%d, PathSum=%.9g. Storing locally.\n", steps, current_node, path_sum_out);
            // Store the node and the path sum *before* traversing from it
            path_compression_visited_local.emplace_back(current_node, path_sum_out);

            ValueType step_sum;
            IndexType next_node;
            bool used_skip = false;

            // Check bounds before accessing cluster data
             assert(static_cast<size_t>(current_node) < clusters_.size() && "Current node index out of bounds in find_rep");
            const Cluster& curr_cluster_ref = clusters_[current_node];

            // Use skip pointers if available for faster traversal
            if (curr_cluster_ref.skip_up >= 0) {
                step_sum = curr_cluster_ref.skip_up_sum;
                next_node = curr_cluster_ref.skip_up;
                used_skip = true;
            } else { // Otherwise, take one step up using merged_into
                step_sum = curr_cluster_ref.moat; // Cost accumulated in this cluster
                next_node = curr_cluster_ref.merged_into; // Move to parent
            }
            format_and_log(4, "  Step %d: NextNode=%d, StepSum=%.9g (UsedSkip=%d)\n", steps, next_node, step_sum, used_skip);

            path_sum_out += step_sum; // Accumulate sum
            current_node = next_node; // Move up
        }
        format_and_log(4, "  Traversal finished after %d steps. FinalNode=%d, FinalPathSum=%.9g\n", steps, current_node, path_sum_out);

        // Check if traversal ended successfully at a valid representative
        if (static_cast<size_t>(current_node) >= clusters_.size()) {
            format_and_log(0, "Error: Path compression traversed out of bounds from start node %d.\n", start_node);
            assert(false && "Path compression out of bounds");
            return kInvalidIndex; // Return invalid index on failure
        }

        IndexType representative_cluster = current_node;
        format_and_log(4, "  Representative cluster found: %d\n", representative_cluster);
        format_and_log(4, "  Applying path compression for %zu visited nodes (updating member clusters_)...\n", path_compression_visited_local.size());

        // Apply path compression: Update skip pointers for all nodes visited during traversal
        for (const auto& [visited_node_idx, sum_at_visited] : path_compression_visited_local) {
            // Check bounds before updating member cluster data
             assert(static_cast<size_t>(visited_node_idx) < clusters_.size() && "Invalid visited node index during compression");
             // Calculate the sum accumulated from the visited node to the final representative
            ValueType new_skip_sum = path_sum_out - sum_at_visited;
            format_and_log(4, "    Updating skip pointer for Node %d: SkipUp=%d, SkipSum=%.9g\n",
                        visited_node_idx, representative_cluster, new_skip_sum);
            // Update the skip pointers in the main clusters_ vector
            clusters_[visited_node_idx].skip_up = representative_cluster;
            clusters_[visited_node_idx].skip_up_sum = new_skip_sum;
        }

        format_and_log(4, "find_representative_and_compress Exit. Rep=%d, PathSum=%.9g\n", representative_cluster, path_sum_out);
        return representative_cluster; // Return the found representative index
     }

    // --- Select Initial Good Nodes ---
    void PCSTFast::select_initial_active_clusters(std::vector<uint8_t>& node_good_flag) {
         format_and_log(3, "select_initial_active_clusters Entry.\n");
         // Use a local buffer for the BFS/DFS traversal queue
         std::vector<IndexType> cluster_queue_local;
         cluster_queue_local.reserve(clusters_.size()); // Reserve capacity

        // Determine starting cluster(s) based on rooted or unrooted case
        if (root_ >= 0) {
             format_and_log(3, "  Rooted case (Root=%d).\n", root_);
            IndexType root_cluster_final_index = root_;
             format_and_log(4, "  Finding final representative for root %d...\n", root_);
             int path_len = 0;
             // Follow merged_into pointers to find the final cluster containing the root
             while(static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].merged_into != kInvalidIndex) {
                 assert(path_len < static_cast<int>(clusters_.size()) && "Infinite loop detected in merge path?");
                 root_cluster_final_index = clusters_[root_cluster_final_index].merged_into;
                 path_len++;
             }
              format_and_log(4, "  Final representative for root %d is cluster %d (path length %d).\n", root_, root_cluster_final_index, path_len);

             // Check if the found representative actually contains the root and is valid
             if (static_cast<size_t>(root_cluster_final_index) < clusters_.size() && clusters_[root_cluster_final_index].contains_root) {
                 format_and_log(2, "  Rooted case: Starting traversal from final root cluster %d.\n", root_cluster_final_index);
                 cluster_queue_local.push_back(root_cluster_final_index); // Add representative to local queue
             } else {
                 format_and_log(0, "Warning: Could not find final cluster containing root %d. No nodes selected initially.\n", root_);
             }
        } else { // Unrooted case
             format_and_log(3, "  Unrooted case.\n");
             format_and_log(2, "  Unrooted case: Starting traversal from all final active clusters.\n");
             // Find all clusters that were never merged AND remained active at the end
             for (size_t ii = 0; ii < clusters_.size(); ++ii) {
                 if (clusters_[ii].merged_into == kInvalidIndex && clusters_[ii].active) {
                      format_and_log(3, "  Adding final active cluster %zu to initial queue.\n", ii);
                     cluster_queue_local.push_back(static_cast<IndexType>(ii)); // Add to local queue
                 }
             }
        }

         format_and_log(3, "  Starting traversal from %zu initial clusters to mark good nodes...\n", cluster_queue_local.size());
         size_t queue_index = 0;
         int nodes_marked_good = 0;
         // Perform BFS/DFS using the local queue to traverse down the merge tree
         while (queue_index < cluster_queue_local.size()) {
             IndexType current_cluster_index = cluster_queue_local[queue_index++];
             format_and_log(4, "    Processing cluster %d from queue (index %zu).\n", current_cluster_index, queue_index - 1);

             // Check for invalid indices (safety measure)
             if (static_cast<size_t>(current_cluster_index) >= clusters_.size()) {
                  format_and_log(1, "Warning: Invalid cluster index %d encountered during good node traversal.\n", current_cluster_index);
                   assert(false && "Invalid cluster index in good node traversal");
                  continue;
             }

             const Cluster& current_cluster = clusters_[current_cluster_index];

             // If it's a merged cluster, add its children (representatives at merge time) to the queue
             if (current_cluster.merged_along >= 0) { // Indicates it's a merged cluster
                  format_and_log(4, "    Cluster %d is merged. Adding children %d, %d to queue.\n",
                              current_cluster_index, current_cluster.child_cluster_1, current_cluster.child_cluster_2);
                 if(current_cluster.child_cluster_1 != kInvalidIndex) cluster_queue_local.push_back(current_cluster.child_cluster_1);
                 if(current_cluster.child_cluster_2 != kInvalidIndex) cluster_queue_local.push_back(current_cluster.child_cluster_2);
             } else { // If it's an original node cluster (index < num_nodes), mark the node as good
                  format_and_log(4,"    Cluster %d is original node cluster.\n", current_cluster_index);
                  // Check index validity before marking the flag
                  if (static_cast<size_t>(current_cluster_index) < node_good_flag.size()) {
                      // Mark only if not already marked
                      if(node_good_flag[current_cluster_index] == 0) {
                           format_and_log(3, "    Marking node %d as good.\n", current_cluster_index);
                           node_good_flag[current_cluster_index] = 1; // Set the flag in the member vector
                           nodes_marked_good++;
                      } else {
                           format_and_log(4, "    Node %d already marked good.\n", current_cluster_index);
                      }
                  } else {
                       format_and_log(0, "Error: Original cluster index %d out of bounds for node_good_flag (%zu).\n", current_cluster_index, node_good_flag.size());
                        assert(false && "node_good_flag index out of bounds");
                  }
             }
         }
         format_and_log(3, "select_initial_active_clusters Exit. Marked %d nodes as good.\n", nodes_marked_good);
    }


    // --- Statistics ---
    void PCSTFast::get_statistics(Statistics* s) const {
        format_and_log(3, "get_statistics Entry.\n");
        if(s) {
            *s = stats_; // Copy collected statistics
             format_and_log(3, "  Statistics copied.\n");
        } else {
             format_and_log(1, "Warning: Null pointer passed to get_statistics.\n");
             // Optionally throw or just log and continue
        }
         format_and_log(3, "get_statistics Exit.\n");
    }

} // namespace cluster_approx
