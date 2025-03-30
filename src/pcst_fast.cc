#include "pcst_fast.h"

#include <algorithm> 
#include <cmath>
#include <cctype>
#include <cstdio>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>


namespace cluster_approx {
    PCSTFast::PruningMethod PCSTFast::parse_pruning_method(std::string_view input) {
        std::string lower_input;
        lower_input.reserve(input.length());

        for (char c : input) {
            lower_input += static_cast<char>(::tolower(static_cast<unsigned char>(c)));
        }

        if (lower_input == "none") return PruningMethod::kNoPruning;
        if (lower_input == "simple") return PruningMethod::kSimplePruning;
        if (lower_input == "gw") return PruningMethod::kGWPruning;
        if (lower_input == "strong") return PruningMethod::kStrongPruning;

        return PruningMethod::kUnknownPruning;
    }

    PCSTFast::PCSTFast(const std::vector<std::pair<int, int> >& edges_,
                    const std::vector<double>& prizes_,
                    const std::vector<double>& costs_,
                    int root_,
                    int target_num_active_clusters_,
                    PruningMethod pruning_,
                    int verbosity_level_,
                    void (*output_function_)(const char*))
        : edges(edges_), prizes(prizes_), costs(costs_), root(root_),
        target_num_active_clusters( (root_ >= 0 && target_num_active_clusters_ > 0) ? 0 : target_num_active_clusters_),
        pruning(pruning_), verbosity_level(verbosity_level_),
        output_function(output_function_)
    {
        const size_t num_nodes = prizes.size();
        const size_t num_edges = edges.size();

        if (root >= static_cast<int>(num_nodes)) {
            throw std::invalid_argument("Root node index out of range.");
        }
        if (target_num_active_clusters_ < 0) {
            throw std::invalid_argument("Target number of active clusters cannot be negative.");
        }
        if (root >= 0 && target_num_active_clusters_ > 0) {
            if (output_function) {
                snprintf(output_buffer, kOutputBufferSize, "Warning: target_num_active_clusters > 0 is ignored in the rooted case.\n");
                output_function(output_buffer);
            }
        }

        edge_parts.resize(2 * num_edges);
        edge_info.resize(num_edges);
        node_deleted.resize(num_nodes, 0);
        node_good.resize(num_nodes, 0);
        build_phase1_included_nodes.resize(num_nodes, 0);

        clusters.reserve(2 * num_nodes);
        inactive_merge_events.reserve(num_edges);
        phase1_result.reserve(num_edges);
        phase2_result.reserve(num_edges);
        phase3_result.reserve(num_edges);
        path_compression_visited.reserve(num_nodes);
        cluster_queue.reserve(num_nodes);

        if (pruning == PruningMethod::kStrongPruning || pruning == PruningMethod::kGWPruning) {
            phase3_neighbors.resize(num_nodes);
        }

        if (pruning == PruningMethod::kStrongPruning) {
            final_component_label.resize(num_nodes);
            strong_pruning_parent.resize(num_nodes);
            strong_pruning_payoff.resize(num_nodes);
            final_components.reserve(num_nodes);
            stack.reserve(num_nodes);
            stack2.reserve(num_nodes);
        }

        for (auto& info : edge_info) {
            info.inactive_merge_event = -1;
        }

        current_time = 0.0;
        eps = 1e-9;

        for (int ii = 0; ii < std::ssize(prizes); ++ii) {
            if (prizes[ii] < 0.0) {
                throw std::invalid_argument("Prize negative.");
            }

            clusters.emplace_back(&pairing_heap_buffer);

            Cluster& current_cluster = clusters.back();
            current_cluster.active = (ii != root);
            current_cluster.active_start_time = 0.0;
            current_cluster.active_end_time = (current_cluster.active) ? -1.0 : 0.0;
            current_cluster.prize_sum = prizes[ii];
            current_cluster.contains_root = (ii == root);

            if (current_cluster.active) {
                clusters_deactivation.insert(prizes[ii], ii);
            }
        }

        for (int ii = 0; ii < std::ssize(edges); ++ii) {
            const auto& edge = edges[ii];
            int uu = edge.first;
            int vv = edge.second;

            if (uu < 0 || vv < 0 || uu >= std::ssize(prizes) || vv >= std::ssize(prizes)) {
                throw std::invalid_argument("Edge endpoint index out of range.");
            }

            if (costs[ii] < 0.0) {
                throw std::invalid_argument("Edge cost negative.");
            }
            
            if (uu == vv) {
                continue;
            }

            double cost = costs[ii];
            EdgePart& uu_part = edge_parts[2 * ii];
            EdgePart& vv_part = edge_parts[2 * ii + 1];
            Cluster& uu_cluster = clusters[uu];
            Cluster& vv_cluster = clusters[vv];
            uu_part.deleted = false;
            vv_part.deleted = false;
            double event_time_u, event_time_v;
            
            if (uu_cluster.active && vv_cluster.active) {
                event_time_u = cost * 0.5;
                event_time_v = cost * 0.5;
            } else if (uu_cluster.active) {
                event_time_u = cost;
                event_time_v = 0.0;
            } else if (vv_cluster.active) {
                event_time_u = 0.0;
                event_time_v = cost;
            } else {
                event_time_u = 0.0;
                event_time_v = 0.0;
            }

            uu_part.next_event_val = event_time_u;
            vv_part.next_event_val = event_time_v;
            uu_part.heap_node = uu_cluster.edge_parts.insert(uu_part.next_event_val, 2 * ii);
            vv_part.heap_node = vv_cluster.edge_parts.insert(vv_part.next_event_val, 2 * ii + 1);
        }

        for (int ii = 0; ii < std::ssize(clusters) && ii < std::ssize(prizes) ; ++ii) {
            if (clusters[ii].active) {
                if (!clusters[ii].edge_parts.is_empty()) {
                    double val;
                    int edge_part;
        
                    (void) clusters[ii].edge_parts.get_min(&val, &edge_part);
                    clusters_next_edge_event.insert(val, ii);
                }
            }
        }
    }

    PCSTFast::~PCSTFast() = default;

    void PCSTFast::get_next_edge_event(double* next_time, int* next_cluster_index, int* next_edge_part_index) {
        if (clusters_next_edge_event.is_empty()) {
            *next_time = std::numeric_limits<double>::infinity();
            *next_cluster_index = -1;
            *next_edge_part_index = -1;
        
            return;
        }

        (void) clusters_next_edge_event.get_min(next_time, next_cluster_index);

        if (*next_cluster_index != -1 && static_cast<size_t>(*next_cluster_index) < clusters.size()) {
            (void) clusters[*next_cluster_index].edge_parts.get_min(next_time, next_edge_part_index);
        } else {
            *next_time = std::numeric_limits<double>::infinity();
            *next_edge_part_index = -1;
        }
    }

    void PCSTFast::remove_next_edge_event(int cluster_index) {
        if (static_cast<size_t>(cluster_index) >= clusters.size() || !clusters[cluster_index].active) {
            clusters_next_edge_event.delete_element(cluster_index);
        
            return;
        }

        clusters_next_edge_event.delete_element(cluster_index);

        double tmp_value;
        int tmp_edge_part;
        
        (void) clusters[cluster_index].edge_parts.delete_min(&tmp_value, &tmp_edge_part);

        if (!clusters[cluster_index].edge_parts.is_empty()) {
            (void) clusters[cluster_index].edge_parts.get_min(&tmp_value, &tmp_edge_part);
            clusters_next_edge_event.insert(tmp_value, cluster_index);
        }
    }

    void PCSTFast::get_next_cluster_event(double* next_time, int* next_cluster_index) {
        if (clusters_deactivation.is_empty()) {
            *next_time = std::numeric_limits<double>::infinity();
            *next_cluster_index = -1;
        
            return;
        }

        (void) clusters_deactivation.get_min(next_time, next_cluster_index);
    }

    void PCSTFast::remove_next_cluster_event() {
        if (!clusters_deactivation.is_empty()) {
            double tmp_value;
            int tmp_cluster;
            
            (void) clusters_deactivation.delete_min(&tmp_value, &tmp_cluster);
        }
    }

    void PCSTFast::get_sum_on_edge_part(int edge_part_index, double* total_sum, double* finished_moat_sum, int* current_cluster_index) {
        const auto& edge = edges[edge_part_index / 2];
        int endpoint = (edge_part_index % 2 == 0) ? edge.first : edge.second;
        *total_sum = 0.0;
        *finished_moat_sum = 0.0;
        *current_cluster_index = endpoint;

        path_compression_visited.clear();

        int path_node = endpoint;
        double path_sum = 0.0;

        while (static_cast<size_t>(path_node) < clusters.size() && clusters[path_node].merged_into != -1) {
            path_compression_visited.emplace_back(path_node, path_sum);
        
            int next_node;
            double step_sum;

            if (clusters[path_node].skip_up >= 0) {
                step_sum = clusters[path_node].skip_up_sum;
                next_node = clusters[path_node].skip_up;
            } else {
                step_sum = clusters[path_node].moat;
                next_node = clusters[path_node].merged_into;
            }
        
            path_sum += step_sum;
            path_node = next_node;
        }

        if (static_cast<size_t>(path_node) >= clusters.size()) {
            *current_cluster_index = -1;
            *total_sum = 0.0;
            *finished_moat_sum = 0.0;
        
            return;
        }

        *current_cluster_index = path_node;

        for (const auto& [visited_cluster_index, sum_at_visited] : path_compression_visited) {
            if (static_cast<size_t>(visited_cluster_index) < clusters.size()) {
                clusters[visited_cluster_index].skip_up = *current_cluster_index;
                clusters[visited_cluster_index].skip_up_sum = path_sum - sum_at_visited;
            }
        }

        if (clusters[*current_cluster_index].active) {
            *finished_moat_sum = path_sum;
            *total_sum = path_sum + (current_time - clusters[*current_cluster_index].active_start_time);
        } else {
            *total_sum = path_sum + clusters[*current_cluster_index].moat;
            *finished_moat_sum = *total_sum;
        }
    }

    void PCSTFast::mark_nodes_as_good(int start_cluster_index) {
        cluster_queue.clear();

        if (static_cast<size_t>(start_cluster_index) >= clusters.size()) return;

        cluster_queue.push_back(start_cluster_index);

        size_t queue_index = 0;
        
        while (queue_index < cluster_queue.size()) {
            int current_cluster_index = cluster_queue[queue_index++];

            if (static_cast<size_t>(current_cluster_index) >= clusters.size()) continue;

            const Cluster& current_cluster = clusters[current_cluster_index];

            if (current_cluster.merged_along >= 0) {
                if(current_cluster.child_cluster_1 != -1) cluster_queue.push_back(current_cluster.child_cluster_1);
                if(current_cluster.child_cluster_2 != -1) cluster_queue.push_back(current_cluster.child_cluster_2);
            } else {
                if (static_cast<size_t>(current_cluster_index) < node_good.size()) {
                    node_good[current_cluster_index] = 1;
                }
            }
        }
    }

    void PCSTFast::mark_clusters_as_necessary(int start_cluster_index) {
        int current_cluster_index = start_cluster_index;

        while (current_cluster_index != -1 && static_cast<size_t>(current_cluster_index) < clusters.size()) {
            if (clusters[current_cluster_index].necessary) {
                break;
            }
        
            clusters[current_cluster_index].necessary = true;
            current_cluster_index = clusters[current_cluster_index].merged_into;
        }
    }

    void PCSTFast::mark_nodes_as_deleted(int start_node_index, int parent_node_index) {
        if (static_cast<size_t>(start_node_index) >= node_deleted.size() || node_deleted[start_node_index]) {
            return;
        }

        cluster_queue.clear();
        cluster_queue.push_back(start_node_index);

        node_deleted[start_node_index] = 1;
        size_t queue_index = 0;
        
        while (queue_index < cluster_queue.size()) {
            int current_node_index = cluster_queue[queue_index++];

            if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                for (const auto& [next_node_index, cost] : phase3_neighbors[current_node_index]) {
                    if (next_node_index == parent_node_index && start_node_index != next_node_index) {
                        continue;
                    }

                    if (static_cast<size_t>(next_node_index) < node_deleted.size() && !node_deleted[next_node_index]) {
                        node_deleted[next_node_index] = 1;
                        cluster_queue.push_back(next_node_index);
                    }
                }
            }
        }
    }

    bool PCSTFast::run(std::vector<int>* result_nodes, std::vector<int>* result_edges) {
        if (!result_nodes || !result_edges) {
            throw std::invalid_argument("Result pointers cannot be null.");
        }

        result_nodes->clear();
        result_edges->clear();

        if (root >= 0 && target_num_active_clusters > 0) {
            if (output_function) {
                snprintf(output_buffer, kOutputBufferSize, "Error: target_num_active_clusters must be 0 in the rooted case.\n");
                output_function(output_buffer);
            }

            return false;
        }

        stats = Statistics();
        const size_t num_nodes = prizes.size();
        int current_num_active_clusters = 0;
        
        for(size_t i = 0; i < clusters.size() && i < num_nodes; ++i) {
            if (clusters[i].active) {
                current_num_active_clusters++;
            }
        }

        phase1_result.clear();
        phase2_result.clear();
        phase3_result.clear();

        while (current_num_active_clusters > target_num_active_clusters) {
            if (verbosity_level >= 2 && output_function) {
                snprintf(output_buffer, kOutputBufferSize, "-----------------------------------------\n");
                output_function(output_buffer);
            }

            double next_edge_time;
            int next_edge_cluster_index;
            int next_edge_part_index = -1;
            
            get_next_edge_event(&next_edge_time, &next_edge_cluster_index, &next_edge_part_index);

            double next_cluster_time;
            int next_cluster_deactivation_index;
            
            get_next_cluster_event(&next_cluster_time, &next_cluster_deactivation_index);

            if (next_edge_time == std::numeric_limits<double>::infinity() &&
                next_cluster_time == std::numeric_limits<double>::infinity()) {
                
                if (verbosity_level >= 1 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize, "No more events possible. Terminating loop.\n");
                    output_function(output_buffer);
                }

                break;
            }

            if (verbosity_level >= 2 && output_function) {
                snprintf(output_buffer,
                        kOutputBufferSize,
                        "Next edge event: time %.9g, cluster %d, part %d\n",
                        next_edge_time,
                        next_edge_cluster_index,
                        next_edge_part_index);
                output_function(output_buffer);
                snprintf(output_buffer,
                        kOutputBufferSize,
                        "Next cluster event: time %.9g, cluster %d\n",
                        next_cluster_time,
                        next_cluster_deactivation_index);
                output_function(output_buffer);
            }

            if (next_edge_time <= next_cluster_time + eps * std::max(1.0, std::abs(next_cluster_time))) {
                stats.total_num_edge_events += 1;
                current_time = next_edge_time;

                if (next_edge_cluster_index < 0 || static_cast<size_t>(next_edge_cluster_index) >= clusters.size()) {
                    remove_next_edge_event(next_edge_cluster_index);
                
                    continue;
                }

                remove_next_edge_event(next_edge_cluster_index);

                if (next_edge_part_index < 0 || static_cast<size_t>(next_edge_part_index) >= edge_parts.size()) {
                    continue;
                }

                if (edge_parts[next_edge_part_index].deleted) {
                    stats.num_deleted_edge_events += 1;
                
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer,
                                kOutputBufferSize,
                                "Edge part %d already deleted, skipping event\n",
                                next_edge_part_index);
                        output_function(output_buffer);
                    }
                
                    continue;
                }

                int other_edge_part_index = get_other_edge_part_index(next_edge_part_index);
                int edge_index = next_edge_part_index / 2;

                if (static_cast<size_t>(edge_index) >= costs.size() || static_cast<size_t>(other_edge_part_index) >= edge_parts.size()) {
                    continue;
                }

                double sum_current_edge_part, current_finished_moat_sum;
                int current_cluster_index;
                
                get_sum_on_edge_part(next_edge_part_index, &sum_current_edge_part, &current_finished_moat_sum, &current_cluster_index);

                double sum_other_edge_part, other_finished_moat_sum;
                int other_cluster_index;
                
                get_sum_on_edge_part(other_edge_part_index, &sum_other_edge_part, &other_finished_moat_sum, &other_cluster_index);

                if (current_cluster_index == -1 || other_cluster_index == -1 ||
                    static_cast<size_t>(current_cluster_index) >= clusters.size() ||
                    static_cast<size_t>(other_cluster_index) >= clusters.size()) {
                        if (verbosity_level >= 1 && output_function) {
                            snprintf(output_buffer, kOutputBufferSize, "Warning: Could not find valid representative clusters for edge part %d or %d.\n", next_edge_part_index, other_edge_part_index);
                            output_function(output_buffer);
                        }

                        continue;
                }

                double current_edge_cost = costs[edge_index];
                double remainder = current_edge_cost - sum_current_edge_part - sum_other_edge_part;
                Cluster& current_cluster = clusters[current_cluster_index];
                Cluster& other_cluster = clusters[other_cluster_index];
                EdgePart& next_edge_part = edge_parts[next_edge_part_index];
                EdgePart& other_edge_part = edge_parts[other_edge_part_index];

                if (verbosity_level >= 2 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Edge event at time %.9g, current edge part %d (cluster %d), other edge part %d (cluster %d)\n",
                            current_time, next_edge_part_index, current_cluster_index,
                            other_edge_part_index, other_cluster_index);
                    output_function(output_buffer);
                    snprintf(output_buffer, kOutputBufferSize,
                            "Sum current part %.9g, other part %.9g, total cost %.9g, remainder %.9g\n",
                            sum_current_edge_part, sum_other_edge_part, current_edge_cost,
                            remainder);
                    output_function(output_buffer);
                }

                if (current_cluster_index == other_cluster_index) {
                    stats.num_merged_edge_events += 1;
                    
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize, "Clusters %d already merged, ignoring edge %d\n", current_cluster_index, edge_index);
                        output_function(output_buffer);
                    }
                    
                    other_edge_part.deleted = true;
                    
                    continue;
                }

                if (remainder <= eps * current_edge_cost || remainder <= eps ) {
                    stats.total_num_merge_events += 1;
                    
                    phase1_result.push_back(edge_index);
                    
                    other_edge_part.deleted = true;

                    int new_cluster_index = clusters.size();

                    clusters.emplace_back(&pairing_heap_buffer);
                    
                    Cluster& new_cluster = clusters.back();
                    Cluster& current_cluster_ref = clusters[current_cluster_index];
                    Cluster& other_cluster_ref = clusters[other_cluster_index];

                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "Merge %d and %d into %d via edge %d\n", current_cluster_index,
                                other_cluster_index, new_cluster_index, edge_index);
                        output_function(output_buffer);
                    }

                    new_cluster.prize_sum = current_cluster_ref.prize_sum + other_cluster_ref.prize_sum;
                    new_cluster.subcluster_moat_sum = current_cluster_ref.subcluster_moat_sum + other_cluster_ref.subcluster_moat_sum;
                    new_cluster.contains_root = current_cluster_ref.contains_root || other_cluster_ref.contains_root;
                    new_cluster.active = !new_cluster.contains_root;
                    new_cluster.merged_along = edge_index;
                    new_cluster.child_cluster_1 = current_cluster_index;
                    new_cluster.child_cluster_2 = other_cluster_index;

                    if (current_cluster_ref.active) {
                        current_cluster_ref.active = false;
                        current_cluster_ref.active_end_time = current_time;
                        current_cluster_ref.moat = current_cluster_ref.active_end_time - current_cluster_ref.active_start_time;
                        
                        if (current_cluster_ref.moat < 0.0) current_cluster_ref.moat = 0.0;

                        clusters_deactivation.delete_element(current_cluster_index);
                        
                        if (!current_cluster_ref.edge_parts.is_empty()) {
                            clusters_next_edge_event.delete_element(current_cluster_index);
                        }
                        
                        current_num_active_clusters -= 1;
                    } else {
                        clusters_deactivation.delete_element(current_cluster_index);
                        clusters_next_edge_event.delete_element(current_cluster_index);
                    }
                    
                    current_cluster_ref.merged_into = new_cluster_index;
                    new_cluster.subcluster_moat_sum += current_cluster_ref.moat;

                    if (other_cluster_ref.active) {
                        stats.num_active_active_merge_events += 1;
                        other_cluster_ref.active = false;
                        other_cluster_ref.active_end_time = current_time;
                        other_cluster_ref.moat = other_cluster_ref.active_end_time - other_cluster_ref.active_start_time;
                        
                        if (other_cluster_ref.moat < 0.0) other_cluster_ref.moat = 0.0;

                        clusters_deactivation.delete_element(other_cluster_index);
                        
                        if (!other_cluster_ref.edge_parts.is_empty()) {
                            clusters_next_edge_event.delete_element(other_cluster_index);
                        }
                        
                        current_num_active_clusters -= 1;
                    } else {
                        stats.num_active_inactive_merge_events += 1;
                        
                        clusters_deactivation.delete_element(other_cluster_index);
                        clusters_next_edge_event.delete_element(other_cluster_index);

                        if (!other_cluster_ref.contains_root) {
                            double edge_event_update_shift = current_time - other_cluster_ref.active_end_time;
                            
                            if (edge_event_update_shift > eps) {
                                other_cluster_ref.edge_parts.add_to_heap(edge_event_update_shift);
                            }

                            inactive_merge_events.emplace_back();
                            
                            InactiveMergeEvent& merge_event = inactive_merge_events.back();
                            merge_event.active_cluster_index = current_cluster_index;
                            merge_event.inactive_cluster_index = other_cluster_index;
                            const auto& edge = edges[edge_index];
                            int active_node_part = (next_edge_part_index % 2 == 0) ? edge.first : edge.second;
                            int inactive_node_part = (next_edge_part_index % 2 == 0) ? edge.second : edge.first;
                            merge_event.active_cluster_node = active_node_part;
                            merge_event.inactive_cluster_node = inactive_node_part;

                            if (static_cast<size_t>(edge_index) < edge_info.size()) {
                                edge_info[edge_index].inactive_merge_event = std::ssize(inactive_merge_events) - 1;
                            }
                        }
                    }

                    other_cluster_ref.merged_into = new_cluster_index;
                    new_cluster.subcluster_moat_sum += other_cluster_ref.moat;
                    new_cluster.edge_parts = PairingHeapType::meld(std::move(current_cluster_ref.edge_parts), std::move(other_cluster_ref.edge_parts));

                    if (new_cluster.active) {
                        new_cluster.active_start_time = current_time;
                        double potential_deactivation_time = new_cluster.active_start_time + new_cluster.prize_sum - new_cluster.subcluster_moat_sum;

                        if (potential_deactivation_time > current_time + eps) {
                            clusters_deactivation.insert(potential_deactivation_time, new_cluster_index);
                        } else {
                            new_cluster.active = false;
                            new_cluster.active_end_time = current_time;
                            new_cluster.moat = 0.0;
                            new_cluster.subcluster_moat_sum += new_cluster.moat;
                        }

                        if (new_cluster.active && !new_cluster.edge_parts.is_empty()) {
                            double tmp_val;
                            int tmp_index;
                    
                            (void) new_cluster.edge_parts.get_min(&tmp_val, &tmp_index);
                            clusters_next_edge_event.insert(tmp_val, new_cluster_index);
                        }

                        if (new_cluster.active) {
                            current_num_active_clusters += 1;
                        }
                    } else {
                        new_cluster.active_start_time = current_time;
                        new_cluster.active_end_time = current_time;
                        new_cluster.moat = 0.0;
                        new_cluster.subcluster_moat_sum += new_cluster.moat;
                    }

                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "Merged %d and %d into %d. New cluster active: %d. Total active: %d\n", current_cluster_index,
                                other_cluster_index, new_cluster_index, new_cluster.active, current_num_active_clusters);
                        output_function(output_buffer);
                    }

                } else {
                    stats.total_num_edge_growth_events += 1;

                    if (other_cluster.active) {
                        stats.num_active_active_edge_growth_events += 1;
                        double next_event_time = current_time + remainder * 0.5;
                        next_edge_part.next_event_val = next_event_time;
                        next_edge_part.heap_node = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);

                        if (!current_cluster.edge_parts.is_empty()) {
                            double current_min_val; 
                            int current_min_idx;
                            
                            (void) current_cluster.edge_parts.get_min(&current_min_val, &current_min_idx);
                            clusters_next_edge_event.delete_element(current_cluster_index);
                            clusters_next_edge_event.insert(current_min_val, current_cluster_index);
                        }

                        other_edge_part.next_event_val = next_event_time;
                        
                        clusters_next_edge_event.delete_element(other_cluster_index);
                        other_cluster.edge_parts.decrease_key(other_edge_part.heap_node, next_event_time);

                        if (!other_cluster.edge_parts.is_empty()) {
                            double other_min_val; 
                            int other_min_idx;
                            
                            (void) other_cluster.edge_parts.get_min(&other_min_val, &other_min_idx);
                            clusters_next_edge_event.insert(other_min_val, other_cluster_index);
                        }

                        if (verbosity_level >= 2 && output_function) {
                            snprintf(output_buffer, kOutputBufferSize, "Active-Active growth for edge %d. New event at time %.9g\n", edge_index, next_event_time);
                            output_function(output_buffer);
                        }
                    } else {
                        stats.num_active_inactive_edge_growth_events += 1;
                        double next_event_time = current_time + remainder;
                        next_edge_part.next_event_val = next_event_time;
                        next_edge_part.heap_node = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);

                        if (!current_cluster.edge_parts.is_empty()) {
                            double current_min_val; 
                            int current_min_idx;
                        
                            (void) current_cluster.edge_parts.get_min(&current_min_val, &current_min_idx);
                        
                            clusters_next_edge_event.delete_element(current_cluster_index);
                            clusters_next_edge_event.insert(current_min_val, current_cluster_index);
                        }

                        other_edge_part.next_event_val = other_cluster.active_end_time;
                        
                        other_cluster.edge_parts.decrease_key(other_edge_part.heap_node, other_cluster.active_end_time);

                        if (verbosity_level >= 2 && output_function) {
                            snprintf(output_buffer, kOutputBufferSize,
                                    "Active-Inactive growth edge %d. New event for active side at %.9g. Inactive part reset.\n",
                                    edge_index, next_event_time);
                            output_function(output_buffer);
                        }
                    }
                }
            } else {
                stats.num_cluster_events += 1;
                current_time = next_cluster_time;

                if (next_cluster_deactivation_index < 0 || static_cast<size_t>(next_cluster_deactivation_index) >= clusters.size()) {
                    remove_next_cluster_event();
                
                    continue;
                }

                Cluster& cluster_to_deactivate = clusters[next_cluster_deactivation_index];
                
                remove_next_cluster_event();

                if (!cluster_to_deactivate.active) {
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize, "Cluster %d already inactive, skipping deactivation event.\n", next_cluster_deactivation_index);
                        output_function(output_buffer);
                    }
                
                    continue;
                }

                cluster_to_deactivate.active = false;
                cluster_to_deactivate.active_end_time = current_time;
                cluster_to_deactivate.moat = cluster_to_deactivate.active_end_time - cluster_to_deactivate.active_start_time;
                
                if (cluster_to_deactivate.moat < 0.0) cluster_to_deactivate.moat = 0.0;

                if (!cluster_to_deactivate.edge_parts.is_empty()) {
                    clusters_next_edge_event.delete_element(next_cluster_deactivation_index);
                }

                current_num_active_clusters -= 1;

                if (verbosity_level >= 2 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Cluster deactivation: cluster %d at time %.9g (moat size %.9g). Total active: %d\n",
                            next_cluster_deactivation_index, current_time, cluster_to_deactivate.moat, current_num_active_clusters);
                    output_function(output_buffer);
                }
            }
        }

        if (verbosity_level >= 1 && output_function) {
            snprintf(output_buffer, kOutputBufferSize,
                    "Finished GW main loop: final event time %.9g, number of edge events %lld, active clusters %d\n",
                    current_time, stats.total_num_edge_events, current_num_active_clusters);
            output_function(output_buffer);
        }

        node_good.assign(num_nodes, 0);

        if (root >= 0) {
            int root_cluster_final_index = root;
            
            if (static_cast<size_t>(root) < clusters.size()) {
                root_cluster_final_index = root;
                
                while(static_cast<size_t>(root_cluster_final_index) < clusters.size() && clusters[root_cluster_final_index].merged_into != -1) {
                    root_cluster_final_index = clusters[root_cluster_final_index].merged_into;
                }

                if (static_cast<size_t>(root_cluster_final_index) < clusters.size() && clusters[root_cluster_final_index].contains_root) {
                        mark_nodes_as_good(root_cluster_final_index);
                } else if (verbosity_level >= 0 && output_function){
                        snprintf(output_buffer, kOutputBufferSize, "Warning: Could not find final cluster containing root %d.\n", root);
                        output_function(output_buffer);
                }
            } else if (verbosity_level >= 0 && output_function) {
                snprintf(output_buffer, kOutputBufferSize, "Warning: Root index %d out of initial cluster bounds.\n", root);
                output_function(output_buffer);
            }
        } else {
            for (size_t ii = 0; ii < clusters.size(); ++ii) {
                if (clusters[ii].merged_into == -1 && clusters[ii].active) {
                    mark_nodes_as_good(ii);
                }
            }
        }

        if (pruning == PruningMethod::kNoPruning) {
            build_phase1_node_set(phase1_result, result_nodes);
        
            *result_edges = phase1_result;
        
            return true;
        }

        if (verbosity_level >= 2 && output_function) {
            snprintf(output_buffer, kOutputBufferSize, "------------------------------------------\n");
            output_function(output_buffer);
            snprintf(output_buffer, kOutputBufferSize, "Starting pruning phase...\n");
            output_function(output_buffer);
        }

        phase2_result.clear();
        phase2_result.reserve(phase1_result.size());
        
        for (int edge_idx : phase1_result) {
            if (static_cast<size_t>(edge_idx) < edges.size()) {
                const auto& edge = edges[edge_idx];
                
                if (static_cast<size_t>(edge.first) < node_good.size() && node_good[edge.first] &&
                    static_cast<size_t>(edge.second) < node_good.size() && node_good[edge.second]) {
                        phase2_result.push_back(edge_idx);
                }
            }
        }

        if (pruning == PruningMethod::kSimplePruning) {
            build_phase2_node_set(result_nodes);
        
            *result_edges = phase2_result;
        
            return true;
        }

        phase3_result.clear();
        phase3_result.reserve(phase2_result.size());

        if (phase3_neighbors.size() != num_nodes) phase3_neighbors.resize(num_nodes);
        
        for (auto& neighbors : phase3_neighbors) { neighbors.clear(); }

        for (int current_edge_index : phase2_result) {
            if (static_cast<size_t>(current_edge_index) < edges.size() && static_cast<size_t>(current_edge_index) < costs.size()) {
                const auto& edge = edges[current_edge_index];
                double current_cost = costs[current_edge_index];
        
                if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                    phase3_neighbors[edge.first].emplace_back(edge.second, current_cost);
                    phase3_neighbors[edge.second].emplace_back(edge.first, current_cost);
                }
            }
        }

        if (node_deleted.size() != num_nodes) node_deleted.resize(num_nodes);
        
        std::fill(node_deleted.begin(), node_deleted.end(), 0);

        if (pruning == PruningMethod::kGWPruning) {
            if (verbosity_level >= 2 && output_function) {
                snprintf(output_buffer, kOutputBufferSize, "Starting GW pruning, phase 2 result size: %zu\n", phase2_result.size());
                output_function(output_buffer);
            }

            for(auto& cluster : clusters) { cluster.necessary = false; }

            for (int ii = std::ssize(phase2_result) - 1; ii >= 0; --ii) {
                int current_edge_index = phase2_result[ii];

                if(static_cast<size_t>(current_edge_index) >= edges.size() || static_cast<size_t>(current_edge_index) >= edge_info.size()) continue;

                const auto& edge = edges[current_edge_index];
                int uu = edge.first;
                int vv = edge.second;
                bool u_deleted = (static_cast<size_t>(uu) >= node_deleted.size() || node_deleted[uu]);
                bool v_deleted = (static_cast<size_t>(vv) >= node_deleted.size() || node_deleted[vv]);

                if (u_deleted && v_deleted) {
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "GW: Edge %d (%d, %d) endpoints deleted, skipping.\n", current_edge_index, uu, vv);
                        output_function(output_buffer);
                    }
                    
                    continue;
                }

                int inactive_merge_idx = edge_info[current_edge_index].inactive_merge_event;

                if (inactive_merge_idx < 0) {
                    if (!u_deleted) mark_clusters_as_necessary(uu);
                    if (!v_deleted) mark_clusters_as_necessary(vv);
                    
                    phase3_result.push_back(current_edge_index);
                    
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "GW: Active-Active edge %d (%d, %d), keeping.\n", current_edge_index, uu, vv);
                        output_function(output_buffer);
                    }
                } else {
                    if (static_cast<size_t>(inactive_merge_idx) >= inactive_merge_events.size()) continue;

                    const InactiveMergeEvent& current_merge_event = inactive_merge_events[inactive_merge_idx];
                    int inactive_cluster_start_node = current_merge_event.inactive_cluster_node;
                    int active_cluster_start_node = current_merge_event.active_cluster_node;
                    int inactive_cluster_rep_at_merge = current_merge_event.inactive_cluster_index;
                    bool is_inactive_side_necessary = (static_cast<size_t>(inactive_cluster_rep_at_merge) < clusters.size() &&
                                                    clusters[inactive_cluster_rep_at_merge].necessary);

                    if (is_inactive_side_necessary) {
                        phase3_result.push_back(current_edge_index);
                        mark_clusters_as_necessary(inactive_cluster_start_node);
                        mark_clusters_as_necessary(active_cluster_start_node);
            
                        if (verbosity_level >= 2 && output_function) {
                            snprintf(output_buffer, kOutputBufferSize,
                                    "GW: Inactive side cluster %d necessary for edge %d (%d, %d), keeping.\n",
                                    inactive_cluster_rep_at_merge, current_edge_index, uu, vv);
                            output_function(output_buffer);
                        }
                    } else {
                        mark_nodes_as_deleted(inactive_cluster_start_node, active_cluster_start_node);
            
                        if (verbosity_level >= 2 && output_function) {
                            snprintf(output_buffer, kOutputBufferSize,
                                    "GW: Inactive side cluster %d not necessary for edge %d (%d, %d), discarding inactive side starting at node %d.\n",
                                    inactive_cluster_rep_at_merge, current_edge_index, uu, vv, inactive_cluster_start_node);
                            output_function(output_buffer);
                        }
                    }
                }
            }

            std::reverse(phase3_result.begin(), phase3_result.end());
            build_phase3_node_set(result_nodes);
            
            *result_edges = phase3_result;
            
            return true;

        } else if (pruning == PruningMethod::kStrongPruning) {
            if (verbosity_level >= 2 && output_function) {
                snprintf(output_buffer, kOutputBufferSize,
                        "Starting Strong pruning, phase 2 result size: %zu\n", phase2_result.size());
                output_function(output_buffer);
            }

            if (final_component_label.size() != num_nodes) final_component_label.resize(num_nodes);
            
            std::fill(final_component_label.begin(), final_component_label.end(), -1);
            final_components.clear();
            
            root_component_index = -1;

            if (strong_pruning_parent.size() != num_nodes) strong_pruning_parent.resize(num_nodes);
            
            std::fill(strong_pruning_parent.begin(), strong_pruning_parent.end(), std::make_pair(-1, 0.0));
            
            if (strong_pruning_payoff.size() != num_nodes) strong_pruning_payoff.resize(num_nodes);
            
            std::fill(strong_pruning_payoff.begin(), strong_pruning_payoff.end(), 0.0);

            for (int ii = 0; ii < std::ssize(prizes); ++ii) {
                if (static_cast<size_t>(ii) < phase3_neighbors.size() && !phase3_neighbors[ii].empty() && final_component_label[ii] == -1) {
                    final_components.emplace_back();
                    label_final_component(ii, std::ssize(final_components) - 1);
                } else if (static_cast<size_t>(ii) < phase3_neighbors.size() && phase3_neighbors[ii].empty() && static_cast<size_t>(ii) < node_good.size() && node_good[ii]) {
                    if (final_component_label[ii] == -1) {
                        final_components.emplace_back();
                        label_final_component(ii, std::ssize(final_components) - 1);
                    }
                }
            }

            for (size_t ii = 0; ii < final_components.size(); ++ii) {
                if (final_components[ii].empty()) continue;

                if (verbosity_level >= 2 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize, "Strong pruning on final component %zu (size %zu):\n",
                            ii, final_components[ii].size());
                    output_function(output_buffer);
                }

                int component_root_node;
            
                if (static_cast<int>(ii) == root_component_index) {
                    component_root_node = root;
            
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize, "Component contains root %d, pruning from there.\n", root);
                        output_function(output_buffer);
                    }
                } else {
                    component_root_node = find_best_component_root(ii);
            
                    if (verbosity_level >= 2 && output_function) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "Best start node for component %zu: %d, pruning from there.\n",
                                ii, component_root_node);
                        output_function(output_buffer);
                    }
                }

                if (component_root_node != -1 && static_cast<size_t>(component_root_node) < num_nodes) {
                    strong_pruning_from(component_root_node, true);
                } else if (verbosity_level >= 1 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize, "Warning: Could not determine valid root for component %zu, skipping pruning.\n", ii);
                    output_function(output_buffer);
                }
            }

            phase3_result.clear();

            for (int current_edge_index : phase2_result) {
                if (static_cast<size_t>(current_edge_index) >= edges.size()) continue;

                const auto& edge = edges[current_edge_index];
                int uu = edge.first;
                int vv = edge.second;
                bool u_deleted = (static_cast<size_t>(uu) >= node_deleted.size() || node_deleted[uu]);
                bool v_deleted = (static_cast<size_t>(vv) >= node_deleted.size() || node_deleted[vv]);

                if (!u_deleted && !v_deleted) {
                    phase3_result.push_back(current_edge_index);
                } else if (verbosity_level >= 2 && output_function) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Strong: Not keeping edge %d (%d, %d) due to deleted endpoint(s).\n",
                            current_edge_index, uu, vv);
                    output_function(output_buffer);
                }
            }

            build_phase3_node_set(result_nodes);
        
            *result_edges = phase3_result;
        
            return true;
        }

        if (output_function) {
            snprintf(output_buffer, kOutputBufferSize, "Error: unknown pruning scheme specified.\n");
            output_function(output_buffer);
        }
        
        return false;
    }

    void PCSTFast::label_final_component(int start_node_index, int new_component_index) {
        cluster_queue.clear();

        if (static_cast<size_t>(start_node_index) >= final_component_label.size() ||
            static_cast<size_t>(new_component_index) >= final_components.size()) return;

        if (final_component_label[start_node_index] != -1) return;

        cluster_queue.push_back(start_node_index);
        
        final_component_label[start_node_index] = new_component_index;
        
        final_components[new_component_index].push_back(start_node_index);

        if (start_node_index == root) {
            root_component_index = new_component_index;
        }

        size_t queue_next = 0;

        while (queue_next < cluster_queue.size()) {
            int current_node_index = cluster_queue[queue_next++];

            if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                for (const auto& [next_node_index, cost] : phase3_neighbors[current_node_index]) {
                    if (static_cast<size_t>(next_node_index) < final_component_label.size() &&
                        final_component_label[next_node_index] == -1)
                    {
                        final_component_label[next_node_index] = new_component_index;
                        
                        final_components[new_component_index].push_back(next_node_index);
                        cluster_queue.push_back(next_node_index);
                        
                        if (next_node_index == root) {
                            root_component_index = new_component_index;
                        }
                    }
                }
            }
        }
    }

    void PCSTFast::strong_pruning_from(int start_node_index, bool mark_as_deleted) {
        stack.clear();
        
        const size_t num_nodes = prizes.size();

        if (static_cast<size_t>(start_node_index) >= num_nodes ||
            static_cast<size_t>(start_node_index) >= strong_pruning_parent.size() ) return;

        strong_pruning_parent[start_node_index] = std::make_pair(-1, 0.0);
        
        stack.emplace_back(true, start_node_index);

        while (!stack.empty()) {
            bool is_pre_order = stack.back().first;
            int current_node_index = stack.back().second;
            
            stack.pop_back();

            if (static_cast<size_t>(current_node_index) >= num_nodes) continue;

            if (is_pre_order) {
                stack.emplace_back(false, current_node_index);
                strong_pruning_payoff[current_node_index] = prizes[current_node_index];

                if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                    for (const auto& [next_node_index, next_cost] : phase3_neighbors[current_node_index]) {
                        if (next_node_index == strong_pruning_parent[current_node_index].first) {
                            continue;
                        }
                        if (static_cast<size_t>(next_node_index) >= num_nodes) continue;
                        
                        strong_pruning_parent[next_node_index] = std::make_pair(current_node_index, next_cost);
                        
                        stack.emplace_back(true, next_node_index);
                    }
                }
            } else {
                if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                    for (const auto& [next_node_index, cost] : phase3_neighbors[current_node_index]) {
                        if (static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                            strong_pruning_parent[next_node_index].first == current_node_index)
                        {
                            double child_edge_cost = strong_pruning_parent[next_node_index].second;
                            double child_subtree_payoff = strong_pruning_payoff[next_node_index] - child_edge_cost;

                            if (child_subtree_payoff <= 0.0) {
                                if (mark_as_deleted) {
                                    if (verbosity_level >= 2 && output_function) {
                                        snprintf(output_buffer, kOutputBufferSize,
                                                "Strong: Subtree at %d non-positive (%.9g), pruning.\n",
                                                next_node_index, child_subtree_payoff);
                                        output_function(output_buffer);
                                    }
                                    
                                    mark_nodes_as_deleted(next_node_index, current_node_index);
                                }
                            } else {
                                strong_pruning_payoff[current_node_index] += child_subtree_payoff;
                            }
                        }
                    }
                }
            }
        }
    }

    int PCSTFast::find_best_component_root(int component_index) {
        if (static_cast<size_t>(component_index) >= final_components.size() || final_components[component_index].empty()) {
            return -1;
        }

        const size_t num_nodes = prizes.size();
        const auto& component_nodes = final_components[component_index];
        int initial_root = component_nodes[0];

        if (static_cast<size_t>(initial_root) >= num_nodes) return -1;

        for (int node_idx : component_nodes) {
            if (static_cast<size_t>(node_idx) < strong_pruning_parent.size()) {
                strong_pruning_parent[node_idx] = std::make_pair(-1, 0.0);
            }
            if (static_cast<size_t>(node_idx) < strong_pruning_payoff.size()) {
                strong_pruning_payoff[node_idx] = 0.0;
            }
        }

        strong_pruning_from(initial_root, false);

        if (static_cast<size_t>(initial_root) >= strong_pruning_payoff.size()) return -1;

        int current_best_root = initial_root;
        double current_best_value = strong_pruning_payoff[initial_root];

        stack2.clear();

        if (static_cast<size_t>(initial_root) < phase3_neighbors.size()) {
            for (const auto& [next_node_index, cost] : phase3_neighbors[initial_root]) {
                if (static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                    strong_pruning_parent[next_node_index].first == initial_root)
                {
                    stack2.push_back(next_node_index);
                }
            }
        }

        while (!stack2.empty()) {
            int current_node_index = stack2.back();
            
            stack2.pop_back();

            if (static_cast<size_t>(current_node_index) >= num_nodes ||
                static_cast<size_t>(current_node_index) >= strong_pruning_parent.size() ||
                static_cast<size_t>(current_node_index) >= strong_pruning_payoff.size()) continue;

            int parent_index = strong_pruning_parent[current_node_index].first;
            double parent_edge_cost = strong_pruning_parent[current_node_index].second;

            if (parent_index == -1 || static_cast<size_t>(parent_index) >= num_nodes ||
                static_cast<size_t>(parent_index) >= strong_pruning_payoff.size()) continue;

            double payoff_subtree_current = strong_pruning_payoff[current_node_index];
            double payoff_parent_rooted = strong_pruning_payoff[parent_index];
            double current_contribution_to_parent = 0.0;
            
            if (payoff_subtree_current > parent_edge_cost) {
                current_contribution_to_parent = payoff_subtree_current - parent_edge_cost;
            }
            
            double payoff_parent_without_current = payoff_parent_rooted - current_contribution_to_parent;
            double parent_contribution_to_current = 0.0;
            
            if (payoff_parent_without_current > parent_edge_cost) {
                parent_contribution_to_current = payoff_parent_without_current - parent_edge_cost;
            }
            
            double current_node_total_payoff_if_root = payoff_subtree_current + parent_contribution_to_current;

            if (current_node_total_payoff_if_root > current_best_value) {
                current_best_root = current_node_index;
                current_best_value = current_node_total_payoff_if_root;
            }

            strong_pruning_payoff[current_node_index] = current_node_total_payoff_if_root;

            if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                for (const auto& [next_node_index, cost] : phase3_neighbors[current_node_index]) {
                    if (static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                        strong_pruning_parent[next_node_index].first == current_node_index)
                    {
                        stack2.push_back(next_node_index);
                    }
                }
            }
        }
        return current_best_root;
    }

    void PCSTFast::build_phase3_node_set(std::vector<int>* node_set) {
        if (!node_set) return;
        
        node_set->clear();
        node_set->reserve(prizes.size());

        for (int ii = 0; ii < std::ssize(prizes); ++ii) {
            bool is_good = (static_cast<size_t>(ii) < node_good.size() && node_good[ii]);
            bool is_deleted = (static_cast<size_t>(ii) < node_deleted.size() && node_deleted[ii]);
            
            if (is_good && !is_deleted) {
                    node_set->push_back(ii);
            }
        }
    }

    void PCSTFast::build_phase2_node_set(std::vector<int>* node_set) {
        if (!node_set) return;
        
        node_set->clear();
        node_set->reserve(prizes.size());

        for (int ii = 0; ii < std::ssize(prizes); ++ii) {
            bool is_good = (static_cast<size_t>(ii) < node_good.size() && node_good[ii]);
            
            if (is_good) {
                node_set->push_back(ii);
            }
        }
    }

    void PCSTFast::build_phase1_node_set(const std::vector<int>& edge_set, std::vector<int>* node_set) {
        if (!node_set) return;
        
        node_set->clear();
        node_set->reserve(prizes.size());

        if (build_phase1_included_nodes.size() != prizes.size()) {
            build_phase1_included_nodes.assign(prizes.size(), 0);
        } else {
            std::fill(build_phase1_included_nodes.begin(), build_phase1_included_nodes.end(), 0);
        }

        for (int edge_idx : edge_set) {
            if (static_cast<size_t>(edge_idx) >= edges.size()) continue;

            const auto& edge = edges[edge_idx];
            int uu = edge.first;
            int vv = edge.second;

            if (static_cast<size_t>(uu) < build_phase1_included_nodes.size() && !build_phase1_included_nodes[uu]) {
                build_phase1_included_nodes[uu] = 1;
                
                node_set->push_back(uu);
            }
            
            if (static_cast<size_t>(vv) < build_phase1_included_nodes.size() && !build_phase1_included_nodes[vv]) {
                build_phase1_included_nodes[vv] = 1;
                
                node_set->push_back(vv);
            }
        }

        for (int ii = 0; ii < std::ssize(prizes); ++ii) {
            bool is_good = (static_cast<size_t>(ii) < node_good.size() && node_good[ii]);
            bool is_included = (static_cast<size_t>(ii) < build_phase1_included_nodes.size() && build_phase1_included_nodes[ii]);
            
            if (is_good && !is_included) {
                    node_set->push_back(ii);
            }
        }
    }

    void PCSTFast::get_statistics(Statistics* s) const {
        if(s) {
            *s = stats;
        }
    }
}
