#include <algorithm>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>
#include <cctype>

#include "pcst_fast.h"

using cluster_approx::PCSTFast;
using std::make_pair;
using std::vector;

PCSTFast::Statistics::Statistics() : total_num_edge_events(0),
                                    num_deleted_edge_events(0),
                                    num_merged_edge_events(0),
                                    total_num_merge_events(0),
                                    num_active_active_merge_events(0),
                                    num_active_inactive_merge_events(0),
                                    total_num_edge_growth_events(0),
                                    num_active_active_edge_growth_events(0),
                                    num_active_inactive_edge_growth_events(0),
                                    num_cluster_events(0) { };

PCSTFast::PruningMethod PCSTFast::parse_pruning_method(const std::string& input) {
    PruningMethod result = kUnknownPruning;
  
    if (input.length() == 4 && tolower(input[0]) == 'n' && tolower(input[1]) == 'o' && tolower(input[2]) == 'n' && tolower(input[3]) == 'e') {
        result = kNoPruning;
    } else if (input.length() == 6 && tolower(input[0]) == 's' && tolower(input[1]) == 'i' && tolower(input[2]) == 'm' && tolower(input[3]) == 'p' && tolower(input[4]) == 'l' && tolower(input[5]) == 'e') {
        result = kSimplePruning;
    } else if (input.length() == 2 && tolower(input[0]) == 'g' && tolower(input[1]) == 'w') {
        result = kGWPruning;
    } else if (input.length() == 6 && tolower(input[0]) == 's' && tolower(input[1]) == 't' && tolower(input[2]) == 'r' && tolower(input[3]) == 'o' && tolower(input[4]) == 'n' && tolower(input[5]) == 'g') {
        result = kStrongPruning;
    }
  
    return result;
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
      target_num_active_clusters(target_num_active_clusters_),
      pruning(pruning_), verbosity_level(verbosity_level_),
      output_function(output_function_) {

    const size_t num_nodes = prizes.size();
    const size_t num_edges = edges.size();

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

    if (pruning == kStrongPruning || pruning == kGWPruning) {
         phase3_neighbors.resize(num_nodes);
    }
    if (pruning == kStrongPruning) {
        final_component_label.resize(num_nodes);
        strong_pruning_parent.resize(num_nodes);
        strong_pruning_payoff.resize(num_nodes);
        final_components.reserve(num_nodes);
        stack.reserve(num_nodes);
        stack2.reserve(num_nodes);
    }

    for (size_t ii = 0; ii < edge_info.size(); ++ii) {
        edge_info[ii].inactive_merge_event = -1;
    }

    current_time = 0.0;
    eps = 1e-6;

    for (int ii = 0; ii < static_cast<int>(num_nodes); ++ii) {
        if (prizes[ii] < 0.0) {
            throw std::invalid_argument("Prize negative.");
        }
      
        clusters.push_back(Cluster(&pairing_heap_buffer)); // Keep original push_back for Cluster
        
        clusters[ii].active = (ii != root);
        clusters[ii].active_start_time = 0.0;
        clusters[ii].active_end_time = (ii == root) ? 0.0 : -1.0;
        clusters[ii].prize_sum = prizes[ii];
        clusters[ii].contains_root = (ii == root);

        if (clusters[ii].active) {
            clusters_deactivation.insert(prizes[ii], ii);
        }
    }

    for (int ii = 0; ii < static_cast<int>(num_edges); ++ii) {
        const auto& edge = edges[ii];
        int uu = edge.first;
        int vv = edge.second;
    
        if (uu < 0 || vv < 0) {
            throw std::invalid_argument("Edge endpoint negative.");
        }
        
        if (uu >= static_cast<int>(prizes.size()) || vv >= static_cast<int>(prizes.size())) {
            throw std::invalid_argument("Edge endpoint out of range (too large).");
        }
    
        double cost = costs[ii];
    
        if (cost < 0.0) {
            throw std::invalid_argument("Edge cost negative.");
        }
    
        EdgePart& uu_part = edge_parts[2 * ii];
        EdgePart& vv_part = edge_parts[2 * ii + 1];
        Cluster& uu_cluster = clusters[uu];
        Cluster& vv_cluster = clusters[vv];
        uu_part.deleted = false;
        vv_part.deleted = false;

        if (uu_cluster.active && vv_cluster.active) {
            double event_time = cost * 0.5;
            uu_part.next_event_val = event_time;
            vv_part.next_event_val = event_time;
        } else if (uu_cluster.active) {
            uu_part.next_event_val = cost;
            vv_part.next_event_val = 0.0;
        } else if (vv_cluster.active) {
            uu_part.next_event_val = 0.0;
            vv_part.next_event_val = cost;
        } else {
            uu_part.next_event_val = 0.0;
            vv_part.next_event_val = 0.0;
        }

        uu_part.heap_node = uu_cluster.edge_parts.insert(uu_part.next_event_val, 2 * ii);
        vv_part.heap_node = vv_cluster.edge_parts.insert(vv_part.next_event_val, 2 * ii + 1);
    }

    for (int ii = 0; ii < static_cast<int>(num_nodes); ++ii) {
        if (clusters[ii].active) {
            if (!clusters[ii].edge_parts.is_empty()) {
                double val;
                int edge_part;
        
                clusters[ii].edge_parts.get_min(&val, &edge_part);
                clusters_next_edge_event.insert(val, ii);
            }
        }
    }
}

PCSTFast::~PCSTFast() {
    for (size_t ii = 0; ii < clusters.size(); ++ii) {
        clusters[ii].edge_parts.release_memory();
    }
}

void PCSTFast::get_next_edge_event(double* next_time, int* next_cluster_index, int* next_edge_part_index) {
    if (clusters_next_edge_event.is_empty()) {
        *next_time = std::numeric_limits<double>::infinity();
        *next_cluster_index = -1;
        *next_edge_part_index = -1;
        
        return;
    }

    clusters_next_edge_event.get_min(next_time, next_cluster_index);
  
    if (static_cast<size_t>(*next_cluster_index) < clusters.size()) {
        clusters[*next_cluster_index].edge_parts.get_min(next_time, next_edge_part_index);
    } else {
        *next_time = std::numeric_limits<double>::infinity();
        *next_edge_part_index = -1;
    }
}

void PCSTFast::remove_next_edge_event(int next_cluster_index) {
    if (static_cast<size_t>(next_cluster_index) >= clusters.size()) {
        clusters_next_edge_event.delete_element(next_cluster_index);
        
        return;
    }

    clusters_next_edge_event.delete_element(next_cluster_index);
  
    double tmp_value;
    int tmp_edge_part;
    
    clusters[next_cluster_index].edge_parts.delete_min(&tmp_value, &tmp_edge_part);

    if (!clusters[next_cluster_index].edge_parts.is_empty()) {
        clusters[next_cluster_index].edge_parts.get_min(&tmp_value, &tmp_edge_part);
        clusters_next_edge_event.insert(tmp_value, next_cluster_index);
    }
}

void PCSTFast::get_next_cluster_event(double* next_time, int* next_cluster_index) {
    if (clusters_deactivation.is_empty()) {
        *next_time = std::numeric_limits<double>::infinity();
        *next_cluster_index = -1;
        
        return;
    }

    clusters_deactivation.get_min(next_time, next_cluster_index);
}

void PCSTFast::remove_next_cluster_event() {
    double tmp_value;
    int tmp_cluster;
    
    clusters_deactivation.delete_min(&tmp_value, &tmp_cluster);
}

void PCSTFast::get_sum_on_edge_part(int edge_part_index, double* total_sum, double* finished_moat_sum, int* current_cluster_index) {
    const auto& edge = edges[edge_part_index / 2];
    int endpoint = (edge_part_index % 2 == 0) ? edge.first : edge.second;
    *total_sum = 0.0;
    *current_cluster_index = endpoint;
  
    path_compression_visited.clear();

    int current_node = endpoint;
  
    while (static_cast<size_t>(current_node) < clusters.size() && clusters[current_node].merged_into != -1) {
        path_compression_visited.emplace_back(current_node, *total_sum);
        
        int next_node;
    
        if (clusters[current_node].skip_up >= 0) {
            *total_sum += clusters[current_node].skip_up_sum;
            next_node = clusters[current_node].skip_up;
        } else {
            *total_sum += clusters[current_node].moat;
            next_node = clusters[current_node].merged_into;
        }
        
        current_node = next_node;
    }

    if (static_cast<size_t>(current_node) >= clusters.size()) {
        *current_cluster_index = -1;
        *total_sum = 0.0;
        *finished_moat_sum = 0.0;
      
        return;
    }
    
    *current_cluster_index = current_node;
    
    for (const auto& visited_pair : path_compression_visited) {
        int visited_cluster_index = visited_pair.first;
        double visited_sum = visited_pair.second;
    
        if (static_cast<size_t>(visited_cluster_index) < clusters.size()) {
            clusters[visited_cluster_index].skip_up = *current_cluster_index;
            clusters[visited_cluster_index].skip_up_sum = *total_sum - visited_sum;
        }
    }

    if (static_cast<size_t>(*current_cluster_index) < clusters.size()) {
        if (clusters[*current_cluster_index].active) {
            *finished_moat_sum = *total_sum;
            *total_sum += current_time - clusters[*current_cluster_index].active_start_time;
        } else {
            *total_sum += clusters[*current_cluster_index].moat;
            *finished_moat_sum = *total_sum;
        }
    } else {
       *total_sum = 0.0;
       *finished_moat_sum = 0.0;
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

    if (clusters[current_cluster_index].merged_along >= 0) {
        if(clusters[current_cluster_index].child_cluster_1 != -1) cluster_queue.push_back(clusters[current_cluster_index].child_cluster_1);
        if(clusters[current_cluster_index].child_cluster_2 != -1) cluster_queue.push_back(clusters[current_cluster_index].child_cluster_2);
    } else {
       if (static_cast<size_t>(current_cluster_index) < node_good.size()) {
            node_good[current_cluster_index] = 1;
       }
    }
  }
}

void PCSTFast::mark_clusters_as_necessary(int start_cluster_index) {
    int current_cluster_index = start_cluster_index;
    
    while (current_cluster_index != -1 && static_cast<size_t>(current_cluster_index) < clusters.size() && !clusters[current_cluster_index].necessary) {
        clusters[current_cluster_index].necessary = true;
        current_cluster_index = clusters[current_cluster_index].merged_into;
    }
}

void PCSTFast::mark_nodes_as_deleted(int start_node_index, int parent_node_index) {
    if (static_cast<size_t>(start_node_index) >= node_deleted.size() || node_deleted[start_node_index]) {
        return;
    }
    
    node_deleted[start_node_index] = 1;
    
    cluster_queue.clear();
    cluster_queue.push_back(start_node_index);
    
    size_t queue_index = 0;
    
    while (queue_index < cluster_queue.size()) {
        int current_node_index = cluster_queue[queue_index++];
        
        if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
            for (const auto& neighbor_pair : phase3_neighbors[current_node_index]) {
                int next_node_index = neighbor_pair.first;
                
                if (next_node_index == parent_node_index) {
                    continue;
                }
            
                if (static_cast<size_t>(next_node_index) >= node_deleted.size()) continue;
                
                if (node_deleted[next_node_index]) {
                    continue;
                }
            
                node_deleted[next_node_index] = 1;
                cluster_queue.push_back(next_node_index);
            }
        }
    }
}

bool PCSTFast::run(std::vector<int>* result_nodes, std::vector<int>* result_edges) {
    result_nodes->clear();
    result_edges->clear();

    if (root >= 0 && target_num_active_clusters > 0) {
        snprintf(output_buffer, kOutputBufferSize, "Error: target_num_active_clusters must be 0 in the rooted case.\n");
        output_function(output_buffer);
        
        return false;
    }

    stats = Statistics();

    size_t num_nodes = prizes.size();
    int num_active_clusters = static_cast<int>(num_nodes);
    
    if (root >= 0 && static_cast<size_t>(root) < num_nodes) {
        num_active_clusters -= 1;
    }

    phase1_result.clear();
    phase2_result.clear();
    phase3_result.clear();

    while (num_active_clusters > target_num_active_clusters) {
        if (verbosity_level >= 2) {
            snprintf(output_buffer, kOutputBufferSize, "-----------------------------------------\n");
            output_function(output_buffer);
        }

        double next_edge_time;
        int next_edge_cluster_index;
        int next_edge_part_index = -1;
    
        get_next_edge_event(&next_edge_time, &next_edge_cluster_index, &next_edge_part_index);
        
        double next_cluster_time;
        int next_cluster_index;
    
        get_next_cluster_event(&next_cluster_time, &next_cluster_index);

        if (verbosity_level >= 2) {
            snprintf(output_buffer, 
                    kOutputBufferSize,
                    "Next edge event: time %e, cluster %d, part %d\n", 
                    next_edge_time, 
                    next_edge_cluster_index, 
                    next_edge_part_index);
            output_function(output_buffer);
            snprintf(output_buffer,
                    kOutputBufferSize,
                    "Next cluster event: time %e, cluster %d\n", 
                    next_cluster_time, 
                    next_cluster_index);
            output_function(output_buffer);
        }

        if (next_edge_time < next_cluster_time) {
            stats.total_num_edge_events += 1;
            current_time = next_edge_time;
            
            if (static_cast<size_t>(next_edge_cluster_index) >= clusters.size()) continue;
      
            remove_next_edge_event(next_edge_cluster_index);

            if (next_edge_part_index < 0 || static_cast<size_t>(next_edge_part_index) >= edge_parts.size()) continue;

            if (edge_parts[next_edge_part_index].deleted) {
                stats.num_deleted_edge_events += 1;
        
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, 
                            kOutputBufferSize, 
                            "Edge part %d already deleted, nothing to do\n", 
                            next_edge_part_index);
                    output_function(output_buffer);
                }
            
                continue;
            }

            int other_edge_part_index = get_other_edge_part_index(next_edge_part_index);
            int edge_index = next_edge_part_index / 2;
      
            if (static_cast<size_t>(edge_index) >= costs.size() || static_cast<size_t>(other_edge_part_index) >= edge_parts.size()) continue;
            
            double current_edge_cost = costs[edge_index];
            double sum_current_edge_part;
            int current_cluster_index;
            double current_finished_moat_sum;
      
            get_sum_on_edge_part(next_edge_part_index,
                                &sum_current_edge_part,
                                &current_finished_moat_sum,
                                &current_cluster_index);
      
            double sum_other_edge_part;
            int other_cluster_index;
            double other_finished_moat_sum;
            
            get_sum_on_edge_part(other_edge_part_index,
                                &sum_other_edge_part,
                                &other_finished_moat_sum,
                                &other_cluster_index);

            if (current_cluster_index == -1 || other_cluster_index == -1 ||
                static_cast<size_t>(current_cluster_index) >= clusters.size() ||
                static_cast<size_t>(other_cluster_index) >= clusters.size()) {
                    continue;
            }

            double remainder = current_edge_cost - sum_current_edge_part - sum_other_edge_part;

            Cluster& current_cluster = clusters[current_cluster_index];
            Cluster& other_cluster = clusters[other_cluster_index];
            EdgePart& next_edge_part = edge_parts[next_edge_part_index];
            EdgePart& other_edge_part = edge_parts[other_edge_part_index];

            if (verbosity_level >= 2) {
                snprintf(output_buffer, kOutputBufferSize,
                        "Edge event at time %e, current edge part %d (cluster %d), "
                        "other edge part %d (cluster %d)\n",
                        current_time, next_edge_part_index, current_cluster_index,
                        other_edge_part_index, other_cluster_index);
                output_function(output_buffer);
                snprintf(output_buffer, kOutputBufferSize,
                        "Sum current part %e, other part %e, total length %e, "
                        "remainder %e\n",
                        sum_current_edge_part, sum_other_edge_part, current_edge_cost,
                        remainder);
                output_function(output_buffer);
            }

            if (current_cluster_index == other_cluster_index) {
                stats.num_merged_edge_events += 1;
        
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize, "Clusters already merged, ignoring edge\n");
                    output_function(output_buffer);
                }
                
                other_edge_part.deleted = true;
        
                continue;
            }

            if (remainder < eps * current_edge_cost || remainder == 0.0) {
                stats.total_num_merge_events += 1;

                phase1_result.push_back(next_edge_part_index / 2);
                
                other_edge_part.deleted = true;
                int new_cluster_index = clusters.size();
        
                clusters.push_back(Cluster(&pairing_heap_buffer));
                
                Cluster& new_cluster = clusters[new_cluster_index];
                Cluster& current_cluster_ref = clusters[current_cluster_index];
                Cluster& other_cluster_ref = clusters[other_cluster_index];

                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Merge %d and %d into %d\n", current_cluster_index,
                            other_cluster_index, new_cluster_index);
                    output_function(output_buffer);
                }

                new_cluster.prize_sum = current_cluster_ref.prize_sum + other_cluster_ref.prize_sum;
                new_cluster.subcluster_moat_sum = current_cluster_ref.subcluster_moat_sum + other_cluster_ref.subcluster_moat_sum;
                new_cluster.contains_root = current_cluster_ref.contains_root || other_cluster_ref.contains_root;
                new_cluster.active = !new_cluster.contains_root;
                new_cluster.merged_along = next_edge_part_index / 2;
                new_cluster.child_cluster_1 = current_cluster_index;
                new_cluster.child_cluster_2 = other_cluster_index;
                current_cluster_ref.active = false;
                current_cluster_ref.active_end_time = current_time + remainder;
                current_cluster_ref.merged_into = new_cluster_index;
                current_cluster_ref.moat = current_cluster_ref.active_end_time - current_cluster_ref.active_start_time;
                
                if (current_cluster_ref.moat < 0.0) current_cluster_ref.moat = 0.0;

                clusters_deactivation.delete_element(current_cluster_index);
        
                num_active_clusters -= 1;
        
                if (!current_cluster_ref.edge_parts.is_empty()) {
                    clusters_next_edge_event.delete_element(current_cluster_index);
                }

                if (other_cluster_ref.active) {
                    stats.num_active_active_merge_events += 1;
                    other_cluster_ref.active = false;
                    other_cluster_ref.active_end_time = current_time + remainder;
                    other_cluster_ref.moat = other_cluster_ref.active_end_time - other_cluster_ref.active_start_time;
                    if (other_cluster_ref.moat < 0.0) other_cluster_ref.moat = 0.0;

                    clusters_deactivation.delete_element(other_cluster_index);
          
                    if (!other_cluster_ref.edge_parts.is_empty()) {
                        clusters_next_edge_event.delete_element(other_cluster_index);
                    }

                    num_active_clusters -= 1;
                } else {
                    stats.num_active_inactive_merge_events += 1;
          
                    if (!other_cluster_ref.contains_root) {
                        double edge_event_update_time = current_time + remainder - other_cluster_ref.active_end_time;
            
                        other_cluster_ref.edge_parts.add_to_heap(edge_event_update_time);
                        inactive_merge_events.emplace_back();
            
                        InactiveMergeEvent& merge_event = inactive_merge_events.back();
                        merge_event.active_cluster_index = current_cluster_index;
                        merge_event.inactive_cluster_index = other_cluster_index;
                        const auto& edge = edges[edge_index];
                        int active_node_part = edge.first;
                        int inactive_node_part = edge.second;
                        
                        if (next_edge_part_index % 2 == 1) {
                            std::swap(active_node_part, inactive_node_part);
                        }
            
                        merge_event.active_cluster_node = active_node_part;
                        merge_event.inactive_cluster_node = inactive_node_part;

                        if (static_cast<size_t>(edge_index) < edge_info.size()) {
                            edge_info[edge_index].inactive_merge_event = inactive_merge_events.size() - 1;
                        }
                    }
                }
        
                other_cluster_ref.merged_into = new_cluster_index;
                new_cluster.edge_parts = PairingHeapType::meld(&(current_cluster_ref.edge_parts), &(other_cluster_ref.edge_parts));
                new_cluster.subcluster_moat_sum += current_cluster_ref.moat;
                new_cluster.subcluster_moat_sum += other_cluster_ref.moat;

                if (new_cluster.active) {
                    new_cluster.active_start_time = current_time + remainder;
                    double becoming_inactive_time = new_cluster.active_start_time + new_cluster.prize_sum - new_cluster.subcluster_moat_sum;
                    
                    clusters_deactivation.insert(becoming_inactive_time, new_cluster_index);
                    
                    if (!new_cluster.edge_parts.is_empty()) {
                        double tmp_val;
                        int tmp_index;
                        
                        new_cluster.edge_parts.get_min(&tmp_val, &tmp_index);
                        clusters_next_edge_event.insert(tmp_val, new_cluster_index);
                    }
          
                    num_active_clusters += 1;
                }

                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Merged %d and %d into %d. New active: %d\n", current_cluster_index,
                            other_cluster_index, new_cluster_index, new_cluster.active);
                    output_function(output_buffer);
                }

            } else if (other_cluster.active) {
                stats.total_num_edge_growth_events += 1;
                stats.num_active_active_edge_growth_events += 1;
                double next_event_time = current_time + remainder * 0.5;

                if (!current_cluster.edge_parts.is_empty()) {
                    clusters_next_edge_event.delete_element(current_cluster_index);
                }
        
                next_edge_part.heap_node = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);
                
                if (!current_cluster.edge_parts.is_empty()) {
                    double tmp_val;
                    int tmp_idx;
                    
                    current_cluster.edge_parts.get_min(&tmp_val, &tmp_idx);
                    clusters_next_edge_event.insert(tmp_val, current_cluster_index);
                }

                clusters_next_edge_event.delete_element(other_cluster_index);
                
                double original_from_value_aa = other_cluster.active_start_time + other_edge_part.next_event_val - other_finished_moat_sum;
        
                other_cluster.edge_parts.decrease_key(other_edge_part.heap_node, original_from_value_aa, next_event_time);
        
                if (!other_cluster.edge_parts.is_empty()) {
                    double tmp_val; 
                    int tmp_idx;
          
                    other_cluster.edge_parts.get_min(&tmp_val, &tmp_idx);
                    clusters_next_edge_event.insert(tmp_val, other_cluster_index);
                }
        
                next_edge_part.next_event_val = sum_current_edge_part + remainder * 0.5;
                other_edge_part.next_event_val = sum_other_edge_part + remainder * 0.5;

                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize, "Added new event at time %e\n", next_event_time);
                    output_function(output_buffer);
                }
            } else {
                stats.total_num_edge_growth_events += 1;
                stats.num_active_inactive_edge_growth_events += 1;
                double next_event_time = current_time + remainder;

                if (!current_cluster.edge_parts.is_empty()) {
                    clusters_next_edge_event.delete_element(current_cluster_index);
                }

                next_edge_part.heap_node = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);
                
                if (!current_cluster.edge_parts.is_empty()) {
                    double tmp_val; 
                    int tmp_idx;
                    
                    current_cluster.edge_parts.get_min(&tmp_val, &tmp_idx);
                    clusters_next_edge_event.insert(tmp_val, current_cluster_index);
                }

                double original_from_value_ai = other_cluster.active_end_time + other_edge_part.next_event_val - other_finished_moat_sum;
                
                other_cluster.edge_parts.decrease_key(other_edge_part.heap_node, original_from_value_ai, other_cluster.active_end_time);
                
                next_edge_part.next_event_val = current_edge_cost - other_finished_moat_sum;
                other_edge_part.next_event_val = other_finished_moat_sum;

                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Added new event at time %e and and event for inactive edge "
                            "part\n", next_event_time);
                    output_function(output_buffer);
                }
            }
        } else {
            stats.num_cluster_events += 1;
            current_time = next_cluster_time;
        
            if (static_cast<size_t>(next_cluster_index) >= clusters.size()) {
                remove_next_cluster_event();
                
                continue;
            }
      
            remove_next_cluster_event();

            Cluster& current_cluster = clusters[next_cluster_index];
            
            if (!current_cluster.active) continue;

            current_cluster.active = false;
            current_cluster.active_end_time = current_time;
            current_cluster.moat = current_cluster.active_end_time - current_cluster.active_start_time;
            
            if (current_cluster.moat < 0.0) current_cluster.moat = 0.0;

            if (!current_cluster.edge_parts.is_empty()) {
                clusters_next_edge_event.delete_element(next_cluster_index);
            }
      
            num_active_clusters -= 1;

            if (verbosity_level >= 2) {
                snprintf(output_buffer, kOutputBufferSize,
                        "Cluster deactivation: cluster %d at time %e (moat size %e)\n",
                        next_cluster_index, current_time, current_cluster.moat);
                output_function(output_buffer);
            }

        }
    }

    if (verbosity_level >= 1) {
        snprintf(output_buffer, kOutputBufferSize,
                "Finished GW clustering: final event time %e, number of edge events "
                "%lld\n", current_time, stats.total_num_edge_events);
        output_function(output_buffer);
    }

    node_good.assign(num_nodes, 0);

    if (root >= 0) {
        for (size_t ii = 0; ii < clusters.size(); ++ii) {
            if (clusters[ii].contains_root && clusters[ii].merged_into == -1) {
                mark_nodes_as_good(ii);
                
                break;
            }
        }
    } else {
        for (size_t ii = 0; ii < clusters.size(); ++ii) {
            if (clusters[ii].active) {
                mark_nodes_as_good(ii);
            }
        }
    }

    if (pruning == kNoPruning) {
        build_phase1_node_set(phase1_result, result_nodes);
    
        *result_edges = phase1_result;
    
        return true;
    }

    if (verbosity_level >= 2) {
        snprintf(output_buffer, kOutputBufferSize, "------------------------------------------\n");
        output_function(output_buffer);
        snprintf(output_buffer, kOutputBufferSize, "Starting pruning\n");
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

    if (pruning == kSimplePruning) {
        build_phase2_node_set(result_nodes);
    
        *result_edges = phase2_result;
    
        return true;
    }

    phase3_result.clear();
    phase3_result.reserve(phase2_result.size());
    
    for (auto& neighbors : phase3_neighbors) { neighbors.clear(); }
  
    if (phase3_neighbors.size() != num_nodes) phase3_neighbors.resize(num_nodes);
  
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

    if (pruning == kGWPruning) {
        if (verbosity_level >= 2) {
            snprintf(output_buffer, kOutputBufferSize, "Starting GW pruning, phase 2 result size: %zu\n", phase2_result.size());
        }

        for(auto& cluster : clusters) { cluster.necessary = false; }

        for (int ii = static_cast<int>(phase2_result.size()) - 1; ii >= 0; --ii) {
            int current_edge_index = phase2_result[ii];
            
            if(static_cast<size_t>(current_edge_index) >= edges.size() || static_cast<size_t>(current_edge_index) >= edge_info.size()) continue;

            const auto& edge = edges[current_edge_index];
            int uu = edge.first;
            int vv = edge.second;
            bool u_deleted = (static_cast<size_t>(uu) >= node_deleted.size() || node_deleted[uu]);
            bool v_deleted = (static_cast<size_t>(vv) >= node_deleted.size() || node_deleted[vv]);

            if (u_deleted && v_deleted) {
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "GW: Edge %d (%d, %d) endpoints deleted, skipping.\n", current_edge_index, uu, vv);
                    output_function(output_buffer);
                }
                
                continue;
            }

            int inactive_merge_idx = edge_info[current_edge_index].inactive_merge_event;
      
            if (inactive_merge_idx < 0) {
                mark_clusters_as_necessary(uu);
                mark_clusters_as_necessary(vv);
                phase3_result.push_back(current_edge_index);

                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "GW: Active-Active edge %d (%d, %d), keeping.\n", current_edge_index, uu, vv);
                    output_function(output_buffer);
                }
            } else {
                if (static_cast<size_t>(inactive_merge_idx) >= inactive_merge_events.size()) continue;

                const InactiveMergeEvent& current_merge_event = inactive_merge_events[inactive_merge_idx];
                int inactive_cluster_index = current_merge_event.inactive_cluster_index;
                bool is_necessary = (static_cast<size_t>(inactive_cluster_index) < clusters.size() && clusters[inactive_cluster_index].necessary);

                if (is_necessary) {
                    phase3_result.push_back(current_edge_index);
                    mark_clusters_as_necessary(current_merge_event.inactive_cluster_node);
                    mark_clusters_as_necessary(current_merge_event.active_cluster_node);

                    if (verbosity_level >= 2) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "GW: Inactive side cluster %d necessary for edge %d (%d, %d), keeping.\n",
                                inactive_cluster_index, current_edge_index, uu, vv);
                        output_function(output_buffer);
                    }
                } else {
                    mark_nodes_as_deleted(current_merge_event.inactive_cluster_node, current_merge_event.active_cluster_node);

                    if (verbosity_level >= 2) {
                        snprintf(output_buffer, kOutputBufferSize,
                                "GW: Inactive side cluster %d not necessary for edge %d (%d, %d), discarding inactive side.\n",
                                inactive_cluster_index, current_edge_index, uu, vv);
                        output_function(output_buffer);
                    }
                }
            }
        }

        build_phase3_node_set(result_nodes);
        
        *result_edges = phase3_result;
    
        return true;
    } else if (pruning == kStrongPruning) {
        if (verbosity_level >= 2) {
            snprintf(output_buffer, kOutputBufferSize,
                    "Starting Strong pruning, phase 2 result size: %zu\n", phase2_result.size());
        }

        if (final_component_label.size() != num_nodes) final_component_label.resize(num_nodes);
        
        std::fill(final_component_label.begin(), final_component_label.end(), -1);
        final_components.clear();
        
        root_component_index = -1;
        
        if (strong_pruning_parent.size() != num_nodes) strong_pruning_parent.resize(num_nodes);
    
        std::fill(strong_pruning_parent.begin(), strong_pruning_parent.end(), make_pair(-1,-1.0));
        
        if (strong_pruning_payoff.size() != num_nodes) strong_pruning_payoff.resize(num_nodes);
        
        std::fill(strong_pruning_payoff.begin(), strong_pruning_payoff.end(), -1.0);

        for (int edge_idx : phase2_result) {
            if (static_cast<size_t>(edge_idx) >= edges.size()) continue;
            
            const auto& edge = edges[edge_idx];
            int uu = edge.first;
            int vv = edge.second;
            
            if (static_cast<size_t>(uu) < num_nodes && final_component_label[uu] == -1) {
                final_components.emplace_back();
                label_final_component(uu, static_cast<int>(final_components.size() - 1));
            }
    
            if (static_cast<size_t>(vv) < num_nodes && final_component_label[vv] == -1) {
                final_components.emplace_back();
                label_final_component(vv, static_cast<int>(final_components.size() - 1));
            }
        }

        for (size_t ii = 0; ii < final_components.size(); ++ii) {
            if (final_components[ii].empty()) continue;

            if (verbosity_level >= 2) {
                snprintf(output_buffer, kOutputBufferSize, "Strong pruning on final component %zu (size %zu):\n",
                        ii, final_components[ii].size());
                output_function(output_buffer);
            }

            if (static_cast<int>(ii) == root_component_index) {
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize, "Component contains root, pruning starting at %d\n", root);
                    output_function(output_buffer);
                }
         
                if (static_cast<size_t>(root) < num_nodes) {
                    strong_pruning_from(root, true);
                }
            } else {
                int best_component_root = find_best_component_root(ii);
         
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Best start node for component %zu: %d, pruning from there\n",
                            ii, best_component_root);
                    output_function(output_buffer);
                }
         
                if (best_component_root != -1 && static_cast<size_t>(best_component_root) < num_nodes) {
                    strong_pruning_from(best_component_root, true);
                }
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

            if (u_deleted || v_deleted) {
                if (verbosity_level >= 2) {
                    snprintf(output_buffer, kOutputBufferSize,
                            "Strong: Not keeping edge %d (%d, %d) due to deleted endpoint(s).\n",
                            current_edge_index, uu, vv);
                    output_function(output_buffer);
                }
            } else {
                phase3_result.push_back(current_edge_index);
            }
        }

        build_phase3_node_set(result_nodes);
        
        *result_edges = phase3_result;
    
        return true;
    }

    snprintf(output_buffer, kOutputBufferSize, "Error: unknown pruning scheme.\n");
    output_function(output_buffer);
  
    return false;
}

void PCSTFast::label_final_component(int start_node_index,  int new_component_index) {
    cluster_queue.clear();
    
    if (static_cast<size_t>(start_node_index) >= final_component_label.size() || final_component_label[start_node_index] != -1) return;
  
    if (static_cast<size_t>(new_component_index) >= final_components.size()) return;

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
            for (const auto& neighbor_pair : phase3_neighbors[current_node_index]) {
                int next_node_index = neighbor_pair.first;
                
                if (static_cast<size_t>(next_node_index) < final_component_label.size() && final_component_label[next_node_index] == -1) {
                    final_component_label[next_node_index] = new_component_index;
                    
                    cluster_queue.push_back(next_node_index);
                    
                    if (static_cast<size_t>(new_component_index) < final_components.size()) {
                        final_components[new_component_index].push_back(next_node_index);
                    }

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
  
    if (static_cast<size_t>(start_node_index) >= prizes.size() || static_cast<size_t>(start_node_index) >= strong_pruning_parent.size()) return;

    strong_pruning_parent[start_node_index] = make_pair(-1, 0.0);
    
    stack.emplace_back(true, start_node_index);

    while (!stack.empty()) {
        bool begin = stack.back().first;
        int current_node_index = stack.back().second;
        
        stack.pop_back();

        if (static_cast<size_t>(current_node_index) >= prizes.size() ||
            static_cast<size_t>(current_node_index) >= strong_pruning_parent.size() ||
            static_cast<size_t>(current_node_index) >= strong_pruning_payoff.size()) continue;
        
        if (begin) {
            stack.emplace_back(false, current_node_index);
            
            if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                for (const auto& neighbor_pair : phase3_neighbors[current_node_index]) {
                    int next_node_index = neighbor_pair.first;
                    double next_cost = neighbor_pair.second;

                    if (next_node_index == strong_pruning_parent[current_node_index].first) {
                        continue;
                    }

                    if (static_cast<size_t>(next_node_index) < prizes.size() && 
                        static_cast<size_t>(next_node_index) < strong_pruning_parent.size()) {
                            strong_pruning_parent[next_node_index] = make_pair(current_node_index, next_cost);
                            
                            stack.emplace_back(true, next_node_index);
                    }
                }
            }
        } else {
            strong_pruning_payoff[current_node_index] = prizes[current_node_index];
            
            if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
                for (const auto& neighbor_pair : phase3_neighbors[current_node_index]) {
                    int next_node_index = neighbor_pair.first;
                    double next_cost = neighbor_pair.second;

                    if (next_node_index == strong_pruning_parent[current_node_index].first) {
                        continue;
                    }

                    if (static_cast<size_t>(next_node_index) < prizes.size() &&
                        static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                        static_cast<size_t>(next_node_index) < strong_pruning_payoff.size() &&
                        strong_pruning_parent[next_node_index].first == current_node_index) {
                            
                            double next_payoff = strong_pruning_payoff[next_node_index] - next_cost;

                            if (next_payoff <= 0.0) {
                                if (mark_as_deleted) {
                                    if (verbosity_level >= 2) {
                                        snprintf(output_buffer, kOutputBufferSize,
                                                "Strong: Subtree at %d non-positive (%e), pruning.\n",
                                                next_node_index, next_payoff);
                                        output_function(output_buffer);
                                    }
                                    
                                    mark_nodes_as_deleted(next_node_index, current_node_index);
                                }
                            } else {
                                strong_pruning_payoff[current_node_index] += next_payoff;
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

    const auto& component_nodes = final_components[component_index];
    int current_best_root = component_nodes[0];

    if (static_cast<size_t>(current_best_root) >= prizes.size()) return -1;

    strong_pruning_from(current_best_root, false);
  
    if (static_cast<size_t>(current_best_root) >= strong_pruning_payoff.size()) return -1;
  
    double current_best_value = strong_pruning_payoff[current_best_root];

    stack2.clear();
   
    if (static_cast<size_t>(current_best_root) < phase3_neighbors.size()) {
        for (const auto& neighbor_pair : phase3_neighbors[current_best_root]) {
            int next_node_index = neighbor_pair.first;
            
            if (static_cast<size_t>(next_node_index) < final_component_label.size() &&
                final_component_label[next_node_index] == component_index &&
                static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                strong_pruning_parent[next_node_index].first == current_best_root) {
                    stack2.emplace_back(next_node_index);
            }
        }
    }

    while (!stack2.empty()) {
        int current_node_index = stack2.back();
    
        stack2.pop_back();

        if (static_cast<size_t>(current_node_index) >= prizes.size() ||
            static_cast<size_t>(current_node_index) >= strong_pruning_parent.size() ||
            static_cast<size_t>(current_node_index) >= strong_pruning_payoff.size()) continue;

        int current_parent_index = strong_pruning_parent[current_node_index].first;
        double parent_edge_cost = strong_pruning_parent[current_node_index].second;

        if (current_parent_index == -1 ||
            static_cast<size_t>(current_parent_index) >= prizes.size() ||
            static_cast<size_t>(current_parent_index) >= strong_pruning_payoff.size()) continue;

        double current_payoff_parent = strong_pruning_payoff[current_parent_index];
        double current_payoff_self = strong_pruning_payoff[current_node_index];
        double current_node_contribution_to_parent = 0.0;
    
        if (current_payoff_self > parent_edge_cost) {
            current_node_contribution_to_parent = current_payoff_self - parent_edge_cost;
        }

        double parent_payoff_without_current_node = current_payoff_parent - current_node_contribution_to_parent;
        double parent_contribution_to_current_node = 0.0;

        if (parent_payoff_without_current_node > parent_edge_cost) {
            parent_contribution_to_current_node = parent_payoff_without_current_node - parent_edge_cost;
        }

        double current_node_total_payoff = current_payoff_self + parent_contribution_to_current_node;

        if (current_node_total_payoff > current_best_value) {
            current_best_root = current_node_index;
            current_best_value = current_node_total_payoff;
        }
        
        strong_pruning_payoff[current_node_index] = current_node_total_payoff;

        if (static_cast<size_t>(current_node_index) < phase3_neighbors.size()) {
            for (const auto& neighbor_pair : phase3_neighbors[current_node_index]) {
                int next_node_index = neighbor_pair.first;
                
                if (static_cast<size_t>(next_node_index) < prizes.size() &&
                    static_cast<size_t>(next_node_index) < strong_pruning_parent.size() &&
                    strong_pruning_parent[next_node_index].first == current_node_index) {
                        stack2.emplace_back(next_node_index);
                }
            }
        }
    }

    return current_best_root;
}

void PCSTFast::build_phase3_node_set(std::vector<int>* node_set) {
    node_set->clear();
    node_set->reserve(prizes.size());
    
    for (int ii = 0; ii < static_cast<int>(prizes.size()); ++ii) {
        if (static_cast<size_t>(ii) < node_good.size() && node_good[ii] &&
            static_cast<size_t>(ii) < node_deleted.size() && !node_deleted[ii]) {
                node_set->push_back(ii);
        }
    }
}

void PCSTFast::build_phase2_node_set(std::vector<int>* node_set) {
    node_set->clear();
    node_set->reserve(prizes.size());

    for (int ii = 0; ii < static_cast<int>(prizes.size()); ++ii) {
        if (static_cast<size_t>(ii) < node_good.size() && node_good[ii]) {
            node_set->push_back(ii);
        }
    }
}

void PCSTFast::build_phase1_node_set(const std::vector<int>& edge_set, std::vector<int>* node_set) {
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
    for (int ii = 0; ii < static_cast<int>(prizes.size()); ++ii) {
        if (static_cast<size_t>(ii) < node_good.size() && node_good[ii] &&
            static_cast<size_t>(ii) < build_phase1_included_nodes.size() && !build_phase1_included_nodes[ii]) {
                node_set->push_back(ii);
        }
    }
}

void PCSTFast::get_statistics(PCSTFast::Statistics* s) const {
    *s = stats;
}
