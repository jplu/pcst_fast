#pragma once

#include <string_view>
#include <utility>
#include <vector>
#include <cstdint>

#include "pairing_heap.h"
#include "priority_queue.h"

namespace cluster_approx {

    class PCSTFast {
    public:
        enum class PruningMethod {
            kNoPruning = 0,
            kSimplePruning,
            kGWPruning,
            kStrongPruning,
            kUnknownPruning,
        };

        struct Statistics {
            long long total_num_edge_events = 0;
            long long num_deleted_edge_events = 0;
            long long num_merged_edge_events = 0;
            long long total_num_merge_events = 0;
            long long num_active_active_merge_events = 0;
            long long num_active_inactive_merge_events = 0;
            long long total_num_edge_growth_events = 0;
            long long num_active_active_edge_growth_events = 0;
            long long num_active_inactive_edge_growth_events = 0;
            long long num_cluster_events = 0;

            Statistics() = default;
        };

        static constexpr int kNoRoot = -1;

        [[nodiscard]] static PruningMethod parse_pruning_method(std::string_view input);

        PCSTFast(const std::vector<std::pair<int, int> >& edges_,
                 const std::vector<double>& prizes_,
                 const std::vector<double>& costs_,
                 int root_,
                 int target_num_active_clusters_,
                 PruningMethod pruning_,
                 int verbosity_level_,
                 void (*output_function_)(const char*));

        ~PCSTFast();

        [[nodiscard]] bool run(std::vector<int>* result_nodes, std::vector<int>* result_edges);

        void get_statistics(Statistics* s) const;

    private:
        using PairingHeapType = PairingHeap<double, int>;
        using PriorityQueueType = PriorityQueue<double, int>;

        struct EdgeInfo {
            int inactive_merge_event = -1;
        };

        struct EdgePart {
            double next_event_val = 0.0;
            bool deleted = false;
            PairingHeapType::ItemHandle heap_node = nullptr;
        };

        struct InactiveMergeEvent {
            int active_cluster_index;
            int inactive_cluster_index;
            int active_cluster_node;
            int inactive_cluster_node;
        };

        struct Cluster {
            PairingHeapType edge_parts;
            bool active = false;
            double active_start_time = 0.0;
            double active_end_time = -1.0;
            int merged_into = -1;
            double prize_sum = 0.0;
            double subcluster_moat_sum = 0.0;
            double moat = 0.0;
            bool contains_root = false;
            int skip_up = -1;
            double skip_up_sum = 0.0;
            int merged_along = -1;
            int child_cluster_1 = -1;
            int child_cluster_2 = -1;
            bool necessary = false;

            explicit Cluster(std::vector<PairingHeapType::ItemHandle>* heap_buffer)
                : edge_parts(*heap_buffer)
                {}
        };

        const std::vector<std::pair<int, int> >& edges;
        const std::vector<double>& prizes;
        const std::vector<double>& costs;
        const int root;
        const int target_num_active_clusters;
        const PruningMethod pruning;
        const int verbosity_level;
        void (*output_function)(const char*);
        Statistics stats;
        std::vector<PairingHeapType::ItemHandle> pairing_heap_buffer;
        std::vector<EdgePart> edge_parts;
        std::vector<EdgeInfo> edge_info;
        std::vector<Cluster> clusters;
        std::vector<InactiveMergeEvent> inactive_merge_events;
        PriorityQueueType clusters_deactivation;
        PriorityQueueType clusters_next_edge_event;
        double current_time = 0.0;
        double eps = 1e-9;
        std::vector<uint8_t> node_good;
        std::vector<uint8_t> node_deleted;
        std::vector<int> phase1_result;
        std::vector<int> phase2_result;
        std::vector<int> phase3_result;
        std::vector<std::pair<int, double>> path_compression_visited;
        std::vector<int> cluster_queue;
        std::vector<std::vector<std::pair<int, double>>> phase3_neighbors;
        std::vector<uint8_t> build_phase1_included_nodes;
        std::vector<int> final_component_label;
        std::vector<std::vector<int>> final_components;
        int root_component_index = -1;
        std::vector<std::pair<int, double>> strong_pruning_parent;
        std::vector<double> strong_pruning_payoff;
        std::vector<std::pair<bool, int>> stack;
        std::vector<int> stack2;
        static constexpr int kOutputBufferSize = 10000;
        char output_buffer[kOutputBufferSize];

        void get_next_edge_event(double* next_time, int* next_cluster_index, int* next_edge_part_index);
        void remove_next_edge_event(int next_cluster_index);
        void get_next_cluster_event(double* next_time, int* next_cluster_index);
        void remove_next_cluster_event();
        void get_sum_on_edge_part(int edge_part_index, double* total_sum, double* finished_moat_sum, int* current_cluster_index);
        void mark_nodes_as_good(int start_cluster_index);
        void mark_clusters_as_necessary(int start_cluster_index);
        void mark_nodes_as_deleted(int start_node_index, int parent_node_index);
        void label_final_component(int start_node_index, int new_component_index);
        void strong_pruning_from(int start_node_index, bool mark_as_deleted);
        int find_best_component_root(int component_index);
        void build_phase1_node_set(const std::vector<int>& edge_set, std::vector<int>* node_set);
        void build_phase3_node_set(std::vector<int>* node_set);
        void build_phase2_node_set(std::vector<int>* node_set);
        [[nodiscard]] inline int get_other_edge_part_index(int edge_part_index) const {
            return edge_part_index ^ 1;
        }
    };
}
