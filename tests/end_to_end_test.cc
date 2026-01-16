#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "pcst_fast/pruning/no_pruner.h"
#include "pcst_fast/pruning/simple_pruner.h"
#include "pcst_fast/pruning/gw_pruner.h"
#include "pcst_fast/pruning/strong_pruner.h"

#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <span>
#include <memory>
#include <algorithm>
#include <iostream>

#include "gtest/gtest.h"
#include "test_helpers.h"

using namespace cluster_approx;

const int kVerbosityLevel = 0;

LogLevel get_test_log_level() {
    if (kVerbosityLevel <= 0) return LogLevel::FATAL;
    if (kVerbosityLevel == 1) return LogLevel::ERROR;
    if (kVerbosityLevel == 2) return LogLevel::WARNING;
    if (kVerbosityLevel == 3) return LogLevel::INFO;
    if (kVerbosityLevel == 4) return LogLevel::DEBUG;
    return LogLevel::TRACE;
}

// Main implementation taking vectors
void RunAlgo(const std::vector<std::pair<NodeId, NodeId>>& edges_vec,
             const std::vector<double>& prizes_vec,
             const std::vector<double>& costs_vec,
             NodeId root,
             int target_num_clusters_input,
             PruningMethod pruning_method,
             const std::vector<NodeId>& expected_node_result,
             const std::vector<EdgeId>& expected_edge_result) {

    StderrLogger logger(get_test_log_level());

    auto edges_span = std::span<const std::pair<NodeId, NodeId>>(edges_vec);
    auto prizes_span = std::span<const double>(prizes_vec);
    auto costs_span = std::span<const double>(costs_vec);

    GraphData graph {
        .edges = edges_span,
        .prizes = prizes_span,
        .costs = costs_span,
        .root = root
    };

    int internal_target_clusters = target_num_clusters_input;
    if (root != kInvalidNodeId) {
        ASSERT_EQ(target_num_clusters_input, 1) << "Target clusters must be 1 for rooted problems.";
        internal_target_clusters = 0;
    } else {
        ASSERT_GE(target_num_clusters_input, 1) << "Target clusters must be >= 1 for unrooted problems.";
    }

    std::unique_ptr<PCSTCoreAlgorithm> core_algo;
    CoreAlgorithmResult core_result;
    try {
        core_algo = std::make_unique<PCSTCoreAlgorithm>(graph, internal_target_clusters, &logger);
        core_result = core_algo->run();
    } catch (const std::exception& e) {
        FAIL() << "Core algorithm execution failed: " << e.what();
        return;
    }

    std::unique_ptr<IPruner> pruner;
    switch (pruning_method) {
    case PruningMethod::kNone:
        pruner = std::make_unique<pruning::NoPruner>();
        break;
    case PruningMethod::kSimple:
        pruner = std::make_unique<pruning::SimplePruner>();
        break;
    case PruningMethod::kGW:
        pruner = std::make_unique<pruning::GWPruner>();
        break;
    case PruningMethod::kStrong:
        pruner = std::make_unique<pruning::StrongPruner>();
        break;
    default:
        FAIL() << "Invalid pruning method enum value in test.";
        return;
    }

    PruningInput pruning_input {
        .graph = graph,
        .core_result = core_result,
        .logger = &logger
    };

    PruningResult final_result;
    try {
        final_result = pruner->prune(pruning_input);
    } catch (const std::exception& e) {
        FAIL() << "Pruner execution failed: " << e.what();
        return;
    }

    CheckResult(expected_node_result, final_result.nodes);
    CheckResult(expected_edge_result, final_result.edges);
}

// Helper template for C-style arrays
template <size_t N1, size_t N2, size_t N3, size_t N4>
void RunAlgo(const std::vector<std::pair<NodeId, NodeId>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             NodeId root,
             int target_num_clusters,
             PruningMethod pruning,
             const int (&expected_node_result)[N3],
             const int (&expected_edge_result)[N4]) {
    std::vector<double> prizes_vec(begin(prizes), end(prizes));
    std::vector<double> costs_vec(begin(costs), end(costs));
    std::vector<NodeId> expected_nodes_vec(begin(expected_node_result), end(expected_node_result));
    std::vector<EdgeId> expected_edges_vec(begin(expected_edge_result), end(expected_edge_result));
    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_clusters, pruning,
            expected_nodes_vec, expected_edges_vec);
}

// Helper template for C-style arrays (Node result only)
template <size_t N1, size_t N2, size_t N3>
void RunAlgo(const std::vector<std::pair<NodeId, NodeId>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             NodeId root,
             int target_num_clusters,
             PruningMethod pruning,
             const int (&expected_node_result)[N3]) {
    std::vector<double> prizes_vec(begin(prizes), end(prizes));
    std::vector<double> costs_vec(begin(costs), end(costs));
    std::vector<NodeId> expected_nodes_vec(begin(expected_node_result), end(expected_node_result));
    std::vector<EdgeId> expected_edges_vec; // Empty edges
    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_clusters, pruning,
            expected_nodes_vec, expected_edges_vec);
}

TEST(EndToEndTest, SimpleTestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = 0;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {0, 1};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kNone, node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kNone, node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kStrong, node_result, edge_result);
}

TEST(EndToEndTest, Simple2TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    int root = 0;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {1, 2}; // 0-1 is not selected because 1-2-3 component dies before reaching 0

    // IMPORTANT: edge_result here is checked against {1, 2}. 
    // Manual analysis suggests 1-2-3 dies. If it dies, how does NoPruner keep nodes 0..3?
    // Because NoPruner returns `phase1_edges`. 
    // If {1,2,3} dies, phase1 edges only contains {1,2} and {2,3}. 
    // Node 0 is root (good). Nodes 1,2,3 are involved in edges. 
    // So nodes {0,1,2,3} is correct. Edge {0} is NOT in phase1.
    
    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kNone, node_result, edge_result);
}

TEST(EndToEndTest, Simple2TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    int root = 0;
    
    // As analyzed, the component {1,2,3} dies due to costs > prizes.
    // It never merges with Root {0}.
    // Result: Root {0} is only valid cluster.
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result);
}

TEST(EndToEndTest, Simple3TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    int root = 0;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kNone, node_result, edge_result);
}

TEST(EndToEndTest, Simple3TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    int root = 0;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, Simple4TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = 0;
    // Edge (0,1) too expensive (100). (1,2) merges (cost 2, prize 6).
    // Result: 0 (root), 1, 2. Edge (1,2).
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kNone, node_result, edge_result);
}

TEST(EndToEndTest, Simple4TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = 0;
    // (1,2) merged but never connects to 0. 
    // GW pruning with root 0 -> 1-2 branch is disconnected from root, so pruned?
    // If not connected to root, GW pruner discards?
    // Yes, node_good[0]=true, others false.
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result);
}

TEST(EndToEndTest, Simple4TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = kInvalidNodeId;
    int target_num_clusters = 2; // Keep 2 clusters: {0} and {1,2}
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, Simple4bTestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    std::vector<double> prizes = {10, 3, 3};
    std::vector<double> costs = {100, 2};
    int root = kInvalidNodeId;
    
    // Same as above but target=1.
    // Cost 100 is too high for prizes 10 and 6. They will die before merging.
    // Core algo stops when everything inactive.
    // Largest cluster {0} (prize 10) vs {1,2} (prize 6).
    // Actually GW pruning behavior on disconnected forest:
    // Core marks all surviving Active clusters as good?
    // If they died, they are not "active".
    // "Unrooted case: Marking nodes from {} remaining active clusters."
    // If they all died, nothing marked good?
    // Wait, core algo marks nodes from `final_root_cluster` OR active clusters.
    // If {0} died and {1,2} died, nothing is active. Nothing marked good.
    // BUT PCST logic usually marks the largest component or relies on target clusters logic.
    // Here, expectation is {0}. Probably because 0 is the "first" or largest dead component?
    // Actually, implementation says: if unrooted, mark from active. If none active, `node_good` is all false.
    // Test expectation {0} implies 0 was somehow selected.
    // Re-reading `pcst_core_algorithm.cc`: "Marking nodes from {} remaining active clusters."
    // If target=1, but all die, num_active=0. 
    // Exception: If unrooted and everything dies, implementation might be returning empty.
    // Let's assume the provided test expectation is ground truth for desired behavior: {0}.
    // This implies {0} survived or was picked?
    // With Prize 10, Cost 100. 0 dies at t=10.
    // With Prize 6, Cost 100. {1,2} dies at t=3.
    // So 0 lives longer. Maybe that's why?
    // We stick to the expectation.
    
    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, 
            std::vector<NodeId>{0}, std::vector<EdgeId>{});
}

TEST(EndToEndTest, Simple5TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 0, 6, 6};
    const double costs[] = {100, 2, 5};
    int root = kInvalidNodeId;
    int target_num_clusters = 2;
    // {0} independent. {1,2,3}?
    // (2,3) cost 5, prizes 6,6. Merge.
    // (1,2) cost 2. Merge.
    // {1,2,3} sum prize 12.
    // 0 prize 10.
    // Connection cost 100.
    // Target 2 -> {0} and {1,2,3}.
    // GW pruning on {1,2,3}:
    // 1 has prize 0. Leaf 1 connected to 2 (cost 2).
    // Subtree 1 cost 2 > prize 0. Should be pruned?
    // Result {0, 2, 3}. Edge {2}. (2,3 is index 2).
    const int node_result[] = {0, 2, 3};
    const int edge_result[] = {2};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, PruningMethod::kGW,
            node_result, edge_result);
}

TEST(EndToEndTest, Medium1TestRootedGWPruning) {
    // Larger random-ish graph test
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {0, 9}, {0, 2}, {0, 3}, {0, 5}, {1, 9},
        {1, 3}, {1, 5}, {1, 7}, {2, 8}, {2, 3}, {3, 4}, {3, 5}, {3, 6},
        {3, 7}, {3, 8}, {3, 9}, {4, 5}, {4, 6}, {4, 7}, {5, 8}, {6, 8}
    };
    std::vector<double> prizes = {0.032052554364677466, 0.32473378289799926,
                                  0.069699345546302638, 0,
                                  0.74867253235151754, 0.19804330340026255,
                                  0.85430521133171622, 0.83819939651391351,
                                  0.71744625276884877, 0.016798567754083948};
    std::vector<double> costs = {
        0.8, 0.8, 0.88, 0.8, 0.8, 0.88, 0.8, 0.8, 0.88, 0.8, 0.8, 0.8,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.8, 0.8, 0.8, 0.8, 0.8
    };
    int root = 3;
    const std::vector<int> node_result = {3, 4, 6, 7, 8};
    const std::vector<int> edge_result = {17, 20, 21, 23};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, Simple6TestUnRootedGWPruning) {
    // Line graph 0..7. 0 and 7 have 100 prize. Middle 3 has 1. Others 0.
    // Costs 0.9.
    // 0 and 7 will essentially eat everything towards the center.
    // Merges will happen until 1 big component.
    const std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}
    };
    const double prizes[] = {100.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 100.0};
    const double costs[] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    int root = kInvalidNodeId;
    const int node_result[] = {0, 1, 2, 3, 4, 5, 6, 7};
    const int edge_result[] = {0, 1, 2, 3, 4, 5, 6};

    RunAlgo(edges, prizes, costs, root, 1, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, Simple7TestUnrootedStrongPruning) {
    // Star-like structure or branching.
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 4});
    std::vector<double> prizes = {0, 2.2, 0, 0, 2.1};
    std::vector<double> costs = {1, 1, 1, 1};
    
    // 0 is hub (prize 0).
    // 1 is leaf (2.2). Cost 1. Net 1.2.
    // 2-3-4 branch. 2(0), 3(0), 4(2.1). Path cost 1+1+1=3.
    // 4 prize 2.1 < Cost 3.
    // Branch 2-3-4 should be pruned by Strong Pruner.
    // Result: {0, 1}. But Strong Pruner marks nodes deleted. 
    // If 0 is kept, 1 is kept.
    // Expected result: {1} ?
    // If {0} has prize 0 and degree 1 (after pruning 2), does it survive?
    // Strong pruner maximizes profit.
    // Component {0,1}: Prize 2.2, Cost 1. Profit 1.2.
    // Component {1}: Prize 2.2. Cost 0. Profit 2.2.
    // Strong pruner finds best root.
    // If root is 1, neighbor 0 (cost 1). Subtree at 0 has prize 0. Net -1. Prune 0.
    // Result {1}.
    
    // Use vector implementation to pass empty edge list safely
    RunAlgo(edges, prizes, costs, kInvalidNodeId, 1, PruningMethod::kStrong, 
            std::vector<NodeId>{1}, std::vector<EdgeId>{});
}

TEST(EndToEndTest, Simple7TestUnrootedGWPruning) {
    // Same graph. GW pruning is less aggressive.
    // 1 merges with 0 (cost 1 < prize 2.2).
    // 4 merges with 3 (cost 1 < prize 2.1).
    // {3,4} merges with 2 (cost 1 < prize 2.1).
    // {2,3,4} merges with {0,1} (cost 1). Prize 2.2 vs 2.1.
    // Total prize 4.3. Total cost 4.
    // Everything merges.
    // GW pruning usually keeps the structure formed.
    const std::vector<std::pair<int, int>> edges = {
        {0, 1}, {0, 2}, {2, 3}, {3, 4}
    };
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    const int node_result[] = {0, 1, 2, 3, 4};
    const int edge_result[] = {0, 1, 2, 3};

    RunAlgo(edges, prizes, costs, kInvalidNodeId, 1, PruningMethod::kGW, node_result, edge_result);
}

TEST(EndToEndTest, Simple8TestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {2, 2, 2};
    const double costs[] = {0, 5}; // Edge 0 (0-1) cost 0. Edge 1 (1-2) cost 5.
    // {0,1} merge immediately. Prize 4.
    // {0,1} vs 2. Edge cost 5. Prize 4 vs 2.
    // Merged component {0,1,2} prize 6, cost 5. Net 1.
    // Strong Pruner optimization:
    // Best root?
    // If root 1:
    //   Branch 0 (cost 0, prize 2) -> Net 2. Keep.
    //   Branch 2 (cost 5, prize 2) -> Net -3. Prune.
    // Result {0,1}.
    const int node_result[] = {0, 1};
    const int edge_result[] = {0};

    RunAlgo(edges, prizes, costs, kInvalidNodeId, 1, PruningMethod::kStrong, node_result, edge_result);
}