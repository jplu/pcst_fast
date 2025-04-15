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
    std::vector<EdgeId> expected_edges_vec;
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
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kNone;

    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {0, 1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kNone;

    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, SimpleTestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kStrong;

    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple2TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kNone;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple2TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;

    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning, node_result);
}

TEST(EndToEndTest, Simple3TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kNone;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple3TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple4TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kNone;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple4TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = 0;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;

    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning, node_result);
}

TEST(EndToEndTest, Simple4TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = kInvalidNodeId;
    int target_num_clusters = 2;
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple4bTestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning, node_result);
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
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {0, 2, 3};
    const int edge_result[] = {2};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Medium1TestRootedGWPruning) {

    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {0, 9}, {0, 2}, {0, 3}, {0, 5}, {1, 9},
        {1, 3}, {1, 5}, {1, 7}, {2, 8}, {2, 3}, {3, 4}, {3, 5}, {3, 6},
        {3, 7}, {3, 8}, {3, 9}, {4, 5}, {4, 6}, {4, 7}, {5, 8}, {6, 8}
    };

    std::vector<double> prizes = {0.032052554364677466, 0.32473378289799926,
                                  0.069699345546302638, 0,
                                  0.74867253235151754, 0.19804330340026255,
                                  0.85430521133171622, 0.83819939651391351,
                                  0.71744625276884877, 0.016798567754083948
                                 };
    std::vector<double> costs = {
        0.8, 0.8, 0.88, 0.8, 0.8, 0.88, 0.8, 0.8, 0.88, 0.8, 0.8, 0.8,
        0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.8, 0.8, 0.8, 0.8, 0.8
    };
    int root = 3;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;

    const std::vector<int> node_result = {3, 4, 6, 7, 8};
    const std::vector<int> edge_result = {17, 20, 21, 23};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple6TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 4});
    edges.push_back({4, 5});
    edges.push_back({5, 6});
    edges.push_back({6, 7});
    const double prizes[] = {100.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 100.0};
    const double costs[] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;

    const int node_result[] = {0, 1, 2, 3, 4, 5, 6, 7};
    const int edge_result[] = {0, 1, 2, 3, 4, 5, 6};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple7TestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 4});
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kStrong;

    const int node_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning, node_result);
}

TEST(EndToEndTest, Simple7TestUnrootedGWPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 4});
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kGW;
    const int node_result[] = {0, 1, 2, 3, 4};
    const int edge_result[] = {0, 1, 2, 3};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}

TEST(EndToEndTest, Simple8TestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    const double prizes[] = {2, 2, 2};
    const double costs[] = {0, 5};
    int root = kInvalidNodeId;
    int target_num_clusters = 1;
    PruningMethod pruning = PruningMethod::kStrong;
    const int node_result[] = {0, 1};
    const int edge_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_clusters, pruning,
            node_result, edge_result);
}