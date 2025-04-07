// pcst_fast_test.cc
#include "pcst_fast.h"
#include "test_helpers.h"
#include "logger.h" // Include the new logger header

#include <vector>
#include <iterator>
#include <utility>
#include <ranges>
#include <string>
#include <iostream>
#include <memory> // For shared_ptr if needed, or stack allocation

#include "gtest/gtest.h"

namespace ca = cluster_approx;
namespace ca_internal = cluster_approx::internal;

using test_helpers::WriteToStderr;
using test_helpers::CheckResult;

namespace {
    // Use the LogLevel enum for verbosity control in tests
    static constexpr ca_internal::LogLevel kTestLogLevel = ca_internal::LogLevel::ERROR; // Or WARNING, INFO, etc.
}

// LogSink implementation for tests, writing to stderr via test helper
void TestLogSink(ca_internal::LogLevel level, const std::string& message) {
    // Optional: Add level prefix for clarity in test output
    std::string prefix;
    switch (level) {
        case ca_internal::LogLevel::FATAL:   prefix = "[FATAL] ";   break;
        case ca_internal::LogLevel::ERROR:   prefix = "[ERROR] ";   break;
        case ca_internal::LogLevel::WARNING: prefix = "[WARNING] "; break;
        case ca_internal::LogLevel::INFO:    prefix = "[INFO] ";    break;
        case ca_internal::LogLevel::DEBUG:   prefix = "[DEBUG] ";   break;
        case ca_internal::LogLevel::TRACE:   prefix = "[TRACE] ";   break;
        case ca_internal::LogLevel::NONE:    return;
    }

    // Use the WriteToStderr helper if available and noexcept, otherwise fallback
    if constexpr (noexcept(WriteToStderr(nullptr))) {
        WriteToStderr((prefix + message).c_str());
    } else {
        std::cerr << prefix << message << std::flush;
    }
}


void RunAlgo(const std::vector<std::pair<int, int>>& edges,
             const std::vector<double>& prizes,
             const std::vector<double>& costs,
             ca::PCSTFast::IndexType root,
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const std::vector<int>& expected_node_result,
             const std::vector<int>& expected_edge_result) {
    std::vector<ca::PCSTFast::IndexType> node_result;
    std::vector<ca::PCSTFast::IndexType> edge_result;

    // Create the logger instance for this test run
    ca_internal::Logger test_logger(TestLogSink, kTestLogLevel);

    // Instantiate PCSTFast, passing the logger reference
    ca::PCSTFast algo(edges, prizes, costs, root, target_num_active_clusters,
                      pruning, test_logger); // Pass logger reference

    ASSERT_TRUE(algo.run(&node_result, &edge_result));

    // Sort results for comparison
    std::ranges::sort(node_result);
    std::ranges::sort(edge_result);

    std::vector<ca::PCSTFast::IndexType> sorted_expected_node_result(
        expected_node_result.begin(), expected_node_result.end());
    std::ranges::sort(sorted_expected_node_result);

    std::vector<ca::PCSTFast::IndexType> sorted_expected_edge_result(
        expected_edge_result.begin(), expected_edge_result.end());
    std::ranges::sort(sorted_expected_edge_result);

    // Check results
    CheckResult(sorted_expected_node_result, node_result);
    CheckResult(sorted_expected_edge_result, edge_result);
}

// Template overloads remain the same, just call the base RunAlgo
template <size_t N1, size_t N2, size_t N3, size_t N4>
void RunAlgo(const std::vector<std::pair<int, int>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             ca::PCSTFast::IndexType root,
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3],
             const int (&expected_edge_result)[N4]) {
    std::vector<double> prizes_vec(std::begin(prizes), std::end(prizes));
    std::vector<double> costs_vec(std::begin(costs), std::end(costs));
    std::vector<int> expected_node_result_vec(std::begin(expected_node_result), std::end(expected_node_result));
    std::vector<int> expected_edge_result_vec(std::begin(expected_edge_result), std::end(expected_edge_result));

    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_active_clusters, pruning,
            expected_node_result_vec, expected_edge_result_vec);
}

template <size_t N1, size_t N2, size_t N3>
void RunAlgo(const std::vector<std::pair<int, int>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             ca::PCSTFast::IndexType root,
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3]) {
    std::vector<double> prizes_vec(std::begin(prizes), std::end(prizes));
    std::vector<double> costs_vec(std::begin(costs), std::end(costs));
    std::vector<int> expected_node_result_vec(std::begin(expected_node_result), std::end(expected_node_result));
    std::vector<int> expected_edge_result_vec{}; // Empty edge result

    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_active_clusters, pruning,
            expected_node_result_vec, expected_edge_result_vec);
}

// --- All TEST cases remain exactly the same ---

TEST(PCSTFastTest, SimpleTestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kNoPruning;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {0, 1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, SimpleTestUnrootedNoPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kNoPruning;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, SimpleTestUnrootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, SimpleTestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {0, 5, 6};
    const double costs[] = {3, 4};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kStrongPruning;
    const int node_result[] = {1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple2TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}};
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kNoPruning;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple2TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}};
    const double prizes[] = {10, 0, 1, 10};
    const double costs[] = {10, 4, 3};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple3TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}};
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kNoPruning;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple3TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}};
    const double prizes[] = {10, 10, 1, 10};
    const double costs[] = {10, 6, 5};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0, 1, 2, 3};
    const int edge_result[] = {0, 1, 2};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple4TestRootedNoPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kNoPruning;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple4TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    ca::PCSTFast::IndexType root = 0;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple4TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 2;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0, 1, 2};
    const int edge_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple4bTestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {10, 3, 3};
    const double costs[] = {100, 2};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple5TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}};
    const double prizes[] = {10, 0, 6, 6};
    const double costs[] = {100, 2, 5};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 2;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {0, 2, 3};
    const int edge_result[] = {2};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Medium1TestRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}, {0, 9}, {0, 2},
                                    {0, 3}, {0, 5}, {1, 9}, {1, 3}, {1, 5},
                                    {1, 7}, {2, 8}, {2, 3}, {3, 4}, {3, 5},
                                    {3, 6}, {3, 7}, {3, 8}, {3, 9}, {4, 5},
                                    {4, 6}, {4, 7}, {5, 8}, {6, 8}};
    const double prizes[] = {0.032052554364677466,
                            0.32473378289799926,
                            0.069699345546302638,
                            0,
                            0.74867253235151754,
                            0.19804330340026255,
                            0.85430521133171622,
                            0.83819939651391351,
                            0.71744625276884877,
                            0.016798567754083948};
    const double costs[] = {0.8, 0.8, 0.8800000000000001, 0.8, 0.8,
                            0.8800000000000001, 0.8, 0.8, 0.8800000000000001, 0.8,
                            0.8, 0.8, 0.8800000000000001, 0.8800000000000001, 0.8800000000000001,
                            0.8800000000000001, 0.8800000000000001, 0.8800000000000001, 0.8800000000000001, 0.8,
                            0.8, 0.8, 0.8, 0.8};
    ca::PCSTFast::IndexType root = 3;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {3, 4, 6, 7, 8};
    const int edge_result[] = {16, 20, 21, 23};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Medium1TestRootedStrongPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}, {0, 9}, {0, 2},
                                    {0, 3}, {0, 5}, {1, 9}, {1, 3}, {1, 5},
                                    {1, 7}, {2, 8}, {2, 3}, {3, 4}, {3, 5},
                                    {3, 6}, {3, 7}, {3, 8}, {3, 9}, {4, 5},
                                    {4, 6}, {4, 7}, {5, 8}, {6, 8}};
    const double prizes[] = {0.032052554364677466,
                            0.32473378289799926,
                            0.069699345546302638,
                            0,
                            0.74867253235151754,
                            0.19804330340026255,
                            0.85430521133171622,
                            0.83819939651391351,
                            0.71744625276884877,
                            0.016798567754083948};
    const double costs[] = {0.8, 0.8, 0.8800000000000001, 0.8, 0.8,
                            0.8800000000000001, 0.8, 0.8, 0.8800000000000001, 0.8,
                            0.8, 0.8, 0.8800000000000001, 0.8800000000000001, 0.8800000000000001,
                            0.8800000000000001, 0.8800000000000001, 0.8800000000000001, 0.8800000000000001, 0.8,
                            0.8, 0.8, 0.8, 0.8};
    ca::PCSTFast::IndexType root = 3;
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kStrongPruning;
    const int node_result[] = {3};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple6TestUnRootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const double prizes[] = {100.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 100.0};
    const double costs[] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    int root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {7};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple6TestUnRootedCFCPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const double prizes[] = {100.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 100.0};
    const double costs[] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    int root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kConnectFinalComponents;
    const int node_result[] = {0, 1, 2, 3, 4, 5, 6, 7};
    const int edge_result[] = {0, 1, 2, 3, 4, 5, 6};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple7TestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {2, 3}, {3, 4}};
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kStrongPruning;
    const int node_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple7TestUnrootedGWPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {2, 3}, {3, 4}};
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    int root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;
    const int node_result[] = {1};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
}

TEST(PCSTFastTest, Simple7TestUnrootedCFCPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {2, 3}, {3, 4}};
    const double prizes[] = {0, 2.2, 0, 0, 2.1};
    const double costs[] = {1, 1, 1, 1};
    int root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kConnectFinalComponents;
    const int node_result[] = {0, 1, 2, 3, 4};
    const int edge_result[] = {0, 1, 2, 3};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}

TEST(PCSTFastTest, Simple8TestUnrootedStrongPruning) {
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}};
    const double prizes[] = {2, 2, 2};
    const double costs[] = {0, 5};
    ca::PCSTFast::IndexType root = ca::PCSTFast::kNoRoot;
    int target_num_active_clusters = 1;
    auto pruning = ca::PCSTFast::PruningMethod::kStrongPruning;
    const int node_result[] = {0, 1};
    const int edge_result[] = {0};

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result, edge_result);
}
