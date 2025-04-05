#include "pcst_fast.h"
#include "test_helpers.h" // Assumes WriteToStderr and CheckResult are defined here

#include <algorithm>
#include <vector>
#include <iterator> // For std::begin, std::end
#include <utility>  // For std::pair
#include <ranges>   // For std::ranges::sort
#include <string>   // For std::string in logger lambda
#include <functional> // For std::function

#include "gtest/gtest.h"

namespace ca = cluster_approx;

// Assume test_helpers::WriteToStderr is defined as:
// namespace test_helpers {
//     inline void WriteToStderr(const char* s) { /* ... fprintf ... */ }
// }
using test_helpers::WriteToStderr;
// Assume test_helpers::CheckResult is defined as:
// namespace test_helpers {
//     template <typename T>
//     void CheckResult(const std::vector<T>& expected, const std::vector<T>& actual) { /* ... gtest checks ... */ }
// }
using test_helpers::CheckResult;

namespace {
    // Set verbosity for tests (0 = silent unless error, higher for more logs)
    // Keep low for standard runs to avoid clutter. Set higher for debugging.
    static constexpr int kVerbosityLevel = 0;
} // namespace

// Main helper function to run the algorithm and check results
void RunAlgo(const std::vector<std::pair<int, int>>& edges, // Input edges use int, compatible with IndexType=int
             const std::vector<double>& prizes,           // Input prizes use double, compatible with ValueType=double
             const std::vector<double>& costs,            // Input costs use double, compatible with ValueType=double
             ca::PCSTFast::IndexType root,                 // Use IndexType for root
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const std::vector<int>& expected_node_result, // Expected results still use int, will be cast/copied
             const std::vector<int>& expected_edge_result) {

    // ***** CHANGE: Declare result vectors using IndexType *****
    std::vector<ca::PCSTFast::IndexType> node_result;
    std::vector<ca::PCSTFast::IndexType> edge_result;

    // ***** CHANGE: Adapt logger signature *****
    // Create a lambda that captures the old WriteToStderr and matches the new Logger type.
    // This lambda ignores the level parameter, mimicking the old behavior.
    auto logger_lambda = [](int /*level*/, const std::string& message) {
        // Check if the old WriteToStderr function pointer is valid before calling
        if constexpr (noexcept(WriteToStderr(nullptr))) { // Basic check if callable
             WriteToStderr(message.c_str());
        } else {
             // Handle case where WriteToStderr might not be available or suitable
             // For tests, maybe just print to std::cerr directly?
             std::cerr << message;
             std::cerr.flush();
        }
    };
    // Or, simply use the default logger provided in pcst_fast.cc if suitable:
    // ca::PCSTFast::Logger test_logger = DefaultStderrLogger;
    ca::PCSTFast::Logger test_logger = logger_lambda; // Use the lambda wrapper


    // Instantiate the refactored PCSTFast class, passing the compatible logger
    // Input vectors (edges, prizes, costs) are passed by const&, compatible types.
    ca::PCSTFast algo(edges, prizes, costs, root, target_num_active_clusters,
                      pruning, kVerbosityLevel, test_logger); // Pass the logger

    // Run the algorithm - pointers to IndexType vectors are passed
    ASSERT_TRUE(algo.run(&node_result, &edge_result));

    // Sort actual results for comparison
    std::ranges::sort(node_result);
    std::ranges::sort(edge_result);

    // Create copies and sort expected results, ensuring type consistency
    // Copy from std::vector<int> to std::vector<IndexType>
    std::vector<ca::PCSTFast::IndexType> sorted_expected_node_result(
        expected_node_result.begin(), expected_node_result.end());
    std::ranges::sort(sorted_expected_node_result);

    std::vector<ca::PCSTFast::IndexType> sorted_expected_edge_result(
        expected_edge_result.begin(), expected_edge_result.end());
    std::ranges::sort(sorted_expected_edge_result);

    // Check results using the provided helper (assumes CheckResult works with IndexType)
    CheckResult(sorted_expected_node_result, node_result);
    CheckResult(sorted_expected_edge_result, edge_result);
}

// Template helper to allow using C-style arrays for input/expected results
// No changes needed here as it converts to std::vector<int> which are then handled by the main RunAlgo
template <size_t N1, size_t N2, size_t N3, size_t N4>
void RunAlgo(const std::vector<std::pair<int, int>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             ca::PCSTFast::IndexType root, // Use IndexType
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3],
             const int (&expected_edge_result)[N4]) {
    // Convert C-arrays to std::vectors
    std::vector<double> prizes_vec(std::begin(prizes), std::end(prizes));
    std::vector<double> costs_vec(std::begin(costs), std::end(costs));
    std::vector<int> expected_node_result_vec(std::begin(expected_node_result), std::end(expected_node_result));
    std::vector<int> expected_edge_result_vec(std::begin(expected_edge_result), std::end(expected_edge_result));

    // Call the main helper function
    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_active_clusters, pruning,
            expected_node_result_vec, expected_edge_result_vec);
}

// Template helper overload for cases where no edges are expected in the result
// No changes needed here
template <size_t N1, size_t N2, size_t N3>
void RunAlgo(const std::vector<std::pair<int, int>>& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             ca::PCSTFast::IndexType root, // Use IndexType
             int target_num_active_clusters,
             ca::PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3]) {
    // Convert C-arrays to std::vectors
    std::vector<double> prizes_vec(std::begin(prizes), std::end(prizes));
    std::vector<double> costs_vec(std::begin(costs), std::end(costs));
    std::vector<int> expected_node_result_vec(std::begin(expected_node_result), std::end(expected_node_result));
    std::vector<int> expected_edge_result_vec{}; // Empty vector for edges

    // Call the main helper function
    RunAlgo(edges, prizes_vec, costs_vec, root, target_num_active_clusters, pruning,
            expected_node_result_vec, expected_edge_result_vec);
}


// --- Test Cases ---
// The bodies of the individual TEST macros remain unchanged.
// They define the specific inputs and expected outputs, relying on the
// RunAlgo helpers to handle the interaction with the PCSTFast class.

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
    const int node_result[] = {0}; // Expected result only contains the root

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result); // Use overload expecting no edges
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
                            0, // Node 3 has 0 prize
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
    ca::PCSTFast::IndexType root = 3; // Root is node 3
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kGWPruning;

    const int node_result[] = {3, 4, 6, 7, 8};
    const int edge_result[] = {16, 20, 21, 23}; // Edges: (3,7), (4,6), (4,7), (6,8)

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
                            0, // Node 3 has 0 prize
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
    ca::PCSTFast::IndexType root = 3; // Root is node 3
    int target_num_active_clusters = 0;
    auto pruning = ca::PCSTFast::PruningMethod::kStrongPruning;
    const int node_result[] = {3}; // Only the root expected after strong pruning

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
    const int node_result[] = {7}; // User's updated expectation

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
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
    const int node_result[] = {1}; // User's updated expectation

    RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
            node_result);
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
