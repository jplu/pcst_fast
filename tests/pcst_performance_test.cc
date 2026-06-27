#include <gtest/gtest.h>
#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pruning/no_pruner.h"
#include "pcst_fast/pruning/simple_pruner.h"
#include "pcst_fast/pruning/gw_pruner.h"
#include "pcst_fast/pruning/strong_pruner.h"

#include <chrono>
#include <random>
#include <iostream>
#include <vector>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cluster_approx {
namespace {

// Helper class to track and display benchmark timings
class Timer {
public:
    Timer(const std::string& name) : name_(name), start_(std::chrono::high_resolution_clock::now()) {}
    
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start_;
        std::cout << "[ BENCHMARK ] " << name_ << " took: " << duration.count() << " ms\n";
    }
private:
    std::string name_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

// Generates a 2D grid graph with randomized node prizes and edge costs
GraphData generate_grid_graph(
    int width, int height, 
    std::vector<std::pair<NodeId, NodeId>>& edges_storage,
    std::vector<double>& prizes_storage,
    std::vector<double>& costs_storage,
    unsigned int seed = 42) {
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prize_dist(0.0, 100.0);
    std::uniform_real_distribution<double> cost_dist(1.0, 50.0);
    std::bernoulli_distribution terminal_dist(0.3); // 30% of nodes are terminals

    int num_nodes = width * height;
    prizes_storage.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        prizes_storage[i] = terminal_dist(rng) ? prize_dist(rng) : 0.0;
    }

    edges_storage.clear();
    costs_storage.clear();
    edges_storage.reserve(2 * num_nodes);
    costs_storage.reserve(2 * num_nodes);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int u = y * width + x;
            if (x + 1 < width) {
                int v = y * width + (x + 1);
                edges_storage.push_back({u, v});
                costs_storage.push_back(cost_dist(rng));
            }
            if (y + 1 < height) {
                int v = (y + 1) * width + x;
                edges_storage.push_back({u, v});
                costs_storage.push_back(cost_dist(rng));
            }
        }
    }

    GraphData graph;
    graph.edges = edges_storage;
    graph.prizes = prizes_storage;
    graph.costs = costs_storage;
    graph.root = kInvalidNodeId;
    return graph;
}

} // namespace

class PCSTPerformanceTest : public ::testing::Test {
protected:
    void SetUp() override {
#ifdef _OPENMP
        int max_threads = omp_get_max_threads();
        std::cout << "[   INFO    ] OpenMP is ENABLED. Max threads available: " << max_threads << "\n";
#else
        std::cout << "[   INFO    ] OpenMP is DISABLED. Pruning will run sequentially.\n";
#endif
    }
};

TEST_F(PCSTPerformanceTest, UnrootedGridBenchmark_Large) {
    // 400x400 grid yields 160,000 nodes and roughly 319,200 edges
    const int width = 400;
    const int height = 400;

    std::vector<std::pair<NodeId, NodeId>> edges_storage;
    std::vector<double> prizes_storage;
    std::vector<double> costs_storage;

    std::cout << "[   INFO    ] Generating " << width << "x" << height << " grid graph...\n";
    GraphData graph = generate_grid_graph(width, height, edges_storage, prizes_storage, costs_storage);
    graph.root = kInvalidNodeId; // Unrooted forest setting

    StderrLogger logger(LogLevel::FATAL); // Disable log I/O overhead

    CoreAlgorithmResult core_result;
    {
        Timer timer("PCST Core Algorithm (Unrooted, 160k Nodes)");
        PCSTCoreAlgorithm core_algo(graph, 1, &logger);
        core_result = core_algo.run();
    }

    PruningInput pruning_input{
        .graph = graph,
        .core_result = core_result,
        .logger = &logger
    };

    // Evaluate NoPruner
    {
        Timer timer("No Pruning");
        pruning::NoPruner pruner;
        auto result = pruner.prune(pruning_input);
        EXPECT_FALSE(result.nodes.empty());
    }

    // Evaluate SimplePruner
    {
        Timer timer("Simple Pruning");
        pruning::SimplePruner pruner;
        auto result = pruner.prune(pruning_input);
        EXPECT_GE(result.nodes.size(), 0);
    }

    // Evaluate GWPruner
    {
        Timer timer("GW Pruning");
        pruning::GWPruner pruner;
        auto result = pruner.prune(pruning_input);
        EXPECT_GE(result.nodes.size(), 0);
    }

    // Evaluate StrongPruner (Exercises multi-threaded disjoint component traversal)
    {
        Timer timer("Strong Pruning (Parallel Component DFS)");
        pruning::StrongPruner pruner;
        auto result = pruner.prune(pruning_input);
        EXPECT_GE(result.nodes.size(), 0);
    }
}

TEST_F(PCSTPerformanceTest, RootedGridBenchmark_Large) {
    const int width = 400;
    const int height = 400;

    std::vector<std::pair<NodeId, NodeId>> edges_storage;
    std::vector<double> prizes_storage;
    std::vector<double> costs_storage;

    GraphData graph = generate_grid_graph(width, height, edges_storage, prizes_storage, costs_storage);
    graph.root = (width * height) / 2; // Set root to the approximate center

    StderrLogger logger(LogLevel::FATAL);

    CoreAlgorithmResult core_result;
    {
        Timer timer("PCST Core Algorithm (Rooted, 160k Nodes)");
        PCSTCoreAlgorithm core_algo(graph, 0, &logger);
        core_result = core_algo.run();
    }

    PruningInput pruning_input{
        .graph = graph,
        .core_result = core_result,
        .logger = &logger
    };

    {
        Timer timer("Strong Pruning (Rooted Tree)");
        pruning::StrongPruner pruner;
        auto result = pruner.prune(pruning_input);
        EXPECT_GE(result.nodes.size(), 0);
    }
}

}