#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

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
#include <map>


namespace py = pybind11;
using namespace cluster_approx;

LogLevel map_verbosity_to_log_level(int verbosity_level) {
    if (verbosity_level < 0) verbosity_level = 0;
    if (verbosity_level > 5) verbosity_level = 5;

    static const std::map<int, LogLevel> level_map = {
        {0, LogLevel::FATAL},
        {1, LogLevel::ERROR},
        {2, LogLevel::WARNING},
        {3, LogLevel::INFO},
        {4, LogLevel::DEBUG},
        {5, LogLevel::TRACE}
    };
    return level_map.at(verbosity_level);
}


std::pair<py::array_t<NodeId>, py::array_t<EdgeId>> pcst_fast(
    py::array_t<NodeId, py::array::c_style | py::array::forcecast> edges,
    py::array_t<double, py::array::c_style | py::array::forcecast> prizes,
    py::array_t<double, py::array::c_style | py::array::forcecast> costs,
    NodeId root,
    int num_clusters,
    const std::string& pruning_method_str,
    int verbosity_level) {
    py::buffer_info edges_info = edges.request();
    py::buffer_info prizes_info = prizes.request();
    py::buffer_info costs_info = costs.request();

    if (edges_info.ndim != 2 || edges_info.shape[1] != 2) {
        throw std::invalid_argument("Edges array must be a 2D array with shape (num_edges, 2).");
    }
    if (prizes_info.ndim != 1) {
        throw std::invalid_argument("Prizes array must be a 1D array.");
    }
    if (costs_info.ndim != 1) {
        throw std::invalid_argument("Costs array must be a 1D array.");
    }

    size_t num_edges = edges_info.shape[0];
    size_t num_nodes = prizes_info.shape[0];

    if (static_cast<size_t>(costs_info.shape[0]) != num_edges) {
        throw std::invalid_argument("Number of costs must match the number of edges.");
    }

    if (root != kInvalidNodeId && (root < 0 || static_cast<size_t>(root) >= num_nodes)) {
        throw std::out_of_range("Root node index " + std::to_string(root) +
                                " is out of range [0, " + std::to_string(num_nodes) + ").");
    }

    int target_num_active_clusters = num_clusters;
    if (root != kInvalidNodeId) {
        if (num_clusters != 1) {
            throw std::invalid_argument("For rooted problems (root != -1), num_clusters must be 1.");
        }
        target_num_active_clusters = 0;
    } else {
        if (num_clusters < 1) {
            throw std::invalid_argument("For unrooted problems (root = -1), num_clusters must be at least 1.");
        }
    }

    PruningMethod pruning_method = parse_pruning_method(pruning_method_str);
    if (pruning_method == PruningMethod::kUnknown) {
        throw std::invalid_argument("Unknown pruning method: " + pruning_method_str +
                                    ". Valid options are: 'none', 'simple', 'gw', 'strong'.");
    }


    StderrLogger logger(map_verbosity_to_log_level(verbosity_level));
    logger.log(LogLevel::INFO, "pcst_fast pybind called. Root: {}, Target Clusters: {}, Pruning: {}, Verbosity: {}",
               root, num_clusters, pruning_method_str, verbosity_level);


    auto* edges_ptr = static_cast<NodeId*>(edges_info.ptr);
    auto* prizes_ptr = static_cast<double*>(prizes_info.ptr);
    auto* costs_ptr = static_cast<double*>(costs_info.ptr);

    auto edges_span = std::span<const std::pair<NodeId, NodeId>>(
                          reinterpret_cast<const std::pair<NodeId, NodeId>*>(edges_ptr), num_edges);
    auto prizes_span = std::span<const double>(prizes_ptr, num_nodes);
    auto costs_span = std::span<const double>(costs_ptr, num_edges);

    GraphData graph {
        .edges = edges_span,
        .prizes = prizes_span,
        .costs = costs_span,
        .root = root
    };


    PCSTCoreAlgorithm core_algo(graph, target_num_active_clusters, &logger);
    CoreAlgorithmResult core_result = core_algo.run();


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
    case PruningMethod::kUnknown:
    default:
        throw std::logic_error("Invalid pruning method reached switch statement.");
    }

    PruningInput pruning_input {
        .graph = graph,
        .core_result = core_result,
        .logger = &logger
    };

    logger.log(LogLevel::INFO, "Core algorithm finished. Running {} pruner.", pruning_method_str);
    PruningResult final_result = pruner->prune(pruning_input);


    logger.log(LogLevel::INFO, "Pruning finished. Result: {} nodes, {} edges.",
               final_result.nodes.size(), final_result.edges.size());

    py::array_t<NodeId> result_nodes_array(final_result.nodes.size());
    py::buffer_info result_nodes_info = result_nodes_array.request();
    NodeId* result_nodes_ptr = static_cast<NodeId*>(result_nodes_info.ptr);
    std::copy(final_result.nodes.begin(), final_result.nodes.end(), result_nodes_ptr);

    py::array_t<EdgeId> result_edges_array(final_result.edges.size());
    py::buffer_info result_edges_info = result_edges_array.request();
    EdgeId* result_edges_ptr = static_cast<EdgeId*>(result_edges_info.ptr);
    std::copy(final_result.edges.begin(), final_result.edges.end(), result_edges_ptr);

    return std::make_pair(result_nodes_array, result_edges_array);
}


PYBIND11_MODULE(pcst_fast, m) {
    m.doc() = R"pbdoc(
        Pybind11 bindings for the pcst_fast C++ library.

        Provides a fast implementation for solving the Prize-Collecting Steiner
        Forest (PCSF) problem, also known as the Prize-Collecting Steiner Tree (PCST)
        problem when rooted or seeking a single tree. Uses a growth-based algorithm
        with different pruning strategies.
    )pbdoc";

    m.def("pcst_fast",
          &pcst_fast,
          py::arg("edges"),
          py::arg("prizes"),
          py::arg("costs"),
          py::arg("root"),
          py::arg("num_clusters"),
          py::arg("pruning"),
          py::arg("verbosity_level") = -1,
          R"pbdoc(
        Runs the Prize-Collecting Steiner Forest algorithm (Fast C++ Implementation).
        Releases the Python Global Interpreter Lock (GIL) during C++ execution.

        Finds a forest (or tree if rooted) connecting subsets of terminals (nodes with
        positive prizes) that maximizes the total prize of connected terminals minus
        the total cost of edges used, subject to constraints on the number of trees
        and potentially requiring a specific root node.

        Args:
            edges (numpy.ndarray[int64]): Array of shape (num_edges, 2) listing undirected edges
                                           using 0-based node indices. Must be C-contiguous.
                                           Indices must fit within a 32-bit signed integer.
            prizes (numpy.ndarray[float64]): Array of shape (num_nodes,) listing non-negative node prizes.
                                             Must be C-contiguous.
            costs (numpy.ndarray[float64]): Array of shape (num_edges,) listing non-negative edge costs.
                                            Must be C-contiguous.
            root (int): The root node index for the rooted variant (PCSTree).
                        Use -1 or any negative value for the unrooted variant (PCSForest).
                        Must be less than num_nodes if non-negative.
            num_clusters (int): The target number of trees (connected components) in the output forest.
                                For the rooted variant (root >= 0), this is typically 1 (or 0 internally),
                                resulting in a single tree containing the root. If set > 1 for rooted,
                                a warning is printed and it's treated as 1 (0 internally).
                                For the unrooted variant (root < 0), this must be positive.
            pruning (str): The pruning method to apply after the main algorithm phase to potentially
                           improve the solution quality or enforce structure.
                           Options: "none", "simple", "gw", "strong", "connectfinal".
            verbosity_level (int, optional): Controls the maximum level of messages printed by the C++ module.
                                             Higher values show more detail. Defaults to 0 (FATAL).
                                              0: FATAL
                                              1: ERROR
                                              2: WARNING
                                              3: INFO
                                              4: DEBUG
                                              5: TRACE
                                             Output goes to Python's stdout/stderr via print().

        Returns:
            tuple[numpy.ndarray[int64], numpy.ndarray[int64]]: A pair containing:
                - nodes: A 1D numpy array of selected node indices (int64) present in the solution forest. Sorted.
                - edges: A 1D numpy array of selected edge indices (int64), corresponding to the
                         indices in the input "costs" and "edges" arrays, forming the solution forest. Sorted.

        Raises:
            ValueError: If input arrays have incorrect shapes, dimensions, or types; if num_clusters
                        is invalid for the rooted/unrooted case; if costs/prizes are negative;
                        if the pruning string is not recognized; or if input
                        indices are too large for internal 32-bit representation.
            IndexError: If node indices in "edges" or the "root" index are out of the valid range [0, num_nodes).
            RuntimeError: If the underlying C++ algorithm encounters an internal error or fails to run successfully. Check logged output.
        )pbdoc"
         );

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const std::invalid_argument &e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const std::runtime_error &e) {
            PyErr_SetString(PyExc_RuntimeError, e.what());
        } catch (const std::out_of_range &e) {
            PyErr_SetString(PyExc_IndexError, e.what());
        } catch (const std::exception &e) {
            PyErr_SetString(PyExc_Exception, e.what());
        }
    });
}