#include "pcst_fast.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>


namespace py = pybind11;
namespace ca = cluster_approx;

namespace {
    void output_function(const char* output) {
        py::print(output, py::arg("flush") = true);
    }

    int safe_cast_int64_to_int(int64_t val) {
        if (std::cmp_less(val, std::numeric_limits<int>::min()) ||
            std::cmp_greater(val, std::numeric_limits<int>::max())) {
                throw py::value_error(
                    "Input index value " + std::to_string(val) +
                    " cannot fit into the internal 32-bit integer type."
                );
        }

        return static_cast<int>(val);
    }
}


std::pair<py::array_t<int64_t>, py::array_t<int64_t>> pcst_fast_pybind(
    py::array_t<int64_t, py::array::c_style | py::array::forcecast> edges_py,
    py::array_t<double, py::array::c_style | py::array::forcecast> prizes_py,
    py::array_t<double, py::array::c_style | py::array::forcecast> costs_py,
    int root,
    int num_clusters,
    std::string_view pruning_sv,
    int verbosity_level) {
        int internal_target_num_active_clusters = num_clusters;

        if (root >= 0) {
            if (num_clusters != 1 && num_clusters != 0) {
                if (verbosity_level >= 0) {
                    py::print("Warning: num_clusters parameter is typically 1 or 0 for rooted PCST.", py::arg("flush")=true);
                }
            }

            internal_target_num_active_clusters = 0;
        } else {
            if (num_clusters <= 0) {
                throw py::value_error("In the unrooted case, num_clusters must be positive.");
            }
            
            internal_target_num_active_clusters = num_clusters;
        }

        py::buffer_info edges_info = edges_py.request();
        py::buffer_info prizes_info = prizes_py.request();
        py::buffer_info costs_info = costs_py.request();

        if (edges_info.ndim != 2 || edges_info.shape[1] != 2) {
            throw py::value_error("Edges must be a 2D array with shape (num_edges, 2).");
        }

        const auto num_edges_ssize = edges_info.shape[0];
        
        if (std::cmp_greater(num_edges_ssize, std::numeric_limits<int>::max())) {
             throw py::value_error("Number of edges exceeds internal integer limits.");
        }
        
        const int num_edges = static_cast<int>(num_edges_ssize);

        if (prizes_info.ndim != 1) {
            throw py::value_error("Prizes must be a 1D array.");
        }
        
        const auto num_nodes_ssize = prizes_info.shape[0];
        
        if (std::cmp_greater(num_nodes_ssize, std::numeric_limits<int>::max())) {
            throw py::value_error("Number of nodes exceeds internal integer limits.");
        }
        
        const int num_nodes = static_cast<int>(num_nodes_ssize);

        if (costs_info.ndim != 1) {
            throw py::value_error("Costs must be a 1D array.");
        }
        
        if (static_cast<int>(costs_info.shape[0]) != num_edges) {
             throw py::value_error("Number of costs must equal number of edges.");
        }
        
        if (root >= num_nodes) {
            throw py::index_error("Root node index out of range [0, num_nodes).");
        }
        
        int actual_root = (root < 0) ? ca::PCSTFast::kNoRoot : root;
        std::vector<std::pair<int, int>> tmp_edges(num_edges);
        const int64_t* edges_ptr = static_cast<const int64_t*>(edges_info.ptr);

        for (int i = 0; i < num_edges; ++i) {
            int u = safe_cast_int64_to_int(edges_ptr[2 * i]);
            int v = safe_cast_int64_to_int(edges_ptr[2 * i + 1]);
        
            if (u < 0 || u >= num_nodes || v < 0 || v >= num_nodes) {
                 throw py::index_error("Edge index out of valid node range [0, num_nodes).");
            }
        
            tmp_edges[i] = {u, v};
        }

        std::vector<double> tmp_prizes(num_nodes);
        const double* prizes_ptr = static_cast<const double*>(prizes_info.ptr);
        std::copy(prizes_ptr, prizes_ptr + num_nodes, tmp_prizes.begin());
        std::vector<double> tmp_costs(num_edges);
        const double* costs_ptr = static_cast<const double*>(costs_info.ptr);
        std::copy(costs_ptr, costs_ptr + num_edges, tmp_costs.begin());
        ca::PCSTFast::PruningMethod pruning_method = ca::PCSTFast::parse_pruning_method(pruning_sv);
        
        if (pruning_method == ca::PCSTFast::PruningMethod::kUnknownPruning) {
             throw py::value_error(std::string("Invalid pruning method string: ") + std::string(pruning_sv));
        }

        ca::PCSTFast algo(tmp_edges, tmp_prizes, tmp_costs, actual_root,
                          internal_target_num_active_clusters, pruning_method, verbosity_level,
                          output_function);

        std::vector<int> result_nodes;
        std::vector<int> result_edges;

        result_nodes.reserve(num_nodes);
        result_edges.reserve(num_edges);

        bool success = algo.run(&result_nodes, &result_edges);

        if (!success) {
            throw std::runtime_error("PCSTFast algorithm run failed. Check logs or inputs.");
        }

        py::array_t<int64_t> result_nodes_array(static_cast<py::ssize_t>(result_nodes.size()));
        py::array_t<int64_t> result_edges_array(static_cast<py::ssize_t>(result_edges.size()));
        py::buffer_info result_nodes_info = result_nodes_array.request();
        py::buffer_info result_edges_info = result_edges_array.request();
        int64_t* result_nodes_ptr = static_cast<int64_t*>(result_nodes_info.ptr);
        int64_t* result_edges_ptr = static_cast<int64_t*>(result_edges_info.ptr);

        std::transform(result_nodes.begin(), result_nodes.end(), result_nodes_ptr,
                       [](int val) { return static_cast<int64_t>(val); });
        std::transform(result_edges.begin(), result_edges.end(), result_edges_ptr,
                       [](int val) { return static_cast<int64_t>(val); });

        return std::make_pair(std::move(result_nodes_array), std::move(result_edges_array));
}

PYBIND11_MODULE(pcst_fast, m) {
    m.doc() = "Fast C++ implementation for the Prize-Collecting Steiner Forest (PCST/PCSF) problem.";

    m.def("pcst_fast",
          &pcst_fast_pybind,
          py::arg("edges"),
          py::arg("prizes"),
          py::arg("costs"),
          py::arg("root"),
          py::arg("num_clusters"),
          py::arg("pruning"),
          py::arg("verbosity_level") = 0,
        R"pbdoc(
        Runs the Prize-Collecting Steiner Forest algorithm (Fast C++ Implementation).

        Finds a forest connecting subsets of terminals that maximizes the sum of prizes
        of connected terminals minus the sum of costs of edges in the forest, subject
        to constraints on the number of trees and potentially a root node.

        Args:
            edges (numpy.ndarray[int64]): Array of shape (num_edges, 2) listing undirected edges
                                           using 0-based node indices.
            prizes (numpy.ndarray[float64]): Array of shape (num_nodes,) listing non-negative node prizes.
            costs (numpy.ndarray[float64]): Array of shape (num_edges,) listing non-negative edge costs.
            root (int): The root node index for the rooted variant (PCSF Tree).
                        Use -1 or any negative value for the unrooted variant (PCSF Forest).
            num_clusters (int): The target number of trees in the output forest.
                                For the rooted variant (root >= 0), this should be 1 (or 0).
                                For the unrooted variant (root < 0), this must be positive.
            pruning (str): The pruning method to apply after the main algorithm phase.
                           Options: 'none', 'simple', 'gw', 'strong'.
            verbosity_level (int): Controls the amount of console output printed by the C++ module
                                   (0: silent, 1: basic info, 2: detailed debug). Defaults to 0.

        Returns:
            tuple[numpy.ndarray[int64], numpy.ndarray[int64]]: A pair containing:
                - nodes: A 1D numpy array of selected node indices (int64).
                - edges: A 1D numpy array of selected edge indices (int64), corresponding to the
                         indices in the input 'costs' and 'edges' arrays.

        Raises:
            ValueError: If input array shapes, dimensions, or parameters are invalid
                        (e.g., negative costs/prizes, inconsistent sizes, invalid pruning string).
            IndexError: If node indices in `edges` or `root` are out of the valid range [0, num_nodes).
            OverflowError: If input indices or counts exceed internal C++ int limits (usually 2^31 - 1).
            RuntimeError: If the underlying C++ algorithm encounters an internal error during execution.
        )pbdoc"
       );
}
