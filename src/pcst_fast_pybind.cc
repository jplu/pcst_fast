#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>
#include <limits>
#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pcst_fast.h"

namespace py = pybind11;

using cluster_approx::PCSTFast;

namespace {
    void output_function(const char* output) {
        py::print(output, py::arg("flush") = true);
    }

    int safe_cast_int64_to_int(int64_t val) {
        if (val < std::numeric_limits<int>::min() || val > std::numeric_limits<int>::max()) {
            throw std::overflow_error("Edge index value cannot fit into an int.");
        }
        
        return static_cast<int>(val);
    }
}

std::pair<py::array_t<int64_t>, py::array_t<int64_t>> pcst_fast(
    py::array_t<int64_t, py::array::c_style | py::array::forcecast> edges,
    py::array_t<double, py::array::c_style | py::array::forcecast> prizes,
    py::array_t<double, py::array::c_style | py::array::forcecast> costs,
    int root,
    int num_clusters,
    const std::string& pruning,
    int verbosity_level) {
        int target_num_active_clusters = num_clusters;
        
        if (root >= 0) {
            if (num_clusters != 1 && num_clusters != 0) {
                throw std::invalid_argument("In the rooted case, num_clusters must be 1 (or 0 if internally handled).");
            }
            
            target_num_active_clusters = 0;
        } else {
            if (num_clusters <= 0) {
                throw std::invalid_argument("In the unrooted case, num_clusters must be positive.");
            }
        }

        py::buffer_info edges_info = edges.request();
        
        if (edges_info.ndim != 2 || edges_info.shape[1] != 2) {
            throw std::invalid_argument("Edges must be a 2D array with shape (num_edges, 2).");
        }
  
        const int64_t num_edges_long = edges_info.shape[0];
        
        if (num_edges_long > std::numeric_limits<int>::max()) {
            throw std::overflow_error("Number of edges exceeds INT_MAX.");
        }

        const int num_edges = static_cast<int>(num_edges_long);
        py::buffer_info prizes_info = prizes.request();
  
        if (prizes_info.ndim != 1) {
            throw std::invalid_argument("Prizes must be a 1D array.");
        }
   
        const int64_t num_nodes_long = prizes_info.shape[0];
        
        if (num_nodes_long > std::numeric_limits<int>::max()) {
            throw std::overflow_error("Number of nodes exceeds INT_MAX.");
        }
   
        const int num_nodes = static_cast<int>(num_nodes_long);
        py::buffer_info costs_info = costs.request();
  
        if (costs_info.ndim != 1) {
            throw std::invalid_argument("Costs must be a 1D array.");
        }
  
        if (static_cast<int>(costs_info.shape[0]) != num_edges) {
            throw std::invalid_argument("Number of costs must equal number of edges.");
        }

        if (root >= num_nodes) {
            throw std::out_of_range("Root node index out of range [0, num_nodes).");
        }
  
        int actual_root = (root < 0) ? -1 : root;
        std::vector<std::pair<int, int>> tmp_edges(num_edges);
        const int64_t* edges_ptr = static_cast<const int64_t*>(edges_info.ptr);
  
        for (int i = 0; i < num_edges; ++i) {
            int u = safe_cast_int64_to_int(edges_ptr[2 * i]);
            int v = safe_cast_int64_to_int(edges_ptr[2 * i + 1]);
            
            if (u < 0 || u >= num_nodes || v < 0 || v >= num_nodes) {
                throw std::out_of_range("Edge index out of valid node range [0, num_nodes).");
            }
            
            tmp_edges[i] = {u, v};
        }

        std::vector<double> tmp_prizes(num_nodes);
        const double* prizes_ptr = static_cast<const double*>(prizes_info.ptr);
        std::copy(prizes_ptr, prizes_ptr + num_nodes, tmp_prizes.begin());
        std::vector<double> tmp_costs(num_edges);
        const double* costs_ptr = static_cast<const double*>(costs_info.ptr);
        std::copy(costs_ptr, costs_ptr + num_edges, tmp_costs.begin());
        PCSTFast::PruningMethod pruning_method;
  
        try {
            pruning_method = PCSTFast::parse_pruning_method(pruning);
        } catch (const std::invalid_argument& e) {
            throw std::invalid_argument("Invalid pruning method string: " + pruning);
        }

        PCSTFast algo(tmp_edges, tmp_prizes, tmp_costs, actual_root,
                    target_num_active_clusters, pruning_method, verbosity_level,
                    output_function);

        std::vector<int> result_nodes;
        std::vector<int> result_edges;
        
        result_nodes.reserve(num_nodes);
        result_edges.reserve(num_edges);
        
        algo.run(&result_nodes, &result_edges);

        py::array_t<int64_t> result_nodes_array(result_nodes.size());
        py::array_t<int64_t> result_edges_array(result_edges.size());
        py::buffer_info result_nodes_info = result_nodes_array.request();
        py::buffer_info result_edges_info = result_edges_array.request();
        int64_t* result_nodes_ptr = static_cast<int64_t*>(result_nodes_info.ptr);
        int64_t* result_edges_ptr = static_cast<int64_t*>(result_edges_info.ptr);
        
        std::copy(result_nodes.begin(), result_nodes.end(), result_nodes_ptr);
        std::copy(result_edges.begin(), result_edges.end(), result_edges_ptr);

        return std::make_pair(std::move(result_nodes_array), std::move(result_edges_array));
}

PYBIND11_MODULE(pcst_fast, m) {
    m.doc() = "A fast algorithm for the PCSF problem.";

    m.def("pcst_fast", &pcst_fast,
        py::arg("edges"), py::arg("prizes"), py::arg("costs"),
        py::arg("root"), py::arg("num_clusters"), py::arg("pruning"),
        py::arg("verbosity_level") = 0,
        "Runs the Prize-Collecting Steiner Forest algorithm (Fast Implementation).\n\n"
        "Args:\n"
        "    edges (numpy.ndarray[int64]): Array of shape (num_edges, 2) listing edges.\n"
        "    prizes (numpy.ndarray[float64]): Array of shape (num_nodes,) listing node prizes.\n"
        "    costs (numpy.ndarray[float64]): Array of shape (num_edges,) listing edge costs.\n"
        "    root (int): The root node index, or -1 for unrooted.\n"
        "    num_clusters (int): The number of trees desired in the output forest.\n"
        "    pruning (str): The pruning method ('none', 'simple', 'gw', 'strong').\n"
        "    verbosity_level (int): Verbosity level (default: 0).\n\n"
        "Returns:\n"
        "    tuple[numpy.ndarray[int64], numpy.ndarray[int64]]: A pair containing:\n"
        "        - Array of selected node indices.\n"
        "        - Array of selected edge indices (indices into the input 'costs'/'edges' array)."
       );
}
