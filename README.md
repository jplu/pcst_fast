# pcst_fast (Fork)

*This is a fork of the original [fraenkel-lab/pcst_fast](https://github.com/fraenkel-lab/pcst_fast).*

A fast C++ implementation of the Prize-Collecting Steiner Forest (PCSF) algorithm with Python bindings. This library allows finding optimal subgraphs (trees or forests) connecting nodes with prizes while minimizing edge costs.

## Features

*   Fast C++ backend for high performance.
*   Implements both rooted (Prize-Collecting Steiner Tree - PCST) and unrooted (Prize-Collecting Steiner Forest - PCSF) variants.
*   Supports various pruning methods to refine solutions.
*   Releases the Python Global Interpreter Lock (GIL) during C++ execution, allowing for better concurrency in threaded applications.
*   Simple Python API using NumPy arrays.
*	See the [NOTES](NOTES.md) file to have a better understanding of what has changed and still remaining to update.

## Installation

### Prerequisites

*   **C++ Compiler:** A C++23 compatible compiler:
    *   GCC version 13+
    *   Clang version 15+
    *   MSVC (Visual Studio 2022 v17.4+) - *Required if installing via `pip` on Windows*
*   **Python:** Version 3.9+
*   **Build System:**
    *   `make` (if building manually using the Makefile)
    *   Standard Python packaging tools (`pip`, `setuptools`, `build`)

### Methods

1.  **Install using pip (Recommended):**
    This command will compile the C++ extension and install the Python package.
    ```bash
    pip install .
    ```
    *Alternatively, install directly from GitHub:*
    ```bash
    pip install git+https://github.com/jplu/pcst_fast.git
    ```

2.  **Build from Source (Manual):**
    Clone the repository and use the Makefile to build the Python bindings.
    ```bash
    git clone https://github.com/jplu/pcst_fast.git
    cd pcst_fast
    make python_binding
    # The package can now be imported from within this directory
    # or installed into your environment using: pip install .
    ```

### Importing the Package

```python
import pcst_fast
```

## Usage

The core functionality is provided by the `pcst_fast` function.

### Function Definition

```python
vertices, edges = pcst_fast.pcst_fast(
    edges,
    prizes,
    costs,
    root,
    num_clusters,
    pruning,
    verbosity_level=0
)
```

### Description

This function executes the Prize-Collecting Steiner Forest algorithm. It aims to find a subgraph (a forest, or a single tree if rooted) that connects a subset of nodes (terminals, those with positive prizes). The objective is to maximize the total prize of the nodes included in the subgraph minus the total cost of the edges used to connect them.

Constraints can be applied, such as requiring a specific `root` node (PCST variant) or targeting a specific `num_clusters` (number of connected components/trees) in the resulting forest (PCSF variant).

### Parameters

*   `edges` (`numpy.ndarray[np.int64]`):
    *   Shape: `(num_edges, 2)`
    *   Content: An array listing the undirected edges. Each row represents an edge defined by the 0-based indices of the two nodes it connects.
    *   Requirements: Must be C-contiguous. Node indices must be non-negative and fit within a 32-bit signed integer (`< 2^31`).
*   `prizes` (`numpy.ndarray[np.float64]`):
    *   Shape: `(num_nodes,)`
    *   Content: An array where each element `prizes[i]` is the non-negative prize associated with node `i`.
    *   Requirements: Must be C-contiguous. Length determines the total number of nodes (`num_nodes`) in the graph.
*   `costs` (`numpy.ndarray[np.float64]`):
    *   Shape: `(num_edges,)`
    *   Content: An array where each element `costs[j]` is the non-negative cost associated with the `j`-th edge listed in the `edges` array.
    *   Requirements: Must be C-contiguous. Must have the same length as the first dimension of `edges`.
*   `root` (`int`):
    *   Specifies the root node for the **rooted** variant (PCST). The resulting subgraph will be a single tree containing this node.
    *   Use `-1` or any negative value to run the **unrooted** variant (PCSF). The result will be a forest.
    *   Requirements: If non-negative, must be a valid node index (`0 <= root < num_nodes`).
*   `num_clusters` (`int`):
    *   Target number of trees (connected components) in the output forest.
    *   **Rooted case (`root >= 0`):** This parameter is effectively ignored (internally treated as 1), as the algorithm seeks a single tree containing the `root`. A warning may be issued if `num_clusters > 1` is provided.
    *   **Unrooted case (`root < 0`):** Must be a positive integer (`>= 1`). Specifies the desired number of trees in the resulting forest.
*   `pruning` (`str`):
    *   Specifies the pruning method applied after the main algorithm phase to potentially improve the solution or enforce structural properties.
    *   Available options: `"none"`, `"simple"`, `"gw"` (Goemans-Williamson based), `"strong"`. *(Note: Based on original README, confirm if "connectfinal" is a valid option)*
*   `verbosity_level` (`int`, optional):
    *   Controls the verbosity of the C++ backend logging (output to Python's stdout/stderr via `print()`).
    *   Levels: `0` (FATAL), `1` (ERROR), `2` (WARNING), `3` (INFO), `4` (DEBUG), `5` (TRACE).
    *   Default: `0`.

### Returns

*   `tuple[numpy.ndarray[np.int64], numpy.ndarray[np.int64]]`:
    A tuple containing two NumPy arrays:
    1.  `nodes`: A 1D array (`np.int64`) containing the sorted indices of the nodes selected for the optimal forest/tree.
    2.  `edges`: A 1D array (`np.int64`) containing the sorted indices of the edges selected for the optimal forest/tree. These indices correspond to the rows in the input `edges` array and elements in the `costs` array.

### Raises

*   `ValueError`: If input arguments are invalid (e.g., incorrect array shapes, dimensions, types; negative costs/prizes; invalid `num_clusters` for the unrooted case; unrecognized `pruning` string; node indices too large for internal 32-bit representation).
*   `IndexError`: If node indices provided in `edges` or the `root` index are outside the valid range `[0, num_nodes)`.
*   `RuntimeError`: If the underlying C++ algorithm encounters an internal error. Check the logged output for more details based on the `verbosity_level`.

## Example

```python
import numpy as np
import pcst_fast

# --- Define Graph Data ---
# Graph:
#    0 --(cost=1)-- 1 --(cost=1)-- 2
#    |              |              |
# (cost=10)      (cost=3)          |
#    |------------- 3 --(cost=1)---|
#
# Nodes: 0, 1, 2, 3
# Edges: (0,1, c=1), (0,3, c=10), (1,2, c=1), (1,3, c=3), (2,3, c=1)
# Prizes: P(0)=5, P(1)=1, P(2)=5, P(3)=6

num_nodes = 4
edges = np.array([
    [0, 1], [0, 3], [1, 2], [1, 3], [2, 3]
], dtype=np.int64)

costs = np.array([
    1.0, 10.0, 1.0, 3.0, 1.0
], dtype=np.float64)

prizes = np.array([5.0, 1.0, 5.0, 6.0], dtype=np.float64) # Node prizes

# --- Run Rooted PCST ---
# Find a tree rooted at node 0, aiming to maximize (sum of prizes) - (sum of costs)
root_node = 0
clusters = 1 # Ignored for rooted, must be 1 tree
pruning_method = "strong"
verbose = 0 # Set higher (e.g., 3 for INFO) to see C++ logs

print("--- Running Rooted PCST (Root = 0) ---")
selected_nodes_rooted, selected_edge_indices_rooted = pcst_fast.pcst_fast(
    edges, prizes, costs, root_node, clusters, pruning_method, verbose
)

print(f"Selected Nodes: {selected_nodes_rooted}")
print(f"Selected Edge Indices: {selected_edge_indices_rooted}")

# --- Run Unrooted PCSF ---
# Find a forest with potentially multiple trees (e.g., target 2 clusters)
root_node = -1 # Unrooted
clusters = 2   # Target 2 trees
pruning_method = "strong"

print("\n--- Running Unrooted PCSF (Target Clusters = 2) ---")
selected_nodes_forest, selected_edge_indices_forest = pcst_fast.pcst_fast(
    edges, prizes, costs, root_node, clusters, pruning_method, verbose
)

print(f"Selected Nodes: {selected_nodes_forest}")
print(f"Selected Edge Indices: {selected_edge_indices_forest}")

# --- Example of Expected Error ---
print("\n--- Example: Triggering ValueError (Negative Prize) ---")
invalid_prizes = np.array([5.0, -1.0, 5.0, 6.0], dtype=np.float64)
try:
    pcst_fast.pcst_fast(edges, invalid_prizes, costs, -1, 1, "strong")
except ValueError as e:
    print(f"Caught expected error: {e}")

```
