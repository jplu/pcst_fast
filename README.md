[![Build Status](https://api.travis-ci.org/fraenkel-lab/pcst_fast.svg)](https://travis-ci.org/fraenkel-lab/pcst_fast)

pcst_fast
=========

A library for solving the **prize-collecting Steiner forest (PCSF)** problem on graphs.
The underlying algorithm is based on the classical Goemans-Williamson approximation scheme.
Our variant provably runs in nearly-linear time and has a factor-2 approximation guarantee.
The following paper contains details about the algorithm:

> [A Nearly-Linear Time Framework for Graph-Structured Sparsity](http://people.csail.mit.edu/ludwigs/papers/icml15_graphsparsity.pdf)
> Chinmay Hegde, Piotr Indyk, Ludwig Schmidt
> ICML 2015

Comparaison with the original version
-------------------------------------
Here is a detailed report comparing the old and new versions of the codebase, focusing on the improvements made in the newer version:

**Analysis of Improvements: Old vs. New Codebase**

The new codebase represents a significant refactoring and modernization effort compared to the old version. The improvements span performance, code structure, maintainability, robustness, and API design.

**I. Performance Improvements**

1.  **`PriorityQueue` Implementation:**
    *   **Old:** Used `std::set<std::pair<ValueType, IndexType>>` as the underlying priority queue implementation, with a `std::vector` mapping indices to set iterators.
    *   **New:** Implements a proper **binary heap** using `std::vector` for heap data and another `std::vector` for direct index-to-position mapping (`index_to_heap_pos_`).
    *   **Improvement:** This is a *major* performance enhancement. The `std::set`-based approach has O(log N) complexity for most operations due to the balanced tree structure *and* significant overhead due to non-contiguous memory allocation (poor cache locality). The new binary heap provides the same theoretical O(log N) for inserts/updates/deletes but with much lower constant factors and significantly better cache performance due to the contiguous `std::vector` storage. O(1) lookup for element positions is maintained.

2.  **`ConnectFinalPruner` Dijkstra Implementation:**
    *   **New Pruning**: `ConnectFinalPruner` that replaces the old `GWPruning` that had an issue into its implementation. So to reproduce the exact same approach, this new pruning has been created.
    *   **Implementation:** Use the Dijkstra search with an expected complexity (e.g., O(E' log V) or O(E' + V log V) depending on PQ implementation, where E' is the number of edges in the relevant graph component), instead of potentially much slower complexity. This is crucial for the `ConnectFinalComponents` pruning strategy.

3.  **`PairingHeap` Memory Management:**
    *   **Old:** Relied on raw `new`/`delete` and passed an external `std::vector* buffer` for managing nodes during deletion and memory release. This buffer management was complex and potentially inefficient.
    *   **New:** Implements an **internal free list** (`allocate_node`, `reuse_node`, `deallocate_node`). Nodes are recycled, reducing the overhead of calls to the global `new` and `delete` operators and potentially improving cache locality. The external buffer dependency is removed, making the class self-contained.
    *   **Improvement:** Reduces memory allocation overhead and potentially improves cache performance during heap operations, especially during intense merging phases. Increases encapsulation.

**II. Code Structure & Maintainability Improvements**

1.  **Centralized Logging (`logger.h`, `logger.cc`):**
    *   **Old:** Used a raw C-style function pointer (`void (*output_function_)(const char*)`) passed through constructors. Logging calls involved manual `snprintf` into a fixed buffer within `PCSTFast`. Log level checking was manual (`if (verbosity_level >= X)`).
    *   **New:** Introduces a dedicated `Logger` class with a `LogLevel` enum. Uses `std::function<void(LogLevel, const std::string&)>` (`LogSink`) for flexible output redirection. Centralizes formatting and level checking. Provides basic thread safety via a mutex.
    *   **Improvement:** Dramatically improves structure, maintainability, and flexibility. Logging is type-safe, decoupled from `PCSTFast`, easily configurable (e.g., directing logs to a file, network, or Python as shown in pybind), and less error-prone than manual `snprintf`.

2.  **Pruning Strategy Pattern (`ipruner.h`, `pruner_factory.h`, specific pruners):**
    *   **Old:** Pruning logic (`Simple`, `GW`, `Strong`) was mixed within `pcst_fast.cc` with conditional logic based on the `PruningMethod` enum.
    *   **New:** Defines an `IPruner` interface. Each pruning strategy is implemented as a separate class inheriting from `IPruner`. A factory function (`create_pruner`) instantiates the correct strategy. `AdvancedPrunerBase` encapsulates common logic for GW/Strong pruning.
    *   **Improvement:** Excellent separation of concerns (Strategy Pattern). Makes the code modular, easier to understand, maintain, test, and extend with new pruning methods without modifying `PCSTFast::run`.

3.  **`PruningContext` (`pruning_context.h`):**
    *   **Old:** Pruning functions within `PCSTFast` likely accessed many member variables directly, or functions took numerous arguments.
    *   **New:** Introduces a `PruningContext` struct passed to `IPruner::prune`. This struct bundles all necessary data (references to edges, costs, prizes, clusters, node flags, logger, etc.).
    *   **Improvement:** Encapsulates data dependencies, simplifies pruner function signatures, and improves readability and maintainability by making data flow explicit.

4.  **`PairingHeap` Encapsulation:**
    *   **Old:** Dependency on the external `buffer` vector for core operations.
    *   **New:** Fully self-contained due to the internal free list.
    *   **Improvement:** Better encapsulation and easier use of the `PairingHeap` class.

5.  **Removal of Helper Members from `PCSTFast`:**
    *   **Old:** `PCSTFast` contained many intermediate vectors used by specific pruning methods (e.g., `path_compression_visited`, `cluster_queue`, `phase3_neighbors`, strong pruning vectors).
    *   **New:** These vectors are now typically members of the specific pruner classes (`AdvancedPrunerBase`, `StrongPruner`) or are local variables within methods, reducing the state managed directly by `PCSTFast`.
    *   **Improvement:** Better encapsulation, reduces the size and complexity of the main `PCSTFast` class.

**III. Correctness & Robustness Improvements**

1.  **Input Validation (`pcst_fast_pybind.cc`):**
    *   **Old:** Pybind wrapper had minimal validation mentioned.
    *   **New:** The pybind wrapper includes extensive checks for array shapes, dimensions, data types, index ranges (nodes in edges, root index), and non-negative prizes/costs. It also uses `safe_cast_int64_to_int` for index conversion.
    *   **Improvement:** Makes the Python interface much more robust, providing clearer error messages to the user instead of potentially causing crashes or incorrect results in the C++ core.

2.  **Error Handling:**
    *   **Old:** Error reporting primarily through the output function pointer, potentially mixing errors with verbose logging. Relied on `std::invalid_argument` for some setup errors.
    *   **New:** Uses the structured `Logger` with distinct `LogLevels` (FATAL, ERROR, WARNING). Input validation in pybind throws specific Python exceptions (`ValueError`, `IndexError`). `PCSTFast::run` returns `bool` and relies on logging for failure details, while the pybind wrapper converts this to `RuntimeError`.
    *   **Improvement:** More structured and informative error handling and reporting.

3.  **Type Safety:**
    *   **Old:** Use of `int` for indices and potentially raw pointers could be less safe. Logging used C-style variadic functions.
    *   **New:** Consistent use of defined types (`IndexType`, `ValueType`). Introduction of `LogLevel` enum for logging. Use of `std::unique_ptr` for `pruner_`. Use of references (`Logger&`).
    *   **Improvement:** Increased type safety and reduced risk of certain classes of errors.

4.  **RAII and Resource Management:**
    *   **Old:** `PairingHeap` had complex manual memory management logic in `release_memory` tied to the external buffer. `PCSTFast` destructor manually called `release_memory` on clusters.
    *   **New:** `PairingHeap` manages its own nodes via the free list and its destructor handles cleanup. `PCSTFast` destructor is simpler due to `unique_ptr<IPruner>` and self-managing `Cluster` members (assuming `Cluster`'s heap handles its own memory).
    *   **Improvement:** Better adherence to RAII principles, making resource management more automatic and less error-prone.

**IV. API & Usability Improvements**

1.  **Python API (`pcst_fast_pybind.cc`):**
    *   **Old:** Basic wrapper.
    *   **New:** Clearer docstrings, explicit argument names (`py::arg`), better error handling and validation (as mentioned above), consistent use of `int64_t` for Python array indices. Handles rooted/unrooted `num_clusters` logic more cleanly.
    *   **Improvement:** More user-friendly and robust Python interface.

2.  **`PairingHeap` API:**
    *   **Old:** `decrease_key` required `from_value`. Constructor required external buffer.
    *   **New:** `decrease_key` only requires the target node and the `new_value`. Default constructor works. Uses `std::optional` for `get_min`. Defines `kInvalidHandle`.
    *   **Improvement:** Cleaner and more modern API.

3.  **Logging Configuration:**
    *   **Old:** Configuration limited to a single function pointer and an integer verbosity level.
    *   **New:** Configuration via `LogSink` allows arbitrary output mechanisms. Log level managed by the `Logger` object.
    *   **Improvement:** More flexible and powerful logging configuration.

**V. Modern C++ Features**

1.  **Concepts (`pairing_heap.h`, `priority_queue.h`):** Used to constrain template parameters (`PairingHeapValue`, `Indexable`), improving template error messages and code clarity.
2.  **`noexcept`:** Applied more liberally to methods where applicable, aiding compiler optimization and indicating intent.
3.  **`[[nodiscard]]`, `[[no_unique_address]]`:** Used to provide hints to the compiler and users about function usage and potential memory layout optimizations.
4.  **RAII:** Better adherence as noted in Robustness.
5.  **`std::optional`:** Used for return types where a value might not exist (`PairingHeap::get_min`, `PriorityQueue::get_min`).
6.  **Move Semantics:** Implemented for `PairingHeap` (new version). `PriorityQueue` (new version) is move-enabled by default due to `std::vector`.
7.  **Range-based for loops, `std::ranges::sort`:** Used where appropriate (e.g., in tests).
8.  **`enum class`:** Used for `LogLevel` and `PruningMethod` for better type safety.

Installation
------------

- **manual**

	The core C++ library has no dependencies besides a basic build system for C++11.
	Both g++ and clang are currently supported.
	The Python wrapper requires a functioning Python build system.

	Clone the repository and compile the python wrapper with the supplied makefile:

	    make pcst_fast_py

	You can then import the package via `import pcst_fast`.

Usage
-----

The `pcst_fast` package contains the following function:

    vertices, edges = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)

The parameters are:
* `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
* `prizes`: the node prizes as a 1D float64 array.
* `costs`: the edge costs as a 1D float64 array.
* `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
* `num_clusters`: the number of connected components in the output.
* `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive). `'none'` and `'simple'` return intermediate stages of the algorithm and do not have approximation guarantees. They are only intended for development. The standard GW pruning method is `'gw'`, which is also the default. `'strong'` uses "strong pruning", which was introduced in [\[JMP00\]](http://dl.acm.org/citation.cfm?id=338637). It has the same theoretical guarantees as GW pruning but better empirical performance in some cases. For the PCSF problem, the output of strong pruning is at least as good as the output of GW pruning.
* `verbosity_level`: an integer indicating how much debug output the function should produce.

The output variables are:
* `vertices`: the vertices in the solution as a 1D int64 array.
* `edges`: the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.

Performance
-----------

The following paper contains many results on standard PCST benchmark instances:

> [A Fast, Adaptive Variant of the Goemans-Williamson Scheme for the Prize-Collecting Steiner Tree Problem](http://people.csail.mit.edu/ludwigs/papers/dimacs14_fastpcst.pdf)
> Chinmay Hegde, Piotr Indyk, Ludwig Schmidt
> Workshop of the 11th DIMACS Implementation Challenge: Steiner Tree Problems, 2014

On instances with up to 350,000 edges, the algorithm typically runs in under 2 seconds on a standard laptop computer from 2010.

