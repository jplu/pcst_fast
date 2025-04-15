# NOTES

## **New Design & Structure Improvements:**

1.  **Separation of Concerns:**
    *   **Core vs. Pruning:** The core Goemans-Williamson algorithm is now isolated in `PCSTCoreAlgorithm`. Pruning strategies are implemented as separate classes (`NoPruner`, `SimplePruner`, `GWPruner`, `StrongPruner`) inheriting from a common `IPruner` interface.
    *   **Clear Interfaces:** Interactions between the core algorithm and pruners are defined by explicit data structures (`CoreAlgorithmResult`, `PruningInput`, `PruningResult`).

2.  **Enhanced Modularity:**
    *   Each component (core algorithm, individual pruners, logger, data structures) resides in its own dedicated files (`.h`/`.cc`).
    *   Adding new pruning strategies is significantly easier – just implement the `IPruner` interface without modifying the core algorithm or other pruners.

3.  **Localized State Management:**
    *   `PCSTCoreAlgorithm` manages only the state relevant to the core GW process.
    *   Each pruner class manages its *own* internal state (e.g., `GWPruner` handles `necessary` flags, `StrongPruner` handles its DFS stacks and payoff calculations). This reduces unnecessary memory usage when simpler pruners are selected.

4.  **Robust Error Handling:**
    *   Replaced boolean returns and manual printing with standard C++ exceptions (`std::invalid_argument`, `std::runtime_error`) for clearer error signaling and propagation, especially useful for the Python bindings.

5.  **Structured Logging:**
    *   Introduced a dedicated `Logger` class with defined severity levels (`FATAL` to `TRACE`).
    *   Logging is centralized through the `Logger` instance (created in Python bindings and passed down), allowing fine-grained control and separation from program logic.
    *   Extensive logging messages added at appropriate levels throughout the code provide deep visibility into the algorithm's execution.

6.  **Modern C++ Practices:**
    *   Systematically uses modern C++ features (up to C++23) like `#pragma once`, `nullptr`, `auto`, range-based for loops, `enum class`, `std::optional`, `std::span`, `std::format`, `[[nodiscard]]`, etc., improving code safety, readability, and expressiveness.

7.  **Organized File Structure:**
    *   A clear, hierarchical directory structure (`include/`, `src/`, `tests/`, `bindings/`) separates the public API, implementation details, tests, and language bindings, improving navigability and build system integration.

8.  **Granular Testing:**
    *   Includes dedicated test files for the core algorithm and *each* individual pruning strategy, alongside end-to-end integration tests, facilitating more targeted testing and easier debugging.

9.  **Explicit Invariants:**
    *   Added `assert()` statements to check crucial internal assumptions and invariants during development and debugging.


## **Encountered Issues:**

1. Core Algorithm Output Logic (Fixed)
    *   Issue: The initial refactoring produced different node/edge sets for NoPruning compared to the original code in some rooted cases (Simple2TestRootedNoPruning, Simple4TestRootedNoPruning).
    *   Root Cause: The original code's "No Pruning" option returned the unfiltered set of edges selected by the core algorithm (phase1_result) and constructed the node set by taking all endpoints of those edges plus any isolated nodes reachable from the root's final cluster (node_good). The refactored design initially filtered the edges based on endpoint node_good status before returning them, and the NoPruner simply returned these filtered results.
    *   Resolution:
        *   CoreAlgorithmResult was modified to store the unfiltered phase1_edges.
        *   NoPruner::prune was updated to replicate the original node set construction logic (endpoints + isolated good nodes) and return the unfiltered phase1_edges.
        *   SimplePruner, GWPruner, and StrongPruner were updated to perform the initial filtering of phase1_edges based on initial_node_filter (the node_good vector) as their starting point, matching the original code's use of phase2_result.

2. GW Pruning Logic (Partially Addressed)
    *   Issue 1: The GW pruning results differed in rooted cases where an edge connected to the root node via an Active-Inactive merge (Simple3TestRootedGWPruning).
        *   Root Cause 1: The original code appeared to always keep such edges connected to the root, deviating from the standard GW pruning rule (which would discard it if the root wasn't marked necessary by another path).
        *   Resolution 1: The refactored GWPruner was modified with a special case to not prune Active-Inactive merge edges where the inactive node is the specified root, mimicking the original behavior.
    *   Issue 2: The GW pruning results differed significantly in an unrooted case (Simple6TestUnRootedGWPruning). The original kept all edges, while the refactored code pruned several based on the standard necessity check.
        *   Root Cause 2: The original code's log/behavior suggested it was marking nodes/clusters as "necessary" more liberally than the standard GW pruning algorithm dictates (possibly keeping any AI edge where either endpoint was necessary, or using a different propagation).
        *   Resolution 2: This discrepancy was not "fixed" by altering the refactored code to match the original's seemingly non-standard behavior. The refactored GWPruner implements the standard logic (keep AI edge only if inactive side is necessary at the time of checking). The test failure for Simple6TestUnRootedGWPruning persists, reflecting this difference, likely indicating an issue or variation in the original implementation/test.

3. Core Algorithm Execution Path (Tie-Breaking) (Not "Fixed")
    *   Issue: In a complex rooted test (Medium1TestRootedGWPruning), the core algorithm produced a different set of selected phase1_edges compared to the original execution trace.
    *   Root Cause: This difference stemmed from the order in which events occurring at the exact same time (t=0.4) were processed. The standard C++ data structures used in the refactoring (std::set for PriorityQueue, PairingHeap's internal logic) resolved ties differently than the original implementation's data structures.
    *   Resolution: This was not "fixed" by forcing arbitrary tie-breaking. It's accepted as a potential, valid outcome of using different underlying data structure implementations while preserving the algorithm's high-level logic. The test still fails against the original expectation because the input to the pruning phase differs.



## **Remaining Improvements**

**I. Performance & Latency**

1.  **Data Structures:**
    *   **Priority Queues:** The `PriorityQueue` class (used for `clusters_deactivation_` and `clusters_next_edge_event_`) is based on `std::set`, providing O(log N) for insert, delete-min, and find (needed for delete/decrease-key). While flexible, specialized heap structures could be faster *if* certain operations dominate:
        *   If only insert/delete-min are critical, `std::priority_queue` (binary heap) offers better constant factors (O(log N) insert, O(log N) delete-min, but no efficient decrease-key/delete-element).
        *   If decrease-key is frequent and performance-critical, investigate alternatives like Fibonacci heaps or specialized d-ary heaps, though their complexity and constant factors can be higher in practice than simpler heaps.
    *   **Pairing Heap (`PairingHeap`)**: Known for good *amortized* performance, especially with decrease-key. However, `delete_min` can have higher worst-case cost due to the multi-pass merge. For specific usage patterns, a simpler binary heap or d-ary heap might offer more consistent performance, *if* decrease-key isn't the bottleneck. Replacing it would require careful reimplementation to match the offset logic.
    *   **Profiling Needed:** The actual bottleneck (core algorithm loop vs. pruning, heap vs. queue vs. traversals) should be determined through profiling on representative datasets before optimizing data structures.

2.  **Memory Allocation & Locality:**
    *   **`PairingHeap::Node` Allocation:** Nodes are allocated individually using `new`. For large heaps, this causes fragmentation and cache misses. Using a **pool allocator** or **arena allocator** for `Node` objects could significantly improve performance by reducing allocation overhead and improving data locality.
    *   **`std::vector` Reallocations:** Vectors like `clusters_`, `edge_parts_`, `inactive_merge_events_`, etc., can reallocate. While `reserve` is used in `initialize`, accurately predicting the final size (especially for `clusters_` which grows during merges) is hard. Over-reserving wastes memory, under-reserving causes copies. Custom allocators or alternative containers (like `std::deque` if locality isn't paramount but reallocations are costly) could be considered, but likely offer minor gains compared to node allocation.

3.  **Algorithm Internals:**
    *   **`get_sum_on_edge_part`:** Path compression helps, but it's still a traversal up the merge tree. If this becomes a bottleneck (unlikely?), further optimization is hard without changing the core cluster representation.
    *   **Pruning Traversals (DFS/BFS):** Strong and GW pruning rely on graph traversals. Standard library algorithms or adjacency list representations are generally efficient (O(N+M)). Minor gains might be possible via loop unrolling or cache-aware traversal patterns, but likely complex and low-impact.

4.  **Python Bindings:**
    *   **Output Conversion:** Copying results from `std::vector` to `py::array_t` incurs overhead. For very large results, investigate pybind11 techniques to return NumPy arrays that directly reference or take ownership of the C++ data (e.g., using `py::capsule` or custom buffer protocols), though this adds complexity.

**II. Memory Usage**

1.  **Data Structure Overheads:**
    *   `std::set` (in `PriorityQueue`) has per-node memory overhead (pointers, color bits). A vector-based heap uses less memory per element.
    *   `PairingHeap::Node` stores multiple pointers.
    *   `std::vector<PairingHeapType::ItemHandle>` in `PriorityQueue` adds memory for storing iterators/pointers.
    *   `Cluster` struct containing a `PairingHeap` member means significant memory per cluster.

2.  **Potential Reductions:**
    *   **Alternative Priority Queue:** Using a `std::priority_queue` (if feasible based on required operations) would reduce per-element memory overhead compared to `std::set`.
    *   **Heap Node Pooling:** As mentioned in performance, pool allocators also reduce memory fragmentation and potentially overall usage.
    *   **Graph Representation during Pruning:** Pruners currently build adjacency lists. If memory is extremely constrained, consider algorithms that work directly on the edge list, although this often increases algorithmic complexity.
    *   **Sparse Data Structures:** If the input graphs are typically very sparse, ensure underlying structures don't allocate excessive memory for non-existent nodes/edges (current `std::vector` usage based on max node ID is generally okay but could be refined).

**III. Software Design & Maintainability**

1.  **Current State (Good):** The refactoring significantly improved design by separating Core, Pruners, Data Structures, Logging, and Interfaces. Modularity is much better.
2.  **Potential Refinements:**
    *   **`PairingHeap` Buffer Dependency:** The external `shared_buffer` is awkward. The heap could potentially manage its own internal buffer (e.g., as a `thread_local` static member or a regular member) to simplify its API, at the cost of potentially more memory use if multiple heaps exist simultaneously without sharing.
    *   **Interface Header Coupling:** `pcst_interfaces.h` currently needs `#include "pcst_fast/pcst_core_internals.h"` because `CoreAlgorithmResult` contains `std::vector<Cluster>`. This exposes `Cluster` (and potentially other internals) more widely. A possible refinement (at the cost of complexity) would be to define an opaque `ClusterStateHandle` in the interface and have `CoreAlgorithmResult` return a container of handles, with the concrete `Cluster` objects managed internally by the Core/Pruner.
    *   **Pruner State Management:** Each pruner now manages its state. This is good, but ensure proper cleanup and minimal state duplication.
    *   **`pruning_utils.cc`:** Ensure this only contains genuinely shared utilities. Functions used by only one pruner should live in that pruner's `.cc` file.
3.  **GW Pruner Coupling & Testability:**
    *   **Observation:** The GW pruning logic in the original implementation appears tightly coupled to the internal state (`Cluster` objects and their `necessary` flags) generated during the core algorithm's execution. Replicating the exact output of the original code for all provided test cases required simulating this coupling (e.g., by passing and potentially modifying cluster state information between the core and the pruner).
    *   **Design Impact:** This coupling, inherited to maintain identical output, compromises the ideal separation of concerns between the core algorithm and the pruner. It makes unit-testing the `GWPruner` in isolation more difficult than other pruners, as it requires more complex state setup mimicking the core algorithm's output.
    *   **Potential Improvement:** A future redesign could implement a *cleaner* `GWPruner` that relies *only* on the edge list (`phase1_edges`), node filter (`initial_node_filter`), and merge event list (`inactive_merge_events`) provided in `CoreAlgorithmResult`. This pruner would manage its own necessity tracking based purely on node IDs. This would improve modularity and testability but might produce different (though potentially more theoretically standard) results on specific edge cases compared to the original implementation, requiring re-validation or adjustment of test expectations. The current implementation prioritizes output consistency over design purity for the GW pruner component.


**IV. Correctness & Robustness**

1.  **Test Discrepancies:** The known differences between refactored output and original test expectations (e.g., `Medium1Test`, `Simple6Test`) need resolution. Either the refactored code has a subtle bug deviating from the original algorithm, or the original code/tests had issues. This is the highest priority for correctness. Focus on:
    *   Tie-breaking in heaps/queues.
    *   Floating point comparisons (`eps_`).
    *   Specific logic in GW pruning (necessity check).
    *   Logic in `mark_nodes_as_good` (especially rooted case).
2.  **Floating Point Precision:** Using `double` and `eps_` comparisons works generally but can be fragile. For ultimate robustness in edge cases, consider using scaled integers or dedicated exact arithmetic libraries if precision issues are suspected.
3.  **Error Handling:** Use of `std::invalid_argument` and `std::runtime_error` is appropriate. Ensure all potential error conditions (invalid inputs, unexpected states) are covered. Constructor validation is now better.
4.  **Assertions:** Assertions help catch internal logic errors during development/debugging. Maintain them for critical invariants.
5.  **Resource Management:** Ensure RAII is used correctly (e.g., `std::vector` handles its memory). The `PairingHeap` cleanup relies on its owner (the `Cluster` vector) calling `release_memory` implicitly via the destructor, which seems okay.

**V. Usability & API**

1.  **Python Bindings:**
    *   The API seems clear.
    *   Logger integration could be enhanced by allowing users to pass a Python logger object instead of just using `StderrLogger`.
    *   Consider exposing `Statistics` back to Python.
2.  **C++ Library API:** The separation into Core + Pruners provides a flexible C++ API if used directly. Header structure is reasonable.

**Updated Summary Table:**

| Category        | Improvement Area                                     | Potential Benefit           | Difficulty/Trade-off                     |
| :-------------- | :--------------------------------------------------- | :-------------------------- | :--------------------------------------- |
| **Performance** | Profile & Optimize Data Structures (Heaps/Queues)    | Latency Reduction           | Requires profiling, specific gains vary  |
|                 | Pool/Arena Allocation for Heap Nodes                 | Latency, Cache Locality     | Moderate complexity                      |
|                 | Optimize Pruning Graph Traversals                    | Latency                     | Minor gains likely, complex              |
|                 | Python Output Views (avoid copy)                     | Latency (large results)     | Higher pybind11 complexity             |
| **Memory**      | Alternative Priority Queue (e.g., binary heap)       | Lower Overhead              | Lose decrease-key efficiency           |
|                 | Pool/Arena Allocation                                | Fragmentation, Usage        | Moderate complexity                      |
|                 | Pruning w/o Adjacency List                           | Lower Memory (Pruning)      | Higher algorithm complexity            |
| **Design**      | Internalize `PairingHeap` Buffer                     | Cleaner API                 | Potentially higher memory if no share    |
|                 | Decouple Interfaces (Opaque Handles?)                | Less Header Coupling        | Higher complexity                        |
|                 | **Decouple GWPruner from Core State**                | **Modularity, Testability** | **May change output vs. original test** |
| **Correctness** | **Resolve Test Discrepancies vs Original**           | **CRITICAL**                | Requires deep debugging/analysis         |
|                 | Floating Point Robustness (Scaled Int?)              | Higher Robustness           | Significant code changes                 |
|                 | Ensure Full Exception Coverage                       | Robustness                  | Review code paths                      |
| **Usability**   | Python Logger Integration                            | Better Python Experience    | Moderate pybind11 complexity             |
|                 | Expose Statistics to Python                          | More Insight for User       | Simple                                   |
| **Testing**     | Add Core Algorithm Unit Tests                        | Better Coverage             | Complex state setup needed             |
|                 | **Improve GWPruner Unit Test Isolation**             | Better Test Quality         | Tied to decoupling GWPruner state        |
