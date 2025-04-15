#pragma once

#include <cstdint>
#include <string>

namespace cluster_approx {

using NodeId = int;
using EdgeId = int;
using EdgePartId = int;
using ClusterId = int;
using EventId = int;

const NodeId kInvalidNodeId = -1;
const EdgeId kInvalidEdgeId = -1;
const EdgePartId kInvalidEdgePartId = -1;
const ClusterId kInvalidClusterId = -1;
const EventId kInvalidEventId = -1;


/**
 * @brief Enumerates the available pruning methods for the PCST algorithm.
 */
enum class PruningMethod {
    kNone = 0,        // No pruning, return the raw GW forest.
    kSimple,          // Simple pruning: remove nodes not connected to the main component(s).
    kGW,              // Goemans-Williamson style pruning based on merge events.
    kStrong,          // Strong pruning based on subtree contribution.
    kUnknown          // Represents an invalid or unparsed pruning method.
};

/**
 * @brief Parses a string representation of a pruning method. Case-insensitive.
 * @param input The string to parse (e.g., "none", "simple", "gw", "strong").
 * @return The corresponding PruningMethod enum value, or PruningMethod::kUnknown if not recognized.
 */
[[nodiscard]] PruningMethod parse_pruning_method(const std::string& input);

}