#pragma once

#include <string>

namespace cluster_approx {

using NodeId = int;
using EdgeId = int;
using EdgePartId = int;
using ClusterId = int;
using EventId = int;

constexpr NodeId kInvalidNodeId = -1;
constexpr EdgeId kInvalidEdgeId = -1;
constexpr EdgePartId kInvalidEdgePartId = -1;
constexpr ClusterId kInvalidClusterId = -1;
constexpr EventId kInvalidEventId = -1;

enum class PruningMethod {
    kNone = 0,
    kSimple,
    kGW,
    kStrong,
    kUnknown
};

[[nodiscard]] PruningMethod parse_pruning_method(const std::string& input);

}