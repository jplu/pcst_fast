#pragma once

#include <cstdio>
#include <vector>
#include <concepts>

#include "gtest/gtest.h"

namespace test_helpers {
    template <typename T, size_t N> requires std::equality_comparable<T>
    constexpr T* begin(T(&arr)[N]) noexcept { return &arr[0]; }
    template <typename T, size_t N> requires std::equality_comparable<T>
    constexpr T* end(T(&arr)[N]) noexcept { return &arr[0] + N; }

    inline void WriteToStderr(const char* s) {
      std::fprintf(stderr, "%s", s);
      std::fflush(stderr);
    }

    template <typename T> requires std::equality_comparable_with<T, T>
    void CheckResult(const std::vector<T>& expected_result,
                     const std::vector<T>& result) {
        ASSERT_EQ(expected_result.size(), result.size())
            << "Result size mismatch. Expected: " << expected_result.size()
            << ", Got: " << result.size();

        for (size_t ii = 0; ii < expected_result.size(); ++ii) {
            EXPECT_EQ(expected_result[ii], result[ii])
                << "Element mismatch at index " << ii
                << ". Expected: " << expected_result[ii]
                << ", Got: " << result[ii];
        }
    }

}
