#include <doctest/doctest.h>
#include <iostream>
#include <psol/psol.hpp>
TEST_CASE("degree_1") {
  SUBCASE("1 2") {
    std::array<float, 2> coe{ 1, 2 };
    auto result = psol::solve_boundless<float, 1>(coe, 0.01f);
    CHECK_EQ(result[0], 0.5f);
  }

  SUBCASE("1 0") {
    std::array<float, 2> coe{ 1, 0 };
    auto result = psol::solve_boundless<float, 1>(coe, 0.01f);
    CHECK(isnan(result[0]));
  }
}
