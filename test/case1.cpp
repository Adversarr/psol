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

TEST_CASE("degree_2") {
  SUBCASE("1 2 1") {
    std::array<float, 3> coe{ 1, 2, 1 };
    auto result = psol::solve_boundless<float, 2>(coe, 0.01f);
    CHECK_EQ(result[0], -1.0f);
    CHECK(isnan(result[1]));
  }

  SUBCASE("1 1 1") {
    std::array<float, 3> coe{ 1, 1, 1 };
    auto result = psol::solve_boundless<float, 2>(coe, 0.01f);
    CHECK(isnan(result[0]));
  }


  SUBCASE("1 0 -2") {
    std::array<float, 3> coe{ -2, 0, 1 };
    auto result = psol::solve_boundless<float, 2>(coe, 0.0001f);
    CHECK(psol::details::is_small_enough(result[0] - (-1.414f), 0.001f));
    CHECK(psol::details::is_small_enough(result[1] - 1.414f, 0.001f));
  }
}
