#include <doctest/doctest.h>
#include <iostream>
#include <psol.hpp>
TEST_CASE("degree_1") {
  std::array<float, 2> coe{1, 2};
  auto result = psol::solve<float, 1>(coe);
  CHECK_EQ(result[0], 0.5f);
}
