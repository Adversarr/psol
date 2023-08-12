#pragma once
#include "common.hpp"

namespace psol::internal {

template <typename _Fp, std::size_t _degree> class Poly {
public:
  using Fp = _Fp;

  static constexpr std::size_t degree = _degree;

  // TODO: add more constructors.
  Poly(std::array<Fp, degree + 1> coe) : c_(coe) {}

  constexpr Fp Get(std::size_t deg) const { return c_.at(deg); }

  template <std::size_t deg> constexpr Fp Get() const noexcept {
    static_assert(deg <= degree && deg >= 0, "Invalid deg.");
    return c_[deg];
  }

  constexpr Fp Eval(Fp x) const noexcept {
    Fp r = c_[0];
    Fp x0 = x;
    for (std::size_t i = 1; i <= degree; ++i) {
      r += x * c_[i];
      x *= x0;
    }

    return r;
  }

  constexpr Fp Grad(Fp x) const noexcept {
    Fp r = c_[1];
    Fp x0 = static_cast<Fp>(1);
    for (std::size_t i = 2; i <= degree; ++i) {
      r += x * c_[i] * i;
      x *= x0;
    }
    return r;
  }

private:
  std::array<Fp, degree + 1> c_;
};

} // namespace psol::internal
