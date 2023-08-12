#include "common.hpp"
#include <cmath>
#include <limits>

namespace psol {

namespace internal {

template <typename Poly, typename Fp, std::size_t degree> struct Algo;

template <typename Fp> bool is_small_enough(Fp value, Fp epsilon) {
  return value > -epsilon && value < epsilon;
}

template <typename Fp> Fp nan(Fp placeholder = 0) {
  return std::numeric_limits<Fp>::quiet_NaN();
}

template <typename Fp, std::size_t low_deg, std::size_t ext>
std::array<Fp, low_deg + ext> inline compose_nan(
    const std::array<Fp, low_deg> &lower) {
  std::array<Fp, low_deg + ext> result;
  std::copy(lower.begin(), lower.end(), result.begin());
  std::fill_n(std::advance(result.begin(), low_deg), nan<Fp>());
  return result;
}

template <typename Poly, typename Fp> struct Algo<Poly, Fp, 1> {
  std::array<Fp, 1> Boundless(const Poly &c, Fp eps) const noexcept {
    Fp bias = c.template Get<0>();
    Fp k = c.template Get<1>();
    Fp base_sol = bias / k;

    // TODO: Deal with |k| < eps
    return std::array<Fp, 1>{base_sol};
  }
};

template <typename Poly, typename Fp = typename Poly::Fp>
Fp do_newton(const Poly &poly, Fp x0, Fp eps, std::size_t max_iter) {
  // Newton's formula:
  //   $$ X_n+1 = X_n - F / F' $$
  Fp f = poly.Eval(x0);
  if (is_small_enough(f, eps)) {
    return x0;
  }

  if (max_iter == 0) {
    // Not converged.
    return nan<Fp>();
  }

  Fp g = poly.Grad(x0);

  // TODO: `g` is small?
  return do_newton(poly, x0 - f / g, eps, max_iter - 1);
}

template <typename Poly, typename Fp> struct Algo<Poly, Fp, 2> {
  std::array<Fp, 2> Boundless(const Poly &coe, Fp eps) const noexcept {
    // c + b x + a x2 = 0
    Fp c = coe.template Get<0>();
    Fp b = coe.template Get<1>();
    Fp a = coe.template Get<2>();

    // TODO: if a is small enough, use solver with degree=1
    if (is_small_enough(a, eps)) {
      return compose_nan<Fp, 1, 1>(Algo<Poly, Fp, 1>::Boundless(coe, eps));
    }

    // test whether the poly does have any root in R.
    Fp delta = b * b - 4 * a * c;
    if (delta < -eps) {
      // no real root.
      return std::array<Fp, 2>{nan<Fp>(), nan<Fp>()};
    }

    Fp mid_point = -b / a * 0.5;
    if (delta < eps) {
      // one root.
      return std::array<Fp, 2>{mid_point, nan<Fp>()};
    }
    // otherwise, compute using Newton formula.

    // guess the solution.
    Fp dt = 0.25 * delta;
    Fp left = mid_point - dt;
    Fp right = mid_point + dt;

    // compute the solution
    return std::array<Fp, 2>{do_newton(coe, left, eps, 10),
                             do_newton(coe, right, eps, 10)};
  }
};
} // namespace internal

} // namespace psol
