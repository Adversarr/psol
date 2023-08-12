
#include "internal/algo.hpp"
#include "internal/poly.hpp"

namespace psol {
template <typename Fp, std::size_t degree>
Solution<Fp, degree> solve(std::array<Fp, degree+1> coefficients) {
  using P = internal::Poly<Fp, degree>;
  P poly(coefficients);

  return internal::Algo<P, Fp, degree>{}.Boundless(poly, 0.001);
}

} // namespace psol
