#ifndef TYPES_HEADER
#define TYPES_HEADER

#include <complex>
#include <vector>
#include <xtensor/xview.hpp>

namespace osiris {
using real = double;
using cmpl = std::complex<real>;

template <typename T> bool are_equal(T a, T b, T epsilon = 1.0e-6) {
  return fabs(a - b) < epsilon;
}

template <class T> struct is_std_vector {
  static const bool value = false;
};

template <class T> struct is_std_vector<std::vector<T>> {
  static const bool value = true;
};

namespace detail {
template <class T> constexpr static auto get_range() {
  auto x = T{};
  return xt::view(x, xt::range(0, 0));
}
} // namespace detail

template <class T> using ViewType = decltype(detail::get_range<T>());

} // namespace osiris

#endif
