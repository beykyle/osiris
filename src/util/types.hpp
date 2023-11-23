#ifndef TYPES_HEADER
#define TYPES_HEADER

#include <complex>

namespace osiris {
using real = double;
using cmpl = std::complex<real>;

template <typename T> bool are_equal(T a, T b, T epsilon = 1.0e-6) {
  return fabs(a - b) < epsilon;
}

} // namespace osiris

#endif
