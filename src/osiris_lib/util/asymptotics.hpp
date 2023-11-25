
#ifndef ASYMPTOTICS_HEADER
#define ASYMPTOTICS_HEADER

#include "util/constants.hpp"
#include "util/types.hpp"

#include <cmath>
#include <complex>

namespace osiris {

namespace asymptotics {

///@brief reduced spherical Bessel function of 1st kind
struct F {
  int l;
  F(int l) : l(l){};
  real operator()(real z) { return z * std::sph_bessel(l, z); }
};

///@brief reduced spherical Bessel function of 2nd kind
struct G {
  int l;
  G(int l) : l(l){};
  real operator()(real z) { return -z * std::sph_neumann(l, z); }
};

///@brief reduced spherical Hankel function; outgoing
struct h_plus {
  int l;
  h_plus(int l) : l(l){};
  cmpl operator()(real z) {
    using constants::i;
    return G{l}(z) + i * F{l}(z);
  }
};

///@brief reduced spherical Hankel function; incoming
struct h_minus {
  int l;
  h_minus(int l) : l(l){};
  cmpl operator()(real z) {
    using constants::i;
    return G{l}(z)-i * F{l}(z);
  }
};

/// @tparam ReducedphBessel a reduced spherical bessel function or linear combo;
/// one of: { F, G, h_minus, h_plus}
template <class ReducedSphBessel>
/// @brief evaluates d/ds s*f(s)
struct d_dz {
  ReducedSphBessel f;
  d_dz(ReducedSphBessel f) : f(f){};
  cmpl operator()(real z) {
    const auto l = f.l;
    return f(z) * static_cast<cmpl>(1 + l) / z - ReducedSphBessel{l + 1}(z);
  }
};
} // namespace asymptotics
} // namespace osiris
#endif
