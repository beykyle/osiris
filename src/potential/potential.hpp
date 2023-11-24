#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include "potential/params.hpp"
#include "potential/params_base.hpp"
#include "solver/channel.hpp"
#include "util/types.hpp"

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <complex>
#include <memory>
#include <type_traits>
#include <vector>

namespace osiris {

using namespace xt::placeholders;

/// @brief Abstract interface for a local potential.
/// Local implies a potential V is diagonal in coordinate space,
/// e.g. <r|V|rp> = V(r) delta(r-rp)
template <class T,
          typename std::enable_if_t<xt::is_xexpression<T>::value, bool> = true>
struct Potential {
  virtual cmpl operator()(real r, T params) const = 0;
};

/// @brief Abstract interface for a non-local potential that is symmetric in the
/// arguments e.g. V(r1,r2) = V(r2,r1)
template <class T,
          typename std::enable_if_t<xt::is_xexpression<T>::value, bool> = true>
struct NonlocalPotential {
  constexpr static bool is_symmetric = true;
  virtual cmpl operator()(real r, real rp, T params) const = 0;
};

/// @brief Common phenomenological potential form used for central potentials
template <class T> struct WoodsSaxon : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    assert(params.size() == 3);
    auto V = params(0);
    auto R = params(1);
    auto a = params(2);
    if ( fabs(a) < 1e-12  ) return 0;
    return V / (1. + exp((r - R) / a));
  };
};

/// @brief Common phenomenological potential form used for surface peakes
/// potentials. The derivative in r of a Woods-Saxon
template <class T> struct DerivWoodsSaxon : public Potential<T> {

  cmpl operator()(real r, T params) const final {
    assert(params.size() == 3);
    auto V = params(0);
    auto R = params(1);
    auto a = params(2);
    if ( fabs(a) < 1e-12  ) return 0;
    const auto y = exp((r - R) / a);
    return -V / a * (y / ((1. + y) * (1. + y)));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
template <class T> struct Thomas : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    return DerivWoodsSaxon<T>{}(r, params) / r;
  };
};

template <class T> struct Gaussian : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    assert(params.size() == 3);
    auto V = params(0);
    auto R = params(1);
    auto sigma = params(2);
    return exp((r - R) * (r - R) / (sigma * sigma));
  };
};

template <size_t N, class T> struct NGaussian : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    assert(params.size() == N * 3);
    cmpl v = 0;
    for (unsigned int i = 0; i < N; ++i) {
      auto p = xt::view(params, xt::range(0, 3));
      v += Gaussian<decltype(p)>(r, p);
    }
    return v;
  };
};

template <class T> struct DoubleGaussian : NGaussian<2, T> {

  ///@brief Thompson, D. R., LeMere, M., & Tang, Y. C. (1977).
  /// Systematic investigation of scattering problems with the resonating-group
  /// method. Nuclear Physics A, 286(1), 53-66.
  static constexpr T minnesota_params_1s0 = {200, -98.15, 1.487, 0.465};
  static constexpr T minnesota_params_3s1 = {200, -178, 1.487, 0.639};
};

template <class T> struct Yukawa : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    assert(params.size() == 2);
    real mass_coupling = params(0);
    real force_coupling = params(1);
    return -force_coupling * exp(-r / mass_coupling) / r;
  };
};

template <class T> struct SphereWell : public Potential<T> {
  cmpl operator()(real r, T params) const final {
    assert(params.size() == 2);
    real R = params(0);
    cmpl depth = params(1);
    if (r < R)
      return -depth;
    return 0;
  };
};

/// @brief Common optical model potential form
template <class T> struct OMP : public Potential<T> {
  using View = ViewType<T>;
  real l_dot_s{};
  
  OMP(real l_dot_s): l_dot_s(l_dot_s) {};
  
  cmpl operator()(real r, T params) const final {
    assert(params.size() == 18);

    return cmpl{
       WoodsSaxon<View>{}(r, xt::view(params, xt::range(0, 3))) 
    +  DerivWoodsSaxon<View>{}(r, xt::view(params, xt::range(6, 9)))
    +  Thomas<View>{}(r, xt::view(params, xt::range(12, 15))) * l_dot_s
    + WoodsSaxon<View>{}(r, xt::view(params, xt::range(3, 6))) * constants::i
    +  DerivWoodsSaxon<View>{}(r, xt::view(params, xt::range(9, 12))) * constants::i
    +  Thomas<View>{}(r, xt::view(params, xt::range(15, _))) * l_dot_s * constants::i
    };
  }
};

template <class GlobalParamsOMP>
xt::xarray<real> get_global_terms(Isotope iso, real erg_cms,
                                  const GlobalParamsOMP &p) {
  const auto A = iso.A;
  const auto Z = iso.Z;
  return {p.real_cent_V(Z, A, erg_cms), p.real_cent_r(Z, A, erg_cms),
          p.real_cent_a(Z, A, erg_cms), p.cmpl_cent_V(Z, A, erg_cms),
          p.cmpl_cent_r(Z, A, erg_cms), p.cmpl_cent_a(Z, A, erg_cms),
          p.real_surf_V(Z, A, erg_cms), p.real_surf_r(Z, A, erg_cms),
          p.real_surf_a(Z, A, erg_cms), p.cmpl_surf_V(Z, A, erg_cms),
          p.cmpl_surf_r(Z, A, erg_cms), p.cmpl_surf_a(Z, A, erg_cms),
          p.real_spin_V(Z, A, erg_cms), p.real_spin_r(Z, A, erg_cms),
          p.real_spin_a(Z, A, erg_cms), p.cmpl_spin_V(Z, A, erg_cms),
          p.cmpl_spin_r(Z, A, erg_cms), p.cmpl_spin_a(Z, A, erg_cms)};
}

/// @brief An arbitrary local potential in r smeared into the off-diagonal
/// by a Gaussian factor in (r-rp), from:
/// Perey, F., and B. Buck.
/// "A non-local potential model for the scattering of neutrons by nuclei."
/// Nuclear Physics 32 (1962): 353-380.
template <class T>
struct PereyBuck : public NonlocalPotential<T> {
  std::unique_ptr<Potential<ViewType<T>>> local_potential;
  Gaussian<T> non_local_factor;

  PereyBuck(std::unique_ptr<Potential<ViewType<T>>> local_potential)
    : local_potential(std::move(local_potential)) {};

  cmpl operator()(real r, real rp, T params) const final {
    return local_potential->operator()(r, xt::view(params, 0)) *
           non_local_factor((r - rp), xt::view(params, xt::range(1, _)));
  };
};

/// @brief Yamaguchi, Yoshio.
/// "Two-nucleon problem when the potential is nonlocal but separable. I."
/// Physical Review 95.6 (1954): 1628.
/// @param mu reduced mass [amu]
/// @param alpha [fm]^-1
/// @param beta [fm]^-1
template <class T> struct Yamaguchi : public NonlocalPotential<T> {

  ///@brief parameters chosen to reproduce bound state of deuteron and
  /// and neutron-proton triplet scattering length
  static constexpr T dn_triplet_params = {0.2316053, 1.3918324, 41.472};

  cmpl operator()(real r, real rp, T params) const final {
    const auto alpha = params(0);
    const auto beta = params(1);
    const auto f = params(2);
    return f * 2 * beta * (alpha + beta) * (alpha + beta) *
           exp(-beta * (r + rp));
  };

  real analytic_swave_kmatrix(real k, T params) const {
    assert(params.size() == 2);
    const auto a = params(0);
    const auto b = params(1);
    const auto d = 2 * (a + b) * (a + b);
    real cot_delta = (a * b * (a + 2 * b) / d +
                      k * k * (a * a + 2 * a * b + 3 * b * b) / (b * d) +
                      k * k * k * k / (b * d)) /
                     k;
    return 1. / cot_delta;
  }

  static real mu_2_f(real mu) {
    return constants::hbar * constants::hbar * constants::c * constants::c /
           mu * constants::MeV_per_amu;
  }
};

} // namespace osiris

#endif
