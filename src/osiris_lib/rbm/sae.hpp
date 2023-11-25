#ifndef SAE_HEADER
#define SAE_HEADER

#include "potential/potential.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"
#include "xtensor/xstrided_view.hpp"

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xslice.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace osiris {
using namespace xt::placeholders;

class EnergizedEIMInteractionSpace {
public:
  using params_t = xt::xtensor<real, 1>;

  const int mesh_size;
  const int nbasis;
  const int lmax;

  const xt::xtensor<real, 1> r_matches;
  const std::array<int, 3> eim_matrices_shape = {lmax, mesh_size, nbasis};

  xt::xtensor<cmpl, 3> Ainv_matrices;

  cmpl tilde(real s, params_t alpha, int l) const {
    return eval_potential(s, alpha, l);
  };

  xt::xtensor<cmpl, 1> coefficients(params_t alpha, int l) const {
    auto Ainv = xt::strided_view(Ainv_matrices, {l, xt::ellipsis()});
    auto u_real = xt::xtensor<cmpl, 1>({nbasis});
    for (int i = 0; i < nbasis; ++i)
      u_real(i) = tilde(r_matches(i), alpha, l);
    return xt::linalg::dot(Ainv, u_real);
  }

  EnergizedEIMInteractionSpace(
      int mesh_size, int nbasis, int lmax,
      xt::xtensor<std::unique_ptr<Potential<xt::xtensor<real, 1>>>, 1>
          &&potential,
      xt::xtensor<cmpl, 3> Ainv_matrices, xt::xtensor<real, 1> r_matches)
      : mesh_size(mesh_size), nbasis(nbasis), lmax(lmax),
        potential(std::move(potential)), Ainv_matrices(Ainv_matrices),
        r_matches(r_matches){};

  real E(params_t alpha) const { return alpha(0); }

  real reduced_mass(params_t alpha) const { return alpha(1); }

  real momentum(params_t alpha) const {
    const auto energy = E(alpha);
    const auto mu = reduced_mass(alpha);
    return sqrt(2 * mu * energy) / constants::hbarc;
  }

private:
  const xt::xtensor<std::unique_ptr<Potential<params_t>>, 1> potential;

  cmpl eval_potential(real s, params_t alpha, int l) const {
    const auto k = momentum(alpha);
    const auto energy = E(alpha);
    const auto r = s / k;
    return potential(l)->operator()(r, xt::view(alpha, xt::range(2, _))) /
           energy;
  }
};

} // namespace osiris

#endif
