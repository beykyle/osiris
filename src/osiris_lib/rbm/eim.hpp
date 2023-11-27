#ifndef EIM_HEADER
#define EIM_HEADER

#include "potential/potential.hpp"
#include "solver/channel.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

namespace osiris {
using namespace xt::placeholders;

/// @brief  Handles the affine-decomposition of a potential where energy and
/// reduced mass are treated as parameters. By default, energy is the first and
/// reduced mass the second elements of the parameter arrays passed in to an
/// instance of this class; the rest are passed into the call to potential
template <
    class array1d_t = xt::xtensor<cmpl, 1>,
    class array2d_t = xt::xtensor<cmpl, 2>,
    class array3d_t = xt::xtensor<cmpl, 3>,
    class array4d_t = xt::xtensor<cmpl, 4>,
    class params_t = xt::xtensor<real, 1>,
    typename std::enable_if_t<xt::is_xexpression<array1d_t>::value, bool> =
        true,
    typename std::enable_if_t<xt::is_xexpression<array2d_t>::value, bool> =
        true,
    typename std::enable_if_t<xt::is_xexpression<array3d_t>::value, bool> =
        true,
    typename std::enable_if_t<xt::is_xexpression<array4d_t>::value, bool> =
        true,
    typename std::enable_if_t<xt::is_xexpression<params_t>::value, bool> = true>
class EnergizedEIMInteractionSpace {
public:
  using potential_t = std::shared_ptr<Potential<params_t>>;

  const int mesh_size;
  const int nbasis;
  const int lmax;

  const std::array<int, 3> eim_matrices_shape = {lmax, nbasis, mesh_size};
  const array1d_t r_matches;
  const array3d_t Ainv_matrices;
  const std::vector<potential_t> potentials;

  EnergizedEIMInteractionSpace(int mesh_size, int nbasis, int lmax,
                               const std::vector<potential_t> &potentials,
                               array3d_t Ainv_matrices, array1d_t r_matches)
      : mesh_size(mesh_size), nbasis(nbasis), lmax(lmax), r_matches(r_matches),
        Ainv_matrices(Ainv_matrices), potentials(potentials) {
    assert(Ainv_matrices.shape() == eim_matrices_shape);
  };

  cmpl tilde(real s, params_t alpha, int l) const {
    return eval_potential(s, alpha, l);
  }

  array1d_t coefficients(params_t alpha, int l) const {
    auto Ainv = xt::strided_view(Ainv_matrices, {l, xt::ellipsis()});
    auto u_real = array1d_t({nbasis});
    for (int i = 0; i < nbasis; ++i)
      u_real(i) = tilde(r_matches(i), alpha, l);
    return xt::linalg::dot(Ainv, u_real);
  }

  array2d_t coefficients(params_t alpha) const {
    auto Ainv = Ainv_matrices;
    auto u_real = array2d_t({lmax, nbasis});
    for (int l = 0; l < lmax; ++l) {
      for (int i = 0; i < nbasis; ++i)
        u_real(l, i) = tilde(r_matches(i), alpha, l);
    }
    return xt::linalg::tensordot(Ainv, u_real, 1);
  }

  real E(params_t alpha) const { return alpha(0); }
  real reduced_mass(params_t alpha) const { return alpha(1); }
  real momentum(params_t alpha) const {
    const auto energy = E(alpha);
    const auto mu = reduced_mass(alpha);
    return sqrt(2 * mu * energy) / constants::hbarc;
  }

private:
  cmpl eval_potential(real s, params_t alpha, int l) const {
    const auto k = momentum(alpha);
    const auto energy = E(alpha);
    const auto r = s / k;
    return potentials(l)->operator()(r, xt::view(alpha, xt::range(2, _))) /
           energy;
  }
};
} // namespace osiris

#endif
