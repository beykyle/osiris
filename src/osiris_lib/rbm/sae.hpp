#ifndef SAE_HEADER
#define SAE_HEADER

#include "rbm/eim.hpp"
#include "solver/channel.hpp"

#include <xtensor/xstrided_view.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xslice.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

namespace osiris {


template <class array1d_t = xt::xtensor<cmpl,1>,
          class array2d_t = xt::xtensor<cmpl,2>,
          class array3d_t = xt::xtensor<cmpl,3>,
          class array4d_t = xt::xtensor<cmpl,4>,
          class params_t = xt::xtensor<real,1>,
          typename std::enable_if_t<xt::is_xexpression<array1d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array2d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array3d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array4d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<params_t>::value, bool> = true>
struct Basis {
  using interaction_t = EnergizedEIMInteractionSpace<array1d_t, array2d_t, array3d_t, array4d_t, params_t>;
  
  const int nbasis;
  const std::array<int, 3> basis_vectors_shape = {interaction.lmax, nbasis, interaction.mesh_size,};
  
  const interaction_t interaction;
  
  const array3d_t phi_l_0;
  const array3d_t vectors;
  
  const array4d_t A2_l;
  const array3d_t A_13_l;
  const array2d_t b2_l;
  const array2d_t b13_l;

  Basis(
      int nbasis, 
      array3d_t phi_l_0,
      array3d_t vectors,
      array4d_t A2_l,
      array3d_t A_13_l,
      array2d_t b2_l,
      array2d_t b13_l, 
      interaction_t interaction):
    nbasis(nbasis),
    basis_vectors_shape({interaction.lmax, interaction.mesh_size, nbasis}),
    interaction(interaction),
    phi_l_0(phi_l_0),
    vectors(vectors),
    A2_l(A2_l),
    A_13_l(A_13_l),
    b2_l(b2_l),
    b13_l(b13_l)
    {}

};

template <class array1d_t = xt::xtensor<cmpl,1>,
          class array2d_t = xt::xtensor<cmpl,2>,
          class array3d_t = xt::xtensor<cmpl,3>,
          class array4d_t = xt::xtensor<cmpl,4>,
          class params_t = xt::xtensor<real,1>,
          typename std::enable_if_t<xt::is_xexpression<array1d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array2d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array3d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<array4d_t>::value, bool> = true,
          typename std::enable_if_t<xt::is_xexpression<params_t>::value, bool> = true>
class ReducedBasisEmulator {
public:
  using basis_t = Basis<array1d_t, array2d_t, array3d_t, array4d_t, params_t>;
  using interaction_t  = typename basis_t::interaction_t;
  
  basis_t basis_up;
  basis_t basis_down;
  interaction_t interactions_up;
  interaction_t interactions_down;
  const std::vector<Channel::FermionSpinOrbitCoupling> couplings_up;
  const std::vector<Channel::FermionSpinOrbitCoupling> couplings_down;
  const std::vector<Channel::Asymptotics> asymptotics;

  template<Polarization p>
  array2d_t coefficients(params_t alpha) {

    auto get_coeffs = [](
        auto shape, auto alpha, auto interaction, auto basis)   { 

      auto result = array2d_t(shape);
      const auto beta = interaction.coefficients(alpha);
      const auto A_utilde = xt::linalg::tensordot(beta, basis.A2_l, 1);
      const auto A = A_utilde + basis.A_13_l;
      const auto b_utilde = xt::linalg::tensordot(beta, basis.b2_l, 1);
      const auto b = b_utilde + basis.b13_l;

      for (int l = 0; l < interaction.lmax; ++l) {
        auto x = xt::linalg::solve(
            xt::view(A, l, xt::all(), xt::all()),
            xt::view(b, l, xt::all(), xt::all())
          );
        xt::view(result, l, xt::all()) = x;
      }
      return result;

    };
    
    if constexpr (p == Polarization::down)  {
      auto shape = std::array<int,2> {interactions_up.lmax - 1, interactions_up.nbasis};
      return get_coeffs(shape, alpha, interactions_down, basis_up);
    }
    else if constexpr (p == Polarization::up)  {
      auto shape = std::array<int,2> {interactions_up.lmax, interactions_up.nbasis};
      return get_coeffs(shape, alpha, interactions_up, basis_down);
    }
  }

};

} // namespace osiris

#endif
