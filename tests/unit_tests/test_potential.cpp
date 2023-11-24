
#include "potential/params.hpp"
#include "potential/potential.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

using namespace osiris;
using namespace osiris::constants;

constexpr Isotope Xe144{144, 54, 143.93851};
const auto erg_cms = 14.001860504200815;

TEST_CASE("test KD neutron potential") {

  auto kd_params = KD03Params<Proj::neutron>();
  auto omp_params = get_global_terms(Xe144, erg_cms, kd_params);
  auto V = OMP<xt::xarray<real>>(1./2.);
  
  REQUIRE( V(0.8, omp_params).real() == Approx(-43.363118795927214) );
  REQUIRE( V(0.8, omp_params).imag() == Approx(-0.8760279215437633) );
}

TEST_CASE("test CH neutron potential") {

  auto ch_params = CH89Params<Proj::neutron>();
  auto omp_params = get_global_terms(Xe144, erg_cms,ch_params);
  auto V = OMP<xt::xarray<real>>(1./2.);
}

TEST_CASE("test WLH neutron potential") {

  auto wlh_params = WLH21Params<Proj::neutron>();
  auto omp_params = get_global_terms(Xe144, erg_cms, wlh_params);
  auto V = OMP<xt::xarray<real>>(1./2.);
}
