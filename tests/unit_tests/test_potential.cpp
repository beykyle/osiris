
#include "potential/params.hpp"
#include "potential/potential.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

TEST_CASE("test KD neutron potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto kd_params = KD03Params<Proj::neutron>();
  auto V = OMP{Xe144, kd_params};
}

TEST_CASE("test CH neutron potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto ch_params = CH89Params<Proj::neutron>();
  auto V = OMP{Xe144, ch_params};
}

TEST_CASE("test WLH neutron potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto wlh_params = WLH21Params<Proj::neutron>();
  auto V = OMP{Xe144, wlh_params};
}

TEST_CASE("test KD proton potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto kd_params = KD03Params<Proj::proton>();
  auto V = OMP{Xe144, kd_params};
}

TEST_CASE("test CH proton potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto ch_params = CH89Params<Proj::proton>();
  auto V = OMP{Xe144, ch_params};
}

TEST_CASE("test WLH proton potential") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr Isotope Xe144{54, 144, 143.93851};

  auto wlh_params = WLH21Params<Proj::proton>();
  auto V = OMP{Xe144, wlh_params};
}
