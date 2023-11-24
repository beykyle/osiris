#include "rbm/bsp.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <xtensor/xarray.hpp>

using namespace osiris;

TEST_CASE("BSP with depth < dimensions") {
  // divides 3D volume into 4 quadrants along x=0 and y=0
  // as long as z in [-1,1], it doesn't affect which partition we're in
  auto bsp = BinarySPTree<int>(2, xt::xarray<real>{-1., -1, -1},
                               xt::xarray<real>{1, 1, 1}, {0, 1, 2, 3});

  // x<0, y<0
  REQUIRE(bsp[xt::xarray<real>{-1. / 2, -1. / 2, 1. / 2}] == 0);
  REQUIRE(bsp[xt::xarray<real>{-1. / 2, -1. / 2, -1. / 2}] == 0);

  // x<0, y>0
  REQUIRE(bsp[xt::xarray<real>{-1. / 2, 1. / 2, 1. / 2}] == 1);
  REQUIRE(bsp[xt::xarray<real>{-1. / 2, 1. / 2, -1. / 2}] == 1);

  // x>0, y<0
  REQUIRE(bsp[xt::xarray<real>{1. / 2, -1. / 2, 1. / 2}] == 2);
  REQUIRE(bsp[xt::xarray<real>{1. / 2, -1. / 2, -1. / 2}] == 2);

  // x>0, y>0
  REQUIRE(bsp[xt::xarray<real>{1. / 2, 1. / 2, 1. / 2}] == 3);
  REQUIRE(bsp[xt::xarray<real>{1. / 2, 1. / 2, -1. / 2}] == 3);
}

TEST_CASE("BSP with depth > dimensions") {
  // divides 2D box into 8 sub-boxes using planes y=1/2 and x=1/4, x=1/2, x=3/4
  //
  //     y
  //     |
  //  1  -------------
  //     |2 |3 |6 |7 |
  // 1/2 -------------
  //     |0 |1 |4 |5 |
  //   0 ---------------x
  //     0    1/2    1
  //
  auto bsp =
      BinarySPTree<int>(3, xt::xarray<real>{0, 0}, xt::xarray<real>{1, 1},
                        {0, 1, 2, 3, 4, 5, 6, 7});

  REQUIRE(bsp[xt::xarray<real>{0.1, 0.1}] == 0);
  REQUIRE(bsp[xt::xarray<real>{0., 0.}] == 0);
  REQUIRE(bsp[xt::xarray<real>{0.25, 0.1}] == 1);
  REQUIRE(bsp[xt::xarray<real>{0.0, 0.999}] == 2);
  REQUIRE(bsp[xt::xarray<real>{0.250, 0.5}] == 3);
  REQUIRE(bsp[xt::xarray<real>{0.4999, 0.99}] == 3);
  REQUIRE(bsp[xt::xarray<real>{0.4999, 0.0}] == 1);
  REQUIRE(bsp[xt::xarray<real>{0.5, 0.0}] == 4);
  REQUIRE(bsp[xt::xarray<real>{0.5, 0.9}] == 6);
  REQUIRE(bsp[xt::xarray<real>{0.9, 0.4999}] == 5);
  REQUIRE(bsp[xt::xarray<real>{0.999, 0.99}] == 7);
}

TEST_CASE("BSP with custom point") {
  // divides 2D box into 8 sub-boxes using planes y=1/2 and x=1/4, x=1/2, x=3/4
  //
  //     y
  //     |
  //  1  -------------
  //     |2 |3 |6 |7 |
  // 1/2 -------------
  //     |0 |1 |4 |5 |
  //   0 ---------------x
  //     0    1/2    1
  //
  struct Point {
    double x{}, y{};
    double &operator[](int i) {
      assert(i >= 0 and i <= 1);
      return i == 0 ? x : y;
    }
    double operator()(int i) const {
      assert(i >= 0 and i <= 1);
      return i == 0 ? x : y;
    }
    auto size() const { return 2; }
  };

  auto bsp =
      BinarySPTree<int, Point>(3, Point{.x = 0, .y = 0}, Point{.x = 1, .y = 1},
                               {0, 1, 2, 3, 4, 5, 6, 7});

  REQUIRE(bsp[{0.1, 0.1}] == 0);
  REQUIRE(bsp[{0., 0.}] == 0);
  REQUIRE(bsp[{0.25, 0.1}] == 1);
  REQUIRE(bsp[{0.0, 0.999}] == 2);
  REQUIRE(bsp[{0.250, 0.5}] == 3);
  REQUIRE(bsp[{0.4999, 0.99}] == 3);
  REQUIRE(bsp[{0.4999, 0.0}] == 1);
  REQUIRE(bsp[{0.5, 0.0}] == 4);
  REQUIRE(bsp[{0.5, 0.9}] == 6);
  REQUIRE(bsp[{0.9, 0.4999}] == 5);
  REQUIRE(bsp[{0.999, 0.99}] == 7);
}
