#include "potential/ch_params.hpp"
#include "potential/kd_params.hpp"
#include "potential/params.hpp"
#include "potential/potential.hpp"
#include "potential/wlh_params.hpp"
#include "solver/scatter.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>

// this file is a loose collection of examples of things you can do with this library

// nice shorthand
constexpr auto n = osiris::Proj::neutron;

void print_potential_vals() {

  using namespace osiris;
  auto kdn_uq = KD03Params<n>::build_KDUQ();
  auto wlh_mean = WLH21Params<n>();

  // let's look at potential values for various isotopes on an energy grid
  constexpr auto erg_min = 0.01;
  constexpr auto erg_max = 10.;
  constexpr auto e_range = erg_max - erg_min;
  constexpr auto e_grid_sz = 500;
  auto e_grid = std::array<real, e_grid_sz>{};
  for (int i = 0; i < e_grid_sz; ++i) {
    e_grid[i] =
        erg_min + e_range * static_cast<real>(i) / static_cast<real>(e_grid_sz);
  }

  // mass 144 isotopes
  using Isotope = std::pair<int, int>;
  constexpr auto niso = 6;
  constexpr auto isotopes = std::array<Isotope, niso>{
      Isotope{58, 144}, Isotope{57, 144}, Isotope{56, 144},
      Isotope{55, 144}, Isotope{54, 144}, Isotope{53, 144}};

  std::cout << "Calculating surface potential depths for " << niso
            << " isotopes.\n"
            << std::flush;

  auto out_wlh = std::ofstream("./cmpl_vs_wlh.csv");
  auto out_kd = std::ofstream("./cmpl_vs_kd.csv");

  out_wlh << "E"
          << "\t";
  out_kd << "E"
         << "\t";
  for (const auto& [Z, A] : isotopes) {
    out_wlh << Z << "_" << A << "\t";
    out_kd << Z << "_" << A << "\t";
  }
  out_wlh << "\n";
  out_kd << "\n";

  for (int j = 0; j < e_grid_sz; ++j) {

    const auto erg = e_grid[j];
    out_wlh << std::scientific << std::setprecision(5) << erg << "\t";
    out_kd << std::scientific << std::setprecision(5) << erg << "\t";

    for (int i = 0; i < niso; ++i) {
      const auto& [Z, A] = isotopes[i];

      out_wlh << std::scientific << std::setprecision(5)
              << wlh_mean.cmpl_surf_V(Z, A, erg) << "\t";
      out_kd << std::scientific << std::setprecision(5)
             << kdn_uq.cmpl_surf_V(Z, A, erg) << "\t";
    }
    out_wlh << "\n";
    out_kd << "\n";
  }
}

int main(int, char**) { return 0; };
