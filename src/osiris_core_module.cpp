#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

#include "potential/wlh_params.hpp"
#include "pybind11/pybind11.h"
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include "potential/potential.hpp"
#include "rbm/bsp.hpp"
#include "solver/scatter.hpp"

namespace py = pybind11;
using namespace osiris;

constexpr auto n = osiris::Proj::neutron;

template <Proj p> auto wlh_n_from_json(std::string fpath) {
  std::ifstream i(fpath);
  json j;
  i >> j;
  return osiris::WLH21Params<p>(j);
}

template <Proj p> auto kduq_n_from_json(std::string fpath) {
  std::ifstream i(fpath);
  json j;
  i >> j;
  return osiris::KD03Params<p>(j);
}

template <typename T>
void declare_bsp_tree(py::module &m, std::string &&typestr) {
  using Class = BinarySPTree<T, xt::pyarray<real>>;
  std::string pyclass_name = std::string("BinarySPTree_") + typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(),
                    py::dynamic_attr())
      .def(
          py::init<int, xt::pyarray<real>, xt::pyarray<real>, xt::pyarray<T>>())
      .def_readonly("depth", &Class::depth)
      .def_readonly("dimensions", &Class::dimensions)
      .def_readonly("bounds_left", &Class::bounds_left)
      .def_readonly("bounds_right", &Class::bounds_right)
      .def("__call__", &Class::operator[])
      .def("at", &Class::at)
      .def("get_bounds", &Class::get_bounds);
}

// Python Module and Docstrings
PYBIND11_MODULE(osiris_core, m) {
  xt::import_numpy();

  m.doc() = R"pbdoc(
        Optical model ScatterIng & ReactIon Software (OSIRIS)

        .. currentmodule:: osiris

        .. autosummary::
           :toctree: _generate

           Isotope
           KD03ParamsNeutron
           KD03ParamsProton
           BinarySPTree
    )pbdoc";

  py::class_<Isotope>(m, "Isotope")
      .def(py::init<int, int, real>(), py::arg("A") = 0, py::arg("Z") = 1,
           py::arg("mass") = constants::p_mass_amu)
      .def_readwrite("A", &Isotope::A)
      .def_readwrite("Z", &Isotope::Z)
      .def_readwrite("mass", &Isotope::mass);

  py::class_<KD03Params<Proj::neutron>>(m, "KD03ParamsNeutron")
      .def(py::init<>())
      .def("from_json", &kduq_n_from_json<Proj::neutron>)
      .def("build_KDUQ", &KD03Params<Proj::neutron>::build_KDUQ)
      .def("real_cent_r", &KD03Params<Proj::neutron>::real_cent_r)
      .def("cmpl_cent_r", &KD03Params<Proj::neutron>::cmpl_cent_r)
      .def("cmpl_surf_r", &KD03Params<Proj::neutron>::cmpl_surf_r)
      .def("real_spin_r", &KD03Params<Proj::neutron>::real_spin_r)
      .def("cmpl_spin_r", &KD03Params<Proj::neutron>::cmpl_spin_r)
      .def("real_cent_a", &KD03Params<Proj::neutron>::real_cent_a)
      .def("cmpl_cent_a", &KD03Params<Proj::neutron>::cmpl_cent_a)
      .def("cmpl_surf_a", &KD03Params<Proj::neutron>::cmpl_surf_a)
      .def("real_spin_a", &KD03Params<Proj::neutron>::real_spin_a)
      .def("cmpl_spin_a", &KD03Params<Proj::neutron>::cmpl_spin_a)
      .def("real_cent_V", &KD03Params<Proj::neutron>::real_cent_V)
      .def("cmpl_cent_V", &KD03Params<Proj::neutron>::cmpl_cent_V)
      .def("cmpl_surf_V", &KD03Params<Proj::neutron>::cmpl_surf_V)
      .def("real_spin_V", &KD03Params<Proj::neutron>::real_spin_V)
      .def("cmpl_spin_V", &KD03Params<Proj::neutron>::cmpl_spin_V)
      .def_readwrite("e_fermi_0", &KD03Params<Proj::neutron>::e_fermi_0)
      .def_readwrite("e_fermi_A", &KD03Params<Proj::neutron>::e_fermi_A)
      .def_readwrite("rv_0", &KD03Params<Proj::neutron>::rv_0)
      .def_readwrite("rv_A", &KD03Params<Proj::neutron>::rv_A)
      .def_readwrite("av_A", &KD03Params<Proj::neutron>::av_A)
      .def_readwrite("rd_0", &KD03Params<Proj::neutron>::rd_0)
      .def_readwrite("rd_A", &KD03Params<Proj::neutron>::rd_A)
      .def_readwrite("ad_0", &KD03Params<Proj::neutron>::ad_0)
      .def_readwrite("ad_A", &KD03Params<Proj::neutron>::ad_A)
      .def_readwrite("rso_0", &KD03Params<Proj::neutron>::rso_0)
      .def_readwrite("rso_A", &KD03Params<Proj::neutron>::rso_A)
      .def_readwrite("aso_0", &KD03Params<Proj::neutron>::aso_0)
      .def_readwrite("v1_0", &KD03Params<Proj::neutron>::v1_0)
      .def_readwrite("v1_asym", &KD03Params<Proj::neutron>::v1_asym)
      .def_readwrite("v1_A", &KD03Params<Proj::neutron>::v1_A)
      .def_readwrite("v2_0", &KD03Params<Proj::neutron>::v2_0)
      .def_readwrite("v2_A", &KD03Params<Proj::neutron>::v2_A)
      .def_readwrite("v3_0", &KD03Params<Proj::neutron>::v3_0)
      .def_readwrite("v3_A", &KD03Params<Proj::neutron>::v3_A)
      .def_readwrite("v4_0", &KD03Params<Proj::neutron>::v4_0)
      .def_readwrite("w1_0", &KD03Params<Proj::neutron>::w1_0)
      .def_readwrite("w1_A", &KD03Params<Proj::neutron>::w1_A)
      .def_readwrite("w2_0", &KD03Params<Proj::neutron>::w2_0)
      .def_readwrite("w2_A", &KD03Params<Proj::neutron>::w2_A)
      .def_readwrite("d1_0", &KD03Params<Proj::neutron>::d1_0)
      .def_readwrite("d1_asym", &KD03Params<Proj::neutron>::d1_asym)
      .def_readwrite("d2_0", &KD03Params<Proj::neutron>::d2_0)
      .def_readwrite("d2_A", &KD03Params<Proj::neutron>::d2_A)
      .def_readwrite("d2_A2", &KD03Params<Proj::neutron>::d2_A2)
      .def_readwrite("d2_A3", &KD03Params<Proj::neutron>::d2_A3)
      .def_readwrite("d3_0", &KD03Params<Proj::neutron>::d3_0)
      .def_readwrite("vso1_0", &KD03Params<Proj::neutron>::vso1_0)
      .def_readwrite("vso1_A", &KD03Params<Proj::neutron>::vso1_A)
      .def_readwrite("vso2_0", &KD03Params<Proj::neutron>::vso2_0)
      .def_readwrite("wso1_0", &KD03Params<Proj::neutron>::wso1_0)
      .def_readwrite("wso2_0", &KD03Params<Proj::neutron>::wso2_0);

  py::class_<KD03Params<Proj::proton>>(m, "KD03ParamsProton")
      .def(py::init<>())
      .def("from_json", &kduq_n_from_json<Proj::proton>)
      .def("build_KDUQ", &KD03Params<Proj::proton>::build_KDUQ)
      .def("real_cent_r", &KD03Params<Proj::proton>::real_cent_r)
      .def("cmpl_cent_r", &KD03Params<Proj::proton>::cmpl_cent_r)
      .def("cmpl_surf_r", &KD03Params<Proj::proton>::cmpl_surf_r)
      .def("real_spin_r", &KD03Params<Proj::proton>::real_spin_r)
      .def("cmpl_spin_r", &KD03Params<Proj::proton>::cmpl_spin_r)
      .def("real_cent_a", &KD03Params<Proj::proton>::real_cent_a)
      .def("cmpl_cent_a", &KD03Params<Proj::proton>::cmpl_cent_a)
      .def("cmpl_surf_a", &KD03Params<Proj::proton>::cmpl_surf_a)
      .def("real_spin_a", &KD03Params<Proj::proton>::real_spin_a)
      .def("cmpl_spin_a", &KD03Params<Proj::proton>::cmpl_spin_a)
      .def("real_cent_V", &KD03Params<Proj::proton>::real_cent_V)
      .def("cmpl_cent_V", &KD03Params<Proj::proton>::cmpl_cent_V)
      .def("cmpl_surf_V", &KD03Params<Proj::proton>::cmpl_surf_V)
      .def("real_spin_V", &KD03Params<Proj::proton>::real_spin_V)
      .def("cmpl_spin_V", &KD03Params<Proj::proton>::cmpl_spin_V)
      .def_readwrite("e_fermi_0", &KD03Params<Proj::proton>::e_fermi_0)
      .def_readwrite("e_fermi_A", &KD03Params<Proj::proton>::e_fermi_A)
      .def_readwrite("rv_0", &KD03Params<Proj::proton>::rv_0)
      .def_readwrite("rv_A", &KD03Params<Proj::proton>::rv_A)
      .def_readwrite("av_A", &KD03Params<Proj::proton>::av_A)
      .def_readwrite("rd_0", &KD03Params<Proj::proton>::rd_0)
      .def_readwrite("rd_A", &KD03Params<Proj::proton>::rd_A)
      .def_readwrite("ad_0", &KD03Params<Proj::proton>::ad_0)
      .def_readwrite("ad_A", &KD03Params<Proj::proton>::ad_A)
      .def_readwrite("rso_0", &KD03Params<Proj::proton>::rso_0)
      .def_readwrite("rso_A", &KD03Params<Proj::proton>::rso_A)
      .def_readwrite("aso_0", &KD03Params<Proj::proton>::aso_0)
      .def_readwrite("v1_0", &KD03Params<Proj::proton>::v1_0)
      .def_readwrite("v1_asym", &KD03Params<Proj::proton>::v1_asym)
      .def_readwrite("v1_A", &KD03Params<Proj::proton>::v1_A)
      .def_readwrite("v2_0", &KD03Params<Proj::proton>::v2_0)
      .def_readwrite("v2_A", &KD03Params<Proj::proton>::v2_A)
      .def_readwrite("v3_0", &KD03Params<Proj::proton>::v3_0)
      .def_readwrite("v3_A", &KD03Params<Proj::proton>::v3_A)
      .def_readwrite("v4_0", &KD03Params<Proj::proton>::v4_0)
      .def_readwrite("w1_0", &KD03Params<Proj::proton>::w1_0)
      .def_readwrite("w1_A", &KD03Params<Proj::proton>::w1_A)
      .def_readwrite("w2_0", &KD03Params<Proj::proton>::w2_0)
      .def_readwrite("w2_A", &KD03Params<Proj::proton>::w2_A)
      .def_readwrite("d1_0", &KD03Params<Proj::proton>::d1_0)
      .def_readwrite("d1_asym", &KD03Params<Proj::proton>::d1_asym)
      .def_readwrite("d2_0", &KD03Params<Proj::proton>::d2_0)
      .def_readwrite("d2_A", &KD03Params<Proj::proton>::d2_A)
      .def_readwrite("d2_A2", &KD03Params<Proj::proton>::d2_A2)
      .def_readwrite("d2_A3", &KD03Params<Proj::proton>::d2_A3)
      .def_readwrite("d3_0", &KD03Params<Proj::proton>::d3_0)
      .def_readwrite("vso1_0", &KD03Params<Proj::proton>::vso1_0)
      .def_readwrite("vso1_A", &KD03Params<Proj::proton>::vso1_A)
      .def_readwrite("vso2_0", &KD03Params<Proj::proton>::vso2_0)
      .def_readwrite("wso1_0", &KD03Params<Proj::proton>::wso1_0)
      .def_readwrite("wso2_0", &KD03Params<Proj::proton>::wso2_0);
  
  py::class_<WLH21Params<Proj::neutron>>(m, "WLH21ParamsNeutron")
      .def(py::init<>())
      .def("from_json", &wlh_n_from_json<Proj::neutron>)
      .def("real_cent_r", &WLH21Params<Proj::neutron>::real_cent_r)
      .def("cmpl_cent_r", &WLH21Params<Proj::neutron>::cmpl_cent_r)
      .def("cmpl_surf_r", &WLH21Params<Proj::neutron>::cmpl_surf_r)
      .def("real_spin_r", &WLH21Params<Proj::neutron>::real_spin_r)
      .def("cmpl_spin_r", &WLH21Params<Proj::neutron>::cmpl_spin_r)
      .def("real_cent_a", &WLH21Params<Proj::neutron>::real_cent_a)
      .def("cmpl_cent_a", &WLH21Params<Proj::neutron>::cmpl_cent_a)
      .def("cmpl_surf_a", &WLH21Params<Proj::neutron>::cmpl_surf_a)
      .def("real_spin_a", &WLH21Params<Proj::neutron>::real_spin_a)
      .def("cmpl_spin_a", &WLH21Params<Proj::neutron>::cmpl_spin_a)
      .def("real_cent_V", &WLH21Params<Proj::neutron>::real_cent_V)
      .def("cmpl_cent_V", &WLH21Params<Proj::neutron>::cmpl_cent_V)
      .def("cmpl_surf_V", &WLH21Params<Proj::neutron>::cmpl_surf_V)
      .def("real_spin_V", &WLH21Params<Proj::neutron>::real_spin_V)
      .def("cmpl_spin_V", &WLH21Params<Proj::neutron>::cmpl_spin_V)
      .def_readwrite("v0", &WLH21Params<Proj::neutron>::v0)
      .def_readwrite("v1", &WLH21Params<Proj::neutron>::v1)
      .def_readwrite("v2", &WLH21Params<Proj::neutron>::v2)
      .def_readwrite("v3", &WLH21Params<Proj::neutron>::v3)
      .def_readwrite("v4", &WLH21Params<Proj::neutron>::v4)
      .def_readwrite("v5", &WLH21Params<Proj::neutron>::v5)
      .def_readwrite("v6", &WLH21Params<Proj::neutron>::v6)
      .def_readwrite("r1", &WLH21Params<Proj::neutron>::r1)
      .def_readwrite("r2", &WLH21Params<Proj::neutron>::r2)
      .def_readwrite("r3", &WLH21Params<Proj::neutron>::r3)
      .def_readwrite("a0", &WLH21Params<Proj::neutron>::a0)
      .def_readwrite("a1", &WLH21Params<Proj::neutron>::a1)
      .def_readwrite("a2", &WLH21Params<Proj::neutron>::a2)
      .def_readwrite("a3", &WLH21Params<Proj::neutron>::a3)
      .def_readwrite("a4", &WLH21Params<Proj::neutron>::a4)
      .def_readwrite("w0", &WLH21Params<Proj::neutron>::w0)
      .def_readwrite("w1", &WLH21Params<Proj::neutron>::w1)
      .def_readwrite("w2", &WLH21Params<Proj::neutron>::w2)
      .def_readwrite("w3", &WLH21Params<Proj::neutron>::w3)
      .def_readwrite("w4", &WLH21Params<Proj::neutron>::w4)
      .def_readwrite("rw0", &WLH21Params<Proj::neutron>::rw0)
      .def_readwrite("rw1", &WLH21Params<Proj::neutron>::rw1)
      .def_readwrite("rw2", &WLH21Params<Proj::neutron>::rw2)
      .def_readwrite("rw3", &WLH21Params<Proj::neutron>::rw3)
      .def_readwrite("rw4", &WLH21Params<Proj::neutron>::rw4)
      .def_readwrite("rw5", &WLH21Params<Proj::neutron>::rw5)
      .def_readwrite("aw0", &WLH21Params<Proj::neutron>::aw0)
      .def_readwrite("aw1", &WLH21Params<Proj::neutron>::aw1)
      .def_readwrite("aw2", &WLH21Params<Proj::neutron>::aw2)
      .def_readwrite("aw3", &WLH21Params<Proj::neutron>::aw3)
      .def_readwrite("aw4", &WLH21Params<Proj::neutron>::aw4)
      .def_readwrite("d0", &WLH21Params<Proj::neutron>::d0)
      .def_readwrite("d1", &WLH21Params<Proj::neutron>::d1)
      .def_readwrite("d2", &WLH21Params<Proj::neutron>::d2)
      .def_readwrite("d3", &WLH21Params<Proj::neutron>::d3)
      .def_readwrite("rs0", &WLH21Params<Proj::neutron>::rs0)
      .def_readwrite("rs1", &WLH21Params<Proj::neutron>::rs1)
      .def_readwrite("rs2", &WLH21Params<Proj::neutron>::rs2)
      .def_readwrite("as0", &WLH21Params<Proj::neutron>::as0)
      .def_readwrite("vso_0", &WLH21Params<Proj::neutron>::vso_0)
      .def_readwrite("vso_1", &WLH21Params<Proj::neutron>::vso_1)
      .def_readwrite("rso_0", &WLH21Params<Proj::neutron>::rso_0)
      .def_readwrite("rso_1", &WLH21Params<Proj::neutron>::rso_1)
      .def_readwrite("aso_0", &WLH21Params<Proj::neutron>::aso_0)
      .def_readwrite("aso_1", &WLH21Params<Proj::neutron>::aso_1);
  
  py::class_<WLH21Params<Proj::proton>>(m, "WLH21ParamsProton")
      .def(py::init<>())
      .def("from_json", &wlh_n_from_json<Proj::proton>)
      .def("real_cent_r", &WLH21Params<Proj::proton>::real_cent_r)
      .def("cmpl_cent_r", &WLH21Params<Proj::proton>::cmpl_cent_r)
      .def("cmpl_surf_r", &WLH21Params<Proj::proton>::cmpl_surf_r)
      .def("real_spin_r", &WLH21Params<Proj::proton>::real_spin_r)
      .def("cmpl_spin_r", &WLH21Params<Proj::proton>::cmpl_spin_r)
      .def("real_cent_a", &WLH21Params<Proj::proton>::real_cent_a)
      .def("cmpl_cent_a", &WLH21Params<Proj::proton>::cmpl_cent_a)
      .def("cmpl_surf_a", &WLH21Params<Proj::proton>::cmpl_surf_a)
      .def("real_spin_a", &WLH21Params<Proj::proton>::real_spin_a)
      .def("cmpl_spin_a", &WLH21Params<Proj::proton>::cmpl_spin_a)
      .def("real_cent_V", &WLH21Params<Proj::proton>::real_cent_V)
      .def("cmpl_cent_V", &WLH21Params<Proj::proton>::cmpl_cent_V)
      .def("cmpl_surf_V", &WLH21Params<Proj::proton>::cmpl_surf_V)
      .def("real_spin_V", &WLH21Params<Proj::proton>::real_spin_V)
      .def("cmpl_spin_V", &WLH21Params<Proj::proton>::cmpl_spin_V)
      .def_readwrite("v0", &WLH21Params<Proj::proton>::v0)
      .def_readwrite("v1", &WLH21Params<Proj::proton>::v1)
      .def_readwrite("v2", &WLH21Params<Proj::proton>::v2)
      .def_readwrite("v3", &WLH21Params<Proj::proton>::v3)
      .def_readwrite("v4", &WLH21Params<Proj::proton>::v4)
      .def_readwrite("v5", &WLH21Params<Proj::proton>::v5)
      .def_readwrite("v6", &WLH21Params<Proj::proton>::v6)
      .def_readwrite("r1", &WLH21Params<Proj::proton>::r1)
      .def_readwrite("r2", &WLH21Params<Proj::proton>::r2)
      .def_readwrite("r3", &WLH21Params<Proj::proton>::r3)
      .def_readwrite("a0", &WLH21Params<Proj::proton>::a0)
      .def_readwrite("a1", &WLH21Params<Proj::proton>::a1)
      .def_readwrite("a2", &WLH21Params<Proj::proton>::a2)
      .def_readwrite("a3", &WLH21Params<Proj::proton>::a3)
      .def_readwrite("a4", &WLH21Params<Proj::proton>::a4)
      .def_readwrite("w0", &WLH21Params<Proj::proton>::w0)
      .def_readwrite("w1", &WLH21Params<Proj::proton>::w1)
      .def_readwrite("w2", &WLH21Params<Proj::proton>::w2)
      .def_readwrite("w3", &WLH21Params<Proj::proton>::w3)
      .def_readwrite("w4", &WLH21Params<Proj::proton>::w4)
      .def_readwrite("rw0", &WLH21Params<Proj::proton>::rw0)
      .def_readwrite("rw1", &WLH21Params<Proj::proton>::rw1)
      .def_readwrite("rw2", &WLH21Params<Proj::proton>::rw2)
      .def_readwrite("rw3", &WLH21Params<Proj::proton>::rw3)
      .def_readwrite("rw4", &WLH21Params<Proj::proton>::rw4)
      .def_readwrite("rw5", &WLH21Params<Proj::proton>::rw5)
      .def_readwrite("aw0", &WLH21Params<Proj::proton>::aw0)
      .def_readwrite("aw1", &WLH21Params<Proj::proton>::aw1)
      .def_readwrite("aw2", &WLH21Params<Proj::proton>::aw2)
      .def_readwrite("aw3", &WLH21Params<Proj::proton>::aw3)
      .def_readwrite("aw4", &WLH21Params<Proj::proton>::aw4)
      .def_readwrite("d0", &WLH21Params<Proj::proton>::d0)
      .def_readwrite("d1", &WLH21Params<Proj::proton>::d1)
      .def_readwrite("d2", &WLH21Params<Proj::proton>::d2)
      .def_readwrite("d3", &WLH21Params<Proj::proton>::d3)
      .def_readwrite("rs0", &WLH21Params<Proj::proton>::rs0)
      .def_readwrite("rs1", &WLH21Params<Proj::proton>::rs1)
      .def_readwrite("rs2", &WLH21Params<Proj::proton>::rs2)
      .def_readwrite("as0", &WLH21Params<Proj::proton>::as0)
      .def_readwrite("vso_0", &WLH21Params<Proj::proton>::vso_0)
      .def_readwrite("vso_1", &WLH21Params<Proj::proton>::vso_1)
      .def_readwrite("rso_0", &WLH21Params<Proj::proton>::rso_0)
      .def_readwrite("rso_1", &WLH21Params<Proj::proton>::rso_1)
      .def_readwrite("aso_0", &WLH21Params<Proj::proton>::aso_0)
      .def_readwrite("aso_1", &WLH21Params<Proj::proton>::aso_1);

  declare_bsp_tree<int>(m, std::string{"int"});

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
