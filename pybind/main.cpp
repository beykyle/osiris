#include "pybind11/pybind11.h"
#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <fstream>
#include <iostream>
#include <numeric>
#include <cmath>

#include "solver/scatter.hpp"
#include "potential/potential.hpp"

namespace py = pybind11;
using namespace osiris;

constexpr auto n = osiris::Proj::neutron;

template<Proj p>
auto kduq_n_from_json(std::string fpath) {
  std::ifstream i(fpath);
  json j;
  i >> j;
  return osiris::KD03Params<p>(j); 
}

// Examples

inline double example1(xt::pyarray<double> &m)
{
    return m(0);
}

inline xt::pyarray<double> example2(xt::pyarray<double> &m)
{
    return m + 2;
}

// Readme Examples

inline double readme_example1(xt::pyarray<double> &m)
{
    auto sines = xt::sin(m);
    return std::accumulate(sines.cbegin(), sines.cend(), 0.0);
}

// Vectorize Examples

inline double scalar_func(double i, double j)
{
    return std::sin(i) + std::cos(j);
}

// Python Module and Docstrings

PYBIND11_MODULE(osiris, m)
{
    xt::import_numpy();

    m.doc() = R"pbdoc(
        Optical model ScatterIng & ReactIon Software (OSIRIS)

        .. currentmodule:: osiris

        .. autosummary::
           :toctree: _generate

           example1
           example2
           readme_example1
           vectorize_example1
    )pbdoc";

    m.def("example1", example1, "Return the first element of an array, of dimension at least one");
    m.def("example2", example2, "Return the the specified array plus 2");

    m.def("readme_example1", readme_example1, "Accumulate the sines of all the values of the specified array");

    m.def("vectorize_example1", xt::pyvectorize(scalar_func), "Add the sine and and cosine of the two specified values");
  

    py::class_<Isotope>(m, "Isotope")
    .def(
        py::init<int, int, real>(),
        py::arg("A") = 0,
        py::arg("Z") = 1,
        py::arg("mass") = constants::p_mass_amu
        )
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
}
