osiris
==============

Optical model ScatterIng & ReactIon Software (OSIRIS)

A nuclear reaction code with a `python` interface.

This is a library of potentials relevant for nuclear reactions, with a primary focus on nucleon-nuclear scattering using optical potentials. It is equipped with an R-Matrix solver on a Lagrange-Legendre mesh for calculating scattering matrix elements and cross sections. Applications include fitting and uncertainty quantification of optical model parameters, and Hauser-Feshbach inelastic scattering calculations. Also included is a module for running the online stage of a trained reduced basis emulator using tangent subspace projection with a binary partition of parameter space. 

The goals of this project are to provide a fast, flexible and high-fidelity solver for nuclear reactions, with a focus on model order reduction and uncertainty quantification. It is meant to be readable, user-friendly, and developed with modern software tools. It can be included as a ``c++`` library using `CMake`, or it can be built as a python module with `pip`.


Python Installation
------------

`osiris` is available on [pypi](pypi.org) at []():

**On Unix (Linux, OS X)**

 - `pip install osiris`

**On Windows (Requires Visual Studio 2015)**

 - For Python 3.7+:
     - `pip install osiris`

Windows runtime requirements
----------------------------

On Windows, the Visual C++ 2015 redistributable packages are a runtime
requirement for this project. It can be found [here](https://www.microsoft.com/en-us/download/details.aspx?id=48145).

If you use the Anaconda python distribution, you may require the Visual Studio
runtime as a platform-dependent runtime requirement for you package:

```yaml
requirements:
  build:
    - python
    - setuptools
    - pybind11

  run:
   - python
   - vs2015_runtime  # [win]
```

Standalone executable 
----------------------------

`osiris` can also be used as a standalone application, with an example `main` function living in`exec/osiris_example_app.cpp`. To build:

```zsh
git clone git@github.com:beykyle/osiris.git
cd omplib
mkdir build
cd build 
cmake -DCMAKE_BUILD_TYPE=Release .. 
make 
make test
make docs
```


Building the documentation
--------------------------

### python documentation

Documentation for the example project is generated using Sphinx. Sphinx has the
ability to automatically inspect the signatures and documentation strings in
the extension module to generate beautiful documentation in a variety formats.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `osiris/docs`
 - `make html`

### c++ documentation

Following a build, run:

```
make docs
```

from the `build/` directory. 

CMake integration
-----------------

`osiris` supports CMake integration using [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) with static linkage. Add the following to your `CMakeLists.txt`:

```cmake
FetchContent_Declare(
  omplib GIT_REPOSITORY https://github.com/beykyle/osiris.git
  GIT_TAG "origin/main"
  )
FetchContent_MakeAvailable(osiris_lib)
```

Now you can `#include` files like `"potential/params.hpp"` into your project, as long as you make `osiris_lib` a dependency of the relevant target.


Running the tests
-----------------

### python unit tests

Running the tests requires `pytest`.

```bash
py.test .
```

### c++ unit tests

Following a build, run:

```
make test
```

from the `build/` directory. 

Dependencies
-----------------

- [CMake](https://cmake.org/) >= 3.18
- a modern compiler. This was tested and built with:

```
clang version 15.0.7 (https://github.com/conda-forge/clangdev-feedstock fc523913ae327dfa0a91bb2b45a36c810e0f55d0)
Target: x86_64-unknown-linux-gnu
Thread model: posix
```

### handled by `CMake` (you don't have to do anything):
- [nlohmann/json](https://github.com/nlohmann/json)
- [Catch2](https://catch2.docsforge.com/)
- [xtensor](https://github.com/xtensor-stack/xtensor)
- [xtl](https://github.com/xtensor-stack/xtl)

### install yourself if you want to use the `python` module, `osiris`:
- [python](https://www.python.org/) 3.7+
- [numpy](https://numpy.org/)
- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)
- [xtensor-python](https://github.com/xtensor-stack/xtensor-python)

It is highly recomended to use use a package, dependency and environment manager like [mamba](https://mamba.readthedocs.io/en/latest/) or [conda](https://docs.conda.io/en/latest/). Then, setting up an environment to run `osiris` with `python` is as easy as (e.g. using `mamba`):

```zsh
mamba create -n osirenv python cmake compilers pybind11 numpy xtensor-python
mamba activate osirenv
git clone git@github.com:beykyle/osiris.git
cd osiris
py setup.py build -j{nproc}
pip install .
```

Citation
-----------------

