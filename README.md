[![CMake](https://github.com/beykyle/osiris/actions/workflows/cmake.yml/badge.svg)](https://github.com/beykyle/osiris/actions/workflows/cmake.yml) [![Python package](https://github.com/beykyle/osiris/actions/workflows/python-package.yml/badge.svg)](https://github.com/beykyle/osiris/actions/workflows/python-package.yml)

OSIRIS   (/oʊˈsaɪərɨs/)
==============
Optical model Scattering & Reaction Software 

A fast, modern nuclear reaction code with a `python` interface using [xtensor-python](https://github.com/xtensor-stack/xtensor-python), built for uncertainty quantificaiton and model order reduction of parametric nucleon-nuclear interactions. Contains two main components; a `c++` library `osiris_lib` containing the core solvers and physics functionality, and a `python` front end; `osiris`.

Python Installation
------------

`osiris` is available at [pypi.org/project/osiris-py](https://pypi.org/project/osiris-py):

 - `pip install osiris-py`

Example use
--------------

```python
import osiris

```

For more in-depth examples and tutorials, see [`examples/`](https://github.com/beykyle/osiris/tree/main/examples)

CMake integration
-----------------

`osiris` supports CMake integration using [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html). Add the following to your `CMakeLists.txt`:

```cmake
FetchContent_Declare(
  osiris 
  GIT_REPOSITORY https://github.com/beykyle/osiris.git
  GIT_TAG "origin/main"
  )
FetchContent_MakeAvailable(osiris)
```

Now you can `#include` files like `"potential/params.hpp"` into your project, as long as you link `osiris_lib` to the relevant target in your `CMakeLists.txt`, e.g.

```cmake
target_link_libraries(<MY_TARGET> osiris_lib)
```

Standalone executable 
----------------------------

`osiris` can also be used as a standalone application, with an example `main` function living in`exec/osiris_example_app.cpp`. To build:

```zsh
git clone git@github.com:beykyle/osiris.git
cd osiris
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
From the `docs/` directory:

```
make html
```

### `c++` documentation

Following a build, run:

```
make docs
```

from the `build/` directory. 

Running the tests
-----------------

### `osiris` (`python`) unit tests

Running the tests requires `pytest`.

```bash
py.test .
```

### `osiris_lib` (`c++`) unit tests

Following a build, run:

```
make test
```

from the `build/` directory. 

`osiris_lib` dependencies
-----------------

- [CMake](https://cmake.org/) >= 3.18
- a modern `c++` compiler and tool chain. This was tested and built primarily with the following compilers, running on Ubuntu, using a dependency stack from `conda-forge`:

```
clang version 15.0.7 (https://github.com/conda-forge/clangdev-feedstock fc523913ae327dfa0a91bb2b45a36c810e0f55d0)
Target: x86_64-unknown-linux-gnu
Thread model: posix
```

and 

```
g++ (conda-forge gcc 10.4.0-19) 10.4.0
```

### handled by `CMake` (you don't have to do anything):
- [`nlohmann/json`](https://github.com/nlohmann/json)
- [`Catch2`](https://catch2.docsforge.com/)
- [`xtensor`](https://github.com/xtensor-stack/xtensor)
- [`xtl`](https://github.com/xtensor-stack/xtl)
- [`xtensor-blas`](https://github.com/xtensor-stack/xtensor-blas)

`osiris` dependencies for `python` bindings
-----------------

If you only want to use `osiris_lib` as a library for your project, these are the only dependencies. To use the `python` bindings, you will also need:

### handled by `CMake` (you don't have to do anything):
- [`xtensor-python`](https://github.com/xtensor-stack/xtensor-python)
- [`pybind11`](https://pybind11.readthedocs.io/en/stable/index.html)

### install yourself
- [`python`](https://www.python.org/) 3.7+
- [`numpy`](https://numpy.org/)
- [`scikit-build-core`](https://github.com/scikit-build/scikit-build-core)
- [`setuptools-scm`](https://pypi.org/project/setuptools-scm/)

### install yourself, optional:
- [`pytest`](https://docs.pytest.org/en/7.4.x/) to run the `python` unit tests

It is recomended to use use a package, dependency and environment manager like [mamba](https://mamba.readthedocs.io/en/latest/) or [conda](https://docs.conda.io/en/latest/). Then, setting up an environment to run `osiris` with `python` is as easy as (e.g. using `mamba`), from `pypi`:

```zsh
mamba create -n osirenv python cmake compilers numpy pytest
mamba activate osirenv
pip install osiris
```
Alternatively, to install for development purposes, clone the repository and create an editable install: 

```zsh
mamba create -n osirenv python cmake compilers pybind11 numpy xtensor-python
mamba activate osirenv
git clone git@github.com:beykyle/osiris.git
cd osiris
py setup.py build -j{nproc}
pip install -e .
```

Parallel execution with MPI
----------------------------
```python
import osiris
import mpi4py
```

Using [`ipyparallel`](https://ipyparallel.readthedocs.io/en/latest/), notebooks can be run with a parallel backend.

Windows runtime requirements
----------------------------

On Windows, the Visual `c++` 2015 redistributable packages are a runtime
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


Citation
-----------------
```latex
@software{Beyer_OSIRIS,
author = {Beyer, Kyle},
license = {BSD-3-Clause},
title = {{OSIRIS}},
url = {https://github.com/beykyle/osiris},
version = {0.1}
}
```
