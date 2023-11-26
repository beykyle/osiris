[![CMake](https://github.com/beykyle/osiris/actions/workflows/cmake.yml/badge.svg)](https://github.com/beykyle/osiris/actions/workflows/cmake.yml) [![Python package](https://github.com/beykyle/osiris/actions/workflows/python-package.yml/badge.svg)](https://github.com/beykyle/osiris/actions/workflows/python-package.yml)

OSIRIS   (/oʊˈsaɪərɨs/)
==============
Optical model Scatterong & Reaction Software 

A fast, modern nuclear reaction code with a `python` interface using [xtensor-python](https://github.com/xtensor-stack/xtensor-python), built for uncertainty quantificaiton and model order reduction of parametric nucleon-nuclear interactions.

Python Installation
------------

`osiris` is available at [pypi.org/project/osiris-py](https://pypi.org/project/osiris-py):

 - `pip install osiris-pu`

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
FetchContent_MakeAvailable(osiris_lib)
```

Now you can `#include` files like `"potential/params.hpp"` into your project, as long as you make `osiris_lib` a dependency of the relevant target in your `CMakeLists.txt`.


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

### c++ documentation

Following a build, run:

```
make docs
```

from the `build/` directory. 

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
- a modern `c++` compiler and tool chain. This was tested and built with:

```
clang version 15.0.7 (https://github.com/conda-forge/clangdev-feedstock fc523913ae327dfa0a91bb2b45a36c810e0f55d0)
Target: x86_64-unknown-linux-gnu
Thread model: posix
```

### handled by `CMake` (you don't have to do anything):
- [`nlohmann/json`](https://github.com/nlohmann/json)
- [`Catch2`](https://catch2.docsforge.com/)
- [`xtensor`](https://github.com/xtensor-stack/xtensor)
- [`xtl`](https://github.com/xtensor-stack/xtl)

### install yourself if you want to `import osiris` in `python`:
- [`python`](https://www.python.org/) 3.7+
- [`numpy`](https://numpy.org/)
- [`pybind11`](https://pybind11.readthedocs.io/en/stable/index.html)
- [`xtensor-python`](https://github.com/xtensor-stack/xtensor-python)
- [`pytest`](https://docs.pytest.org/en/7.4.x/) to run the `python` unit tests

It is highly recomended to use use a package, dependency and environment manager like [mamba](https://mamba.readthedocs.io/en/latest/) or [conda](https://docs.conda.io/en/latest/). Then, setting up an environment to run `osiris` with `python` is as easy as (e.g. using `mamba`), from `pypi`:

```zsh
mamba create -n osirenv python cmake compilers numpy pybind11 xtensor-python pytest
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


Citation
-----------------

