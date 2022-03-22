# diptest

[![Linux Build](https://github.com/RUrlus/diptest/actions/workflows/linux.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/linux.yml)
[![Windows Build](https://github.com/RUrlus/diptest/actions/workflows/windows.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/windows.yml)
[![MacOS build](https://github.com/RUrlus/diptest/actions/workflows/macos.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/macos.yml)

A Python/C implementation of Hartigan & Hartigan's dip test for unimodality.

The dip test measures multimodality in a sample by the maximum difference, over
all sample points, between the empirical distribution function, and the
unimodal distribution function that minimizes that maximum difference. Other
than unimodality, it makes no further assumptions about the form of the null
distribution.

## Dependencies
* `numpy`
* [Optional] `OpenMP`

Parallelisation of the p-value computation using bootstrapping is offered using OpenMP.
OpenMP is disabled by default but can be enabled, see installation section
below.
 Multi-threading can be turned off by setting the number of threads equal to 1. See the docstring of `diptest` for details.


## Installation
```bash
    cd diptest
    pip install .
```

#### Enable OpenMP

To enable OpenMP use:
```bash
    SKBUILD_CONFIGURE_OPTIONS="-DDIPTEST_ENABLE_OPENMP=ON" pip3 install . -v
```

#### Debug installation

To enable a debug build use:
```bash
    SKBUILD_CONFIGURE_OPTIONS="-DCMAKE_BUILD_TYPE=Debug" pip3 install . -v
```

#### Debug printing

To enable the debug print statements use:
```bash
    SKBUILD_CONFIGURE_OPTIONS="-DDIPTEST_ENABLE_DEBUG=ON" pip3 install . -v
```
then call the function with debug argument set to a value greater than zero:
```python3
    diptest(x, debug=1)
```

## Usage

This library provides two functions:
* `dipstat`
* `diptest`

The first only computes Hartigan's dip statistic. `diptest` computes both the
statistic and the p-value. The p-value can be computed using interpolation of a
critical value table (default) or by bootstrapping the null hypothesis.
Note that for larger samples (N > 1e5) this is quite compute and memory intensive.

## References

Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality. The
Annals of Statistics.

Hartigan, P. M. (1985). Computation of the Dip Statistic to Test for
Unimodality. Journal of the Royal Statistical Society. Series C (Applied
Statistics), 34(3), 320-325.

## Acknowledgement

`diptest` is just a Python port of [Martin Maechler's R module of the same
name](http://cran.r-project.org/web/packages/diptest/index.html).
The package wrapping the C implementation was originally written by [Alistair Muldal](https://github.com/alimuldal/diptest).
The fork is an update with a number of changes:
* Fixes a buffer overrun issue in `_dip.c` by reverting to the original C implementation
* Python bindings using Pybind11 (C++) instead of Cython
* P-value computation using bootstrapping has been moved down to C++ with optional parallelisation support through OpenMP
* Removed overhead caused by debug branching statements by placing them under a compile-time definition
* Added tests and wheel support

## License

`diptest` is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
