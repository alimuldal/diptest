diptest
======

A Python/C implementation of Hartigan & Hartigan's dip test for unimodality.

The dip test measures multimodality in a sample by the maximum difference, over
all sample points, between the empirical distribution function, and the
unimodal distribution function that minimizes that maximum difference. Other
than unimodality, it makes no further assumptions about the form of the null
distribution.

### WARNING: This code is no longer actively maintained.

Dependencies
----
* `numpy`
* `Cython`

Installation
----
    $ python setup.py install

References
----

Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality. The
Annals of Statistics.

Hartigan, P. M. (1985). Computation of the Dip Statistic to Test for
Unimodality. Journal of the Royal Statistical Society. Series C (Applied
Statistics), 34(3), 320-325.

Acknowledgement
---

`diptest` is just a Python port of [Martin Maechler's R module of the same
name](http://cran.r-project.org/web/packages/diptest/index.html), and uses a
slightly modified version of his C function for computing the dip statistic.

License
---

`diptest` is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

[![No Maintenance Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)
