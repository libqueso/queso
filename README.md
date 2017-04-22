The `QUESO` Library [![Build Status](https://travis-ci.org/libqueso/queso.svg?branch=dev)](https://travis-ci.org/libqueso/queso) [![Coverage](https://codecov.io/gh/libqueso/queso/coverage.svg?branch=dev)](https://codecov.io/gh/libqueso/queso)
===================

`QUESO` stands for Quantification of Uncertainty for Estimation,
Simulation and Optimization.

`QUESO` is a collection of algorithms and other functionalities aimed
for the solution of statistical inverse problems, for the solution of
statistical forward problems, for the validation of a model and for
the prediction of quantities of interest from such model along with
the quantification of their uncertainties.

`QUESO` is designed for flexibility, portability, easiness of use and
easiness of extension. Its software design follows an object-oriented
approach and its code is written on C++ and over MPI. It can run over
uniprocessor or multiprocessor environments.

Installation
------------

You can obtain `QUESO` tarballs
[here](https://github.com/libqueso/queso/releases).

If you do not have a `configure` script in the top level directory,
run `bootstrap` to generate a configure script using `autotools`.

Before compiling, you must run the `configure` script.  To run, type
`./configure`.  Additional options may be provided if desired.  Run
`./configure --help` for details.

After successfully running `configure`, type `make` to build the
`QUESO` library.

Then type `make install` to install it in the directory previously
specified by the `--prefix` option of the `configure` script.

Documentation
-------------

`QUESO` documentation is available
[here](http://libqueso.github.io/queso/html/).

Documentation for older versions:

-  [v0.56.1](http://libqueso.github.io/queso/v0.56.1/html/)
-  [v0.56.0](http://libqueso.github.io/queso/v0.56.0/html/)
-  [v0.55.0](http://libqueso.github.io/queso/v0.55.0/html/)
-  [v0.54.0](http://libqueso.github.io/queso/v0.54.0/html/)
-  [v0.53.0](http://libqueso.github.io/queso/v0.53.0/html/)
-  [v0.52.0](http://libqueso.github.io/queso/v0.52.0/html/)
-  [v0.51.1](http://libqueso.github.io/queso/v0.51.1/html/)
-  [v0.50.1](http://libqueso.github.io/queso/v0.50.1/html/)

Dependencies
------------

At a minimum, `QUESO` compilation requires MPI and linkage against two
external libraries:

- [Boost](http://www.boost.org/) v1.35+
- [GSL](https://www.gnu.org/software/gsl/) v1.10+

`QUESO` also has several optional dependencies that enable additional functionality:

- [Teuchos](http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html) v11.0.2+
- [GRVY](https://red.ices.utexas.edu/projects/hpct/files) v0.29+
- [HDF5](http://www.hdfgroup.org/HDF5/) v1.8.0+
- [GLPK](https://www.gnu.org/software/glpk/) v4.35+
- [PETSc](http://www.mcs.anl.gov/petsc/)

Should you be interested in using the optional infinite dimensional
capabilities of `QUESO`, then you also need the following dependencies:

- [libMesh](http://libmesh.sourceforge.net) v0.9.1+ compiled with [SLEPc](http://www.grycap.upv.es/slepc/)
- [HDF5](http://www.hdfgroup.org/HDF5/) v1.8.0+

License
-------

See [`LICENSE`](https://github.com/libqueso/queso/blob/dev/LICENSE) file
distributed with `QUESO` for more information.

Contributing
------------

Contributions are very welcome.  If you wish to contribute, please
take a few moments to review the [branching
model](http://nvie.com/posts/a-successful-git-branching-model/)
`QUESO` utilizes.

Support
-------

If you have questions or need help with using or contributing to `QUESO`,
feel free to ask questions on one of the mailing lists:

- [queso-users](https://groups.google.com/forum/#!forum/queso-users) mailing
  list for questions regarding usage and reporting bugs
- [queso-dev](https://groups.google.com/forum/#!forum/queso-dev) mailing list
  for discussion regarding development of `QUESO`

Citing QUESO
-------
Please add the following citation to any paper, technical report or
article describing the use of the `QUESO` library:

```bibtex
@inproceedings{Prudencio2012,
  author = {Prudencio, Ernesto E and Schulz, Karl W},
  booktitle = {Euro-Par 2011: Parallel Processing Workshops},
  pages = {398--407},
  publisher = {Springer},
  title = {{The parallel C++ statistical library ‘QUESO’: Quantification of
    Uncertainty for Estimation, Simulation and Optimization}},
  url = {http://dx.doi.org/10.1007/978-3-642-29737-3\_44},
  year = {2012}
}
```
