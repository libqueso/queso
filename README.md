The `QUESO` Library
=================

`QUESO` stands for Quantification of Uncertainty for Estimation,
Simulation and Optimization.

`QUESO` is a collection of algorithms and other functionalities aimed
for the solution of statiscal inverse problems, for the solution of
statistical forward problems, for the validation of a model and for
the prediction of quantities of interest from such model along with
the quantification of their uncertainties.

`QUESO` is designed for flexibility, portability, easiness of use and
easiness of extension. Its software design follows an object-oriented
approach and its code is written on C++ and over MPI. It can run over
uniprocessor or multiprocessor environments.

Installation
------------

If you do not have a `configure` script in the top level directory,
run `bootstrap` to generate a configure script using autotools.

Before compiling, you must run the `configure` script.  To run, type
`./configure`.  Additional options may be provided if desired.  Run
`./configure --help` for details.

After successfully running `configure`, type `make` to build the
`QUESO` library

Then type `make install` to install it in the directory previously
specified by the `--prefix` option of the `configure` script.

Dependencies
------------

`QUESO` has three required dependencies:

- [boost](http://www.boost.org/)
- mpi (either [mpich](http://www.mpich.org/) or
  [openmpi](http://www.open-mpi.org/))
- [gsl](https://www.gnu.org/software/gsl/)

`QUESO` also has several optional dependencies:

- [teuchos](http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html)
- [grvy](https://red.ices.utexas.edu/projects/hpct/files)
- [hdf5](http://www.hdfgroup.org/HDF5/)
- [glpk](https://www.gnu.org/software/glpk/)
- [petsc](http://www.mcs.anl.gov/petsc/)

License
-------

See `LICENSE` file distributed with `QUESO` for more information.

Contributing
------------

Contributions are very welcome.  If you wish to contribute, please take a few
moments to review the [branching
model](http://nvie.com/posts/a-successful-git-branching-model/) `QUESO`
utilises.

Support
-------

If you have questions or need help with using or contributing to `QUESO`,
feel free to ask questions on one of the mailing lists:

- [queso-users](https://groups.google.com/forum/#!forum/queso-users) mailing
  list for questions regarding usage and reporting bugs
- [queso-dev](https://groups.google.com/forum/#!forum/queso-dev) mailing list
  for discussion regarding development of `QUESO`
