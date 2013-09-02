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

The project homepage is [here](http://pecos.ices.utexas.edu).

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

License
-------

See `LICENSE` file distributed with `QUESO` for more information.

Contributing
------------

Contributions are very welcome.  If you wish to contribute, please take a few
moments to review the [branching
model](http://nvie.com/posts/a-successful-git-branching-model/) `QUESO`
utilises.
