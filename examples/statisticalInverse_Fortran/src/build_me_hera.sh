#!/bin/bash

mpif90 -c -g likelihood_function.f90

mpic++ -g -I./ \
    -I$QUESO_DIR/include \
    -I$TRILINOS_DIR/include \
    -I$BOOST_DIR/include \
    -I$HDF5_DIR/include  \
    -I$GLPK_DIR/include  \
    -DMPICH_IGNORE_CXX_SEEK \
    exStatisticalInverseProblem1_gsl.C likelihood_function.o -L$QUESO_DIR/lib -Wl,-rpath,$QUESO_DIR/lib -lqueso